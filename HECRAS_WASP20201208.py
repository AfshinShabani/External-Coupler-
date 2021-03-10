# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 11:58:31 2020

@author: afshin
"""

import os
import sys 
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
import pyproj
from pyproj import _datadir, datadir
from fiona import _shim,schema
import geopandas as gpd
import pandas as pd
import h5py
import numpy as np
from scipy import interpolate 
from random import randrange,uniform
import time 
from datetime import date
import itertools
from numpy import inf
from shapely.geometry import shape
from shapely.geometry import Polygon
from shapely.geometry import LineString,Point
from shapely.geometry.polygon import LinearRing
from shapely.geometry.multipolygon import MultiPolygon
import csv
#import psutil 
#pyproj.datadir.set_data_dir(r'C:\Users\afshin\AppData\Local\Temp\_MEI315602\pyproj\proj_dir\share\proj')

# Memory shortage warning 
# memory=psutil.virtual_memory()
# def memory_swamp():
#     if memory[2]>98:
#        print("Warning:program is running out of memroy on this machine (<2% of memory is availabel)")
#     if memory[2]>=99:
#        print("Memory error: the program was determined due to insffuicient memory")
#        sys.exis()

print ("HEC-RAS 2D and WASP External Coupler")
 
print ("Please enter directory for External Coupler configuration file")
Config_File=input()
with open(os.path.join(r'',Config_File,'Configuration.txt'),'r+') as address:
#with open(r'D:\Example\External Coupler\Configuration.txt','r+') as address:    
    lines=address.readlines()
    address.close()
PATH_HECRAS=os.path.join(r'',lines[14].rstrip())
PATH_WASP=os.path.join(r'',lines[15].rstrip())
HDF=os.path.join(r'',lines[16].rstrip())
PostProcess=os.path.join(r'',lines[17].rstrip())
PATH=os.path.join(r'',lines[18].rstrip())
outlet=os.path.join(r'',lines[19].rstrip())
# timestep=int(os.path.join(r'',lines[21].rstrip()))
MinVolume=float(lines[20].rstrip())
maxVolume=float(lines[21].rstrip())
Num_layer=int(lines[22].rstrip())
Sediment_depth=float(lines[23].rstrip())
HDT=int(lines[24].rstrip())
Year,Month,Day,Hour,Minute,Second=[int(i) for i in (lines[25].rstrip().split(','))]
Hydro_name=os.path.join(lines[26].rstrip())

    
print ("Reading HEC-RAS 2D output")

HRAS=h5py.File(HDF,'r') # Read HEC-RAS 2D HDF output 
Floodplain=HRAS['/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas']
Name=list(Floodplain.keys())[0]
Odir='/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/%s/' %Name # Directory for the HEC-RAS output file 
Gdir='/Geometry/2D Flow Areas/%s/' %Name                                                               # Directory for the HEC-RAS geometery file 
BoundDir='/Geometry/Boundary Condition Lines'                                                          # Directory for boundaries numbers in HEC-RAS model 
BoundFlow='/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/%s/Boundary Conditions/' %Name# Directory for boundaries's flow in HEC-RAS model 
Odir1='/Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/%s/Processed Data/Profile (Horizontal)/'%Name # directory for the HEC-RAS output file
Odir2=r'/Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/2D Flow Areas/%s/Boundary Conditions/' %Name

# HEC-RAS output 
Output=HRAS.get(Odir)                       
CDepth=np.array(Output.get('Depth'))             # Cell depth 
WSurfElev=np.array(Output.get('Water Surface'))  # Water surface elevation  
rows,columns=WSurfElev.shape                     # Number of cell and time step 
Flow=np.array(Output.get('Face Q'))              # Net flow at the cell face

#### HEC-RAS Geometery
Geo=HRAS.get(Gdir)
CelVEInfo=np.array(Geo.get('Cells Volume Elevation Info')) # Read Cell Volume Elevation Information 
CelVE=np.array(Geo.get('Cells Volume Elevation Values'))   # Read Cell Volume Elevation Information 
FAEInfo=np.array(Geo.get('Faces Area Elevation Info'))     # Read Face Area Elevation Information 
FCIndex=np.array(Geo.get('Faces Cell Indexes'))            # Read Faces and Their Corresponding Cell Number 
CArea=np.array(Geo.get('Cells Surface Area'))              # Read Cell Face Area 


## HEC-RAS Post process
FVelocity=np.array(Output.get('Face Velocity')) # Face velocity 
Post=h5py.File(PostProcess,'r')
Output=Post.get(Odir1)
Farea=np.array(Output.get('Face Hydraulic Area'))           # Read Face Area  
FDepth=np.array(Output.get('Face Hydraulic Depth'))   
# Flow=Farea*FVelocity

print ("Generating WASP segments")
#create WASP segments
WASP=gpd.read_file(PATH_WASP)                                                  #WASP preliminary segments                                                                   # Get coordinate system
HECRAS_CELL=gpd.read_file(PATH_HECRAS)
crr=HECRAS_CELL.crs  
WASP.crs=crr
def WASP_SEGMENT():
    if 'Segment' in WASP.columns:
        del WASP['Segment']
    
    if 'Cell Index' in WASP.columns:
        WASP_SEG=WASP.rename(columns={'Cell Index':'Segment'})
        WASP_SEG.Segment=WASP_SEG.index+1
    else:
        WASP.insert(0,'Segment',WASP.index+1)
        WASP_SEG=WASP
    HECRAS_POINTS=HECRAS_CELL.copy()                                               # create a copy to convert HEC-RAS cells to points
    HECRAS_POINTS['geometry']=HECRAS_POINTS['geometry'].centroid                   # convert 
    HEWASP_join = gpd.sjoin(WASP_SEG, HECRAS_POINTS,how='inner',op="intersects")   # Join HEC-RAS cell to WASP SEGMENT
    HEWASP=HEWASP_join[['Segment','Cell Index']]                                   # HEC-RAS and WASP segments 
    HEWASP1= HEWASP.drop_duplicates(subset='Cell Index', keep="first")             # Remove points that are sheared by WASP boundary lines
    HECW=HECRAS_CELL.merge(HEWASP1,on=['Cell Index'])
    HECW.FID=HECW['Segment']
    WSEGMENT=HECW.dissolve(by='FID',aggfunc='first')
    WSEGMENT=WSEGMENT.explode() # convert multipart polygon to single part 
    WSEGMENT.Segment=[x for x in range(1,len(WSEGMENT)+1)]
    WSEGMENT.index=[x for x in range(1,len(WSEGMENT)+1)] # reset WSEGMENT multipindex to single index 
    Segment_Cells=gpd.sjoin(WSEGMENT,HECRAS_POINTS,how='inner',op="intersects")
    DSC=Segment_Cells[['Segment','Cell Index_right']]
    DSC.columns=['Segment','CellIndex']
    WSEGMENT.Area=WSEGMENT.geometry.area
    WSEGMENT.to_crs(crr)
    MD=WSEGMENT
    WSEGMENT.to_file(os.path.join(PATH,'WASP-Segments.shp'))

# Creating flow path by converting segments to lines
    def polygon_polyline (poly,crr):
        Line_coordinate=[]  # coordinate for boundaries
        for rw in WSEGMENT.Segment:
            r=WSEGMENT[WSEGMENT.Segment==rw]
            lines=list(r.geometry.boundary.iloc[0].coords) # convert polygons to points
            for i in range(len(lines)-1):
                lcord=[lines[i],lines[i+1]]
                if not any(lcord) in Line_coordinate:
                    Line_coordinate.append(lcord)
        line=[] # Boundaries
        for ln in Line_coordinate:
            line.append(LineString(ln))
        polylin=gpd.GeoDataFrame({'geometry':line},crs=crr)
        polylin.insert(0,'name',[x for x in range(1,len(polylin)+1)]) 
        Cordlin=polylin.geometry.bounds
        Cordlin.insert(0,'name',polylin.name)
        Group=Cordlin.groupby(['minx','miny','maxx','maxy']).agg(['first','last','count']) # this will help to delete the duplicated lines 
        Subset=Group['name']['first'].tolist()
        polylin=polylin[polylin.name.isin(Subset)] # remove duplicated lines
        polylin.to_file(os.path.join(PATH,'Boundaries.shp'))
        return(line,polylin) 
    
    line,polyline=polygon_polyline(WSEGMENT,WASP.crs) # creat polyline and its shapefile

## Create HEC-RAS 2D and WASP Catalog File
    point=gpd.GeoDataFrame({'geometry':polyline.geometry.centroid.buffer(0.1),'name':polyline.name},crs=crr)
    #point.to_file(os.path.join(PATH,'points.shp'))
    Fmap=gpd.sjoin(point,WSEGMENT,how='left',op='intersects')                            # intersect flow path and WASP Segments
    Fmap=Fmap.rename(columns={'name':'OBJECTID'})
    #Fmap.to_file(os.path.join(PATH,'pointline.shp')) 
    COMSC=gpd.sjoin(WSEGMENT, HECRAS_POINTS,how='inner',op="contains")  # intersecting WASP and HEC-RAS computional girds
    HECCELL=pd.DataFrame({'Cell Index':HECRAS_CELL['Cell Index'],'geometry':HECRAS_CELL['geometry']})
    WSHEC=pd.DataFrame({'Segment':COMSC['Segment'],'Cell Index':COMSC['index_right']}) # creating dataframe for WASP segments and HEC-RAS Cells 
    CSJ=pd.merge(HECCELL,WSHEC,on='Cell Index')                                    # Join HEC-RAS Cell and WASP segment 
    HECB=gpd.GeoDataFrame({'geometry':CSJ['geometry'],'CellIndex':CSJ['Cell Index'],'Segment':CSJ['Segment']},crs=crr)
    Cdmap=gpd.sjoin(point,HECB,how='left',op='intersects')
    del Cdmap['index_right']
    Cdmap.columns=['geometry','OBJECTID','CellIndex','Segment']
    Cdmap.to_file(os.path.join(PATH,'HEC-RAS2D_WASP_catalog.shp')) # HEC-RAS and WASP flow path catalog file 
    return(MD,Fmap,Cdmap,DSC) 

MD,FMAP,Cdmap,DSC=WASP_SEGMENT()

print ("Calculating initial volume and depth for WASP segments")
def FC_Index():
    FIndex=pd.DataFrame(FCIndex,columns=['Cell1','Cell2'])   # This will trun the face cell indedx to the a dataframe 
    FIndex.insert(0,'Index',FIndex.index)
    return (FIndex)
FCIndex1=FC_Index()

# def DEPTH_VOLUME():
TVE=pd.DataFrame(columns=['Index','Volume','Elevation']) ### Dataframe that store volme and elevation of cells 
for i in range(columns):
    b=CelVEInfo[i,0]
    NRow= CelVEInfo[i,1:2][0] # number of data for each cell 
    VE=CelVE[b:b+NRow,:]
    indx=[]
    for rw in range(NRow):
        indx.append(i)
    dataset=pd.DataFrame({'Index':indx,'Elevation':VE[:,0],'Volume':VE[:,1]})
    TVE=pd.concat([TVE,dataset],ignore_index=True,sort=True) 
CellVolume=[]  # Volume for cells within a segment
CellDepth=[]   # Depth for cels within a segment
CellArea=[]  
Cell_Index=[]
Cn=0
for j in DSC.CellIndex:
    TVS=TVE[TVE.Index==int(j)]
    Cell_Index.append(j)
    for i in range (0,1): #rows
        CE=WSurfElev[i,j]
        CellArea.append(CArea[j])         # Get the cell area
        CellDepth.append(CDepth[i,j])     # Get cell depth 
        if len(TVS)>1:                    # calcualte volume for non-flat cells 
            Extrapolate=interpolate.interp1d(TVS.Elevation,TVS.Volume, fill_value='extrapolate') # this code will interaploate-extrapolate the cell volume 
            CV=max(Extrapolate(CE),0) # to make sure the negative value are not included in cell volume  
            CellVolume.append(np.round(CV,4))
        else:                            # Volume for flat cells is calculated by depth*Area
            CellVolume.append(np.round(CDepth[i,j]*CArea[j],4))
    for i in range (0,rows): #rows
        CellArea.append(CArea[j])         # Get cell area
        CellDepth.append(CDepth[i,j])     # Get water depth 
    Cn=Cn+1
def DEPTH_VOLUME():
    CellV=pd.DataFrame(np.reshape(np.asarray(CellVolume),(len(DSC),len(CellVolume)//len(DSC))))
    CellV.index=Cell_Index
    DSCV=pd.merge(DSC,CellV,left_on='CellIndex', right_on=CellV.index) # Join cell volume to segments 
    SegVolume=DSCV.groupby(['Segment']).sum()
    SegmentVolume=SegVolume.iloc[:,1:].values  # This value should be 2 
    CellD=pd.DataFrame((np.reshape(np.asarray(CellDepth),(len(DSC),len(CellDepth)//len(DSC)))))
    CellA=pd.DataFrame((np.reshape(np.asarray(CellArea),(len(DSC),len(CellArea)//len(DSC)))))
    multiply=CellD*CellA  # cell area multiply by depth 
    multiply['CellIndex']=Cell_Index   # insert cell index 
    multiply['CellArea']=CellA.iloc[:,0]
    Smultiply=pd.merge(DSC,multiply,left_on='CellIndex', right_on='CellIndex')
    DASum=Smultiply.groupby(['Segment']).sum()
    WDASum=DASum.div(DASum['CellArea'],axis=0)
    SegmentDepth=WDASum.iloc[:,1:-1].values
    # CellA['Segment']=DSC['Segment'].values
    # multiply=CellD*CellA
    # multiply['Segment']=DSC['Segment'].values
    # DASum=multiply.groupby(['Segment']).sum()
    # Asum=CellA.groupby(['Segment']).sum()
    # SegmentDepth=(DASum/Asum).values
    return(SegmentVolume,SegmentDepth)   
SegmentVolume,SegmentDepth=DEPTH_VOLUME()

#Delete some of data to clean up the memory 
# del CelVEInfo
# del CelVE
# del FAEInfo
# del CDepth
# del WSurfElev


print ("Calculating flow for WASP segments")

def SEGMENT_FLOW():
    Bound=HRAS.get(BoundDir)
    BouNames=np.array(Bound.get('Attributes'))
    BEXNames=[]                                               # Boundary External names
    BINNames=[]                                               # Boundary Internal names
    BEXID=[]                                                  # Boundary External ID
    BINID=[]                                                  # Boundary Internal ID
    for i in range(0,len(BouNames)):
        Name=BouNames[i][0]
        if outlet in BouNames[i][0].decode('utf-8'):  #Modified
           BEXNames.append(Name)
           BEXID.append(i)
        else:
           BINNames.append(Name) 
           BINID.append(i) 
    
    ExternalFaces=np.array(Bound.get('External Faces'))        # Downstream Boundary Faces
    Bounflow=HRAS.get(Odir2)
    DBFlow= np.array(Bounflow.get('%s - Flow'%BEXNames[0].decode('utf-8')))                                              # Upstream boundary name and its flow valeus are stored in this dictionary 
    
    def Boundary_Flow(BINNames,BINID,FCIDndex1):
        UPB_dict={} 
        for name in BINNames:
            UPB_dict.update({'%s'%name:np.array(Bounflow.get('%s - Flow' %name.decode('utf-8')))})
        DExternalFaces=pd.DataFrame(ExternalFaces)
        UF=DExternalFaces[DExternalFaces['BC Line ID'].isin(BINID)]
        Bfaces=UF['Face Index'].tolist()
        #UF=ExternalFaces[ExternalFaces['BC Line ID']==BINID]   #  faces # Modified 
        #Bfaces=[UF[i][1] for i in range(len(UF))]              # select Boundary faces    
        BSEGMENT=[]                                   
        for j in Bfaces:
            SB=FCIndex1[FCIndex1.Index==j]
            Cell=SB.Cell1
            BSEG=DSC[DSC.CellIndex==int(Cell)]
            BSEGMENT.append(int(BSEG.Segment))           #Modified Name
        BUNSEGMENT=pd.DataFrame({'Faces':Bfaces,'Segment':BSEGMENT})
        SEGNT=[]                                                   # name of the segment 
        SEGOUT=[]                                                  # Upstream boundary or 0  
        LENSEG=[]                                                   # Get lenght of cell in the segment           
        for i in BUNSEGMENT.Segment.unique():                      
            SEGNT.append(i)
            SEGOUT.append(0)
            SBUNSEGMENT=BUNSEGMENT[BUNSEGMENT.Segment==i]
            LENSEG.append(len(SBUNSEGMENT))
        SEGFLOW=[]                                                # outflow from boundary segments
        UPB=pd.DataFrame()
        for key in UPB_dict.keys():
            DUPB=pd.DataFrame.from_dict(UPB_dict[key])
            UPB=pd.concat([UPB,DUPB],axis=1)
        UBFlow=UPB.to_numpy()
        for index,j in enumerate(LENSEG):
            r=int(np.sum(LENSEG[:index]))
            SUBFlow=list(np.sum(UBFlow[:,r:r+j],axis=1))
            SEGFLOW=SEGFLOW+SUBFlow
        return(SEGNT,SEGOUT,SEGFLOW)
    
    UPSEGN,UUPSEGN,UPSEGFLOW=Boundary_Flow(BINNames,BINID,FCIndex1) # upstream flow 
    BNSEGMENT,BSEGOUT,BSEGFlow=Boundary_Flow(BEXNames,BEXID,FCIndex1) # upstream flow 
    
    Fgroup=Cdmap.groupby(['OBJECTID']).count() # group the faces with number of thier cells 
    Direction=Fgroup[Fgroup.Segment>1]         # Face between two segments 
    SFgroup=Cdmap[Cdmap['OBJECTID'].isin(Direction.index)]
    From=SFgroup.groupby('OBJECTID').first()  # get first cell for a segment 
    To=SFgroup.groupby('OBJECTID').last()     # get second cell for a segment 
    From_To=pd.merge(From,To,on='OBJECTID')   # create flow direction 
    From_To=From_To[['CellIndex_x','CellIndex_y','Segment_x','Segment_y']] # join boundary, cell, and segment names
    From_To.columns=['Cell1','Cell2','From','To'] # segment boundaries and corresponding HEC-RAS cells 
    From_To.columns=['Cell1','Cell2','From','To']
    To_From=From_To.copy() # inverse direciton 
    To_From.columns=['Cell2','Cell1','From','To']
    DFlow=pd.DataFrame(np.transpose(Flow)) # create a dataframe for HEC-RAS flow
    Inverse_DFlow=pd.DataFrame(np.transpose(Flow*-1))     # Inverse flow 
    Face_Flow=pd.merge(FCIndex1,DFlow,left_on='Index',right_on=DFlow.index) #face flow for direction 
    Face_Flow_Inve=pd.merge(FCIndex1,Inverse_DFlow,left_on='Index',right_on=Inverse_DFlow.index) # face flow for opposite direction
    Sub1=pd.merge(From_To,Face_Flow,left_on=['Cell1','Cell2'],right_on=['Cell1','Cell2'])
    Sub1=Sub1.drop_duplicates()
    Sub2=pd.merge(To_From,Face_Flow_Inve,left_on=['Cell1','Cell2'],right_on=['Cell1','Cell2']) 
    Sub2=Sub2.drop_duplicates()
    Sub1_Sub2=pd.concat([Sub1,Sub2])  # combine both flow datasets 
    Total_Flow=(Sub1_Sub2.drop(['Cell1','Cell2','Index'],axis=1)).groupby(['From','To']).sum()
    FlowDirection=np.array(list(Total_Flow.index))
    TFlow=Total_Flow.to_numpy()
    SegFlow=np.concatenate([FlowDirection,TFlow],axis=1) 
    
    ###remove duplicate flow path
    SegDFlow=pd.DataFrame(SegFlow,columns=[str(x) for x in range(rows+2)])
    t=pd.DataFrame(SegFlow[:,:2],columns=['0','1'])
    lst=[]
    for i in t.iterrows():
        a=(i[1][0],i[1][1])
        b=(i[1][1],i[1][0])
        if a and b not in lst:
            lst.append(a)
    t1=pd.DataFrame(lst,columns=['0','1'])
    t2=pd.merge(t1,SegDFlow,left_on=['0','1'],right_on=['0','1'])
    SegFlow=t2.to_numpy()
    ##
    BUpstream=np.array(list(zip(UUPSEGN,UPSEGN)))
    DUpstream=np.asarray(UPSEGFLOW).reshape(len(UUPSEGN),len(UPSEGFLOW)//len(UPSEGN)) # fix this rows
    Upstream=np.concatenate([BUpstream,DUpstream],axis=1)
    BDownstream=np.array(list(zip(BNSEGMENT,BSEGOUT)))
    DDownstream=np.asarray(BSEGFlow).reshape(len(BNSEGMENT),len(BSEGFlow)//len(BNSEGMENT)) # fix this rows 
    Downstream=np.concatenate([BDownstream,DDownstream],axis=1)
    TFlow=np.concatenate([Upstream,SegFlow,Downstream])
    FlowDirection=TFlow[:,0:2] # direction of flow 
    FLOW=TFlow[:,2:]
    np.savetxt(os.path.join(PATH,'Flow.txt'),TFlow) # Save flow and direction for segments 
    np.savetxt(os.path.join(PATH,'Direction1.txt'),FlowDirection)
    # convert instantaneous flow to integrated flow 
    T=[]
    timestep=2
    for i in range (0,rows,timestep):
        if i < rows-timestep:
            a=i+timestep
            for j in range(0,timestep):
                T.append(((FLOW[:,i]+FLOW[:,i+timestep])/2).tolist())
        else:
            for d in range(0,1):
                T.append(((FLOW[:,rows-timestep]+FLOW[:,-1])/2).tolist())
    FLOW=np.transpose(np.array(T))
    DFLOW=np.concatenate([FlowDirection,FLOW],axis=1)
    return(DFLOW,FLOW,FlowDirection)
DFLOW,FLOW,FlowDirection=SEGMENT_FLOW()
#DFLOW1=DFLOW

print("Calculating volume for WASP segments")


# def SEGMENT_VOLUME(DFLOW):
FlowD=pd.DataFrame(DFLOW)    
SVolume=[]
Segment=[]
minvl=[]
TimeInd=[]
BoundaryInd=[]
# MinVolume=(MD.Area*0.1).tolist()
for index,name in MD.iterrows(): # Put number of segments here
    i=int(name['Segment'])
    index=i-1
    FlowD1=FlowD[(FlowD.iloc[:,0]==i)]*-1  # multiple by negaive -1 because outflow decreases the volume
    FlowD2=FlowD[(FlowD.iloc[:,1]==i)]
    FlowD3=pd.concat([FlowD1,FlowD2])
    # FD1C=FlowD3.cumsum(axis=1)
    # FD2C=FD1C.iloc[:,3:671]
    FlowD3.insert(0,'Segment',i)
    FlowD4=FlowD3.groupby(['Segment']).sum()
    FlowD4=FlowD4.drop([0,1],axis=1)
    FlowD4.columns=[str(x) for x in range(len(FlowD4.columns))]
    FlowD5=FlowD4*60
    if SegmentVolume[index,0]<MinVolume:    
        FlowD5['0']=int(MinVolume)
    else:
        FlowD5['0']=SegmentVolume[index,MinVolume]    
    FlowD6=FlowD5.cumsum(axis=1)
    FlowD7=FlowD6.T.loc[FlowD6.T[i]<=0]  # subset data were flow become negative
    if len(FlowD7)>0:
        print(i)
        INDF=int(FlowD7.index[0])                        # get index of negaitve volume
        INP=[]
        #INPvalue=[]
        INPID=FlowD6.iloc[:,INDF-1:]                      
        FlowD8=FlowD6.iloc[0,INDF-1]                       # Get location of first value < minvolume
        #FlowD8=MinVolume
        for NAME, values in INPID.iteritems():             # subset negative volume
            #FlowD9=MinVolume+values.tolist()[0]
            FlowD9=values.tolist()[0]
            if FlowD9<MinVolume:                        # add flow index to the list if the volume is < minvolume
               INP.append(int(NAME))
               #INPvalue.append(0)
            else:
                FlowD8=FlowD8+values.tolist()[0]           # update the minimum volume for segment
        FlowSHDW=FlowD3.index.tolist()                      # Index for boudnaries to make thier flow zero
        INDEX=np.argmin(FlowD3.iloc[:,INDF])
        # for k in FlowSHDW[INDEX]:
        K=FlowSHDW[INDEX]                                 # find the largest negative value
        for l in INP:
            FlowD.iloc[k,l+2]=0
        break
    return(FlowD)            
FLOWD=SEGMENT_VOLUME(DFLOW)           
FL=SEGMENT_VOLUME(FLOWD)            
                
    V=FlowD6.iloc[0,:].tolist()
    if np.min(V)<MinVolume and np.min(V)>(MinVolume+maxVolume)*-1:  
        FlowD5['0']=FlowD5['0']+abs(np.min(V))+MinVolume+maxVolume
        FlowD6=FlowD5.cumsum(axis=1)
        V=FlowD6.iloc[0,:].tolist() 
    if np.min(V)<(MinVolume+maxVolume)*-1:         
       Segment.append(i)
       minvl.append(np.min(V))
    SVolume.append(V)
LSVolume=list(itertools.chain.from_iterable(SVolume))  
SegVolume=np.reshape(np.asarray(LSVolume),(len(MD),rows))
if len(Segment)>0:
    Modify_Segments=pd.DataFrame({'Segment':Segment,'mminimum volume': minvl})
    #Modify_Segments.to_csv(os.path.join(PATH,'Segments.csv'),chunksize=1000,mode='a')
    Modify_Segments.to_csv(os.path.join(PATH,'Segments.csv'))
    # print ('Error: Please revise following segments in segments.csv\
    #        \n or increase the maximum initial volume for dry segments. \
    #            \n The maximm initial volume must be > %s' %(np.min(minvl)*-1))
    return(SegVolume,Segment,minvl,FlowD)

SegVolume,Segment,minvl,FlowD=SEGMENT_VOLUME()
def terminate ():
     if len(Segment)>0:
         # print ('Error: Please revise segments in segments.csv')
         print ('Error: Please set additional volume for wet/dry segmets to %s' %(np.min(minvl)*-1))
         time.sleep(3)
         sys.exis()

if len (HECRAS_CELL)> len (MD):
    Itr=1
    while len(Segment)>0 and Itr<30:
        if Itr==1:
            print('Volume Error: %s \n program is modifying WASP segments\n this may take several minutes' %(np.min(minvl)))
        # print('Iteration number:%s'%Itr) 
        print ('Volume Error: %s'%np.min(minvl))
        print('Modifying WASP segment numbers:',Segment)  
        Nsegment=len(Segment)
        WD=DSC[DSC.Segment.isin(Segment)] # Wet\Dry segments
        WDCells=list(WD.CellIndex.unique())
        WDSegments=list(WD.Segment.unique())
        HEC_Sub=HECRAS_CELL[HECRAS_CELL['Cell Index'].isin(WDCells)] # select HEC-RAS cells for Wet/Dry segments
        SUB_WSEG=MD[~MD.Segment.isin(WDSegments)] # Excluding Wet/Dry Segment
        New_Segment=SUB_WSEG.append(HEC_Sub)
        WASP=New_Segment['geometry']
        WASP=gpd.GeoDataFrame(WASP)
        print('Regenerating WASP segments')
        MD,FMAP,Cdmap,DSC=WASP_SEGMENT()
        print('Recalculating inital volume and depth')
        SegmentVolume,SegmentDepth=DEPTH_VOLUME()
        print("Recalculating flow")
        DFLOW,FLOW,FlowDirection=SEGMENT_FLOW()
        print ("Recalculating volume")
        SegVolume,Segment,minvl,FlowD=SEGMENT_VOLUME()
        if Nsegment==len(Segment):
           NSegment=[]
           for i in Segment:
               Rflow=FlowD[(FlowD.iloc[:,0]==i)|(FlowD.iloc[:,1]==i)]
               NSegment.extend([x for x in (Rflow.iloc[:,0].unique().tolist())]) # add neighbors segment for revision 
               NSegment.extend([x for x in (Rflow.iloc[:,1].unique().tolist())])
           # Segment.extend([x for x in NSegment if x not in Segment])
           Segment=NSegment
        Itr=Itr+1          
    terminate()
    print("Volume error was resolved")
terminate()
    
SSVolume=(MD.Area*Sediment_depth).tolist()


print ("Calculating segment velocity")
### HEC-RAS post process file 


def Cell_Velocity(FCIndex1,FVA): # function to join face velocity and face area for cells 
    CFVA=pd.merge(FCIndex1,FVA,left_on='Index',right_on=FVA.index)  # Join VA to the cells 
    Cindexes=pd.concat([CFVA.groupby('Cell1').sum(),CFVA.groupby('Cell2').sum()],sort=True) # Join cells to their faces
    CFVAI=Cindexes.groupby(Cindexes.index).sum()
    CFVAI=CFVAI.iloc[:,2:rows+2]#rows
    return(CFVAI)

FA=pd.DataFrame(np.transpose(abs(Farea)))
FVA=pd.DataFrame(np.transpose(abs(FVelocity*Farea)))       # velocity * hydraulic area for faces
CellVelocity=pd.DataFrame(np.nan_to_num(Cell_Velocity(FCIndex1,FVA)/Cell_Velocity(FCIndex1,FA)))#*np.reshape(CArea,(columns,1)
Cell_Seg=pd.merge(CellVelocity,DSC,left_on=CellVelocity.index,right_on='CellIndex')  # Join cells velocities to segments
DCArea=pd.DataFrame(CArea,columns=['Area'])
Cell_Seg=pd.merge(Cell_Seg,DCArea,left_on='CellIndex',right_on=DCArea.index)
CVA=Cell_Seg.mul(Cell_Seg['Area'],axis=0)
CVA.Segment=Cell_Seg.Segment
SVelocity=CVA.groupby(['Segment']).sum()
SegmentVelocity=SVelocity.iloc[:,0:rows].to_numpy()/np.reshape(np.asarray(MD.Area),(len(MD),1)) # calculate segment veolcity 

# Write hydrodynamic data to output 
np.savetxt(os.path.join(PATH,'Velocity.txt'),SegmentVelocity)
np.savetxt(os.path.join(PATH,'Volume.txt'),SegVolume)
np.savetxt(os.path.join(PATH,'Segment_Depth.txt'),SegmentDepth)    # save segment depth 
np.savetxt(os.path.join(PATH,'Flow.txt'),FlowD)


print ("Generating input file for WASP Hydrolink API")
#Hydro_output=r'C:\Users\afshin\Documents\WASP20200627\WASP20200620\wasp\bin'
# Hydrolink_name='Hydrodynamic.txt'      # Name 
Input=open(os.path.join(PATH,'%s.txt'%Hydro_name),'w+')
Nseg=len(SegVolume)*Num_layer   # number of segments
NFLOW=len(FLOW)       # flow 
From=FlowDirection[:,0] # Flow direction from segments 
To=FlowDirection[:,1]    # flow direction to segments
Drow,Dcolumns=SegVolume.shape
Frow,Fcolumns=FLOW.shape
Header=[Nseg,NFLOW,HDT,Month,Day,Year,Hour,Minute,Num_layer] # inforamtion that will be saved in header file and will be read by the WASP Hydrolink API
for i in Header:
  Input.write(str(i)+',')  
Input.write('\n')  
if Num_layer==1:
    for j in range(len(From)):
        Input.write(str(int(From[j]))+','+str(int(To[j]))+'\n') 
    for r in range(0,Dcolumns):
        for d in range(0,Drow):
            Input.write(str(SegVolume[d,r])+','+str(max(SegmentDepth[d,r],0.05))+','+str( SegmentVelocity[d,r])+str('\n'))
        for m in range(0,Frow):
            Input.write(str(FLOW[m,r]))
            Input.write('\n')    
if Num_layer==2:
    for j in range(len(From)):
        Input.write(str(int(From[j]))+','+str(int(To[j]))+'\n') 
    for r in range(0,Dcolumns):
        for d in range(0,Drow):
            Input.write(str(SegVolume[d,r])+','+str(max(SegmentDepth[d,r],0.05))+','+str( SegmentVelocity[d,r])+str('\n'))
        for d in range(0,Drow):
            Input.write(str(SSVolume[d])+','+str(Sediment_depth)+','+str(0)+','+'\n') 
        for m in range(0,Frow):
            Input.write(str(FLOW[m,r]))
            Input.write('\n')
Input.close()        
  
SegName=[str(x+1) for x in range (0,(len(SegVolume))*2)]
with open (os.path.join(PATH,'Segment-name.txt'),'w') as SNAME:
    writer=csv.writer(SNAME,delimiter='\n')
    writer.writerow(SegName)
        

with open(os.path.join(PATH,'hydrolink.ctl'),'w+') as Ctl: # creating new output
     Ctl.write('a'+'\n') 
     Ctl.write('H'+'\n')
     Ctl.write('%s.txt'%Hydro_name+'\n')
     Ctl.write('%s2D.hyd'%Hydro_name+'\n')
     Ctl.write('   %s'%Day+'\n')
     Ctl.write('   %s'%Month+'\n')
     Ctl.write('   %s'%Year+'\n')
     Ctl.write('   %s'%Hour+'\n')
     Ctl.write('   %s'%Minute+'\n')
     Ctl.write('   %s'%Second+'\n')
     Ctl.write('Segment-name.txt'+'\n')
     Ctl.write('   0')
Ctl.close()

print ("The Hydrolink API hydrodynamic input file was sucessfully generated")








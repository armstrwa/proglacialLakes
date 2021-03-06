# -*- coding: utf-8 -*-
"""
Script to convert Hannah's GEEDIT lake polyline outlines to polygons

Created on Wed Mar 27 10:03:04 2019

@author: armstrongwh
"""

## IMPORT MODULES

#import modules - the order you do this in is very important due to bugs in shapely and fiona. Make sure you import in this order.
from shapely import geometry
from shapely.ops import transform
import fiona
from osgeo import ogr, osr, gdal
from fiona.crs import from_epsg
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from pyproj import Proj, transform
import shutil
import utm

## DEFINING FUNCTIONS

def readLineCoordsFromMultiLineShp(shapefileFn,transNum):
    '''
    	# TAKEN FROM M. FAHNESTOCK'S L8_sample_frames... code. 08 sep 2016
    	# Unfinished as of 08 sep 2016 wha
    	# 22 sep 2016 think it's finished wha
    	# Mark's notes	
    	# use fiona to read in points along centerline in local projection (called PS here)
    	# save the PS original points in sample_pts_PS, then reproject to 4326, push back into
    	# coordinates of input Point object (f) and make a list of these (out_pts_ll) and when done save
    	# them as a Shapely MultiPoint object
    	'''
    with fiona.open(shapefileFn,'r') as c:
    	
        # Initialize
        sample_pts_local=[]
        sample_pts_lon_lat=[]
        xList = []
        yList = []
        lonList = []
        latList = []		
        	
        # Get coordinate reference system from original shapefile
        original = Proj(c.crs)
        
        destination = Proj(init='EPSG:4326') # dest is WGS84 = EPSG 4326
        
        f = c[transNum] # open specified feature (transect)
        props = f['properties']

        for j in f['geometry']['coordinates']:
            x = j[0]
            y = j[1]
            xList.append(j[0])
            yList.append(j[1])
            sample_pts_local.append( (x,y) )
            lon,lat = transform(original, destination,x,y)
            lonList.append(lon)
            latList.append(lat)
            sample_pts_lon_lat.append( (lon,lat) )
                
        medLon = np.median(lonList) # get meadian longitude
        utmZone = utm.from_latlon(51.2,medLon)[2] # get utm zone

        epsgLocalUtm = 32600+utmZone
        destination2 = Proj(init='EPSG:'+str(epsgLocalUtm))
        easting,northing = transform(destination,destination2,lonList,latList)
        # Stick it together in a dict for storage
        dataOut = {'inProps':props,'shpFn':shapefileFn,'xLocal':xList,'yLocal':yList,'lat':latList,'lon':lonList,'samplePtsXy':sample_pts_local,'samplePtsLatLon':sample_pts_lon_lat,'utmEasting':easting,'utmNorthing':northing,'utmZone':utmZone}
        		
        return dataOut

def makePolygonShpFromDictList(dictList,shpFnOut):
    '''
    	# Function uses shapely and fiona to make a polygon from readLineCoordsFromMultiLineShp dict list
    '''

    nDicts = len(dictList) 
    polyList = []
    uPolyList = [] # doing utm polygons to calculate physically meangingful area
     
    for j in range(0,nDicts):
        dictNow = dictList[j]
        x = dictNow['lon']
        y = dictNow['lat']
        e = dictNow['utmEasting']
        n = dictNow['utmNorthing']
        
        coordPairs = []
        utmPairs = []
        
        for i in range(0,len(x)):
            coordPairs.append( (x[i],y[i]) )
            utmPairs.append( (e[i],n[i]) )
            
        polyList.append(geometry.Polygon(coordPairs))
        uPolyList.append(geometry.Polygon(utmPairs))
    
    mPoly = geometry.MultiPolygon(polyList)
     
    # Define a polygon feature geometry with one attribute
    schema = {
    'geometry': 'Polygon',
    'properties': {'id': 'int',
                   'imageDate':'str',
                   'imagePath':'str',
                   'shpName':'str',
                   'area_m2':'float',
                   'cent_lon':'float',
                   'cent_lat':'float'
                   }
    }

    # Write a new Shapefile
    with fiona.open(shpFnOut, 'w', 'ESRI Shapefile', schema,crs=from_epsg(4326)) as c:
    ## If there are multiple geometries, put the "for" loop here
        for k in range(0,len(mPoly)):
            dictNow = dictList[k]
            shpFnNow = dictNow['shpFn']
            inPropsNow = dictNow['inProps']
            imagePathNow = inPropsNow['ImagePath']
            dateNow = inPropsNow['Date']
            # polygon area in m^2
            areaNow = uPolyList[k].area
            # polygon centroid in lon, lat
            centNow = polyList[k].centroid
            cX = centNow.x
            cY = centNow.y
            
            c.write({
                    'geometry': geometry.mapping(mPoly[k]),
                    'properties': {'id':k,'imageDate':dateNow,'imagePath':imagePathNow,'shpName':shpFnNow,'area_m2':areaNow,'cent_lon':cX,'cent_lat':cY},
                    })
    		
    print "Sucessfully wrote shapefile: " + shpFnOut

def makeSinglePolygonFromSingleLine(lineGeom):
    '''
    This is incredibly simple at barebones. Input needs to be a shapely line geometry
    '''
    geom = geometry.asShape(lineGeom)
    poly = geometry.asPolygon(geom)
    
    return poly

def getPolygonVertexCoords(polyGeom):
    '''
    get vertex coordinates from a shapely polygon geometry
    returns xList, yList lists of vertex locations
    '''
    
    ext = polyGeom.exterior
    coords = ext.coords.xy
    xArr = coords[0]
    yArr = coords[1]
    xList = xArr.tolist()
    yList = yArr.tolist()    
    
    return xList, yList

def makePolyGeomFromCoordLists(xList,yList):
    '''
    make a shapely polygon geometry from a list of x (lon or easting) and y (lat or nothing) vertex coordinates
    outputs polyOut, a shapely polygon geometry
    '''
    
    coordPairs = []
    
    # make list of coordinate pair tuples
    for i in range(0,len(xList)):   
        coordPairs.append( (xList[i],yList[i]) )
    
    polyOut = geometry.Polygon(coordPairs)
    
    return polyOut

def getPolyAream2(polyGeom):
    '''
    This assumes the polygon has lat/lon coordinates (epsg 4326, wgs84)
    input = shapely polygon geometry with lat/lon coordinates
    output = polygon area in m^2
    '''
    
    src_crs = Proj(init='EPSG:4326') # get wgs84 projection (assumes this input geometry projection)
    
    # get centroid
    c = polyGeom.centroid

        
    utmZone = utm.from_latlon(c.y,c.x)[2] # get utm zone
    epsgLocalUtm = 32600+utmZone # epsg code for utm zone             
    dst_crs = Proj(init='EPSG:'+str(epsgLocalUtm))
    
    lonList, latList = getPolygonVertexCoords(polyGeom) # get vertex coordinates in lat/lon
    easting,northing = transform(src_crs,dst_crs,lonList,latList) # transform these coordiantes to local utm
    
    utmPoly = makePolyGeomFromCoordLists(easting,northing) # make polygon in utm coordinate system
    
    area_km2 = utmPoly.area # polygon area in m^2
    
    return area_km2

def fixTupleList(allLakeCentList):
    '''
    clean up function to deal with fact that I can't get the outline lists to work quite right
    '''    
    
    listNum = len(allLakeCentList) # how many entries are in this list?
    
    fixedList = [] # empty list for filling
    
    for i in range(0,listNum):
        entryNow = allLakeCentList[i]
        
        if type(entryNow) is tuple:
            fixedList.append(entryNow)
        elif type(entryNow) is list:
            subListNum = len(entryNow)
            
            for j in range(0,subListNum):
                fixedList.append(entryNow[j])
    
    return fixedList
            
# Function to sample raster at specified coordinates
# Returns values of raster at x,y
def sampleRasterAtXY(rasterIn,x,y):
	
    z = []
	
    xyTfm = createTransformation(rasterIn)
    imageIj = xyToImageIj(x,y,xyTfm)
	
    rast = gdal.Open(rasterIn)
    band = rast.GetRasterBand(1) # assumes single band raster
    bandArr = band.ReadAsArray() # convert band into array
	
    samplePtNum = int(imageIj.shape[0]) # number of sample points
	
    for i in range(0,samplePtNum):
        z.append( bandArr[ int(np.round(imageIj[i,1])), int(np.round(imageIj[i,0] )) ] )
	
    return z

# Function to convert pixel location to ground location
def imageToXy(i,j,xyTfm):
	
	outX = []
	outY = []
	
	it = np.nditer([i,j])
	
	for ii,jj in it:
		outX.append( xyTfm[0] + ( (ii+0.5) * xyTfm[2] ) )
		outY.append( xyTfm[1] + ( (jj+0.5) * xyTfm[3] ) )
	
	outArray = np.array((outX,outY))
	outArrayTranspose = outArray.transpose()
	return outArrayTranspose

# Function to convert ground coordinates to image pixel location
# Returns i,j pairs. I.e., outArrayTranspose[:,0] = i (col) coord,
# outArrayTransport[:,1] = j (row) coord
def xyToImageIj(x,y,xyTfm):
	
	outI = []
	outJ = []
	  
	it = np.nditer([x,y])
	
	for xx,yy in it:
		outI.append( (xx - xyTfm[0])/xyTfm[2] - 0.5 )
		outJ.append( (yy - xyTfm[1])/xyTfm[3] - 0.5 )
	
	outArray = np.array((outI,outJ))
	outArrayTranspose = outArray.transpose()
	return outArrayTranspose

# Creates transformation parameters for future work
# Outputs xyTfm vector that is [minX,maxY,pxSizeX,pxSizeY]
def createTransformation(rasterIn):
	rast = gdal.Open(rasterIn) # open raster
	gt = rast.GetGeoTransform() # get image to ground geotransformation params
	numPxX = rast.RasterXSize # number of pixels in x direction
	numPxY = rast.RasterYSize # number of pixels in y direction
	pxSizeX = gt[1] # pixel size in x
	pxSizeY = gt[5] # pixel size in y

	# upper left coordinates
	minX = gt[0] # minimum x value
	maxY = gt[3] # maximum y value

	maxX = minX + numPxX*pxSizeX # maximum x value
	minY = maxY + numPxY*pxSizeY # minimum y value

	xyTfm = np.array([minX,maxY,pxSizeX,pxSizeY]).astype('float')

	return xyTfm
	    

## USER INPUTS
    
# identify folder with shapefiles to merge
folderPath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\allLakes_25apr2019'
folderPathOut = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\allLakes_25apr2019'

fullFnCentOut = folderPathOut + '\\allLakeCentroids.shp'
fullFnTableOut = folderPathOut + '\\allLakeCentroids.csv'
fullFnPolyOut = folderPathOut + '\\allLakePolygons.shp'

# raster with elevation data
demFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\dem\\dem_wgs84.tif'

# rasters with winter precipiation
pptOldFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\pr_decadal_summaries_AK_CAN_2km_CRU_TS31_historical\\decadal_mean\\pr_decadal_mean_DJF_total_mm_cru_TS31_historical_1960_1969.tif'
pptMidFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\pr_decadal_summaries_AK_CAN_2km_CRU_TS31_historical\\decadal_mean\\pr_decadal_mean_DJF_total_mm_cru_TS31_historical_1980_1989.tif'
pptNewFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\pr_decadal_summaries_AK_CAN_2km_CRU_TS31_historical\\decadal_mean\\pr_decadal_mean_DJF_total_mm_cru_TS31_historical_2000_2009.tif'
pptChangeFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\pr_diff_djf_2000s-1960s.tif'

# rasters with summer air temperature
tasOldFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\tas_decadal_summaries_AK_CAN_2km_CRU_TS31_historical\\decadal_mean\\tas_decadal_mean_JJA_mean_c_cru_TS31_historical_1960_1969.tif'
tasMidFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\tas_decadal_summaries_AK_CAN_2km_CRU_TS31_historical\\decadal_mean\\tas_decadal_mean_JJA_mean_c_cru_TS31_historical_1980_1989.tif'
tasNewFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\tas_decadal_summaries_AK_CAN_2km_CRU_TS31_historical\\decadal_mean\\tas_decadal_mean_JJA_mean_c_cru_TS31_historical_2000_2009.tif'
tasChangeFn = 'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\climate\\tas_diff_jja_2000s-1960s.tif'


## PROCESSING

# delete output files if they already exist
if  os.path.exists(fullFnCentOut) or os.path.exists(fullFnCentOut):
    os.remove(fullFnCentOut)   
    os.remove(fullFnTableOut)

shpList = glob.glob(folderPath + '\\T*.shp') # this should produce a list of filenames for lake outline shapefiles

# make empty storage containers
allLakeDict = {}

shpIterCount = 1 # iteration counter

for fn in shpList:
    
    
    shpNow = fiona.open(fn) # open current lake (all outlines)
    
    shortFn = str(os.path.split(fn)[1]) # get shapefile name without path
    lakeCode = shortFn[5:-4] # isolate name code name
    
    print('Processing outline ' + str(shpIterCount) + ' of ' + str(len(shpList)) + ': ' + shortFn)
    
    thisLakeDict = {}
    thisLakeList = []
    
    nLines = len(shpNow)
    
    try:
        for i in range(0,nLines):
            
            lineNow = shpNow[i] # select current lake outline
            lineDate = lineNow['properties']['Date']
            imageName = lineNow['properties']['ImagePath']
            
            lineGeom = geometry.asShape(lineNow['geometry'])
            polyGeom = makeSinglePolygonFromSingleLine(lineGeom)
            
            polyArea = getPolyAream2(polyGeom) # polygon area [m^2]
            
            # get centroid coordinates
            cent = lineGeom.centroid
            cx = cent.x
            cy = cent.y
            centTuple = (cx,cy)
            
            ## SAMPLE CLIMATE DATA         
            akAlbers_epsg = 3338 # climate data is in alaska albers
            wgs_epsg = 4326 # lake centroids are in wgs84
            
            dst_crs = Proj(init='EPSG:'+str(akAlbers_epsg))
            src_crs = Proj(init='EPSG:'+str(wgs_epsg))
            
            akAlbersCoords = transform(src_crs,dst_crs,cx,cy) # get albers coordinates
            
            # get winter precip in 1960s, 1980s, 2000s, and change 2000s-1960s
            pptOldNum = sampleRasterAtXY(pptOldFn,akAlbersCoords[0],akAlbersCoords[1])
            pptMidNum = sampleRasterAtXY(pptMidFn,akAlbersCoords[0],akAlbersCoords[1])
            pptNewNum = sampleRasterAtXY(pptNewFn,akAlbersCoords[0],akAlbersCoords[1])
            pptChangeNum = sampleRasterAtXY(pptChangeFn,akAlbersCoords[0],akAlbersCoords[1])
       
            # get summer air temp in 1960s, 1980s, 2000s, and change 2000s-1960s
            tasOldNum = sampleRasterAtXY(tasOldFn,akAlbersCoords[0],akAlbersCoords[1])
            tasMidNum = sampleRasterAtXY(tasMidFn,akAlbersCoords[0],akAlbersCoords[1])
            tasNewNum = sampleRasterAtXY(tasNewFn,akAlbersCoords[0],akAlbersCoords[1])
            tasChangeNum = sampleRasterAtXY(tasChangeFn,akAlbersCoords[0],akAlbersCoords[1])
            
            # get elevation data
            elevationNum = sampleRasterAtXY(demFn,cx,cy)    
            
            # MAKE DICTIONARY OUTPUT
            dictEntryName = lakeCode + '-' + lineDate
            
            # make output dictionary entry
            dictEntry = {'area_m2':polyArea,'cent_x':cx,'cent_y':cy \
                         ,'date':lineDate,'imageName':imageName \
                         ,'polygonGeometry':polyGeom, 'tempJJA_1960':tasOldNum[0] \
                         ,'tempJJA_1980':tasMidNum[0], 'tempJJA_2000':tasNewNum[0] \
                         ,'tempJJA_2000-1960_degC':tasChangeNum[0], 'prcpDJF_1960':pptOldNum[0] \
                         ,'prcpJJA_1980':pptMidNum[0], 'prcpJJA_2000':pptNewNum[0] \
                         ,'prcpDJF_2000-1960_mm':pptChangeNum[0],'elevation_m':elevationNum[0]
                         }
            
            # store it
            thisLakeDict[dictEntryName] = dictEntry
        
            # MAKE LAKE OUTPUT
            listEntry = [lakeCode, dictEntryName, lineDate, cx, cy, polyArea \
                        ,elevationNum[0],tasOldNum[0],tasMidNum[0],tasNewNum[0] \
                        ,tasChangeNum[0],pptOldNum[0],pptMidNum[0],pptNewNum[0] \
                        , pptChangeNum[0]
                        ]
            
            arrayEntry = np.array(listEntry)
            
            if i == 0:
                allOutlines = arrayEntry
                lakeCentList = [centTuple]
            else:
                allOutlines = np.vstack((allOutlines,arrayEntry))
                lakeCentList.append(centTuple)
        
        # add outline data to master list    
        if fn == shpList[0]:
            allLakeArray = allOutlines
            allLakeCentList = lakeCentList
        else:
            allLakeArray = np.vstack((allLakeArray,allOutlines))
            allLakeCentList.append(lakeCentList)
            
        # add outline to master dictionary    
        allLakeDict[lakeCode] = thisLakeDict
    except:
        print "Skipping lake above"    
    
    shpIterCount+=1
    
    shpNow.close()


## SAVE OUT

# Define a feature geometry with attributes
schema = {
'geometry': 'Point',
'properties': {'id': 'int',
               'imageDate':'str',
               'cent_lon':'float',
               'cent_lat':'float',
               'lakeCode':'str',
               'area_m2':'float'
               }
}

# Define another feature geometry with attributes
schemaPoly = {
'geometry': 'Polygon',
'properties': {'id': 'int',
               'imageDate':'str',
               'cent_lon':'float',
               'cent_lat':'float',
               'lakeCode':'str',
               'area_m2':'float'
               }
}

# Write a new polygon shapefile
with fiona.open(fullFnPolyOut, 'w', 'ESRI Shapefile', schemaPoly,crs=from_epsg(4326)) as c:
    
    idCount = 1 # iteration counter for naming shape id
    
    keyIter = allLakeDict.iterkeys()
    
    for n in range(0,len(allLakeDict)):
        
        keyNow = keyIter.next()
        
        lakeNow = allLakeDict[keyNow]
        
        otherKeyIter = lakeNow.iterkeys()
        
        for m in range(0,len(lakeNow)):
            
            otherKeyNow = otherKeyIter.next()
            
            dataNow = lakeNow[otherKeyNow]
            
            c.write({
                    'geometry': geometry.mapping(dataNow['polygonGeometry']),
                    'properties': {'id':idCount,'lakeCode':keyNow,'imageDate':dataNow['date'],'cent_lon':dataNow['cent_x'],'cent_lat':dataNow['cent_y'],'area_m2':dataNow['area_m2']},
                    })
    
            idCount+=1 # update feature count
    

# make shapely geometry of all lake centroids
fixedLakeCentList = fixTupleList(allLakeCentList) # make the list look like how shapely expects it
mPoint = geometry.MultiPoint(fixedLakeCentList) # make a multipoint geometry 

# Write a new Shapefile
with fiona.open(fullFnCentOut, 'w', 'ESRI Shapefile', schema,crs=from_epsg(4326)) as c:

    for k in range(0,len(mPoint)):
        
        linkedArray = allLakeArray[k]
        
        c.write({
                'geometry': geometry.mapping(mPoint[k]),
                'properties': {'id':k,'lakeCode':linkedArray[0],'imageDate':linkedArray[2],'cent_lon':linkedArray[3],'cent_lat':linkedArray[4],'area_m2':linkedArray[5]},
                })


# overwrite projection file
# very hacky; copies source projection file because fiona having issues outputting one
srcPrjFile = fn[:-3] + 'prj' # pull last shapefile for projection info
dstPrjFile = fullFnCentOut[:-3] + 'prj'
dstPrjFile2 = fullFnPolyOut[:-3] + 'prj'
shutil.copyfile(srcPrjFile, dstPrjFile)
shutil.copyfile(srcPrjFile, dstPrjFile2)

# save master file out as a csv
np.savetxt(fullFnTableOut,allLakeArray,delimiter=',',fmt='%s')
# output = [lakeCode,lakeCodeWithDate,date,longitude_degrees,latitude_degrees,area_m2,elevation_m,tempJJA_1960s_degC,tempJJA_1980s_degC,tempJJA_2000s_degC,tempChangeJJA_2000-1960_degC,prcpDJF_1960s_mm,prcpDJF_1980s_mm,prcpDJF_2000s_mm,prcpDJF_2000-1960_mm]
# prcp = precipitation, temp = air temperature
# negative precip change indicates drying in recent times. negative temp chang eindicates cooling in recent time.
# DJF = december, january february (winter)
# JJA = june, july, august (summer)

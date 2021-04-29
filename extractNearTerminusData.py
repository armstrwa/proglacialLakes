# -*- coding: utf-8 -*-
"""
Script to extract near-terminus data for Hannah's study lakes

Created on Thu Nov 07 12:37:32 2019

@author: armstrongwh
"""

import os
# import sys
import numpy as np
# import scipy as sci 
# import csv
import glob
from functools import partial
import matplotlib.pyplot as plt
import json
from shapely import geometry
from shapely.ops import transform
import fiona
from osgeo import ogr, osr, gdal
from pyproj import Proj, transform
from fiona.crs import from_epsg
from gdalconst import *
import datetime as dt
import time
import matplotlib.gridspec as gridspec
import utm

# these needed for querying intersection with rgi (in fast way)
import geopandas as gpd
from geopandas.tools import sjoin

# needed for sampling raster by mask, geotiff reprojection, merging rasters
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.crs import CRS
from rasterio.merge import merge

# below needed for finding epsg.csv
os.environ['GDAL_DATA'] = u'C:\\Users\\armstrongwh\\AppData\\Local\\Continuum\\anaconda2\\Library\\share\\gdal'


'''
Based on the website below, I got my GDAL_DATA path right, which was causing an error CRSError: Invalid CRS on the line     transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds)
https://stackoverflow.com/questions/45883445/how-to-fix-the-enviroment-variable-gdal-data-path-set
Used to say 
os.environ['GDAL_DATA'] and outputs 'C:\\cygwin64\\bin'

then ran
(website said line below, did not work)
#os.environ['GDAL_DATA'] = os.environ['CONDA_PREFIX'] + r'\Library\share\gdal'
(line below DOES WORK)
os.environ['GDAL_DATA'] = u'C:\\Users\\armstrongwh\\AppData\\Local\\Continuum\\anaconda2\\Library\\share\\gdal'

proj lib was that shown below, and did not change this
'C:\\Users\\armstrongwh\\AppData\\Local\\Continuum\\anaconda2\\lib\\site-packages\\rasterio\\proj_data'
'''

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

# Function to get coordinates from shapefile from study glaciers shapefile
# Outputs array of x, y coordinates of shapefile vertices
def getShapefileCoordinatesFromMultiline(shapefileIn,transName=None):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	ds = driver.Open(shapefileIn,0)
	lyr = ds.GetLayer(0)
	
	numFeats = lyr.GetFeatureCount()
	
	for i in range(0,numFeats):
		featNow = lyr.GetFeature(i)
		transNow = featNow.GetField(1)
		
		if transName is None or transName == transNow:
			
			geom = featNow.geometry()
			numPoints = geom.GetPointCount()
			
			x = []
			y = []
			
			for j in range(0,numPoints):
				x.append( geom.GetX(j) )
				y.append( geom.GetY(j) )
		
	coordArray = np.array((x,y))
	coordArrayTranspose = coordArray.transpose()
		
	return coordArrayTranspose


# Function to get coordinates from shapefile from evenly space shapefile
# Outputs array of x, y coordinates of shapefile vertices
def getShapefileCoordinatesFromMultipoint(shapefileIn):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	ds = driver.Open(shapefileIn,0)
	lyr = ds.GetLayer(0)
	
	numPoints = lyr.GetFeatureCount()
			
	x = []
	y = []
	
	for i in range(0,numPoints):
		pointNow = lyr.GetFeature(i)
		geom = pointNow.geometry()
		x.append( geom.GetX() )
		y.append( geom.GetY() )
		
	coordArray = np.array((x,y))
	coordArrayTranspose = coordArray.transpose()
		
	return coordArrayTranspose

# Function to sample raster at specified coordinates
# Returns values of raster at x,y
def sampleRasterAtXY(rasterIn,x,y):
	
	z = []
	
	imageIj = xyToImageIj(x,y)
	
	rast = gdal.Open(rasterIn)
	band = rast.GetRasterBand(1) # assumes single band raster
	bandArr = band.ReadAsArray() # convert band into array
	
	samplePtNum = imageIj.shape[0] # number of sample points
	
	for i in range(0,samplePtNum):
		z.append( bandArr[ np.round(imageIj[i,1]), np.round(imageIj[i,0] ) ] )
	
	return z
	
# Calculate the distance along a line, given vertex coordinates
def distanceAlongLine(x,y):
	
	dist = [0] # starts at 0 distance
	
	numPoints = len(x)
	
	for i in range(1,numPoints):
		oldEasting = x[i-1]
		oldNorthing = y[i-1]
		nowEasting = x[i]
		nowNorthing = y[i]
		dE = nowEasting - oldEasting
		dN = nowNorthing - oldNorthing
		distChange = np.sqrt(dE**2 + dN**2)
		dist.append( distChange + dist[i-1] )
		
	return dist

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


def convert4326toUtm(polyGeom):
    '''
    This assumes the polygon has lat/lon coordinates (epsg 4326, wgs84)
    input = shapely polygon geometry with lat/lon coordinates
    output = same polygon, but in utm coordinates
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
   
    return utmPoly


def sampleRasterByShpGeom(shapelyGeometry,rasterFn):
    # extract the raster values values within the polygon 
    with rasterio.open(rasterFn) as src:
        src = rasterio.open(rasterFn)
        out_image, out_transform = mask(src,[geometry.mapping(shapelyGeometry)], crop=True)
    
    return out_image


def reprojectGeotiff(fnIn,fnOut,dst_epsg_num):
    '''
    Write a new reprojected geotiff to file
    input: fnIn = input full file path (file to reproject)
        fnOut = output full file path
        dst_epsg_num = epsg number (integer) for output coordinate system
    output: reprojected geotiff located at fnOut
    
    slightly modified from https://rasterio.readthedocs.io/en/latest/topics/reproject.html#reprojecting-a-geotiff-dataset
    '''
    
    dst_crs = 'EPSG:' + str(dst_epsg_num)
    
    with rasterio.open(pathNow) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
    
        with rasterio.open(outFn, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)


# near-terminus shapefile path
sampleShpPath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\Glacier Ablation Area Shape File and Valley Section Line-20191111T174812Z-001\\'
nearTermPath = sampleShpPath + 'AblationZoneSampleArea.shp'
nearTermPath4326 = sampleShpPath + 'AblationZoneSampleArea_epsg4326.shp'

# lake centroid shapefile path
lakeCentShpFn = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\allLakes_TYrYrNameLaLon_23sep2019\/allLakeOneCentroidPerLake.shp'

# joined region 01-02 RGI path
rgiPath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\gis_datasets\\00_rgi60\\0102_rgi60_alaska_westernCandaa_merged.shp'
rgiPathHussSubset = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\data\\huss\\huss_rgi_coverage.shp'


# ice thickness data folder
iceThicknessPath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\gis_datasets\\farinotti_ice_thickness\\'
iceThicknessSubsetPath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\data\\iceThicknessSubset\\'


# alaska-canada merged dem
akCanMergeDem = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\gis_datasets\\ALASKA_GEOSPATIAL\\mergedAkCanAndPnw_dem_epsg4326.tif'

# path to shapefile showing extent of hydrosheds dem data in study area
hydroshedExtentShpFn = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\gis_datasets\\hydrosheds_pnw_extent_epsg4326.shp'
# now AK-Can merged DEM extent
akCanMergeExtentShpFn = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\gis_datasets\\alaska_canada_demMerge_googDataExtent_epsg4326.shp'

# open near terminus sampling area shapefile
nearTermShp = fiona.open(nearTermPath4326)

#shapefile = gpd.read_file(nearTermPath4326)


rgiShp = fiona.open(rgiPath)
rgiHussShp = fiona.open(rgiPathHussSubset)

# find overlapping features
samplePoly = gpd.GeoDataFrame.from_file(nearTermPath4326) 
rgiPoly = gpd.GeoDataFrame.from_file(rgiPathHussSubset) 
pointInPolys = sjoin(samplePoly, rgiPoly, how='left')
# output is a dictionary (geopandas geodataframe), pointInPolys.keys() to see what fields are called

# testing use geopandas
rgiIsectList = pointInPolys['rgiId'] # this is list of RGI IDs that overlap with sample shapefile
rgiIndexList = pointInPolys['id'] # this is index, use like rgiHussShp[int(rgiIndexList[3])]
test = pointInPolys['geometry'][4] # this is a near-terminus sample geometry

listLen = len(rgiIsectList) # hard coded range in loop below to go to 100 b/c issue below

'''

Start getting index error at j = 101. Dataframe size is 108, but rows beyond 101 don't exist?
Maybe I'm just not understanding how pandas works?
KeyError: u'the label [101] is not in the [index]'

'''

noDataValThickness = 0 # no data value in farinotti ice thickness data set, mask out for later analysis


if 0: # run this to reproject ice thickness rasters into WGS84
    # loop through overlapping rgi outlines, reproject ice thickness data if it exists
    for j in range(0,101): # see above (listLen) for why this hard coded
        rgiIdNow = str(pointInPolys.loc[j,'rgiId']) # get current rgi id
        pathNow = iceThicknessPath + rgiIdNow + '_thickness.tif'
        outPath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\data\\iceThicknessSubset\\'
        outFn = outPath  + rgiIdNow + '_thickness_epsg4326.tif'
       
    
        if rgiIdNow == 'nan':
            print "Skipping, no ice thickness data on iteration: " + str(j)
        else:
            if type(pointInPolys.loc[j,'rgiId']) == unicode:
                reprojectGeotiff(pathNow,outFn,4326)
                print "Reprojected geotiff: " + rgiIdNow + '_thickness_epsg4326.tif'     
            else:
                rgiMatchList = pointInPolys.loc[j,'rgiId']                
                numMatches = len(rgiMatchList)
                print str(numMatches) + " RGI matches on iteration " + str(j)
                
                for k in range(0,numMatches):
                    matchNow = str(rgiMatchList.iloc[k])
                    
                    pathNow = iceThicknessPath + matchNow + '_thickness.tif'
                    outFn = outPath  + matchNow + '_thickness_epsg4326.tif'
                   
                    reprojectGeotiff(pathNow,outFn,4326)
                    print "Reprojected geotiff: " + matchNow + '_thickness_epsg4326.tif'     
                        


# loop over overlapping ice thickness rasters and sample
outArr = None
for i in range(0,101): # iterate over overlapping rgi and near-terminus sample polygons
    print i
    rgiNow = rgiIsectList[i] # get current RGI ID
    
    if type(rgiNow) == float: # this true when nan
        print "Skipping, no Farinotti ice thickness data on iteration: " + str(i)
    else:        

        if type(pointInPolys.loc[i,'rgiId']) == unicode: # if there is only 1 overlapping ice thickness dataset
            fnNow4326 = iceThicknessSubsetPath + str(rgiNow) + '_thickness_epsg4326.tif' # relevant ice thickness filename
            # get geometry (which will be WGS84 for near-terminus shapefile) and convert to UTM (to match Farinotti dataset)
            sampleGeomNow_wgs84 = pointInPolys['geometry'][i]
            sampleGeomNow_utm = convert4326toUtm(sampleGeomNow_wgs84)            
            tmp = sampleRasterByShpGeom(sampleGeomNow_wgs84,fnNow4326) # subset ice thickness dataset by near-terminus polygon
            thicknessWithNoData = tmp.data[0] # just pull data, has lots of other stuff in tmp (which is a numpy masked array)
        else:
            rgiMatchList = pointInPolys.loc[i,'rgiId']    
            sampleGeomNow_wgs84 = pointInPolys['geometry'][i]    
            
            numMatches = len(rgiMatchList)
            print str(numMatches) + " RGI matches on iteration " + str(i)
            
            mergeList = [] # filenames to merge
            
            for k in range(0,numMatches):
                matchNow = str(rgiMatchList.iloc[k])
                
                pathNow = iceThicknessSubsetPath + matchNow + '_thickness_epsg4326.tif'
                
                rastNow = rasterio.open(pathNow)
                
                mergeList.append(rastNow)
            
            mosaic, out_trans = merge(mergeList)
            
            #newProfile = rastNow.profile # copy profile from ice thickness map
            newMeta = rastNow.meta
            #newProfile['affine'] = out_trans # change transformation to reflect new image extent
            
            out_meta = rastNow.meta.copy()
            
            # reformat out_trans to align with newMeta['transform']
            transFmt = (out_trans[2], out_trans[0], out_trans[1],  out_trans[5], out_trans[3], out_trans[4])
            
            out_meta.update({"driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": transFmt
            #"crs": "+proj=utm +zone=35 +ellps=GRS80 +units=m +no_defs "
            }
            )
            
            # Write your the ndvi raster object
            with rasterio.open(iceThicknessSubsetPath + 'tempRast.tif', 'w', **out_meta) as dst:
                dst.write(mosaic)
            
            # mosaic explanation here https://automating-gis-processes.github.io/CSC18/lessons/L6/raster-mosaic.html
            tmp = sampleRasterByShpGeom(sampleGeomNow_wgs84.iloc[0],iceThicknessSubsetPath + 'tempRast.tif') # subset ice thickness dataset by near-terminus polygon
            thicknessWithNoData = tmp.data[0] # just pull data, has lots of other stuff in tmp (which is a numpy masked array)

        # the line below removes no data from above to allow computation of stastics, but loses spatial reference
        thicknessNoNoData = np.extract(thicknessWithNoData != noDataValThickness, thicknessWithNoData)
        
        sampleAreaM2 = sampleGeomNow_utm.area # near terminus sample area is m^2
        sampleAreaPixels = len(thicknessNoNoData) # how many data points used for computing near-terminus thickness statistics
        
        # computing ice thickness statistics for output
        p25,p50,p75 = np.percentile(thicknessNoNoData,[25,50,75]) # 25, 50 (median), 75th percentile near terminus ice thickness
        iqr = p75-p25 # interquartile range of near-terminus ice thickness
        outMean = np.nanmean(thicknessNoNoData)
        outMin = np.nanmin(thicknessNoNoData)
        outMax = np.nanmax(thicknessNoNoData)
        outStd = np.nanstd(thicknessNoNoData)
    
        if type(rgiNow) == unicode:
            nameOut1 = str(rgiNow)
            nameOut2 = str(-999)
        else:
            nameOut1 = str(rgiNow.iloc[0])
            nameOut2 = str(rgiNow.iloc[1])
            
        tmpArr = [i,nameOut1,nameOut2,sampleAreaM2,sampleAreaPixels,outMin,p25,p50,outMean,p75,outMax,iqr,outStd]
        
        if outArr is None:
            outArr = tmpArr
        else:
            outArr = np.vstack((outArr,tmpArr))

    

np.savetxt(iceThicknessSubsetPath+'nearTerminus_iceThicknessStates.csv',outArr,delimiter=',',fmt="%s")


if 0:
    '''
    stuff below may be needed for auto dem downloading, having gotten it working as of 17 jan 2020, 12:30
    '''
    from bounds import RasterBounds, GeoBounds
    def makePolyGeomFromRasterBounds(rastFn):
        rast = rasterio.open(rastFn)
        bb = rast.bounds
        l = bb[0]
        b = bb[1]
        r = bb[2]
        t = bb[3]
        
        ul = (l,t)
        ll = (l,b)
        ur = (r,t)
        lr = (r,l)
        
        mpGeo = geometry.asPolygon((ul,ur,lr,ll))
        
    bb = RasterBounds(affine_transform=profile['affine'],
                      profile=profile, latlon=True)
    
    
    




# open lake centriod shapefile
lakeCentShp = fiona.open(lakeCentShpFn)

nLake = len(lakeCentShp)


# get array of lake centroid coordinates
centArr = getShapefileCoordinatesFromMultipoint(lakeCentShpFn)


for i in range(0,nLake):
    lakeNow = lakeCentShp[i]
    lakeGeom = geometry.asShape(lakeNow['geometry'])
    
    

# not testing for intsersection after merging DEM
if 0:
    if lakeGeom.intersects(hShdGeom): # if data within hydrosheds region
        # print(i,lakeGeom.intersects(hShdGeom))
        
        # extract elevation data from hydrosheds
    #elif lakeGeom.intersects(akCanMergeGeom):
        # extract elevation data from merged AK/CAN DEM
    else:
        print('i = ',str(i),'lake code = ',str(lakeNow['properties']['lakeCode']),'no DEM intersection')



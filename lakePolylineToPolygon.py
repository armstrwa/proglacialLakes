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


## USER INPUTS
    
# identify shapefile to process
path = u"C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\sampleshapefile\\"
fn = 'T1987Mala59140.shp'

fnOut = path + 'T1987Mala59140_polygons.shp' # output filename


## PROCESSING

fullFn = path + fn

# open shapefile
with fiona.open(fullFn,'r') as c:        
        lineNum = len(c) # how many outlines are in collection?
        
        dataDicts = [] # empty container
        
        # iterate over transects
        for i in range(0,lineNum):
            dataDicts.append(readLineCoordsFromMultiLineShp(fullFn,i))
            
        # make multipolygon once all dictionaries have been collected
        makePolygonShpFromDictList(dataDicts,fnOut)
        
# overwrite projection file
# very hacky; copies source projection file because fiona having issues outputting one
srcPrjFile = fullFn[:-3] + 'prj'
dstPrjFile = fnOut[:-3] + 'prj'
shutil.copyfile(srcPrjFile, dstPrjFile)


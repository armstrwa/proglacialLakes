# -*- coding: utf-8 -*-
"""
Script to pull Farinotti ice thickness data for lake area change analysis

Created on Mon Jul 22 09:51:59 2019

@author: armstrongwh
"""

## import modules - the order you do this in is very important due to bugs in shapely and fiona. Make sure you import in this order.
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

## user defined parameters

lakeChangeFn = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\GlacierDataWithLakeID_23sep2019.csv'
iceThicknessPath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\gis_datasets\\farinotti_ice_thickness\\'
outFolder = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\data\\iceThicknessSubset\\'

## ANALYSIS

outArr = [] # storage for writing out

lakeData = np.genfromtxt(lakeChangeFn,delimiter=',',skip_header=1,dtype=None)
lakeNum = len(lakeData)

for i in range(0,lakeNum): # loop over lakes
    lakeNow = lakeData[i] # pull current lake
    rgiNow = lakeNow[1] # RGI ID number for current lake
    
    rastFn = glob.glob(iceThicknessPath + rgiNow + '*.tif') # find ice thickness tiff for upstream glacier
    
    ds = gdal.Open(rastFn[0]) # open file
    arr = ds.ReadAsArray() # turn into numpy raster
    
    gt = ds.GetGeoTransform() # vector with raster origin (x,y at upper left), pixel size, and rotation param
    pxArea = gt[1]**2 # pixel area (m^2 for farinotti dataset, 50x50 m)
    
    # make 0s (off glacier) into nan (so doesn't change stats)
    noDataInd = arr == 0   
    arr[noDataInd] = np.nan
    
    # calculate statistics
    meanH = np.nanmean(arr) # mean ice thickness
    maxH = np.nanmax(arr) # max ice thickness
    stdH = np.nanstd(arr) # standard deviation of ice thickness
    p25,p50,p75 = np.nanpercentile(arr,[25,50,75]) # median and interquartile range of ice thickness
    vol_m3 = pxArea*np.nansum(arr) # glacier volume [m^3]
    
    # copy raster to subset folder
    fnOnly = str(rastFn[0].split('\\')[-1]) # just get raster filename (no path)
    shutil.copyfile(rastFn[0],outFolder+fnOnly) # copy file

    outLine = [rgiNow, p25,p50,p75,meanH,maxH,stdH,vol_m3]
    
    if i == 0:
        outArr = outLine
    else:
        outArr = np.vstack((outArr,outLine))

np.savetxt(outFolder+'iceThicknessExtractedData.csv',outArr,fmt="%s",delimiter=',')


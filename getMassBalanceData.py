# -*- coding: utf-8 -*-
"""
Script to subset and compile Huss & Hock mass balance data to our study glaciers


Created on Thu Jul 25 10:25:55 2019

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
import time

## user defined parameters

lakeChangeFn = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\gis\\GlacierDataWithLakeID_23sep2019.csv'
massBalancePath = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\data\\huss\\'
outFolder = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\projects\\lakes\\data\\massBalanceSubset\\'

akTauFn = massBalancePath + 'dBdz_alaska.csv'
canTauFn = massBalancePath + 'dBdz_westerncanada.csv'
akNetBalFn = massBalancePath + 'alaska_annual_balance_sfc_r1.csv'
canNetBalFn = massBalancePath + 'westerncanada_Annual_Balance_sfc_r1.csv'

# shapefile for merged RGI regions 01 and 02
rgiMergedShpFn = u'C:\\Users\\armstrongwh\\Google Drive\\appstate\\gis_datasets\\00_rgi60\\0102_rgi60_alaska_westernCandaa_merged.shp'
mergedDbDzFn = massBalancePath + 'merged_dBdz_alaskaAndWCanada.csv'
mPolyShpFnOut = massBalancePath + 'huss_rgi_coverage.shp'

## ANALYSIS

outArr = [] # storage for writing out

lakeData = np.genfromtxt(lakeChangeFn,delimiter=',',skip_header=1,dtype=None)
lakeNum = len(lakeData)

akTauData = np.genfromtxt(akTauFn,delimiter=',',skip_header=1,dtype=None)
canTauData = np.genfromtxt(canTauFn,delimiter=',',skip_header=1,dtype=None)
akNetBalData = np.genfromtxt(akNetBalFn,delimiter=',',skip_header=1,dtype=None)
canNetBalData = np.genfromtxt(canNetBalFn,delimiter=',',skip_header=1,dtype=None)
mergedTauData = np.genfromtxt(mergedDbDzFn,delimiter=',',skip_header=1,dtype=None)

## Get data from Matthias' dbdz files
for i in range(0,lakeNum): # loop over lakes
    lakeNow = lakeData[i] # pull current lake
    rgiNow = lakeNow[1] # RGI ID number for current lake
    
    rgiZone = int(rgiNow[6:8]) # identify RGI zone from id
    rgiNum = int(rgiNow[-5:]) # just RGI number
    
    # point to proper dataset whether ak or western canada
    if rgiZone == 1:
        tauData = akTauData
    elif rgiZone == 2:
        tauData = canTauData
    
    # loop over all glaciers in relevant RGI zone
    for j in range(0,len(tauData)):
        glacierId = tauData[j][0] # get current glacier's RGI id
        glacierNum = int(glacierId[-5:]) # just get rgi id number
        if rgiNum == glacierNum: # is this the glacier we're interesting in?
            switch = 1 # switch to 1 to show we've found a match for this lake
            #print "Found match at i, j: "
            #print i, j
            
            tauDataNow = list(tauData[j]) # highlight matching data, make it a list for proper indexing
            outLine = [rgiNow] + tauDataNow[3:] # first part appends RGI code from lake
            
            if i == 0:
                outArr = outLine
            else:
                outArr = np.vstack((outArr,outLine))
        
    if j == len(tauData)-1 and switch == 0:
        print "No match found on i = " + str(i)
        print "RGI ID = " + rgiNow
        print "Lake ID = " + lakeNow[23]
        
        # make nan's to fill in array
        nanArr =  np.nan*np.ones((12))
        nanList = nanArr.tolist()
        outLine = [rgiNow] + nanList
        
        if i == 0:
            outArr = outLine
        else:
            outArr = np.vstack((outArr,outLine))
            
    switch = 0 # reset the switch that says we've found a match
        
# save out
np.savetxt(outFolder+'dbdzExtractedData.csv',outArr,fmt="%s",delimiter=',')


## Get data from Matthias' annual balance time series data
#  Need to do this twice b/c files are not ordered the same so j in one ~= j in other
outArr = [] # storage for writing out
for i in range(0,lakeNum): # loop over lakes
    lakeNow = lakeData[i] # pull current lake
    rgiNow = lakeNow[1] # RGI ID number for current lake
    
    rgiZone = int(rgiNow[6:8]) # identify RGI zone from id
    rgiNum = int(rgiNow[-5:]) # just RGI number
    
    # point to proper dataset whether ak or western canada
    if rgiZone == 1:
        netData = akNetBalData
    elif rgiZone == 2:
        netData = canNetBalData
    
    # loop over all glaciers in relevant RGI zone
    for j in range(0,len(netData)):
        glacierNum = netData[j][0] # get current glacier's RGI number
        if rgiNum == glacierNum: # is this the glacier we're interesting in?
            switch = 1 # switch to 1 to show we've found a match for this lake
            #print "Found match at i, j: "
            #print i, j
            
            netDataNow = list(netData[j]) # highlight matching data, make it a list for proper indexing
            meanBal80s = np.mean(netDataNow[1:11]) # average net annual balance in the 1980s
            meanBal90s = np.mean(netDataNow[11:21]) # average net annual balance in the 1990s
            meanBal00s = np.mean(netDataNow[21:31]) # average net annual balance in the 2000s
            meanBal10s = np.mean(netDataNow[31:]) # average net annual balance in the 2010s
            meanBalAll = np.mean(netDataNow[1:]) # mean balance over whole period
            sumBalAll = np.sum(netDataNow[1:],axis=0)
            
            
            outLine = [rgiNow, meanBal80s, meanBal90s, meanBal00s, meanBal10s, meanBalAll, sumBalAll] # first part appends RGI code from lake
            
            if i == 0:
                outArr = outLine
            else:
                outArr = np.vstack((outArr,outLine))
        
    if j == len(tauData)-1 and switch == 0:
        print "No match found on i = " + str(i)
        print "RGI ID = " + rgiNow
        print "Lake ID = " + lakeNow[23]
        
        # make nan's to fill in array
        nanArr =  np.nan*np.ones((5))
        nanList = nanArr.tolist()
        outLine = [rgiNow] + nanList
        
        if i == 0:
            outArr = outLine
        else:
            outArr = np.vstack((outArr,outLine))
            
    switch = 0 # reset the switch that says we've found a match
          
# save out
np.savetxt(outFolder+'annualBalanceExtractedData.csv',outArr,fmt="%s",delimiter=',')


### testing for RGI glaciers in Huss dataset
#   this is very slow, only run if needed
shp = fiona.open(rgiMergedShpFn) # open merged RGI shapefile
nShp = len(shp) # how many outlines do we have

idList = []

for m in range(0,nShp): # loop over merged RGI to make a list of all RGI ids for later indexing
    f = shp[m]
    rgiNow = str(f['properties']['RGIId'])
    idList.append(rgiNow[-8:])

idArr = np.array(idList) # needs to be an array for indexing to work right

polyList = [] # create empty geometry list
rgiList = []

tic = time.clock()
for k in range(0,len(mergedTauData)):
    glacierId = mergedTauData[k][0] # get current glacier's RGI id
    #print "Processing glaciers % : " + str(float(k)/float(len(mergedTauData)))
    
    try: # get errors if no match found
        ind = int(np.where(idArr==glacierId[-8:])[0][0]) # where does this rgi id appear in the list, use this as index
        f = shp[ind]
      
        geom = f['geometry']
        polyList.append(geometry.shape(geom))
        rgiList.append(str(f['properties']['RGIId']))
    except:
        print "No match for: " + str(glacierId)
#    for m in range(0,nShp): # loop over merged RGI
#        f = shp[m]
#        rgiNow = str(f['properties']['RGIId'])
#        
#        if rgiNow[-8:] == glacierId[-8:]: # is this the glacier we're interested in? need to do last 8 b/c different period placement in 6.0 between RGI and Huss

toc = time.clock()
toc - tic

mPoly = geometry.MultiPolygon(polyList)

# Define a polygon feature geometry with one attribute
schema = {
'geometry': 'Polygon',
'properties': {'id': 'int',
               'rgiId':'str',
               }
}


# Write a new Shapefile
with fiona.open(mPolyShpFnOut, 'w', 'ESRI Shapefile', schema,crs=from_epsg(4326)) as c:
## If there are multiple geometries, put the "for" loop here
    for k in range(0,len(mPoly)):
        idNow = rgiList[k]
        
        c.write({
                'geometry': geometry.mapping(mPoly[k]),
                'properties': {'id':k,'rgiId':idNow},
                })
		
print "Sucessfully wrote shapefile: " + mPolyShpFnOut

# overwrite projection file
# very hacky; copies source projection file because fiona having issues outputting one
srcPrjFile = rgiMergedShpFn[:-3] + 'prj' # pull last shapefile for projection info
dstPrjFile = mPolyShpFnOut[:-3] + 'prj'
shutil.copyfile(srcPrjFile, dstPrjFile)
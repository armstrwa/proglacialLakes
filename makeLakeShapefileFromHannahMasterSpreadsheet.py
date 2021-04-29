'''
Make shapefile with lots of lake attributes from Hannah's master spreadsheet
10 Apr 2020, wha

NOTE: BUILT FOR FORMAT OF 'MASTER_LAKE_DATA_20200205...reformat.csv' - if column numbers, number of header rows, nan representation, etc. changes, this may break
'''

import numpy as np
import geopandas as gp
from shapely import geometry
from shapely.ops import transform
import fiona
from fiona.crs import from_epsg

# point ot master spreadsheet
fn = '/Users/williamarmstrong/Google Drive/appstate/projects/lakes/data/Master_Lake_Data_20200205_reformat.csv'
fullFnOut = '/Users/williamarmstrong/Google Drive/appstate/projects/lakes/gis/HannahMasterLake_20200205_tableData_epsg4326.shp'

# read in data
data = np.genfromtxt(fn,delimiter=',',dtype=None)

# initialize storage matrices
if 1:
	lakeCodes = []
	lons = []
	lats = []
	elevs = []
	jjaTemps = []
	jjaDTemps = []
	djfPrcps = []
	djfDPrcps = []
	rgiLinks = []
	glAreas = []
	glSlopes = []
	glZRanges = []
	glSurges = []
	glWidths = []
	cumMassBals = []
	respTimes = []
	glNames = []
	dAreas = []
	lakeTypes = []
	glTermThicks = []
	growthStyles = []


# Define a feature geometry with attributes
if 1:
	schema = {
	'geometry': 'Point',
	'properties': {'id': 'int',
				   'lakeCode':'str',
				   'rgiLink':'str',
				   'glName':'str',
				   'lakeType':'str',
				   'lon':'float',
				   'lat':'float',
				   'elev':'int',
				   'jjaTempDegC':'float',
				   'jjaTempChagne':'float',
				   'djfPrcpMm':'float',
				   'dfjPrcpChange':'float',
				   'glAreaKm2':'float',
				   'glSlope':'float',
				   'glWidth':'float',
				   'cumMassBal':'float',
				   'glResponseTimeYr':'int',
				   'glMeanThickM':'float',
				   'glElevRange':'float',
				   'glSurge':'int',
				   'areaChangeM2':'float',
				   'growthStyle':'int',
				   'startArea':'float',
				   'logAreaChange':'float'
				   }
	}

# Write a new shapefile, hard-coded to WGS84
with fiona.open(fullFnOut, 'w', 'ESRI Shapefile', schema,crs=from_epsg(4326)) as c:

	for i in range(0,len(data)):
		d = data[i]
	
		#lakeCode = str(d[0]) # lake code
		lakeCode = str(d[1][5:14]) # lake code, above errored on some b/c ascii issue, just stripping from shp name
		lon = d[2] # lake centroid longitude [deg]
		lat = d[3] # lake centroid latitude
		elev = d[4] # lake elevation [m]
		jjaTemp = d[7] # modeled jja temp in 2000s, degrees C
		jjaDTemp = d[8] # modeled jja temp change, 2000s-1960s, deg C
		djfPrcp = d[11] # modeled winter precip, mm
		djfDPrcp = d[12] # modeled winter precip change ,mm
		rgiLink = str(d[55])
		glArea = d[63] # glacier area, km2
		glSlope = d[67] # glacier slope, degrees
		glZRange = d[65] - d[64] # glacier elevation range, m
		glSurge = d[74] # surge type glacier status
		glWidth = d[80] # glacier width, m
		glTermThick = d[113] # near terminus glacier thickness [m]
		cumMassBal = d[88] # cumulative mass balance from Huss & Hock [m we]
		respTime = d[99] # response time from huss [yr]
		glName = str(d[76])
		dArea = d[78] # lake area change, m^2
		logDArea = np.log10(dArea) # log10 area change
		startArea = d[33] # starting lake area [m^2]
		lakeType = str(d[103]) # lake setting classification [see Hannah's notes]
		growthStyle = d[104] # growth style of timeseries (accelerating, decelerating, etc.)

		# make data lists for later plotting
		if 1:
			lakeCodes.append(lakeCode)
			lons.append(lon)
			lats.append(lat)
			elevs.append(elev)
			jjaTemps.append(jjaTemp)
			jjaDTemps.append(jjaDTemp)
			djfPrcps.append(djfPrcp)
			djfDPrcps.append(djfDPrcp)
			rgiLinks.append(rgiLink)
			glAreas.append(glArea)
			glSlopes.append(glSlope)
			glZRanges.append(glZRange)
			glSurges.append(glSurge)
			glWidths.append(glWidth)
			cumMassBals.append(cumMassBal)
			respTimes.append(respTime)
			glNames.append(glName)
			dAreas.append(dArea)
			lakeTypes.append(lakeType)
			glTermThicks.append(glTermThick)
			growthStyles.append(growthStyle)
	
		# make lake centroid point geometry
		pointNow = geometry.Point( float(lon) , float(lat) ) # make it a geometry		

		# write it to the shapefile, and save its hash with it
		c.write({
				'geometry': geometry.mapping(pointNow),
				'properties': {'id': 'int',
							   'lakeCode':lakeCode,
							   'rgiLink':rgiLink,
							   'glName':glName,
							   'lakeType':lakeType,
							   'lon':lon,
							   'lat':lat,
							   'elev':elev,
							   'jjaTempDegC':jjaTemp,
							   'jjaTempChagne':jjaDTemp,
							   'djfPrcpMm':djfPrcp,
							   'dfjPrcpChange':djfDPrcp,
							   'glAreaKm2':glArea,
							   'glSlope':glSlope,
							   'glWidth':glWidth,
							   'cumMassBal':cumMassBal,
							   'glResponseTimeYr':respTime,
							   'glMeanThickM':glTermThick,
							   'glElevRange':glZRange,
							   'glSurge':glSurge,
							   'areaChangeM2':dArea,
							   'growthStyle':growthStyle,
							   'startArea':startArea,
							   'logAreaChange':logDArea							   
							   }
				})



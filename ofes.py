## Get transects from regular grid
# Coordinates and Indexes: Predominantly Meridional
def get_transect_meridional(lato,lono,latt,lont,interpolation='linear'):
	"""
		Parameters: 
		lato,lono = original grid
		latt, lont = transect latitudes and longitudes
		interpolation = interpolation type (default: linear)

		Returns: 
		lat, lon: grid coordinates at the transect location
		ilat, jlon: grid indexes of the transects 

	"""

	from scipy.interpolate import interp1d
	import numpy as np

	lat = lato[(lato >= latt.min()) & (lato <= latt.max())] 
	
	# Longitude axis
	fi = interp1d(latt, lont, kind=interpolation)
	lon = fi(lat)

	# Find longitudes from the ARGO grid
	ilat = []
	jlon = []

	for i in xrange(lat.shape[0]):

		ilat.append(np.abs(lato-lat[i]).argmin())
		jlon.append(np.abs(lono-lon[i]).argmin())


	return lon,lat,jlon,ilat


# Coordinates and Indexes: Predominantly Zonal
def get_transect_zonal(lato,lono,latt,lont,interpolation='linear'):
	"""
		Parameters: 
		lato,lono = original grid
		latt, lont = transect latitudes and longitudes
		interpolation = interpolation type (default: linear)

		Returns: 
		lat, lon: grid coordinates at the transect location
		ilat, jlon: grid indexes of the transects 

	"""

	from scipy.interpolate import interp1d
	import numpy as np


	lon = lono[(lono >= lont.min()) & (lono <= lont.max())] 

	# Longitude axis
	fi = interp1d(lont, latt, kind=interpolation)
	lat = fi(lon)


	# Find longitudes from the ARGO grid
	ilat = []
	jlon = []

	for i in xrange(lon.shape[0]):

		ilat.append(np.abs(lato-lat[i]).argmin())
		jlon.append(np.abs(lono-lon[i]).argmin())

	return lon,lat,jlon,ilat



## Derivatives
def derivative_4d(variable,dc,coordinate_name):
	r"""
		Calculate the central difference derivative 
		between two consecutive points of a variable.

		Parameters: 
		variable: indexes order [time,z,lat,lon]
		dc = resolution (in OFES is constant everywhere)

		Returns: 
		dvdc: variable derivative		
	"""

	import numpy as np
	import gsw


	if coordinate_name == 'time':
		dvdc = (variable[1:,:,:,:] - variable[:-1,:,:,:])/dc

	elif coordinate_name == 'z':
		dvdc = np.ones(variable.shape)*np.nan
		dvdc = dvdc[:,1:,:,:]

		for k in xrange(dc.shape[0]):
			dvdc[:,k,:,:] = (variable[:,k+1:,:,:] - variable[:,k,:,:])/dc[k]	


	elif coordinate_name == 'lon':
		dvdc = np.ones(variable.shape)*np.nan
		dvdc = dvdc[:,:,:,1:]

		for k in xrange(dc.shape[0]):
			dvdc[:,:,k,:] = (variable[:,:,k,1:] - variable[:,:,k,:-1])/dc[k]

	elif coordinate_name == 'lat':
		dvdc = (variable[:,:,1:,:] - variable[:,:,:-1,:])/dc

	return dvdc


def derivative_space3d(variable,dc,coordinate_name):
	r"""
		Calculate the central difference derivative 
		between two consecutive points of a variable.

		Parameters: 
		variable: indexes order [z,lat,lon]
		dc = resolution (in OFES is constant everywhere)

		Returns: 
		dvdc: variable derivative		
	"""

	import numpy as np
	import gsw


	if coordinate_name == 'z':
		dvdc = np.ones(variable.shape)*np.nan
		dvdc = dvdc[1:,:,:]

		for k in xrange(dc.shape[0]):
			dvdc[k,:,:] = (variable[k+1:,:,:] - variable[k,:,:])/dc[k]		

	elif coordinate_name == 'lon':
		dvdc = np.ones(variable.shape)*np.nan
		dvdc = dvdc[:,:,1:]

		for k in xrange(dc.shape[0]):
			dvdc[:,k,:] = (variable[:,k,1:] - variable[:,k,:-1])/dc[k]

	elif coordinate_name == 'lat':
		dvdc = (variable[:,1:,:] - variable[:,:-1,:])/dc

	return dvdc


def derivative_space2d(variable,dc,coordinate_name):
	r"""
		Calculate the central difference derivative 
		between two consecutive points of a variable.

		Parameters: 
		variable: indexes order [lat,lon]
		dc = resolution (in OFES is constant everywhere)

		Returns: 
		dvdc: variable derivative		
	"""

	import numpy as np
	import gsw
		
	if coordinate_name == 'lon':
		dvdc = np.ones(variable.shape)*np.nan
		dvdc = dvdc[:,1:]

		for k in xrange(dc.shape[0]):
			dvdc[k,:] = (variable[k,1:] - variable[k,:-1])/dc[k]

	elif coordinate_name == 'lat':
		dvdc = (variable[1:,:] - variable[:-1,:])/dc

	return dvdc



def derivative_time_space2d(variable,dc,coordinate_name):
	r"""
		Calculate the central difference derivative 
		between two consecutive points of a variable.

		Parameters: 
		variable: indexes order [time,lat,lon]
		dc = resolution (in OFES is constant everywhere)

		Returns: 
		dvdc: variable derivative		
	"""

	import numpy as np
	import gsw

	if coordinate_name == 'time':
		dvdc = (variable[1:,:,:] - variable[:-1,:,:])/dc	

	elif coordinate_name == 'lon':
		dvdc = np.ones(variable.shape)*np.nan
		dvdc = dvdc[:,:,1:]

		for k in xrange(dc.shape[0]):
			dvdc[:,k,:] = (variable[:,k,1:] - variable[:,k,:-1])/dc[k]

	elif coordinate_name == 'lat':
		dvdc = (variable[:,1:,:] - variable[:,:-1,:])/dc

	return dvdc



## Calculate distances from a rectangular lat/lon grid
def dx_from_dlon(lon,lat):
	"""
	Calculate zonal distance at the Earth's surface in m from a 
	longitude and latitude rectangular grid

	Input:
			lon [M,N]: Longitude in degrees
			lat [M,N]: Latitude in degrees

	Output:
			dx   [M,N-1]: Distance in meters

	"""

	## Required packages
	import numpy as np

	# Earth Radius in [m]
	earth_radius = 6371.0e3


	# Convert angles to radians
	lat = np.radians(lat)
	lon = np.radians(lon)
	
	# Zonal distance in radians
	dlon = np.diff(lon,axis=1)
	lat = (lat[:,1:]+lat[:,:-1])/2.0


	# Zonal distance arc 
	a = np.cos(lat)*np.cos(lat)*(np.sin(dlon/2.0))**2.0
	angles = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

	# Distance in meters
	dx = earth_radius * angles

	return dx


def dy_from_dlat(lat):
	"""
	Calculate meridional distance at the Earth's surface in m from a 
	longitude and latitude rectangular grid

	Input:
			lat [M,N]: Latitude in degrees

	Output:
			dy   [M-1,N]: Distance in meters

	"""

	## Required packages
	import numpy as np

	## Meridional resolution (m)
	dy = np.diff(lat,axis=0)
	dy = dy*111194.928

	return dy
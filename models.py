### Auxiliary classes and functions 
class class_builder(object):
	def __init__(self,var,varnames):
		for i in xrange(len(var)):

			exec("self."+varnames[i]+" = var[i]")
			self.varnames = varnames



###Primary functions 
def write_netcdf(filename,varnames,data,title):
	r"""
	Parameters
	----------
	filename: string
		path + filename of netCDF file e. g., '~/Desktop/netcdf_file.nc'

	varnames: list
		names of the variables names e.g., ['var1','var2','var3'] (names should be exactly the same of 
		the data.varnames)
	
	data: class 
		class containing the variables to be saved

	title: string 
		Title or any brief description of the data set 
	"""
	from netCDF4 import Dataset
	import numpy as np 

	# Opening file
	filenc = Dataset(filename,'w',clobber=True)


	# Global Attributes of the file
	filenc.title = title
	filenc.institution = 'The Rosenstiel School of Marine and Atmospheric Science, University of Miami'
	filenc.author = 'MSc Tiago Carrilho Bilo'

	## Creating variables
	for v in xrange(len(varnames)):
		exec("tp = data."+varnames[v]+".dtype")
		exec("d = len(data."+varnames[v]+".shape)")

		if d == 1:
			exec("filenc.createDimension('i"+np.str(v)+"', None)")
			exec("var"+np.str(v)+" = filenc.createVariable(varnames[v], tp.name, ('i"+np.str(v)+"'), zlib=True)")
			exec("var"+np.str(v)+"[:] = data."+varnames[v])			
		elif d == 2:		
			exec("filenc.createDimension('i"+np.str(v)+"', None)")
			exec("filenc.createDimension('j"+np.str(v)+"', None)")
			exec("var"+np.str(v)+" = filenc.createVariable(varnames[v], tp.name, ('i"+np.str(v)+"','j"+np.str(v)+"'), zlib=True)")
			exec("var"+np.str(v)+"[:,:] = data."+varnames[v])			
		elif d == 3:		
			exec("filenc.createDimension('i"+np.str(v)+"', None)")
			exec("filenc.createDimension('j"+np.str(v)+"', None)")
			exec("filenc.createDimension('k"+np.str(v)+"', None)")
			exec("var"+np.str(v)+" = filenc.createVariable(varnames[v], tp.name, ('i"+np.str(v)+"','j"+np.str(v)+"','k"+np.str(v)+"'), zlib=True)")
			exec("var"+np.str(v)+"[:,:,:] = data."+varnames[v])			
		elif d == 4:
			exec("filenc.createDimension('i"+np.str(v)+"', None)")
			exec("filenc.createDimension('j"+np.str(v)+"', None)")
			exec("filenc.createDimension('k"+np.str(v)+"', None)")
			exec("filenc.createDimension('t"+np.str(v)+"', None)")
			exec("var"+np.str(v)+" = filenc.createVariable(varnames[v], tp.name, ('i"+np.str(v)+"','j"+np.str(v)+"','k"+np.str(v)+"','t"+np.str(v)+"'), zlib=True)")			
			exec("var"+np.str(v)+"[:,:,:,:] = data."+varnames[v])

	filenc.close()
	print "File "+filename+" created"

	return




def ECCO2(varnames,latrange,lonrange,timerange):
	r"""
	Download a subset from the Estimating the Circulation and Climate of 
	the Ocean, Phase II (ECCO2): High-Resolution Global-Ocean and Sea-Ice 
	Data Synthesis (ECCO2), from the NASA OpenDap server 



	Parameters
	----------
	varnames: list
		names of the variables to retrieve e.g., ['var1','var2','var3']

	latrange: list
		latitudes of the northern and southern limits of the domain e. g., [90,-90]

	lonrange: list
		longitudes of the eastern and western limits of the domain e. g., [-180,180]	

	timerange: list
		choose from monthly (yyyymm) or daily/3-day (yyyymmdd) resolution products and the time range
		(IMPORTANT: DOWNLOAD MONTHLY, DAILY OR 3-DAY VARIABLES SET SEPARATLY)

		e.g., [initial_time,final_time] or [initial_time,final_time]


	Returns
	-------
	subset: class 
		subset.lat; subset.lon; subset.time; subset.depth
		subset.var1; subset.var2; subset.var3


	Variables 			Long Name 							Product
	--------- 			--------- 							-------
	MXLDTH 				mix layer depth						daily 

	PHIBOT				Bottom pressure anomaly 			daily

	SALT (SALT_monthly) Salinity 							3-day/monthly

	SIarea 				SeaIce fractional ice coverage		daily 

	SIheff				SeaIce effective ice thickness 		daily

	SIsnow				SeaIce snow coverage				daily 

	SSH 				Sea Surface Height 					daily

	SSS 				Sea Surface Salinity 				daily 

	SST 				Sea Surface Temperature 			daily 

	THETA (THETA_monthly)	Potential Temperature 				3-day/monthly

	UVEL (UVEL_monthly) 	Zonal velocity component 			3-day/monthly

	VVEL (VVEL_monthly)		Meridional velocity component 		3-day/monthly

	WVEL (WVEL_monthly)		Vertical velocity component 		3-day/monthly

	oceFWflux_daily		Fresh water flux into the ocean 	daily

	oceQnet_daily 		Net heat flux 						daily 

	oceQsw_daily		Short-wave radiation flux 			daily 

	oceTAUX_daily		Zonal wind stress 					daily 

	oceTAUY_daily 		Meridional wind stress 				daily 

	"""

	# Auxiliary packages
	from netCDF4 import Dataset, date2num 
	import datetime 
	import numpy as np 
	import pylab as py
	from dateutil.parser import parse



	## Time axes and depth axes
	if varnames[0] == 'SALT' or varnames[0] == 'THETA' or varnames[0] == 'UVEL'  or varnames[0] == 'VVEL'  or varnames[0] == 'WVEL':
		start  = datetime.datetime.strptime(timerange[0], "%Y%m%d")
		end  = datetime.datetime.strptime(timerange[1], "%Y%m%d")
		time = [start + datetime.timedelta(days=t) for t in xrange(0, (end-start).days,3)]

		for t in xrange(len(time)):
			time[t] = time[t].strftime('%Y%m%d')


		depth = np.array([0])



	elif varnames[0] == 'SALT_monthly' or varnames[0] == 'THETA_monthly' or varnames[0] == 'UVEL_monthly' or varnames[0] == 'VVEL_monthly'  or varnames[0] == 'WVEL_monthly':
		start  = datetime.datetime.strptime(timerange[0], "%Y%m")
		end  = datetime.datetime.strptime(timerange[1], "%Y%m")
		time = [start + datetime.timedelta(days=t) for t in xrange(0, (end-start).days,30)]

		for t in xrange(len(time)):
			time[t] = time[t].strftime('%Y%m')

		time = list(np.unique(time))
		
		depth = np.array([5.00000000e+00,   1.50000000e+01,   2.50000000e+01,
			3.50000000e+01,   4.50000000e+01,   5.50000000e+01,
			6.50000000e+01,   7.50049973e+01,   8.50250015e+01,
			9.50950012e+01,   1.05309998e+02,   1.15870003e+02,
			1.27150002e+02,   1.39740005e+02,   1.54470001e+02,
			1.72399994e+02,   1.94735001e+02,   2.22710007e+02,
			2.57470001e+02,   2.99929993e+02,   3.50679993e+02,
			4.09929993e+02,   4.77470001e+02,   5.52710022e+02,
			6.34734985e+02,   7.22400024e+02,   8.14469971e+02,
			9.09739990e+02,   1.00715503e+03,   1.10590503e+03,
			1.20553503e+03,   1.30620496e+03,   1.40915002e+03,
			1.51709497e+03,   1.63417505e+03,   1.76513501e+03,
			1.91415002e+03,   2.08403491e+03,   2.27622510e+03,
			2.49125000e+03,   2.72925000e+03,   2.99025000e+03,
			3.27425000e+03,   3.58125000e+03,   3.91125000e+03,
			4.26425000e+03,   4.64025000e+03,   5.03925000e+03,
			5.46125000e+03,   5.90625000e+03])

	else:
		start  = datetime.datetime.strptime(timerange[0], "%Y%m%d")
		end  = datetime.datetime.strptime(timerange[1], "%Y%m%d")
		time = [start + datetime.timedelta(days=t) for t in xrange(0, (end-start).days,1)]

		for t in xrange(len(time)):
			time[t] = time[t].strftime('%Y%m%d')

		depth = np.array([0])


	## Domain subset
	# Longitude
	longitude = np.arange(0.125,359.875+0.25,0.25)
	latitude  = np.arange(-89.875,89.875+0.25,0.25)

	if lonrange[0] < 0. and lonrange[1] < 0.:
		loni = lonrange[0]+360.
		lonf = lonrange[1]+360.

		jgood = (longitude>=loni) & (longitude<=lonf)
		jgood = py.find(jgood)

	elif lonrange[0] < 0. and lonrange[1] >= 0.:
		loni = lonrange[0]+360.
		lonf = lonrange[1]		

		jgood = (longitude>=loni) | (longitude<=lonf)
		jgood = py.find(jgood)

		shift = len(longitude[longitude>=loni])
		jgood = np.roll(jgood,shift)

	else:
		loni = lonrange[0]
		lonf = lonrange[1]			
		jgood = (longitude>=loni) | (longitude<=lonf)
		jgood = py.find(jgood)	


	# Latitude 
	igood = (latitude>=latrange[0]) & (latitude<=latrange[1])
	igood = py.find(igood)


	# Latitude and longitude axis
	lat = latitude[igood]
	lon = longitude[jgood]; lon[lon>180] = lon[lon>180]-360.


	igood,jgood = np.meshgrid(igood,jgood)

	# OpenDap server 
	url  = 'http://ecco2.jpl.nasa.gov:80/opendap/data1/cube/cube92/lat_lon/quart_90S_90N/'

	# Retrieving the variables from data
	var = []
	def download(name,time):

		print time 
		if name == 'SALT' or name == 'THETA' or name == 'UVEL' or name == 'VVEL' or name == 'WTHETA' or name == 'SALT_monthly' or name == 'THETA_monthly' or name == 'UVEL_monthly' or name == 'VVEL_monthly'  or name == 'WVEL_monthly':
			if len(name)>6:
				filename = url+name+'.nc/'+name[:-8]+'.1440x720x50.'+time+'.nc'
				name = name[:-8]
			else:
				filename = url+name+'.nc/'+name+'.1440x720x50.'+time+'.nc'

			data = Dataset(filename)
			data1 = data[name][:].squeeze()[:,igood,jgood]

		else:
			filename = url+name+'.nc/'+name+'.1440x720.'+time+'.nc'
			data = Dataset(filename,'r')
			data1 = data[name][:].squeeze()[igood,jgood]

		data.close()

		return data1 

	for k in xrange(len(varnames)):
		print varnames[k]
		rvar = []
		for t in xrange(len(time)):
			rvar.append(download(varnames[k],time[t]))
		
		var.append(np.array(rvar))


	# Adding variables to the subset
	for t in xrange(len(time)):

		if varnames[0][-7:] == 'monthly':
			time[t] = parse(time[t]+'01')
		else:
			time[t] = parse(time[t])

	time = np.array(time) 
	time = date2num(time, 'days since 1992-1-1 00:00:00')
	var.append(time); var.append(lat); var.append(lon); var.append(depth)
	varnames.append('time'); varnames.append('lat'); varnames.append('lon'); varnames.append('depth') 
	subset = class_builder(var,varnames)

	return subset   



def OFES(varnames,latrange,lonrange,timerange,product='OFES_NCEP_RUN/MONTHLY_3D/'):
	r"""
	Download a subset from the EOGCM for the Earth Simulator, from the 
	JAMSTEC OpenDap server 

	Please check the products and variables names in the http://www.jamstec.go.jp
	website. Since this is a beta version I can not guarantee that it will work in 
	all variables  

	Parameters
	----------
	varnames: list
		names of the variables to retrieve e.g., ['var1','var2','var3']

	latrange: list
		latitudes of the northern and southern limits of the domain e. g., [90,-90]

	lonrange: list
		longitudes of the eastern and western limits of the domain e. g., [0,360]	

	timerange: list
		time period of interest (yyyymmdd) e. g., [19920101,20100510]


	Returns
	-------
	subset: class 
		subset.lat; subset.lon; subset.time; subset.topo
		subset.var1; subset.var2; subset.var3

	"""
	from netCDF4 import Dataset, date2num 
	import datetime 
	import numpy as np 
	import pylab as py
	from dateutil.parser import parse


	path = 'http://www.jamstec.go.jp/esc/fes/dods/OFES/'

	# local depth (topography) and depth of the deepest layer
	filename1 = 'TOPOG/ht' 

	data1 = Dataset(path+filename1,'r')

	longitude = data1.variables.pop('lon'); longitude = longitude[:].squeeze()
	latitude  = data1.variables.pop('lat'); latitude = latitude[:].squeeze()

	loni = lonrange[0]				
	lonf = lonrange[-1]				

	# Coordinates for topography
	latgood1 = (latitude>=latrange[0]) & (latitude<=latrange[1])
	longood1 = (longitude>=loni) & (longitude<=lonf)


	# Topography 
	print "Retrieving model topography"
	topo = data1.variables.pop('ht')
	topo = topo[0,0,latgood1,:]; topo = topo[:,longood1]
	data1.close()



	## Retrieving variables
	# Time range 
	tstart  = datetime.datetime.strptime(timerange[0], "%Y%m%d")
	tend  = datetime.datetime.strptime(timerange[-1], "%Y%m%d")

	tstart = date2num(tstart, 'days since 1-1-1 00:00:00')
	tend = date2num(tend, 'days since 1-1-1 00:00:00')


	var = []
	for k in xrange(len(varnames)):
		print 'Retrieving '+varnames[k]

		filename = path+product+varnames[k]
		data = Dataset(filename,'r')

		# Time 
		if k == 0:
			time = data.variables.pop('time'); time = time[:].squeeze()
			pressure = data.variables.pop('lev'); pressure = pressure[:].squeeze()

			latgg = data.variables.pop('lat'); longg = data.variables.pop('lon')

			tgood = py.find((time >= tstart) & (time <= tend))
			pgood = py.find(pressure<20000.)
			latgood = py.find((latitude>=latrange[0]) & (latitude<=latrange[1]))
			longood = py.find((longitude>=loni) & (longitude<=lonf))

			time = time[(time >= tstart) & (time <= tend)]
		else:
			latgg = data.variables.pop('lat'); longg = data.variables.pop('lon')
			timeg = data.variables.pop('time'); 
			pressureg = data.variables.pop('lev');


		rvar = np.ones((tgood.shape[0],pgood.shape[0],latgood.shape[0],longood.shape[0]))*np.nan
		rvar1 = data.variables.pop(varnames[k])

		count = 0.
		for j in xrange(longood.shape[0]):
			for i in xrange(latgood.shape[0]):
				rvar[:,:,i,j] = rvar1[tgood,:,latgood[i],longood[j]]
				count+=1.


		if varnames[k] == 'temp':
			varnames[k] = 'theta'
		
		var.append(rvar)
		data.close()


	lat = latitude[(latitude>=latrange[0]) & (latitude<=latrange[1])]
	lon = longitude[(longitude>=loni) & (longitude<=lonf)]


	var.append(time); var.append(lat); var.append(lon); var.append(topo); var.append(pressure)
	varnames.append('time'); varnames.append('lat'); varnames.append('lon'); varnames.append('topo'); varnames.append('pressure')

	subset = class_builder(var,varnames)
	return subset



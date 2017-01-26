### Auxiliary classes and functions 
class class_builder(object):
	def __init__(self,var,varnames):
		for i in xrange(len(var)):

			exec("self."+varnames[i]+" = var[i]")
			
		self.varnames = varnames


class load_netcdf(object):
	def __init__(self,filename):

		from netCDF4 import Dataset 

		data = Dataset(filename,'r')
		varnames = data.variables.keys()

		for i in xrange(len(varnames)):

			exec("self."+varnames[i]+" = data[varnames[i]][:]")


		data.close()



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

		IMPORTANT: FOR MONTHLY DATA THE TIME RANGE MUST BE MAXIMUM OF 4 YEARS

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
		time = [start + datetime.timedelta(days=t) for t in xrange(0, (end-start).days+3,3)]

		for t in xrange(len(time)):
			time[t] = time[t].strftime('%Y%m%d')


		depth = np.array([0])



	elif varnames[0] == 'SALT_monthly' or varnames[0] == 'THETA_monthly' or varnames[0] == 'UVEL_monthly' or varnames[0] == 'VVEL_monthly'  or varnames[0] == 'WVEL_monthly':
		start  = datetime.datetime.strptime(timerange[0], "%Y%m")
		end  = datetime.datetime.strptime(timerange[1], "%Y%m")
		time = [start + datetime.timedelta(days=t) for t in xrange(0, (end-start).days+31,31)]

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
		time = [start + datetime.timedelta(days=t) for t in xrange(0, (end-start).days+1,1)]

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



def OFES(filename,user,passw,varnames,latrange,lonrange,zrange,timerange,product='ncep_0.1_global_3day',nz=1):
	r"""
	Download and save a subset from the OGCM for the Earth Simulator, from JAMSTEC 
	The data is accessed through the APDRC OpenDap server 
	(http://apdrc.soest.hawaii.edu/dods/esc_only/OfES/) 

	This version only will correctly retrieve the 3 days snapshots of the 0.1 degree resolution
	experiments  

	Parameters
	----------
	filename: char
		name of the file containing the OFES output

	varnames: list
		names of the variables to retrieve e.g., ['var1','var2','var3']

	latrange: list
		latitudes of the southern and northern limits of the domain [j_south,j_north]
		i. e., position of the array np.arange(LATmin,LATmax,0.1)

	lonrange: list
		longitude of the western and eastern  limits of the domain [i_east,i_west]
		i. e., position of the array np.arange(LONmin,LONmax,0.1)	

	timerange: list
		time period of interest (yyyymmdd) e. g., [19920101,20100510]

	nz: integer
		number of z-levels to be retrieved from 3D variables

	zrange: list
		Z-levels of interested. IMPORTANT: THE CONNECTION WITH THE SERVER MIGHT BE INTERRUPTED WITH
		IF YOU TRY TO RETRIEVE TOO MUCH DATA AT ONCE. THEREFORE IF YOU INTEND TO GET MORE THAN 10 
		VERTICAL LEVELS CONSIDER DO THE FOLLOWING: 

			e.g., zrange = [[0,10],[11,21],[22,32],[33,43],[44,53]] or zrange = [[0,10],[48,53]], or ... 

		if z-levels < 15: e. g., zrange = [[0,4]]

	Returns
	-------
	subset: class 
		subset.var1; subset.var2; subset.var3; subset.time (if applicable)
	"""

	from netCDF4 import Dataset, date2num, num2date 
	import datetime 
	import numpy as np 
	import pylab as py

	# URL 
	url = 'http://'+user+':'+passw+'@apdrc.soest.hawaii.edu/dods/esc_only/OfES/'+product+'/'


	## Finding the time indexes
	# 2d
	start_2d  = datetime.datetime.strptime('19500101', "%Y%m%d")
	start_2d = date2num(start_2d, 'days since 1-1-1 00:00:00')

	end_2d  = datetime.datetime.strptime('20131229', "%Y%m%d")
	end_2d = date2num(end_2d, 'days since 1-1-1 00:00:00')

	# 2d time axis
	time_2d = np.arange(start_2d,end_2d+3,3)


	# 3d
	start_3d  = datetime.datetime.strptime('19800103', "%Y%m%d")
	start_3d = date2num(start_3d, 'days since 1-1-1 00:00:00')

	end_3d  = datetime.datetime.strptime('20131229', "%Y%m%d")
	end_3d = date2num(end_3d, 'days since 1-1-1 00:00:00')

	# 3d time axis
	time_3d = np.arange(start_3d,end_3d+3,3)


	# time range conversion
	trange = [0,0] 
	trange[0] = datetime.datetime.strptime(timerange[0], "%Y%m%d")
	trange[-1] = datetime.datetime.strptime(timerange[-1], "%Y%m%d")

	trange[0] = date2num(trange[0], 'days since 1-1-1 00:00:00')
	trange[-1] = date2num(trange[-1], 'days since 1-1-1 00:00:00')


	# Indexes
	t_2d = [py.find(time_2d >= trange[0])[0]]
	t_2d.append(py.find(time_2d < trange[-1])[-1])

	t_3d = [py.find(time_3d >= trange[0])[0]]
	t_3d.append(py.find(time_3d < trange[-1])[-1])


	# Initializing variables
	var = []

	# Retrieving the outputs
	for varname in varnames:
		VARNAMES = varname+' '

		if varname == u'eta' or varname == u'hflx' or varname == u'sflx' or varname == 'taux' or varname == 'tauy' or varname == 'convU' or varname == 'hblt' or varname == 'hmxl':
			data  = Dataset(url+varname,'r')
			var.append(data[varname][t_2d[0]:t_2d[-1]+1,0,latrange[0]:latrange[-1]+1,lonrange[0]:lonrange[-1]+1].squeeze())
			data.close()
		elif varname == u'ht' or varname == 'kmt':
			data  = Dataset(url+varname,'r')
			var.append(data[varname][0,0,latrange[0]:latrange[-1]+1,lonrange[0]:lonrange[-1]+1].squeeze())
			data.close()
		elif varname == 'salt':
			for ZR in zrange:
				data  = Dataset(url+varname,'r')
				VAR = data['salinity'][t_3d[0]:t_3d[-1]+1,ZR[0]:ZR[-1]+1,latrange[0]:latrange[-1]+1,lonrange[0]:lonrange[-1]+1].squeeze() 
				data.close()

				if ZR[0] == zrange[0][0]:
					NVAR = np.ones((VAR.shape[0],nz,VAR.shape[2],VAR.shape[3]))*np.nan
					z = np.ones(54)*np.nan


				NVAR[:,ZR[0]:ZR[-1]+1,:,:] = VAR.copy()
				z[ZR[0]:ZR[-1]+1] = 1

				del(VAR)

			igood = ~np.isnan(z) 
			var.append(NVAR[:,igood,:,:].copy())
			del(NVAR)
			del(z)
			del(igood)

			# var.append(data['salinity'][t_3d[0]:t_3d[-1]+1,zrange[0]:zrange[-1]+1,latrange[0]:latrange[-1]+1,lonrange[0]:lonrange[-1]+1].squeeze())
		else:
			for ZR in zrange:
				data  = Dataset(url+varname,'r')
				VAR = data[varname][t_3d[0]:t_3d[-1]+1,ZR[0]:ZR[-1]+1,latrange[0]:latrange[-1]+1,lonrange[0]:lonrange[-1]+1].squeeze() 
				data.close()

				if ZR[0] == zrange[0][0]:
					NVAR = np.ones((VAR.shape[0],nz,VAR.shape[2],VAR.shape[3]))*np.nan
					z = np.ones(54)*np.nan


				NVAR[:,ZR[0]:ZR[-1]+1,:,:] = VAR.copy()
				z[ZR[0]:ZR[-1]+1] = 1

				del(VAR)

			igood = ~np.isnan(z) 
			var.append(NVAR[:,igood,:,:].copy())
			del(NVAR)
			del(z)
			del(igood)

			# var.append(data[varname][t_3d[0]:t_3d[-1]+1,zrange[0]:zrange[-1]+1,latrange[0]:latrange[-1]+1,lonrange[0]:lonrange[-1]+1].squeeze())

	# Organize outputs in a class
	subset = class_builder(var,varnames)

	write_netcdf(filename,
		subset.varnames,subset,'OFES '+VARNAMES+' from'+timerange[0]+' to '+timerange[-1])



## Extract metadata from output.b ASCII files
def get_metadata(filename):
  r"""
    Get variables names, order (i.e., variable index), 
    vertical position, and dimensions from HYCOM outputs.b files

    Obs: vertical position = 0: mix layer or surface
         vertical position > 0: subsurface layer 

    Function based on Rafael Goncalves readHYCOM.py

    Author: Tiago Bilo
  """

  import numpy as np

  # Open ASCII outputs.b file
  f = open(filename,'r')

  # Initialize variables
  varnames = []
  k = []
  order = []

  ivar = 0
  idm = 0

  # Extract information from ASCII file
  for line in f.readlines():

    # Zonal axis dimension  
    if 'idm' in line: 
      idm = np.int(line.split()[0])
    
    # Meridional axis dimension
    elif 'jdm' in line:
      jdm = np.int(line.split()[0])

    # Variables names, order and vertical position
    elif ('=' in line) & (idm>0):
      ivar += 1  

      elements = line.split()

      # Info
      varnames.append(elements[0])
      k.append(np.int(elements[4]))
      order.append(ivar)
            
  # Close file
  f.close()

  return varnames,k,order,idm,jdm



## Extract variables from outputs.a binary files 
def get_variables(filename,var):

  r"""
    Retrieve variables from HYCOM outputs.a binary files. 

    Input:
      filename = 'path/outputs.a'
      var = ['var_name1','var_name2', ... 'var_nameN']

    
    For retrieving all variables

      var = 'all' 

    Function based on Rafael Goncalves readHYCOM.py

    Author: Tiago Bilo
  """

  import numpy as np
  import math

  # Extract metadata
  varnames,k,order,idm,jdm = get_metadata(filename[:-2]+'.b')


  # Total number of bytes offset for each variables. Each HYCOM record is padded to 16 Kb, 
  # or 4096*32-bit words.
  nn  = 4*np.int(math.ceil(float(jdm*idm) / 4096) * 4096)


  # Open file
  fi  = open(filename,'rb')



  # In case the user wants all variable
  if var == 'all':
    var = np.unique(varnames)


  # Variables loop
  variables = []

  for i in xrange(len(var)):

    # Find indexes of the variable of interest
    indexes = []

    for ii in xrange(len(varnames)):
      if varnames[ii] == var[i]:
        indexes.append(ii)


    # Extract variables
    var_temp = np.ones((len(indexes),jdm,idm))*np.nan

    for ii in xrange(len(indexes)):

      # Bytes offset
      off = (order[indexes[ii]]-1)*nn

      # Data stream position relative to the beginning of the file
      fi.seek(off,0)

      # Retrieve and reorganize data
      data = np.fromfile(fi,'>f',jdm*idm)
      data[data>1e10] = np.nan
      
      var_temp[ii,:,:] = data.reshape((jdm,idm))

      del(data)


    variables.append(var_temp.squeeze())

  # Close file
  fi.close()

  return variables


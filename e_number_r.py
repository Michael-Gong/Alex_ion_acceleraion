import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec

import multiprocessing as mp


######## Constant defined here ########
pi        =     3.1415926535897932384626
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 20,  
        }  


if __name__ == '__main__':
  start   =  1 # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
  
  x_start   = 15.0 # micron 
  x_stop    = 35.0 # micron 
  R_size    = 2.0  # micron
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  ######### Script code drawing figure ################
  #dir = ['cannon_a190_n15', 'cannon_a190_n20', 'cannon_a190', 'cannon_a190_n40', 'cannon_a190_n50', 'cannon_a190_n60']
  dir = ['uniform_a190_n15', 'uniform_a190_n20', 'uniform_a190_n30', 'uniform_a190_n40', 'uniform_a190_n50', 'uniform_a190_n60'] 

  for i in range(np.size(dir)):
      from_path = './'+dir[i]+'/' 
      to_path = './'+dir[i]+'/' 

      e_number  = np.zeros((stop-start+step)/step)
      e_time    = np.zeros((stop-start+step)/step)
      name = 'Electron'
      for n in range(start,stop+step,step):
          data = sdf.read(from_path+'q'+str(n).zfill(4)+".sdf",dict=True)
          header=data['Header']
          time=header['time']/1.0e-15
          x  = data['Grid/Grid_mid'].data[0]/1.0e-6
          y  = data['Grid/Grid_mid'].data[1]/1.0e-6
          z  = data['Grid/Grid_mid'].data[2]/1.0e-6
    
          X,Y,Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
          RR = (Y**2+Z**2)**0.5
      
          eexx = data['Derived/Number_Density/'+str.capitalize(name)].data/denunit
          x_index = ((X[:,0,0] > x_start) & (X[:,0,0]<x_stop))
          eexx = eexx[:, RR[0,:,:]<R_size]
          eexx = eexx[x_index,:]
          print(eexx.shape)
          ex = np.sum(eexx)*V
          e_time[(n-start)/step] = time
          e_number[(n-start)/step] = ex
          print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
    
    
      print('finised '+from_path)
      np.savetxt(to_path+'e_number_cylinder_r20.txt',e_number)
      np.savetxt(to_path+'e_time_cylinder_r20.txt',e_time)


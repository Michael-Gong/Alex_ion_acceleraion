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
  ######### Parameter you should set ###########
  #start   =  23  # start time
  #stop    =  30  # end time
  #step    =  1  # the interval or step
  #youwant =  ['ey','ex','ey_averaged','bz','bz_averaged'] #,'electron_en','electron_ekbar','electron_density']
  #youwant.append('Ion_ekbar')
  #youwant.append('positron_ekbar')
  #youwant.append('electron_en')
  #youwant.append('photon_en')
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px...
  
  to_path='./uniform_a190_n30/'
  from_path = to_path
 
  e_number = np.zeros(stop-start+step)
  e_series = np.zeros(stop-start+step)
  e_time   = np.zeros(stop-start+step)
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
    #### header data ####
    #data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    #header=data['Header']
    #time=header['time']
    #x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    #y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    #X, Y = np.meshgrid(x, y)
    name = 'Electron'
    data = sdf.read(from_path+'q'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    z  = data['Grid/Grid_mid'].data[2]/1.0e-6
  
    X,Y,Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
    RR = (Y**2+Z**2)**0.5
    
    eexx = data['Derived/Number_Density/'+str.capitalize(name)].data#/denunit
    print(eexx.shape)
    eexx = eexx[:,RR[0,:,:]<5.0]
    print(eexx.shape)
    eexx = eexx[abs(X[:,0,0]-30)<15.0,:]
    print(eexx.shape)
    ex = np.sum(eexx)*30.0*np.pi*5**2*1e-18
    e_number[n-start] = ex
    e_series[n-start] = n
    e_time[n-start]   = time
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
  np.savetxt(to_path+'eden_volume_box.txt',e_number)
  np.savetxt(to_path+'eden_volume_box_time.txt',e_time)
  np.savetxt(to_path+'eden_volume_box_n.txt',e_series)

    

import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec

import multiprocessing as mp


######## Constant defined here ########


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
        'style'  : 'normal',
        'color'  : 'black',  
	    'weight' : 'normal',  
        'size'   : 20,  
       }  
######### Parameter you should set ###########

def processplot(n): 
    from_path_2='./cannon_a190_e_div/'
    data = sdf.read(from_path_2+'q'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    z  = data['Grid/Grid_mid'].data[2]/1.0e-6
    X,Y,Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
    RR = (Y**2+Z**2)**0.5

    data_ne_in = data['Derived/Number_Density/E_in'].data/denunit+data['Derived/Number_Density/E_in_1'].data/denunit
    data_ne_out= data['Derived/Number_Density/E_out'].data/denunit+data['Derived/Number_Density/E_out_1'].data/denunit
    for R_lim in range(5,45,5): 
      eexx = data_ne_in[:,RR[0,:,:]<R_lim/10.]
      n_size = eexx[-1,:].size
      value_n_e_in = np.sum(eexx,axis=1)/n_size

      eexx = data_ne_out[:,RR[0,:,:]<R_lim/10.]
      n_size = eexx[-1,:].size
      value_n_e_out = np.sum(eexx,axis=1)/n_size

      value_axisx = np.linspace(15.25,44.25,60)
      value_axisy = np.linspace(15.25,44.25,60)
      value_grid = np.linspace(15,45,61)
    
      value_total_in  = np.zeros_like(value_axisy)
      value_total_out = np.zeros_like(value_axisy)
      value_num   = np.zeros_like(value_axisy)
    
      for i in range(60):
        value_total_in[i] = np.sum(value_n_e_in[(value_grid[i]<=x) & (value_grid[i+1]>x)],0)
        value_total_out[i]= np.sum(value_n_e_out[(value_grid[i]<=x) & (value_grid[i+1]>x)],0)
#        value_num[i] = np.size(work_y[(value_grid[i]<=gamma) & (value_grid[i+1]>gamma)])
        print('in-:',value_total_in[i]/(value_total_in[i]+value_total_out[i]),'; out-:',value_total_out[i]/(value_total_in[i]+value_total_out[i]))
    
    #    plt.subplot()
      in_ = value_total_in/(value_total_in+value_total_out)
      out_ = value_total_out/(value_total_in+value_total_out)
      width=0.4
      pl=plt.bar(value_axisx, in_*100, width, color='orangered',edgecolor='black',linewidth=2)
      pt=plt.bar(value_axisx, out_*100, width, bottom=in_*100, color='dodgerblue',edgecolor='black',linewidth=2)
    
      plt.xlim(15,45)
      plt.ylim(0,104)
      plt.xlabel('X [$\mu m$]',fontdict=font)
      plt.ylabel('Percentage [%]',fontdict=font)
      plt.xticks(fontsize=20); plt.yticks(fontsize=20);
      plt.legend(['n$_e$-in','n$_e$-out'],loc='best',fontsize=18)
    #plt.text(200,650,' t=400fs',fontdict=font)
    
    #plt.show()
    #lt.figure(figsize=(100,100))
      fig = plt.gcf()
      fig.set_size_inches(10.2, 8.4)
      fig.savefig(from_path_2+'e_in_out_ratio_'+str(n).zfill(4)+'_r'+str(R_lim).zfill(2)+'.png',format='png',dpi=160)
      plt.close("all")
      print('finised e_in_out_ratio_'+str(n).zfill(4)+'_r'+str(R_lim).zfill(2)+'.png')
    return 0

if __name__ == '__main__':
  start   =  15 # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=4)
  results = pool.map(processplot,inputs)
  print(results)

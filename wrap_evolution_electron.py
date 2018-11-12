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
from functools import partial

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
        'size'   : 25,  
        }  

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

  
def processdata_1(n, R_limit, L_length, dir_string): 
  from_path = dir_string
 
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
  eexx = eexx[:,RR[0,:,:]<R_limit]
  eexx = eexx[abs(X[:,0,0]-(15.+L_length/2.0))<L_length,:]
  ex = np.sum(eexx)*30.0*np.pi*5**2*1e-18
  return(n,ex)
  print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
  


if __name__ == '__main__':
  start   =  1  # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
    
  e_number = np.zeros(stop-start+step)
  e_series = np.zeros(stop-start+step)

  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=8)
  for l0 in range(5,35,5):
    for r0 in range(1,6):
      func_partial = partial(processdata_1,R_limit=r0,L_length=l0,dir_string='./cannon_a190/')  
      results = pool.map(func_partial,inputs)
      #print(results)
      for i in range(stop-start+step):
          e_number[i] = results[i][1]
          e_series[i] = results[i][0]
      plt.plot(e_series, e_number/1e19,'-r',marker='o',markersize=12, linewidth=4, label='structured',zorder=0)
    
      func_partial = partial(processdata_1,R_limit=r0,L_length=l0,dir_string='./uniform_a190_n30/')  
      results = pool.map(func_partial,inputs)
      print(results)
      for i in range(stop-start+step):
          e_number[i] = results[i][1]
          e_series[i] = results[i][0]
      plt.plot(e_series, e_number/1e19,'-k',marker='o',markersize=12, linewidth=4, label='uniform',zorder=0)
    
      plt.xlim(0,32)
    #  plt.ylim(0.,1.5)
      plt.xlabel('time snapshot [10fs]',fontdict=font)
      plt.ylabel('Number [$10^{19}$]',fontdict=font)
    #  plt.xticks([15,25,35,45],fontsize=25); plt.yticks([0,0.5,1.0,1.5],fontsize=25);
      plt.xticks(fontsize=25); plt.yticks(fontsize=25);
    #  plt.text(-100,650,' t = '++' fs',fontdict=font)
      plt.legend(loc='upper left',fontsize=20,framealpha=1.0)
      plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                    wspace=None, hspace=None)
    #  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
    
      fig = plt.gcf()
      fig.set_size_inches(10, 8.5)
      fig.savefig('./'+'wrap_evolution_electron_r'+str(r0)+'_l'+str(l0).zfill(2)+'.png',format='png',dpi=160)
      print('./'+'wrap_evolution_electron_r'+str(r0)+'_l'+str(l0).zfill(2)+'.png')
      plt.close("all")



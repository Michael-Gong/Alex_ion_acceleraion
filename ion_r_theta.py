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

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def processplot(n): 
  
  from_path = './txt_cannon/'
  to_path   = './txt_cannon/'
  
  px = np.loadtxt(from_path+'px2d_x.txt')
  py = np.loadtxt(from_path+'py2d_x.txt')
  pz = np.loadtxt(from_path+'pz2d_x.txt')
  x  = np.loadtxt(from_path+'xx2d_x.txt')
  y  = np.loadtxt(from_path+'yy2d_x.txt')
  z  = np.loadtxt(from_path+'zz2d_x.txt')

  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
  r     = (y**2+z**2)**0.5  

  gg = (px**2+py**2+pz**2+1)**0.5
  Ek = (gg-1)*1836*0.51
  
  fig = plt.figure()
  ax  = fig.add_subplot(111)
#  ax.set_ylim(0, 90)
#  ax.set_yticks([0,30,60])
#  ax.set_rlabel_position(0)

  img = ax.scatter(r[:,n-10], theta[:,n-10], c='red',  s=10., edgecolors='None', alpha=0.5)
#  fig.colorbar(img,cax=cax,label='time [fs]', ticks=[200,240,280,320,360])
  ax.set_xlim(0,5.5)
  ax.set_ylim(0.,45.)
  ax.set_xlabel(r'$r\ [\mu m]$',fontdict=font)
  ax.set_ylabel(r'$\theta\ [^o]$',fontdict=font)
  ax.tick_params(axis='x',labelsize=20)
  ax.tick_params(axis='y',labelsize=20)
  ax.set_title(r'$r-\theta_r$'+'_'+str(n), va='bottom', y=1., fontsize=20)
#  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.10, bottom=None, right=0.86, top=None,
                wspace=None, hspace=None)
#  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
#plt.show()
#lt.figure(figsize=(100,100))

  #par1.set_ylabel(r'$E_x\ [m_ec\omega/|e|]$',fontdict=font)
  #par1.yaxis.label.set_color('red')
  #par1.tick_params(axis='y', colors='red', labelsize=20)
  #par1.set_ylim(-5,12)

  fig = plt.gcf()
  fig.set_size_inches(10, 7.5)
  fig.savefig(to_path+'r_theta_'+str(n).zfill(4)+'.png',format='png',dpi=160)
  plt.close("all")
  print('finised ')

if __name__ == '__main__':
  start   =  10 # start time
  stop    =  27  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=10)
  results = pool.map(processplot,inputs)
  print(results)

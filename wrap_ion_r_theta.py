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

if __name__ == '__main__':
  start   =  10 # start time
  stop    =  27  # end time
  step    =  1  # the interval or step
  
  from_path = './cannon_a190_v484/'
  to_path   = './cannon_a190_v484_fig/'
  
  data = sdf.read(from_path+"i_tot_loc0027.sdf",dict=True)
  #grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
  px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
  gg = (px**2+py**2+pz**2+1)**0.5
  Ek = (gg-1)*1836*0.51

  part13_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
  part13_id = part13_id[ (Ek>225) & (abs(theta)<10) & (Ek<245)]

  print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))
  

  fig, ax = plt.subplots()
  serial_color = ['red','gold','lawngreen','lightskyblue','darkviolet']
  serial_n = [27,23,19,15,11]
  serial_range = np.array([ [34,41], [28,34], [21,26], [16.5,19.5], [13,16.5] ])

  serial_color = ['red']*10+['gold']*10+['lawngreen']*10+['lightskyblue']*10+['darkviolet']*10
  cmap=matplotlib.colors.ListedColormap(serial_color)
  norm = matplotlib.colors.Normalize(vmin=180, vmax=380)

  for i in range(5):
      data = sdf.read(from_path+"i_tot_loc"+str(serial_n[i]).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']/1.0e-15
      px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      grid_y = data['Grid/Particles/subset_Only_Ions0/Ion'].data[1]/1.0e-6      
      grid_z = data['Grid/Particles/subset_Only_Ions0/Ion'].data[2]/1.0e-6      
      temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data

      theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
      grid_r = (grid_y**2+grid_z**2)**0.5

      theta = theta[np.in1d(temp_id,part13_id)]
      grid_r = grid_r[np.in1d(temp_id,part13_id)]

      time_c = np.zeros_like(grid_r) + time1

      img = ax.scatter(grid_r, theta, c=time_c, norm=norm, cmap=cmap, s=0.04, edgecolors='None', alpha=0.66)
#  fig.colorbar(img,cax=cax,label='time [fs]', ticks=[200,240,280,320,360])
  ax.set_xlim(0,5.5)
  ax.set_ylim(0.,42.)
  ax.set_xlabel(r'$r\ [\mu m]$',fontdict=font)
  ax.set_ylabel(r'$\theta\ [^o]$',fontdict=font)
  ax.tick_params(axis='x',labelsize=20)
  ax.tick_params(axis='y',labelsize=20)
#  ax.set_title(r'$r-\theta_r$'+'_'+str(n), va='bottom', y=1., fontsize=20)
#  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.12, bottom=0.12, right=0.95, top=0.95,
                wspace=None, hspace=None)
#  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
#plt.show()
#lt.figure(figsize=(100,100))

  #par1.set_ylabel(r'$E_x\ [m_ec\omega/|e|]$',fontdict=font)
  #par1.yaxis.label.set_color('red')
  #par1.tick_params(axis='y', colors='red', labelsize=20)
  #par1.set_ylim(-5,12)

  fig = plt.gcf()
  fig.set_size_inches(7, 6)
  fig.savefig('./wrap_r_theta.png',format='png',dpi=160)
  plt.close("all")
  print('finised ')


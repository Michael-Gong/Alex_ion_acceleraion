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
from mpl_toolkits.mplot3d import Axes3D

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
        'size'   : 26,  
        }  

font_size = 20


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def processplot(n): 
  from_path='./cannon_a190_e_div/'  
  to_path=from_path  

#  data = sdf.read(from_path+"new_loc0027.sdf",dict=True)
#  part_id = data['Particles/ID/subset_New_particles/E_in_1'].data
#  choice = np.random.choice(range(part_id.size), int(part_id.size/1000), replace=False)
#  part_id_in = part_id[choice]
#  print('in part_id size is ',part_id_in.size,' max ',np.max(part_id_in),' min ',np.min(part_id_in))
#  
#  data = sdf.read(from_path+"new_loc0027.sdf",dict=True)
#  part_id = data['Particles/ID/subset_New_particles/E_out_1'].data
#  choice = np.random.choice(range(part_id.size), int(part_id.size/1000), replace=False)
#  part_id_out = part_id[choice]
#  print('out part_id size is ',part_id_out.size,' max ',np.max(part_id_out),' min ',np.min(part_id_out))
#  
  data = sdf.read(from_path+"new_loc"+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time1=header['time']

  grid_x_in = data['Grid/Particles/subset_New_particles/E_in_1'].data[0]/1.0e-6      
  grid_y_in = data['Grid/Particles/subset_New_particles/E_in_1'].data[1]/1.0e-6      
  grid_z_in = data['Grid/Particles/subset_New_particles/E_in_1'].data[2]/1.0e-6      

  grid_x_out = data['Grid/Particles/subset_New_particles/E_out_1'].data[0]/1.0e-6      
  grid_y_out = data['Grid/Particles/subset_New_particles/E_out_1'].data[1]/1.0e-6      
  grid_z_out = data['Grid/Particles/subset_New_particles/E_out_1'].data[2]/1.0e-6      

  fig = plt.figure()
  ax = plt.axes(projection='3d')

  reduce_fac = 100

  makersize = 0.5
  print(grid_x_in[::reduce_fac].shape,grid_y_in[::reduce_fac].shape,grid_z_in[::reduce_fac].shape)
  print(grid_x_out[::reduce_fac].shape,grid_y_out[::reduce_fac].shape,grid_z_out[::reduce_fac].shape)
#    plt.subplot()
  #normalize = matplotlib.colors.Normalize(vmin=0, vmax=20, clip=True)
  pt3d=ax.scatter(grid_x_in[::reduce_fac], grid_y_in[::reduce_fac], grid_z_in[::reduce_fac], c='salmon', s=makersize*5, edgecolors='face', alpha=0.7, marker='.')
  pt3d=ax.scatter(grid_x_out[::reduce_fac], grid_y_out[::reduce_fac], grid_z_out[::reduce_fac], c='skyblue', s=makersize, edgecolors='face', alpha=0.7, marker='.')
 #   plt.legend(loc='upper right')
  print('hellow')
  ax.set_xlim([0,55])
  ax.set_ylim([-19.,19])
  ax.set_zlim([-19.,19])
  ax.set_xlabel('\n\nX [$\mu m$]',fontdict=font)
  ax.set_ylabel('\n\nY [$\mu m$]',fontdict=font)
  ax.set_zlabel('\n\nZ [$\mu m$]',fontdict=font)
  for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(font_size)
  for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(font_size)
  for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(font_size)

  ax.grid(linestyle='--', linewidth='0.5', color='grey')
  ax.view_init(elev=45, azim=-45)

  ax.xaxis.pane.set_edgecolor('black')
  ax.yaxis.pane.set_edgecolor('black')
  ax.zaxis.pane.set_edgecolor('black')
  # Set the background color of the pane YZ
  ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
  ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
  ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

  ax.scatter(grid_x_in[::reduce_fac],grid_z_in[::reduce_fac],c='salmon',s=makersize*5, alpha=0.5,zdir='y',zs=19,marker='.')
  ax.scatter(grid_x_in[::reduce_fac],grid_y_in[::reduce_fac],c='salmon',s=makersize*5, alpha=0.5,zdir='z',zs=-19,marker='.')
  ax.scatter(grid_y_in[::reduce_fac],grid_z_in[::reduce_fac],c='salmon',s=makersize*5, alpha=0.5,zdir='x',zs=0,marker='.')

  ax.scatter(grid_x_out[::reduce_fac],grid_z_out[::reduce_fac],c='skyblue',s=makersize, alpha=0.5,zdir='y',zs=19,marker='.')
  ax.scatter(grid_x_out[::reduce_fac],grid_y_out[::reduce_fac],c='skyblue',s=makersize, alpha=0.5,zdir='z',zs=-19,marker='.')
  ax.scatter(grid_y_out[::reduce_fac],grid_z_out[::reduce_fac],c='skyblue',s=makersize, alpha=0.5,zdir='x',zs=0,marker='.')

#  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                wspace=None, hspace=None)
  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
#plt.show()
#lt.figure(figsize=(100,100))


  fig = plt.gcf()
  fig.set_size_inches(12, 10.5)
  fig.savefig(to_path+'e_in_out_3d_scatter'+str(n).zfill(4)+'.png',format='png',dpi=80)
  plt.close("all")
  print('finised '+str(n).zfill(4))
#  return 0

if __name__ == '__main__':
  start   =  3  # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=6)
  results = pool.map(processplot,inputs)
  print(results)

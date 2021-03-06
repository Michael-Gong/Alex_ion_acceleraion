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
        'size'   : 20,  
        }  

font_size = 20


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def processplot(n): 
  
  from_path = './'
  to_path   = './'
  
  data = sdf.read(from_path+"i_tot_loc0027.sdf",dict=True)
  #grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
  px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
  gg = (px**2+py**2+pz**2+1)**0.5
  Ek = (gg-1)*1836*0.51

  part13_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
  part13_id = part13_id[ (Ek>5) ]
#  part13_id = part13_id[ (Ek>225) & (abs(theta)<20) & (Ek<245)]

#  choice = np.random.choice(range(part13_id.size), 2000, replace=False)
#  part13_id = part13_id[choice]

  choice = np.random.choice(range(part13_id.size), 500000, replace=False)
  part13_id = part13_id[choice]
  print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))
  
  ######### Script code drawing figure ################
  #for n in range(start,stop+step,step):
  #### header data ####
  #data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  #header=data['Header']
  #time=header['time']
  #x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  #y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  #X, Y = np.meshgrid(x, y)
  
  data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time1=header['time']
  px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  gg = (px**2+py**2+pz**2+1)**0.5
  Ek = (gg-1)*1836*0.51
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi

  grid_x = data['Grid/Particles/subset_Only_Ions0/Ion'].data[0]/1.0e-6      
  grid_y = data['Grid/Particles/subset_Only_Ions0/Ion'].data[1]/1.0e-6      
  grid_z = data['Grid/Particles/subset_Only_Ions0/Ion'].data[2]/1.0e-6      
  temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
  grid_r =  (grid_y**2+grid_z**2)**0.5

  px = px[np.in1d(temp_id,part13_id)]
  grid_x = grid_x[np.in1d(temp_id,part13_id)]
  grid_y = grid_y[np.in1d(temp_id,part13_id)]
  grid_z = grid_z[np.in1d(temp_id,part13_id)]
  theta = theta[np.in1d(temp_id,part13_id)]
  grid_r =  grid_r[np.in1d(temp_id,part13_id)]
  Ek = Ek[np.in1d(temp_id,part13_id)]
 
  

 

  if np.size(px) == 0:
    return 0;
  Ek[Ek > 500] = 500
  Ek[0] = 0.0
  Ek[-1] = 500.0
  color_index = Ek
  fig = plt.figure()
  ax = plt.axes(projection='3d')

  makersize = 0.01
#    plt.subplot()
  #normalize = matplotlib.colors.Normalize(vmin=0, vmax=20, clip=True)
  pt3d=ax.scatter(grid_x, grid_y, grid_z, c=color_index, s=makersize*10, cmap='nipy_spectral', edgecolors='face', alpha=1.0)

  cbar=plt.colorbar(pt3d, ticks=np.linspace(np.min(color_index), np.max(color_index), 5) ,pad=0.05)
  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
  cbar.set_label(r'$E_k$'+' [MeV]',fontdict=font)
  cbar.set_clim(0,500)
#plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
#plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
#plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
#plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
 #   plt.legend(loc='upper right')
  ax.set_xlim([10,50])
  ax.set_ylim([-12.,12])
  ax.set_zlim([-12.,12])
  ax.set_xlabel('\n\nX [$\mu m$]',fontdict=font)
  ax.set_ylabel('\n\nY [$\mu m$]',fontdict=font)
  ax.set_zlabel('\n\nZ [$\mu m$]',fontdict=font)
  for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(font_size)
  for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(font_size)
  for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(font_size)

  ax.grid(linestyle='--', linewidth='0.5', color='grey')
  #ax.grid()
  ax.view_init(elev=45, azim=-45)

  ax.xaxis.pane.set_edgecolor('black')
  ax.yaxis.pane.set_edgecolor('black')
  ax.zaxis.pane.set_edgecolor('black')
  # Set the background color of the pane YZ
  ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
  ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))
  ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 1.0))

  ax.scatter(grid_x,grid_z,c=color_index,s=makersize*0.5, alpha=0.5,zdir='y',zs=12,cmap='nipy_spectral')
  ax.scatter(grid_x,grid_y,c=color_index,s=makersize*0.5, alpha=0.5,zdir='z',zs=-12,cmap='nipy_spectral')
  ax.scatter(grid_y,grid_z,c=color_index,s=makersize*0.5, alpha=0.5,zdir='x',zs=10,cmap='nipy_spectral')


#  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                wspace=None, hspace=None)
  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font, y=1.08)
#plt.show()
#lt.figure(figsize=(100,100))


  fig = plt.gcf()
  fig.set_size_inches(12, 8.5)
  fig.savefig(to_path+'proton_3d_scatter_Ek'+str(n).zfill(4)+'.png',format='png',dpi=320)
  plt.close("all")
  print('finised '+str(n).zfill(4))
  return 0

if __name__ == '__main__':
  start   =  1  # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=10)
  results = pool.map(processplot,inputs)
  print(results)

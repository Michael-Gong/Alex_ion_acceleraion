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

c_red = matplotlib.colors.colorConverter.to_rgba('red')
c_yellow= matplotlib.colors.colorConverter.to_rgba('yellow')
c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
cmap_wyr = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_white_trans,c_yellow,c_red],128)

c_blue = matplotlib.colors.colorConverter.to_rgba('dodgerblue')
c_green= matplotlib.colors.colorConverter.to_rgba('lime')
c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
cmap_wbg = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_white_trans,c_blue,c_green],128)

if __name__ == '__main__':

  n0=30.0
  R=1.8e-6
  L=15e-6
  Ntot = np.pi*R*R*L*n0*denunit
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  
  start   =  3 # start time
  stop    =  51  # end time
  step    =  2  # the interval or step

  for dir_n in ['v484_2']:
    to_path   = './cannon_a190_'+dir_n+'_fig/'
    from_path = './cannon_a190_'+dir_n+'/'  
    if not os.path.exists(to_path):
          os.mkdir(to_path)
    n0=100.0
#    n0=float(dir_n[-2:])
    R=1.8e-6
    L=15e-6
    Ntot = np.pi*R*R*L*n0*denunit
    V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  ######### Script code drawing figure ################
    for n in range(start,stop+step,step):
      weight = V*denunit*n0/10.0 
      data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']
      if 'Particles/Px/subset_Only_Energetic_Electrons/Electron' not in data:
          continue
      px = data['Particles/Px/subset_Only_Energetic_Electrons/Electron'].data/(m0*v0)
      py = data['Particles/Py/subset_Only_Energetic_Electrons/Electron'].data/(m0*v0)
      pz = data['Particles/Pz/subset_Only_Energetic_Electrons/Electron'].data/(m0*v0)
      gg = (px**2+py**2+pz**2+1)**0.5
      theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
    
      grid_x = data['Grid/Particles/subset_Only_Energetic_Electrons/Electron'].data[0]/1.0e-6      
      grid_y = data['Grid/Particles/subset_Only_Energetic_Electrons/Electron'].data[1]/1.0e-6      
      grid_z = data['Grid/Particles/subset_Only_Energetic_Electrons/Electron'].data[2]/1.0e-6      
      temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data

      grid_r = (grid_y**2+grid_z**2)**0.5    
      px = px[grid_r<1.5]
      py = py[grid_r<1.5]
      pz = pz[grid_r<1.5]
      grid_x = grid_x[grid_r<1.5]
    
      if np.size(px) == 0:
        continue;
    
      color_index = abs(theta)
    
#      plt.subplot(2,1,1)
#      print('hehe1')
#      weight = np.zeros_like(grid_x)+weight
#      print(max(weight),min(weight))
#    #    plt.subplot()
#      plt.hist2d(grid_x, px, bins=(100, 100), range=[[10,55],[0,2000]], cmap='jet', weights=weight, normed=False, norm=colors.LogNorm(vmin=1e3, vmax=1e8))
#    #  cbar=plt.colorbar(pad=0.3)
#    #  cbar.set_label(r'$dN/(dxdp_x)$'+' [A.U.]',fontdict=font)
#    #  cbar.ax.tick_params(labelsize=20)
#    #  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)
#    #plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
#    #plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
#    #plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
#    #plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
#     #   plt.legend(loc='upper right')
#      plt.xlim(0,55)
##      plt.ylim(0.,1.5)
#      plt.xlabel('X [$\mu m$]',fontdict=font)
#      plt.ylabel('$p_x$ [m$_e$c]',fontdict=font)
#      plt.xticks([10,20,30,40,50],fontsize=25); 
#      plt.yticks(fontsize=25);
#    #  plt.text(-100,650,' t = '++' fs',fontdict=font)
 
      plt.subplot(3,1,1)
      plt.scatter(grid_x, px, c='crimson', s=0.2, edgecolors='None', alpha=0.2)
      plt.xlabel('X [$\mu m$]',fontdict=font)
      plt.ylabel('$p_x$ [m$_e$c]',fontdict=font)
      plt.xlim(0,55)
      plt.ylim(0,1850)
      plt.xticks([0,10,20,30,40,50],fontsize=25); 
      plt.yticks([500,1000,1500],fontsize=25);

      plt.subplot(3,1,2)
      plt.scatter(grid_x, py, c='limegreen', s=0.2, edgecolors='None', alpha=0.2)
      plt.xlabel('X [$\mu m$]',fontdict=font)
      plt.ylabel('$p_y$ [m$_e$c]',fontdict=font)
      plt.xlim(0,55)
      plt.ylim(-300,300)
      plt.xticks([0,10,20,30,40,50],fontsize=25); 
      plt.yticks([-200,0,200],fontsize=25);

      plt.subplot(3,1,3)
      plt.scatter(grid_x, pz, c='dodgerblue', s=0.2, edgecolors='None', alpha=0.2)
      plt.xlabel('X [$\mu m$]',fontdict=font)
      plt.ylabel('$p_x$ [m$_e$c]',fontdict=font)
      plt.xlim(0,55)
      plt.ylim(-300,300)
      plt.xticks([0,10,20,30,40,50],fontsize=25); 
      plt.yticks([-200,0,200],fontsize=25);


      plt.subplots_adjust(left=0.15, bottom=0.15, right=0.96, top=0.99,
                    wspace=0.2, hspace=0.2)
    #plt.show()
    #lt.figure(figsize=(100,100))
   


 
    
      fig = plt.gcf()
      fig.set_size_inches(8.5, 20.)
      fig.savefig(to_path+'pxpypz_x_e_'+str(n).zfill(4)+'.png',format='png',dpi=160)
      plt.close("all")
      print(to_path+'pxpypz_x_e_'+str(n).zfill(4)+'.png')

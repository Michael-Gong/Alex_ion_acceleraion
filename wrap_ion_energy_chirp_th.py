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
        'size'   : 25,  
        }  

font2 = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 18,  
        }  

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

if __name__ == '__main__':
  start   =  3  # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
    
  
  from_path = './cannon_a190/'
  to_path   = './cannon_a190/'
  
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
  

  fig, ax1 = plt.subplots()
  serial_color = ['red','gold','lawngreen','lightskyblue','darkviolet']
  serial_n = [27,23,19,15,11]
  serial_range = np.array([ [34,41], [28,34], [21,26], [16.5,19.5], [13,16.5] ])
  ax2 = ax1.twinx()

  serial_color = ['red']*10+['gold']*10+['lawngreen']*10+['lightskyblue']*10+['darkviolet']*10
  cmap=matplotlib.colors.ListedColormap(serial_color)
  norm = matplotlib.colors.Normalize(vmin=180, vmax=380)

  for i in range(5):
      data = sdf.read(from_path+"i_tot_loc"+str(serial_n[i]).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']/1.0e-15
      px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      grid_x = data['Grid/Particles/subset_Only_Ions0/Ion'].data[0]/1.0e-6      
      temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
     
      px = px[np.in1d(temp_id,part13_id)]
      grid_x = grid_x[np.in1d(temp_id,part13_id)]

      time_c = np.zeros_like(grid_x) + time1

      img = ax1.scatter(grid_x, px, c=time_c, norm=norm, cmap=cmap, s=0.03, edgecolors='None', alpha=0.66)

      x  = np.loadtxt(from_path+'ex_lineout_x.txt')
      ex = np.loadtxt(from_path+'ex_lineout_r15_'+str(serial_n[i]).zfill(4)+'.txt')
      ex = ex[ (x>= serial_range[i][0]) & (x<=serial_range[i][1])] 
      x = x[ (x>= serial_range[i][0]) & (x<=serial_range[i][1])] 
      ax2.plot(x,ex, '--', linewidth=3, color=serial_color[45-i*10], label="Ex")

#plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
#plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
#plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
#plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
 #   plt.legend(loc='upper right')
  cax = fig.add_axes([0.25,0.8,0.5,0.02])
  cbar = fig.colorbar(img,cax=cax,label='time [fs]', ticks=[200,240,280,320,360], orientation='horizontal')
  cbar.set_label('time [fs]',fontdict=font2)
  cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),fontsize=18)

  ax1.set_xlim(10,45)
  ax1.set_ylim(0.,1.)
  ax1.set_xlabel('X [$\mu m$]',fontdict=font)
  ax1.set_ylabel('$p_x$ [m$_i$c$^2$]',fontdict=font)
  ax1.tick_params(axis='x',labelsize=25)
  ax1.tick_params(axis='y',labelsize=25)
  ax1.grid(linestyle='--', linewidth=0.4, color='grey')

  ax2.set_ylim(0.,25)
  ax2.set_ylabel('$E_x$ [$m_ec\omega/|e|$]',fontdict=font,color='r')
  ax2.tick_params(axis='y',labelsize=25,colors='r')
#  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.12, bottom=None, right=0.90, top=None,
                wspace=None, hspace=None)
#  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
#plt.show()
#lt.figure(figsize=(100,100))

  #par1.set_ylabel(r'$E_x\ [m_ec\omega/|e|]$',fontdict=font)
  #par1.yaxis.label.set_color('red')
  #par1.tick_params(axis='y', colors='red', labelsize=20)
  #par1.set_ylim(-5,12)

  fig = plt.gcf()
  fig.set_size_inches(10, 8.0)
  #fig.savefig(to_path+'ion_energy_chirp.png',format='png',dpi=160)
  fig.savefig('./wrap_ion_energy_chirp_th.png',format='png',dpi=160)
  plt.close("all")
  print('finised ')

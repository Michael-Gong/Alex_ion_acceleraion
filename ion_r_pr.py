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
  
  from_path = './cannon_a190/'
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
#  part13_id = part13_id[ (Ek>225) & (abs(theta)<10) & (Ek<245)]
  part13_id = part13_id[ (Ek>234.8) & (abs(theta)<10) & (Ek<235.2)]
  print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))
  
  fig = plt.figure()
  ax  = fig.add_subplot(111,polar=True)
  ax.set_ylim(0, 90)
  ax.set_yticks([0,30,60])
  ax.set_rlabel_position(0)

  serial_color = ['red']*10+['gold']*10+['lawngreen']*10+['lightskyblue']*10+['darkviolet']*10
  cmap=matplotlib.colors.ListedColormap(serial_color)
  norm = matplotlib.colors.Normalize(vmin=100, vmax=400)

  serial_n = [11,15,19,23,27] # [27,23,19,15,11]
#  serial_range = np.array([ [34,42], [28,34], [21,26], [16,21], [12,17] ])
  for i in ['ok']:
      i = n
      print('start ',i)
      data = sdf.read(from_path+"i_tot_loc"+str(i).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']/1.e-15
      px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
      print('in ',i)
      px = px[np.in1d(temp_id,part13_id)]
      py = py[np.in1d(temp_id,part13_id)]
      pz = pz[np.in1d(temp_id,part13_id)]
      print('finish ',i,' px,py,pz')
      phi    = (np.arctan2(py,pz) + np.pi)/np.pi*180
      theta  = np.arctan2((py**2+pz**2)**0.5,px)/np.pi*180
      time_c = np.zeros_like(phi) + time1

      img = ax.scatter(phi, theta, c=time_c, norm=norm, cmap=cm.nipy_spectral, s=5., edgecolors='None', alpha=0.5)
    #plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
    #plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
    #plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
    #plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
     #   plt.legend(loc='upper right')
      cax = fig.add_axes([0.8,0.1,0.02,0.3])
      fig.colorbar(img,cax=cax,label='time [fs]', ticks=[200,240,280,320,360])
      #ax.set_xlim(10,50)
      #ax.set_ylim(0.,1.)
      ax.set_xlabel(r'$\phi\ [^o]$',fontdict=font)
    #  ax.set_ylabel(r'$\theta\ [^o]$',fontdict=font)
      ax.tick_params(axis='x',labelsize=20)
      ax.tick_params(axis='y',labelsize=20)
      ax.set_title('proton_angular_time='+str(time1), va='bottom', y=1., fontsize=20)
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
      fig.savefig(to_path+'single+proton_theta_phi_'+str(i).zfill(4)+'.png',format='png',dpi=160)
      plt.close("all")
      print('finised ')

if __name__ == '__main__':
  start   =  10 # start time
  stop    =  30  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=10)
  results = pool.map(processplot,inputs)
  print(results)

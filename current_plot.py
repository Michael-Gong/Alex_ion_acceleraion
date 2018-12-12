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
jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',  
        'color'  : 'black',  
        'weight' : 'normal',  
        'size'   : 20,  
        }  


if __name__ == '__main__':
  start   =  23 # start time
  stop    =  23  # end time
  step    =  1  # the interval or step
    
  
  youwant = ['jx_averaged','jx','jy_averaged','jy','jz_averaged','jz']
  #youwant =  ['ey','ex','ey_averaged','bz','bz_averaged'] #,'electron_en','electron_ekbar','electron_density']
  #youwant.append('Ion_ekbar')
  #youwant.append('positron_ekbar')
  #youwant.append('electron_en')
  #youwant.append('photon_en')
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px...
  
  from_path = './uniform_a190_n30/'
  to_path   = './uniform_a190_n30/'
  
  
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
  #### header data ####
  #data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  #header=data['Header']
  #time=header['time']
  #x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  #y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  #X, Y = np.meshgrid(x, y)
  
    for name in youwant:
      if (name[0:2] == 'jx') or (name[0:2] == 'jy') or (name[0:2] == 'jz'):
              data = sdf.read(from_path+'current'+str(n).zfill(4)+".sdf",dict=True)
              header=data['Header']
              time=header['time']
              x  = data['Grid/Grid_mid'].data[0]/1.0e-6
              y  = data['Grid/Grid_mid'].data[1]/1.0e-6
              X, Y = np.meshgrid(x, y)
              eexx = data['Current/'+str.capitalize(name)].data/jalf
              n3d=len(eexx[0,0,:])
              ex = (eexx[:,:,n3d//2-1]+eexx[:,:,n3d//2])/2 
              if np.min(ex.T) == np.max(ex.T):
                  continue
              eee=30. #np.max([-np.min(ex.T),np.max(ex.T)])
              #if (name == 'ex'):
              #    eee = 50
              #elif  (name == 'ex_averaged'):
              #    eee = 30
              #elif (name == 'ey') or (name == 'bz'):
              #    eee = 380 
              #elif (name == 'ey_averaged') or (name == 'ez_averaged'):
              #    eee = 30
              levels = np.linspace(-eee, eee, 40)
              plt.contourf(X, Y, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-eee, vmax=eee), cmap=cm.jet)
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
              cbar.set_label('j [$J_{alfven}/\lambda^2$]',fontdict=font)
              cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=15)        
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Y [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);
              plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
              fig = plt.gcf()
              fig.set_size_inches(12, 7)
              fig.savefig(to_path+name+'cut_z=0'+str(n).zfill(4)+'.png',format='png',dpi=100)
              plt.close("all")
      if (name[0:2] == 'jx') or (name[0:2] == 'jy') or (name[0:2] == 'jz'):
              data = sdf.read(from_path+'current'+str(n).zfill(4)+".sdf",dict=True)
              header=data['Header']
              time=header['time']
              x  = data['Grid/Grid_mid'].data[0]/1.0e-6
              z  = data['Grid/Grid_mid'].data[1]/1.0e-6
              X, Z = np.meshgrid(x, z)
              eexx = data['Current/'+str.capitalize(name)].data/jalf
              n3d=len(eexx[0,:,0])
              ex = (eexx[:,n3d//2-1,:]+eexx[:,n3d//2,:])/2 
              if np.min(ex.T) == np.max(ex.T):
                  continue
              eee=30. #np.max([-np.min(ex.T),np.max(ex.T)])
              #if (name == 'ex'):
              #    eee = 50
              #elif  (name == 'ex_averaged'):
              #    eee = 30
              #elif (name == 'ey') or (name == 'bz'):
              #    eee = 380 
              #elif (name == 'ey_averaged') or (name == 'ez_averaged'):
              #    eee = 30
              levels = np.linspace(-eee, eee, 40)
              plt.contourf(X, Z, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-eee, vmax=eee), cmap=cm.jet)
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
              cbar.set_label('j [$J_{alfven}/\lambda^2$]',fontdict=font)
              cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=15)        
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Z [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);
              plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
              fig = plt.gcf()
              fig.set_size_inches(12, 7)
              fig.savefig(to_path+name+'cut_y=0'+str(n).zfill(4)+'.png',format='png',dpi=100)
              plt.close("all")
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')


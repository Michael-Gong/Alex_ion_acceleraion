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


upper = matplotlib.cm.jet(np.arange(256))
#upper = matplotlib.cm.nipy_spectral_r(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolormap = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

    
def processplot(n): 
  
  youwant = ['jx_averaged','jy_averaged','jz_averaged']
  #youwant =  ['ey','ex','ey_averaged','bz','bz_averaged'] #,'electron_en','electron_ekbar','electron_density']
  #youwant.append('Ion_ekbar')
  #youwant.append('positron_ekbar')
  #youwant.append('electron_en')
  #youwant.append('photon_en')
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px...
  
  from_path = './cannon_a190/'
  to_path   = './cannon_a190/'
  
  
  ######### Script code drawing figure ################
  if 3>2:
  #### header data ####
  #data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  #header=data['Header']
  #time=header['time']
  #x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  #y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  #X, Y = np.meshgrid(x, y)
              from_path = './cannon_a190/'
              plt.subplot(2,2,1)
              data = sdf.read(from_path+'current'+str(n).zfill(4)+".sdf",dict=True)
              header=data['Header']
              time=header['time']
              x  = data['Grid/Grid_mid'].data[0]/1.0e-6
              y  = data['Grid/Grid_mid'].data[1]/1.0e-6
              X, Y = np.meshgrid(x, y)
              jx = data['Current/Jx_averaged'].data/jalf
              jy = data['Current/Jy_averaged'].data/jalf
              j_xy = (jx**2+jy**2)**0.5
              n3d=len(j_xy[0,0,:])
              j_xy = (j_xy[:,:,n3d//2-1]+j_xy[:,:,n3d//2])/2 
              jx   = (jx[:,:,n3d//2-1]+jx[:,:,n3d//2])/2 
              jy   = (jy[:,:,n3d//2-1]+jy[:,:,n3d//2])/2 
              if np.min(j_xy) == np.max(j_xy):
                  return 0
              tickld = [-10,-1,0,1,10]
              tick_r = np.logspace(-0.5,1.5,20); tick_r = tick_r.tolist()
              tick_l = -np.logspace(-0.5,1.5,20); tick_l = tick_l.tolist(); tick_l.reverse()
              tick_l = tick_l + tick_r
              #print(tick_l)
              plt.contourf(X, Y, jx.T, levels=tick_l, norm=mcolors.SymLogNorm(linthresh=10**(-0.5), linscale=0.2, vmin=-10**1.5,  vmax=10**1.5), cmap='RdBu')
              eee=100. #np.max([-np.min(ex.T),np.max(ex.T)])
              #levels = np.logspace(-1, 2, 40)
              #plt.contourf(X, Y, jx.T, norm=mcolors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-30.0,  vmax=30.0), cmap='RdBu_r')
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(pad=0.01,ticks=tickld)
              cbar.set_label('$j_{x}$ [$I_A/\lambda^2$]',fontdict=font)
              cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=15)
##              s_j=1
##              X = X[::s_j,::s_j]
##              Y = Y[::s_j,::s_j]
##              j_xy = j_xy.T; jx = jx.T; jy=jy.T
##              j_xy = j_xy[::s_j,::s_j]
##              jx = jx[::s_j,::s_j]
##              jy = jy[::s_j,::s_j]
##              threshold_quiver = eee/4.0
##              X = X[j_xy > threshold_quiver] 
##              Y = Y[j_xy > threshold_quiver] 
##              jx = jx[j_xy > threshold_quiver] 
##              jy = jy[j_xy > threshold_quiver] 
##              plt.quiver(X, Y, 0.1*jx, 0.1*jy, color='red')
              plt.text(35,8,'cannon cut_z=0',fontsize=24)
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Y [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);

              plt.subplot(2,2,3)
              x  = data['Grid/Grid_mid'].data[0]/1.0e-6
              z  = data['Grid/Grid_mid'].data[2]/1.0e-6
              X, Z = np.meshgrid(x, z)
              jx = data['Current/Jx_averaged'].data/jalf
              jz = data['Current/Jz_averaged'].data/jalf
              j_xz = (jx**2+jz**2)**0.5
              n3d=len(j_xz[0,:,0])
              j_xz = (j_xz[:,n3d//2-1,:]+j_xz[:,n3d//2,:])/2 
              jx   = (jx[:,n3d//2-1,:]+jx[:,n3d//2,:])/2 
              jz   = (jz[:,n3d//2-1,:]+jz[:,n3d//2,:])/2 
              if np.min(j_xz) == np.max(j_xz):
                  return 0
              tickld = [-10,-1,0,1,10]
              tick_r = np.logspace(-0.5,1.5,20); tick_r = tick_r.tolist()
              tick_l = -np.logspace(-0.5,1.5,20); tick_l = tick_l.tolist(); tick_l.reverse()
              tick_l = tick_l + tick_r
              #print(tick_l)
              plt.contourf(X, Y, jx.T, levels=tick_l, norm=mcolors.SymLogNorm(linthresh=10**(-0.5), linscale=0.2, vmin=-10**1.5,  vmax=10**1.5), cmap='RdBu')
              #plt.contourf(X, Y, jx.T, levels=tick_l, cmap='RdBu')
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(pad=0.01,ticks=tickld)
              cbar.set_label('$j_{x}$ [$I_A/\lambda^2$]',fontdict=font)
              cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=15)        
              plt.text(35,8,'cannon cut_y=0',fontsize=24)
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Z [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);



              from_path = './uniform_a190_n30/'
              plt.subplot(2,2,2)
              data = sdf.read(from_path+'current'+str(n).zfill(4)+".sdf",dict=True)
              header=data['Header']
              time=header['time']
              x  = data['Grid/Grid_mid'].data[0]/1.0e-6
              y  = data['Grid/Grid_mid'].data[1]/1.0e-6
              X, Y = np.meshgrid(x, y)
              jx = data['Current/Jx_averaged'].data/jalf
              jy = data['Current/Jy_averaged'].data/jalf
              j_xy = (jx**2+jy**2)**0.5
              n3d=len(j_xy[0,0,:])
              j_xy = (j_xy[:,:,n3d//2-1]+j_xy[:,:,n3d//2])/2 
              jx   = (jx[:,:,n3d//2-1]+jx[:,:,n3d//2])/2 
              jy   = (jy[:,:,n3d//2-1]+jy[:,:,n3d//2])/2 
              if np.min(j_xy) == np.max(j_xy):
                  return 0
              tickld = [-10,-1,0,1,10]
              tick_r = np.logspace(-0.5,1.5,20); tick_r = tick_r.tolist()
              tick_l = -np.logspace(-0.5,1.5,20); tick_l = tick_l.tolist(); tick_l.reverse()
              tick_l = tick_l + tick_r
              #print(tick_l)
              plt.contourf(X, Y, jx.T, levels=tick_l, norm=mcolors.SymLogNorm(linthresh=10**(-0.5), linscale=0.2, vmin=-10**1.5,  vmax=10**1.5), cmap='RdBu')
              #plt.contourf(X, Y, jx.T, norm=mcolors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-30.0,  vmax=30.0), cmap='RdBu_r')
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(pad=0.01,ticks=tickld)
              cbar.set_label('$j_{x}$ [$I_A/\lambda^2$]',fontdict=font)
              cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=15)        
              plt.text(35,8,'uniform cut_z=0',fontsize=24)
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Y [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);

              plt.subplot(2,2,4)
              x  = data['Grid/Grid_mid'].data[0]/1.0e-6
              z  = data['Grid/Grid_mid'].data[2]/1.0e-6
              X, Z = np.meshgrid(x, z)
              jx = data['Current/Jx_averaged'].data/jalf
              jz = data['Current/Jz_averaged'].data/jalf
              j_xz = (jx**2+jz**2)**0.5
              n3d=len(j_xz[0,:,0])
              j_xz = (j_xz[:,n3d//2-1,:]+j_xz[:,n3d//2,:])/2 
              jx   = (jx[:,n3d//2-1,:]+jx[:,n3d//2,:])/2 
              jz   = (jz[:,n3d//2-1,:]+jz[:,n3d//2,:])/2 
              if np.min(j_xz) == np.max(j_xz):
                  return 0
              tickld = [-10,-1,0,1,10]
              tick_r = np.logspace(-0.5,1.5,20); tick_r = tick_r.tolist()
              tick_l = -np.logspace(-0.5,1.5,20); tick_l = tick_l.tolist(); tick_l.reverse()
              tick_l = tick_l + tick_r
              #print(tick_l)
              plt.contourf(X, Y, jx.T, levels=tick_l, norm=mcolors.SymLogNorm(linthresh=10**(-0.5), linscale=0.2, vmin=-10**1.5,  vmax=10**1.5), cmap='RdBu')
              #plt.contourf(X, Y, jx.T, norm=mcolors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=-30.0,  vmax=30.0), cmap='RdBu_r')
              #### manifesting colorbar, changing label and axis properties ####
              cbar=plt.colorbar(pad=0.01,ticks=tickld)
              cbar.set_label('$j_{x}$ [$I_A/\lambda^2$]',fontdict=font)
              cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=15)        
              plt.text(35,8,'uniform cut_y=0',fontsize=24)
              plt.xlabel('X [$\mu m$]',fontdict=font)
              plt.ylabel('Z [$\mu m$]',fontdict=font)
              plt.xticks(fontsize=20); plt.yticks(fontsize=20);


              plt.subplots_adjust(left=0.05, bottom=0.06, right=0.98, top=0.98,
                wspace=0.091, hspace=0.151)

              fig = plt.gcf()
              fig.set_size_inches(24, 14)
              fig.savefig('./J_cut=0_bilog'+str(n).zfill(4)+'.png',format='png',dpi=100)
              plt.close("all")
              print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
              return 0


if __name__ == '__main__':
  start   =  12 # start time
  stop    =  12  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=1)
  results = pool.map(processplot,inputs)
  print(results)

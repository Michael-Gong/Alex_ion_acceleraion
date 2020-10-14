#%matplotlib inline
import sdf
import matplotlib
import matplotlib as mpl
#mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
import sys
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_jet = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.viridis(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_viridis = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = matplotlib.cm.rainbow(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolor_rainbow = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

def pxpy_to_energy(gamma, weight):
    binsize = 200
    en_grid = np.linspace(50,19950,200)
    en_bin  = np.linspace(0,20000.0,201)
    en_value = np.zeros_like(en_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])
    return (en_grid, en_value)


def theta_to_grid(theta, weight):
    binsize = 240
    theta_grid = np.linspace(-119.5,119.5,240)
    theta_bin  = np.linspace(-120,120,241)
    theta_value = np.zeros_like(theta_grid) 
    for i in range(binsize):
#        if i == binsize-1:
#            en_value[i] = sum(weight[en_bin[i]<=gamma])
#        else:
            theta_value[i] = sum(weight[ (theta_bin[i]<=theta) & (theta<theta_bin[i+1]) ])
    return (theta_grid, theta_value)


if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
  ######## Constant defined here ########
  pi        =     3.1415926535897932384626
  q0        =     1.602176565e-19 # C
  m0        =     9.10938291e-31  # kg
  v0        =     2.99792458e8    # m/s^2
  kb        =     1.3806488e-23   # J/K
  mu0       =     4.0e-7*np.pi       # N/A^2
  epsilon0  =     8.8541878176203899e-12 # F/m
  h_planck  =     6.62606957e-34  # J s
  wavelength=     1.06e-6
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
  font_size = 20
  font_size_2 = 16
  from_path = './cannon_a190_v484/'
  n0=30.0
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
#  weight_1 = V*denunit*n0/50.0 
  weight_1 = V*denunit*n0/50.0 
 # nsteps      = int(sum(1 for line in open(from_path+'t_tot_s.txt'))/part_number)
  nphi  = 180
  ntheta= 120
  to_path   = './cannon_a190_v484_fig/'
  for n in range(15,33,12):
      data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time1=header['time']
      px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      py = py[px>0]
      pz = pz[px>0]
      px = px[px>0]
      gg = (px**2+py**2+pz**2+1)**0.5
      theta = np.arctan((py**2+pz**2)**0.5,px)*180.0/np.pi
      ek_1  = (gg-1.0)*0.51*1836.0
      ww = np.zeros_like(ek_1) + weight_1
      phi   = np.arctan2(pz,py)*180.0/np.pi
   
      for thresh_en in [160,200,240,280,320]:
          px_con = px[(ek_1>(thresh_en-20)) & (ek_1<(thresh_en+20))]  
          py_con = py[(ek_1>(thresh_en-20)) & (ek_1<(thresh_en+20))]  
          pz_con = pz[(ek_1>(thresh_en-20)) & (ek_1<(thresh_en+20))]  
          ww_con = np.zeros_like(px_con) + weight_1
          px_th  = ((thresh_en/918.0+1.0)**2-1)**0.5
          py_con = py_con/px_th
          pz_con = pz_con/px_th
          r1=np.tan(10.0/180*np.pi)
          r2=np.tan(20.0/180*np.pi)
          py_edges = np.linspace(-r2*1.05,r2*1.05,320)
          pz_edges = np.linspace(-r2*1.05,r2*1.05,319)
          py_edges_1 = np.linspace(-r2*1.05,r2*1.05,319)
          pz_edges_1 = np.linspace(-r2*1.05,r2*1.05,318)
          print(py_con,py_edges)
          H, _, _   = np.histogram2d(pz_con, py_con, [pz_edges, py_edges], weights=ww_con)
          PY, PZ = np.meshgrid(py_edges_1,pz_edges_1)
          print()
          fig = plt.figure()
          ax  = fig.add_subplot(111)
          ax.set_facecolor('black')
          levels = np.logspace(1e5,1e7, 41)
          print(H)
          print('Max H:',np.max(H))
          print('Min H:',np.min(H[H>0]))
#          H[H<1e4] = np.nan
          print('Max H:',np.max(H))
          img=ax.pcolormesh(PY,  PZ,  H, norm=colors.LogNorm(vmin=1e5, vmax=1e7), cmap='magma')
          #cax = fig.add_axes([0.65,0.94,0.25,0.02])
#          cbar=fig.colorbar(img,cax=cax, ticks=[1e5,1e6,1e7],orientation='horizontal')
#          cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), fontsize=font_size)
#          cbar.set_label(r'dN/d$p_x$d$p_y$ [A.U.]',fontdict=font)

          x = np.linspace(-r1,r1,320)
          y = (r1**2-x**2)**0.5
          ax.plot(x, y,'--',color='lime',linewidth=2.5)
          ax.plot(x, -y,'--',color='lime',linewidth=2.5)
          x = np.linspace(-r2,r2,320)
          y = (r2**2-x**2)**0.5
          ax.plot(x, y,'--y',color='lime',linewidth=2.5)
          ax.plot(x, -y,'--y',color='lime',linewidth=2.5)

          ax.set_xlabel(r'$p_y/p_x$',fontdict=font)
          ax.set_ylabel(r'$p_z/p_x$',fontdict=font)
  #ax.set_yscale('log')
          ax.set_ylim(-1.05*r2,1.05*r2)
          ax.set_xlim(-1.05*r2,1.05*r2)
          ax.set_xticks([-0.2,0,0.2])
          ax.set_yticks([-0.2,0,0.2])
#          ax.set_xticklabels([-0.2,0,0.2])
#          ax.set_yticklabels([-0.2,0.,0.2])

  #ax.set_theta_zero_location('N')
    #  ax.set_ylabel(r'$\theta\ [^o]$',fontdict=font)
          ax.tick_params(axis='x',labelsize=font_size) 
          ax.tick_params(axis='y',labelsize=font_size)
  #ax.set_title('proton_angular_time='+str(time1), va='bottom', y=1., fontsize=20)
    #  plt.text(-100,650,' t = '++' fs',fontdict=font)

#plt.pcolormesh(x, y, ex.T, norm=mpl.colors.Normalize(vmin=0,vmax=100,clip=True), cmap=cm.cubehelix_r)
#  plt.axis([x.min(), x.max(), y.min(), y.max()])
#### manifesting colorbar, changing label and axis properties ####
#  cbar=plt.colorbar(pad=0.01)#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=font_size)
#  cbar.set_label('dN/dE [A.U.]',fontdict=font)
#  a0=200.0
#  alpha=np.linspace(-3.5,0.5,501)
#  plt.xlabel(r'$\theta$'+' [degree]',fontdict=font)
#  plt.ylabel('time [fs]',fontdict=font)
#  plt.xticks([-135,-90,-45,0,45,90,135],fontsize=font_size); plt.yticks([0,500,1000,1500],fontsize=font_size);
#  plt.title(r'$dN/d\theta$'+' for no RR', fontsize=font_size)
#  plt.xlim(-120,120)
#  plt.ylim(0,1650)
#plt.title('electron at y='+str(round(y[n,0]/2/np.pi,4)),fontdict=font)

          plt.subplots_adjust(top=0.96, bottom=0.01, left=0.01, right=0.96, hspace=0.10, wspace=0.05)


          fig = plt.gcf()
          fig.set_size_inches(5.5, 5.5)
#fig.set_size_inches(5, 4.5)
          fig.savefig(to_path+'px_py'+str(n).zfill(4)+'enth_'+str(thresh_en).zfill(4)+'_'+'.png',format='png',dpi=160)
          plt.close("all")
          print(to_path+'px_py'+str(n).zfill(4)+'enth_'+str(thresh_en).zfill(4)+'_'+'.png')


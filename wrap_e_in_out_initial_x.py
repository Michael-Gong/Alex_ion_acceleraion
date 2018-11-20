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

upper = matplotlib.cm.viridis(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolormap = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

if __name__ == '__main__':
  from_path='./cannon_a190_e_div_2/2_'
  to_path  =from_path
  
  x_in = np.loadtxt(from_path+'xx3d_in.txt')
  y_in = np.loadtxt(from_path+'yy3d_in.txt')
  z_in = np.loadtxt(from_path+'zz3d_in.txt')
  w_in = np.loadtxt(from_path+'ww3d_in.txt')

  x_out= np.loadtxt(from_path+'xx3d_out.txt')
  y_out= np.loadtxt(from_path+'yy3d_out.txt')
  z_out= np.loadtxt(from_path+'zz3d_out.txt')
  w_out= np.loadtxt(from_path+'ww3d_out.txt')

  xx   = np.zeros([np.size(x_in[:,0])+np.size(x_out[:,0]),2])
  yy   = np.zeros_like(xx)
  zz   = np.zeros_like(xx)
  weight = np.zeros_like(xx)
  xx[:np.size(x_in[:,0]),:]=x_in
  xx[np.size(x_in[:,0]):,:]=x_out
  yy[:np.size(y_in[:,0]),:]=y_in
  yy[np.size(y_in[:,0]):,:]=y_out
  zz[:np.size(z_in[:,0]),:]=z_in
  zz[np.size(z_in[:,0]):,:]=z_out
  weight[:np.size(z_in[:,0]),:]=w_in
  weight[np.size(z_in[:,0]):,:]=w_out

  fig,host = plt.subplots()
  #plt.hist2d(xx[:,0], (yy[:,0]**2+zz[:,0]**2)**0.5, bins=(100, 50), range=[[0,15],[0,8]], cmap=mycolormap, weights=weight[:,0], normed=False, norm=colors.Normalize(vmin=0, vmax=10))
  plt.hist2d(xx[:,0], (yy[:,0]**2+zz[:,0]**2)**0.5, bins=(100, 50), range=[[0,15],[0,8]], cmap=mycolormap, weights=weight[:,0]/1.5e5, normed=False,norm=colors.Normalize(vmin=0, vmax=100))
  cbar=plt.colorbar(pad=0.01)
  cbar.set_label(r'$dN/(dxdr)$'+' [A.U.]',fontdict=font)
  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)

  x=np.linspace(0,15,101)
  y=np.zeros_like(x)+1.8
  plt.plot(x,y,':r',linewidth=4.5)
    #plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
    #plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
    #plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
     #   plt.legend(loc='upper right')
#  plt.Circle((0,0),1.8, color='k',fill=False,linestyle=':',linewidth=4,edgecolor='black')
  plt.xlim(0, 15)
  plt.ylim(0, 8)
  plt.xlabel('X [$\mu m$]',fontdict=font)
  plt.ylabel('r [$\mu m$]',fontdict=font)
  plt.xticks([0,5,10,15],fontsize=25); plt.yticks([0,2,4,6,8],fontsize=25);
    #  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.12, bottom=0.15, right=0.99, top=0.96,
                    wspace=None, hspace=None)
  #plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
  print('hehe')
    #plt.show()
    #lt.figure(figsize=(100,100))
  to_path= './cannon_a190_e_div_2/' 
    
  fig = plt.gcf()
  fig.set_size_inches(9.5, 5.5)
  fig.savefig(to_path+'wrap_e_in_out_initial_x.png',format='png',dpi=160)
  plt.close("all")
  print('wrap_e_in_out_initial.png')

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
        'size'   : 24,  
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


def processplot(n): 
  
  from_path = './cannon_a190_v484/'
  to_path   = './cannon_a190_v484_fig/'
  
  n0=30.0
  R=1.8e-6
  L=15e-6
  Ntot = np.pi*R*R*L*n0*denunit
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  weight = V*denunit*n0/50.0 
  
  data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time1=header['time']
  px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  gg = (px**2+py**2+pz**2+1)**0.5
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi

  grid_x = data['Grid/Particles/subset_Only_Ions0/Ion'].data[0]/1.0e-6      
  temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data

  px = px[theta < 10]
  grid_x = grid_x[theta < 10]
  theta = theta[theta < 10]

  if np.size(px) == 0:
    return 0;
  theta[theta < -10] = -10
  theta[theta >  10] =  10

  color_index = abs(theta)

  fig,host = plt.subplots()
  print('hehe1')
  weight = np.zeros_like(grid_x)+weight
  print(max(weight),min(weight))
#    plt.subplot()
#  plt.scatter(grid_x, px, c=color_index, s=0.03, cmap='rainbow_r', edgecolors='None', alpha=0.66)
  #plt.hist2d(grid_x, px, bins=(50, 50), range=[[10,50],[0,1.5]], cmap=plt.cm.jet, weights=weight, norm=mcolors.Normalize(vmin=-2.0, vmax=2.0))
#  plt.hist2d(grid_x, px, bins=(100, 100), range=[[15,45],[0,1.5]], cmap='cubehelix_r', weights=weight, normed=False, norm=colors.Normalize(vmin=0, vmax=1e9))
  plt.hist2d(grid_x, px, bins=(100, 100), range=[[10,45],[0,1.5]], cmap=cmap_wyr, weights=weight, normed=False, norm=colors.LogNorm(vmin=5e7, vmax=1e9))
#  cbar=plt.colorbar(pad=0.3)
#  cbar.set_label(r'$dN/(dxdp_x)$'+' [A.U.]',fontdict=font)
#  cbar.ax.tick_params(labelsize=20)
#  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)
#plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
#plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
#plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
#plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
 #   plt.legend(loc='upper right')
  plt.xlim(10,43)
  plt.ylim(0.,1.2)
  plt.xlabel('X [$\mu m$]',fontdict=font)
  plt.ylabel('$p_x$ [m$_p$c$^2$]',fontdict=font)
  plt.xticks([15,25,35],fontsize=25); 
  plt.yticks([0,0.5,1.0],fontsize=25);
#  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.15, bottom=0.15, right=0.76, top=None,
                wspace=None, hspace=None)
  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
  print('hehe')
#plt.show()
#lt.figure(figsize=(100,100))

  par1 = host.twinx()
  par2 = host.twinx()
  #par3 = host.twinx()

  par2.spines["right"].set_position(("axes", 1.2))
  make_patch_spines_invisible(par2)
  par2.spines["right"].set_visible(True)

  #par3.spines["right"].set_position(("axes", 1.1))
  #make_patch_spines_invisible(par3)
  #par3.spines["right"].set_visible(True)

  tkw = dict(size=25, width=2.)

  x  = np.loadtxt(from_path+'ex_lineout_x.txt')
  ex = np.loadtxt(from_path+'ex_lineout_r15_'+str(n).zfill(4)+'.txt')*3.5
  #p1, = par1.plot(x,ex, "-",color='lime', label="$E_x$",linewidth=3)

  x  = np.loadtxt(from_path+'iden_lineout_x.txt')
  iden = np.loadtxt(from_path+'iden_lineout_r15_'+str(n).zfill(4)+'.txt')#*exunit/denunit
  x  = np.loadtxt(from_path+'eden_lineout_x.txt')
  eden = np.loadtxt(from_path+'eden_lineout_r15_'+str(n).zfill(4)+'.txt')#*exunit/denunit
  p2, = par2.plot(x,eden-iden, "-b", label="$n_e$",linewidth=2.5,zorder=2)
  p1, = par1.plot(x,ex, "-",color='lime', label="$E_x$",linewidth=3,zorder=3)
  

#  p1, = par1.plot(x,iden, "-r", label="$Z_in_i$",linewidth=3)
#  p3, = par3.plot(x,iden, "-r", label="Ion")
#  p3, = par3.plot(x,iden, "-y", label="Ion")
#  p3, = par3.plot(x,cden*6, "-g", label="Ion")
#  par3.set_ylabel('$n^+\ [n_c]$')

#  par1.legend(loc='upper center',fontsize=20,framealpha=1.0)
  par1.set_ylim(0,28)
  par1.set_ylabel('$E_x$ [TV/m]',fontdict=font,color='lime')
  par1.tick_params(axis='y',labelsize=22,colors='lime')

#  par2.legend(loc='upper center',fontsize=20,framealpha=1.0)
  par2.set_ylim(0,2.2)
#  par2.set_ylabel('$Z_i\overline{n}_i-\overline{n}_e$ [$n_c$]',fontdict=font,color='b')
  par2.set_ylabel('$\overline{n}_e-Z_i\overline{n}_i$ [$n_c$]',fontdict=font,color='b')
  par2.tick_params(axis='y',labelsize=22,colors='b')


  fig = plt.gcf()
  fig.set_size_inches(8.5, 6)
  fig.savefig(to_path+'wrap_px_x_'+str(n).zfill(4)+'.png',format='png',dpi=160)
  plt.close("all")
  print('finised '+str(n).zfill(4))
  return 0

if __name__ == '__main__':
  start   =  3  # start time
  stop    =  51  # end time
  step    =  2  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=1)
  results = pool.map(processplot,inputs)
  print(results)


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
  
  to_path='./cannon_a190_bulk200/'
  from_path = './cannon_a190_bulk200/'
  
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

  #choice = np.random.choice(range(part13_id.size), 20000, replace=False)
  #part13_id = part13_id[choice]
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
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi

  grid_x = data['Grid/Particles/subset_Only_Ions0/Ion'].data[0]/1.0e-6      
  temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data

  px = px[theta < 20]
  grid_x = grid_x[theta < 20]
  theta = theta[theta < 20]

  if np.size(px) == 0:
    return 0;
  theta[theta < -20] = -20
  theta[theta >  20] =  20

  color_index = abs(theta)

  fig,host = plt.subplots()
#    plt.subplot()
  plt.scatter(grid_x, px, c=color_index, s=0.03, cmap='rainbow_r', edgecolors='None', alpha=0.66)
  cbar=plt.colorbar( ticks=np.linspace(np.min(color_index), np.max(color_index), 5) ,pad=0.01)
  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
  cbar.set_label(r'$|\theta|$'+' [degree]',fontdict=font)

#plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
#plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
#plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
#plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
 #   plt.legend(loc='upper right')
  plt.xlim(-5,55)
  plt.ylim(0.,1.6)
  plt.xlabel('X [$\mu m$]',fontdict=font)
  plt.ylabel('$p_x$ [m$_i$c$^2$]',fontdict=font)
  plt.xticks(fontsize=20); plt.yticks(fontsize=20);
#  plt.text(-100,650,' t = '++' fs',fontdict=font)
  plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                wspace=None, hspace=None)
  plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
#plt.show()
#lt.figure(figsize=(100,100))

  par1 = host.twinx()
  par2 = host.twinx()
  par3 = host.twinx()

  par2.spines["right"].set_position(("axes", 1.05))
  make_patch_spines_invisible(par2)
  par2.spines["right"].set_visible(True)

  par3.spines["right"].set_position(("axes", 1.1))
  make_patch_spines_invisible(par3)
  par3.spines["right"].set_visible(True)

  tkw = dict(size=20, width=1.)

  x  = np.loadtxt(from_path+'ex_lineout_x.txt')
  ex = np.loadtxt(from_path+'ex_lineout_r15_'+str(n).zfill(4)+'.txt')
  p1, = par1.plot(x,ex, "-k", label="Ex")
  par1.set_ylabel(r'$E_x\ [m_ec\omega/|e|]$')
  par1.yaxis.label.set_color(p1.get_color())
  par1.tick_params(axis='y', colors=p1.get_color(), **tkw)
  par1.set_ylim(-10,15)

  x  = np.loadtxt(from_path+'eden_lineout_x.txt')
  eden = np.loadtxt(from_path+'eden_lineout_r15_'+str(n).zfill(4)+'.txt')#*exunit/denunit
  p2, = par2.plot(x,eden, "-b", label="Electron")
  par2.set_ylabel('$n_e\ [n_c]$')
  par2.yaxis.label.set_color(p2.get_color())
  par2.tick_params(axis='y', colors=p2.get_color(), **tkw)
  par2.set_ylim(0,30)
  

  x  = np.loadtxt(from_path+'iden_lineout_x.txt')
  iden = np.loadtxt(from_path+'iden_lineout_r15_'+str(n).zfill(4)+'.txt')#*exunit/denunit
  p3, = par3.plot(x,iden, "-r", label="Ion")
  par3.set_ylabel('$n_i\ [n_c]$')
  par3.yaxis.label.set_color(p3.get_color())
  par3.tick_params(axis='y', colors=p3.get_color(), **tkw)
  par3.set_ylim(0,30)

  fig = plt.gcf()
  fig.set_size_inches(12, 7.5)
  fig.savefig(to_path+'r15_comb_proton'+str(n).zfill(4)+'.png',format='png',dpi=80)
  plt.close("all")
  print('finised '+str(n).zfill(4))
  return 0

if __name__ == '__main__':
  start   =  3  # start time
  stop    =  31  # end time
  step    =  2  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=5)
  results = pool.map(processplot,inputs)
  print(results)

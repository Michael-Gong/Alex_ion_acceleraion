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


if __name__ == '__main__':

  n0=30.0
  R=1.8e-6
  L=15e-6
  Ntot = np.pi*R*R*L*n0*denunit
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  
  start   =  3 # start time
  stop    =  51  # end time
  step    =  2  # the interval or step

  for dir_n in ['t20','t40','t50','t60','t80','t100','t120']:
    to_path   = './cannon_a190_'+dir_n+'_fig/'
    from_path = './cannon_a190_'+dir_n+'/'  
    if not os.path.exists(to_path):
          os.mkdir(to_path)
    n0=30.0
#    n0=float(dir_n[-2:])
    R=1.8e-6
    L=15e-6
    Ntot = np.pi*R*R*L*n0*denunit
    V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  ######### Script code drawing figure ################
    fig,host = plt.subplots()
  
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
  
    for n in range(start,stop+step,step):
      if not os.path.exists(from_path+'ex_lineout_r15_'+str(n).zfill(4)+'.txt'):
          continue
      x  = np.loadtxt(from_path+'ex_lineout_x.txt')
      ex = np.loadtxt(from_path+'ex_lineout_r15_'+str(n).zfill(4)+'.txt')
      p1, = par1.plot(x,ex, "-",color='lime', label="$E_x$",linewidth=1)
    
      x  = np.loadtxt(from_path+'iden_lineout_x.txt')
      iden = np.loadtxt(from_path+'iden_lineout_r15_'+str(n).zfill(4)+'.txt')#*exunit/denunit
      x  = np.loadtxt(from_path+'eden_lineout_x.txt')
      eden = np.loadtxt(from_path+'eden_lineout_r15_'+str(n).zfill(4)+'.txt')#*exunit/denunit
      plt.plot(x,eden-iden, "-b", label="$n_e-Z_in_i$",linewidth=1)
    #  p3, = par3.plot(x,iden, "-r", label="Ion")
    #  p3, = par3.plot(x,iden, "-y", label="Ion")
    #  p3, = par3.plot(x,cden*6, "-g", label="Ion")
    #  par3.set_ylabel('$n^+\ [n_c]$')
    
  #  par1.legend(loc='upper center',fontsize=20,framealpha=1.0)
#    plt.xlim(0,55)
#    plt.ylim(0,15)
#    plt.xlabel('X [$\mu m$]',fontdict=font)
#    plt.ylabel('$\overline{n}_e-Z_i\overline{n}_i$ [$n_c$]',fontdict=font)

#    plt.xticks(fontsize=25); 
#    plt.yticks(fontsize=25);

    par1.set_xlim(0,50)
    par1.set_xlabel('$x$ [$\mu m$]',fontdict=font,color='k')
    par1.set_ylim(0,15)
    par1.set_ylabel('$E_x$ [$m_ec\omega/|e|$]',fontdict=font,color='lime')
    par1.tick_params(axis='y',labelsize=25,colors='lime')
  
    par2.set_ylim(0,15)
#    par2.legend(loc='upper center',fontsize=20,framealpha=1.0)
    par2.set_ylabel('$Z_i\overline{n}_i-\overline{n}_e$ [$n_c$]',fontdict=font,color='b')
  
  
    fig = plt.gcf()
    fig.set_size_inches(8.5, 6)
    fig.savefig(to_path+'system_ex_x_'+dir_n+'.png',format='png',dpi=160)
    plt.close("all")
    print(to_path+'system_ex_x_'+dir_n+'.png')

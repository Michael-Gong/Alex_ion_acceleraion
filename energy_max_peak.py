#import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
  
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
          'size'   : 23,  
          }  
  
##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####



  to_path='./'

  n0=30.0
  R=1.8e-6
  L=15e-6
  Ntot = np.pi*R*R*L*n0*denunit
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  weight = V*denunit*n0/50.0 

  density = np.array([15,  20,  30,  40,  50,  60])
  en_max  = np.array([435, 500, 550, 605, 630, 600])
  en_peak = np.array([245, 225, 235, 245, 205, 185])


# plt.plot(density,upper,'--k',linewidth=4, label='Eqs.(4,5)',zorder=0)
# plt.plot(power,lower,'--g',linewidth=4, label='With depletion',zorder=1)
  plt.scatter(density, en_max, c='deepskyblue',marker='o',s=200, label='3D PIC (cut-off energy)', edgecolors='black', linewidth='3',alpha=1,zorder=2)
  plt.scatter(density, en_peak,c='tomato',marker='^',s=200, label='3D PIC (peak energy)', edgecolors='black',linewidth='3',alpha=1,zorder=3)


  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('$n_e$ [$n_c$]',fontdict=font)
  plt.ylabel('E$_k$ [MeV]',fontdict=font)
  plt.xticks(fontsize=23); plt.yticks([0,250,500,750,1000],fontsize=23);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(10,70) 
  plt.ylim(0,1000) 
  plt.legend(loc='upper right',fontsize=18,framealpha=0.5)
  plt.subplots_adjust(left=0.16, bottom=0.15, right=0.95, top=0.95,
          wspace=None, hspace=None)
  #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(8.4, 7.0)
  fig.savefig('./en_max_peak.png',format='png',dpi=160)
  plt.close("all")

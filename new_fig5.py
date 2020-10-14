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

  ######### Parameter you should set ###########
  start   =  12  # start time
  stop    =  30  # end time
  step    =  1  # the interval or step
  
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  youwant =  ['electron_en']#,'electron_no_en']#,'ey','ex','ey_averaged','bz','bz_averaged','Subset_high_e_density','Subset_high_e_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  #if (os.path.isdir('jpg') == False):
  #  os.mkdir('jpg')
  ######### Script code drawing figure ################
  n0=30.0
  R=1.8e-6
  L=15e-6
  Ntot = np.pi*R*R*L*n0*denunit
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  weight = V*denunit*n0/50.0 

  set_relativistic =0 
  for i in [1]:
      dataf = np.loadtxt('./fig5_data.txt')
      power =  dataf[:,0]
      upper =  dataf[:,1]
      lower =  dataf[:,2]
     
      power_c = np.array([0.5,1.0,2.0,4.0])
      result_c = np.array([0.27,0.5,0.55,0.65])
      power_p = np.array([0.5,1.0,2.0,4.0])
      result_p = np.array([0.149,0.183,0.221,0.230])

      plt.plot(power,upper,'--k',linewidth=4, label='Eqs.(4,5)',zorder=0)
      plt.plot(power,lower,'--g',linewidth=4, label='With depletion',zorder=1)
      plt.scatter(power_c,result_c,c='deepskyblue',marker='o',s=200, label='3D PIC (cut-off energy)', edgecolors='black', linewidth='3',alpha=1,zorder=2)
      plt.scatter(power_p,result_p,c='tomato',marker='^',s=200, label='3D PIC (peak energy)', edgecolors='black',linewidth='3',alpha=1,zorder=3)


      #### manifesting colorbar, changing label and axis properties ####
      plt.xlabel('P [PW]',fontdict=font)
      plt.ylabel('E$_k$ [GeV]',fontdict=font)
      plt.xticks(fontsize=23); plt.yticks([0.0,0.25,0.5,0.75,1.0,1.25],fontsize=23);
      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
      plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
      plt.xlim(0,5) 
      plt.ylim(0,1.25) 
      plt.legend(loc='upper right',fontsize=18,framealpha=0.5)
      plt.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95,
                wspace=None, hspace=None)
    #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
      fig = plt.gcf()
      fig.set_size_inches(8.4, 7.0)
      fig.savefig('./new_fig5.png',format='png',dpi=160)
      plt.close("all")

#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
  
if __name__ == "__main__":
  print ('This is main of module "*.py"')
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


  from_path = './'
  to_path='./'

  ######### Parameter you should set ###########
  start   =  1  # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
  
  #  os.mkdir('jpg')
  ######### Script code drawing figure ################
#  dir = ['cannon_a190_n15', 'cannon_a190_n20', 'cannon_a190', 'cannon_a190_n40', 'cannon_a190_n50', 'cannon_a190_n60'] 
  dir = ['uniform_a190_n15', 'uniform_a190_n20', 'uniform_a190_n30', 'uniform_a190_n40', 'uniform_a190_n50', 'uniform_a190_n60'] 
 
  label = ['n15', 'n20', 'n30', 'n40', 'n50', 'n60'] 

  for i in range(6):
      data_y = np.loadtxt(from_path+dir[i]+'/e_number_cylinder_r20.txt')*denunit
      data_x = np.loadtxt(from_path+dir[i]+'/e_time_cylinder.txt')
     

      plt.plot(data_x,data_y,linewidth=4, label=label[i])
      print('finished in '+dir[i]) 
  #### manifesting colorbar, changing label and axis properties ####
  plt.xlabel('time [fs]',fontdict=font)
  plt.ylabel('Number',fontdict=font)
  plt.xticks(fontsize=20);
  plt.yticks(fontsize=20);
#  plt.yticks([0.0,0.5,1.0,1.5,2.0],fontsize=20);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)

  plt.legend(loc='best',fontsize=18,framealpha=1.0)
  plt.subplots_adjust(left=0.16, bottom=0.15, right=0.99, top=0.95,
            wspace=None, hspace=None)
#        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(10.0, 8.5)
  fig.savefig(to_path+'./time_electron_num_uniform_r20.png',format='png',dpi=160)
  plt.close("all")

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


  def pxpy_to_energy(gamma, weight):
      binsize = 300
      en_grid = np.linspace(0.5,599.5,300)
      en_bin  = np.linspace(0,600.0,301)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(600.0/binsize)
      return (en_grid, en_value)

  to_path='./'

  ######### Parameter you should set ###########
  start   =  1  # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
 
  x  = np.loadtxt('ex_lineout_x.txt')

  for i in range(start,stop+step,step):
      ex = np.loadtxt('ex_lineout_'+str(i).zfill(4)+'.txt')
      
      plt.plot(x,ex,'-b',linewidth=4)
    
      #### manifesting colorbar, changing label and axis properties ####
      plt.xlabel('X [$\mu m$]',fontdict=font)
      plt.ylabel('$E_x$ [$m_ec\omega/|e|$]',fontdict=font)
      plt.xticks(fontsize=20); plt.yticks(fontsize=20);
      #plt.yscale('log')
      plt.ylim(-40,40)
      plt.xlim(-5,55)
      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
      plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)

      plt.legend(loc='best',fontsize=20,framealpha=1.0)
      plt.subplots_adjust(left=None, bottom=0.15, right=0.95, top=0.95,
                wspace=None, hspace=None)
    #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
      fig = plt.gcf()
      fig.set_size_inches(10.0, 6.0)
      fig.savefig('./ex_line_'+str(i).zfill(4)+'.png',format='png',dpi=80)
      plt.close("all")
      print('finished!'+str(i).zfill(4))

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
import matplotlib.colors as mcolors
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec

import multiprocessing as mp

  
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
      binsize = 500
      en_grid = np.linspace(1,999,500)
      en_bin  = np.linspace(0,1000.0,501)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(1000.0/binsize)
      return (en_grid, en_value)

def processplot(n):
    to_path='./uniform_a190_n60/'
    from_path = './uniform_a190_n60/'
    ######### Parameter you should set ###########
    n0=30.0
    R=1.8e-6
    L=15e-6
    Ntot = np.pi*R*R*L*n0*denunit
    V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
    weight = V*denunit*n0/20.0 
  #  weight = Ntot/(1200*360*360*50)
  
    set_relativistic =1 
    if 1 > 0:
     
        data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
        header=data['Header']
        time1=header['time']
        if set_relativistic == 1:
            px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
            py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
            pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
            gg = (px**2+py**2+pz**2+1)**0.5
            theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
      
            ek_1 = (gg - 1.0)*0.51*1836
            ww_1 = np.zeros_like(ek_1) + weight
            dist_x1, den1 = pxpy_to_energy(ek_1,ww_1)
      
            ek_2 = (gg[abs(theta)<10.0] - 1.0)*0.51*1836
            ww_2 = np.zeros_like(ek_2) + weight
            dist_x2, den2 = pxpy_to_energy(ek_2,ww_2)
            print('set_relativistic 1'+str(n).zfill(4))
        elif set_relativistic == 0: 
            px = data['Particles/Px/subset_Only_Ions0/Ion'].data
            py = data['Particles/Py/subset_Only_Ions0/Ion'].data
            pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data
            pp = (px**2+py**2+pz**2)**0.5
            theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
      
            ek_1 = pp**2/2.0/1836/m0/(1.6e-13)
            ww_1 = np.zeros_like(ek_1) + weight
            dist_x1, den1 = pxpy_to_energy(ek_1,ww_1)
      
            ek_2 = (pp[abs(theta)<10.0])**2/2.0/1836/m0/(1.6e-13)
            ww_2 = np.zeros_like(ek_2) + weight
            dist_x2, den2 = pxpy_to_energy(ek_2,ww_2)
            print('set_relativistic 0'+str(n).zfill(4))
            
      
        plt.plot(dist_x1,den1,':b',linewidth=4, label=str(round(time1/1e-15,0))+'; total')
        plt.plot(dist_x2,den2,':r',linewidth=4, label=str(round(time1/1e-15,0))+'; '+r'$\theta$'+'< 10$^o$')
      
        #### manifesting colorbar, changing label and axis properties ####
        plt.xlabel('Energy [MeV]',fontdict=font)
        plt.ylabel('dN/dE [per MeV]',fontdict=font)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);
        plt.yscale('log')
        #plt.ylim(2e7,8e9)
        plt.xlim(5,1000)
        plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
        plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  
        plt.legend(loc='best',fontsize=20,framealpha=1.0)
        plt.subplots_adjust(left=None, bottom=0.15, right=0.95, top=0.95,
                  wspace=None, hspace=None)
      #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
        fig = plt.gcf()
        fig.set_size_inches(10.0, 6.0)
        fig.savefig(to_path+'proton_spectrum_'+str(n).zfill(4)+'.png',format='png',dpi=80)
        plt.close("all")
        print('finished!'+str(n).zfill(4))      

if __name__ == '__main__':
  start   =  3 # start time
  stop    =  31  # end time
  step    =  2  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=10)
  results = pool.map(processplot,inputs)
  print(results)

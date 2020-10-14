#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from os import path
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
  wavelength=     1.0e-6
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
      binsize = 1000
      en_grid = np.linspace(0.5,999.5,1000)
      en_bin  = np.linspace(0,1000.0,1001)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(1000.0/binsize)
      return (en_grid, en_value)


  dir_list = ['no40','no60','no140','no180','no220']
  dir_n    = [11,     11,    11,     11,     11]
  dir_limit= [[205,245], [207,247], [206,246], [168,208],  [156,196]]


  for i in range(5):
      from_path='./cannon_a190_'+dir_list[i]+'/'
      n = dir_n[i]
#      to_path  ='./cannon_a190_'+dir_n+'_fig/'
#      if not os.path.exists(to_path):
#          os.mkdir(to_path)
      ######### Parameter you should set ###########
      
	  
	  #  os.mkdir('jpg')
      ######### Script code drawing figure ################
      n0=30.0
#      n0=float(dir_list[i][-2:]) #30.0
      R=1.8e-6
      L=15e-6
      Ntot = np.pi*R*R*L*n0*denunit
      V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
      weight = V*denunit*n0/50.0 
    #  weight = Ntot/(1200*360*360*50)
    
      set_relativistic =1
      if 1>0: 
          if not os.path.exists(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf"):
              continue
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
              condition = (ek_1>dir_limit[i][0]) & (ek_1<dir_limit[i][1])
              charge_1 = np.sum(ww_1[condition])*1.6e-19              

              ek_2 = (gg[abs(theta)<10.0] - 1.0)*0.51*1836
              ww_2 = np.zeros_like(ek_2) + weight
              dist_x2, den2 = pxpy_to_energy(ek_2,ww_2)
              condition = (ek_2>dir_limit[i][0]) & (ek_2<dir_limit[i][1])
              charge_2 = np.sum(ww_2[condition])*1.6e-19              
       
              print(from_path+str(n).zfill(4)+'<10 degree:',charge_2, '  all degree:',charge_1, ' range',dir_limit[i][0],dir_limit[i][1]) 

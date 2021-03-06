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


  def pxpy_to_energy(gamma, weight):
      binsize = 500
      en_grid = np.linspace(0.5,499.5,binsize)
      en_bin  = np.linspace(0,500.0,binsize+1)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(500.0/binsize)
      return (en_grid, en_value)

  to_path='./cannon_a190_v484_fig/'
  from_path = './cannon_a190_v484/' 

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
  weight2 = V*denunit*n0/10.0 

  set_relativistic =1 
  for n in np.arange(15,16):
    for m in np.arange(27,28): 
      data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time=header['time']
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
           
      data = sdf.read(from_path+"/i_tot_loc"+str(m).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time=header['time']
      if set_relativistic == 1:
          px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
          py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
          pz = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
          gg = (px**2+py**2+pz**2+1)**0.5
          theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
    
          ek_3 = (gg - 1.0)*0.51*1836
          ww_3 = np.zeros_like(ek_3) + weight
          dist_x3, den3 = pxpy_to_energy(ek_3,ww_3)
    
          ek_4 = (gg[abs(theta)<10.0] - 1.0)*0.51*1836
          ww_4 = np.zeros_like(ek_4) + weight
          dist_x4, den4 = pxpy_to_energy(ek_4,ww_4)
          print('set_relativistic 1'+str(m).zfill(4))

          
#      plt.plot(dist_x1,den1,':b',linewidth=4, label='before; total')
#      plt.plot(dist_x2,den2,':r',linewidth=4, label='before; '+r'$\theta$'+'< 10$^o$')
#      plt.plot(dist_x3,den3,'-b',linewidth=4, label='after;   total')
#      plt.plot(dist_x4,den4,'-r',linewidth=4, label='after; '+r'  $\theta$'+'< 10$^o$')
   
#      plt.plot(dist_x1,den1/1e9,':b',linewidth=4, label='before; total')
      plt.plot(dist_x2,den2/1e9,':r',linewidth=4, label='before; '+r'$\theta$'+'< 10$^o$')
#      plt.plot(dist_x3,den3/1e9,'-b',linewidth=4, label='after;   total')
      plt.plot(dist_x4,den4/1e9,'-r',linewidth=4, label='after; '+r'  $\theta$'+'< 10$^o$')

      start_index = 220
      stop_index  = 240
      y1 = den4[start_index:stop_index]/1e9
      y2 = np.zeros_like(y1)
      print(y1)
      print(y2)
      print(dist_x4)
      plt.fill_between(dist_x4[start_index:stop_index],y2,y1,color='green',alpha=0.15)


      from_p = '../alex_ion_acc/uniform_a190_n30/'
      data = sdf.read(from_p+"i_tot_loc"+str(27).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time=header['time']
      px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
      gg = (px**2+py**2+pz**2+1)**0.5
      theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
    
      ek_2 = (gg[abs(theta)<10.0] - 1.0)*0.51*1836
      ww_2 = np.zeros_like(ek_2) + weight2
      dist_x2, den2 = pxpy_to_energy(ek_2,ww_2)
      plt.plot(dist_x2,den2/1e9,linestyle='-',color='k',linewidth=4, label='uniform; '+r'$\theta$'+'< 10$^o$')
	  
      #### manifesting colorbar, changing label and axis properties ####
      plt.xlabel(r'$\varepsilon_p$ [MeV]',fontdict=font)
      plt.ylabel(r'dN/d$\varepsilon_p$ [$10^9$MeV$^{-1}$]',fontdict=font)
      plt.xticks([100,200,300],fontsize=23); 
      plt.yticks([0.5,1,1.5,2],fontsize=23);
#      plt.yscale('log')
#      plt.ylim(0.02e9,3e9)
      plt.ylim(0.02,2.45)
      plt.xlim(100,350)
#      plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
#      plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)

      plt.legend(loc='best',fontsize=21,framealpha=0.5)
      plt.subplots_adjust(left=0.14, bottom=0.15, right=0.95, top=0.95,
                wspace=None, hspace=None)
    #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
      fig = plt.gcf()
      fig.set_size_inches(5.5, 7.0)
      fig.savefig(to_path+'relativistic_new'+str(n).zfill(4)+'_'+str(m).zfill(4)+'.png',format='png',dpi=160)
      plt.close("all")
      print(to_path+'relativistic_new'+str(n).zfill(4)+'_'+str(m).zfill(4)) 

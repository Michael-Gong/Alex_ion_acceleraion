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
          'size'   : 25,  
          }  
  font_size =25
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
    

  ##scan for temporal duration
  t_duration = np.array([20,   40,  50,  60,  80, 100, 120])
  peaken_t   = np.array([132,  183, 216, 230, 264, 278, 290])
  charge_t   = np.array([1.65e-9, 3.43e-9, 5.56e-9, 7.82e-9, 8.31e-9, 1.14e-8, 1.24e-8]) #+-20 MeV
#  charge_t   = np.array([8.55e-10, 1.92e-9, 3.59e-9, 5.03e-9, 5.16e-9, 6.65e-9, 7.33e-9]) #+-10 MeV
#  charge_t   = np.array([1.17e-9, 3.10e-9, 6.85e-9, 1.07e-8, 1.42e-8, 2.04e-8, 1.95e-8]) 200~300 MeV
  ##scan for channel density
  n_inner    = np.array([10, 20, 30, 40, 50, 60])
  peaken_n_in= np.array([259, 272, 230, 219, 195, 161])
  charge_n_in= np.array([6.98e-9, 6.95e-9, 7.82e-9, 4.95e-9, 5.66e-9, 6.66e-9])  # +-20 MeV
#  charge_n_in= np.array([9.98e-9, 1.16e-8, 1.07e-8, 6.71e-9, 2.64e-9, 6.72e-10]) 200~300 MeV
  ##scan for wall density
  n_outer    = np.array([40, 60, 100, 140, 180, 220])
  peaken_n_ou= np.array([225, 227, 230, 226, 188, 176])
  charge_n_ou= np.array([1.95e-9, 4.41e-9, 7.82e-9, 4.58e-9, 3.07e-9, 2.13e-9]) #+-20 MeV
#  charge_n_ou= np.array([3.72e-9, 1.36e-8, 1.07e-8, 9.56e-9, 1.61e-9, 4.46e-12]) 200~300 MeV

  a0_        = np.array([60, 85, 110, 150, 190, 240])
  peaken_a0  = np.array([36, 82, 140, 173, 230, 338])
  charge_a0  = np.array([1.59e-9, 3.93e-9, 2.45e-9, 4.71e-9, 7.82e-9, 3.76e-9]) #+-20 MeV
  

  plt.subplot(2,4,1)
#  plt.plot(power,upper,'--k',linewidth=4, label='Eqs.(4,5)',zorder=0)
#      plt.plot(power,lower,'--g',linewidth=4, label='With depletion',zorder=1)
  plt.scatter(t_duration,peaken_t,c='deepskyblue',marker='o',s=200, label='3D PIC (peak energy)', edgecolors='black', linewidth='3',alpha=1,zorder=2)
  t= np.linspace(10,130,300)
#  y0= t**0.5*28
#  plt.plot(t,y0,linestyle='-',linewidth=2.0,zorder=1)
  y1= (1-30.0/(0.5*190))**0.5*0.35*190*(t/3.33333333333)**0.5
  plt.plot(t,y1,c='blue',linestyle='-',linewidth=2.0,zorder=1)
#  y2= t**0.5*30
#  plt.plot(t,y2,linestyle='-',linewidth=2.0,zorder=1)
#  y3= t**0.5*27
#  plt.plot(t,y3,linestyle='-',linewidth=2.0,zorder=1)
  plt.xlabel(r'$\tau$ [fs]',fontdict=font)
  plt.ylabel(r'$\varepsilon_p$ [MeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(0,140) 
  plt.ylim(0,380) 
  plt.legend(loc='best',fontsize=18,framealpha=0.5)


  plt.subplot(2,4,5)
  plt.scatter(t_duration,charge_t/1e-9,c='deepskyblue',marker='o',s=200, label=r'charge with $\theta<10^\circ$ (200~300 MeV)', edgecolors='black',linewidth='3',alpha=1,zorder=3)
  plt.xlabel(r'$\tau$ [fs]',fontdict=font)
  plt.ylabel(r'$\mathcal{Q}_p$ [nC]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(0,140) 
  plt.ylim(0,14) 



  plt.subplot(2,4,2)
#  plt.plot(power,upper,'--k',linewidth=4, label='Eqs.(4,5)',zorder=0)
  plt.scatter(n_inner,peaken_n_in,c='tomato',marker='s',s=200, label='3D PIC (peak energy)', edgecolors='black', linewidth='3',alpha=1,zorder=2)
  n_in = np.linspace(5,65,300)
  y1= (1-n_in/(0.5*190))**0.5*0.35*190*18.0**0.5 #(1-n_in/95.0)**0.5*285  # k*a_0*[1-n_e/(0.5*a_0n_cr)]**0.5*(tau)**0.5; k= 
  plt.plot(n_in,y1,c='red',linestyle='-',linewidth=2.0,zorder=1)
  plt.xlabel(r'$n_{e,in}$ [$n_{cr}$]',fontdict=font)
  plt.ylabel(r'$\varepsilon_p$ [MeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(0,70) 
  plt.ylim(0,400) 
  plt.legend(loc='best',fontsize=18,framealpha=0.5)

  plt.subplot(2,4,6)
  plt.scatter(n_inner,charge_n_in/1e-9,c='tomato',marker='s',s=200, label=r'charge with $\theta<10^\circ$ (200~300 MeV)', edgecolors='black',linewidth='3',alpha=1,zorder=3)
  plt.xlabel(r'$n_{e,in}$ [$n_{cr}$]',fontdict=font)
  plt.ylabel(r'$\mathcal{Q}_p$ [nC]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(0,70) 
  plt.ylim(0,14) 


  plt.subplot(2,4,3)
#  plt.plot(power,upper,'--k',linewidth=4, label='Eqs.(4,5)',zorder=0)
  no_line = np.linspace(25,235,200)
  y1= (1-30.0/(0.5*190))**0.5*0.35*190*18.0**0.5+np.zeros_like(no_line)  # k*a_0*[1-n_e/(0.5*a_0n_cr)]**0.5*(tau)**0.5; k=0.354
  plt.plot(no_line,y1,c='green',linestyle='-',linewidth=2.0,zorder=1)
  plt.scatter(n_outer,peaken_n_ou,c='limegreen',marker='^',s=250, label='3D PIC (peak energy)', edgecolors='black', linewidth='3',alpha=1,zorder=2)
  plt.xlabel(r'$n_{e,out}$ [$n_{cr}$]',fontdict=font)
  plt.ylabel(r'$\varepsilon_p$ [MeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(0,240) 
  plt.ylim(0,400) 
  plt.legend(loc='best',fontsize=18,framealpha=0.5)

  plt.subplot(2,4,7)
  plt.scatter(n_outer,charge_n_ou/1e-9,c='limegreen',marker='^',s=250, label=r'charge with $\theta<10^\circ$ (200~300 MeV)', edgecolors='black',linewidth='3',alpha=1,zorder=3)
  plt.xlabel(r'$n_{e,out}$ [$n_{cr}$]',fontdict=font)
  plt.ylabel(r'$\mathcal{Q}_p$ [nC]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(0,240) 
  plt.ylim(0,14) 


  plt.subplot(2,4,4)
#  plt.plot(power,upper,'--k',linewidth=4, label='Eqs.(4,5)',zorder=0)
  line_a0 = np.linspace(60,255,400)
  y1= (1-30.0/(0.5*line_a0))**0.5*0.35*line_a0*18.0**0.5  # k*a_0*[1-n_e/(0.5*a_0n_cr)]**0.5*(tau)**0.5; k=0.354
  plt.plot(line_a0,y1,c='purple',linestyle='-',linewidth=2.0,zorder=1)
  plt.scatter(a0_,peaken_a0,c='purple',marker='D',s=200, label='3D PIC (peak energy)', edgecolors='black', linewidth='3',alpha=1,zorder=2)
  plt.xlabel(r'$a_0$',fontdict=font)
  plt.ylabel(r'$\varepsilon_p$ [MeV]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(20,280) 
  plt.ylim(0,400) 
  plt.legend(loc='best',fontsize=18,framealpha=0.5)

  plt.subplot(2,4,8)
  plt.scatter(a0_,charge_a0/1e-9,c='purple',marker='D',s=200, label=r'charge with $\theta<10^\circ$ (200~300 MeV)', edgecolors='black',linewidth='3',alpha=1,zorder=3)
  plt.xlabel(r'$a_0$',fontdict=font)
  plt.ylabel(r'$\mathcal{Q}_p$ [nC]',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.grid(which='major',color='k', linestyle='--', linewidth=0.1)
#  plt.grid(which='minor',color='k', linestyle='--', linewidth=0.1)
  plt.xlim(20,280) 
  plt.ylim(0,14) 


  plt.subplots_adjust(left=0.1, bottom=0.1, right=0.98, top=0.98,
                wspace=0.24, hspace=None)
    #        plt.text(250,6e9,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
  fig = plt.gcf()
  fig.set_size_inches(29.33, 12.0)
  fig.savefig('./wrap_peaken_scan_20MeV.png',format='png',dpi=160)
  plt.close("all")

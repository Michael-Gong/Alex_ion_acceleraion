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


def processplot(n): 
  ######### Parameter you should set ###########
  #start   =  23  # start time
  #stop    =  30  # end time
  #step    =  1  # the interval or step
  #youwant =  ['ey','ex','ey_averaged','bz','bz_averaged'] #,'electron_en','electron_ekbar','electron_density']
  #youwant.append('Ion_ekbar')
  #youwant.append('positron_ekbar')
  #youwant.append('electron_en')
  #youwant.append('photon_en')
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px...
  
  from_path='./cannon_a190/'
  to_path  =from_path
  
  
  n0=30.0
  R=1.8e-6
  L=15e-6
  Ntot = np.pi*R*R*L*n0*denunit
  V=(1.0/20.0)*(1.0/15.0)*(1.0/15.0)*1.0e-18
  weight = V*denunit*n0/50.0 
  
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

  px = px[theta < 10]
  grid_x = grid_x[theta < 10]
  theta = theta[theta < 10]

  if np.size(px) == 0:
    return 0;
  theta[theta < -10] = -10
  theta[theta >  10] =  10

  color_index = abs(theta)

  print('hehe1')
  weight = np.zeros_like(grid_x)+weight
  print(max(weight),min(weight))
  ######### Script code drawing figure ################
  #for n in range(start,stop+step,step):
  #### header data ####
  #data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
  #header=data['Header']
  #time=header['time']
  #x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  #y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  #X, Y = np.meshgrid(x, y)
  from_path_2='./cannon_a190_e_div/'
  data = sdf.read(from_path_2+'q'+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time=header['time']
  x  = data['Grid/Grid_mid'].data[0]/1.0e-6
  y  = data['Grid/Grid_mid'].data[1]/1.0e-6
  z  = data['Grid/Grid_mid'].data[2]/1.0e-6

  X,Y,Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
  RR = (Y**2+Z**2)**0.5

  data_ne_in = data['Derived/Number_Density/E_in'].data/denunit+data['Derived/Number_Density/E_in_1'].data/denunit
  data_ne_out= data['Derived/Number_Density/E_out'].data/denunit+data['Derived/Number_Density/E_out_1'].data/denunit
  data_np= data['Derived/Number_Density/Ion'].data/denunit+data['Derived/Number_Density/Ion_1'].data/denunit
  data_nc= data['Derived/Number_Density/Carbon'].data/denunit+data['Derived/Number_Density/Carbon_1'].data/denunit


  data_1 = sdf.read(from_path+'e_fields'+str(n).zfill(4)+".sdf",dict=True)
  data_ex= data_1['Electric Field/'+str.capitalize('ex_averaged')].data/exunit

  for R_lim in range(5,45,5): 
      eexx = data_ex[:,RR[0,:,:]<R_lim/10.]
    #  print(eexx.shape)
      n_size = eexx[-1,:].size
      value_ex = np.sum(eexx,axis=1)/n_size

      eexx = data_ne_in[:,RR[0,:,:]<R_lim/10.]
#      print(eexx.shape)
      n_size = eexx[-1,:].size
      value_n_e_in = np.sum(eexx,axis=1)/n_size

      eexx = data_ne_out[:,RR[0,:,:]<R_lim/10.]
#      print(eexx.shape)
      n_size = eexx[-1,:].size
      value_n_e_out = np.sum(eexx,axis=1)/n_size

      eexx = data_np[:,RR[0,:,:]<R_lim/10.]
#      print(eexx.shape)
      n_size = eexx[-1,:].size
      value_n_p = np.sum(eexx,axis=1)/n_size
   
      eexx = data_nc[:,RR[0,:,:]<R_lim/10.]
#      print(eexx.shape)
      n_size = eexx[-1,:].size
      value_n_c = np.sum(eexx,axis=1)/n_size
    #    plt.subplot()
    #  plt.scatter(grid_x, px, c=color_index, s=0.03, cmap='rainbow_r', edgecolors='None', alpha=0.66)
      #plt.hist2d(grid_x, px, bins=(50, 50), range=[[10,50],[0,1.5]], cmap=plt.cm.jet, weights=weight, norm=mcolors.Normalize(vmin=-2.0, vmax=2.0))
      fig,host = plt.subplots()
      plt.hist2d(grid_x, px, bins=(100, 100), range=[[15,45],[0,1.5]], cmap='cubehelix_r', weights=weight, normed=False, norm=colors.Normalize(vmin=0, vmax=1e9))
      cbar=plt.colorbar(pad=0.1)
      cbar.set_label(r'$dN/(dxdp_x)$'+' [A.U.]',fontdict=font)
      cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(),fontsize=20)
    #plt.plot(np.linspace(-500,900,1001), np.zeros([1001]),':k',linewidth=2.5)
    #plt.plot(np.zeros([1001]), np.linspace(-500,900,1001),':k',linewidth=2.5)
    #plt.plot(np.linspace(-500,900,1001), np.linspace(-500,900,1001),'-g',linewidth=3)
    #plt.plot(np.linspace(-500,900,1001), 200-np.linspace(-500,900,1001),'-',color='grey',linewidth=3)
     #   plt.legend(loc='upper right')
      plt.xlim(15,45)
      plt.ylim(0.,1.5)
      plt.xlabel('X [$\mu m$]',fontdict=font)
      plt.ylabel('$p_x$ [m$_i$c$^2$]',fontdict=font)
      plt.xticks([15,25,35,45],fontsize=25); plt.yticks([0,0.5,1.0,1.5],fontsize=25);
    #  plt.text(-100,650,' t = '++' fs',fontdict=font)
      plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                    wspace=None, hspace=None)
      plt.title('At '+str(round(time1/1.0e-15,2))+' fs',fontdict=font)
      print('hehe')
    #plt.show()
    #lt.figure(figsize=(100,100))
    
      par1 = host.twinx()
      #par2 = host.twinx()
      #par3 = host.twinx()
    
      #par2.spines["right"].set_position(("axes", 1.05))
      #make_patch_spines_invisible(par2)
      #par2.spines["right"].set_visible(True)
    
      #par3.spines["right"].set_position(("axes", 1.1))
      #make_patch_spines_invisible(par3)
      #par3.spines["right"].set_visible(True)
    
      tkw = dict(size=25, width=2.)
    
      p1, = par1.plot(x,value_ex, "-",color='lime', label="$E_x$",linewidth=3)
      par1.set_ylim(0,5)
    
      p1, = par1.plot(x,value_n_e_in+value_n_e_out, "-b", label="$n_e$",linewidth=3)
      p1, = par1.plot(x,value_n_e_in, "-c", label="$n_e-in$",linewidth=3)
      p1, = par1.plot(x,value_n_e_out, "-",color='orange', label="$n_e-out$",linewidth=3)
    
      p1, = par1.plot(x,value_n_p+value_n_c*6, "-r", label="$Z_in_i$",linewidth=3)
    #  p3, = par3.plot(x,iden, "-r", label="Ion")
    #  p3, = par3.plot(x,iden, "-y", label="Ion")
    #  p3, = par3.plot(x,cden*6, "-g", label="Ion")
    #  par3.set_ylabel('$n^+\ [n_c]$')
    
      par1.legend(loc='middle right',fontsize=20,framealpha=1.0)
      par1.set_ylim(0,5)
      par1.set_ylabel('$E_x$ [$m_ec\omega/|e|$] & $n_e$($Z_in_i$) [$n_c$]',fontdict=font,color='r')
      par1.tick_params(axis='y',labelsize=25,colors='r')
    
      fig = plt.gcf()
      fig.set_size_inches(12, 7.5)
      fig.savefig(from_path_2+'new_cut_px_x_phase_'+str(n).zfill(4)+'_r'+str(R_lim).zfill(2)+'.png',format='png',dpi=160)
      plt.close("all")
      print('finised cut_px_x_phase_'+str(n).zfill(4)+'_r'+str(R_lim).zfill(2)+'.png')
  return 0

if __name__ == '__main__':
  start   =  1 # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=3)
  results = pool.map(processplot,inputs)
  print(results)

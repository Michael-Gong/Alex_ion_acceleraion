import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
  
def one_procedure_scatter(n,part13_id):
    data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time1=header['time']/1.0e-15
    px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
    grid_x = data['Grid/Particles/subset_Only_Ions0/Ion'].data[0]/1.0e-6      
    grid_y = data['Grid/Particles/subset_Only_Ions0/Ion'].data[1]/1.0e-6 
    grid_z = data['Grid/Particles/subset_Only_Ions0/Ion'].data[2]/1.0e-6      
    temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
    px = px[np.in1d(temp_id,part13_id)]
    grid_x = grid_x[np.in1d(temp_id,part13_id)]
    grid_y = grid_y[np.in1d(temp_id,part13_id)]
    grid_z = grid_z[np.in1d(temp_id,part13_id)]
    plt.scatter(grid_y, grid_z, c='lime', s=0.3, edgecolors='None', alpha=0.8,zorder=2)


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
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2  # j_aflven/wavelength^2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font  = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 25,  
          }  
 
  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 14,  
          }  
 
  font_size=25
  font_size2=17
  marker_size = 0.05 
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
      binsize = 200
      en_grid = np.linspace(0,1000,200)
      en_bin  = np.linspace(0,1000,201)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(en_bin[-1]-en_bin[-2])
      return (en_grid, en_value)

  c_green= matplotlib.colors.colorConverter.to_rgba('limegreen')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_wg = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_white_trans,c_green],128)

  upper = matplotlib.cm.Reds(np.arange(256))
  lower = np.ones((int(256/4),4))
  for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
  cmap = np.vstack(( lower, upper ))
  mycolor_reds = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

  upper = matplotlib.cm.binary(np.arange(256))
  lower = np.ones((int(256/4),4))
  for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
  cmap = np.vstack(( lower, upper ))
  mycolor_binary = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

  upper = matplotlib.cm.Oranges(np.arange(256))
  lower = np.ones((int(256/4),4))
  for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
  cmap = np.vstack(( lower, upper ))
  mycolor_orange = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])


  ######### Parameter you should set ###########
  from_path = './cannon_a190_v484/'
  to_path   = './cannon_a190_v484_fig/' 
  if not os.path.exists(to_path):
         os.mkdir(to_path)  
  data = sdf.read(from_path+"i_tot_loc0027.sdf",dict=True)
  px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
  theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
  gg = (px**2+py**2+pz**2+1)**0.5
  Ek = (gg-1)*1836*0.51
  part13_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
  part13_id = part13_id[ (Ek>220) & (abs(theta)<10) & (Ek<240)]
  print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))

  x_start = 460
  x_stop  = 520
  quiver_space = 5
  ######### Script code drawing figure ################
  n=17        
#### header data ####

  plt.subplot(1,3,1)
  data = sdf.read(from_path+'b_fields'+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time =header['time']
  x    = data['Grid/Grid_mid'].data[0]/1.e-6
  y    = data['Grid/Grid_mid'].data[1]/1.e-6
  z    = data['Grid/Grid_mid'].data[2]/1.e-6
  var1  = data['Magnetic Field/By_averaged'].data/bxunit
  var2  = data['Magnetic Field/Bz_averaged'].data/bxunit
  Y,Z  = np.meshgrid(y,z)
  eexx = (var1**2+var2**2)**0.5
  ex = np.sum(eexx[x_start:x_stop,:,:],axis=0)/np.size(eexx[x_start:x_stop,0,0])
  var1 = np.sum(var1[x_start:x_stop,:,:],axis=0)/np.size(var1[x_start:x_stop,0,0])
  var2 = np.sum(var2[x_start:x_stop,:,:],axis=0)/np.size(var2[x_start:x_stop,0,0])
  eee = 8#np.max([np.max(ex),abs(np.min(ex))])
  ex[ex>eee] =eee
  levels = np.linspace(0, eee, 40)
  plt.contourf(Y, Z, ex.T, levels=levels, norm=colors.Normalize(vmin=0, vmax=eee), cmap=mycolor_orange)
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_xlim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_ylim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])
  cbar = plt.colorbar(pad=0.01, ticks=np.linspace(0, eee, 3))
  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
  cbar.set_label(r'$\overline{B}_\phi\ [m_e\omega/|e|]$',fontdict=font)
  Q=plt.quiver(Y[::quiver_space,::quiver_space], Z[::quiver_space,::quiver_space], (var1.T)[::quiver_space,::quiver_space], (var2.T)[::quiver_space,::quiver_space], pivot='mid', units='x',scale=15,headwidth=4.5,width=0.01)
  plt.xlim(-2.5,2.5)
  plt.ylim(-2.5,2.5)
  plt.xlabel('$Y\ [\mu m]$',fontdict=font)
  plt.ylabel('$Z\ [\mu m]$',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size);
  plt.text(0.2,2,str(round(time/1.0e-15,1))+' fs',fontdict=font2,color='k')


  plt.subplot(1,3,2)
  data = sdf.read(from_path+'e_fields'+str(n).zfill(4)+".sdf",dict=True)
  var1  = data['Electric Field/Ey_averaged'].data/exunit
  var2  = data['Electric Field/Ez_averaged'].data/exunit
  Y,Z  = np.meshgrid(y,z)
  eexx = (var1**2+var2**2)**0.5
  ex = np.sum(eexx[x_start:x_stop,:,:],axis=0)/np.size(eexx[x_start:x_stop,0,0])
  var1 = np.sum(var1[x_start:x_stop,:,:],axis=0)/np.size(var1[x_start:x_stop,0,0])
  var2 = np.sum(var2[x_start:x_stop,:,:],axis=0)/np.size(var2[x_start:x_stop,0,0])
  eee = 8#np.max([np.max(ex),abs(np.min(ex))])
  ex[ex>eee] =eee
  levels = np.linspace(0, eee, 40)
  plt.contourf(Y, Z, ex.T, levels=levels, norm=colors.Normalize(vmin=0, vmax=eee), cmap=mycolor_reds)
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_xlim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_ylim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])
  cbar = plt.colorbar(pad=0.01, ticks=np.linspace(0, eee, 3))
  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
  cbar.set_label(r'$\overline{E}_r\ [m_ec\omega/|e|]$',fontdict=font)
#  one_procedure_scatter(n=n,part13_id=part13_id)
  Q=plt.quiver(Y[::quiver_space,::quiver_space], Z[::quiver_space,::quiver_space], (var1.T)[::quiver_space,::quiver_space], (var2.T)[::quiver_space,::quiver_space], pivot='mid', units='x',scale=15,headwidth=8)
  plt.xlim(-2.5,2.5)
  plt.ylim(-2.5,2.5)
  plt.xlabel('$Y\ [\mu m]$',fontdict=font)
  plt.ylabel('$Z\ [\mu m]$',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size); 


  plt.subplot(1,3,3)
  data = sdf.read(from_path+'current'+str(n).zfill(4)+".sdf",dict=True)
  header=data['Header']
  time =header['time']
  x    = data['Grid/Grid_mid'].data[0]/1.e-6
  y    = data['Grid/Grid_mid'].data[1]/1.e-6
  z    = data['Grid/Grid_mid'].data[2]/1.e-6
  var  = data['Current/Jx_averaged'].data/jalf
  Y,Z  = np.meshgrid(y,z)
  eexx = var #(var1**2+var2**2)**0.5
  ex = np.sum(eexx[x_start:x_stop,:,:],axis=0)/np.size(eexx[x_start:x_stop,0,0])
#    var1 = np.sum(var1[x_start:x_stop,:,:],axis=0)/np.size(var1[x_start:x_stop,0,0])
#    var2 = np.sum(var2[x_start:x_stop,:,:],axis=0)/np.size(var2[x_start:x_stop,0,0])
  eee = 10#np.max([np.max(ex),abs(np.min(ex))])
  ex[ex>eee] =eee
  ex[ex<-eee] =-eee
  levels = np.linspace(-eee, eee, 40)
  plt.contourf(Y, Z, ex.T, levels=levels, norm=colors.Normalize(vmin=-eee, vmax=eee), cmap='bwr')
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_xlim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_ylim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])
  cbar = plt.colorbar(pad=0.01, ticks=np.linspace(-eee, eee, 3))
  cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
  cbar.set_label(r'$\overline{j}_x\ [J_\alpha/\lambda^2]$',fontdict=font)
#    Q=plt.quiver(Y[::quiver_space,::quiver_space], Z[::quiver_space,::quiver_space], (var1.T)[::quiver_space,::quiver_space], (var2.T)[::quiver_space,::quiver_space], pivot='mid', units='x',scale=100,headwidth=4.5)
  plt.xlim(-2.5,2.5)
  plt.ylim(-2.5,2.5)
  plt.xlabel('$Y\ [\mu m]$',fontdict=font)
  plt.ylabel('$Z\ [\mu m]$',fontdict=font)
  plt.xticks(fontsize=font_size); 
  plt.yticks(fontsize=font_size); 
  plt.text(0.2,2,'at x='+str(round(x[(x_start+x_stop)//2],1))+' micron',fontdict=font2,color='k')


  plt.subplots_adjust(left=0.08, bottom=0.10, right=0.95, top=0.97, wspace=0.3, hspace=0.3)
  fig = plt.gcf()
  fig.set_size_inches(33, 8.5)
  fig.savefig('./er_ave_section_'+str(n).zfill(4)+'.png',format='png',dpi=160)
  plt.close("all")
  print('./er_ave_section_'+str(n).zfill(4)+'.png')


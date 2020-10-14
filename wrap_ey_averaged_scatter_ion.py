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
          'size'   : 32,  
          }  
  font2 = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 16,  
          } 
  space_1 = 5
  space_2 = 5
  font_size = 32
  marker_size=0.05
##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('crimson')
  c_blue= matplotlib.colors.colorConverter.to_rgba('mediumblue')
  c_yellow= matplotlib.colors.colorConverter.to_rgba('yellow')
  c_cyan= matplotlib.colors.colorConverter.to_rgba('cyan')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128)
  cmap_ryb= matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_yellow,c_blue],128)
  cmap_yc = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_yellow,c_white_trans,c_cyan],128) 

  c_brown = matplotlib.colors.colorConverter.to_rgba('gold')
  c_green = matplotlib.colors.colorConverter.to_rgba('springgreen')
  cmap_bg = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_brown,c_white_trans,c_white_trans,c_green],128)
  c_black = matplotlib.colors.colorConverter.to_rgba('black')
  c_white = matplotlib.colors.colorConverter.to_rgba('white')
  cmap_bw = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_white_trans,c_white,c_black],128)
   
##end for transparent colorbar##
  upper = matplotlib.cm.gist_rainbow_r(np.arange(256))
  lower = np.ones((int(256/4),4))
  for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
  cmap = np.vstack(( lower, upper ))
  mycolor_gist = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

  upper = matplotlib.cm.rainbow(np.arange(256))
  lower = np.ones((int(256/4),4))
  for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
  cmap = np.vstack(( lower, upper ))
  mycolor_rainbow = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

  upper = matplotlib.cm.Greens(np.arange(256))
  lower = np.ones((int(256/4),4))
  for i in range(3):
      lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
  cmap = np.vstack(( lower, upper ))
  mycolor_Greens = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
 
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

  def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

  c = mcolors.ColorConverter().to_rgb
  rvb_test = make_colormap([c('white'),0.03,c('white'), c('royalblue'), 0.23, c('royalblue'), c('springgreen'), 0.43, c('springgreen'),c('yellow'),0.6, c('yellow'),c('red'), 0.8, c('red'),c('firebrick')])


  color_list = ['blue','limegreen','red'] 
 
def one_procedure(n):
    data = sdf.read(from_path+'e_fields'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time  =header['time']
    x     = data['Grid/Grid_mid'].data[0]/1.0e-6
    y     = data['Grid/Grid_mid'].data[1]/1.0e-6
    z     = data['Grid/Grid_mid'].data[2]/1.0e-6
    X, Z  = np.meshgrid(x, z) 
    ey1  = data['Electric Field/Ez_averaged'].data/exunit*3.5
    ey    = np.sum(ey1[:,179:181,:],axis=1)/2
    levels = np.linspace(-20,20, 41)
    ey[ey>np.max(levels)]=np.max(levels); 
    ey[ey<np.min(levels)]=np.min(levels);

    plt.contourf(X, Z, ey.T, levels=levels, cmap='bwr', norm=colors.Normalize(vmin=np.min(levels),vmax=np.max(levels)),zorder=0)
#    cbar=plt.colorbar(pad=0.01,ticks=[-20,0,20])
#    cbar.ax.tick_params(labelsize=font2['size']) 
#    cbar.set_label('$\overline{E}_z$ [TV/m]',fontdict=font2)        
    data  = sdf.read(from_path+'q'+str(n).zfill(4)+".sdf",dict=True)
    den_ea= data['Derived/Number_Density_averaged/Electron'].data/denunit # should be averaged density
    den_e = np.sum(den_ea[:,179:181,:],axis=1)/2
    plt.contour(X, Z, den_e.T,levels=[1.0], colors='black', linewidths=2,origin='lower',zorder=1)
    plt.contour(X, Z, den_e.T,levels=[30.0], colors='blue', linewidths=2,origin='lower',zorder=1)
#         plt.contour(X, Z, den_e.T,levels=[1.0], colors='black', linewidths=2.,origin='lower')

def one_procedure_scatter(n,part13_id):
    data = sdf.read(from_path+"i_tot_loc"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time1=header['time']/1.0e-15
    px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
    grid_x = data['Grid/Particles/subset_Only_Ions0/Ion'].data[0]/1.0e-6      
    grid_z = data['Grid/Particles/subset_Only_Ions0/Ion'].data[2]/1.0e-6      
    temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
    px = px[np.in1d(temp_id,part13_id)]
    grid_x = grid_x[np.in1d(temp_id,part13_id)]
    grid_z = grid_z[np.in1d(temp_id,part13_id)]
    plt.scatter(grid_x, grid_z, c='lime', s=0.03, edgecolors='None', alpha=0.8,zorder=2)

    
 
if __name__ == '__main__':
    from_path = './cannon_a190_v484/'
    to_path   = './cannon_a190_v484_fig/'
    ######### Script code drawing figure ################
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

    plt.subplot(3,1,1)
    one_procedure(n=15)
    one_procedure_scatter(n=15,part13_id=part13_id)
#         plt.contour(X, Z, den_e.T,levels=[30.0], colors='white', linewidths=2.,origin='lower')
    plt.xlim(8,28)
    plt.ylim(-4.5,4.5)
#    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Z [$\mu m$]',fontdict=font)
    plt.xticks([10,15,20,25],fontsize=0.01); 
    plt.yticks([-3,0,3],fontsize=font_size);
#    plt.title('t='+str(round(time/1e-15,2))+' fs',fontdict=font)

    plt.subplot(3,1,2)
    one_procedure(n=17)
    one_procedure_scatter(n=17,part13_id=part13_id)
#         plt.contour(X, Z, den_e.T,levels=[30.0], colors='white', linewidths=2.,origin='lower')
    plt.xlim(8,28)
    plt.ylim(-4.5,4.5)
#    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Z [$\mu m$]',fontdict=font)
    plt.xticks([10,15,20,25],fontsize=0.01); 
    plt.yticks([-3,0,3],fontsize=font_size);
#    plt.title('t='+str(round(time/1e-15,2))+' fs',fontdict=font)

    plt.subplot(3,1,3)
    one_procedure(n=19)
    one_procedure_scatter(n=19,part13_id=part13_id)
    plt.xlim(8,28)
    plt.ylim(-4.5,4.5)
    plt.xlabel('X [$\mu m$]',fontdict=font)
    plt.ylabel('Z [$\mu m$]',fontdict=font)
    plt.xticks([10,15,20,25],fontsize=font_size); 
    plt.yticks([-3,0,3],fontsize=font_size);
#    plt.title('t='+str(round(time/1e-15,2))+' fs',fontdict=font)

    
    plt.subplots_adjust(left=0.14, bottom=0.12, right=0.99, top=0.99, wspace=0.02, hspace=0.05)
    fig = plt.gcf()
    fig.set_size_inches(10, 10)
    fig.savefig(to_path+'wrap_ey_scatter_ion.png',format='png',dpi=160)
    plt.close("all")
    print(to_path+'wrap_ey_scatter_ion.png')

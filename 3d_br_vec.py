import numpy as np
import sdf 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from functools import reduce
import multiprocessing as mp
import sys, getopt
import os, time
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import gc

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
        'size'   : 25 ,  
        }  
font_size = 20

upper = matplotlib.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
  lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
mycolormap = matplotlib.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

def rebin3d(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1],shape[2],a.shape[2]//shape[2]
    return a.reshape(sh).mean(-1).mean(3).mean(1)

def reg_cmap_transparent(iname,alpha):
    oname = iname + '_transparent'
    cmap = plt.get_cmap(iname)
    values = np.linspace(0,1,256)
    colors = cmap(values)
    for i in range(256):
        colors[i][3] = alpha[i]
    colorlist = [(values[i],colors[i]) for i in range(256)]
    cmap = plt.cm.colors.LinearSegmentedColormap.from_list(oname,colorlist)
    plt.cm.register_cmap(cmap=cmap)
    return cmap

def create_alpha(func):
    return [ 1 if func(i)>1 else 0 if func(i)<0 else func(i) for i in range(256)]


def processplot(n): 
  #for n in range(start,stop+step,step):
    from_path='./'
    to_path='./'
    x_start=100; x_stop=700; y_start=60; y_stop=300; z_start=60; z_stop=300;
    x_size = x_stop-x_start; y_size = y_stop-y_start; z_size = z_stop-z_start
    name = 'Br_averaged'

    data = sdf.read(from_path+'b_fields'+str(n).zfill(4)+'.sdf',dict=True)
    header=data['Header']
    time =header['time']
    x    = data['Grid/Grid_mid'].data[0]/1.e-6
    y    = data['Grid/Grid_mid'].data[1]/1.e-6
    z    = data['Grid/Grid_mid'].data[2]/1.e-6
    var1  = data['Magnetic Field/By_averaged'].data/bxunit
    var2  = data['Magnetic Field/Bz_averaged'].data/bxunit
    var0  = data['Magnetic Field/Bx_averaged'].data/bxunit

    var   = (var1**2+var2**2+var0**2)**0.5

    X, Y, Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
    var  = var[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    var1  = var1[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    var2  = var2[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    var0  = var0[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    X    =  X[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    Y    =  Y[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    Z    =  Z[x_start:x_stop,y_start:y_stop,z_start:z_stop]

    plotkws = {'marker':'.','edgecolors':'none'}
    norm = None

    index = 6.0
    _abs  = False # True is for ex; False is for density
    log   = False
    elev  = None
    azim  = None

    print('here1')
    fig = plt.figure()
    ax  = plt.axes(projection='3d')
    ax.view_init(elev=elev, azim=azim)
    print('here2')
    jump_space = 20
    vec_X = X[::jump_space,::jump_space,::jump_space]
    vec_Y = Y[::jump_space,::jump_space,::jump_space]
    vec_Z = Z[::jump_space,::jump_space,::jump_space]
    vec_var0 = var0[::jump_space,::jump_space,::jump_space]
    vec_var1 = var1[::jump_space,::jump_space,::jump_space]
    vec_var2 = var2[::jump_space,::jump_space,::jump_space]
    vec_var  = (vec_var0**2+vec_var1**2+vec_var2**2)**0.5
    vec_index = vec_var > 4.0
    # Color by azimuthal angle
    # Flatten and normalize
    c = (vec_var[vec_index] - 4.0) / 16
    c = c.ravel()
    # Repeat for each body line and two head lines
    c = np.concatenate((c, np.repeat(c, 2)))
    # Colormap
    c = plt.cm.jet(c)
    im = ax.quiver(vec_X[vec_index], vec_Y[vec_index], vec_Z[vec_index], vec_var0[vec_index], vec_var1[vec_index], vec_var2[vec_index], colors=c, length=1.3, lw=2, arrow_length_ratio=0.5, alpha=0.8, normalize=True)#, norm=mcolors.Normalize(vmin=8, vmax=20))
    print('here3')
    ax.set_xlabel('\n\nX'+ '[$\mu m$]',fontdict=font)
    ax.set_ylabel('\n\nY'+ '[$\mu m$]',fontdict=font)
    ax.set_zlabel('\n\nZ'+ '[$\mu m$]',fontdict=font)

    ax.set_xlim([x_start/20-5,x_stop/20-5])
    ax.set_ylim([-(y_stop-y_start)/2/15-5,(y_stop-y_start)/2/15+5])
    ax.set_zlim([-(z_stop-z_start)/2/15-5,(z_stop-z_start)/2/15+5])

    
    #cbar=plt.colorbar(im, ticks=np.linspace(0,30, 5) ,pad=0.01)
    #cbar=plt.colorbar(im, ticks=np.linspace(np.min(color_index), np.max(color_index), 5) ,pad=0.01)
    #cbar=plt.colorbar(im, pad=0.01)
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    #cbar.set_label(name+r'$[m_e\omega/|e|]$',fontdict=font)
    #cbar.set_clim(300,600)
    print('here4')

    #print('here4')

    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(font_size)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(font_size)
    for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(font_size)



#    ax.scatter(X,Z,c=var, zdir='y',zs=(y_stop-y_start)/2/12, marker='.', edgecolors='none', cmap=cmap)
#    ax.scatter(X,Y,c=var, zdir='z',zs=-(z_stop-z_start)/2/12,  marker='.', edgecolors='none', cmap=cmap)
#    ax.scatter(Y,Z,c=var, zdir='x',zs=x_start/20-5,  marker='.', edgecolors='none', cmap=cmap)

    #plot for y_z plane
    var1  = data['Magnetic Field/By_averaged'].data/bxunit
    var2  = data['Magnetic Field/Bz_averaged'].data/bxunit

    Y,Z  = np.meshgrid(y,z,indexing='ij')
    ex = (var1**2.+var2**2.)**0.5
    ex = (ex[420-1,:,:]+ex[420,:,:])/2
    ex = ex[y_start:y_stop,z_start:z_stop]
    eee = 15#np.max([np.max(ex),abs(np.min(ex))])
    ex[ex>eee] =eee
    Y  = Y[y_start:y_stop,z_start:z_stop]
    Z  = Z[y_start:y_stop,z_start:z_stop]
    levels = np.linspace(0, eee, 40)
    im2=ax.contourf(ex.T, Y.T, Z.T, levels=levels, norm=mcolors.Normalize(vmin=0, vmax=eee), cmap=mycolormap, zdir='x', offset=x_start/20-5)
#    im2=ax.quiver(Y[ex>8], Z[ex>8], var1[ex>8], var2[ex>8], colors='black', length=1., lw=1, arrow_length_ratio=0.5, alpha=1, normalize=False)#, norm=mcolors.Normalize(vmin=8, vmax=20))
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_xlim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_ylim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])
    cbar = plt.colorbar(im2,  ticks=np.linspace(0, eee, 3))
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label(name+r'$[m_e\omega/|e|]$',fontdict=font)
    

    #plot for x_z plane
    X,Z = np.meshgrid(x,z,indexing='ij')
    eexx = var1
    ex = (eexx[:,(y_start+y_stop)//2-1,:]+eexx[:,(y_start+y_stop)//2,:])/2
    ex = ex[x_start:x_stop,z_start:z_stop]
    X  = X[x_start:x_stop,z_start:z_stop]
    Z  = Z[x_start:x_stop,z_start:z_stop]
    if np.min(ex.T) == np.max(ex.T):
         #continue
         return
    eee = 30
    levels = np.linspace(-eee, eee, 40)
    ex[ex>eee] = eee
    ex[ex<-eee] = -eee
    ax.contourf(X.T, ex.T, Z.T, levels=levels, norm=mcolors.Normalize(vmin=-eee, vmax=eee), cmap=cm.bwr, zdir='y', offset=(y_stop-y_start)/2/15+5)
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_ylim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_ylim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])

    #plot for x_y plane
    X,Y = np.meshgrid(x,y,indexing='ij')
    ex = var2
    ex = (ex[:,:,(z_start+z_stop)//2-1]+ex[:,:,(z_start+z_stop)//2])/2
    ex = ex[x_start:x_stop,y_start:y_stop]
    X  = X[x_start:x_stop,y_start:y_stop]
    Y  = Y[x_start:x_stop,y_start:y_stop]
    if np.min(ex.T) == np.max(ex.T):
         #continue
         return
    eee = 30
    levels = np.linspace(-eee, eee, 40)
    ex[ex>eee] = eee
    ex[ex<-eee] = -eee
    im2=ax.contourf(X.T, Y.T, ex.T, levels=levels, norm=mcolors.Normalize(vmin=-eee, vmax=eee), cmap=cm.bwr, zdir='z', offset=-(z_stop-z_start)/2/15-5)
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_ylim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_zlim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])
    cbar = plt.colorbar(im2,  ticks=np.linspace(-eee, eee, 5))
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label(r'$B_\theta\ [m_e\omega/|e|]$',fontdict=font)

    #ax.scatter(X,Z,c=var,**plotkws ,zdir='y',zs=4)
    #ax.scatter(X,Y,c=var,**plotkws, zdir='z',zs=-4)
    #ax.scatter(Y,Z,c=var,**plotkws, zdir='x',zs=15)

    plt.show()
    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    #ax.grid(linestyle='None', linewidth='0.5', color='white')
    plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                    wspace=None, hspace=None)
  #  plt.title('At '+str(round(time/1.0e-15,2))+' fs',fontdict=font)


    fig = plt.gcf()
    fig.set_size_inches(20, 10.5)
    fig.savefig(to_path+'3d_'+name+'_vec_'+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    print('finised '+str(n).zfill(4))
    print('here5')
    gc.collect()

if __name__ == '__main__':
  start   =  10 # start time
  stop    =  18  # end time
  step    =  2  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=1)
  results = pool.map(processplot,inputs)
  print(results)

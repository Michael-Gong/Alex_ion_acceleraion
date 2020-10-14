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
        'size'   : 25,  
        }  
font_size = 20

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
    from_path='./cannon_a190/'
    to_path='./cannon_a190/'
    x_start=100; x_stop=700; y_start=60; y_stop=300; z_start=60; z_stop=300;
    x_size = x_stop-x_start; y_size = y_stop-y_start; z_size = z_stop-z_start
    name = 'Electron_density'

    data = sdf.read(from_path+'q'+str(n).zfill(4)+'.sdf',dict=True)
    header=data['Header']
    time =header['time']
    x    = data['Grid/Grid_mid'].data[0]/1.e-6
    y    = data['Grid/Grid_mid'].data[1]/1.e-6
    z    = data['Grid/Grid_mid'].data[2]/1.e-6
    var  = data['Derived/Number_Density/Electron'].data/denunit

    X, Y, Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
    var  = var[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    X    =  X[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    Y    =  Y[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    Z    =  Z[x_start:x_stop,y_start:y_stop,z_start:z_stop]
   
    var = rebin3d(var, (x_size//2, y_size//2, z_size//2))
    X = rebin3d(X, (x_size//2, y_size//2, z_size//2))
    Y = rebin3d(Y, (x_size//2, y_size//2, z_size//2))
    Z = rebin3d(Z, (x_size//2, y_size//2, z_size//2))

    var  = var.reshape(np.size(var))
    var[var > 30] =0
    X    = X.reshape(np.size(X))
    Y    = Y.reshape(np.size(Y))
    Z    = Z.reshape(np.size(Z))

    plotkws = {'marker':'.','edgecolors':'none'}
    norm = None

    index = 6.0
    _abs  = False # True is for ex; Flase is for density
    log   = False
    elev  = None
    azim  = None

    if _abs == True:
        norm = 0
        _min = max(np.max(var),np.min(var))**(0.002**(1.0/index)) if log else max(np.max(var),np.min(var))*0.002**(1.0/index)
        plt.set_cmap(reg_cmap_transparent('bwr',create_alpha(lambda x:abs(x/127.5-1)**index)))
    else:
        _min = np.max(var)**(0.002**(1.0/index)) if log else np.max(var)*0.002**(1.0/index)
        plt.set_cmap(reg_cmap_transparent('hsv_r',create_alpha(lambda x:abs(x/255.0)**index)))

        #special code
        _min = max(_min,1.1e27*0.8)
    #    var._cutrange(lambda x : x[1] < 3)

    if log:
        plotkws['norm'] = matplotlib.colors.LogNorm()

    #var.cutrange(_min=_min,_abs=_abs)
    #point_scatter3D(var,norm=norm,plotkws=plotkws)
    #def point_scatter3D(var,elev=None,azim=None,hold=False,iso=False,norm=None,plotkws={}):
    cmap = plt.get_cmap()
    if norm is not None:
        v0 = np.min(var) - norm
        v1 = np.max(var) - norm
        if abs(v0/v1) > 1:
            low = 0
            high = 0.5 * (1 - v1/v0)
        else:
            low = 0.5 * (1 + v0/v1)
            high = 1.0

        cmap = plt.cm.colors.LinearSegmentedColormap.from_list('tr',
                cmap(np.linspace(low,high,256)))

#    print('here1')
    fig = plt.figure()
    ax  = plt.axes(projection='3d')
    ax.view_init(elev=elev, azim=azim)
#    print('here2')
    im = ax.scatter(X, Y, Z, c=var, cmap=cmap, **plotkws)
#    print('here3')
    ax.set_xlabel('\n\nX'+ '[$\mu m$]',fontdict=font)
    ax.set_ylabel('\n\nY'+ '[$\mu m$]',fontdict=font)
    ax.set_zlabel('\n\nZ'+ '[$\mu m$]',fontdict=font)

    ax.set_xlim([x_start/20-5,x_stop/20-5])
    ax.set_ylim([-(y_stop-y_start)/2/15-5,(y_stop-y_start)/2/15+5])
    ax.set_zlim([-(z_stop-z_start)/2/15-5,(z_stop-z_start)/2/15+5])

    
    #cbar=plt.colorbar(im, ticks=np.linspace(np.min(color_index), np.max(color_index), 5) ,pad=0.01)
    cbar=plt.colorbar(im, pad=0.01)
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label(r'$n_e\ [n_c]$',fontdict=font)
    #cbar.set_clim(300,600)

    #print('here4')

    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(font_size)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(font_size)
    for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(font_size)

    #plot for y_z plane
    Y,Z  = np.meshgrid(y,z,indexing='ij')
    eexx  = data['Derived/Number_Density/Electron'].data/denunit
    ex = (eexx[420-1,:,:]+eexx[420,:,:])/2
    ex = ex[y_start:y_stop,z_start:z_stop]
    eee = 30 # np.max([np.max(ex),abs(np.min(ex))])
    ex[ex>30] =30
    Y  = Y[y_start:y_stop,z_start:z_stop]
    Z  = Z[y_start:y_stop,z_start:z_stop]
    levels = np.linspace(0, eee, 40)
    im2=ax.contourf(ex.T, Y.T, Z.T, levels=levels, norm=mcolors.Normalize(vmin=0, vmax=eee), cmap=cm.pink_r, zdir='x', offset=x_start/20-5)
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_xlim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_ylim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])
    cbar = plt.colorbar(im2,  ticks=np.linspace(0, eee, 3))
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label(r'$n_e\ [n_c]$',fontdict=font)
    

    #plot for x_z plane
    X,Z = np.meshgrid(x,z,indexing='ij')
    eexx  = data['Derived/Number_Density/Electron'].data/denunit
    ex = (eexx[:,(y_start+y_stop)//2-1,:]+eexx[:,(y_start+y_stop)//2,:])/2
    ex = ex[x_start:x_stop,z_start:z_stop]
    X  = X[x_start:x_stop,z_start:z_stop]
    Z  = Z[x_start:x_stop,z_start:z_stop]
    if np.min(ex.T) == np.max(ex.T):
         #continue
         return
    eee = 30
    ex[ex>eee] = eee
    levels = np.linspace(0, eee, 40)
    ax.contourf(X.T, ex.T, Z.T, levels=levels, norm=mcolors.Normalize(vmin=0, vmax=eee), cmap=cm.pink_r, zdir='y', offset=(y_stop-y_start)/2/15+5)
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_ylim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_ylim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])

    #plot for x_y plane
    X,Y = np.meshgrid(x,y,indexing='ij')
    eexx  = data['Derived/Number_Density/Electron'].data/denunit
    ex = (eexx[:,:,(z_start+z_stop)//2-1]+eexx[:,:,(z_start+z_stop)//2])/2
    ex = ex[x_start:x_stop,y_start:y_stop]
    X  = X[x_start:x_stop,y_start:y_stop]
    Y  = Y[x_start:x_stop,y_start:y_stop]
    if np.min(ex.T) == np.max(ex.T):
         #continue
         return
    eee = 30
    ex[ex>eee] = eee
    levels = np.linspace(0, eee, 40)
    im2=ax.contourf(X.T, Y.T, ex.T, levels=levels, norm=mcolors.Normalize(vmin=0, vmax=eee), cmap=cm.pink_r, zdir='z', offset=-(z_stop-z_start)/2/15-5)
#    ax.set_xlim([x_start/20-5,x_stop/20-5])
#    ax.set_ylim([-(y_stop-y_start)/2/12,(y_stop-y_start)/2/12])
#    ax.set_zlim([-(z_stop-z_start)/2/12,(z_stop-z_start)/2/12])
    cbar = plt.colorbar(im2,  ticks=np.linspace(-eee, eee, 5))
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label(r'$n_e\ [n_c]$',fontdict=font)

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
    fig.savefig(to_path+'3d_'+name+str(n).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    print('finised '+str(n).zfill(4))
    #print('here5')

if __name__ == '__main__':
  start   =  1 # start time
  stop    =  31  # end time
  step    =  1  # the interval or step
    
  inputs = range(start,stop+step,step)
  pool = mp.Pool(processes=5)
  results = pool.map(processplot,inputs)
  print(results)

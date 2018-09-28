%matplotlib inline
import numpy as np
import sdf 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from __future__ import print_function
from functools import reduce
import multiprocessing
import sys, getopt
import os, time
from plot_class import *
from mpl_toolkits.mplot3d import Axes3D

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
font_size = 16
to_path='./'

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
    """alpha可见的最低值是0.002"""
    return [ 1 if func(i)>1 else 0 if func(i)<0 else func(i) for i in range(256)]

x_start=400; x_stop=500; y_start=120; y_stop=240; z_start=120; z_stop=240;
x_size = x_stop-x_start; y_size = y_stop-y_start; z_size = z_stop-z_start
name = 'Ion'

for n in range(12,30,2):
    data = sdf.read('q'+str(n).zfill(4)+'.sdf',dict=True)
    header=data['Header']
    time =header['time']
    x    = data['Grid/Grid_mid'].data[0]/1.e-6
    y    = data['Grid/Grid_mid'].data[1]/1.e-6
    z    = data['Grid/Grid_mid'].data[2]/1.e-6
    var  = data['Derived/Number_Density/'+name].data/denunit

    X, Y, Z = np.meshgrid(x, y, z, sparse=False, indexing='ij')
    var  = var[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    X    =  X[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    Y    =  Y[x_start:x_stop,y_start:y_stop,z_start:z_stop]
    Z    =  Z[x_start:x_stop,y_start:y_stop,z_start:z_stop]

    #X = rebin(X, (x_size//2, y_size//2, z_size//2))
    #Y = rebin(Y, (x_size//2, y_size//2, z_size//2))
    #Z = rebin(Z, (x_size//2, y_size//2, z_size//2))

    var  = var.reshape(np.size(var))
    X    = X.reshape(np.size(X))
    Y    = Y.reshape(np.size(Y))
    Z    = Z.reshape(np.size(Z))

    plotkws = {'marker':'.','edgecolors':'none'}
    norm = None

    index = 3
    _abs  = False # True is for ex; Flase is for density
    log   = False
    elev  = None
    azim  = None

    if _abs:
        norm = 0
        _min = max(np.max(var),np.min(var))**(0.002**(1.0/index)) if log else max(np.max(var),np.min(var))*0.002**(1.0/index)
        plt.set_cmap(reg_cmap_transparent('bwr',create_alpha(lambda x:abs(x/127.5-1)**index)))
    else:
        _min = np.max(var)**(0.002**(1.0/index)) if log else np.max(var)*0.002**(1.0/index)
        plt.set_cmap(reg_cmap_transparent('plasma_r',create_alpha(lambda x:abs(x/255.0)**index)))

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
        v0 = np.min(var.data) - norm
        v1 = np.max(var.data) - norm
        if abs(v0/v1) > 1:
            low = 0
            high = 0.5 * (1 - v1/v0)
        else:
            low = 0.5 * (1 + v0/v1)
            high = 1.0

        cmap = plt.cm.colors.LinearSegmentedColormap.from_list('tr',
                cmap(np.linspace(low,high,256)))

    print('here1')
    fig = plt.figure()
    ax  = plt.axes(projection='3d')
    ax.view_init(elev=elev, azim=azim)
    print('here2')
    im = ax.scatter(X, Y, Z, c=var, cmap=cmap, **plotkws)
    print('here3')
    ax.set_xlabel('X'+ '[$\mu m$]')
    ax.set_ylabel('y'+ '[$\mu m$]')
    ax.set_zlabel('Z'+ '[$\mu m$]')

    ax.set_xlim([15,20])
    ax.set_ylim([-4.,4])
    ax.set_zlim([-4.,4])

    
    #cbar=plt.colorbar(im, ticks=np.linspace(np.min(color_index), np.max(color_index), 5) ,pad=0.01)
    cbar=plt.colorbar(im, pad=0.01)
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    cbar.set_label(name+r'$[n_c]$',fontdict=font)
    #cbar.set_clim(300,600)

    print('here4')

    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(font_size)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(font_size)
    for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(font_size)

    ax.scatter(X,Z,c=var,**plotkws ,zdir='y',zs=4)
    ax.scatter(X,Y,c=var,**plotkws, zdir='z',zs=-4)
    ax.scatter(Y,Z,c=var,**plotkws, zdir='x',zs=15)


    ax.grid(linestyle='--', linewidth='0.5', color='grey')
    plt.subplots_adjust(left=0.16, bottom=None, right=0.97, top=None,
                    wspace=None, hspace=None)
    plt.title('At '+str(round(time/1.0e-15,2))+' fs',fontdict=font)


    fig = plt.gcf()
    fig.set_size_inches(12, 10.5)
    fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=320)
    plt.close("all")
    #print('finised '+str(n).zfill(4))
    print('here5')

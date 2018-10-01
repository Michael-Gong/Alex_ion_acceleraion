import sdf
#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
#from matplotlib import colors, ticker, cm
#from matplotlib.mlab import bivariate_normal
#from optparse import OptionParser
#import os
#from colour import Color

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
#print 'electric field unit: '+str(exunit)
#print 'magnetic field unit: '+str(bxunit)
#print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
    	'weight' : 'normal',  
        'size'   : 20,  
       }  


from_path='./cannon_a190/'
to_path='./txt_cannon/'

data = sdf.read(from_path+"i_tot_loc0027.sdf",dict=True)
#grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
theta = np.arctan2((py**2+pz**2)**0.5,px)*180.0/np.pi
gg = (px**2+py**2+pz**2+1)**0.5
Ek = (gg-1)*1836*0.51
part13_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
#  part13_id = part13_id[ (Ek>225) & (abs(theta)<10) & (Ek<245)]
part13_id = part13_id[ (Ek>234.8) & (abs(theta)<10) & (Ek<235.2)]
print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))

  

data = sdf.read(from_path+"i_tot_loc0010.sdf",dict=True)
part00_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
part13_id = np.intersect1d(part00_id,part13_id)
print('after intersect i_tot_loc0010.sdf part_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))


######### Parameter you should set ###########
start   =  10  # start time
stop    = 27  # end time
step    =  1  # the interval or step

#  if (os.path.isdir('jpg') == False):
#    os.mkdir('jpg')
######### Script code drawing figure ################
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+'i_tot_loc'+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    if ( n==start ):
        part_id = data['Particles/ID/subset_Only_Ions0/Ion'].data
    else:
        part_id = np.intersect1d(data['Particles/ID/subset_Only_Ions0/Ion'].data, part_id)
    print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))

part_id = np.intersect1d(part_id,part13_id)
print('After intersecting with 0013.sdf')
print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))


######### Parameter you should set ###########

px_2d = np.zeros([part_id.size,stop-start+1])
py_2d = np.zeros([part_id.size,stop-start+1])
pz_2d = np.zeros([part_id.size,stop-start+1])
xx_2d = np.zeros([part_id.size,stop-start+1])
yy_2d = np.zeros([part_id.size,stop-start+1])
zz_2d = np.zeros([part_id.size,stop-start+1])
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+'i_tot_loc'+str(n).zfill(4)+".sdf",dict=True)
    px = data['Particles/Px/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
    py = data['Particles/Py/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
    pz = data['Particles/Pz/subset_Only_Ions0/Ion'].data/(1836*m0*v0)
    grid_x = data['Grid/Particles/subset_Only_Ions0/Ion'].data[0]/wavelength
    grid_y = data['Grid/Particles/subset_Only_Ions0/Ion'].data[1]/wavelength
    grid_z = data['Grid/Particles/subset_Only_Ions0/Ion'].data[2]/wavelength
    temp_id = data['Particles/ID/subset_Only_Ions0/Ion'].data

    px = px[np.in1d(temp_id,part_id)]
    py = py[np.in1d(temp_id,part_id)]
    pz = pz[np.in1d(temp_id,part_id)]
    grid_x = grid_x[np.in1d(temp_id,part_id)]
    grid_y = grid_y[np.in1d(temp_id,part_id)]
    grid_z = grid_z[np.in1d(temp_id,part_id)]

    temp_id = temp_id[np.in1d(temp_id,part_id)]


    for ie in range(part_id.size):
        px_2d[ie,n-start] = px[temp_id==part_id[ie]]
        py_2d[ie,n-start] = py[temp_id==part_id[ie]]
        pz_2d[ie,n-start] = pz[temp_id==part_id[ie]]
        xx_2d[ie,n-start] = grid_x[temp_id==part_id[ie]]
        yy_2d[ie,n-start] = grid_y[temp_id==part_id[ie]]
        zz_2d[ie,n-start] = grid_z[temp_id==part_id[ie]]
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

np.savetxt(to_path+'px2d_x.txt',px_2d)
np.savetxt(to_path+'py2d_x.txt',py_2d)
np.savetxt(to_path+'pz2d_x.txt',pz_2d)
np.savetxt(to_path+'xx2d_x.txt',xx_2d)
np.savetxt(to_path+'yy2d_x.txt',yy_2d)
np.savetxt(to_path+'zz2d_x.txt',zz_2d)





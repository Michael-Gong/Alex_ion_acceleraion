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

from_path = './cannon_a190_e_div_2/'

data = sdf.read(from_path+"new_loc0032.sdf",dict=True)
header=data['Header']
time1=header['time']
grid_x = data['Grid/Particles/subset_New_particles/E_in_1'].data[0]/1.0e-6      
grid_y = data['Grid/Particles/subset_New_particles/E_in_1'].data[1]/1.0e-6      
grid_z = data['Grid/Particles/subset_New_particles/E_in_1'].data[2]/1.0e-6      
px     = data['Particles/Px/subset_New_particles/E_in_1'].data/(m0*v0)
py     = data['Particles/Py/subset_New_particles/E_in_1'].data/(m0*v0)
pz     = data['Particles/Pz/subset_New_particles/E_in_1'].data/(m0*v0)
gg     = (px**2+py**2+pz**2+1.)**0.5
rr     = (grid_y**2+grid_z**2)**0.5
part13_id = data['Particles/ID/subset_New_particles/E_in_1'].data
part13_id = part13_id[ (rr<4.0) & (abs(grid_x-29.5)<2.5) ]
print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))

#start   =  0  # start time
#stop    =  31  # end time
#step    =  1  # the interval or step
#  if (os.path.isdir('jpg') == False):
#    os.mkdir('jpg')
######### Script code drawing figure ################
#for n in range(start,stop+step,step):
#    #### header data ####
#    data = sdf.read(from_path+'new_loc'+str(n).zfill(4)+".sdf",dict=True)
#    header=data['Header']
#    time=header['time']
#    part_id = np.intersect1d(data['Particles/ID/subset_New_particles/E_in_1'].data, part13_id)
#    print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))
#
#part_id = np.intersect1d(part_id,part13_id)
#print('After intersecting with 0023.sdf')
#print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))
part_id = part13_id

######### Parameter you should set ###########
start   =  0  # start time
stop    =  32  # end time
step    =  1  # the interval or step

px_2d = np.zeros([part_id.size,(stop-start)/step+1])
py_2d = np.zeros([part_id.size,(stop-start)/step+1])
pz_2d = np.zeros([part_id.size,(stop-start)/step+1])
xx_2d = np.zeros([part_id.size,(stop-start)/step+1])
yy_2d = np.zeros([part_id.size,(stop-start)/step+1])
zz_2d = np.zeros([part_id.size,(stop-start)/step+1])
ww_2d = np.zeros([part_id.size,(stop-start)/step+1])
for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+"new_loc"+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time1=header['time']
    grid_x = data['Grid/Particles/subset_New_particles/E_in_1'].data[0]/1.0e-6      
    grid_y = data['Grid/Particles/subset_New_particles/E_in_1'].data[1]/1.0e-6      
    grid_z = data['Grid/Particles/subset_New_particles/E_in_1'].data[2]/1.0e-6      
    px     = data['Particles/Px/subset_New_particles/E_in_1'].data/(m0*v0)
    py     = data['Particles/Py/subset_New_particles/E_in_1'].data/(m0*v0)
    pz     = data['Particles/Pz/subset_New_particles/E_in_1'].data/(m0*v0)
    ww     = data['Particles/Weight/subset_New_particles/E_in_1'].data
    temp_id = data['Particles/ID/subset_New_particles/E_in_1'].data

    px = px[np.in1d(temp_id,part_id)]
    py = py[np.in1d(temp_id,part_id)]
    pz = pz[np.in1d(temp_id,part_id)]
    grid_x = grid_x[np.in1d(temp_id,part_id)]
    grid_y = grid_y[np.in1d(temp_id,part_id)]
    grid_z = grid_z[np.in1d(temp_id,part_id)]
    ww     = ww[np.in1d(temp_id,part_id)]
    temp_id = temp_id[np.in1d(temp_id,part_id)]

    for ie in range(part_id.size):
        px_2d[ie,(n-start)/step] = px[temp_id==part_id[ie]]
        py_2d[ie,(n-start)/step] = py[temp_id==part_id[ie]]
        pz_2d[ie,(n-start)/step] = pz[temp_id==part_id[ie]]
        xx_2d[ie,(n-start)/step] = grid_x[temp_id==part_id[ie]]
        yy_2d[ie,(n-start)/step] = grid_y[temp_id==part_id[ie]]
        zz_2d[ie,(n-start)/step] = grid_z[temp_id==part_id[ie]]
        ww_2d[ie,(n-start)/step] = ww[temp_id==part_id[ie]]
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
from_path = './cannon_a190_e_div_2/'

np.savetxt(from_path+'px3d_in.txt',px_2d)
np.savetxt(from_path+'py3d_in.txt',py_2d)
np.savetxt(from_path+'pz3d_in.txt',py_2d)
np.savetxt(from_path+'xx3d_in.txt',xx_2d)
np.savetxt(from_path+'yy3d_in.txt',yy_2d)
np.savetxt(from_path+'zz3d_in.txt',zz_2d)
np.savetxt(from_path+'ww3d_in.txt',ww_2d)






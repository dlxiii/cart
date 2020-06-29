#!/bin/env python
import time
import numpy as np
import matplotlib.tri as tri
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.io import netcdf_file

nf = netcdf_file('tokyobay_0001.nc', 'r', False)
x = nf.variables['lon'][:]
y = nf.variables['lat'][:]
h = nf.variables['h'][:]
u = nf.variables['u'][:]
v = nf.variables['v'][:]
xc = nf.variables['lonc'][:]
yc = nf.variables['latc'][:]
nv = nf.variables['nv'][:]
zeta = nf.variables['zeta'][:]
nele = nf.dimensions['nele']
node = nf.dimensions['node']
nt = zeta.shape[0]
nf.close()

triang = tri.Triangulation(x,y,nv.transpose()-1)

plt.figure()
plt.triplot(triang, color='k', linewidth = 0.2)
plt.tripcolor(triang, h)
ax = plt.gca()
ax.set_aspect('equal')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')
plt.title('depth (m)')
plt.colorbar()

plt.figure()
plt.triplot(triang, color='k', linewidth = 0.2)
plt.tripcolor(triang, zeta[-1,:])
ax = plt.gca()
ax.set_aspect('equal')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')
plt.title('zeta (m)')
plt.colorbar()

plt.figure()
plt.triplot(triang, color='k', linewidth = 0.2)
ind = range(0,nele,5)
plt.quiver(xc[ind], yc[ind], u[-1,0,ind], v[-1,0,ind])
#plt.tripcolor(triang, u[-1,0,:])
ax = plt.gca()
ax.set_aspect('equal')
plt.xlabel('Longitude (degrees)')
plt.ylabel('Latitude (degrees)')

plt.show()


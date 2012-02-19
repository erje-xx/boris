#! /usr/bin/env python

import matplotlib.pyplot as plt
from numpy import pi, deg2rad, linspace, meshgrid, sqrt, arctan2, where, fft 
from boris_lib import vomega, vdwdk

# material parameters
Ms = 139260.5752 
H = 55704.230082
L = 5.1e-06
gamma = 35185.83772
alpha = 3.1e-16
theta = deg2rad(90.0)
trans_n = 0.0
mp = [Ms, H, L, gamma, alpha, theta, trans_n]

# k_max = no of nodes*smallest k for sample width W
W = 0.0013
k0 = pi/W
k_max = 100.0*k0

# freq, epsilon for emission pattern [GHz]
freq = 3.65
epsilon = 0.01
# GHz -> Hz
freq = freq*10**9
epsilon = epsilon*10**9

#plt.subplot(231)
#plt.title('1D Dispersion')
#plt.xlabel('Wavevector (1/cm)')
#plt.ylabel('Frequency (GHz)')
#plt.grid(True)

#nsteps = 100
#k_BV = linspace(0,k_max,nsteps)
#k_DE = linspace(0,k_max,nsteps)
## we plot k in 1/cm and omega in GHz
#plt.plot(10**-2*k_BV,10**-9*vomega(k_BV,0.0,mp),'b-',10**-2*k_DE,10**-9*vomega(k_DE,pi/2,mp),'r-')

plt.subplot(221)
plt.title('Dispersion Surface')
plt.xlabel('kpara (1/m)')
plt.ylabel('kperp (1/m)')

gridpts = 100
kpara, kperp = meshgrid(linspace(-k_max,k_max,gridpts), linspace(-k_max,k_max,gridpts))
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
vo = vomega(kzeta,phi,mp)
disp_surface_limits = [ kpara.min(), kpara.max(), kperp.min(), kperp.max() ]
plt.imshow(vo, extent = disp_surface_limits)
plt.contour(vo, extent = disp_surface_limits)

plt.subplot(222)
plt.title('Group Veloctiy')
plt.xlabel('kz (1/m)')
plt.ylabel('ky (1/m)')
plt.grid(True)

gridpts = 20
kpara, kperp = meshgrid(linspace(-k_max,k_max,gridpts), linspace(-k_max,k_max,gridpts))
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
vd_z = vdwdk(kzeta,phi,mp)[0]
vd_y = vdwdk(kzeta,phi,mp)[1]
plt.quiver(kpara, kperp, vd_z, vd_y)

plt.subplot(223)
plt.title('Slowness Surface, ' + str(round(10**-9*freq,2)) + ' GHz')
plt.xlabel('kpara (1/m)')
plt.ylabel('kperp (1/m)')

lowb = freq - epsilon
highb = freq + epsilon
iso = where(vo < lowb, 0, where(vo > highb, 0, vo))
plt.imshow(iso, extent = disp_surface_limits)

plt.subplot(224)
plt.title('Emission Pattern')
plt.imshow(abs(fft.ifftshift(fft.ifft2(iso)))**2)

plt.tight_layout()
plt.show()

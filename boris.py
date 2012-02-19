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

# k_max = no of nodes*minimum k for sample width W
W = 0.0013
k0 = pi/W
k_max = 100.0*k0

# freq, epsilon for emission pattern [GHz]
freq = 3.65
epsilon = 0.01
# GHz -> Hz
freq = freq*10**9
epsilon = epsilon*10**9

# Number of points for dispersion surface & group velocity
gridpts_ds = 100
gridpts_vg = 10

plt.subplot(221)
plt.title('Dispersion Surface')
plt.xlabel('kpara [1/m]')
plt.ylabel('kperp [1/m]')

kpara, kperp = meshgrid(linspace(-k_max,k_max,gridpts_ds), linspace(-k_max,k_max,gridpts_ds))
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
vo = vomega(kzeta,phi,mp)
disp_surface_limits = [ kpara.min(), kpara.max(), kperp.min(), kperp.max() ]
plt.imshow(vo, extent = disp_surface_limits)
plt.contour(vo, extent = disp_surface_limits)

plt.subplot(222)
plt.title('Group Veloctiy')
plt.xlabel('kz [1/m]')
plt.ylabel('ky [1/m]')
plt.grid(True)

kpara, kperp = meshgrid(linspace(-k_max,k_max,gridpts_vg), linspace(-k_max,k_max,gridpts_vg))
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
vd_z = vdwdk(kzeta,phi,mp)[0]
vd_y = vdwdk(kzeta,phi,mp)[1]
plt.quiver(kpara, kperp, vd_z, vd_y)

plt.subplot(223)
plt.title('Slowness Surface, ' + str(round(10**-9*freq,2)) + ' GHz')
plt.xlabel('kpara [1/m]')
plt.ylabel('kperp [1/m]')

lowb = freq - epsilon
highb = freq + epsilon
iso = where(vo < lowb, 0, where(vo > highb, 0, vo))
plt.imshow(iso, extent = disp_surface_limits)

plt.subplot(224)
plt.title('Emission Pattern')
plt.xlabel('para [m]')
plt.ylabel('perp [m]')
n_values_z = fft.fftfreq(len(linspace(-k_max,k_max,gridpts_ds)), 2*k_max/gridpts_ds)
emis_pat_limits = [min(n_values_z), max(n_values_z), min(n_values_z), max(n_values_z)]
plt.imshow(abs(fft.ifftshift(fft.ifft2(iso)))**2, extent = emis_pat_limits)

plt.tight_layout()
plt.show()

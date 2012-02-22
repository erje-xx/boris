#! /usr/bin/env python

import matplotlib.pyplot as plt
from numpy import pi, deg2rad, linspace, meshgrid, sqrt, arctan2, where, fft, angle
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
nodes = 100.0
k_max = nodes*k0

# freq, epsilon for emission pattern [GHz]
freq = 3.65
epsilon = 0.005
# GHz -> Hz
freq = freq*10**9
epsilon = epsilon*10**9

# Number of points for dispersion surface & group velocity
# Corresponds to sample rate in k-space, gridpts_ds
gridpts_ds = 1000
gridpts_vg = 10

#summary figure number
sfn= 1
plt.figure(sfn)
plt.subplot(121)
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

plt.subplot(122)
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

# number of emission patterns
no_of_ep = 20
ep_freqs= linspace(vo.min(),vo.max(),no_of_ep)

# emission pattern limits
n_values_z = fft.fftfreq(len(linspace(-k_max,k_max,gridpts_ds)), 2*k_max/gridpts_ds)
epl = [min(n_values_z), max(n_values_z), min(n_values_z), max(n_values_z)]
zoom = 0.1

# figure number (don't want to overwrite DS/Vg)
epfn = sfn + 1
for ef in ep_freqs:
	plt.figure(epfn)
	epfn += 1

	plt.subplot(221)
	plt.title('Slowness Surface, ' + str(round(10**-9*ef,2)) + ' GHz')
	plt.xlabel('para [m]')
	plt.ylabel('perp [m]')

	lowb = ef - epsilon
	highb = ef + epsilon
	iso = where(vo < lowb, 0, where(vo > highb, 0, vo))
	plt.imshow(iso, extent = disp_surface_limits)

	plt.subplot(222)
	plt.title('Intensity')
	plt.xlabel('para [m]')
	plt.ylabel('perp [m]')

	plt.xlim(zoom*epl[0],zoom*epl[1])
	plt.ylim(zoom*epl[2],zoom*epl[3])
	plt.imshow(abs(fft.ifftshift(fft.ifft2(iso)))**2, extent = epl)

	plt.subplot(223)
	plt.title('Amplitude')
	plt.xlabel('para [m]')
	plt.ylabel('perp [m]')

	plt.xlim(zoom*epl[0],zoom*epl[1])
	plt.ylim(zoom*epl[2],zoom*epl[3])
	plt.imshow(abs(fft.ifftshift(fft.ifft2(iso))), extent = epl)

	plt.subplot(224)
	plt.title('Phase')
	plt.xlabel('para [m]')
	plt.ylabel('perp [m]')

	plt.xlim(zoom*epl[0],zoom*epl[1])
	plt.ylim(zoom*epl[2],zoom*epl[3])
	plt.imshow(angle(fft.ifftshift(fft.ifft2(iso))), extent = epl)
	
plt.tight_layout()
plt.show()

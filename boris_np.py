#! /usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
from boris_lib import omega, vomega, dwdk, d2wdk2

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

#print('Calculating 1D dispersion...')
plt.subplot(221)
plt.title('1D Dispersion')
plt.xlabel('Wavevector (1/cm)')
plt.ylabel('Frequency (GHz)')
plt.grid(True)

nsteps = 100
k_BV = linspace(0,k_max,nsteps)
k_DE = linspace(0,k_max,nsteps)
# we plot k in 1/cm and omega in GHz
plt.plot(10**-2*k_BV,10**-9*vomega(k_BV,0.0,mp),'b-',10**-2*k_DE,10**-9*vomega(k_DE,pi/2,mp),'r-')

#print('Done calculating 1D dispersion!')

#print('Calculating group velocity vector field...')
#plt.subplot(222)
#plt.title('Group Veloctiy')
#plt.xlabel('kz (1/m)')
#plt.ylabel('ky (1/m)')
#plt.grid(True)

#kzl = []
#kyl = []
#dwdkz = []
#dwdky = []
#ssf = 20
#for i in range(0, int(k_max + 1), int(k_max/ssf)):
	#kz = k0*(2*i-k_max)
	#for l in range(0, int(k_max + 1), int(k_max/ssf)):
		#ky = k0*(2*l-k_max)
		#if kz == 0 and ky == 0:
			#continue
#
		#kzeta = math.sqrt(kz**2 + ky**2)
		#phi = math.atan2(kz,ky)
		#dwdk_temp = dwdk(kzeta, phi, Ms, H, L, gamma, alpha, W, theta, k_max, trans_n)
#
		#kzl.append(kz)
		#kyl.append(ky)
		#dwdkz.append(dwdk_temp[0])
		#dwdky.append(dwdk_temp[1])
#
#plt.quiver(kzl, kyl, dwdkz, dwdky)
#print('Done calculating group velocity vector field!')

#print('Calculating dispersion surface')
plt.subplot(222)
plt.title('Dispersion Surface')
plt.xlabel('kpara (1/m)')
plt.ylabel('kperp (1/m)')

gridpts = 100
kpara, kperp = meshgrid(linspace(-k_max,k_max,gridpts), linspace(-k_max,k_max,gridpts))
disp_surface_limits = [ kpara.min(), kpara.max(), kperp.min(), kperp.max() ]
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
vo = vomega(kzeta,phi,mp)
plt.imshow(vo, extent = disp_surface_limits)
plt.contour(vo, extent = disp_surface_limits)

#print('Minimum frequency: ' + str(10**-9*vo.min()) + ' GHz')
#print('Maximum frequency: ' + str(10**-9*vo.max()) + ' GHz')
#print('Difference: ' + str(10**-9*(vo.max() - vo.min())) + ' GHz')
#print('Done calculating dispersion surface!')

#print('Calculating slowness surface...')
plt.subplot(223)
plt.title('Slowness Surface, ' + str(round(10**-9*freq,2)) + ' GHz')
plt.xlabel('kpara (1/m)')
plt.ylabel('kperp (1/m)')

lowb = freq - epsilon
highb = freq + epsilon
iso = where(vo < lowb, 0, where(vo > highb, 0, vo))
plt.imshow(iso, extent = disp_surface_limits)
#print('Done calculating slowness surface!')

plt.subplot(224)
plt.title('Emission Pattern')
plt.imshow(abs(fft.ifftshift(fft.ifft2(iso)))**2)

print(disp_surface_limits)
plt.tight_layout()
plt.show()

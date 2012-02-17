#! /usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
from boris_lib import omega
from boris_lib import dwdk
from boris_lib import d2wdk2
import math

pi = math.pi

Ms = 139260.5752 
H = 55704.230082
L = 5.1e-06
gamma = 35185.83772
alpha = 3.1e-16
phi_selected = 0.0
theta = 90.0
W = 0.0013
k_max = 1000.0
freq = 4.2
trans_n = 0.0

# Unit conversions from set_the_boris to what the routines need
# degrees -> radians
phi_selected = phi_selected * pi/180
phi = phi_selected
# degrees -> radians
theta = theta * pi/180
# GHz -> Hz
freq = freq*10**9

print('Calculating 1D dispersion...')
k0 = 2*pi/W
k = []
w1 = []
w2 = []
ssf = 100
for i in range(0,int(k_max + 1), int(k_max/ssf)):
	kzeta = k0*i
	# Frequency written to file in GHz, wavevector in inverse centimeters
	# Step size determined by film thickness, equal to 2*pi/W
	#disp_out.write(str(kzeta*10**-2) + '\t' + str(omg_temp*10**-9) + '\n')
	k.append(10**-2*kzeta)
	phi = 0.0
	w1.append(10**-9*omega(kzeta, phi, Ms, H, L, gamma, alpha, theta, trans_n))
	phi = 90*pi/180
	w2.append(10**-9*omega(kzeta, phi, Ms, H, L, gamma, alpha, theta, trans_n))

plt.subplot(221)
plt.title('1D Dispersion')
plt.xlabel('Wavevector (1/cm)')
plt.ylabel('Frequency (GHz)')
plt.grid(True)
plt.plot(k,w1,'b-',k,w2,'r-')
plt.tight_layout()
#plt.show()
print('Done calculating 1D dispersion!')

print('Calculating group velocity vector field...')
kzl = []
kyl = []
dwdkz = []
dwdky = []
ssf = 20
for i in range(0, int(k_max + 1), int(k_max/ssf)):
	kz = k0*(2*i-k_max)
	for l in range(0, int(k_max + 1), int(k_max/ssf)):
		ky = k0*(2*l-k_max)
		if kz == 0 and ky == 0:
			continue

		kzeta = math.sqrt(kz**2 + ky**2)
		phi = math.atan2(kz,ky)
		dwdk_temp = dwdk(kzeta, phi, Ms, H, L, gamma, alpha, W, theta, k_max, trans_n)

		kzl.append(kz)
		kyl.append(ky)
		dwdkz.append(dwdk_temp[0])
		dwdky.append(dwdk_temp[1])

plt.subplot(222)
plt.title('Group Veloctiy')
plt.xlabel('kz (1/m)')
plt.ylabel('ky (1/m)')
plt.grid(True)
plt.quiver(kzl, kyl, dwdkz, dwdky)
print('Done calculating group velocity vector field!')

print('Calculating dispersion surface')
# We want to write a function that calcualtes the
# dispersion surface from -k_max to k_max in both
# direction, with a given step_size
# Step-size factor
ssf = 100
# What follows works, but abomination much?
c = zeros(( len(range(0, int(k_max + 1), int(k_max/ssf))), len(range(0, int(k_max + 1), int(k_max/ssf))) ))
d = zeros(( len(range(0, int(k_max + 1), int(k_max/ssf))), len(range(0, int(k_max + 1), int(k_max/ssf))) ))
#c[::,100] = range(0, int(k_max + 1), 10)

kperp = zeros(( len(range(0, int(k_max + 1), int(k_max/ssf))), len(range(0, int(k_max + 1), int(k_max/ssf))) ))
kpara = zeros(( len(range(0, int(k_max + 1), int(k_max/ssf))), len(range(0, int(k_max + 1), int(k_max/ssf))) ))
testest = 0
for i in range(len(kpara)):
	kpara[i] = range(int(-k_max/2), int(k_max/2 + 1), int(k_max/ssf))
for i in range(len(kperp)):
	kperp[i] = kperp[i] + k_max - testest - k_max/2
	testest += k_max/ssf
disp_surface_limits = [ kpara.min(), kpara.max(), kperp.min(), kperp.max() ]
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
for i in range(len(kzeta)):
	for j in range(len(kzeta[i])):
		c[i,j] = omega(kzeta[i,j],phi[i,j],Ms,H,L,gamma,alpha,theta,trans_n)
		#d[i,j] = linalg.det(d2wdk2(kzeta[i,j],phi[i,j],Ms,H,L,gamma,alpha, W, theta, k_max, trans_n))
plt.subplot(223)
#plt.imshow(c[::-1,:])
#plt.contour(c[::-1,:])
plt.title('Dispersion Surface')
plt.xlabel('kz (1/m)')
plt.ylabel('ky (1/m)')
plt.imshow(c[::-1,:], extent = disp_surface_limits )
plt.contour(c[::-1,:], extent = disp_surface_limits)

print('Done calculating dispersion surface!')

plt.subplot(224)
print('Calculating slowness surface...')
flip = c[::-1,:].copy()
freq = 4.0*10**9
epsilon = 0.01*10**9
for i in range(len(flip)):
	for j in range(len(flip[i])):
		if flip[i,j] < (freq - epsilon) or flip[i,j] > (freq + epsilon):
			#print(flip[i,j])
			flip[i,j] = 0

plt.imshow(abs(fft.ifftshift(fft.ifft2(flip)))**2)
print('Done calculating slowness surface!')

plt.tight_layout()
plt.show()

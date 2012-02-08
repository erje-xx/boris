#! /usr/bin/env python

from numpy import *
import matplotlib.pyplot as plt
from boris_lib import omega
from boris_lib import dwdk
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
w = []
for i in range(0,int(k_max + 1)):
	kzeta = k0*i
	# Frequency written to file in GHz, wavevector in inverse centimeters
	# Step size determined by film thickness, equal to 2*pi/W
	#disp_out.write(str(kzeta*10**-2) + '\t' + str(omg_temp*10**-9) + '\n')
	k.append(10**-2*kzeta)
	w.append(10**-9*omega(kzeta, phi, Ms, H, L, gamma, alpha, theta, trans_n))

plt.subplot(221)
plt.title('1D Dispersion')
plt.xlabel('Wavevector (1/cm)')
plt.ylabel('Frequency (GHz)')
plt.grid(True)
plt.plot(k,w)
plt.tight_layout()
#plt.show()
print('Done calculating 1D dispersion!')

print('Calculating group velocity vector field...')
kzl = []
kyl = []
dwdkz = []
dwdky = []
for i in range(0, int(k_max + 1), 10):
	kz = k0*(2*i-k_max)
	for l in range(0,int(k_max + 1), 10):
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
# 
#disp = []
#disp_line = []
#for i in range(0, int(k_max + 1), 10):
	#ky = k0*(2*i-k_max)
	#for l in range(0,int(k_max + 1), 10):
		#kz = k0*(2*l-k_max)
		#if kz == 0 and ky == 0:
			#continue
		#kzeta = math.sqrt(kz**2 + ky**2)
		#phi = math.atan2(kz,ky)
		#disp_line.append(10**-9*omega(kzeta, phi, Ms, H, L, gamma, alpha, theta, trans_n))
	#disp.insert(0, disp_line)
	#disp_line = []

# What follows works, but is an abomination
a = zeros(( len(range(0, int(k_max + 1), 10)), len(range(0, int(k_max + 1), 10)) ))
b = zeros(( len(range(0, int(k_max + 1), 10)), len(range(0, int(k_max + 1), 10)) ))
c = zeros(( len(range(0, int(k_max + 1), 10)), len(range(0, int(k_max + 1), 10)) ))
#c[::,100] = range(0, int(k_max + 1), 10)
testest = 0
for i in range(len(a)):
	a[i] = a[i] + k_max - testest
	testest += 10
	#print(rowrow)
for i in range(len(b)):
	b[i] = range(0, int(k_max + 1), 10)
#print(b)
#print(a)
kperp = a
kpara = b
#print(a[0,0])
kzeta = sqrt(kperp**2 + kpara**2)
phi = arctan2(kperp,kpara)
for i in range(len(kzeta)):
	for j in range(len(kzeta[i])):
		c[i,j] = omega(kzeta[i,j],phi[i,j],Ms,H,L,gamma,alpha,theta,trans_n)
plt.subplot(223)
plt.imshow(c[::-1,:])
print('Done calculating dispersion surface!')

plt.show()

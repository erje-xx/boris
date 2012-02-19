#! /usr/bin/env python

import matplotlib.pyplot as plt
from numpy import pi, deg2rad, linspace
from boris_lib import vomega

# material parameters
Ms = 139260.5752 
H = 55704.230082
L = 5.1e-06
gamma = 35185.83772
alpha = 3.1e-16
theta = deg2rad(90.0)
trans_n = 0.0
mp = [Ms, H, L, gamma, alpha, theta, trans_n]

phi = deg2rad(45.0)

# k_max = no of nodes*minimum k for sample width W
W = 0.0013
k0 = pi/W
k_max = 100.0*k0

plt.subplot(111)
plt.title('1D Dispersion')
plt.xlabel('Wavevector (1/cm)')
plt.ylabel('Frequency (GHz)')
plt.grid(True)

nsteps = 100
k_BV = linspace(0,k_max,nsteps)
k_DE = linspace(0,k_max,nsteps)
k_user = linspace(0,k_max,nsteps)
# we plot k in 1/cm and omega in GHz
plt.plot(10**-2*k_BV,10**-9*vomega(k_BV,0.0,mp),'b-',10**-2*k_DE,10**-9*vomega(k_DE,pi/2,mp),'r-', 10**-2*k_user,10**-9*vomega(k_user,phi,mp),'k-')

plt.tight_layout()
plt.show()

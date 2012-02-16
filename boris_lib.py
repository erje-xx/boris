import math
# frequency
def omega(kzeta, phi, Ms, H, L, gamma, alpha, theta, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms

	# Exception for FMR
	if kzeta == 0.0:
		return math.sqrt(wh*(wh + wm))
	
	kappa_n = n*pi/L
	kn = math.sqrt(kzeta**2 + kappa_n**2)
	Q = (wh + alpha*wm*kn*kn)
	
	# Fn, same for totally unpinned and totally pinned cases
	Fn = (2/(kzeta*L))*(1 - ((-1)**n)*math.exp(-kzeta*L))

	# Pnn, totally UNPINNED surface spins
	if n == 0:
		B = 0.5
	else:
		B = 1
	
	Pnn = kzeta**2/(kn**2) - Fn*B*kzeta**4/(kn**4)

	A = math.sin(theta)*math.sin(theta)
	C = math.cos(phi)*math.cos(phi)
	D = math.sin(phi)*math.sin(phi)
	Fnn_1 = Pnn
	Fnn_2 = A
	Fnn_3 = -A*Pnn*(1 + C)
	Fnn_4 = A*wm*Pnn*(1 - Pnn)*D/Q
	Fnn = (Fnn_1 + Fnn_2 + Fnn_3 + Fnn_4)

	return (math.sqrt(Q*(Q + wm*Fnn)))

# group velocity
def dwdk(kzeta, phi, Ms, H, L, gamma, alpha, W, theta, k_max, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms

	kn = math.sqrt(kzeta**2 + (n*pi/L)**2)
	dkndk = (kzeta/kn)

	Q = (wh + alpha*wm*kn*kn)
	dQdk = (2*alpha*wm*kzeta)

	Fn = ((2/(kzeta*L))*(1 - ((-1)**n)*math.exp(-kzeta*L)))

	dFndk_1 = -2/(L*kzeta*kzeta)
	dFndk_2 = 2*((-1)**n)*math.exp(-kzeta*L)/(kzeta*kzeta*L)
	dFndk_3 = 2*((-1)**n)*math.exp(-kzeta*L)/kzeta
	dFndk = (dFndk_1 + dFndk_2 + dFndk_3)

	if n==0:
		B = 0.5
	else:
		B = 1.0

	Pnn = (kzeta**2/(kn**2) - Fn*B*kzeta**4/(kn**4))

	A = math.sin(theta)*math.sin(theta)
	C = math.cos(phi)*math.cos(phi)
	D = math.sin(phi)*math.sin(phi)
	Fnn_1 = Pnn
	Fnn_2 = A
	Fnn_3 = -A*Pnn*(1 + C)
	Fnn_4 = A*wm*Pnn*(1 - Pnn)*D/Q
	Fnn = (Fnn_1 + Fnn_2 + Fnn_3 + Fnn_4)

	dPnndk_1 = 2*kzeta*(kn**-2)
	dPnndk_2 = -2*(kzeta**2)*(kn**-3)*dkndk
	dPnndk_3 = -4*(kzeta**3)*(kn**-4)*B*Fn
	dPnndk_4 = 4*(kzeta**4)*(kn**-5)*B*Fn*dkndk
	dPnndk_5 = -(kzeta**4)*(kn**-4)*B*dFndk
	dPnndk = (dPnndk_1 + dPnndk_2 + dPnndk_3 + dPnndk_4 + dPnndk_5)

	E = wm*math.sin(theta)*math.sin(theta)*math.sin(phi)*math.sin(phi)
	dFnndk_1 = dPnndk
	dFnndk_2 = -dPnndk*math.sin(theta)*math.sin(theta)*(1 + math.cos(phi)*math.cos(phi))
	dFnndk_3 = - (1/Q)*(1/Q)*E*Pnn*dQdk
	dFnndk_4 = (1/Q)*E*dPnndk
	dFnndk_5 = (1/Q)*(1/Q)*E*Pnn*Pnn*dQdk
	dFnndk_6 = -2*(1/Q)*E*Pnn*dPnndk
	dFnndk = (dFnndk_1 + dFnndk_2 + dFnndk_3 + dFnndk_4 + dFnndk_5 + dFnndk_6)

	# We use our own omega here, since it's simple and calling omega_n costs time
	omega = math.sqrt(Q*(Q + wm*Fnn))

	dbetadk = dQdk*(2*Q + wm*Fnn) + wm*Q*dFnndk
	M = wm*math.sin(theta)*math.sin(theta)*Pnn*(1 - Pnn)/Q
	
	dwdkzeta = (1/(2*omega))*dbetadk
	# note change due to henning's check, compare with previous results
	dwdphi = (1/(2*omega))*Q*wm*math.sin(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + M)

	dwdkz = dwdkzeta*(math.cos(phi)) - (1/kzeta)*dwdphi*math.sin(phi)
	dwdky = dwdkzeta*(math.sin(phi)) + (1/kzeta)*dwdphi*math.cos(phi)

	# Have to multiply by 2pi, since we work with the gyromagnetic ratio
	# in Hz/Oe...i.e. we work with normal gamma divided by 2pi
	# This 2pi yields real group velocity
	return [2*pi*dwdkz, 2*pi*dwdky, 2*pi*dwdkzeta, 2*pi*dwdphi]

# 2nd deriv
def d2wdk2(kzeta, phi, Ms, H, L, gamma, alpha, W, theta, k_max, n):
	if kzeta == 0:
		return [ [0 , 0] , [0 , 0] ]
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms

	kn = math.sqrt(kzeta**2 + (n*pi/L)**2)
	dkndk = (kzeta/kn)
	d2kndk2 = (1/kn) - (kzeta**2/kn**3)

	dphidky = math.cos(phi)/kzeta
	dphidkz = -math.sin(phi)/kzeta
	d2phidky2 = -math.sin(2*phi)/kzeta**2
	d2phidkz2 = math.sin(2*phi)/kzeta**2
	d2phidkydkz = (1 - 2*math.cos(phi)*math.cos(phi))/(kzeta**2)

	dkzetadky = math.sin(phi)
	dkzetadkz = math.cos(phi)
	d2kzetadky2 = math.cos(phi)*math.cos(phi)/kzeta
	d2kzetadkz2 = math.sin(phi)*math.sin(phi)/kzeta
	d2kzetadkydkz = -math.cos(phi)*math.sin(phi)/kzeta

	R = (wh + alpha*wm*kn*kn)
	dRdk = (2*alpha*wm*kzeta)
	d2Rdk2 = (2*alpha*wm)

	Fn = ((2/(kzeta*L))*(1 - ((-1)**n)*math.exp(-kzeta*L)))

	dFndk_1 = -2/(L*kzeta*kzeta)
	dFndk_2 = 2*((-1)**n)*math.exp(-kzeta*L)/(kzeta*kzeta*L)
	dFndk_3 = 2*((-1)**n)*math.exp(-kzeta*L)/kzeta
	dFndk = (dFndk_1 + dFndk_2 + dFndk_3)

	d2Fndk2_1 = 4*(1/(L*kzeta**3))
	d2Fndk2_2 = -4*((-1)**n)*(1/(L*kzeta**3))*math.exp(-kzeta*L) - 2*((-1)**n)*(1/kzeta**2)*math.exp(-kzeta*L)
	d2Fndk2_3 = -2*((-1)**n)*(1/kzeta**2)*math.exp(-kzeta*L) - 2*L*((-1)**n)*(1/kzeta)*math.exp(-kzeta*L)
	d2Fndk2 = (d2Fndk2_1 + d2Fndk2_2 + d2Fndk2_3)

	if n==0:
		B = 0.5
	else:
		B = 1.0

	Pnn = (kzeta**2/(kn**2) - Fn*B*kzeta**4/(kn**4))

	A = math.sin(theta)*math.sin(theta)
	C = math.cos(phi)*math.cos(phi)
	D = math.sin(phi)*math.sin(phi)
	Fnn_1 = Pnn
	Fnn_2 = A
	Fnn_3 = -A*Pnn*(1 + C)
	Fnn_4 = A*wm*Pnn*(1 - Pnn)*D/R
	Fnn = (Fnn_1 + Fnn_2 + Fnn_3 + Fnn_4)

	# We use our own omega here, since it's simple and calling omega_n costs time
	omega = math.sqrt(R*(R + wm*Fnn))

	dPnndk_1 = 2*kzeta*(kn**-2)
	dPnndk_2 = -2*(kzeta**2)*(kn**-3)*dkndk
	dPnndk_3 = -4*(kzeta**3)*(kn**-4)*B*Fn
	dPnndk_4 = 4*(kzeta**4)*(kn**-5)*B*Fn*dkndk
	dPnndk_5 = -(kzeta**4)*(kn**-4)*B*dFndk
	dPnndk = (dPnndk_1 + dPnndk_2 + dPnndk_3 + dPnndk_4 + dPnndk_5)

	d2Pnndk2_1 = 2/kn**2 - 4*(kzeta/kn**3)*dkndk
	d2Pnndk2_2 = -4*(kzeta/kn**3)*dkndk + 6*(kzeta**2/kn**4)*dkndk**2 - 2*(kzeta**2/kn**3)*d2kndk2
	d2Pnndk2_3 = -12*(kzeta**2/kn**4)*Fn*B + 16*(kzeta**3/kn**5)*Fn*B*dkndk - 4*(kzeta**3/kn**4)*B*dFndk
	#d2Pnndk2_4 = 16*(kzeta**3/kn**5)*Fn*B*dkndk - 20*(kzeta**4/kn**6)*Fn*B*dkndk**2 + 4*(kzeta**4/kn**5)*B*dkndk*dFndk + 4*(kzeta**4/kn**5)*Fn*B*d2kndk2
	d2Pnndk2_4_1 = 16*(kzeta**3/kn**5)*Fn*B*dkndk - 20*(kzeta**4/kn**6)*Fn*B*dkndk**2  
	d2Pnndk2_4_2 = 4*(kzeta**4/kn**5)*B*dkndk*dFndk + 4*(kzeta**4/kn**5)*Fn*B*d2kndk2
	d2Pnndk2_4 = d2Pnndk2_4_1 + d2Pnndk2_4_2
	d2Pnndk2_5 = -4*(kzeta**3/kn**4)*B*dFndk + 4*(kzeta**4/kn**5)*B*dFndk*dkndk - (kzeta**4/kn**4)*B*d2Fndk2
	d2Pnndk2 = (d2Pnndk2_1 + d2Pnndk2_2 + d2Pnndk2_3 + d2Pnndk2_4 + d2Pnndk2_5)

	E = wm*math.sin(theta)*math.sin(theta)*math.sin(phi)*math.sin(phi)
	dFnndk_1 = dPnndk
	dFnndk_2 = -dPnndk*math.sin(theta)*math.sin(theta)*(1 + math.cos(phi)*math.cos(phi))
	dFnndk_3 = - (1/R)*(1/R)*E*Pnn*dRdk
	dFnndk_4 = (1/R)*E*dPnndk
	dFnndk_5 = (1/R)*(1/R)*E*Pnn*Pnn*dRdk
	dFnndk_6 = -2*(1/R)*E*Pnn*dPnndk
	dFnndk = (dFnndk_1 + dFnndk_2 + dFnndk_3 + dFnndk_4 + dFnndk_5 + dFnndk_6)

	d2Fnndk2_1 = d2Pnndk2
	d2Fnndk2_2 = -d2Pnndk2*math.sin(theta)*math.sin(theta)*(1 + math.cos(phi)*math.cos(phi))
	d2Fnndk2_3 = 2*E*Pnn*dRdk**2*(1/R**3) - (1/R**2)*E*dRdk*dPnndk - (1/R**2)*E*Pnn*d2Rdk2
	d2Fnndk2_4 = (-1/R**2)*E*dPnndk*dRdk + (1/R)*E*d2Pnndk2
	d2Fnndk2_5 = -2*(1/R**3)*E*Pnn**2*dRdk**2 + 2*(1/R**2)*E*Pnn*dRdk*dPnndk + (1/R**2)*E*Pnn**2*d2Rdk2
	d2Fnndk2_6 = 2*(1/R**2)*E*Pnn*dRdk*dPnndk - 2*(1/R)*E*dPnndk**2 - 2*(1/R)*E*Pnn*d2Pnndk2
	d2Fnndk2 = (d2Fnndk2_1 + d2Fnndk2_2 + d2Fnndk2_3 + d2Fnndk2_4 + d2Fnndk2_5 + d2Fnndk2_6)

	dbetadk = dRdk*(2*R + wm*Fnn) + wm*R*dFnndk
	M = wm*math.sin(theta)*math.sin(theta)*Pnn*(1 - Pnn)/R
	
	# note change in phi due to henning's check, compare with previous results
	dwdkzeta = (1/(2*omega))*dbetadk
	dwdphi = (1/(2*omega))*R*wm*math.sin(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + M)

	# first order derivs, cartesian
	dwdkz = dwdkzeta*(math.cos(phi)) - (1/kzeta)*dwdphi*math.sin(phi)
	dwdky = dwdkzeta*(math.sin(phi)) + (1/kzeta)*dwdphi*math.cos(phi)

	# kzeta 2nd deriv
	d2w2dk2 = (2*R+wm*Fnn)*d2Rdk2 + R*wm*d2Fnndk2 + (2*dRdk + 2*wm*dFnndk)*dRdk
	d2wdkzeta2 = (1/(2*omega))*d2w2dk2 - (1/omega)*(dwdkzeta)**2

	# phi 2nd deriv
	d2Fnndphi2 = 2*math.cos(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + wm*Pnn*(1-Pnn)*math.sin(theta)*math.sin(theta)/R)
	d2w2dphi2 = R*wm*d2Fnndphi2
	d2wdphi2 = (1/(2*omega))*d2w2dphi2 - (1/omega)*(dwdphi)**2

	# mixed 2nd deriv
	T = wm*math.sin(2*phi)*math.sin(theta)*math.sin(theta)
	U = wm*wm*math.sin(2*phi)*math.sin(theta)*math.sin(theta)
	d2wdkdphi_1 = (-T/(2*omega**2))*R*Pnn*dwdkzeta + (T/(2*omega))*Pnn*dRdk + (T/(2*omega))*R*dPnndk
	d2wdkdphi_2 = (-U/(2*omega**2))*(Pnn - Pnn**2)*dwdkzeta + (U/(2*omega))*dPnndk*(1 - 2*Pnn)
	d2wdkdphi = d2wdkdphi_1 + d2wdkdphi_2

	# 2nd derivs to cartesian, conversion factors
	d2wdkzdkzeta = dkzetadkz*d2wdkzeta2 + dphidkz*d2wdkdphi
	d2wdkzdphi = dkzetadkz*d2wdkdphi + dphidkz*d2wdphi2
	d2wdkydkzeta = dkzetadky*d2wdkzeta2 + dphidky*d2wdkdphi
	d2wdkydphi = dkzetadky*d2wdkdphi + dphidky*d2wdphi2
	
	# 2nd derivs Cartesian
	d2wdkz2 = d2wdkzdkzeta*dkzetadkz + d2wdkzdphi*dphidkz + dwdkzeta*d2kzetadkz2 + dwdphi*d2phidkz2
	d2wdky2 = d2wdkydkzeta*dkzetadky + d2wdkydphi*dphidky + dwdkzeta*d2kzetadky2 + dwdphi*d2phidky2
	d2wdkydkz = d2wdkzdkzeta*dkzetadky + d2wdkzdphi*dphidky + dwdkzeta*d2kzetadkydkz + dwdphi*d2phidkydkz

	# Can be conveniently reshaped into 2x2 numpy array :)
	return [ [d2wdkz2, d2wdkydkz] , [d2wdkydkz, d2wdky2] ]

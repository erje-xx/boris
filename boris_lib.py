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
	
	# Have to multiply by 2pi, since we work with the gyromagnetic ratio
	# in Hz/Oe...i.e. we work with normal gamma divided by 2pi
	dwdkzeta = (1/(2*omega))*dbetadk
	dwdphi = (1/(2*omega))*Q*wm*math.sin(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + M)

	dwdkz = dwdkzeta*(math.cos(phi)) - (1/kzeta)*dwdphi*math.sin(phi)
	dwdky = dwdkzeta*(math.sin(phi)) + (1/kzeta)*dwdphi*math.cos(phi)

	d2wdkzeta2 = (1/(2*omega))*d2w2dkzeta2 - (1/omega)*(dwdkzeta)**2
	return [dwdkz, dwdky, dwdkzeta, dwdphi]


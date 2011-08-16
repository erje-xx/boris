import math
def sc():
	pi = 3.14159265
	L = 5.1*(10**-6)
	kzl = 0
	alpha = 3*(10**-16)

	k_z = 0
	n = 1
	delta = 0
	kappa_n = n*pi/L
	phi = 0
	M = 1750
	H = 700

	for x in range(0,11):
		phi = pi*x/20
		output = open('/home/erje/output_test_'+str(int(math.degrees(phi)))+'.dat','w')
#		output = open('/home/erje/output_test_'+str(n)+'.dat','w')
		for y in range(0,10000):
			kzl=y/100
			k_n = math.sqrt((1/L**2)*(n*n*pi*pi + kzl**2))
			if kzl == 0:
				Fn = 0.0
			else:
				Fn = (2/kzl)*(1 + math.exp(-kzl))
			Pnn = (kzl**2/(n*n*pi*pi+kzl**2))*(1 + (n*n*pi*pi/(n*n*pi*pi + kzl**2))*Fn)
			omega = math.sqrt((1.4 + alpha*k_n**2 - Pnn)*(0.4 + Pnn*math.sin(phi)*math.sin(phi) + alpha*k_n**2))
			output.write(str(kzl) + ' ' + str(omega) + '\n')
		output.close()

def sc2():
	Ms = 1750
	H = 700
	A = 3*(10**-16)
	d= 5.1*(10**-6)
	L = 0.0013
	pi = 3.14159265

	qy = 2*pi/L
	qz= 2*pi/L

	qp = math.sqrt(qz**2 + qy**2)

	P = 1 + ((qp*d)**-1)*(1-math.exp(-qp*d))
	F = 1 + P*(1-P)*( (4*pi*Ms)/(H + 2*A*qp**2/Ms) )*(qy**2/qp**2) - P*(qz**2/qp**2)
	D = H + 2*A*qp**2/Ms
	G = D + 4*pi*Ms*F

	omega_bar = math.sqrt((H + 2*A*qp**2/Ms)*(H + 2*A*qp**2/Ms + 4*pi*F))


	dqpdqy = qy*1/(math.sqrt(qy**2 + qz**2))
	dqpdqz = qz*1/(math.sqrt(qy**2 + qz**2))
	
	dPdqy = dqpdqy*(d*math.exp(-qp*d)*(qp*d)**-1 - (1-math.exp(-qp*d))*d*(qp*d)**-2)
	dPdqz = dqpdqy*(d*math.exp(-qp*d)*(qp*d)**-1 - (1-math.exp(-qp*d))*d*(qp*d)**-2)

	dFdqy_1 = dPdqy -2*P*dPdqy
	dFdqy_2 = -16*pi*A*dqpdqy*D**-2
	dFdqy_3 = -2*qy*qp**-2 - 2*qy**2*qp**-3*dqpdqy
	dFdqy_4 = qz**2*qp**-2*dPdqy -2*P*qz**-2*qp**-3*dqpdqy
	dFdqy = dFdqy_1 + dFdqy_2 + dFdqy_3 + dFdqy_4
	dDdqy = 4*A*qp*dqpdqy/Ms
	dGdqy = dDdqy + 4*pi*Ms*dFdqy
	dwdqy = 0.5*(G*D)**(-0.5)*(D*dGdqy + G*dDdqy)

	dFdqz_1 = (1 - 2*P)*dPdqz
	dFdqz_2 = -4*pi*Ms*(H + 2*A*qp**2/Ms)**(-2)*(4*A*qp*dqpdqz/Ms)
	dFdqz_3 = qy*qy*-2*dqpdqz*qp**(-3)
	dFdqz_4 = qz**2*dPdqz*qp**(-2) + P*(2*qz*qp**(-2) + -2*qz**2*qp**(-3)*dqpdqz)
	dFdqz = dFdqz_1 + dFdqz_2 + dFdqz_3 + dFdqz_4
	dDdqz = 4*A*qp*dqpdqz/Ms
	dGdqz = dDdqz + 4*pi*Ms*dFdqz
	dwdqz = 0.5*(D*G)**(-0.5)*(G*dDdqz + D*dGdqz)

	print(str(qz) + ' ' + str(dwdqz))
	print(str(qy) + ' ' + str(dwdqy))

def sc3():
	WORKING_DIRECTORY = '/home/erje/boris/'
	pi = math.pi

	# Ms in Gauss
	Ms = 1750
	# H in Oe
	H = 700
	# A in m**2
	alpha = 3*(10**-16)
	# Film thickness d in meters
	d= 5.1*(10**-6)
	# Film in-plane extent L in meters
	L = 1.3*(10**-3)
	# gamma in Hz/Oe
	gamma = 2.8 * 10**6

	wh = gamma*H
	wm = gamma*Ms

	# angle between film normal and Ms in radians
	theta = pi/2
	# in-plane angle between propagation and Ms, [0,pi/2]
	phi = 0

	# number of nodes in the transverse (parallel to film normal) direction
	n = 0

	# Transverse (along film thickness) wavevector
	kappa_n = n*pi/d

	for i in range(0,641):
		phi = i*pi/1280
		disp_out = open(WORKING_DIRECTORY + 'disp_' + str(math.degrees(phi)) + '.dat', 'w')
		#dispdiv_out = open(WORKING_DIRECTORY + 'dispdiv_' + str(math.degrees(phi)) + '.dat', 'w')
		for l in range(1,20):
			# in-plane wavevector, in direction of propagation
			kzeta = 2*pi*l/L
	
			kn = math.sqrt(kzeta**2 + kappa_n**2)
	
			# Fn, same for totally unpinned and totally pinned cases
			Fn = (2/(kzeta*d))*(1 - ((-1)**n)*math.exp(-kzeta*d))
	
			# Pnn, totally UNPINNED surface spins
			if n == 0:
				Pnn = kzeta**2/(kn**2) - (kzeta**4/(kn**4))*Fn/2
				B = 0.5
			else:
				Pnn = kzeta**2/(kn**2) - (kzeta**4/(kn**4))*Fn
				B = 1
	
			Fnn = Pnn + math.sin(theta)*math.sin(theta)*(1 - Pnn*(1 + math.cos(phi)*math.cos(phi)) + wm*(Pnn*(1 - Pnn)*math.sin(phi)*math.sin(phi))/(wh + alpha*wm*kn**2))
	
			omega_n = math.sqrt((wh + alpha*wm*kn**2)*(wh + alpha*wm*kn**2 + wm*Fnn))
	
			kzeta_reform = kzeta*10**-2
			omega_n_reform = omega_n*10**-9
	
			print('k: ' + str(kzeta_reform) + ' w: ' + str(omega_n_reform))
			#disp_out.write('k: ' + str(kzeta_reform) + ' w: ' + str(omega_n_reform) + '\n')
			disp_out.write(str(kzeta_reform) + '\t' + str(omega_n_reform) + '\n')
	
			kern = wh + alpha*wm*kn**2
			dkndk = kzeta/kn

			dFndk = (2/kzeta)*(-1)**n*math.exp(-kzeta*d) - (2/(kzeta*kzeta*d))*(1 - (-1)**n*math.exp(-kzeta*d))
			dPnndk = 2*kzeta*kn**-2 - 2*kzeta**2*kn**-3*dkndk - 4*kzeta**3*kn**-4*Fn*B + 4*kzeta**4*kn**-5*Fn*B*dkndk - kzeta**4*kn**-4*B*dFndk
			D = math.sin(theta)*math.sin(theta)*math.sin(phi)*math.sin(phi)*wm
			E = (1/kern)*(D*Pnn - 2*D*Pnn*dPnndk) - (D*Pnn - D*Pnn**2)*2*alpha*wm*kzeta/(kern**2)
			dFnndk = dPnndk - dPnndk*math.sin(theta)*math.sin(theta)*(1 - math.cos(phi)*math.cos(phi)) + E
			dbetadk = 2*alpha*wm*kzeta*(kern + wm*Fnn) + (2*alpha*wm*kzeta + wm*dFnndk)*kern
			dwdk = (1/(2*omega_n))*dbetadk
			#dispdiv_out.write(str(kzeta_reform) + '\t' + str(dwdk) + '\n')

		disp_out.close()
		#dispdiv_out.close()

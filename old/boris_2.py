import math
import os
def boris_calculates():
	#WORKING_DIRECTORY = '/home/erje/boris/'
	WORKING_DIRECTORY = os.getcwd() + '/'
	conf_file = open(WORKING_DIRECTORY + 'set_the_boris.conf', 'r')

	pi = math.pi
	
	for line in conf_file:
		form_line = line.split('=')
		if len(form_line) == 0:
			continue
		if form_line[0].strip() == "SATURATION_MAGNETIZATION":
			Ms = float(form_line[1].strip())
		elif form_line[0].strip() == "EXTERNAL_FIELD":
			H = float(form_line[1].strip())
		elif form_line[0].strip() == "FILM_THICKNESS":
			d = float(form_line[1].strip())
		elif form_line[0].strip() == "GYROMAGNETIC_RATIO":
			gamma = float(form_line[1].strip())
		elif form_line[0].strip() == "EXCHANGE_CONSTANT":
			alpha = float(form_line[1].strip())
		elif form_line[0].strip() == "PHI":
			phi_selected = float(form_line[1].strip())
		elif form_line[0].strip() == "THETA":
			theta = float(form_line[1].strip())
		elif form_line[0].strip() == "FILM_WIDTH":
			L = float(form_line[1].strip())
		elif form_line[0].strip() == "MAXIMUM_WAVEVECTOR":
			k_max = float(form_line[1].strip())
		elif form_line[0].strip() == "EXCITATION_FREQUENCY":
			freq = float(form_line[1].strip())
		elif form_line[0].strip() == "TRANSVERSE_NODES":
			trans_n = float(form_line[1].strip())
		elif form_line[0].strip() == "DISPERSION":
			dispersion_test = float(form_line[1].strip())
		elif form_line[0].strip() == "SLOWNESS_SURFACE":
			slowness_surface_test = float(form_line[1].strip())
		else:
			continue

	print('Ms:\t\t' + str(Ms))
	print('H:\t\t' + str(H))
	print('gamma:\t\t' + str(gamma))
	print('alpha:\t\t' + str(alpha))
	print('d:\t\t' + str(d))
	print('L:\t\t' + str(L))
	print('k_max:\t\t' + str(k_max))
	print('n:\t\t' + str(trans_n))
	print('phi:\t\t' + str(phi_selected))
	print('theta:\t\t' + str(theta))
	print('frequency:\t' + str(freq))

	# Convert phi to radians
	phi_selected = phi_selected * pi/180
	# Convert theta to radians
	theta = theta * pi/180
	# Convert GHz to Hz
	freq = freq*10**9

	if abs(dispersion_test) != 0:
		dispersion(Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, trans_n)

	if abs(slowness_surface_test) != 0:
		slowness_surface(Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L)

	print('That\'s all folks!')
	return

def omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms
	
	kz = kzeta*math.cos(phi)
	ky = kzeta*math.sin(phi)
	kappa_n = n*pi/d
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

	Fnn_1 = Pnn
	Fnn_2 = math.sin(theta)*math.sin(theta)
	Fnn_3 = (-1)*math.sin(theta)*math.sin(theta)*Pnn*(1 + math.cos(phi)*math.cos(phi))
	Fnn_4 = math.sin(theta)*math.sin(theta)*wm*(Pnn*(1 - Pnn)*math.sin(phi)*math.sin(phi))/(wh + alpha*wm*kn**2)
	Fnn_new = Fnn_1 + Fnn_2 + Fnn_3 + Fnn_4
	
	omega_n = math.sqrt((wh + alpha*wm*kn**2)*(wh + alpha*wm*kn**2 + wm*Fnn))
	return omega_n

def dwdk(Ms, H, d, gamma, alpha, phi, L, theta, k_max, freq, trans_n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms
	
	kz = kzeta*math.cos(phi)
	ky = kzeta*math.sin(phi)
	kappa_n = n*pi/d
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

	Fnn_1 = Pnn
	Fnn_2 = math.sin(theta)*math.sin(theta)
	Fnn_3 = (-1)*math.sin(theta)*math.sin(theta)*Pnn*(1 + math.cos(phi)*math.cos(phi))
	Fnn_4 = math.sin(theta)*math.sin(theta)*wm*(Pnn*(1 - Pnn)*math.sin(phi)*math.sin(phi))/(wh + alpha*wm*kn**2)
	Fnn_new = Fnn_1 + Fnn_2 + Fnn_3 + Fnn_4
	
	omega_n = math.sqrt((wh + alpha*wm*kn**2)*(wh + alpha*wm*kn**2 + wm*Fnn))

	kern = wh + alpha*wm*kn**2
	M = wm*math.sin(theta)*math.sin(theta)*Pnn*(1 - Pnn)/kern

	dFndk = (2/kzeta)*(-1)**n*math.exp(-kzeta*d) - (2/(kzeta*kzeta*d))*(1 - (-1)**n*math.exp(-kzeta*d))
	dPnndk = 2*kzeta*kn**-2 - 2*kzeta**2*kn**-3*dkndk - 4*kzeta**3*kn**-4*Fn*B + 4*kzeta**4*kn**-5*Fn*B*dkndk - kzeta**4*kn**-4*B*dFndk
	D = math.sin(theta)*math.sin(theta)*math.sin(phi)*math.sin(phi)*wm
	E = (1/kern)*(D*Pnn - 2*D*Pnn*dPnndk) - (D*Pnn - D*Pnn**2)*2*alpha*wm*kzeta/(kern**2)
	dFnndk = dPnndk - dPnndk*math.sin(theta)*math.sin(theta)*(1 - math.cos(phi)*math.cos(phi)) + E
	dbetadk = 2*alpha*wm*kzeta*(kern + wm*Fnn) + (2*alpha*wm*kzeta + wm*dFnndk)*kern

	# Have to multiply by 2pi, since we work with the gyromagnetic ratio
	# in Hz/Oe...i.e. we work with normal gamma divided by 2pi
	dwdkzeta = 2*pi*(1/(2*omega_n))*dbetadk
	dwdphi = 2*pi*(1/(2*omega_n))*kern*wm*math.cos(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + M)

	dwdky = dwdkzeta*(ky/kzeta) + dwdphi*(kz/kzeta**2)
	dwdkz = dwdkzeta*(kz/kzeta) + dwdphi*(-ky/kzeta**2)
	return [dwdky, dwdkz, dwdkzeta, dwdphi]

def dispersion(Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, n):
	#WORKING_DIRECTORY = '/home/erje/boris/'
	WORKING_DIRECTORY = os.getcwd() + '/'
	disp_out = open(WORKING_DIRECTORY + 'dispersion.dat', 'w')
	#dispdiv_out = open(WORKING_DIRECTORY + 'group_velocity.dat', 'w')
	dispersion_surface_out = open(WORKING_DIRECTORY + 'dispersion_surface.dat', 'w')

	print('Calculating dispersion surface...')

	pi = math.pi

	for i in range(0,181):
		phi = i*2*pi/180
		for l in range(1,int(k_max + 1)):
			# in-plane wavevector, in direction of propagation
			kzeta = 2*pi*l/L
			ky = kzeta*math.sin(phi)
			kz = kzeta*math.cos(phi)

			omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)

			# Frequency written to file in GHz, wavevector inverse meters
			dispersion_surface_out.write(str(ky) + '\t' + str(kz) + '\t' + str(omg_temp*10**-9) + '\n')

			if phi == phi_selected:
				if kzeta == 2*pi/L:
					print('Writing 1D dispersion for specified value of PHI: ' + str(math.degrees(phi)) + ' degrees')
				# Frequency written to file in GHz, wavevector in inverse centimeters
				disp_out.write(str(kzeta*10**-2) + '\t' + str(omg_temp*10**-9) + '\n')

			#dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, omega_fixed, n)
	
			#if phi == phi_selected:
			#	# Wavevector in inverse centimeters, group velocity in m/s
			#	dispdiv_out.write(str(kzeta*10**-2) + '\t' + str(dwdkzeta) + '\n')
	
		dispersion_surface_out.write('\n')

	dispersion_surface_out.close()
	disp_out.close()
	#dispdiv_out.close()
	print('Done calculating dispersion surface!')

def group_velocity(Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, n):
	#WORKING_DIRECTORY = '/home/erje/boris/'
	WORKING_DIRECTORY = os.getcwd() + '/'
	dispdiv_out = open(WORKING_DIRECTORY + 'group_velocity.dat', 'w')

	print('Calculating group velocity...')

	pi = math.pi

	for i in range(0,181):
		phi = i*2*pi/180
		for l in range(1,int(k_max + 1)):
			# in-plane wavevector, in direction of propagation
			kzeta = 2*pi*l/L
			ky = kzeta*math.sin(phi)
			kz = kzeta*math.cos(phi)

			omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)


			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, omega_fixed, n)
	
			if phi == phi_selected:
				# Wavevector in inverse centimeters, group velocity in m/s
				dispdiv_out.write(str(kzeta*10**-2) + '\t' + str(dwdkzeta) + '\n')
	
	dispdiv_out.close()
	print('Done calculating group velocity!')

def slowness_surface(Ms, H, d, gamma, alpha, theta, n, omega_fixed, k_max, L):
	#WORKING_DIRECTORY = '/home/erje/boris/'
	WORKING_DIRECTORY = os.getcwd() + '/'
	slowness_surface_out = open(WORKING_DIRECTORY + 'slowness_surface.dat', 'w')

	pi = math.pi

	epsilon = omega_fixed/10**4

	print('Calculating slowness surface at frequency: ' + str(omega_fixed*10**-9) + ' GHz')

	psi_resolution = 360
	for i in range(0, psi_resolution + 1):
		psi = i*2*pi/psi_resolution
		for l in range(1,int(k_max + 1)):
			kzeta = 2*pi*l/L
			ky = kzeta*math.sin(psi)
			kz = kzeta*math.cos(psi)

			omg_temp = omega_n(kzeta, psi, Ms, H, d, gamma, alpha, theta, n)

			if abs(omg_temp - omega_fixed) < epsilon:
				slowness_surface_out.write(str(ky) + '\t' + str(kz) + '\n')

	#ky = 0.0
	#kz = 0.01

	#while True:
		#kzeta = math.sqrt(ky**2 + kz**2)
		#if kz == 0.0 and ky > 0:
			#phi = pi/2
		#elif kz == 0.0 and ky < 0:
			#phi = 3*pi/2
		#else:
			#phi = math.atan(ky/kz)
#
		#omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)
#
		#if abs(omg_temp - omega_fixed) < epsilon:
			#ky_init = ky
			#kz_init = kz
			#print('Found initial ky, kz...starting along slowness surface...')
			#print('ky: ' + str(ky))
			#print('kz: ' + str(kz))
			#print('w: ' + str(omg_temp*10**-9) + ' GHz')
			#break
		#kz = kz + 10.0
#
	#psi = 0
	#first_pass = 0
	#delta_k = kz_init/10
	#print('delta_k: ' + str(delta_k) + ' inverse meters')
	#print('epsilon: ' + str(epsilon*10**-9) + ' GHz')
	#while True:
		#ky = ky_init + delta_k*math.cos(psi)
		#kz = kz_init + delta_k*math.sin(psi)
#
		#kzeta = math.sqrt(ky**2 + kz**2)
		#phi = math.atan(ky/kz)
#
		#omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)
#
		#if abs(omg_temp - omega_fixed) < epsilon:
			#ky_init = ky
			#kz_init = kz
			#slowness_surface_out.write(str(ky) + '\t' + str(kz) + '\n')
			#continue
		#if math.atan(ky/kz) > 3*pi/2 and first_pass == 0:
			#fourth_quadrant_flag = 1
			#first_pass = 1
			#print('Entered quadrant IV for slowness surface---there\'s no going back!')
#
		#if first_pass == 1 and math.atan(ky/kz) < pi/2:
			#print('Finished looping the slowness surface in k-space---returning to base...')
			#break
		#psi = psi + 0.001


	#for i in range(0,181):
		#phi = i*2*pi/180
		#for l in range(1,int(k_max + 1)):
			## in-plane wavevector, in direction of propagation
			#kzeta = 2*pi*l/L
			#kz = kzeta*math.cos(phi)
			#ky = kzeta*math.sin(phi)
	#
			#kn = math.sqrt(kzeta**2 + kappa_n**2)
			#dkndk = kzeta/kn
	#
			## Fn, same for totally unpinned and totally pinned cases
			#Fn = (2/(kzeta*d))*(1 - ((-1)**n)*math.exp(-kzeta*d))
	#
			## Pnn, totally UNPINNED surface spins
			#if n == 0:
				#Pnn = kzeta**2/(kn**2) - (kzeta**4/(kn**4))*Fn/2
				#B = 0.5
			#else:
				#Pnn = kzeta**2/(kn**2) - (kzeta**4/(kn**4))*Fn
				#B = 1
	#
			#Fnn = Pnn + math.sin(theta)*math.sin(theta)*(1 - Pnn*(1 + math.cos(phi)*math.cos(phi)) + wm*(Pnn*(1 - Pnn)*math.sin(phi)*math.sin(phi))/(wh + alpha*wm*kn**2))
	#
			#omega_n = math.sqrt((wh + alpha*wm*kn**2)*(wh + alpha*wm*kn**2 + wm*Fnn))
	#
			#if abs(omega_n - omega_fixed) < epsilon:
				#slowness_surface_out.write(str(ky) + '\t' + str(kz) + '\n')
				#break
#
			## Convert inverse meters to inverse centimeters
			##kzeta_reform = kzeta*10**-2
	#
			##kern = wh + alpha*wm*kn**2
			##M = wm*math.sin(theta)*math.sin(theta)*Pnn*(1 - Pnn)/kern
	#
			##dFndk = (2/kzeta)*(-1)**n*math.exp(-kzeta*d) - (2/(kzeta*kzeta*d))*(1 - (-1)**n*math.exp(-kzeta*d))
			##dPnndk = 2*kzeta*kn**-2 - 2*kzeta**2*kn**-3*dkndk - 4*kzeta**3*kn**-4*Fn*B + 4*kzeta**4*kn**-5*Fn*B*dkndk - kzeta**4*kn**-4*B*dFndk
			##D = math.sin(theta)*math.sin(theta)*math.sin(phi)*math.sin(phi)*wm
			##E = (1/kern)*(D*Pnn - 2*D*Pnn*dPnndk) - (D*Pnn - D*Pnn**2)*2*alpha*wm*kzeta/(kern**2)
			##dFnndk = dPnndk - dPnndk*math.sin(theta)*math.sin(theta)*(1 - math.cos(phi)*math.cos(phi)) + E
			##dbetadk = 2*alpha*wm*kzeta*(kern + wm*Fnn) + (2*alpha*wm*kzeta + wm*dFnndk)*kern
	#
			## Have to multiply by 2pi, since we work with the gyromagnetic ratio
			## in Hz/Oe...i.e. we work with normal gamma divided by 2pi
			##dwdkzeta = 2*pi*(1/(2*omega_n))*dbetadk
			##dwdphi = 2*pi*(1/(2*omega_n))*kern*wm*math.cos(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + M)
	#
			##dwdky = dwdkzeta*(ky/kzeta) + dwdphi*(kz/kzeta**2)
			##dwdkz = dwdkzeta*(kz/kzeta) + dwdphi*(-ky/kzeta**2)
	#
	#
	slowness_surface_out.close()
	print('Done calculating slowness surface!')

if __name__ == "__main__":
	boris_calculates()

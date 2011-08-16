import math
import os
def conf_file_readin():
	#WORKING_DIRECTORY = '/home/erje/boris/'
	WORKING_DIRECTORY = os.getcwd() + '/'
	conf_file = open(WORKING_DIRECTORY + 'set_the_boris.conf', 'r')
	
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
			phi = float(form_line[1].strip())
		elif form_line[0].strip() == "THETA":
			theta = float(form_line[1].strip())
		elif form_line[0].strip() == "FILM_WIDTH":
			L = float(form_line[1].strip())
		elif form_line[0].strip() == "NODE_NUMBER":
			node_number = float(form_line[1].strip())
		elif form_line[0].strip() == "EXCITATION_FREQUENCY":
			freq = float(form_line[1].strip())
		elif form_line[0].strip() == "TRANSVERSE_NODES":
			trans_n = float(form_line[1].strip())
		else:
			continue

	print('Ms:\t\t' + str(Ms))
	print('H:\t\t' + str(H))
	print('d:\t\t' + str(d))
	print('gamma:\t\t' + str(gamma))
	print('alpha:\t\t' + str(alpha))
	print('L:\t\t' + str(L))
	print('node_number:\t' + str(node_number))
	print('phi:\t\t' + str(phi))
	print('theta:\t\t' + str(theta))
	print('n:\t\t' + str(trans_n))
	print('frequency:\t' + str(freq))
	parameter_list = [Ms, H, d, gamma, alpha, phi, L, theta, node_number, freq, trans_n]
	return parameter_list

def dispersion_calc():
	#WORKING_DIRECTORY = '/home/erje/boris/'
	WORKING_DIRECTORY = os.getcwd() + '/'
	disp_out = open(WORKING_DIRECTORY + 'dispersion.dat', 'w')
	dispdiv_out = open(WORKING_DIRECTORY + 'group_velocity.dat', 'w')
	dispersion_surface_out = open(WORKING_DIRECTORY + 'dispersion_surface.dat', 'w')

	pi = math.pi
	parameter_list = conf_file_readin()

	Ms = parameter_list[0]
	H = parameter_list[1]
	d = parameter_list[2]
	gamma = parameter_list[3]
	alpha = parameter_list[4]
	phi_selected = parameter_list[5]
	L = parameter_list[6]
	theta = parameter_list[7]
	node_number = parameter_list[8]
	omega_fixed = parameter_list[9]
	n = parameter_list[10]

	# Ms in Gauss
	#Ms = 1750
	# H in Oe
	#H = 700
	# A in m**2
	#alpha = 3*(10**-16)
	# Film thickness d in meters
	#d= 5.1*(10**-6)
	# Film in-plane extent L in meters
	#L = 1.3*(10**-3)
	# gamma in Hz/Oe
	#gamma = 2.8 * 10**6
	# angle between film normal and Ms in radians
	#theta = pi/2
	# in-plane angle between propagation and Ms, [0,pi/2]
	#phi = 0

	# Convert phi to radians
	phi_selected = phi_selected * pi/180
	# Convert theta to radians
	theta = theta * pi/180

	wh = gamma*H
	wm = gamma*Ms

	# Transverse (along film thickness) wavevector
	kappa_n = n*pi/d

	for i in range(0,181):
		phi = i*2*pi/180
		for l in range(1,int(node_number + 1)):
			# in-plane wavevector, in direction of propagation
			kzeta = 2*pi*l/L
			kz = kzeta*math.cos(phi)
			ky = kzeta*math.sin(phi)
	
			kn = math.sqrt(kzeta**2 + kappa_n**2)
			dkndk = kzeta/kn
	
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
	
			# Convert inverse meters to inverse centimeters
			kzeta_reform = kzeta*10**-2
			# Conver Hz to GHz
			omega_n_reform = omega_n*10**-9
	
			#print('k: ' + str(kzeta_reform) + ' w: ' + str(omega_n_reform))
			#disp_out.write('k: ' + str(kzeta_reform) + ' w: ' + str(omega_n_reform) + '\n')
			if phi == phi_selected:
				disp_out.write(str(kzeta_reform) + '\t' + str(omega_n_reform) + '\n')
	
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
	
			if phi == phi_selected:
				dispdiv_out.write(str(kzeta_reform) + '\t' + str(dwdkzeta) + '\n')
			dispersion_surface_out.write(str(ky) + '\t' + str(kz) + '\t' + str(omega_n_reform) + '\n')
	
		disp_out.close()
		dispdiv_out.close()
		dispersion_surface_out.write('\n')
	dispersion_surface_out.close()

def slowness_surface():
	#WORKING_DIRECTORY = '/home/erje/boris/'
	WORKING_DIRECTORY = os.getcwd() + '/'
	slowness_surface_out = open(WORKING_DIRECTORY + 'slowness_surface.dat', 'w')

	pi = math.pi
	parameter_list = conf_file_readin()

	Ms = parameter_list[0]
	H = parameter_list[1]
	d = parameter_list[2]
	gamma = parameter_list[3]
	alpha = parameter_list[4]
	phi_selected = parameter_list[5]
	L = parameter_list[6]
	theta = parameter_list[7]
	node_number = parameter_list[8]
	omega_fixed = parameter_list[9]
	n = parameter_list[10]

	# Convert phi to radians
	phi_selected = phi_selected * pi/180
	# Convert theta to radians
	theta = theta * pi/180

	phi = phi_selected

	wh = gamma*H
	wm = gamma*Ms

	# Transverse (along film thickness) wavevector
	kappa_n = n*pi/d

	# Convert GHz to Hz
	omega_fixed = omega_fixed*10**9

	epsilon = omega_fixed/10**6
	delta_k = omega_fixed/10**6

	ky = 0
	kz = 0.01

	while True:
		kzeta = math.sqrt(ky**2 + kz**2)
		kn = math.sqrt(kzeta**2 + kappa_n**2)
		#dkndk = kzeta/kn
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

		if abs(omega_n - omega_fixed) < epsilon:
			ky_init = ky
			kz_init = kz
			break
		kz = kz + 10.0

	print('Found initial ky, kz...starting along slowness surface...')
	psi = 0
	while True:
		ky = ky_init + delta_k*math.cos(psi)
		kz = kz_init + delta_k*math.sin(psi)
		kzeta = math.sqrt(ky**2 + kz**2)
		kn = math.sqrt(kzeta**2 + kappa_n**2)
		#dkndk = kzeta/kn
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

		if abs(omega_n - omega_fixed) < epsilon:
			ky_init = ky
			kz_init = kz
			slowness_surface_out.write(str(ky) + '\t' + str(kz) + '\n')
			continue
		if psi > 0 and math.atan(ky/kz) == 0:
			break
		psi = psi + 0.01



	#for i in range(0,181):
		#phi = i*2*pi/180
		#for l in range(1,int(node_number + 1)):
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

if __name__ == "__main__":
    #dispersion_calc()
    slowness_surface()

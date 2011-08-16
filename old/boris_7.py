import math
import os
def boris_calculates():
	conf_file = open(os.getcwd() + '/' + 'set_the_boris.conf', 'r')

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
		elif form_line[0].strip() == "SLOWNESS_SURFACE_EQUAL":
			slowness_surface_equal_test = float(form_line[1].strip())
		elif form_line[0].strip() == "GROUP_VELOCITY":
			group_velocity_test = float(form_line[1].strip())
		elif form_line[0].strip() == "POINT_SOURCE":
			point_source_test = float(form_line[1].strip())
		elif form_line[0].strip() == "POINT_SOURCE_SPECTRUM":
			# We want to remove any leading "/" they might give us
			POINT_SOURCE_SPECTRUM = form_line[1].strip('./\n ')
		elif form_line[0].strip() == "WORKING_DIRECTORY":
			# There are supposed to have given us an absolute directory
			# so we assume they did that, even if they did not lead with
			# a "/". We also add the trailing "/" to the WORKING_DIRECTORY
			# to deal more directly with filenames later.
			# This processing gets us the correct directory from the following:
			# /home/erje/boris
			# /home/erje/boris/
			# home/erje/boris
			# home/erje/boris/
			# These all result in /home/erje/boris/
			WORKING_DIRECTORY_in = form_line[1].strip('/\n ')
			WORKING_DIRECTORY = '/' + WORKING_DIRECTORY_in + '/'
			
		else:
			continue

	print('WORKING_DIRECTORY: ' + WORKING_DIRECTORY)
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
		filename = WORKING_DIRECTORY + 'dispersion.dat'
		surface_filename = WORKING_DIRECTORY + 'dispersion_surface.dat'
		disp_out = open(filename, 'w')
		dispersion_surface_out = open(surface_filename, 'w')
		disp_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		disp_out.write('# Ms:\t\t' + str(Ms) + '\n')
		disp_out.write('# H:\t\t' + str(H) + '\n')
		disp_out.write('# gamma:\t' + str(gamma) + '\n')
		disp_out.write('# alpha:\t' + str(alpha) + '\n')
		disp_out.write('# d:\t\t' + str(d) + '\n')
		disp_out.write('# L:\t\t' + str(L) + '\n')
		disp_out.write('# k_max:\t' + str(k_max) + '\n')
		disp_out.write('# n:\t\t' + str(trans_n) + '\n')
		disp_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		disp_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		disp_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		disp_out.close()
		dispersion_surface_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		dispersion_surface_out.write('# Ms:\t\t' + str(Ms) + '\n')
		dispersion_surface_out.write('# H:\t\t' + str(H) + '\n')
		dispersion_surface_out.write('# gamma:\t' + str(gamma) + '\n')
		dispersion_surface_out.write('# alpha:\t' + str(alpha) + '\n')
		dispersion_surface_out.write('# d:\t\t' + str(d) + '\n')
		dispersion_surface_out.write('# L:\t\t' + str(L) + '\n')
		dispersion_surface_out.write('# k_max:\t' + str(k_max) + '\n')
		dispersion_surface_out.write('# n:\t\t' + str(trans_n) + '\n')
		dispersion_surface_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		dispersion_surface_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		dispersion_surface_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		dispersion_surface_out.close()

		dispersion(filename, surface_filename, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, trans_n)

	if abs(slowness_surface_test) != 0:
		slowness_surface_out = open(WORKING_DIRECTORY + 'slowness_surface.dat', 'w')
		slowness_surface_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		slowness_surface_out.write('# Ms:\t\t' + str(Ms) + '\n')
		slowness_surface_out.write('# H:\t\t' + str(H) + '\n')
		slowness_surface_out.write('# gamma:\t' + str(gamma) + '\n')
		slowness_surface_out.write('# alpha:\t' + str(alpha) + '\n')
		slowness_surface_out.write('# d:\t\t' + str(d) + '\n')
		slowness_surface_out.write('# L:\t\t' + str(L) + '\n')
		slowness_surface_out.write('# k_max:\t' + str(k_max) + '\n')
		slowness_surface_out.write('# n:\t\t' + str(trans_n) + '\n')
		slowness_surface_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		slowness_surface_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		slowness_surface_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		slowness_surface_out.close()

		slowness_surface(WORKING_DIRECTORY + 'slowness_surface.dat', Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L)

	if abs(slowness_surface_equal_test) != 0:
		surface_filename = WORKING_DIRECTORY + 'slowness_surface_equal.dat'
		vector_filename = WORKING_DIRECTORY + 'slowness_surface_equal_vector.dat'

		slowness_surface_equal_out = open(surface_filename, 'w')
		slowness_surface_equal_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		slowness_surface_equal_out.write('# Ms:\t\t' + str(Ms) + '\n')
		slowness_surface_equal_out.write('# H:\t\t' + str(H) + '\n')
		slowness_surface_equal_out.write('# gamma:\t' + str(gamma) + '\n')
		slowness_surface_equal_out.write('# alpha:\t' + str(alpha) + '\n')
		slowness_surface_equal_out.write('# d:\t\t' + str(d) + '\n')
		slowness_surface_equal_out.write('# L:\t\t' + str(L) + '\n')
		slowness_surface_equal_out.write('# k_max:\t' + str(k_max) + '\n')
		slowness_surface_equal_out.write('# n:\t\t' + str(trans_n) + '\n')
		slowness_surface_equal_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		slowness_surface_equal_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		slowness_surface_equal_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		slowness_surface_equal_out.close()

		slowness_surface_equal_vector_out = open(vector_filename, 'w')
		slowness_surface_equal_vector_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		slowness_surface_equal_vector_out.write('# Ms:\t\t' + str(Ms) + '\n')
		slowness_surface_equal_vector_out.write('# H:\t\t' + str(H) + '\n')
		slowness_surface_equal_vector_out.write('# gamma:\t' + str(gamma) + '\n')
		slowness_surface_equal_vector_out.write('# alpha:\t' + str(alpha) + '\n')
		slowness_surface_equal_vector_out.write('# d:\t\t' + str(d) + '\n')
		slowness_surface_equal_vector_out.write('# L:\t\t' + str(L) + '\n')
		slowness_surface_equal_vector_out.write('# k_max:\t' + str(k_max) + '\n')
		slowness_surface_equal_vector_out.write('# n:\t\t' + str(trans_n) + '\n')
		slowness_surface_equal_vector_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		slowness_surface_equal_vector_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		slowness_surface_equal_vector_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		slowness_surface_equal_vector_out.close()

		slowness_surface_equal(surface_filename, vector_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L)

	if abs(group_velocity_test) != 0:
		regular_filename = WORKING_DIRECTORY + 'group_velocity_1D.dat'
		surface_filename = WORKING_DIRECTORY + 'group_velocity.dat'
		dispdiv_out = open(regular_filename, 'w')
		dispdiv_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		dispdiv_out.write('# Ms:\t\t' + str(Ms) + '\n')
		dispdiv_out.write('# H:\t\t' + str(H) + '\n')
		dispdiv_out.write('# gamma:\t' + str(gamma) + '\n')
		dispdiv_out.write('# alpha:\t' + str(alpha) + '\n')
		dispdiv_out.write('# d:\t\t' + str(d) + '\n')
		dispdiv_out.write('# L:\t\t' + str(L) + '\n')
		dispdiv_out.write('# k_max:\t' + str(k_max) + '\n')
		dispdiv_out.write('# n:\t\t' + str(trans_n) + '\n')
		dispdiv_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		dispdiv_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		dispdiv_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		dispdiv_out.close()
		dispdiv_surface_out = open(surface_filename, 'w')
		dispdiv_surface_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		dispdiv_surface_out.write('# Ms:\t\t' + str(Ms) + '\n')
		dispdiv_surface_out.write('# H:\t\t' + str(H) + '\n')
		dispdiv_surface_out.write('# gamma:\t' + str(gamma) + '\n')
		dispdiv_surface_out.write('# alpha:\t' + str(alpha) + '\n')
		dispdiv_surface_out.write('# d:\t\t' + str(d) + '\n')
		dispdiv_surface_out.write('# L:\t\t' + str(L) + '\n')
		dispdiv_surface_out.write('# k_max:\t' + str(k_max) + '\n')
		dispdiv_surface_out.write('# n:\t\t' + str(trans_n) + '\n')
		dispdiv_surface_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		dispdiv_surface_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		dispdiv_surface_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		dispdiv_surface_out.close()

		group_velocity(regular_filename, surface_filename, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, trans_n)

	if abs(point_source_test) != 0:
		data_filename = WORKING_DIRECTORY + POINT_SOURCE_SPECTRUM
		phase_filename = WORKING_DIRECTORY + 'point_source_phase.dat'
		amplitude_filename = WORKING_DIRECTORY + 'point_source_amplitude.dat'

		amplitude_out = open(amplitude_filename, 'w')
		phase_out = open(phase_filename, 'w')

		phase_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		phase_out.write('# Ms:\t\t' + str(Ms) + '\n')
		phase_out.write('# H:\t\t' + str(H) + '\n')
		phase_out.write('# gamma:\t' + str(gamma) + '\n')
		phase_out.write('# alpha:\t' + str(alpha) + '\n')
		phase_out.write('# d:\t\t' + str(d) + '\n')
		phase_out.write('# L:\t\t' + str(L) + '\n')
		phase_out.write('# k_max:\t' + str(k_max) + '\n')
		phase_out.write('# n:\t\t' + str(trans_n) + '\n')
		phase_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		phase_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		phase_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		phase_out.close()

		amplitude_out.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
		amplitude_out.write('# Ms:\t\t' + str(Ms) + '\n')
		amplitude_out.write('# H:\t\t' + str(H) + '\n')
		amplitude_out.write('# gamma:\t' + str(gamma) + '\n')
		amplitude_out.write('# alpha:\t' + str(alpha) + '\n')
		amplitude_out.write('# d:\t\t' + str(d) + '\n')
		amplitude_out.write('# L:\t\t' + str(L) + '\n')
		amplitude_out.write('# k_max:\t' + str(k_max) + '\n')
		amplitude_out.write('# n:\t\t' + str(trans_n) + '\n')
		amplitude_out.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
		amplitude_out.write('# theta:\t' + str(math.degrees(theta)) + '\n')
		amplitude_out.write('# frequency:\t' + str(freq*10**-9) + '\n')
		amplitude_out.close()

		point_source(data_filename, phase_filename, amplitude_filename, L)

	return


def dispersion(filename, surface_filename, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, n):
	disp_out = open(filename, 'a')
	dispersion_surface_out = open(surface_filename, 'a')

	print('Calculating dispersion surface...')

	pi = math.pi
	
	phi_resolution = 5 * 360
	for i in range(0,phi_resolution + 1):
		phi = i*2*pi/(phi_resolution)
		for l in range(0,int(k_max + 1)):
			# in-plane wavevector, in direction of propagation
			kzeta = 2*pi*l/L
			ky = kzeta*math.sin(phi)
			kz = kzeta*math.cos(phi)

			omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)

			# Frequency written to file in GHz, wavevector inverse meters
			dispersion_surface_out.write(str(kz) + '\t' + str(ky) + '\t' + str(omg_temp*10**-9) + '\n')

			if phi == phi_selected:
				if kzeta == 2*pi/L:
					print('Writing 1D dispersion for specified value of PHI: ' + str(math.degrees(phi)) + ' degrees')
				# Frequency written to file in GHz, wavevector in inverse centimeters
				disp_out.write(str(kzeta*10**-2) + '\t' + str(omg_temp*10**-9) + '\n')

		# We need this newline for gnuplot's pm3d function
		dispersion_surface_out.write('\n')

	dispersion_surface_out.close()
	disp_out.close()
	print('Done calculating dispersion surface!')
	return


def group_velocity(filename, surface_filename, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, n):
	dispdiv_out = open(filename, 'a')
	dispdiv_surface_out = open(surface_filename, 'a')

	print('Calculating group velocity vector field...')

	pi = math.pi

	phi_resolution = 180
	for i in range(0,phi_resolution + 1):
		phi = i*2*pi/phi_resolution
		for l in range(1,int(k_max + 1)):
			# in-plane wavevector, in direction of propagation
			kzeta = 2*pi*l/L
			ky = kzeta*math.sin(phi)
			kz = kzeta*math.cos(phi)

			dwdk_temp = dwdk_test(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n)
			#dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n)
			#dwdk_temp = dwdk_new(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n)
	
			dispdiv_surface_out.write(str(kz) + '\t' + str(ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')

			if phi == phi_selected:
				# Wavevector in inverse centimeters, group velocity in m/s
				dispdiv_out.write(str(kzeta*10**-2) + '\t' + str(dwdk_temp[2]) + '\n')

		# We need this newline for gnuplot's pm3d function
		dispdiv_surface_out.write('\n')
	
	dispdiv_out.close()
	dispdiv_surface_out.close()
	print('Done calculating group velocity vector field!')
	return

def slowness_surface(surface_filename, Ms, H, d, gamma, alpha, theta, n, omega_fixed, k_max, L):
	slowness_surface_out = open(filename, 'a')

	pi = math.pi

	epsilon = omega_fixed/10**4

	print('Calculating slowness surface at frequency: ' + str(omega_fixed*10**-9) + ' GHz')

	psi_resolution = 10*360
	for i in range(0, int(psi_resolution/4 + 1)):
		psi = i*2*pi/(psi_resolution)
		for l in range(1,int(k_max + 1)):
			kzeta = 2*pi*l/L
			ky = kzeta*math.sin(psi)
			kz = kzeta*math.cos(psi)

			omg_temp = omega_n(kzeta, psi, Ms, H, d, gamma, alpha, theta, n)

			if abs(omg_temp - omega_fixed) < epsilon:
				slowness_surface_out.write(str(kz) + '\t' + str(ky) + '\n')
				slowness_surface_out.write(str(-kz) + '\t' + str(ky) + '\n')
				slowness_surface_out.write(str(-kz) + '\t' + str(-ky) + '\n')
				slowness_surface_out.write(str(kz) + '\t' + str(-ky) + '\n')

	slowness_surface_out.close()
	print('Done calculating slowness surface!')
	return

def slowness_surface_equal(surface_filename, vector_filename, Ms, H, d, gamma, alpha, theta, n, omega_fixed, k_max, L):
	slowness_surface_equal_out = open(surface_filename, 'a')
	slowness_surface_equal_vector_out = open(vector_filename, 'a')

	pi = math.pi

	epsilon = omega_fixed/10**4

	print('Calculating slowness surface at frequency: ' + str(omega_fixed*10**-9) + ' GHz')

	k_max = (2*pi/L)*k_max
	kz = 0.01
	ky = 0.00
	didnt_find_it = 0

	while True:
		kzeta = math.sqrt(ky**2 + kz**2)
		phi = math.atan(ky/kz)

		omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)

		if abs(omg_temp - omega_fixed) < epsilon:
			ky_init = ky
			kz_init = kz
			print('Found initial ky, kz...starting along slowness surface...')
			print('ky: ' + str(ky))
			print('kz: ' + str(kz))
			print('w: ' + str(omg_temp*10**-9) + ' GHz')
			break
		if kz > k_max:
			didnt_find_it = 1
			break
		kz = kz + 10.0

	if didnt_find_it == 1:
		kz = 0.0
		ky = 0.01
		while True:
			kzeta = math.sqrt(ky**2 + kz**2)
			if kz == 0.0:
				phi = pi/2
			else:
				phi = math.atan(ky/kz)
	
			omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)
	
			if abs(omg_temp - omega_fixed) < epsilon:
				ky_init = ky
				kz_init = kz
				print('Found initial ky, kz...starting along slowness surface...')
				print('ky: ' + str(ky))
				print('kz: ' + str(kz))
				print('w: ' + str(omg_temp*10**-9) + ' GHz')
				break
			if ky > k_max:
				didnt_find_it = 1
				print('We found nothing!!!')
				return
			ky = ky + 10.0

	if ky_init == 0.0:
		psi = 0.0
		delta_k = kz_init/100
	elif kz_init == 0.0:
		psi = 0.0
		delta_k = ky_init/100

	psi = 0.0
	psi_old = 0

	reverse_switch = 0
	print('delta_k: ' + str(delta_k) + ' inverse meters')
	print('epsilon: ' + str(epsilon*10**-9) + ' GHz')
	while True:
		ky = ky_init + delta_k*math.sin(psi)
		kz = kz_init + delta_k*math.cos(psi)

		kzeta = math.sqrt(ky**2 + kz**2)
		if kz == 0.0:
			phi = pi/2
		else:
			phi = math.atan(ky/kz)
		
		if phi < 0 or phi > pi/2:
			print('Finished sweeping first quadrant!')
			break

		omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)

		if abs(omg_temp - omega_fixed) < epsilon:
			ky_init = ky
			kz_init = kz
		
			slowness_surface_equal_out.write(str(kz) + '\t' + str(ky) + '\n')
			slowness_surface_equal_out.write(str(-kz) + '\t' + str(ky) + '\n')
			slowness_surface_equal_out.write(str(-kz) + '\t' + str(-ky) + '\n')
			slowness_surface_equal_out.write(str(kz) + '\t' + str(-ky) + '\n')

			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n)
			slowness_surface_equal_vector_out.write(str(kz) + '\t' + str(ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')

			psi_old = psi
			reverse_switch = 0
			continue
	
		if psi > (psi_old + pi/2):
			reverse_switch = 1
			psi = psi_old

		if reverse_switch == 1:
			psi = psi - math.radians(1)
		else:
			psi = psi + math.radians(1)

	slowness_surface_equal_out.close()
	slowness_surface_equal_vector_out.close()
	print('Done calculating slowness surface!')
	return

def point_source(data_filename, phase_filename, amplitude_filename, L):
	phase_out = open(phase_filename, 'a')
	amplitude_out = open(amplitude_filename,'a')
	#grid_out = open('/home/erje/boris/' + 'grid.dat')
	print('Calculating emission pattern...')

	emission_amplitude_accum = 0.0
	emission_phase_accum = 0.0

	y_steps = 100
	z_steps = 100

	for y_prim in range(0, y_steps + 1):
		slowness_surface_data = open(data_filename, 'r')
		y = L*y_prim/(y_steps) - L/2
		for z_prim in range(0, z_steps + 1):
			z = L*z_prim/z_steps - L/2
			for line in slowness_surface_data:
				if line[0] == '#':
					continue
				kz = float(line.split()[0].strip())
				ky = float(line.split()[1].strip())

				#emission_profile = complex(math.cos(y*ky + z*kz), math.sin(y*ky + z*kz))
				#emission_amplitude = math.sqrt(emission_profile.real**2 + emission_profile.imag**2)
				#emission_phase = math.atan(emssion_profile.imag/emission_profile.real)
				emission_amplitude = 1
				emission_phase = math.sqrt(y**2 + z**2)*math.sqrt(ky**2 + kz**2)
	
				emission_amplitude_accum = emission_amplitude_accum + emission_amplitude
				emission_phase_accum = emission_phase_accum + emission_phase

			phase_out.write(str(z) + '\t' + str(y) + '\t' + str(emission_phase_accum) + '\n')
			amplitude_out.write(str(z) + '\t' + str(y) + '\t' + str(emission_amplitude_accum) + '\n')
		emission_amplitude_accum = 0.0
		emission_phase_accum = 0.0
		# For gnuplot pm3d!
		phase_out.write('\n')
		amplitude_out.write('\n')
		slowness_surface_data.close()

	phase_out.close()
	amplitude_out.close()
	slowness_surface_data.close()
	print('Done calculating emission pattern!')
	return

def omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms

	# Exception for FMR
	if kzeta == 0.0:
		return math.sqrt(wh*(wh + wm))
	
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

def dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms
	
	kz = kzeta*math.cos(phi)
	ky = kzeta*math.sin(phi)
	kappa_n = n*pi/d
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

	#dwdky = dwdkzeta*(ky/kzeta) + dwdphi*(kz/kzeta**2)
	#dwdkz = dwdkzeta*(kz/kzeta) + dwdphi*(-ky/kzeta**2)
	dwdkz = dwdkzeta*(math.cos(phi)) - dwdphi*(math.sin(phi)/kzeta)
	dwdky = dwdkzeta*(math.sin(phi)) + dwdphi*(math.cos(phi)/kzeta)
	return [dwdkz, dwdky, dwdkzeta, dwdphi]

def dwdk_test(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms
	kn = math.sqrt(kzeta**2 + (n*pi/d)**2)

	#A = (wm/wh)*math.sin(phi)*math.sin(phi)
	#Foo = 1 - math.cos(phi)*math.cos(phi) - 3*Fo/2 + A*Fo*Fo/4

	omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)
	dF_temp = dFnndk(kzeta,kn,d,n,wh,wm,alpha,theta,phi)
	dQ_temp = dQdk(kzeta,wm,alpha)
	Q_temp = Q(kn,wh,wm,alpha)
	Fnn_temp = Fnn(kzeta,kn,d,n,theta,phi,alpha,wm,wh)

	dbetadk = dQ_temp*(2*Q_temp + wm*Fnn_temp) + wm*Q_temp*dF_temp
	
	# Have to multiply by 2pi, since we work with the gyromagnetic ratio
	# in Hz/Oe...i.e. we work with normal gamma divided by 2pi
	dwdkzeta = 2*pi*(1/(2*omg_temp))*dbetadk

	return [0.0, 0.0, dwdkzeta, 0.0]

def dQdk(kzeta,wm,alpha):
	return (2*alpha*wm*kzeta)

def dkndk(kzeta,kn):
	return (kzeta/kn)

def dFndk(kzeta,d,n):
	dFndk_1 = -2/(d*kzeta*kzeta)
	dFndk_2 = 2*((-1)**n)*math.exp(-kzeta*d)/(kzeta*kzeta*d)
	dFndk_3 = 2*((-1)**n)*math.exp(-kzeta*d)/kzeta
	return (dFndk_1 + dFndk_2 + dFndk_3)

def dPnndk(kzeta,kn,d,n):
	if n==0:
		B = 0.5
	else:
		B = 1.0

	dPnndk_1 = 2*kzeta*(kn**-2)
	dPnndk_2 = -2*(kzeta**2)*(kn**-3)*dkndk(kzeta,kn)
	dPnndk_3 = -4*(kzeta**3)*(kn**-4)*B*Fn(kzeta,d,n)
	dPnndk_4 = 4*(kzeta**4)*(kn**-5)*B*Fn(kzeta,d,n)*dkndk(kzeta,kn)
	dPnndk_5 = -(kzeta**4)*(kn**-4)*B*dFndk(kzeta,d,n)
	return (dPnndk_1 + dPnndk_2 + dPnndk_3 + dPnndk_4 + dPnndk_5)

def dFnndk(kzeta,kn,d,n,wh,wm,alpha,theta,phi):
	E = wm*math.sin(theta)*math.sin(theta)*math.sin(phi)*math.sin(phi)
	dFnndk_1 = dPnndk(kzeta,kn,d,n)
	dFnndk_2 = -dPnndk(kzeta,kn,d,n)*math.sin(theta)*math.sin(theta)*(1 + math.cos(phi)*math.cos(phi))
	dFnndk_3 = - (1/Q(kn,wh,wm,alpha))*(1/Q(kn,wh,wm,alpha))*E*Pnn(kzeta,kn,d,n)*dQdk(kzeta,wm,alpha)
	dFnndk_4 = (1/Q(kn,wh,wm,alpha))*E*dPnndk(kzeta,kn,d,n)
	dFnndk_5 = (1/Q(kn,wh,wm,alpha))*(1/Q(kn,wh,wm,alpha))*E*Pnn(kzeta,kn,d,n)*Pnn(kzeta,kn,d,n)*dQdk(kzeta,wm,alpha)
	dFnndk_6 = -2*(1/Q(kn,wh,wm,alpha))*E*Pnn(kzeta,kn,d,n)*dPnndk(kzeta,kn,d,n)
	return (dFnndk_1 + dFnndk_2 + dFnndk_3 + dFnndk_4 + dFnndk_5 + dFnndk_6)

def Q(kn,wh,wm,alpha):
	return (wh + alpha*wm*kn*kn)

def Fn(kzeta,d,n):
	return ((2/(kzeta*d))*(1 - ((-1)**n)*math.exp(-kzeta*d)))

def Pnn(kzeta,kn,d,n):
	if n == 0:
		B = 0.5
	else:
		B = 1.0
	return (kzeta**2/(kn**2) - Fn(kzeta,d,n)*B*kzeta**4/(kn**4))

def Fnn(kzeta,kn,d,n,theta,phi,alpha,wm,wh):
	A = math.sin(theta)*math.sin(theta)
	C = math.cos(phi)*math.cos(phi)
	D = math.sin(phi)*math.sin(phi)
	Fnn_1 = Pnn(kzeta,kn,d,n)
	Fnn_2 = A
	Fnn_3 = -A*Pnn(kzeta,kn,d,n)*(1 + C)
	Fnn_4 = A*wm*Pnn(kzeta,kn,d,n)*(1 - Pnn(kzeta,kn,d,n))*D/Q(kn,wh,wm,alpha)
	return (Fnn_1 + Fnn_2 + Fnn_3 + Fnn_4)

def dispersion_rect(filename, surface_filename, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, n):
	disp_out = open(filename, 'a')
	dispersion_surface_out = open(surface_filename, 'a')

	print('Calculating dispersion surface...')

	pi = math.pi
	
	for i in range(0,int(k_max + 1)):
		kz = i*2*pi/L
		for l in range(1,int(k_max + 1)):
			# in-plane wavevector, in direction of propagation
			ky = 2*pi*l/L
			kzeta = math.sqrt(ky**2 + kz**2)
			phi = math.atan2(ky,kz)
			omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n)

			# Frequency written to file in GHz, wavevector inverse meters
			dispersion_surface_out.write(str(kz) + '\t' + str(ky) + '\t' + str(omg_temp*10**-9) + '\n')
			dispersion_surface_out.write(str(-kz) + '\t' + str(ky) + '\t' + str(omg_temp*10**-9) + '\n')
			dispersion_surface_out.write(str(-kz) + '\t' + str(-ky) + '\t' + str(omg_temp*10**-9) + '\n')
			dispersion_surface_out.write(str(kz) + '\t' + str(-ky) + '\t' + str(omg_temp*10**-9) + '\n')

			if phi == phi_selected:
				if kzeta == 2*pi/L:
					print('Writing 1D dispersion for specified value of PHI: ' + str(math.degrees(phi)) + ' degrees')
				# Frequency written to file in GHz, wavevector in inverse centimeters
				disp_out.write(str(kzeta*10**-2) + '\t' + str(omg_temp*10**-9) + '\n')

		# We need this newline for gnuplot's pm3d function
		dispersion_surface_out.write('\n')

	dispersion_surface_out.close()
	disp_out.close()
	print('Done calculating dispersion surface!')
	return

if __name__ == "__main__":
	import sys
	if len(sys.argv) == 1:
		print('Using script version -> ' + str(sys.argv[0]))
		boris_calculates()
	else:
		# Generalize this!
		boris_calculates(sys.argv[1])

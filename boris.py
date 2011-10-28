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
			# IOW, these all result in /home/erje/boris/
			WORKING_DIRECTORY = '/' + form_line[1].strip('/\n ') + '/'
			
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

	# Unit conversions from set_the_boris to what the routines need
	# degrees -> radians
	phi_selected = phi_selected * pi/180
	# degrees -> radians
	theta = theta * pi/180
	# GHz -> Hz
	freq = freq*10**9

	if abs(dispersion_test) != 0:
		disp_filename = WORKING_DIRECTORY + 'dispersion.dat'
		disp_surface_filename = WORKING_DIRECTORY + 'dispersion_surface.dat'

		file_header_prep(WORKING_DIRECTORY, disp_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)
		file_header_prep(WORKING_DIRECTORY, disp_surface_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

		dispersion_rect(disp_filename, disp_surface_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

	if abs(group_velocity_test) != 0:
		gp_filename = WORKING_DIRECTORY + 'group_velocity_1D.dat'
		vf_filename = WORKING_DIRECTORY + 'group_velocity.dat'

		file_header_prep(WORKING_DIRECTORY, gp_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)
		file_header_prep(WORKING_DIRECTORY, vf_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

		group_velocity(gp_filename, vf_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

	if abs(slowness_surface_test) != 0:
		ss_filename = WORKING_DIRECTORY + 'slowness_surface.dat'

		file_header_prep(WORKING_DIRECTORY, ss_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

		slowness_surface(ss_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

	if abs(slowness_surface_equal_test) != 0:
		ssurface_filename = WORKING_DIRECTORY + 'slowness_surface_equal.dat'
		vector_filename = WORKING_DIRECTORY + 'slowness_surface_equal_vector.dat'

		#file_header_prep(WORKING_DIRECTORY, ssurface_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)
		file_header_prep(WORKING_DIRECTORY, vector_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

		slowness_surface_equal(ssurface_filename, vector_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

	if abs(point_source_test) != 0:
		data_filename = WORKING_DIRECTORY + POINT_SOURCE_SPECTRUM
		phase_filename = WORKING_DIRECTORY + 'point_source_phase.dat'
		amplitude_filename = WORKING_DIRECTORY + 'point_source_amplitude.dat'
		vector_filename = WORKING_DIRECTORY + 'slowness_surface_equal_vector.dat'

		file_header_prep(WORKING_DIRECTORY, amplitude_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)
		file_header_prep(WORKING_DIRECTORY, phase_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

		point_source(vector_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected)

	return


#def dispersion(filename, surface_filename, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, n):

def dispersion_polar(disp_filename, disp_surface_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	disp_out = open(disp_filename, 'a')
	dispersion_surface_out = open(disp_surface_filename, 'a')

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

def dispersion_rect(disp_filename, disp_surface_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	disp_out = open(disp_filename, 'a')
	dispersion_surface_out = open(disp_surface_filename, 'a')

	print('Calculating dispersion surface...')

	pi = math.pi
	
	for i in range(0,int(k_max + 1)):
		kz = 2*pi/L*(2*i-k_max)
		for l in range(0,int(k_max + 1)):
			ky = 2*pi/L*(2*l-k_max)

			kzeta= math.sqrt(kz**2 + ky**2)
			phi = math.atan2(kz,ky)
			omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, trans_n)

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

#def group_velocity(gp_filename, vf_filename, Ms, H, d, gamma, alpha, phi_selected, L, theta, k_max, n):
def group_velocity(gp_filename, vf_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	dispdiv_out = open(gp_filename, 'a')
	dispdiv_surface_out = open(vf_filename, 'a')

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

			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n)
	
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

#def slowness_surface(surface_filename, Ms, H, d, gamma, alpha, theta, n, omega_fixed, k_max, L):
def slowness_surface(ss_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	slowness_surface_out = open(ss_filename, 'a')

	pi = math.pi

	epsilon = freq/10**4

	print('Calculating slowness surface at frequency: ' + str(freq*10**-9) + ' GHz')

	psi_resolution = 10*360
	for i in range(0, int(psi_resolution/4 + 1)):
		psi = i*2*pi/(psi_resolution)
		for l in range(1,int(k_max + 1)):
			kzeta = 2*pi*l/L
			ky = kzeta*math.sin(psi)
			kz = kzeta*math.cos(psi)

			omg_temp = omega_n(kzeta, psi, Ms, H, d, gamma, alpha, theta, n)

			if abs(omg_temp - freq) < epsilon:
				slowness_surface_out.write(str(kz) + '\t' + str(ky) + '\n')
				slowness_surface_out.write(str(-kz) + '\t' + str(ky) + '\n')
				slowness_surface_out.write(str(-kz) + '\t' + str(-ky) + '\n')
				slowness_surface_out.write(str(kz) + '\t' + str(-ky) + '\n')

	slowness_surface_out.close()
	print('Done calculating slowness surface!')
	return

def slowness_surface_rect(ss_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	slowness_surface_out = open(ss_filename, 'a')

	pi = math.pi

	epsilon = freq/10**4

	print('Calculating slowness surface at frequency: ' + str(freq*10**-9) + ' GHz')

	for i in range(0, int(k_max + 1)):
		kz = 2*pi/L*(2*i-k_max)
		for l in range(1,int(k_max + 1)):
			ky = 2*pi/L*(2*l-k_max)
			kzeta = math.sqrt(kz**2 + ky**2)
			phi = math.atan2(kz,ky)

			omg_temp = omega_n(kzeta, psi, Ms, H, d, gamma, alpha, theta, trans_n)

			if abs(omg_temp - freq) < epsilon:
				slowness_surface_out.write(str(kz) + '\t' + str(ky) + '\n')
				slowness_surface_out.write(str(-kz) + '\t' + str(ky) + '\n')
				slowness_surface_out.write(str(-kz) + '\t' + str(-ky) + '\n')
				slowness_surface_out.write(str(kz) + '\t' + str(-ky) + '\n')

	slowness_surface_out.close()
	print('Done calculating slowness surface!')
	return

#def slowness_surface_equal(ssurface_filename, vector_filename, Ms, H, d, gamma, alpha, theta, n, omega_fixed, k_max, L):
def slowness_surface_equal(ssurface_filename, vector_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	#slowness_surface_equal_out = open(ssurface_filename, 'a')
	slowness_surface_equal_vector_out = open(vector_filename, 'a')

	pi = math.pi

	epsilon = freq/10**4

	print('Calculating slowness surface at frequency: ' + str(freq*10**-9) + ' GHz')

	k_max = (2*pi/L)*k_max
	kz = 0.01
	ky = 0.00
	didnt_find_it = 1

	while kz < k_max:
		kzeta = math.sqrt(ky**2 + kz**2)
		phi = math.atan2(kz,ky)

		omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, trans_n)

		if abs(omg_temp - freq) < epsilon:
			ky_init = ky
			kz_init = kz
			print('Found initial ky, kz...starting along slowness surface...')
			print('ky: ' + str(ky))
			print('kz: ' + str(kz))
			print('w: ' + str(omg_temp*10**-9) + ' GHz')
			didnt_find_it = 0
			break
		kz = kz + 10.0

	if didnt_find_it == 1:
		kz = 0.0
		ky = 0.01
		while True:
			kzeta = math.sqrt(ky**2 + kz**2)
			phi = math.atan2(kz,ky)
	
			omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, trans_n)
	
			if abs(omg_temp - freq) < epsilon:
				ky_init = ky
				kz_init = kz
				print('Found initial ky, kz...starting along slowness surface...')
				print('ky: ' + str(ky))
				print('kz: ' + str(kz))
				print('w: ' + str(omg_temp*10**-9) + ' GHz')
				break
			if ky > k_max:
				didnt_find_it = 1
				print('No wavevector on the slowness surface up to wavevector k_max has been found')
				return
			ky = ky + 10.0

	if ky_init == 0.0:
		delta_k = kz_init/50
	elif kz_init == 0.0:
		delta_k = ky_init/50

	print('delta_k: ' + str(delta_k) + ' inverse meters')
	print('epsilon: ' + str(epsilon*10**-9) + ' GHz')

	psi = 0.0
	psi_old = 0
	reverse_switch = 0

	while True:
		kz = kz_init + delta_k*math.cos(psi)
		ky = ky_init + delta_k*math.sin(psi)

		kzeta = math.sqrt(ky**2 + kz**2)
		phi = math.atan2(kz,ky)
		
		if phi < 0 or phi > pi/2:
			print('Finished sweeping first quadrant!')
			break

		omg_temp = omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, trans_n)

		if abs(omg_temp - freq) < epsilon:
			ky_init = ky
			kz_init = kz
		
			# For polar plot
			#slowness_surface_equal_out.write(str(kzeta) + '\t' + str(phi) + '\n')
			#slowness_surface_equal_out.write(str(kzeta) + '\t' + str(phi + pi) + '\n')
			#slowness_surface_equal_out.write(str(kzeta) + '\t' + str(pi - phi) + '\n')
			#slowness_surface_equal_out.write(str(kzeta) + '\t' + str(-phi) + '\n')

			phi = math.atan2(kz,ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			slowness_surface_equal_vector_out.write(str(kz) + '\t' + str(ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')

			phi = math.atan2(-kz,ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			slowness_surface_equal_vector_out.write(str(-kz) + '\t' + str(ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')

			phi = math.atan2(-kz,-ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			slowness_surface_equal_vector_out.write(str(-kz) + '\t' + str(-ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')

			phi = math.atan2(kz,-ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			slowness_surface_equal_vector_out.write(str(kz) + '\t' + str(-ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')

			phi = math.atan2(kz,ky)
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

	#slowness_surface_equal_out.close()
	slowness_surface_equal_vector_out.close()
	print('Done calculating slowness surface!')
	return

def point_source(ref_ssurface_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	print('Calculating emission pattern')

	#grid_out = open('/home/erje/boris/grid_out.dat','w')
	phase_out = open('/home/erje/boris/phase_out.dat','w')
	amplitude_out = open('/home/erje/boris/amplitude_out.dat','w')
	reduced_out = open('/home/erje/boris/reduced.dat','w')

	pi = math.pi
	LLG_alpha = 0.0006
	omega_h = gamma*H
	omega_m = gamma*Ms
	omega_r = LLG_alpha*(omega_h + omega_m/2)

	#Fix a way to get out delta_k!
	delta_k = 15347.2002

	z_step_size = 200
	y_step_size = 200
	new_L = L

	ssurface_in = open(ref_ssurface_filename, 'r')
	i=0
	reduced_ss = []
	for line in ssurface_in:
		if line[0] == '#':
			continue
		if i % 96 == 0.0:
			kz = float(line.split()[0].strip())
			ky = float(line.split()[1].strip())
			dwdkz = float(line.split()[2].strip())
			dwdky = float(line.split()[3].strip())
			reduced_ss.append([kz,ky,dwdkz,dwdky])
			#reduced_out.write(str(kz) + '\t' + str(ky) + '\n')

			kzeta = math.sqrt(ky**2 + kz**2)
			phi = math.atan2(kz,ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			reduced_out.write(str(kz) + '\t' + str(ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')
			reduced_ss.append([kz,ky,dwdk_temp[0],dwdk_temp[1]])

			phi = math.atan2(-kz,ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			reduced_out.write(str(-kz) + '\t' + str(ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')
			reduced_ss.append([-kz,ky,dwdk_temp[0],dwdk_temp[1]])

			phi = math.atan2(-kz,-ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			reduced_out.write(str(-kz) + '\t' + str(-ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')
			reduced_ss.append([-kz,-ky,dwdk_temp[0],dwdk_temp[1]])

			phi = math.atan2(kz,-ky)
			dwdk_temp = dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, trans_n)
			reduced_out.write(str(kz) + '\t' + str(-ky) + '\t' + str(dwdk_temp[0]) + '\t' + str(dwdk_temp[1]) + '\n')
			reduced_ss.append([kz,-ky,dwdk_temp[0],dwdk_temp[1]])
			i = i + 1
		else:
			i = i + 1
			continue
	ssurface_in.close()
	reduced_out.close()

	for in_z in range(0,int(z_step_size) + 1):
		z = 2*new_L*in_z/z_step_size - new_L
		if (100*in_z/z_step_size) % 20.0 == 0.0 and (100*in_z/z_step_size) != 0.0:
			print('Finished processing ' + str(100*in_z/z_step_size) + ' percent of grid...')
		for in_y in range(0,int(y_step_size) + 1):
			y = 2*new_L*in_y/y_step_size - new_L
			grid_point = [z,y,0.0]
			amplitude = [z,y,0.0]
			phase = [z,y,0.0]
			ssurface_in = open(ref_ssurface_filename, 'r')
			#for line in reduced_ss:
			for line in ssurface_in:
				if line[0] == '#':
					continue
				kz = float(line.split()[0].strip())
				ky = float(line.split()[1].strip())
				dwdkz = float(line.split()[2].strip())
				dwdky = float(line.split()[3].strip())

				#kz = float(line[0])
				#ky = float(line[1])
				#dwdkz = float(line[2])
				#dwdky = float(line[3])

				r_dot_vg = abs(z*dwdkz + y*dwdky)
				if r_dot_vg == 0.0:
					continue
				else:
					#grid_point[2] = grid_point[2] + (1/(2*pi))*math.exp((-(z**2+y**2)*omega_r)/r_dot_vg)*math.cos(delta_k*(kz*z+ky*y))
					amplitude[2] = amplitude[2] + (1/(2*pi))*math.exp((-(z**2+y**2)*omega_r)/(r_dot_vg))
					#phase[2] = phase[2] + math.cos(delta_k*(kz*z+ky*y))
					phase[2] = phase[2] + math.cos((kz*z+ky*y))

			ssurface_in.close()
			#grid_out.write(str(grid_point[0]) + '\t' + str(grid_point[1]) + '\t' + str(grid_point[2]) + '\n')
			amplitude_out.write(str(amplitude[0]) + '\t' + str(amplitude[1]) + '\t' + str(amplitude[2]) + '\n')
			phase_out.write(str(phase[0]) + '\t' + str(phase[1]) + '\t' + str(phase[2]) + '\n')
		#grid_out.write('\n')
		amplitude_out.write('\n')
		phase_out.write('\n')

	print('Done calculating emission pattern!')
	ssurface_in.close()
	#grid_out.close()
	phase_out.close()
	amplitude_out.close()
	return

# frequency
def omega_n(kzeta, phi, Ms, H, d, gamma, alpha, theta, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms

	# Exception for FMR
	#if kzeta == 0.0:
		#return math.sqrt(wh*(wh + wm))
	
	kappa_n = n*pi/d
	kn = math.sqrt(kzeta**2 + kappa_n**2)
	Q = (wh + alpha*wm*kn*kn)
	
	# Fn, same for totally unpinned and totally pinned cases
	Fn = (2/(kzeta*d))*(1 - ((-1)**n)*math.exp(-kzeta*d))

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
def dwdk(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms

	kn = math.sqrt(kzeta**2 + (n*pi/d)**2)
	dkndk = (kzeta/kn)

	Q = (wh + alpha*wm*kn*kn)
	dQdk = (2*alpha*wm*kzeta)

	Fn = ((2/(kzeta*d))*(1 - ((-1)**n)*math.exp(-kzeta*d)))

	dFndk_1 = -2/(d*kzeta*kzeta)
	dFndk_2 = 2*((-1)**n)*math.exp(-kzeta*d)/(kzeta*kzeta*d)
	dFndk_3 = 2*((-1)**n)*math.exp(-kzeta*d)/kzeta
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
	dwdkzeta = 2*pi*(1/(2*omega))*dbetadk
	dwdphi = 2*pi*(1/(2*omega))*Q*wm*math.cos(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + M)

	dwdkz = dwdkzeta*(math.cos(phi)) - (1/kzeta)*dwdphi*math.sin(phi)
	dwdky = dwdkzeta*(math.sin(phi)) + (1/kzeta)*dwdphi*math.cos(phi)
	return [dwdkz, dwdky, dwdkzeta, dwdphi]

# 2nd deriv
def d2wdk2(kzeta, phi, Ms, H, d, gamma, alpha, L, theta, k_max, n):
	pi = math.pi
	wh = gamma*H
	wm = gamma*Ms

	kn = math.sqrt(kzeta**2 + (n*pi/d)**2)
	dkndk = (kzeta/kn)

	Q = (wh + alpha*wm*kn*kn)
	dQdk = (2*alpha*wm*kzeta)

	Fn = ((2/(kzeta*d))*(1 - ((-1)**n)*math.exp(-kzeta*d)))

	dFndk_1 = -2/(d*kzeta*kzeta)
	dFndk_2 = 2*((-1)**n)*math.exp(-kzeta*d)/(kzeta*kzeta*d)
	dFndk_3 = 2*((-1)**n)*math.exp(-kzeta*d)/kzeta
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
	dwdkzeta = 2*pi*(1/(2*omega))*dbetadk
	dwdphi = 2*pi*(1/(2*omega))*Q*wm*math.cos(2*phi)*(Pnn*math.sin(theta)*math.sin(theta) + M)

	dwdkz = dwdkzeta*(math.cos(phi)) - (1/kzeta)*dwdphi*math.sin(phi)
	dwdky = dwdkzeta*(math.sin(phi)) + (1/kzeta)*dwdphi*math.cos(phi)
	return [dwdkz, dwdky, dwdkzeta, dwdphi]

def file_header_prep(WORKING_DIRECTORY, prep_filename, Ms, H, d, gamma, alpha, theta, trans_n, freq, k_max, L, phi_selected):
	file_to_be_prep = open(prep_filename, 'w')
	file_to_be_prep.write('# WORKING_DIRECTORY: ' + WORKING_DIRECTORY + '\n')
	file_to_be_prep.write('# Ms:\t\t' + str(Ms) + '\n')
	file_to_be_prep.write('# H:\t\t' + str(H) + '\n')
	file_to_be_prep.write('# gamma:\t' + str(gamma) + '\n')
	file_to_be_prep.write('# alpha:\t' + str(alpha) + '\n')
	file_to_be_prep.write('# d:\t\t' + str(d) + '\n')
	file_to_be_prep.write('# L:\t\t' + str(L) + '\n')
	file_to_be_prep.write('# k_max:\t' + str(k_max) + '\n')
	file_to_be_prep.write('# n:\t\t' + str(trans_n) + '\n')
	file_to_be_prep.write('# phi:\t\t' + str(math.degrees(phi_selected)) + '\n')
	file_to_be_prep.write('# theta:\t' + str(math.degrees(theta)) + '\n')
	file_to_be_prep.write('# frequency:\t' + str(freq*10**-9) + '\n')
	file_to_be_prep.close()

if __name__ == "__main__":
	import sys
	if len(sys.argv) == 1:
		print('Using script version -> ' + str(sys.argv[0]))
		boris_calculates()
	else:
		# Fix this!
		boris_calculates(sys.argv[1])

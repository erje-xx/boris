import os
def readin():
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
		elif form_line[0].strip() == "ANGLE":
			phi = float(form_line[1].strip())
		elif form_line[0].strip() == "FILM_WIDTH":
			# Following three lines could be used to handle TWO 
			# equal sines as separators...but isn't worth extra code
			#if form_line[1] == '':
				#L = float(form_line[2].strip())
				#continue
			L = float(form_line[1].strip())
		else:
			continue

	print('Ms:\t' + str(Ms))
	print('H:\t' + str(H))
	print('d:\t' + str(d))
	print('gamma:\t' + str(gamma))
	print('alpha:\t' + str(alpha))
	print('phi:\t' + str(phi))
	print('L:\t' + str(L))

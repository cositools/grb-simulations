def adjust_times(input_file, output_file, start_time):
	'''
	Adjust timestamps of .tra or .sim file.

	Parameters
	----------
	input_file : pathlib.PosixPath
		Path to .sim or .tra file
	output_file : pathlib.PosixPath
		Path to .sim or .tra file to save with adjusted times
	start_time : float
		Start time of adjusted simulation in s
	'''

	initial_start_time = None

	if input_file.suffix == '.tra' or input_file.suffix == '.sim':

		with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:

			for line in in_file:
				line = line.rstrip('\n')

				if 'TI' in line:
					if initial_start_time is None:
						initial_start_time = float(line[3:18])

					new_time = float(line[3:18]) + start_time - initial_start_time
					out_file.write(f'TI {new_time:.15g}\n')

				elif 'TE' in line:
					if initial_start_time is None:
						initial_start_time = float(line[3:18])

					new_time = float(line[3:18]) + start_time - initial_start_time
					out_file.write(f'TE {new_time:.15g}\n')

				else:
					out_file.write(line + '\n')

	else:
		
		raise RuntimeError("File must be .sim or .tra.")

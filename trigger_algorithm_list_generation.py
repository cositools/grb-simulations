import ROOT as root
import numpy as np
import os
import argparse
import yaml
import random
from contextlib import contextmanager,redirect_stderr,redirect_stdout

# Suppress MEGAlib output
class suppress_stdout_stderr(object):
    
  def __init__(self):
    # Open a pair of null files
    self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
    # Save the actual stdout (1) and stderr (2) file descriptors.
    self.save_fds = [os.dup(1), os.dup(2)]

  def __enter__(self):
    # Assign the null pointers to stdout and stderr.
    os.dup2(self.null_fds[0],1)
    os.dup2(self.null_fds[1],2)

  def __exit__(self, *_):
    # Re-assign the real stdout/stderr back to (1) and (2)
    os.dup2(self.save_fds[0],1)
    os.dup2(self.save_fds[1],2)
    # Close all file descriptors
    for fd in self.null_fds + self.save_fds:
      os.close(fd)

# Initialize parser & read arguments from command line
def parse_args():

  parser = argparse.ArgumentParser()
  parser.add_argument("-y", "--yaml", help = "Path to input .yaml file", default='examples/example_trigger_algorithm.yaml')
  args = parser.parse_args()

  return args

# Read yaml file
def read_yaml(file):

  with open(file, "r") as myfile:
    inputs = yaml.safe_load(myfile)

  return inputs

# Load and initialize MEGAlib
def init_megalib():

  with suppress_stdout_stderr():
    root.gSystem.Load("$(MEGALIB)/lib/libMEGAlib.so")
    g = root.MGlobal()
    g.Initialize()

  return g

# Define paths to source files, background file, and output files
def define_paths(inputs):

  source_path = os.path.abspath(inputs['source_path']) + '/'
  if inputs['background_type'] == 'random':
    background_path = os.path.abspath(inputs['background_path']) + '/'
  elif inputs['background_type'] == 'file':
    background_path = os.path.abspath(inputs['background_path'])
  else:
    raise RuntimeError("background_type in yaml file must be either 'random' or 'file'")
  output_path = os.path.abspath(inputs['output_path']) + '/'
  if not os.path.isdir(output_path):
    os.makedirs(output_path)

  return source_path, background_path, output_path

# Load geometry
def load_geo(geometry_path):
  
  with suppress_stdout_stderr():
    geometry = root.MDGeometryQuest()
  if geometry.ScanSetupFile(root.MString(geometry_path)) == True:
    print("Geometry " + geometry_path + " loaded!")
  else:
    raise RuntimeError('Unable to load geometry ' + geometry_path)

  with suppress_stdout_stderr():
    reader = root.MFileEventsSim(geometry)

  return reader

# Open file
def open_file(reader, file_path):

  if reader.Open(root.MString(file_path)) == False:
    raise RuntimeError('Unable to open file ' + file_path)

# Open source file
def open_source_file(reader, file_path, output_dir):

  if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
  open_file(reader, file_path)

# Initialize variables
def init_variables(mm_version):

  bgo_num = 4
  ged_num = 3
  bgo_elim = 80

  if mm_version == 8:
    bgob_pos = [-19.2, 19.2, -17.5, 17.5, 14.6, 16.5]
    bgox1_pos = [17.5, 19.3, -17.5, 17.5, 16.9, 34.7]
    bgox2_pos = [-19.3, -17.5, -17.5, 17.5, 16.9, 34.7]
    bgoy1_pos = [-17.1, 17.1, 15.8, 17.6, 16.9, 34.7]
    bgoy2_pos = [-17.1, 17.1, -17.6, -15.8, 16.9, 34.7]
    ged1_z = [23.1, 24.5]
    ged2_z = [25.6, 27.0]
    ged3_z = [28.2, 29.6]
    ged4_z = [30.7, 32.1]
  elif mm_version == 12:
    bgob_pos = [-21.1, 21.1, -18.3, 18.3, 13.6, 16.]
    bgox1_pos = [19., 21.2, -18.3, 18.3, 16.3, 35.3]
    bgox2_pos = [-21.2, -19., -18.3, 18.3, 16.3, 35.3]
    bgoy1_pos = [-18.6, 18.6, 16.2, 18.4, 16.3, 35.3]
    bgoy2_pos = [-18.6, 18.6, -18.4, -16.2, 16.3, 35.3]
    ged1_z = [22.5, 23.9]
    ged2_z = [25.1, 26.5]
    ged3_z = [27.6, 29.0]
    ged4_z = [30.2, 31.6]
  else:
    raise RuntimeError('Mass model version ' + str(mm_version) + ' is not supported')

  bgo_pos = [bgob_pos, bgox1_pos, bgox2_pos, bgoy1_pos, bgoy2_pos]
  ged_pos = [ged1_z, ged2_z, ged3_z, ged4_z]

  detector_keys = ['ged1', 'ged2', 'ged3', 'ged4', 'bgox1', 'bgox2', 'bgoy1', 'bgoy2', 'bgob']

  return bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, detector_keys

# Create and fill directories for each hit
# Note: If background simulations are longer than a day, this may not work correctly
def make_hit_dict(reader, bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, detector_keys, end_time=None):
  
  times = {}
  energies = {}

  for item in detector_keys:
    times[item] = []
    energies[item] = []

  hrs_init = 19
  reach_end = False
  x=1

  with suppress_stdout_stderr():
    while True: 
      event = reader.GetNextEvent()
      if not event or reach_end:
        break
      root.SetOwnership(event, True)

      for i in range(event.GetNHTs()):
        hit = event.GetHTAt(i)
        position = hit.GetPosition()
        energy = hit.GetEnergy()
        hour = int(event.GetTime().GetHours() - hrs_init)
        hour = int(hour + 24) if hour < 0 else hour
        nanosecond = str(event.GetTime().GetNanoSeconds()).zfill(9)
        seconds = str((hour * 3600) + (int(event.GetTime().GetMinutes()) * 60) + int(event.GetTime().GetSeconds()))
        time = float(seconds + '.' + nanosecond)
  
        if not end_time is None:
          if time > end_time:
            reach_end = True
            break
      
        if hit.GetDetector() == bgo_num and energy >= bgo_elim:
          if bgo_pos[0][0] <= position[0] <= bgo_pos[0][1] and bgo_pos[0][2] <= position[1] <= bgo_pos[0][3] and bgo_pos[0][4] <= position[2] <= bgo_pos[0][5]:
            times['bgob'].append(time)
            energies['bgob'].append(energy)
          elif bgo_pos[1][0] <= position[0] <= bgo_pos[1][1] and bgo_pos[1][2] <= position[1] <= bgo_pos[1][3] and bgo_pos[1][4] <= position[2] <= bgo_pos[1][5]:
            times['bgox1'].append(time)
            energies['bgox1'].append(energy)
          elif bgo_pos[2][0] <= position[0] <= bgo_pos[2][1] and bgo_pos[2][2] <= position[1] <= bgo_pos[2][3] and bgo_pos[2][4] <= position[2] <= bgo_pos[2][5]:
            times['bgox2'].append(time)
            energies['bgox2'].append(energy)
          elif bgo_pos[3][0] <= position[0] <= bgo_pos[3][1] and bgo_pos[3][2] <= position[1] <= bgo_pos[3][3] and bgo_pos[3][4] <= position[2] <= bgo_pos[3][5]:
            times['bgoy1'].append(time)
            energies['bgoy1'].append(energy)
          elif bgo_pos[4][0] <= position[0] <= bgo_pos[4][1] and bgo_pos[4][2] <= position[1] <= bgo_pos[4][3] and bgo_pos[4][4] <= position[2] <= bgo_pos[4][5]:
            times['bgoy2'].append(time)
            energies['bgoy2'].append(energy)
        elif hit.GetDetector() == ged_num:
          if ged_pos[0][0] <= position[2] <= ged_pos[0][1]:
            times['ged1'].append(time)
            energies['ged1'].append(energy)
          elif ged_pos[1][0] <= position[2] <= ged_pos[1][1]:
            times['ged2'].append(time)
            energies['ged2'].append(energy)
          elif ged_pos[2][0] <= position[2] <= ged_pos[2][1]:
            times['ged3'].append(time)
            energies['ged3'].append(energy)
          elif ged_pos[3][0] <= position[2] <= ged_pos[3][1]:
            times['ged4'].append(time)
            energies['ged4'].append(energy)
        
  return times, energies

# Randomly choose background files in specified directory
def choose_background(num, bkg_length, bkg_interval, background_path, components, file_type):

  file_num = str(random.randint(1, num))

  if file_type == 'sequential':
    start_time = (int(file_num) - 1) * bkg_length
    background_start = random.randint(start_time, start_time + bkg_length - bkg_interval)
    background_end = background_start + bkg_interval
  elif file_type == 'simultaneous':
    background_start = random.randint(0, bkg_length - bkg_interval)
    background_end = background_start + bkg_interval
  else:
    raise RuntimeError("background_file_type in yaml file must be 'sequantial' or 'simultaneous'")

  background_paths = []

  for file in os.listdir(background_path):
    for item in components:
      if file.startswith(item + '.') and 'inc' + file_num + '.' in file and (file.split('.')[-1] == 'sim' or (file.split('.')[-1] == 'gz' and file.split('.')[-2] == 'sim')):
        background_paths.append(background_path + os.fsdecode(file))

  return background_paths, background_start, background_end

# Fill hit directories for each background component
def read_background(reader, path, bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, end_time, detector_keys):

  open_file(reader, path)
  print('Reading background file: ' + path.split('/')[-1])
  times, energies = make_hit_dict(reader, bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, detector_keys, end_time)

  return times, energies

# Read in background events
def select_background(reader, background_paths, bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, start_time, end_time, detector_keys):

  times = {}
  energies = {}

  for item in detector_keys:
    times[item] = []
    energies[item] = []

  for i in range(len(background_paths)):
    component_times, component_energies = read_background(reader, background_paths[i], bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, end_time, detector_keys)
    for key in component_times.keys():
      for j in range(len(component_times[key])):
        if component_times[key][j] > start_time:
          times[key].append(component_times[key][j])
          energies[key].append(component_energies[key][j])

  return times, energies

# Combine source and background by placing source in the middle of background
def combine(source_times, source_energies, background_times, background_energies, detector_keys):

  min_times = []
  max_times = []
  times = {}
  energies = {}
  times_sorted = {}
  energies_sorted = {}

  for item in detector_keys:
    times[item] = []
    energies[item] = []
    times_sorted[item] = []
    energies_sorted[item] = []

  for key in background_times.keys():
    try:
      min_times.append(min(background_times[key]))
      max_times.append(max(background_times[key]))
    except:
      print(key, background_times[key])

  for key in source_times.keys():
    for i in range(len(source_times[key])):
      source_times[key][i] += (min(min_times) + max(max_times)) / 2

  for key in times.keys():
    for i in range(len(source_times[key])):
      times[key].append(source_times[key][i])
      energies[key].append(source_energies[key][i])
    for i in range(len(background_times[key])):
      times[key].append(background_times[key][i])
      energies[key].append(background_energies[key][i])
    times_sorted[key], energies_sorted[key] = (list(x) for x in zip(*sorted(zip(times[key], energies[key]))))

  return times_sorted, energies_sorted

# Write times and energies to file
def write_events(file_path, times, energies):

  print('Writing to file')
  with open(file_path, 'w') as f:
    f.write('time (s)           energy (keV)\n')
    for i in range(len(times)):
      f.write("{:.9f}".format(times[i]) + '        ' + str(energies[i]) + '\n')


def main():
  g = init_megalib()
  input_file = parse_args().yaml
  inputs = read_yaml(input_file)
  source_path, background_path, output_path = define_paths(inputs)
  reader = load_geo(inputs['geometry_path'])
  bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, detector_keys = init_variables(inputs['mass_model_version'])

  for file in os.listdir(source_path):
    filename = os.fsdecode(file)
    source_name = file.split('.')[0]
    if not os.path.isdir(output_path + source_name):
      open_source_file(reader, source_path + filename, output_path + source_name + '/')
      print('Reading source file: ' + filename)
      source_times, source_energies = make_hit_dict(reader, bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, detector_keys)

      if inputs['background_type'] == 'random':
        background_paths, start, end = choose_background(inputs['background_number'], inputs['background_file_length'], inputs['background_time'], background_path, inputs['background_components'], inputs['background_file_type'])
        background_times, background_energies = select_background(reader, background_paths, bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, start, end, detector_keys)
      else:
        open_file(reader, background_path)
        print('Reading background file: ' + background_path.split('/')[-1])
        background_times, background_energies = make_hit_dict(reader, bgo_num, ged_num, bgo_elim, bgo_pos, ged_pos, detector_keys)

      times, energies = combine(source_times, source_energies, background_times, background_energies, detector_keys)

      for key in times.keys():
        write_events(output_path + source_name + '/' + key + '.txt', times[key], energies[key])
  

if __name__ == "__main__":
    main()



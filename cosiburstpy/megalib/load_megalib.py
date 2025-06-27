import ROOT as root
from cosiburstpy.utility.utility import SuppressOutput

class LoadMEGAlib():

	def __init__(self, mass_model):
		'''
		Load MEGAlib.

		Parameters
		----------
		mass_model : pathlib.PosixPath
			Path to mass model
		'''

		self.mass_model = str(mass_model)

		self.initialize_megalib()
		self.load_mass_model()

	def initialize_megalib(self):
		'''
		Load and initialize MEGAlib.
		'''

		with SuppressOutput():

			root.gSystem.Load('$(MEGALIB)/lib/libMEGAlib.so')
			g = root.MGlobal()
			g.Initialize()

	def load_mass_model(self):
		'''
		Load mass model.
		'''

		with SuppressOutput():

			mass_model = root.MDGeometryQuest()

			if mass_model.ScanSetupFile(root.MString(self.mass_model)) == True:
				print(f"Mass model {self.mass_model} loaded.")
			else:
				raise RuntimeError(f"Unable to load {self.mass_model}.")

			self.reader = root.MFileEventsSim(mass_model)

	def open_file(self, file):
		'''
		Open MEGAlib file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to file
		'''

		if self.reader.Open(root.MString(str(file))) == False:
			raise RuntimeError(f"Unable to open {file}.")
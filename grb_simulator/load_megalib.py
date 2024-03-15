import ROOT as root
from .config import suppress_output

class load_megalib():

	def __init__(self, geometry_path):
		"""
		Load MEGAlib.

		Parameters
		----------
		geometry_path : str
			Path to mass model
		"""

		self.geometry_path = geometry_path

		self.init_megalib()
		self.load_geo()

	def init_megalib(self):
		"""
		Load and initialize MEGAlib.
		"""

		with suppress_output():
			root.gSystem.Load("$(MEGALIB)/lib/libMEGAlib.so")
			g = root.MGlobal()
			g.Initialize()

	def load_geo(self):
		"""
		Load mass model.
		"""

		with suppress_output():
			geometry = root.MDGeometryQuest()

			if geometry.ScanSetupFile(root.MString(self.geometry_path)) == True:
				print("Geometry " + self.geometry_path + " loaded!")
			else:
				raise RuntimeError('Unable to load geometry ' + self.geometry_path)

			with suppress_output():
				self.reader = root.MFileEventsSim(geometry)

	def open_file(self, file_path):
		"""
		Open MEGAlib file.

		Parameters
		----------
		file_path : str
			Path to file
		"""

		if self.reader.Open(root.MString(file_path)) == False:
			raise RuntimeError('Unable to open file ' + file_path)
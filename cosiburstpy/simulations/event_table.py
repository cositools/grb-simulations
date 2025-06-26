from cosiburstpy.utility.utility import read_csv, write_csv

class EventTable():

	def __init__(self, data):
		'''
		Event table handling.

		Parameters
		----------
		data : dict
			Event table data where keys are column names and entries are lists of values
		'''

		self.data = data

	@classmethod
	def from_file(cls, file, delimiter='\t'):
		'''
		Read in event table.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to event table file
		delimiter : str, optional
			Delimiter for file

		Returns
		-------
		lightcurve : cosiburstpy.simulations.event_list.EventTable
			Event table
		'''

		data = read_csv(file, delimiter=delimiter)

		event_list = cls(data)

		return event_list

	def write_file(self, file, delimiter='\t'):
		'''
		Create event table file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to event table file
		delimiter : str, optional
			Delimiter for file
		'''

		write_csv(file, self.data, delimiter=delimiter)

	def add_columns(self, file, data, delimiter='\t'):
		'''
		Add columns to event table file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to updated event table file
		data : dict
			Event table data for additional columns
		delimiter : str, optional
			Delimiter for file
		'''

		for column in data:
			self.data[column] = data[column]

		self.write_file(file, delimiter)

import os
from gbm.finder import BurstCatalog, TriggerFtp
from .config import fill_dict, write_yaml, define_paths

class download_gbm_data():

	def __init__(self, inputs):
		"""
		Download data from GBM.

		Parameters
		----------
		inputs : dict
			Contents of input yaml file 
		"""

		self.output_path = define_paths([inputs['output_path']], [True])[0]

		if 'filters' in inputs.keys():
			self.filters = inputs['filters']
		else:
			self.filters = None

		self.download = inputs['download']

		self.slice_catalog()

	def slice_catalog(self):
		"""
		Slice GBM event catalog based on filters.
		"""

		filter_list = []

		if self.filters == None:
			print('No filters defined. Downloading all GBM triggers')
			self.sliced_event_catalog = BurstCatalog()
		else:
			for item in self.filters:
				filter_list.append((item, self.filters[item][0], self.filters[item][1]))
			
			self.sliced_event_catalog = BurstCatalog().slices(filter_list)

	def download_events(self):
		"""
		Download TTE data from GBM.

		Parameters
		----------
		download : list
			Names of columns to download from burst catalog
		output_path : str
			Path to directory to store downloaded files
		"""

		event_list = self.sliced_event_catalog.get_table(columns=self.download)
		trig_finder = TriggerFtp(event_list[len(event_list)-1][self.download.index('trigger_name')].replace('bn', ''))

		for i in range(len(event_list)):
			name = event_list[i][self.download.index('trigger_name')]
			print('Downloading data for ' + name + ' (' + str(i+1) + '/' + str(len(event_list)) + ')')
			if not os.path.isdir(self.output_path + name):
				os.mkdir(self.output_path + name)
			mydict = fill_dict(self.download, event_list[i])
			write_yaml(self.output_path + event_list[i][self.download.index('trigger_name')] + '/' + event_list[i][self.download.index('trigger_name')] + '.yaml', mydict)
			detectors = event_list[i][self.download.index('bcat_detector_mask')]
			trig_finder.set_trigger(name.replace('bn', ''))
		
			detector_list = []
			for j in range(len(detectors)):
				if detectors[j] == '1':
					if j <= 9:
						det = 'n' + str(j)
					elif j == 10:
						det = 'na'
					elif j == 11:
						det = 'nb'
					detector_list.append(det)
			bgo1 = 0
			bgo2 = 0
			for item in detectors[0:6:1]:
				bgo1 += float(item)
			for item in detectors[6:13:1]:
				bgo2 += float(item)
			if bgo1 > bgo2:
				detector_list.append('b0')
			elif bgo1 < bgo2:
				detector_list.append('b1')
			else:
				detector_list.append('b0')
				detector_list.append('b1')
			trig_finder.get_tte(self.output_path + name, dets=detector_list)


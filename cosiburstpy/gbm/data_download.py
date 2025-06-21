from gbm.finder import BurstCatalog, TriggerFtp
import logging
from cosiburstpy.utility.utility import write_yaml, make_dict

logger = logging.getLogger(__name__)

class DataDownload():

	def __init__(self, output_dir, columns, filters=None):
		'''
		Download data from GBM.

		Parameters
		----------
		output_dir : pathlib.PosixPath
			Path to directory to save output files
		'''

		self.output = output_dir
		self.output.mkdir(parents=True, exist_ok=True)

		self.download = columns

		if not 'trigger_name' in self.download:
			self.download.append('trigger_name')

		if filters:

			self.filters = filters
			filter_list = []

			for item in self.filters:
				filter_list.append((item, self.filters[item][0], self.filters[item][1]))

			self.catalog = BurstCatalog().slices(filter_list)

		else:

			self.catalog = BurstCatalog()

	def download_burst_catalog(self):
		'''
		Download burst catalog data.
		'''

		event_list = self.catalog.get_table(columns=self.download)

		for i, event in enumerate(event_list):

			name = event[self.download.index('trigger_name')]
			logger.info(f'Downloading burst catalog data for {name} ({i+1}/{len(event_list)})')

			(self.output / name).mkdir(parents=True, exist_ok=True)

			event_data = make_dict(self.download, event)
			write_yaml(self.output / name / f'{name}.yaml', event_data)

	def download_tte_data(self):
		'''
		Download TTE data.
		'''

		if not 'bcat_detector_mask' in self.download:
			self.download.append('bcat_detector_mask')

		event_list = self.catalog.get_table(columns=self.download)
		trigger_finder = TriggerFtp(event_list[len(event_list)-1][self.download.index('trigger_name')].replace('bn', ''))

		for i, event in enumerate(event_list):

			name = event[self.download.index('trigger_name')]
			logger.info(f'Downloading TTE data for {name} ({i+1}/{len(event_list)})')

			(self.output / name).mkdir(parents=True, exist_ok=True)

			detectors = event[self.download.index('bcat_detector_mask')]
			trigger_finder.set_trigger(name.replace('bn', ''))
		
			detector_list = []

			for j in range(len(detectors)):

				if detectors[j] == '1':

					if j <= 9:
						detector_list.append('n' + str(j))

					elif j == 10:
						detector_list.append('na')

					elif j == 11:
						detector_list.append('nb')

			bgo1 = 0
			for item in detectors[0:6]:
				bgo1 += float(item)

			bgo2 = 0
			for item in detectors[6:13]:
				bgo2 += float(item)

			if bgo1 > bgo2:
				detector_list.append('b0')

			elif bgo1 < bgo2:
				detector_list.append('b1')

			else:
				detector_list.append('b0')
				detector_list.append('b1')

			trigger_finder.get_tte(self.output / name, dets=detector_list)
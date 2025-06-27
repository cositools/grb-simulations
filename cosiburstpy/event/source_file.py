import astropy.units as u
import logging

logger = logging.getLogger(__name__)

class SourceFile():

	def __init__(self, event, mass_model, lightcurve, orientation=None):
		'''
		MEGAlib .source file handling.

		Parameters
		----------
		event : cosiburstpy.event.event.Event
			Source model
		mass_model : pathlib.PosixPath
			Path to simulation mass model file
		lightcurve : pathlib.PosixPath
			Path to lightcurve file
		orientation : pathlib.PosixPath, optional
			Path to orientation file
		'''

		self.name = event.name

		self.spectrum = event.spectrum.generate_text(event.energy_range)
		self.flux = event.average_photon_flux().to(u.cm**-2 * u.s**-1).value

		self.end_time = event.time_range[1].to(u.s).value + 10.

		if hasattr(event, 'polarization'):

			if orientation:
				self.polarization = event.polarization
			else:
				raise RuntimeError("Polarization is not supported in spacecraft frame.")

		position = event.position
		self.occulted = event.occulted

		if orientation:

			self.orientation = str(orientation)

			if position.frame.name == 'spacecraftframe':
				position = position.transform_to('galactic')
			
			self.l = position.l.degree
			self.b = position.b.degree

		else:

			if position.frame.name != 'spacecraftframe':
				position = position.transform_to(SpacecraftFrame(attitude=event.orientation.attitudes[0]))

			self.z = 90. - position.lat.degree
			self.a = position.lon.degree

		self.mass_model = str(mass_model)
		self.lightcurve = str(lightcurve)

	def write_file(self, file, save_acs_hits=False, earth_occultation=False, run_name='TransientSim'):
		'''
		Create .source file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .source file
		save_acs_hits : bool, optional
			Whether to save hits in anti-coincidence shields in .sim file
		earth_occultation : bool, optional
			Whether to simulate Earth occultation
		run_name : str, optional
			Name of run in .source file
		'''

		with open(file, 'w') as f:

			f.write('# Global Parameters\n')
			f.write('Version                      1\n')
			f.write(f'Geometry                     {self.mass_model}\n')

			f.write('\n# Physics list\n')  
			f.write('PhysicsListEM                LivermorePol\n')

			f.write('\n# Output formats\n')
			f.write('StoreSimulationInfo          init-only\n')

			if save_acs_hits == True:
				f.write('\n# Store shield counts\n')
				f.write('PreTriggerMode               EveryEventWithHits\n')

			f.write('\n# Run and source parameters\n')
			f.write(f'Run                          {run_name}\n')

			f.write(f'{run_name}.FileName        {self.name}\n')
			f.write(f'{run_name}.Time            {self.end_time}\n')

			if hasattr(self, 'orientation'):
				f.write(f'{run_name}.OrientationSky  Galactic File NoLoop {self.orientation}\n')

			f.write(f'{run_name}.Source          {self.name}\n')
			f.write(f'{self.name}.ParticleType     1\n')

			if hasattr(self, 'orientation'):

				f.write(f'{self.name}.Beam             FarFieldPointSource 0 0\n')
				f.write('\n# Orientation\n')
				f.write(f'{self.name}.Orientation      Galactic Fixed {self.b:.3f} {self.l:.3f}\n')

			else:

				f.write(f'{self.name}.Beam             FarFieldPointSource {self.z} {self.a}\n')

			f.write('\n# Spectrum \n')
			f.write(f'{self.name}.Spectrum         {self.spectrum}\n')

			f.write('\n# Average photon flux in photon/cm2/s\n')
			f.write(f'{self.name}.Flux             {self.flux:.5f}\n')

			if hasattr(self, 'polarization'):
				f.write('\n# Polarization \n')
				f.write(f'{self.name}.Polarization     galactic {self.polarization[0]:.2f} {self.polarization[1].degree:.2f}\n')

			f.write('\n# Lightcurve\n')
			f.write(f'{self.name}.Lightcurve       File false {self.lightcurve}\n')

			if earth_occultation:

				f.write('\n# Earth occultation\n')
				f.write(f'{self.name}.EarthOccultation false')

			elif self.occulted:

				logger.warning(f'{self.name} is at an occulted position, but Earth occultation is not being included in simulation.')

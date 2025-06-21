import astropy.units as u
import numpy as np

def mono(e, energy):
	'''
	Monoenergetic function.

	Parameters
	----------
	e : astropy.units.quantity.Quantity
		Energy
	ebreak : astropy.units.quantity.Quantity
		Energy of delta function

	Returns
	-------
	amplitude : astropy.units.quantity.Quantity
		Amplitude of spectrum at given energy
	'''

	if e == energy:
		amplitude = 1. / e.unit
	else:
		amplitude = 0. / e.unit

	return amplitude

def band(e, alpha, beta, ebreak, piv=100.*u.keV):
	'''
	Band function.

	Parameters
	----------
	e : astropy.units.quantity.Quantity
		Energy
	alpha : float
		Low energy spectral index
	beta : float
		High energy spectral index
	ebreak : astropy.units.quantity.Quantity
		Break energy
	piv : astropy.units.quantity.Quantity, optional
		Pivot energy

	Returns
	-------
	amplitude : astropy.units.quantity.Quantity
		Amplitude of spectrum at given energy
	'''

	if e <= (alpha - beta) * ebreak:
		amplitude = (e / piv)**alpha * np.exp(-e / ebreak) / e.unit
	else:
		amplitude = (e / piv)**beta * np.exp(beta - alpha) * ((alpha - beta) * ebreak / piv)**(alpha - beta) / e.unit
	
	return amplitude

def comp(e, index, epeak, piv=100.*u.keV):
	'''
	Comptonized function.

	Parameters
	----------
	e : astropy.units.quantity.Quantity
		Energy
	index : float
		Spectral index
	epeak : astropy.units.quantity.Quantity
		Peak energy
	piv : astropy.units.quantity.Quantity, optional
		Pivot energy
	
	Returns
	-------
	amplitude : astropy.units.quantity.Quantity
		Amplitude of spectrum at given energy
	'''

	amplitude = (e / piv)**index * np.exp(-(index + 2) * e / epeak) / e.unit

	return amplitude

def pl(e, index, piv=100.*u.keV):
	'''
	Power law function.

	Parameters
	----------
	e : astropy.units.quantity.Quantity
		Energy
	index : float
		Spectral index
	piv : astropy.units.quantity.Quantity, optional
		Pivot energy

	Returns
	-------
	amplitude : astropy.units.quantity.Quantity
		Amplitude of spectrum at given energy
	'''

	amplitude = (e / piv)**(-index) / e.unit

	return amplitude

def bpl(e, ebreak, index_lo, index_hi, e_max, piv=100.*u.keV):
	'''
	Broken power law function.

	Parameters
	----------
	e : astropy.units.quantity.Quantity
		Energy
	ebreak : astropy.units.quantity.Quantity
		Break energy
	index_lo : float
		Low energy spectral index
	index_hi : float
		High energy spectral index
	emax : astropy.units.quantity.Quantity
		Maximum energy
	piv : astropy.units.quantity.Quantity, optional
		Pivot energy

	Returns
	-------
	amplitude : astropy.units.quantity.Quantity
		Amplitude of spectrum at given energy
	'''

	if e <= ebreak:
		amplitude = (e / piv)**(-index_lo) / e.unit
	else:
		amplitude = (e / piv)**(-index_hi) * emax**(index_hi - index_lo) / e.unit

	return amplitude

def sbpl(e, ebreak, index_lo, index_hi, bscale, piv=100.*u.keV):
	'''
	Smoothly broken power law function.

	Parameters
	----------
	e : astropy.units.quantity.Quantity
		Energy
	ebreak : astropy.units.quantity.Quantity
		Break energy
	index_lo : float
		Low energy spectral index
	index_hi : float
		High energy spectral index
	bscale : float
		Break scale
	piv : astropy.units.quantity.Quantity, optional
		Pivot energy

	Returns
	-------
	amplitude : astropy.units.quantity.Quantity
		Amplitude of spectrum at given energy
	'''

	q = np.log10(e / ebreak) / bscale
	q_piv = np.log10(piv / ebreak) / bscale

	m = (index_hi - index_lo) / 2
	b = (index_hi + index_lo) / 2

	a = m * bscale * np.log((np.exp(q) + np.exp(-q)) / 2)
	a_piv = m * bscale * np.log((np.exp(q_piv) + np.exp(-q_piv)) / 2)

	amplitude = (e / piv)**b * 10**(a - a_piv)

	return amplitude

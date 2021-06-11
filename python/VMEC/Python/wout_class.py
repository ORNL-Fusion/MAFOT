# wout_class.py
# Author: A. Wingen
#	but based on wout_c by Aaron Sontag and IDL scripts by John Canik
# Date: Feb. 25. 2013
# ----------------------------------------------------------------------------------------
# Reads wout file form VMEC
#	wout file content stored in Member self.data, keys are names in wout-Namelist
# Provides interpolation functions for all 1-D and 2-D variables with respect to toroidal flux s
#	interpolation functions stored in Member self.spl, same keys as variable-names in self.data
# Evaluates all spectral variables, see Member-Functions self.ev() and others
# Provides all VMEC datasets for plotScripts
#
# Input:
#	Filename (string): wout Filename, netCDF format
# Output:
#	see Member functions
# typical calls:
#	W = wout_class.Wout(Filename)
#	w = W.getData()
#		w (dictionary): wout-data; a real copy of the class Member, not a reference
#	 or
#	w = W.getAll()
#		w (dictionary): wout-data + spline-functions (keys start with 'spl'); a real copy
#
# ----------------------------------------------------------------------------------------
from numpy import *
from netCDF4 import Dataset
import scipy.interpolate as scinter
import scipy.integrate as integ

from Misc.deriv import deriv

#### List of local imports in certain subroutines ####
# from warnings import filterwarnings
# import multiprocessing as mp
# import os
# import pyvisfile.silo as silo
# import VMEC.Python.vmecdb as vdb

# ----------------------------------------------------------------------------------------
# --- Wout Class -------------------------------------------------------------------------
class Wout:

	def __init__(self, Filename = None, shot = None, timeID = None, runnum = 0):
		
		# --- Member Variables ---
		self.data = {}
		self.spl = {}
		
		# --- load wout file into self.data using either file or MDS+ database---
		self.load_data(Filename, shot, timeID, runnum)
			
		# --- load spline functions into self.spl ---
		# 2-D: s mesh
		if self.data['lasym']:	# asymmetric (not stellarator symmetric)
			LIST = ['rmns', 'rmnc', 'zmns', 'zmnc', 'bsubsmns', 'bsubsmnc']
		else:					# up-down symmetric (stellarator symmetric)
			LIST = ['rmnc', 'zmns', 'bsubsmns']

		for key in LIST:
			if self.data.has_key(key):
				mnmax = self.data[key].shape[1]
				self.spl[key] = [0] * mnmax	
				for i in xrange(mnmax):
					self.spl[key][i] = scinter.UnivariateSpline(self.data['s'], self.data[key][:,i], s = 0)
				
		# 2-D: s-ds mesh (so called half-mesh)
		if self.data['lasym']:	# asymmetric (not stellarator symmetric)
			LIST = ['bmns', 'bmnc', 'gmns', 'gmnc', 'lmns', 'lmnc',
					'bsupumns', 'bsupumnc', 'bsupvmns', 'bsupvmnc',
					'bsubumns', 'bsubumnc', 'bsubvmns', 'bsubvmnc',
					'drdsmnc','dzdsmns','drdsmns','dzdsmnc']
		else:					# up-down symmetric (stellarator symmetric)
			LIST = ['bmnc', 'gmnc', 'lmns',
					'bsupumnc', 'bsupvmnc',
					'bsubumnc', 'bsubvmnc',
					'drdsmnc','dzdsmns']

		shalf = append(self.data['shalf'][0],self.data['shalf'])	# insert another grid point at s = 0
		shalf[1] = 0	
		for key in LIST:
			if self.data.has_key(key):
				mnmax = self.data[key].shape[1]
				if self.use_nyq & (not (key in ['lmns','lmnc','drdsmnc','dzdsmns','drdsmns','dzdsmnc'])): xm = 'xm_nyq'
				else: xm = 'xm'
				self.spl[key] = [0] * mnmax
				self.data[key][0,0:mnmax] = self.data[key][1,0:mnmax] * (-1)**self.data[xm] # expand beyond magnetic axis: shalf[0] = -shalf[1] => Amn(-s)*sin(mu-nv) = Amn(s)*sin(m(u+pi)-nv) = (-1)**m * Amn(s)*sin(mu-nv); same for cos
				for i in xrange(mnmax):
					#self.spl[key][i] = scinter.UnivariateSpline(self.data['shalf'], self.data[key][:,i], s = 0)
					dat = append(self.data[key][0,i],self.data[key][:,i]) # insert another data point at s = 0
					if self.data[xm][i] > 0: dat[1] = 0
					else: dat[1] = dat[2] - 0.5*shalf[2] * (dat[3] - dat[2])/(shalf[3] - shalf[2])  # extrapolate as y(s) = as**2 + c 0.5*(dat[2] + dat[0])
					self.spl[key][i] = scinter.UnivariateSpline(shalf, dat, s = 0)

		# 1-D: s mesh
		LIST = ['equif', 'specw', 'bdotgradv',
				'jcuru', 'jcurv', 'jdotb',
				'presf', 'iotaf', 'q_factor', 'qf',
				'chi', 'chipf',
				'phi', 'phipf']
		for key in LIST:
			if self.data.has_key(key):
				self.spl[key] = scinter.UnivariateSpline(self.data['s'], self.data[key], s = 0)

		# 1-D: s-ds mesh (half mesh)
		LIST = ['beta_vol', 'over_r', 'vp', 'mass',
				'DShear', 'DCurr', 'DGeod', 'DMerc', 'DWell',
				'pres', 'iotas', 'qs',
				'buco', 'bvco', 
				'phips']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key][0] = self.data[key][1]	# profiles are axisymmetric; shalf[0] = -shalf[1]
				self.spl[key] = scinter.UnivariateSpline(self.data['shalf'], self.data[key], s = 0)
				
		# identify which machine
		self.data['machine'] = 'unknown'
		Raxis, Zaxis = self.get_axis(0, n0only = True)
		Rlcfs, Zlcfs = self.get_surface()
		if (Raxis < 1.0) & (Rlcfs.min() > 0.43) & (Rlcfs.max() < 0.92) & (Zlcfs.min() > -0.45) & (Zlcfs.max() < 0.45):
			self.data['machine'] = 'cmod'
		elif (Raxis > 1.0) & (Rlcfs.min() > 1.0) & (Rlcfs.max() < 2.21) & (Zlcfs.min() > -1.06) & (Zlcfs.max() < 1.2) & ('asdex' in Filename): 
			self.data['machine'] = 'asdex'
		elif (Raxis > 1.0) & (Rlcfs.min() > 1.0) & (Rlcfs.max() < 2.4) & (Zlcfs.min() > -1.4) & (Zlcfs.max() < 1.4): 
			self.data['machine'] = 'd3d'
		elif (Raxis > 4.0) & (Rlcfs.min() > 4.25) & (Rlcfs.max() < 6.5) & (Zlcfs.min() > -1.3) & (Zlcfs.max() < 1.3): 
			self.data['machine'] = 'w7x'
		elif (Raxis > 4.0) & (Rlcfs.min() > 4.0) & (Rlcfs.max() < 8.5) & (Zlcfs.min() > -4.8) & (Zlcfs.max() < 4.8): 
			self.data['machine'] = 'iter'
			
	# ------------------------------------------------------------------------------------
	# --- Member Functions ---------------------------------------------------------------

	# --- openFile(Filename) -------------------------------------------------------------
	# open wout file and parse its contents into self.data dictionary
	def openFile(self, Filename):
		from warnings import filterwarnings
		filterwarnings('ignore')
		data = Dataset(Filename, 'r')
		varnames = data.variables.keys()
		for n in varnames:
			if(data.variables[n].size > 1):										# arrays (all types)
				self.data[n] = array(data.variables[n][:])
				if(data.variables[n].dtype == 'S1'):							# character array
					if(data.variables[n].size == data.variables[n].shape[0]):	# 1-D						
						self.data[n] = ''.join(data.variables[n][:]).strip()
					else:														# 2-D
						self.data[n] = []
						for i in xrange(data.variables[n].shape[0]):
							self.data[n].append(''.join(data.variables[n][i,:]).strip())
			else:																# single variable
				if(data.variables[n].dtype == 'float'):							# float
					self.data[n] = float64(data.variables[n][:])
				elif(data.variables[n].dtype == 'int32'):						# int
					self.data[n] = int(data.variables[n][:])
				elif(data.variables[n].dtype == 'S1'):							# char
					try: self.data[n] = ''.join(data.variables[n][:])	# at fixed bndy: mgrid_mode exists but is masked -> error
					except: self.data[n] = data.variables[n][:]
				else:															# unknown
					print 'unknown datatype in variable:', n
					self.data[n] = data.variables[n][:]
	
	
	# --- load_data(Filename) -----------------------------------------------------------------
	# open wout file and parse its contents into data dictionary
	# adds some usefull additions to data
	def load_data(self, Filename, shot, timeID, runnum):
		if Filename != None:
			self.openFile(Filename)
		else:
			self.loadMDSplus(shot, timeID, runnum)
			
		if self.data.has_key('lrfp__logical__'): 
			self.data['lrfp'] = bool(self.data['lrfp__logical__'])

		if self.data['lrfp']:
			dchids = self.data['chi'][-1]
			s = abs(self.data['chi']/dchids)
			self.data['dchids'] = dchids
		else:
			dphids = self.data['phi'][-1]
			s = abs(self.data['phi']/dphids)
			self.data['dphids'] = dphids
			if self.data.has_key('chi'):
				dchids = self.data['chi'][-1]
				psi = abs(self.data['chi']/dchids)
				self.data['dchids'] = dchids
			else:
				dchids = dphids
				psi = s
				self.data['dchids'] = dphids
		
		# some variables can have spectrum up to Nyquist frequency instead
		# only rmn, zmn and lmn are always on truncated spectrum
		# all others can be either on truncated or Nyquist spectrum
		self.use_nyq = False
		if self.data.has_key('mnmax_nyq'):
			if self.data['gmnc'].shape[1] > self.data['rmnc'].shape[1]:
				self.use_nyq = True
		
		# add self-defined keys
		self.data['s'] = s
		self.data['shalf'] = self.data['s'] - 0.5*append(self.data['s'][1] - self.data['s'][0], diff(self.data['s'])) #s - 0.5*(s[1] - s[0])
		self.data['psi'] = psi
		self.data['psihalf'] = self.data['psi'] - 0.5*append(self.data['psi'][1] - self.data['psi'][0], diff(self.data['psi'])) #s - 0.5*(s[1] - s[0])
		if self.data.has_key('lasym__logical__'): 
			self.data['lasym'] = bool(self.data['lasym__logical__'])
		if self.data.has_key('lrecon__logical__'): 
			self.data['lrecon'] = bool(self.data['lrecon__logical__'])
		if self.data.has_key('lfreeb__logical__'): 
			self.data['lfreeb'] = bool(self.data['lfreeb__logical__'])
		if isinstance(self.data['raxis_cc'], ndarray) & (self.data['ntor'] > 0):
			if self.data['lasym']:
				self.data['delta_h'] = sqrt(self.data['raxis_cc'][1]**2 + self.data['raxis_cs'][1]**2 	# helical core (n=1) axis displacement in meters
							+ self.data['zaxis_cc'][1]**2 + self.data['zaxis_cs'][1]**2)/sqrt(2.0)
			else:
				self.data['delta_h'] = sqrt(self.data['raxis_cc'][1]**2 + self.data['zaxis_cs'][1]**2)/sqrt(2.0)
		else:
			self.data['delta_h'] = 0
		self.data['qf'] = 1.0/self.data['iotaf']
		self.data['qs'] = 1.0/self.data['iotas']
		
		# add derivatives drds and dzds
		self.make_RZderiv()
						
		# set minor Radius, major Radius, correct aspect ratio
		# 'aspect' in wout is Rmajor_p/aminor_p, which are different from the ones here		
		self.data['Rmajor'] = 0.5*(self.data['rmax_surf'] + self.data['rmin_surf'])	# major radius
		self.data['aminor'] = 0.5*(self.data['rmax_surf'] - self.data['rmin_surf']) # minor radius
		self.data['aspect'] = self.data['Rmajor']/self.data['aminor']  				# correct aspect ratio
		
		# calculate beta normal
		self.data['betaN'] = self.data['betatotal']*self.data['aminor']*abs(self.data['b0']/self.data['ctor'])*1e+6
			
		# compensate for APHI on the variables, which depend on s instead of phi (found by checking each, there may be more...)
		norm_prof_s = self.data['phipf']/self.data['phi'][-1]
		norm_func = scinter.UnivariateSpline(s, norm_prof_s, s = 0)
		norm_prof_shalf = append(norm_func(self.data['shalf'][1]), norm_func(self.data['shalf'][1::]))
		
		# spectral variables s mesh
		if self.data['lasym']:	# asymmetric (not stellarator symmetric)
			LIST = ['bsubsmns', 'bsubsmnc']
		else:					# up-down symmetric (stellarator symmetric)
			LIST = ['bsubsmns']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key + '_raw'] = self.data[key].copy()
				for i in xrange(self.data[key].shape[1]):
					self.data[key][:,i] /= norm_prof_s

		# spectral variables s-ds mesh
		if self.data['lasym']:	# asymmetric (not stellarator symmetric)
			LIST = ['gmns', 'gmnc']
		else:					# up-down symmetric (stellarator symmetric)
			LIST = ['gmnc']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key + '_raw'] = self.data[key].copy()
				for i in xrange(self.data[key].shape[1]):
					self.data[key][:,i] /= norm_prof_shalf

		# profile variables s mesh
		LIST = ['equif', 'jcuru', 'jcurv', 'chipf', 'phipf']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key + '_raw'] = self.data[key].copy()
				self.data[key] /= norm_prof_s

		# profile variables s-ds mesh
		LIST = ['vp', 'phips']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key + '_raw'] = self.data[key].copy()
				self.data[key] /= norm_prof_shalf
				
		# make phips constant
		self.data['phips'] = self.data['phips'][1::].mean() + zeros(self.data['phips'].shape)
		
		# get 3D Plasma Volume, lives on s, needs aphi corrected 'vp'
		self.data['V'] = append(0, append(integ.cumtrapz(self.data['vp'][1::],self.data['shalf'][1::])*2*pi*2*pi, self.data['volume_p'])) 
		
		# profile variables s-ds mesh with funky second and last point (first point is at -ds/2 anyway)
		# use linear extrapolation
		LIST = ['DShear', 'DCurr', 'DGeod', 'DMerc', 'DWell']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key][1]  = ((self.data[key][3] - self.data[key][2])/(self.data['shalf'][3] - self.data['shalf'][2])*self.data['shalf'][1] 
									 + (self.data[key][2]*self.data['shalf'][3] - self.data[key][3]*self.data['shalf'][2])/(self.data['shalf'][3] - self.data['shalf'][2]))
				self.data[key][-1] = ((self.data[key][-2] - self.data[key][-3])/(self.data['shalf'][-2] - self.data['shalf'][-3])*self.data['shalf'][-1]
									 + (self.data[key][-3]*self.data['shalf'][-2] - self.data[key][-2]*self.data['shalf'][-3])/(self.data['shalf'][-2] - self.data['shalf'][-3]))

		# profile variables s mesh with funky first point
		# use linear extrapolation
		LIST = ['chipf']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key][0]  = ((self.data[key][2] - self.data[key][1])/(self.data['s'][2] - self.data['s'][1])*self.data['s'][0] 
									 + (self.data[key][1]*self.data['s'][2] - self.data[key][2]*self.data['s'][1])/(self.data['s'][2] - self.data['s'][1]))

		# cut zeros off the end of aux profiles
		LIST = ['ac_aux', 'am_aux', 'ai_aux']
		for key in LIST:
			if self.data.has_key(key + '_s') & self.data.has_key(key + '_f'):
				idx = where(self.data[key + '_s'] != 0)[0]
				if len(idx) == 0: continue
				if not (idx[0] == 0): idx = append(0,idx)
				self.data[key + '_s'] = self.data[key + '_s'][idx]
				self.data[key + '_f'] = self.data[key + '_f'][idx]


	# --- data = getData() ---------------------------------------------------------------
	# returns a copy of self.data
	def getData(self):
		data = {}
		for key in self.data.keys():
			if isinstance(self.data[key], ndarray):
				data[key] = self.data[key].copy()
			else:
				data[key] = self.data[key]
		return data


	# --- data = getAll() ----------------------------------------------------------------
	# returns a copy of self.data and self.spl as one dict
	def getAll(self):
		data = self.getData()
		spl = self.spl.copy()
		for key in self.spl.keys():
			data['spl' + key] = spl[key]
		return data


	# --- data = getKeys() ---------------------------------------------------------------
	# returns all keys in self.data -> all variables in wout file
	def getKeys(self):
		return self.data.keys()


	# --- data = get(key) ----------------------------------------------------------------
	# returns variable key from wout file
	def get(self, key): 
		return self.data[key]


	# --- X = ev(keys, s, u, v, n0only, factor, ntor, derivs) ----------------------------
	# evaluates all spectral variables in the keys list
	# do NOT specify sin and cos keys separately, ONLY ONE OF THEM!!!
	# returns dictionary of variables
	# if keys has only one key, the output is the respective variable, not a dictionary!
	# s, u, v need to bee either arrays of same size or scalars
	# Input
	#	keys (string or list of strings):	variablenames in wout file; with or without 's' or 'c' at the end
	#	s (array or float):					flux coordiante
	#	u (array or float):					curv-linear poloidal angle
	#	v (array or float):					curv-linear toroidal angle
	#	n0only (bool):						full spectrum (default) or n = 0 only
	#	n_neq0_only (bool):					full spectrum (default) or n != 0 only
	#	factor(float):						amplification factor for n != 0 components, default = 1
	#	ntor (int):							evaluate only one n-component, or all (default)
	#	derivs (bool):						return angular derivatives as well, or not (default)
	#   VMECparity (bool):                  True: mu - nv;  False: mu + nv
	#   m0ony (bool):                       full spectrum (default) or m = 0 only
	# Output
	#	X (array or dictionary):			evaluated variable(s) (and optional: derivatives)
	#										keys are the same as input keys, but without 'mns' or 'mnc' at end
	#										derivs == True: for each key you have: key, dkeydu, dkeydv
	def ev(self, keys, s, u, v, n0only = False, factor = 1, ntor = None, derivs = False, 
			n_neq0_only = False, VMECparity = True, m0only = False):
		# make sure keys is a list even if it has only one element
		if not isinstance(keys, list):	
			keys = [keys]		
		
		# check sanity of keys and remove 'mns' or 'mnc' at the end
		for i, key in enumerate(keys):
			if not ('mn' in key):
				raise  RuntimeError('Wrong key to evaluate: ' + key)
			if (key[-1] == 's') | (key[-1] == 'c'):
				keys[i] = key[0:-1]		# remove 's' or 'c' at the end
		for i, key in enumerate(keys):
			if key[-1] == 'n':
				keys[i] = key[0:-2]		# remove 'mn' at the end
		
		
		# intialize
		if isinstance(s, ndarray): Ns = len(s)
		else: Ns = 1
		if isinstance(u, ndarray): Nu = len(u)
		else: Nu = 1
		if isinstance(v, ndarray): Nv = len(v)
		else: Nv = 1
		N = max([Ns, Nu, Nv])
		
		# set sign of nv
		if VMECparity: p = -1
		else: p = 1
		
		# loop through keys
		X = {}
		for key in keys:
			X[key] = zeros(N)
			if derivs:
				X['d' + key + 'du'] = zeros(N)
				X['d' + key + 'dv'] = zeros(N)
				
			# check for Nyquist spectrum
			if self.use_nyq & (not (key in ['r','z','l','drds','dzds'])):
				xn = 'xn_nyq'
				xm = 'xm_nyq'
				mnmax = 'mnmax_nyq'
			else:
				xn = 'xn'
				xm = 'xm'
				mnmax = 'mnmax'
			
			# which modes to sum
			idx0 = where(self.data[xn] == 0)[0]
			if n0only:
				idx = idx0
			elif n_neq0_only:
				idx = where(self.data[xn] != 0)[0]
			elif not (ntor is None):
				idx = where(self.data[xn] == ntor)[0]
			elif m0only:
				idx = where(self.data[xm] == 0)[0]
			else:
				idx = xrange(self.data[mnmax])
			
			# amplification of perturbation
			amp = factor * ones(self.data[xn].shape)
			amp[idx0] = 1.0
		
			# sum up modes
			for i in idx:
				m = self.data[xm][i]
				n = self.data[xn][i]
				cosuv = amp[i] * cos(m*u + p*n*v)
				sinuv = amp[i] * sin(m*u + p*n*v)
				if self.spl.has_key(key + 'mns'):
					spls = self.spl[key + 'mns'][i](s)
					if isinstance(s, ndarray):
						if(m > 0): spls[s == 0] = 0		# s = 0 -> no u dependence, so use only m=0 modes
					else:
						if((m > 0) & (s == 0)): spls = 0
					X[key] += spls * sinuv
					if derivs:
						X['d' + key + 'du'] += m * spls * cosuv
						X['d' + key + 'dv'] += p*n * spls * cosuv
				if self.spl.has_key(key + 'mnc'):
					splc = self.spl[key + 'mnc'][i](s)
					if isinstance(s, ndarray):
						if(m > 0): splc[s == 0] = 0		# s = 0 -> no u dependence, so use only m=0 modes
					else:
						if((m > 0) & (s == 0)): splc = 0
					X[key] += splc * cosuv
					if derivs:
						X['d' + key + 'du'] += -m * splc * sinuv
						X['d' + key + 'dv'] += -p*n * splc * sinuv
						
		if (not derivs) & (len(keys) == 1): 
			return X.values()[0]
		else:
			return X


	# --- make_RZderiv() -----------------------------------------------------------------
	# derivatives of R and Z with respect to s in spectral form
	# derivatives live on half grid then
	def make_RZderiv(self):
		self.data['drdsmnc'] = zeros(self.data['rmnc'].shape)
		self.data['dzdsmns'] = zeros(self.data['rmnc'].shape)
		
		if self.data['lasym']:
			self.data['drdsmns'] = zeros(self.data['rmnc'].shape)
			self.data['dzdsmnc'] = zeros(self.data['rmnc'].shape)
		
		# drds[i+1,:] = (r[i+1,:] - r[i,:]) / (s[i+1] - s[i])
		# drds is on the half grid and then i+1 is same as i+0.5 on full grid
		for i in xrange(self.data['ns']-1):
			self.data['drdsmnc'][i+1,:] = ((self.data['rmnc'][i+1,:] - self.data['rmnc'][i,:]) 
											/ (self.data['s'][i+1] - self.data['s'][i]))
			self.data['dzdsmns'][i+1,:] = ((self.data['zmns'][i+1,:] - self.data['zmns'][i,:]) 
											/ (self.data['s'][i+1] - self.data['s'][i]))
			if self.data['lasym']:
				self.data['drdsmns'][i+1,:] = ((self.data['rmns'][i+1,:] - self.data['rmns'][i,:]) 
												/ (self.data['s'][i+1] - self.data['s'][i]))
				self.data['dzdsmnc'][i+1,:] = ((self.data['zmnc'][i+1,:] - self.data['zmnc'][i,:]) 
												/ (self.data['s'][i+1] - self.data['s'][i]))
				
			
	# --- Raxis, Zaxis = get_axis(v, n0only) ---------------------------------------------
	# returns (R,Z) position of magnetic axis at toroidal angle v
	# Input
	#	v (float):			curv-linear toroidal angle
	#	 n0only (bool):		 full spectrum (default) or n = 0 only
	# Output
	#	 Raxis (float):		R of magnetic axis
	#	Zaxis (float):		Z of magnetic axis
	def get_axis(self, v, n0only = False, factor = 1):	
		if isinstance(self.data['raxis_cc'], ndarray):
			Raxis = self.data['raxis_cc'][0]
			if self.data['lasym']: 
				Zaxis = self.data['zaxis_cc'][0]
			else: 
				Zaxis = 0	# up-down symmetry
		else:
			Raxis = self.data['raxis_cc']
			if self.data['lasym']: 
				Zaxis = self.data['zaxis_cc']
			else: 
				Zaxis = 0	# up-down symmetry
			n0only = True
	
		if isinstance(v, ndarray):
			Raxis *= ones(len(v))
			Zaxis *= ones(len(v))
	
		if not n0only:
			for i in xrange(1, int(self.data['ntor']) + 1):
				n = self.data['xn'][i]
				cosnv = cos(n*v)
				sinnv = -sin(n*v)
				Raxis += factor * self.data['raxis_cc'][i] * cosnv
				Zaxis += factor * self.data['zaxis_cs'][i] * sinnv
				if self.data['lasym']:
					Raxis += factor * self.data['raxis_cs'][i] * sinnv
					Zaxis += factor * self.data['zaxis_cc'][i] * cosnv

		return Raxis, Zaxis

			
	# --- R, Z = get_surface(s, phi, ntheta, n0only, amp_fac, ntor) ----------------------
	# returns (R,Z) of single flux surface in 2-D poloidal cross-section
	# calls self.ev for given s, evenly spaced u with ntheta points and 
	# if machine == 'd3d': toroidal angle phi, given in DIII-D coordiantes -> left-handed and in degrees
	# else: toroidal angle phi is just v in degrees -> right-handed
	def get_surface(self, s = 1.0, phi = 0, ntheta = 200, n0only = False, amp_fac = 1.0, ntor = None):
		u = linspace(0, 2*pi, ntheta)
		if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
		else: v = (phi/180.0)%2 * pi
		
		dic = self.ev(['rmn','zmn'], s, u, v, n0only = n0only, factor = amp_fac, ntor = ntor)
		return dic['r'], dic['z']


	# --- dic = get_plasmaSurface(ntheta, nphi, n0only, amp_fac, ntor) -------------------
	# returns (x,y,z) cartesian coordinates of one or multiple plasma surface(s) in 3-D space
	# calls self.ev for each given s, evenly spaced u with nu points and 
	# evenly spaced v with nv points
	# Input
	#	s (double/array): normalized flux coordinate of surface(s), default s = 1
	#	nu (int):		number of poloidal points
	#	nv (int):		number of toroidal points
	#	 n0only (bool):	 full spectrum (default) or n = 0 only
	#	factor(float):	amplification factor for n != 0 components, default = 1
	#	ntor (int):		evaluate only one n-component, or all (default)
	# Output
	#	dic (dictionary):	'x','y','z','s'
	def get_plasmaSurface(self, s = 1.0, nu = 200, nv = 200, n0only = False, amp_fac = 1.0, ntor = None):		
		# s grid
		if (isinstance(s, int) | isinstance(s, float)): s = [s]
		
		# u, v grid
		u, v = meshgrid(linspace(0, 2*pi, nu), linspace(0, 2*pi, nv))
		
		x = zeros((nv, nu, len(s)))
		y = zeros((nv, nu, len(s)))
		z = zeros((nv, nu, len(s)))
		
		for i, s0 in enumerate(s):	
			# get surface
			dic = self.ev(['rmn','zmn'], s0, u.flatten(), v.flatten(), n0only = n0only, factor = amp_fac, ntor = ntor)

			R = dic['r'].reshape(u.shape)
			x[:,:,i] = R * cos(v)
			y[:,:,i] = R * sin(v)
			z[:,:,i] = dic['z'].reshape(u.shape)

		return {'x':x, 'y':y, 'z':z, 's':s}

		
	# --- u = midplane(s, v, n0only) ------------------------------------------------------
	# gets the vmec angle u for the geometric angle 0, so that Z(s,u,v) = Zaxis(v) 
	# Input
	#	s (double/array):   normalized flux coordinate of surface(s)
	#	v (float):			curv-linear toroidal angle
	#	n0only (bool):		full spectrum (default) or n = 0 only
	#   HFS (bool):         High Field Side or Low Field Side (default)
	# Output
	#	u (array or float):	curv-linear poloidal angle
	def midplane(self, s, v, n0only = False, HFS = False):
		eps = 1e-12
		_,Zaxis = self.get_axis(v, n0only = n0only)
		if HFS: u = pi*ones(s.shape)
		else: u = zeros(s.shape)
		du = ones(s.shape)
		idx = where(s == 0)
		if size(idx) > 0: fix_s_is_0 = True
		else: fix_s_is_0 = False
		n = 0
		while any(abs(du) > eps):
			Z = self.ev('zmn',s,u,v, n0only = n0only)
			dZ = (self.ev('zmn',s,u + 0.001,v, n0only = n0only) - Z)/(0.001)
			Z -= Zaxis
			du = Z/dZ
			if fix_s_is_0: du[idx] = 0
			u -= du
			n += 1
			if n > 30:
				print 'no convergence, max error: ', abs(du).max()
				break
		return u
	
		
	# --- dic = get_bcovar(s, u, v) ------------------------------------------------------
	# gets covariant magnetic field in curv-linear coordinates
	# calls self.ev for 'bsupumn' and 'bsupvmn' keys
	# !!! s only valid inside positive half-mesh: self.data['shalf'][1::] !!!
	# Output
	#	 dic (dictionary):	'bsupu','bsupv'
	#		bsupu (array):	B^u
	#		bsupv (array):	B^v
	def get_bcovar(self, s, u, v, n0only = False):
		return self.ev(['bsupumn','bsupvmn'], s, u, v, n0only = n0only)
	
	
	# --- dic = get_bcontra(s, u, v) -----------------------------------------------------
	# gets contra-variant magnetic field in curv-linear coordiantes
	# ??? DOES THIS REQUIRE A SIGN FLIP in m or n ??? - here: signs as usual for now !!!
	# calls self.ev for 'bsubsmn', 'bsubumn' and 'bsubvmn' keys
	# !!! s only valid inside positive half-mesh: self.data['shalf'][1::] !!!
	# Output
	#	 dic (dictionary):	'bsubu','bsubv','bsubs'
	#		bsubu (array):	B_u
	#		bsubv (array):	B_v
	#		bsubs (array):	B_s
	def get_bcontra(self, s, u, v, n0only = False):
		return self.ev(['bsubsmn','bsubumn','bsubvmn'], s, u, v, n0only = n0only)


	# --- g = get_jacobian(s, u, v) ----------------------------------------------------
	# gets the jacobian g; calls self.ev for 'gmn' key
	# !!! s only valid inside positive half-mesh: self.data['shalf'][1::] !!!
	# Output
	#	 g (array):	jacobian
	def get_jacobian(self, s, u, v, n0only = False):
		return self.ev('gmn', s, u, v, n0only = n0only)
		

	# --- dic = get_lamda(s, u, v) -------------------------------------------------------
	# gets lambda; calls self.ev for 'lmn' key; with derivatives
	# !!! s only valid inside positive half-mesh: self.data['shalf'][1::] !!!
	# Output
	#	 dic (dictionary):	'l','dldu','dldv'
	def get_lamda(self, s, u, v, n0only = False):
		return self.ev('lmn', s, u, v, derivs = True, n0only = n0only)
			
		
	# --- dic = get_B2D(phiwant, n0only, npts) -------------------------------------------
	# returns poloidal, toroidal and total magnetic fields from vmec 
	# Input
	#	 phi (double):		 toroidal angle in degrees, only for DIII-D: left-handed, default = 0
	#	 n0only (bool):		 full spectrum (default) or n = 0 only
	#	npts (int):			number of theta points, default = 128
	#	s (float/1D array/None), u (float/1D array/None), v (float/1D array/None):
	#		3-D coordinates; if either s, u or v is None (default), then polidal cross-section is used
	#		and s, u, v are ignored; else phi and npts are ignored and grid is given by s, u, v
	#		s only valid inside positive half-mesh: self.data['shalf'][1::]
	# Output
	#	 dict (dictionary):	'R','Z','Bpol2D','Btor2D','B2D','bees','s', 'u', 'v', 'theta'
	#		 R (array):		R coordinates in m
	#		Z (array):		Z coordinates in m
	#		phi (double/array):	toroidal angle coordinate in degrees
	#		Bpol2D (array):	poloidal magnetic field in T
	#		Btor2D (array):	toroidal magnetic field in T
	#		B2D (array):	total magnetic field in T
	#		BR2D (array):	Radial magnetic field in T
	#		BZ2D (array):	vertical magnetic field in T
	#		bees (dict):	output from self.get_bcovar
	#		s (array):		toroidal flux coordiante
	#		u (array):		curv-linear poloidal angle
	#		v (float/array):toroidal angle in rad and right-handed
	def get_B2D(self, phi = 0.0, n0only = False, npts = 128, s = None, u = None, v = None):	
		# set grid
		if (s is None) | (u is None) | (v is None):
			s = linspace(0, self.data['shalf'][-1], npts)
			u = linspace(0, 2*pi, npts)
			s, u = meshgrid(s, u)
			s = s.flatten()
			u = u.flatten()
			if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
			else: v = (phi/180.0)%2 * pi
		else:
			if self.data['machine'] == 'd3d': phi = 360.0 - v * 180.0/pi		# angle in degrees and left-handed (as in DIII-D)
			else: phi = v * 180.0/pi

		# get R,Z
		rz = self.ev(['rmn','zmn'], s, u, v, n0only = n0only, derivs = True)
		R = rz['r']
		Z = rz['z']
	
		# get B
		bees = self.get_bcovar(s, u, v, n0only)
		BR = rz['drdu']*bees['bsupu'] + rz['drdv']*bees['bsupv']
		BZ = rz['dzdu']*bees['bsupu'] + rz['dzdv']*bees['bsupv']
		Btor = R * bees['bsupv']
		Bpol = sqrt(BR**2 + BZ**2)
		Btot = sqrt(Bpol**2 + Btor**2)

		return {'R':R, 'Z':Z, 'phi':phi, 'Bpol2D':Bpol, 'Btor2D':Btor, 'B2D':Btot, 'BR2D':BR, 'BZ2D':BZ,
				'bees':bees, 's':s, 'u':u, 'v':v}
	
	
	# --- make_jcontra() ------------------------------------------------------------------
	# provides spectral representation of contra-variant local current density and adds it to self.data
	# also adds spline functions to self.spl
	# all jsub will be on half grid
	def make_jcontra(self):
		ns = self.data['ns']
		if self.use_nyq:
			mnmax = self.data['mnmax_nyq']
			xn = 'xn_nyq'
			xm = 'xm_nyq'
		else:
			mnmax = self.data['mnmax']
			xn = 'xn'
			xm = 'xm'
		shalf = self.data['shalf'][1::]		# ignore first value, since it is non-physical (negative)
		mu0 = 4*pi*1e-7
	
		self.data['jsubsmns'] = zeros((ns, mnmax))
		self.data['jsubumnc'] = zeros((ns, mnmax))
		self.data['jsubvmnc'] = zeros((ns, mnmax))
		
		if self.data['lasym']:
			self.data['jsubsmnc'] = zeros((ns, mnmax))
			self.data['jsubumns'] = zeros((ns, mnmax))
			self.data['jsubvmns'] = zeros((ns, mnmax))

		# mu0 * j = curl B = nabla x B;  nabla = (d/ds, d/du, d/dv);  B = (Bsubs, Bsubu, Bsubv)
		# jsubs = dBsubv/du - dBsubu/dv;  jsubu = dBsubs/dv - dBsubv/ds;  jsubv = dBsubu/ds - dBsubs/du
		# Bsubu & Bsubv are on half-grid;  Bsubs is on full grid;  Bx = sum(bxmns * sin(mu - nv) + bxmnc * cos(mu - nv))
		for i in xrange(mnmax):
			n = self.data[xn][i]
			m = self.data[xm][i]
			bsubsmns = append(0, self.spl['bsubsmns'][i](shalf))		# this one is on full grid -> move to half-grid
			self.data['jsubsmns'][:,i] = -m * self.data['bsubvmnc'][:,i] - n * self.data['bsubumnc'][:,i]
			self.data['jsubumnc'][:,i] = append(0, -deriv(self.data['bsubvmnc'][1::,i], shalf)) - n * bsubsmns
			self.data['jsubvmnc'][:,i] = append(0, deriv(self.data['bsubumnc'][1::,i], shalf)) - m * bsubsmns
		
			if self.data['lasym']:
				bsubsmnc = append(0, self.spl['bsubsmnc'][i](shalf))	# this one is on full grid -> move to half-grid
				self.data['jsubsmnc'][:,i] = m * self.data['bsubvmns'][:,i] + n * self.data['bsubumns'][:,i]
				self.data['jsubumns'][:,i] = append(0, -deriv(self.data['bsubvmns'][1::,i], shalf)) + n * bsubsmnc
				self.data['jsubvmns'][:,i] = append(0, deriv(self.data['bsubumns'][1::,i], shalf)) + m * bsubsmnc
				
		self.data['jsubsmns'] /= mu0
		self.data['jsubumnc'] /= mu0
		self.data['jsubvmnc'] /= mu0
		
		if self.data['lasym']:
			self.data['jsubsmnc'] /= mu0
			self.data['jsubumns'] /= mu0
			self.data['jsubvmns'] /= mu0
				
		# get spline functions; they all live on the half-grid
		if self.data['lasym']:	# asymmetric (not stellarator symmetric)
			LIST = ['jsubsmns', 'jsubsmnc', 'jsubumns', 'jsubumnc', 'jsubvmns', 'jsubvmnc']
		else:					# up-down symmetric (stellarator symmetric)
			LIST = ['jsubsmns', 'jsubumnc', 'jsubvmnc']

		for key in LIST:
			self.spl[key] = [0] * mnmax
			self.data[key][0,:] = self.data[key][1,:] * (-1)**self.data[xm]
			for i in xrange(mnmax):
				self.spl[key][i] = scinter.UnivariateSpline(self.data['shalf'], self.data[key][:,i], s = 0)


	# --- dic = get_jcontra(s, u, v) ------------------------------------------------------
	# gets contra-variant local current density; calls self.ev for 'jsubsmn','jsubumn' and 'jsubvmn' keys
	# Output
	#	 dic (dictionary):	'jsubs','jsubu','jsubv'
	#		jsubs (array):	j_s
	#		jsubu (array):	j_u
	#		jsubv (array):	j_v
	def get_jcontra(self, s, u, v, n0only = False, n_neq0_only = False):
		if not self.data.has_key('jsubsmns'):
			self.make_jcontra()
		return self.ev(['jsubsmn','jsubumn','jsubvmn'], s, u, v, n0only = n0only, n_neq0_only = n_neq0_only)

	
	# --- dic = get_j2D(phiwant, n0only, npts) -------------------------------------------
	# returns poloidal, toroidal, parallel and total current density from vmec 
	# also returns poloidal, toroidal and total magnetic fields from vmec, same as self.get_B2D()
	# Input
	#	 phi (double):		 toroidal angle in degrees, only for DIII-D: left-handed, default = 0
	#	 n0only (bool):		 full spectrum (default) or n = 0 only
	#	npts (int):			number of theta points, default = 128
	#	s (float/1D array/None), u (float/1D array/None), v (float/1D array/None):
	#		3-D coordinates; if either s, u or v is None (default), then polidal cross-section is used
	#		and s, u, v are ignored; else phi and npts are ignored and grid is given by s, u, v
	#		s only valid inside positive half-mesh: self.data['shalf'][1::]
	# Output
	#	 dict (dictionary):	'R','Z','jpol2D','jtor2D','j2D','jpar2D','jees','s', 'u', 'v', 'theta'
	#		 R (array):		R coordinates in m
	#		Z (array):		Z coordinates in m
	#		phi (double/array):	toroidal angle coordinate in degrees
	#		jpol2D (array):	poloidal current density in A/m^2
	#		jtor2D (array):	toroidal current density in A/m^2
	#		j2D (array):	total current density in A/m^2
	#		jpar2D (array):	parallel current density in A/m^2
	#		jR2D (array):	Radial current density in A/m^2
	#		jZ2D (array):	vertical current density in A/m^2
	#		Bpol2D (array):	poloidal magnetic field in T
	#		Btor2D (array):	toroidal magnetic field in T
	#		B2D (array):	total magnetic field in T
	#		BR2D (array):	Radial magnetic field in T
	#		BZ2D (array):	vertical magnetic field in T
	#		bees (dict):	output from self.get_bcovar
	#		jees (dict):	output from self.get_jcontra
	#		 g (array):		jacobian
	#		s (array):		toroidal flux coordiante
	#		u (array):		curv-linear poloidal angle
	#		v (double):		toroidal angle in rad and right-handed
	#		theta (array):	geometric poloidal angle
	def get_j2D(self, phi = 0.0, n0only = False, npts = 128, s = None, u = None, v = None):	
		# set grid
		if (s is None) | (u is None) | (v is None):
			s0 = linspace(0, self.data['shalf'][-1], npts)
			u0 = linspace(0, 2*pi, npts)
			s, u = meshgrid(s0, u0)
			s = s.flatten()
			u = u.flatten()
			if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
			else: v = (phi/180.0)%2 * pi
		else:
			if self.data['machine'] == 'd3d': phi = 360.0 - v * 180.0/pi		# angle in degrees and left-handed (as in DIII-D)
			else: phi = v * 180.0/pi

		# get R,Z
		rz = self.ev(['rmn','zmn'], s, u, v, n0only = n0only, derivs = True)
		R = rz['r']
		Z = rz['z']
		
		# get jacobian
		g = self.get_jacobian(s, u, v, n0only)
	
		# get jcontra and turn it into covar
		jees = self.get_jcontra(s, u, v, n0only)
		#jsups = jees['jsubs'] / g
		jsupu = jees['jsubu'] / g
		jsupv = jees['jsubv'] / g
		
		# get j
		jR = rz['drdu']*jsupu + rz['drdv']*jsupv
		jZ = rz['dzdu']*jsupu + rz['dzdv']*jsupv
		jtor = R * jsupv
		jpol = sqrt(jR**2 + jZ**2)
		jtot = sqrt(jR**2 + jZ**2 + jtor**2)
	
		# get B
		bees = self.get_bcovar(s, u, v, n0only)
		BR = rz['drdu']*bees['bsupu'] + rz['drdv']*bees['bsupv']
		BZ = rz['dzdu']*bees['bsupu'] + rz['dzdv']*bees['bsupv']
		Btor = R * bees['bsupv']
		Bpol = sqrt(BR**2 + BZ**2)
		Btot = sqrt(Bpol**2 + Btor**2)

		# jpar = (j dot B)/B
		jpar = (jR*BR + jZ*BZ + jtor*Btor)/Btot

		return {'R':R, 'Z':Z, 'phi':phi, 'jpol2D':jpol, 'jtor2D':jtor, 'j2D':jtot, 'jpar2D':jpar, 'jR2D':jR, 'jZ2D':jZ,
				'Bpol2D':Bpol, 'Btor2D':Btor, 'B2D':Btot, 'BR2D':BR, 'BZ2D':BZ, 'bees':bees,
				'jees':jees, 'g':g, 's':s, 'u':u, 'v':v}


	# --- dic = get_j2D_parallel(phiwant, n0only, npts) ----------------------------------
	# same as self.get_j2D, but computes in parallel. 
	# Flat s, u, v arrays are divided into chunks of 100 elements per process
	def get_j2D_parallel(self, phi = 0.0, n0only = False, npts = 128, s = None, u = None, v = None):
		import multiprocessing as mp
		import os

		# set grid
		if (s is None) | (u is None) | (v is None):
			s0 = linspace(0, self.data['shalf'][-1], npts)
			u0 = linspace(0, 2*pi, npts)
			s, u = meshgrid(s0, u0)
			s = s.flatten()
			u = u.flatten()
			if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
			else: v = (phi/180.0)%2 * pi
		else:
			if self.data['machine'] == 'd3d': phi = 360.0 - v * 180.0/pi		# angle in degrees and left-handed (as in DIII-D)
			else: phi = v * 180.0/pi
		
		is_s_array = isinstance(s, ndarray)
		if is_s_array: Ns = len(s)
		else: Ns = 1
		
		is_u_array = isinstance(u, ndarray)
		if is_u_array: Nu = len(u)
		else: Nu = 1
		
		is_v_array = isinstance(v, ndarray)
		if is_v_array: Nv = len(v)
		else: Nv = 1
		
		M = max([Ns, Nu, Nv])		# total number of elements
		Nch = 5000					# number of elements per chunk
		N = int(ceil(M/float(Nch)))	# number of chunks
		
		# Prepare Processes
		if N < 10: numProcs = N
		else: numProcs = 10
		ProcID = range(numProcs)	# start numProcs processes
		Lock = mp.Lock()
		Procs = []
	
		# these are memory shared objects, all procs write into them, they have to be 1D arrays, 'd' for double
		out =  {'R':mp.Array('d', M, lock = False), 'Z':mp.Array('d', M, lock = False), 
				'jR':mp.Array('d', M, lock = False), 'jZ':mp.Array('d', M, lock = False), 
				'jtor':mp.Array('d', M, lock = False), 
				'BR':mp.Array('d', M, lock = False), 'BZ':mp.Array('d', M, lock = False), 
				'Btor':mp.Array('d', M, lock = False)}
				
		# Set the Queue Size: Queue has built-in maximum length of 32767. If Qsize > 32767 -> Error occurs!
		Qsize_limit = 1000			# don't has to be the maximum length
		if (N > Qsize_limit):
			Qsize = Qsize_limit
		else:
			Qsize = N
		
		# Total number of elements to go through the queue -> N
		Qlength = mp.Value('i', N)	# memory shared variable, set to N, 'i' for int

		# Set the Queue and fill it
		Q = mp.Queue(Qsize)				# If Qsize > 32767 -> Error occurs!
		for i in range(Qsize):
			Q.put(i)
	
		# Create all processes and run them
		for ID in ProcID:
			proc = mp.Process(target = self.__get_j2D_loop__, args = (Lock, Q, Qlength, out, N, Nch, 
												s, u, v, is_s_array, is_u_array, is_v_array, n0only))
			Procs.append(proc)
			proc.start()			# runs process
	
		# Continue filling the Queue, put waits until slot opens up in Queue
		for i in range(Qsize, N):	# if Qsize == N, nothing happens here
			Q.put(i)
	
		# Wait for all processes to complete
		for t in Procs:
			t.join()

		# Done multiprocessing -> collect data
		R = frombuffer(out['R'])
		Z = frombuffer(out['Z'])
		
		jR = frombuffer(out['jR'])
		jZ = frombuffer(out['jZ'])
		jtor = frombuffer(out['jtor'])
		jpol = sqrt(jR**2 + jZ**2)
		jtot = sqrt(jR**2 + jZ**2 + jtor**2)
		
		BR = frombuffer(out['BR'])
		BZ = frombuffer(out['BZ'])
		Btor = frombuffer(out['Btor'])
		Bpol = sqrt(BR**2 + BZ**2)
		Btot = sqrt(Bpol**2 + Btor**2)

		# jpar = (j dot B)/B
		jpar = (jR*BR + jZ*BZ + jtor*Btor)/Btot
		
		return {'R':R, 'Z':Z, 'phi':phi, 'jpol2D':jpol, 'jtor2D':jtor, 'j2D':jtot, 'jpar2D':jpar, 'jR2D':jR, 'jZ2D':jZ,
				'Bpol2D':Bpol, 'Btor2D':Btor, 'B2D':Btot, 'BR2D':BR, 'BZ2D':BZ, 'bees':None,
				'jees':None, 'g':None, 's':s, 'u':u, 'v':v}
				
				
	# --- __get_j2D_loop__() -------------------------------------------------------------
	# private function!  loop body, called by self.get_j2D_parallel
	def __get_j2D_loop__(self, Lock, Q, Qlength, out, N, Nch, s0, u0, v0, is_s_array, is_u_array, is_v_array, n0only):
		while not Q.empty():
			# get one element out of the queue
			Lock.acquire()
			# make sure, Queue is really not empty, now within a Lock with Q.get()
			if Q.empty():		 # if it is, exit here
				Lock.release()	# don't forget to unlock!
				break;
			i = Q.get()
			Qlength.value -= 1
			Lock.release()
		
			# Do the work
			iStart = Nch * i
			iEnd = Nch * (i + 1)
			if (i == N - 1): iEnd = None

			if is_s_array: s = s0[iStart:iEnd]
			else: s = s0
			if is_u_array: u = u0[iStart:iEnd]
			else: u = u0
			if is_v_array: v = v0[iStart:iEnd]
			else: v = v0

			# get R,Z
			rz = self.ev(['rmn','zmn'], s, u, v, n0only = n0only, derivs = True)
			R = rz['r']; out['R'][iStart:iEnd] = R
			out['Z'][iStart:iEnd] = rz['z']
			
			# get dR/ds & dZ/ds
			 #ds = min([0.1*(1 - self.data['shalf'][-1]), 1e-5])	# at least 1e-5, at most a tenth of the shalf step to s = 1
			 #drz = self.ev(['rmn','zmn'], s + ds, u, v, n0only = n0only)
			 #drds = (drz['r'] - rz['r'])/ds
			 #dzds = (drz['z'] - rz['z'])/ds
		
			# get jacobian
			g = self.get_jacobian(s, u, v, n0only)
	
			# get jcontra and turn it into covar
			jees = self.get_jcontra(s, u, v, n0only)
			#jsups = jees['jsubs'] / g
			jsupu = jees['jsubu'] / g
			jsupv = jees['jsubv'] / g
		
			# get j
			out['jR'][iStart:iEnd] = rz['drdu']*jsupu + rz['drdv']*jsupv #+ drds*jsups
			out['jZ'][iStart:iEnd] = rz['dzdu']*jsupu + rz['dzdv']*jsupv #+ dzds*jsups
			out['jtor'][iStart:iEnd] = R * jsupv
	
			# get B
			bees = self.get_bcovar(s, u, v, n0only)
			out['BR'][iStart:iEnd] = rz['drdu']*bees['bsupu'] + rz['drdv']*bees['bsupv']
			out['BZ'][iStart:iEnd] = rz['dzdu']*bees['bsupu'] + rz['dzdv']*bees['bsupv']
			out['Btor'][iStart:iEnd] = R * bees['bsupv']


	# --- delta_r = getSurfDisp(s, phi) --------------------------------------------------
	# returns displacement dr of surface s compared to axisymmetric (n = 0) case
	# theta is geometric angle, u is VMEC angle and usfl is a straight field line angle
	def getSurfDisp(self, s, phi = 0, ntheta = 400):
		# get axis for n0only
		Raxis, Zaxis  = self.get_axis(0, n0only = True)

		# get axisymmetric part only
		rsurf, zsurf = self.get_surface(s, phi = phi, n0only = True, ntheta = ntheta)
	
		theta = self.__get_theta__(rsurf, zsurf, None, Raxis, Zaxis)	# poloidal angle theta
		r = self.__get_r__(rsurf, zsurf, None, Raxis, Zaxis)			# minor radius
		
		# first value is often the same as the last one, so remove
		theta = theta[1::]					
		r = r[1::]
	
		# sometimes first values are about 2pi, so sort along theta
		idx = theta.argsort()
		theta = theta[idx]
		r = r[idx]
	
		# make periodic
		theta = append(theta - 2*pi, append(theta, theta + 2*pi))
		r = append(r, append(r, r))
	
		# interpolate r
		f_r0 = scinter.UnivariateSpline(theta, r, s = 0)
	
		# get full surface
		rsurf, zsurf = self.get_surface(s, phi = phi, n0only = False, ntheta = ntheta)
	
		theta = self.__get_theta__(rsurf, zsurf, None, Raxis, Zaxis)	# poloidal angle theta
		r = self.__get_r__(rsurf, zsurf, None, Raxis, Zaxis)			# minor radius
		idx = theta.argsort()							# sort along theta
		theta = theta[idx]
		r = r[idx]

		# subtract axisymmetric part; convert to cm
		dr = (r - f_r0(theta)) * 100
		
		# get straight field line coordinate
		if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
		else: v = (phi/180.0)%2 * pi
		u = linspace(0,2*pi,ntheta)	# default from self.get_surface
		usfl = u + self.ev('lmn', s, u, v)
		
		return {'theta':theta, 'dr':dr, 'u':u, 'usfl':usfl}


	# --- delta_r = getSurfDisp(s, phi) --------------------------------------------------
	# returns displacement dr of surface s compared to axisymmetric (n = 0) case
	# theta is geometric angle, u is VMEC angle and usfl is a straight field line angle
	def getSurfDispTorroidal(self, s, u = 0, nv = 400):
		# get axis for n0only
		Raxis, Zaxis  = self.get_axis(0, n0only = True)
		v = linspace(0,2*pi,nv)

		# get axisymmetric part only
		dic = self.ev(['rmn','zmn'], s, u, v, n0only = True)
		rsurf, zsurf = dic['r'], dic['z']
	
		theta0 = self.__get_theta__(rsurf, zsurf, None, Raxis, Zaxis)	# poloidal angle theta
		r0 = self.__get_r__(rsurf, zsurf, None, Raxis, Zaxis)			# minor radius
		theta = theta0.mean()	# geometric angle that corresponds to u (in axisymmetric n = 0 limit)
	
		# get full surface
		dic = self.ev(['rmn','zmn'], s, u, v, n0only = False)
		rsurf, zsurf = dic['r'], dic['z']	
		r = self.__get_r__(rsurf, zsurf, None, Raxis, Zaxis)			# minor radius

		# subtract axisymmetric part; convert to cm
		dr = (r - r0) * 100
				
		return {'theta':theta, 'dr':dr, 'u':u, 'v':v}

		
	# --------------------------------------------------------------------------------
	# Flux surface enclosed 2D-Volume (Area inside surface) and derivative (s)
	def volume2D(self, phi = 0, n0only = False, ntheta = 200):
		V = zeros(self.data['ns'])
		s = self.data['s']
		
		# get magnetic axis
		if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
		else: v = (phi/180.0)%2 * pi
		Raxis, Zaxis  = self.get_axis(v)

		for i in xrange(1, self.data['ns']):
			Rs, Zs = self.get_surface(s[i], phi = phi, n0only = n0only, ntheta = ntheta)
			rsq = (Rs - Raxis)**2 + (Zs - Zaxis)**2
			theta = self.__get_theta__(Rs, Zs, phi, Raxis, Zaxis)
			idx = where(abs(diff(theta)) > 5)[0]			# locate 2pi jump in theta, if any
			if(len(idx) > 0): theta[idx[0] + 1 ::] += 2*pi		# make theta increase monotonically
			V[i] = 0.5 * integ.simps(rsq, theta)
		
		# dV/ds
		dV = deriv(V, s)

		return {'V':V, 'Vprime':dV}

		
	# --------------------------------------------------------------------------------
	# get Pprime with respect to toroidal flux (phi) and poloidal flux (chi)
	def pprime(self):
		dPdphi = deriv(self.data['presf'],abs(self.data['phi']))
		dPdchi = deriv(self.data['presf'],abs(self.data['chi'])) * 2*pi
		return {'dPdphi':dPdphi, 'dPdchi':dPdchi}
		
		
	# --------------------------------------------------------------------------------
	# get flux surface averaged Fpol, Fprime and FFprime with respect to poloidal flux (chi)
	# axisymmetric only, to compare with Grad-Shafranov solutions
	def Grad_Shafranov_Fpol(self):
		f = scinter.UnivariateSpline(self.data['shalf'][1::], self.data['bvco'][1::],s = 0)
		F = f(self.data['s'][1:-1])
		F = append(append(self.data['rbtor0'], F), self.data['rbtor'])
		Fprime = deriv(F, abs(self.data['chi'])) * 2*pi;
		return {'Fpol':F, 'Fprime':Fprime, 'FFprime':F*Fprime}
		
		
	# --- R, Z = getBoundary(phiwant, ntheta) --------------------------------------------
	# same as self.get_surface(...) with s = 1; kept for back-compatibility
	def getBoundary(self, phiwant = 0, ntheta = 200):
		return self.get_surface(phi = phiwant, ntheta = ntheta)
		
		
	# --- R, Z = getBoundPerturb(phiwant, ntheta, amp_fac) -------------------------------
	# same as self.get_surface(...) with s = 1 and factor = amp_fac; kept for back-compatibility
	def getBoundPerturb(self, phiwant = 0, ntheta = 200, amp_fac = 50):
		return self.get_surface(phi = phiwant, ntheta = ntheta, amp_fac = amp_fac)


	# --- delta_r = getBoundPerturbAmp(phiwant, ntheta, amp_fac) -------------------------
	# returns amplitude of boundary surface deviation
	# calls self.get_surface(...) with s = 1, factor = amp_fac and ntor; kept for back-compatibility		
	def getBoundPerturbAmp(self, phiwant = 0, ntheta = 200, amp_fac = 1.0, ntor = None):
		R, Z = self.get_surface(phi = phiwant, ntheta = ntheta, amp_fac = amp_fac, ntor = ntor)
		if (ntor is None):
			R0, Z0 = self.get_surface(phi = phiwant, ntheta = ntheta, n0only = True)
		else:
			R0, Z0 = 0, 0
		return sqrt((R - R0)**2 + (Z - Z0)**2)
		
		
	# --- dic = bndy_shaping() ---
	# returns all n = 0 plasma shaping parameters, according to paper by T. Luce
	# Plasma Phys. Control. Fusion 57 (2015) 049501
	def bndy_shaping(self):
		from numpy.linalg.linalg import norm
		R, Z = self.get_surface(s = 1.0, n0only = True, ntheta = 2000)
		idx = R.argmax()
		Rmax = R[idx]
		Z_Rmax = Z[idx]
		idx = R.argmin()
		Rmin = R[idx]
		Z_Rmin = Z[idx]
		idx = Z.argmax()
		R_Zmax = R[idx]
		Zmax = Z[idx]
		idx = Z.argmin()
		R_Zmin = R[idx]
		Zmin = Z[idx]
		
		Rgeo = 0.5*(Rmax + Rmin)	# major radius
		a = 0.5*(Rmax - Rmin) 		# minor radius
		eps = a/Rgeo				# inverse aspect ratio
		
		# Elongation
		kappa = 0.5*(Zmax - Zmin)/a	# elongation
		kappa_u = (Zmax - Z_Rmax)/a	# upper elongation
		kappa_l = (Z_Rmax - Zmin)/a	# lower elongation
		
		# Triangularity
		delta_u = (Rgeo - R_Zmax)/a	# upper triangularity
		delta_l = (Rgeo - R_Zmin)/a	# lower triangularity
		
		# squareness
		# 0 = lower outer, 1 = lower inner, 2 = upper outer, 3 = upper inner
		zeta = zeros(4)
		Es = [[Rmax,Zmin], [Rmin,Zmin], [Rmax,Zmax], [Rmin,Zmax]]
		Os = [[R_Zmin,Z_Rmax], [R_Zmin,Z_Rmin], [R_Zmax,Z_Rmax], [R_Zmax,Z_Rmin]]
		
		for i,E,O in zip(arange(4),Es,Os):
			E = array(E); O = array(O)
			OE = E - O	# O -> E
			LOE = norm(OE)
			LOC = LOE/sqrt(2)
			LCE = LOE - LOC
		
			t1 = (R - O[0])/OE[0]
			t2 = (Z - O[1])/OE[1]
			t1[(R > max(E[0],O[0])) | (R < min(E[0],O[0]))] = 1e+10
			idx = abs(t2-t1).argmin()
			D = array([R[idx], Z[idx]])
			OD = D - O
			LOD = norm(OD)		
			zeta[i] = (LOD - LOC)/LCE
				
		return {'Rgeo':Rgeo, 'a':a, 'eps':eps, 'Zoff':Z_Rmax, 'asp_ratio':1.0/eps,
				'kappa':kappa, 'kappa_u':kappa_u, 'kappa_l':kappa_l, 
				'delta_u':delta_u, 'delta_l':delta_l, 
				'zeta_lo':zeta[0], 'zeta_li':zeta[1], 'zeta_uo':zeta[2], 'zeta_ui':zeta[3]}
	
	
	# --- R,Z = ev_bndy_shape(theta, dic) ---
	# calculates the n = 0 plasma boundary shape based on shaping parameters in dic
	# theta is a poloidal angle like parameter, NOT the geometric poloidal angle or u
	# according to paper by T. Luce, Plasma Phys. Control. Fusion 57 (2015) 049501
	def ev_bndy_shape(self, theta, dic = None):
		if dic is None: dic = self.bndy_shaping()
		Roff_u = dic['a']*(dic['asp_ratio'] - dic['delta_u'])
		Roff_l = dic['a']*(dic['asp_ratio'] - dic['delta_l'])
		n1 = log(2)
		n2 = 1.0/sqrt(2)
		n3 = 1 - n2
		
		if isinstance(theta, ndarray):
			a = zeros(theta.shape)
			b = 10*ones(theta.shape)
		else:
			a = 0
			b = 10
		
		N = 40	# (b-a) < 1e-12
		
		for i in xrange(N):
			r = 0.5*(a + b)
			R = r*cos(theta) + dic['Rgeo']
			Z = r*sin(theta) + dic['Zoff']
		
			if isinstance(theta, ndarray):
				Q = ones(theta.shape)
				x = zeros(theta.shape)
				y = zeros(theta.shape)
				A = zeros(theta.shape)
				B = zeros(theta.shape)
				n = zeros(theta.shape)
				Q[(Z > dic['Zoff']) & (R <= Roff_u)] = 2
				Q[(Z <= dic['Zoff']) & (R <= Roff_l)] = 3
				Q[(Z <= dic['Zoff']) & (R > Roff_l)] = 4
				# Quadrant I
				x[Q==1] = R[Q==1] - Roff_u
				y[Q==1] = Z[Q==1] - dic['Zoff']
				A[Q==1] = dic['a']*(1 + dic['delta_u'])
				B[Q==1] = dic['a']*dic['kappa_u']
				n[Q==1] = -n1/log(n2 + dic['zeta_uo']*n3)
				# Quadrant II
				x[Q==2] = Roff_u - R[Q==2]
				y[Q==2] = Z[Q==2] - dic['Zoff']
				A[Q==2] = dic['a']*(1 - dic['delta_u'])
				B[Q==2] = dic['a']*dic['kappa_u']
				n[Q==2] = -n1/log(n2 + dic['zeta_ui']*n3)		
				# Quadrant IV
				x[Q==4] = R[Q==4] - Roff_l
				y[Q==4] = dic['Zoff'] - Z[Q==4]
				A[Q==4] = dic['a']*(1 + dic['delta_l'])
				B[Q==4] = dic['a']*dic['kappa_l']
				n[Q==4] = -n1/log(n2 + dic['zeta_lo']*n3)
				# Quadrant III
				x[Q==3] = Roff_l - R[Q==3]
				y[Q==3] = dic['Zoff'] - Z[Q==3]
				A[Q==3] = dic['a']*(1 - dic['delta_l'])
				B[Q==3] = dic['a']*dic['kappa_l']
				n[Q==3] = -n1/log(n2 + dic['zeta_li']*n3)
				
			else:
				if Z > dic['Zoff']:
					y = Z - dic['Zoff']
					if R > Roff_u:	# Quadrant I
						x = R - Roff_u
						A = dic['a']*(1 + dic['delta_u'])
						B = dic['a']*dic['kappa_u']
						n = -n1/log(n2 + dic['zeta_uo']*n3)
					else:			# Quadrant II
						x = Roff_u - R
						A = dic['a']*(1 - dic['delta_u'])
						B = dic['a']*dic['kappa_u']
						n = -n1/log(n2 + dic['zeta_ui']*n3)		
				else:
					y = dic['Zoff'] - Z
					if R > Roff_l:	# Quadrant IV
						x = R - Roff_l
						A = dic['a']*(1 + dic['delta_l'])
						B = dic['a']*dic['kappa_l']
						n = -n1/log(n2 + dic['zeta_lo']*n3)
					else:			# Quadrant III
						x = Roff_l - R
						A = dic['a']*(1 - dic['delta_l'])
						B = dic['a']*dic['kappa_l']
						n = -n1/log(n2 + dic['zeta_li']*n3)		
			
			f = (x/A)**n + (y/B)**n - 1

			if isinstance(theta, ndarray):
				b[f > 0] = r[f > 0]
				a[f < 0] = r[f < 0]
			else:
				if f > 0: b = r
				else: a = r
		
		return R,Z


	# --- dR = helical_core_inflection_point() ---
	# calculates the radial point where the flux compression is max
	def helical_core_inflection_point(self, full = False,vertical = False):
		from scipy.optimize import bisect
		if self.data['ntor'] == 0: return 0		# no n = 1 mode available
		
		# Move helical core to outer midplane
		ans = self.data['raxis_cs'][1]    # n = 1 component of the raxis_cs variable
		anc = self.data['raxis_cc'][1]
		v = arctan(-ans/anc)
		if anc < 0: v += pi
		u = 0
		if vertical: 
			v += pi/2
			u = pi/2

		rdic = self.ev(['rmn','zmn'],self.data['s'],u,v)
		R,Z = rdic['r'],rdic['z']
		Raxis,Zaxis = self.get_axis(0, n0only = True)
		r = sqrt((R-Raxis)**2 + (Z-Zaxis)**2)

		# locate point where curvature of r(s) changes from compression to expansion
		dr = deriv(r,self.data['s'])
		d2r = deriv(dr,self.data['s'])
		f = scinter.UnivariateSpline(self.data['s'],d2r,s = 0)
		
		# find d^2r/ds^2 = 0
		s0 = 0.005
		a,b = None,None
		while (a is None) | (b is None):
			if s0 > 1: 
				print 'Fail to locate root'
				break
			y0 = f(s0)
			if y0 < 0: a = s0
			else: b = s0
			s0 += 0.015

		s0 = bisect(f,a,b)
		
		R = self.ev('rmn',s0,linspace(0,2*pi,300),v)
		dR = 0.5*(R.max() - R.min())
		Z = self.ev('zmn',s0,linspace(0,2*pi,300),v)
		dZ = 0.5*(Z.max() - Z.min())
		dr = sqrt(dR**2 + dZ**2)/sqrt(2.0)
		da = dr/self.data['aminor']
		
		if full: return {'s0':s0, 'dR':dR, 'dZ':dZ, 'dr':dr, 'da':da, 'r':r, 'd2r':d2r, 'v':v}
		else: return s0


	# --- rho_fourier_coeff() ---
	# calculates spectral compression coefficients for rho (= minor r) conversion
	# mexp = 0 for polar,  = 1 for equal arclength, > 1 for smaller spectral width
	# mexp = 4 corresponds to STELLOPT & VMEC choice
	# WARNING: These coefficients only work for fixed boundary created by DESCUR
	def rho_fourier_coeff(self, mexp = 4):
		mpol = self.data['mpol']+2
		t1m = zeros(mpol)
		t2m = zeros(mpol)
		t2m[0] = 0.5**mexp
		t1m[1] = 1
		t1m[2] = t2m[0]	# BETTER FIT IF AXIS IS FIXED (mimics polar m=1 coupling)
		for m in xrange(1,mpol-1):
			t2m[m] = (float64(m+1)/m)**mexp
		for m in xrange(3,mpol):
			t1m[m] = (float64(m-1)/m)**mexp	
		return t1m,t2m
		
	
	# --- rho_fourier_center() ---
	# calculates center point for rho (= minor r) conversion
	# so that Rrev00c == R00c and Zrev00c == Z00c
	def rho_fourier_center(self):
		off = 2*self.data['ntor'] + 1
		R0 = self.data['rmnc'][:,0] - 0.5*(self.data['rmnc'][:,2*off] + self.data['zmns'][:,2*off]) #*t2m[0]/t1m[2]
		if self.data['lasym']:
			Z0 = self.data['zmnc'][:,0] - 0.5*(self.data['rmns'][:,2*off] - self.data['zmnc'][:,2*off]) #*t2m[0]/t1m[2]
		else:
			Z0 = 0
		return R0,Z0

			
	# --- get_rho_fourier() ---
	# calculates Fourier modes for rho (= minor r) from R,Z modes
	def get_rho_fourier(self):
		if self.data['lfreeb']:
			print 'Warning: In free boundary mode this can only give an approximate fit. Exact fits only in fixed boundary mode.'
		ntor = self.data['ntor']
		off = 2*ntor + 1
		self.data['rhomnc'] = zeros((self.data['ns'],self.data['mnmax']+off))
		t1m,t2m = self.rho_fourier_coeff()
		R0,Z0 = self.rho_fourier_center()
		LIST = ['rhomnc']			
			
		# all n, m >= 2 and n > 0, m = 1
		for i in xrange(off+1, self.data['mnmax']-off):
			m = self.data['xm'][i]
			tnorm = 0.5/(t1m[m+1]**2 + t2m[m-1]**2)
			self.data['rhomnc'][:,i] = tnorm*((self.data['rmnc'][:,i+off] + self.data['zmns'][:,i+off])*t1m[m+1]
												+ (self.data['rmnc'][:,i-off] - self.data['zmns'][:,i-off])*t2m[m-1])
		# all n, m = mpol-1, m = mpol
		for i in xrange(self.data['mnmax']-off, self.data['mnmax']+off):
			if i < self.data['mnmax']: m = self.data['mpol'] - 1
			else: m = self.data['mpol']
			tnorm = 0.5/(t1m[m+1]**2 + t2m[m-1]**2)
			self.data['rhomnc'][:,i] = tnorm*(self.data['rmnc'][:,i-off] - self.data['zmns'][:,i-off])*t2m[m-1]
			
		# m = 0 and m = 1, n <= 0
		for i in xrange(0, off+1):
			m = self.data['xm'][i]
			self.data['rhomnc'][:,i] = 0.5*(self.data['rmnc'][:,i+off] + self.data['zmns'][:,i+off])/t1m[m+1]
			
		if self.data['lasym']:	# asymmetric (not stellarator symmetric)
			LIST = ['rhomnc','rhomns']
			self.data['rhomns'] = zeros((self.data['ns'],self.data['mnmax']+off))

			# all n, m >= 2 and n > 0, m = 1
			for i in xrange(off+1, self.data['mnmax']-off):
				m = self.data['xm'][i]
				tnorm = 0.5/(t1m[m+1]**2 + t2m[m-1]**2)
				self.data['rhomns'][:,i] = tnorm*((self.data['rmns'][:,i+off] - self.data['zmnc'][:,i+off])*t1m[m+1]
													+ (self.data['rmns'][:,i-off] + self.data['zmnc'][:,i-off])*t2m[m-1])
			# all n, m = mpol-1, m = mpol
			for i in xrange(self.data['mnmax']-off, self.data['mnmax']+off):
				if i < self.data['mnmax']: m = self.data['mpol'] - 1
				else: m = self.data['mpol']
				tnorm = 0.5/(t1m[m+1]**2 + t2m[m-1]**2)
				self.data['rhomns'][:,i] = tnorm*(self.data['rmns'][:,i-off] + self.data['zmnc'][:,i-off])*t2m[m-1]

			# m = 0, n > 0 and m = 1, n <= 0
			for i in xrange(1, off+1):
				m = self.data['xm'][i]
				self.data['rhomns'][:,i] = 0.5*(self.data['rmns'][:,i+off] - self.data['zmnc'][:,i+off])/t1m[m+1]

			# n = 0, m = 0
			self.data['rhomns'][:,0] = 0 #0.5*(self.data['rmns'][:,off] - self.data['zmnc'][:,off])/t1m[1]

		for key in LIST:
			if self.data.has_key(key):
				mnmax = self.data[key].shape[1]
				self.spl[key] = [0] * mnmax	
				for i in xrange(mnmax):
					self.spl[key][i] = scinter.UnivariateSpline(self.data['s'], self.data[key][:,i], s = 0)


	# --- reverse_RZmn() ---
	# calculates Fourier modes for RZ from rho (= minor r) modes
	def reverse_RZmn(self):
		if not self.data.has_key('rhomnc'): self.get_rho_fourier()
		self.data['rrevmnc'] = zeros(self.data['rmnc'].shape)
		self.data['zrevmns'] = zeros(self.data['zmns'].shape)
		ntor = self.data['ntor']
		off = 2*ntor + 1
		t1m,t2m = self.rho_fourier_coeff()
		R0,Z0 = self.rho_fourier_center()
		LIST = ['rrevmnc','zrevmns']

		# n = 0, m = 0
		self.data['rrevmnc'][:,0] = R0 + t2m[0]*self.data['rhomnc'][:,off]
		# other modes
		for i in xrange(off, self.data['mnmax']):
			m = self.data['xm'][i]
			if self.data['rhomnc'][-1,m*off] != 0:
				self.data['rrevmnc'][:,i] += t1m[m]*self.data['rhomnc'][:,i-off]
				self.data['zrevmns'][:,i] += t1m[m]*self.data['rhomnc'][:,i-off]
		for i in xrange(1, self.data['mnmax']):
			m = self.data['xm'][i]
			if self.data['rhomnc'][-1,m*off] != 0:
				self.data['rrevmnc'][:,i] += t2m[m]*self.data['rhomnc'][:,i+off]
				self.data['zrevmns'][:,i] -= t2m[m]*self.data['rhomnc'][:,i+off]	
			if(m < self.data['mpol']-2): x = self.data['rhomnc'][-1,(m+2)*off]
			else: x = 0
			if(x == 0):
				tnorm = t1m[m+2]**2/t2m[m]
				self.data['rrevmnc'][:,i] += tnorm*self.data['rhomnc'][:,i+off]
				self.data['zrevmns'][:,i] -= tnorm*self.data['rhomnc'][:,i+off]

		if self.data['lasym']:	# asymmetric (not stellarator symmetric)
			LIST = ['rrevmnc','zrevmns','rrevmns','zrevmnc']
			self.data['rrevmns'] = zeros(self.data['rmns'].shape)
			self.data['zrevmnc'] = zeros(self.data['zmnc'].shape)
			# n = 0, m = 0
			self.data['zrevmnc'][:,0] = Z0 + t2m[0]*self.data['rhomns'][:,off]
			# other modes
			for i in xrange(off, self.data['mnmax']):
				m = self.data['xm'][i]
				if self.data['rhomnc'][-1,m*off] != 0:
					self.data['rrevmns'][:,i] += t1m[m]*self.data['rhomns'][:,i-off]
					self.data['zrevmnc'][:,i] -= t1m[m]*self.data['rhomns'][:,i-off]
			for i in xrange(1, self.data['mnmax']):
				m = self.data['xm'][i]
				if self.data['rhomnc'][-1,m*off] != 0:
					self.data['rrevmns'][:,i] += t2m[m]*self.data['rhomns'][:,i+off]
					self.data['zrevmnc'][:,i] += t2m[m]*self.data['rhomns'][:,i+off]
				if(m < self.data['mpol']-2): x = self.data['rhomnc'][-1,(m+2)*off]
				else: x = 0
				if(x == 0):
					tnorm = t1m[m+2]**2/t2m[m]
					self.data['rrevmns'][:,i] += tnorm*self.data['rhomns'][:,i+off]
					self.data['zrevmnc'][:,i] += tnorm*self.data['rhomns'][:,i+off]    

		for key in LIST:
			if self.data.has_key(key):
				mnmax = self.data[key].shape[1]
				self.spl[key] = [0] * mnmax	
				for i in xrange(mnmax):
					self.spl[key][i] = scinter.UnivariateSpline(self.data['s'], self.data[key][:,i], s = 0)
					

	# --- make_updwn_sym() ---
	# force the equilibrium to be up-down symmetric by setting rmns and zmnc to zero
	def make_updwn_sym(self):
		self.data['rmns'] = zeros(self.data['rmns'].shape)
		self.data['zmnc'] = zeros(self.data['zmnc'].shape)
		self.data['lasym'] = False
		LIST = ['rmns','zmnc']
		for key in LIST:
			if self.data.has_key(key):
				mnmax = self.data[key].shape[1]
				self.spl[key] = [0] * mnmax	
				for i in xrange(mnmax):
					self.spl[key][i] = scinter.UnivariateSpline(self.data['s'], self.data[key][:,i], s = 0)

	
	# --- write(Filename, mesh, vars) ----------------------------------------------------
	# writes Silo-File for 3-D visualization with Visit
	# Input:
	#		IMPORTANT: 3-D input arrays are in typical Python (C) order: shape = (nv, nu, ns)
	#				   but Visit needs Fortran order: shape = (ns, nu, nv), so arrays are transposed
	#				   in this function. Do not flip order yourself!
	#	Filename (string):	full PathName of file to be created; the .silo extension is recommended
	#	mesh (dict):		optional; contains the 'x', 'y' and 'z' arrays for the mesh & the surface coordinate 's'
	#						IMPORTANT: each array must be 3 dimensional
	#						if not given, then the default is self.get_plasmaSurface()
	#	vars (dict):		optional; contains the variables living on the mesh
	#						each key is one array, which gets written into the Silo file
	#						IMPORTANT: each array must be 3 dimensional, same shape as each mesh array
	#						if not given, then the default is self.get_j2D()
	# Output:
	#	Silo file written to disk, readable by Visit
	def write(self, Filename, mesh = None, vars = None):
		try: import pyvisfile.silo as silo
		except:
			print 'pyvisfile and silo not found. You must install them to use this feature. Abort!'
			return
		
		# create Silo file
		f = silo.SiloFile(Filename, mode = silo.DB_CLOBBER, target = silo.DB_LOCAL, 
						  fileinfo = 'B-field and current density on a toroidal surface')
						  
		# default if there is no input for mesh
		if (mesh is None):
			mesh = self.get_plasmaSurface()	# uses same s, u, v by default, as defined for vars below
				
		# write mesh
		coord = [mesh['x'].T, mesh['y'].T, mesh['z'].T]
	
		options = dict()
		options[silo.DBOPT_XLABEL] = 'R'
		options[silo.DBOPT_YLABEL] = 'R'
		options[silo.DBOPT_ZLABEL] = 'Z'
		options[silo.DBOPT_XUNITS] = 'm'
		options[silo.DBOPT_YUNITS] = 'm'
		options[silo.DBOPT_ZUNITS] = 'm'

		f.put_quadmesh('mesh', coord, coordtype = silo.DB_NONCOLLINEAR, optlist = options)
		
		# default if there is no input for vars
		if (vars is None):
			# 3-D mesh:
			nv, nu, ns = mesh['x'].shape
			s = mesh['s']
			u = linspace(0, 2*pi, nu)
			v = linspace(0, 2*pi, nv)
			if (ns == 1):
				if (isinstance(s, ndarray) | isinstance(s, list)): s = s[0]
				u, v = meshgrid(u, v)
			else:
				s, u, v = broadcast_arrays(s, u[:,newaxis], v[:,newaxis,newaxis]) # 3-D meshgrid
				s = s.flatten()
			u = u.flatten()
			v = v.flatten()
			
			# variables on the mesh	
			vars = self.get_j2D(s = s, u = u, v = v)
			
		
		# write vars
		if (vars.has_key('jees') | vars.has_key('bees')):	# vars is output from self.get_j2D() or self.get_B2D()
			options = dict()
				
			# check if vars arrays and mesh arrays are compatible
			if not (mesh['x'].size == vars['B2D'].size):
				raise RuntimeError('Mesh shape does not match vars shape')
			
			# B-field
			options[silo.DBOPT_UNITS] = 'T'
			Translate = {'Bpol':'Bpol2D', 'Btor':'Btor2D', 'Bmod':'B2D', 'BR':'BR2D', 'BZ':'BZ2D'}
			for name, key in Translate.iteritems():
				if vars.has_key(key):
					var = vars[key].reshape(mesh['x'].shape).T
					f.put_quadvar1(name, 'mesh', var, var.shape, centering = silo.DB_NODECENT, optlist = options)

			# current density
			options[silo.DBOPT_UNITS] = 'A/m^2'
			Translate = {'jpol':'jpol2D', 'jtor':'jtor2D', 'jmod':'j2D', 'jR':'jR2D', 'jZ':'jZ2D', 'jpar':'jpar2D'}
			for name, key in Translate.iteritems():
				if vars.has_key(key):
					var = vars[key].reshape(mesh['x'].shape).T
					f.put_quadvar1(name, 'mesh', var, var.shape, centering = silo.DB_NODECENT, optlist = options)
		else:
			# check if vars arrays and mesh arrays are compatible
			if not (mesh['x'].size == vars[vars.keys()[0]].size):
				raise RuntimeError('Mesh shape does not match vars shape')

			# anything
			for key in vars.keys():
				if not vars[key].flags.contiguous: vars[key] = vars[key].copy()
				f.put_quadvar1(key, 'mesh', vars[key], vars[key].shape, centering = silo.DB_NODECENT)
				  
		# close Silo file
		f.close()


	# --------------------------------------------------------------------------------
	# get poloidal angle from (R,Z) coordinates
	def __get_theta__(self, R, Z, phi, Raxis = None, Zaxis = None, n0only = False):
		# get magnetic axis
		if (Raxis is None) | (Zaxis is None):
			if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
			else: v = (phi/180.0)%2 * pi
			Raxis, Zaxis  = self.get_axis(v, n0only = n0only)
		
		Rm = R - Raxis # R relative to magnetic axis
		Zm = Z - Zaxis # Z relative to magnetic axis

		if isinstance(Rm, ndarray):
			Rm[(Rm < 1e-16) & (Rm >= 0)] = 1e-16
			Rm[(Rm > -1e-16) & (Rm < 0)] = -1e-16
		else:
			if (Rm < 1e-16) & (Rm >= 0): Rm = 1e-16
			if (Rm > -1e-16) & (Rm < 0): Rm = -1e-16

		theta = arctan(Zm/Rm)
		
		if isinstance(theta, ndarray):
			theta[Rm < 0] += pi
			theta[(Rm >= 0) & (Zm < 0)] += 2*pi
		else:
			if(Rm < 0): theta += pi
			if((Rm >= 0) & (Zm < 0)): theta += 2*pi

		return theta


	# --------------------------------------------------------------------------------
	# get minor radius from (R,Z) coordinates
	def __get_r__(self, R, Z, phi, Raxis = None, Zaxis = None, n0only = False):
		# get magnetic axis
		if (Raxis is None) | (Zaxis is None):
			if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
			else: v = (phi/180.0)%2 * pi
			Raxis, Zaxis  = self.get_axis(v, n0only = n0only)

		Rm = R - Raxis # R relative to magnetic axis
		Zm = Z - Zaxis # Z relative to magnetic axis

		return sqrt(Rm*Rm + Zm*Zm);
		
		
	# --------------------------------------------------------------------------------
	# get flux coordinates from cylindrical
	# s, u are initial guesses; if none: u = geometric angle, s preconditioned by bisection 
	def get_flux_coord(self, R, phi, Z, s_in = None, u_in = None, quiet = True, n0only = False):
		# setup initial guesses
		if self.data['machine'] == 'd3d': v = (2 - phi/180.0)%2 * pi		# angle back in radiants and right-handed
		else: v = (phi/180.0)%2 * pi
		if (u_in is None): u = self.__get_theta__(R, Z, phi, n0only = n0only)
		else: u = u_in.copy()
		if (s_in is None): s = self.__bisec__(R, Z, v, u, quiet = quiet, n0only = n0only)
		else: s = s_in.copy()
		
		s, u, success = self.__newton2D__(R, Z, s, u, v, quiet, n0only = n0only)
		
		if (len(s) == 1): s = s[0] 
		if (len(u) == 1): u = u[0] 
		
		if success: return s, u, v
		else:
			if not quiet: print 'no convergence' 
			return s, u, v
		
	# --------------------------------------------------------------------------------
	def write_to_mdsplus(self, shot, timeID, runnum_in = -1, cybele_user = None,
						 local_port = 8000, verbose = True, connection = None):
		"""
		Populate the MDS+ tree with the values from the wout file
		If runnum is not set, it will populate the first tree where USED=0
		
		2015/5 R Wilcox
		
		Inputs:
		  shot (int)          shot number
		  timeID (int)        either the time in ms, or if <100,
		                      the actual time-node number
		  runnum_in (int)     RUN-node number: 1-25
		                      or -1 = next free node (default)
		  mode (str)          Mode to open the MDS+ tree in (default = 'READ')
		  cybele_user (str)   User name to remotely log on to Cybele with
		                      to tunnel to r2d2
		  local_port (int)    Port number to use for tunnel on initiating
		                      system end (ie, swim/OLCF)
		Output:
		  datTree             MDS+ tree object
		"""
		import VMEC.Python.vmecdb as vdb
		
		if connection is None: datTree = vdb.connect_vmecdb(shot, 'EDIT', cybele_user, local_port)
		else: datTree = connection
		if datTree is None: return
		timenum, runnum = vdb.get_runnum(shot, timeID, runnum_in, verbose = False, check_used = False, connection = datTree)
		
		if timeID >= 100:
			# Write to TIME nodes
			timePathStr = 'TIME{:02}'.format(timenum)
			try:
				time_node = datTree.getNode(timePathStr + ':TIME')
				if not time_node.data() == timeID:
					if verbose: print "TIME gets overwritten. Was:", time_node.data(), "  Now:", timeID
					time_node.putData(timeID)
			except:
				if verbose: print "TIME did not get written into MDS+ database"
		
		pathStr = 'TIME{:02}'.format(timenum) + '.RUN{:02}'.format(runnum) + '.WOUT'
		if verbose:
			print '--------------------------------------'
			print "Writing WOUT file into VMECDB MDS+ node " + pathStr

		used_bool = datTree.getNode(pathStr + ':USED').getData()
		if used_bool and verbose: print '**Overwriting ' + pathStr + ', unwritten nodes may contain data from previous run'
		
		skip_nodes = ['input_extension','version_','lfreeb__logical__',
					  'lrfp__logical__','lasym__logical__','lrecon__logical__','used']
		
		subnodes = ['flags', 'parameters', '0D', '1D', 'spectral']
		
		for d in self.data.keys():
			if d in skip_nodes: continue
			doneBool = False
			for i in subnodes:
				try:
					node = datTree.getNode(pathStr + '.' + i + ':' +d)
					if d == 'curlabel':
						node.putData(array(self.data[d]))
					else:
						node.putData(self.data[d])
					doneBool = True
					break
				except: pass
			if (not doneBool) and verbose:
				print d + " did not get written into MDS+ database"
		
		# write outputs that do not line up with the node names
		try:
			node_ext = datTree.getNode(pathStr + '.FLAGS:INPUT_EXT')
			node_ext.putData(self.data['input_extension'])
		except:
			if verbose: print "INPUT_EXT did not get written into MDS+ database"
			
		try:
			node_version = datTree.getNode(pathStr + '.FLAGS:VERSION')
			node_version.putData(self.data['version_'])
		except:
			if verbose: print "VERSION did not get written into MDS+ database"
			
		try:
			node_used = datTree.getNode(pathStr[:-5] + ':USED')
			node_used.putData(1)
			node_used = datTree.getNode(pathStr+':USED')
			node_used.putData(1)
		except:
			if verbose: print "USED did not get written into MDS+ database"
			
	
	# --------------------------------------------------------------------------------
	def loadMDSplus(self, shot, timeID, runnum = 0, verbose = False):
		"""
		Load values from MDS+ tree into Wout object
		"""
		import VMEC.Python.vmecdb as vdb
		
		mdsDICT = vdb.read_from_mdsplus(shot, timeID, runnum, verbose)
		
		for k in mdsDICT['wout'].keys():
			self.data[k] = mdsDICT['wout'][k]
				
		# undo APHI adjustments, which are done again during loading
		LIST = ['bsubsmns', 'bsubsmnc', 'gmns', 'gmnc', 'equif', 'jcuru', 
				'jcurv', 'chipf', 'phipf', 'vp', 'phips']
		for key in LIST:
			if self.data.has_key(key):
				self.data[key] = self.data[key + '_raw'].copy()
				
		if self.data.has_key('used'): _=self.data.pop('used')
			
						
	# --------------------------------------------------------------------------------
	# get flux coordinates (s,u) from (R,Z) coordinates	using 2-D Newton method
	def __newton2D__(self, R0, Z0, s_in, u_in, v, quiet = True, n0only = False):
		imax = 100
		delta = 1e-14
		length_old = 1
		s = s_in.copy()
		u = u_in.copy()
	
		for i in xrange(imax):
	
			R, Z, J = self.__get_coordinates_and_jacobian__(s, u, v, n0only = n0only)
		
			fr = R - R0
			fz = Z - Z0
			
			# J[0] = dR/ds,  J[1] = dR/du,  J[2] = dZ/ds,  J[3] = dZ/du
			det = J[0] * J[3] - J[1] * J[2]
			ds = (J[3] * fr - J[1] * fz) / det
			du = (J[0] * fz - J[2] * fr) / det
	
			length = sqrt(ds**2 + du**2)
			if not quiet: print i, array([length]).max()
			if all(length < delta): return s, u, True	# convergence
			
			# if caught in a loop, relax delta somewhat and allow convergence
			if (array([length]).max() >= array([length_old]).max()) & (all(length < 10*delta)):
				return s, u, True	# convergence
	
			s -= ds
			u -= du
			u = u % (2*pi)	# u modulo 2pi
			
			length_old = length
		
		return s, u, False		# no convergence
	
	
	# --------------------------------------------------------------------------------
	# function called inside __newton2__
	# returns R, Z and the derivatives:
	# J[0] = dR/ds,  J[1] = dR/du,  J[2] = dZ/ds,  J[3] = dZ/du
	def __get_coordinates_and_jacobian__(self, s, u, v, n0only = False):
		ds = 0.0001
		rz = self.ev(['rmn','zmn'], s, u, v, derivs = True, n0only = n0only)
		rzds = self.ev(['rmn','zmn'], s+ds, u, v, n0only = n0only)
		
		# J[0] = dR/ds,  J[1] = dR/du,  J[2] = dZ/ds,  J[3] = dZ/du
		J = list(zeros(4))
		J[0] = (rzds['r'] - rz['r']) / ds
		J[1] = rz['drdu']
		J[2] = (rzds['z'] - rz['z']) / ds
		J[3] = rz['dzdu']
		
		return rz['r'], rz['z'], J


	# --------------------------------------------------------------------------------
	# bisection with only a few steps to find a crude estimate of s; use as preconditioner for newton2D
	# set u fixed as the geometric angle
	# assume f(a) < 0 & f(b) > 0 with f = (r(s,u,v) - Raxis)^2 + (z(s,u,v) - Zaxis)^2 - (r0 - Raxis)^2 - (z0 - Zaxis)^2
	def __bisec__(self, R0, Z0, v, u, a = 0, b = 1, quiet = True, Raxis = None, Zaxis = None, n0only = False):
		if (Raxis is None) | (Zaxis is None): Raxis, Zaxis = self.get_axis(v, n0only = n0only)
		r0 = (R0 - Raxis)**2 + (Z0 - Zaxis)**2
	
		if isinstance(u, ndarray):
			a = a*ones(u.shape)
			b = b*ones(u.shape)
	
		for i in xrange(10):	# after N steps is (b-a) = 0.5**N; here N = 10 -> (b-a) < 1e-3
			s = (a + b)/2.0

			rz = self.ev(['rmn','zmn'], s, u, v, n0only = n0only)
			r = (rz['r'] - Raxis)**2 + (rz['z'] - Zaxis)**2
			f = r - r0
		
			if isinstance(u, ndarray):
				b[f > 0] = s[f > 0]
				a[f < 0] = s[f < 0]
			else:
				if(f > 0): b = s
				else: a = s
		
			if not quiet: print i, array([abs(b-a)]).max()
		return s

	
# ----------------------------------------------------------------------------------------
# --- End of Class -----------------------------------------------------------------------

	
	
	
	
	
	
	
	
	
	
	


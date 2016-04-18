import sys,os
import numpy as np
import h5py
import scipy.interpolate as interp
from netCDF4 import Dataset
import pylab as plt

import VMEC.Python.wout_class as WC
from Misc.deriv import deriv5pt

# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------

class chk_inside_class:
	"""
	class to determine if (R,phi,Z) is inside s = 1 surface or outside	
	phi must be a scalar, R,Z can be arrays
	"""
	def __init__(self, W, phi):
		"initialize: get spline of r_s(theta_s) for s = 1"
		# wout file
		self.W = W
		self.data = self.W.data
		self.phi = phi
		
		# s = 1 surface
		Rs, Zs = self.W.get_surface(phi = (2-phi/np.pi)*180) # angle in degree and left-handed
		
		# axis
		self.Raxis, self.Zaxis  = self.W.get_axis(phi)

		r = self.W.__get_r__(Rs, Zs, None, self.Raxis, self.Zaxis)
		theta = self.W.__get_theta__(Rs, Zs, None, self.Raxis, self.Zaxis)
		
		# make theta periodic & monotonic
		idx = np.where(np.abs(np.diff(theta)) > 5)[0]	# locate 2pi jump in theta
		theta[idx + 1 ::] += 2*np.pi				# make theta increase monotonically
		if (theta[0] > 0):
			theta = np.append(theta - 2*np.pi, theta[1:idx + 3])
			r = np.append(r, r[1:idx + 3])
		else: 
			theta = np.append(theta, theta[1:idx + 3] + 2*np.pi)
			r = np.append(r, r[1:idx + 3])
		
		# define function
		self.r_func = interp.UnivariateSpline(theta, r, s = 0)
	
	# R,Z can be scalars or arrays
	def __call__(self, R, Z):
		"""
		return inside
		evaluate point(s) R,Z; returns (list of) boolean
		"""
		r = self.W.__get_r__(R, Z, None, self.Raxis, self.Zaxis)
		theta = self.W.__get_theta__(R, Z, None, self.Raxis, self.Zaxis)
		
		rs = self.r_func(theta)
		
		return r <= rs
		
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------

class xpand_class:
	"""
	Read magnetic field from xpand_mpi on a evenly spaced 3-D grid (R, phi, Z)

	How to: Initialize
	   import use_xpand as ux
	   xpand = ux.xpand_class(wout)

	To write points.dat file for xpand run
		xpand.write_points()

	To start, read results:
	   xpand.load(file)

	Show fields for angle self.phi[k,0,0] and mark miscalculated grid points.
	Run several times to find the proper eps that identifies all the bad grid points and no more.
	Then heal them
	   xpand.view(k, eps = 1)
	   xpand.heal_all(eps = 1)   # Heal all bad grid points

	Write results back to file for distribution
	   xpand.wite(filename)
	"""
	def __init__(self, wout, Rmin = 1, Rmax = 2.4, NR = 2, Zmin = -1.5, Zmax = 1.5, 
				 NZ = 2, pmin = 0, pmax = 2*np.pi, Np = 1):
		"initialize & set 3-D grid"
		self.W = WC.Wout(wout)
		self.Rmin = Rmin
		self.Rmax = Rmax
		self.Zmin = Zmin
		self.Zmax = Zmax
		self.pmin = pmin
		self.pmax = pmax
		self.NR = NR
		self.NZ = NZ
		self.Np = Np
		
		self.R = np.zeros((Np, NZ, NR))
		self.phi = np.zeros(self.R.shape)
		self.Z = np.zeros(self.R.shape)
		self.BR = np.zeros(self.R.shape)
		self.Bphi = np.zeros(self.R.shape)
		self.BZ = np.zeros(self.R.shape)
		self.Pressure = np.zeros(self.R.shape)
		
		self.inside_all = []
		
		self.PresFunc = interp.UnivariateSpline(self.W.data['s'], self.W.data['presf'], s = 0)
		
		# get working dir from wout-file name
		idx = wout[::-1].find('/')	# returns location of last '/' in wout or -1 if not found
		if(idx == -1):
			self.cwd = os.getcwd()
		else:
			idx *= -1
			self.cwd = os.path.abspath(wout[0:idx - 1])	# absolute path without a final '/'

		# set grid
		self._init_grid()
		
		
	def _init_grid(self):
		"""
		return R, phi, Z
		sets equidistand 3-D grid; phi is in radiants and right-handed;
		calls chk_inside_class to id points inside s = 1 or not
		"""
		# set grid
		R = np.linspace(self.Rmin, self.Rmax, self.NR)
		Z = np.linspace(self.Zmin, self.Zmax, self.NZ)
		
		if (self.Np == 1):
			R, Z = np.meshgrid(R, Z)
			phi = self.pmin*np.ones((self.NZ, self.NR))
			self.R[0,:,:] = R
			self.phi[0,:,:] = phi
			self.Z[0,:,:]= Z
		else:
			if(self.pmin % (2*pi) == self.pmax % (2*pi)): # periodic
				phi = np.linspace(self.pmin, self.pmax, self.Np+1)[0:-1]	# last = first -> remove last
			else: phi = np.linspace(self.pmin, self.pmax, self.Np)
			R, Z, phi = np.broadcast_arrays(R, Z[:,np.newaxis], phi[:,np.newaxis,np.newaxis]) # 3-D meshgrid
			self.R = R
			self.phi = phi
			self.Z = Z
	
		for k in xrange(self.Np):
			r, z = self.R[k,:,:].flatten(), self.Z[k,:,:].flatten()
			chk_inside = chk_inside_class(self.W, self.phi[k,0,0])
			inside = chk_inside(r,z)
			self.inside_all.append(inside)			
	
	
	def _heal_spots(self, B, eps = 1, index = None, quiet = True, index_only = False):
		"""
		return B, index
		finds miscalculated grid points and fixes them by averaging over surrounding values 
		"""
		B = B.copy()
		Ny, Nx = B.shape

		inside = lambda i,j: False if ((i < 0) | (i >= Ny) | (j < 0) | (j >= Nx)) else True
		good = lambda i,j: False if ((i,j) in index) else True
		
		LIST = [(1,0), (1,1), (1,-1), (-1,0), (-1,1), (-1,-1), (0,1), (0,-1)]
		
		if (index == None):
			index = []
			for i in xrange(Ny):
				for j in xrange(Nx):
					A = np.array([])
					for m,n in LIST:
						if inside(i+m,j+n): A = np.append(A,B[i+m,j+n])			
					if (abs(B[i,j] - np.mean(A)) > eps):
						index.append((i,j))
						
		if not index_only:
			for i,j in index:
				A = np.array([])
				for m,n in LIST:
					for k in xrange(1,4):
						if good(i+k*m,j+k*n):
							if inside(i+k*m,j+k*n): A = np.append(A,B[i+k*m,j+k*n])
							break;
				B[i,j] = np.mean(A)

# 		if (index == None):
# 			index = []
# 
# 			for i in xrange(Ny):
# 				for j in xrange(Nx):
# 					# bottom bndy
# 					if (i == 0) & (j > 0) & (j < Nx-1):
# 						A = array([B[i+1,j], B[i+1,j+1], B[i+1,j-1], B[i,j+1], B[i,j-1]])
# 					# left bndy
# 					elif (j == 0) & (i > 0) & (i < Ny-1):
# 						A = array([B[i+1,j], B[i+1,j+1], B[i-1,j], B[i-1,j+1], B[i,j+1]])
# 					# top bndy
# 					elif (i == Ny-1) & (j > 0) & (j < Nx-1):
# 						A = array([B[i-1,j], B[i-1,j+1], B[i-1,j-1], B[i,j+1], B[i,j-1]])
# 					# right bndy
# 					elif (j == Nx-1) & (i > 0) & (i < Ny-1):
# 						A = array([B[i+1,j], B[i+1,j-1], B[i-1,j], B[i-1,j-1], B[i,j-1]])
# 					# lower left corner
# 					elif (i == 0) & (j == 0):
# 						A = array([B[i+1,j], B[i+1,j+1], B[i,j+1]])
# 					# lower right corner
# 					elif (i == 0) & (j == Nx-1):
# 						A = array([B[i+1,j], B[i+1,j-1], B[i,j-1]])
# 					# top left corner
# 					elif (i == Ny-1) & (j == 0):
# 						A = array([B[i-1,j], B[i-1,j+1], B[i,j+1]])
# 					# top right corner
# 					elif (i == Ny-1) & (j == Nx-1):
# 						A = array([B[i-1,j], B[i-1,j-1], B[i,j-1]])
# 					# inside
# 					else:
# 						A = array([B[i+1,j], B[i+1,j+1], B[i+1,j-1], B[i-1,j], B[i-1,j+1], B[i-1,j-1], B[i,j+1], B[i,j-1]])
# 					
# 					meanA = mean(A)
# 					if (abs(B[i,j] - meanA) > eps):
# 						index.append((i,j))
# 						B[i,j] = meanA
# 
# 		else:
# 			for i,j in index:
# 				if (i == 0) & (j > 0) & (j < Nx-1):
# 					A = array([B[i+1,j], B[i+1,j+1], B[i+1,j-1], B[i,j+1], B[i,j-1]])
# 				# left bndy
# 				elif (j == 0) & (i > 0) & (i < Ny-1):
# 					A = array([B[i+1,j], B[i+1,j+1], B[i-1,j], B[i-1,j+1], B[i,j+1]])
# 				# top bndy
# 				elif (i == Ny-1) & (j > 0) & (j < Nx-1):
# 					A = array([B[i-1,j], B[i-1,j+1], B[i-1,j-1], B[i,j+1], B[i,j-1]])
# 				# right bndy
# 				elif (j == Nx-1) & (i > 0) & (i < Ny-1):
# 					A = array([B[i+1,j], B[i+1,j-1], B[i-1,j], B[i-1,j-1], B[i,j-1]])
# 				# lower left corner
# 				elif (i == 0) & (j == 0):
# 					A = array([B[i+1,j], B[i+1,j+1], B[i,j+1]])
# 				# lower right corner
# 				elif (i == 0) & (j == Nx-1):
# 					A = array([B[i+1,j], B[i+1,j-1], B[i,j-1]])
# 				# top left corner
# 				elif (i == Ny-1) & (j == 0):
# 					A = array([B[i-1,j], B[i-1,j+1], B[i,j+1]])
# 				# top right corner
# 				elif (i == Ny-1) & (j == Nx-1):
# 					A = array([B[i-1,j], B[i-1,j-1], B[i,j-1]])
# 				# inside
# 				else:
# 					A = array([B[i+1,j], B[i+1,j+1], B[i+1,j-1], B[i-1,j], B[i-1,j+1], B[i-1,j-1], B[i,j+1], B[i,j-1]])
# 				B[i,j] = mean(A)
# 			

		if not quiet: print 'Found:', len(index), 'Points'
		return B, index

			
	def heal_all(self, eps = 1, quiet = False):
		"fix bad grid points for all fields and all angles; see self.heal()"
		for k in xrange(self.Np):
			self.heal(k, eps, quiet)


	def heal(self, k, eps = 1, quiet = False):
		"""
		fix bad grid points for all fields at angle k; eps is the threshold to find bad
		points in the Bphi field (most reliable); use self.view() with various eps to find 
		a threshold that works for all; the default is quite good in generic cases
		"""
		_,index = self._heal_spots(self.Bphi[k,:,:], eps = eps, index_only = True)
		if len(index) > 0:
			if not quiet: print 'k =', k, '\t Angle =', self.phi[k,0,0], ' \t',
			self.BR[k,:,:],_ = self._heal_spots(self.BR[k,:,:], index = index)
			self.Bphi[k,:,:],_ = self._heal_spots(self.Bphi[k,:,:], index = index, quiet = quiet)
			self.BZ[k,:,:],_ = self._heal_spots(self.BZ[k,:,:], index = index)
			
			
	def write_points(self, filename = 'points.dat'):
		"creates points.dat file in cwd"
		R = self.R.flatten()
		phi = self.phi.flatten()
		Z = self.Z.flatten()
		N = len(R)
		
		with open(self.cwd + '/' + filename,'w') as f:
			f.write("# NR = " + str(self.NR) + "\n")
			f.write("# Np = " + str(self.Np) + "\n")
			f.write("# NZ = " + str(self.NZ) + "\n")
			for i in xrange(N): 
				f.write(str(R[i]) + "\t" + str(phi[i]) + "\t" + str(Z[i]) + "\n")
	

	def view_all(self, spots = True, eps = 1, Bmod_only = True, quiet = True):
		"shows contour plots for all angles; see self.view()"
		for k in xrange(self.Np):
			self.view(k, spots, eps, Bmod_only)


	def view(self, k, spots = True, eps = 1, Bmod_only = False, BRrange = 0.6, BZrange = 0.8, Bthrange = 0.7, quiet = True, BRBZ = False):
		"""
		filled contour plots of BR, Bphi & BZ for angle phi[k]; marks bad grid points
		scan through multiple eps for best results; see self.heal_all() for more info on eps
		"""
		title_string = str(k) + ': Angle = ' + str(round(self.phi[k,0,0],3))
		Bmod = np.sqrt(self.BR[k,:,:]**2 + self.Bphi[k,:,:]**2 + self.BZ[k,:,:]**2)
		if spots: 
			_,index = self._heal_spots(self.Bphi[k,:,:], eps = eps, index_only = True, quiet = quiet)
		
		if not Bmod_only:
			plt.rcParams['font.size'] = 18
			plt.rcParams['font.family'] = 'Arial'
			plt.figure(figsize = (7,9)); plt.subplot(111, aspect = 'equal')
			cs = plt.contourf(self.R[k,:,:], self.Z[k,:,:], np.sqrt(self.BR[k,:,:]**2 + self.BZ[k,:,:]**2), np.linspace(0,Bthrange,128))
			plt.xlabel('R [m]'); plt.ylabel('Z [m]'); plt.title(title_string, size = 18)
			plt.xlim(self.Rmin, self.Rmax); plt.ylim(self.Zmin, self.Zmax)
			C = plt.colorbar(cs, pad = 0.01, extend = 'neither', format = '%.2g')
			C.set_label('B$_{\\theta}$ [T]', rotation = 270, va = 'bottom')
			if spots: 
				for i,j in index: 
					plt.plot(self.R[k,i,j], self.Z[k,i,j], 'ko', mfc = 'none')
		
			plt.figure(figsize = (7,9)); plt.subplot(111, aspect = 'equal')
			cs = plt.contourf(self.R[k,:,:], self.Z[k,:,:], self.Bphi[k,:,:], 128)
			plt.xlabel('R [m]'); plt.ylabel('Z [m]'); plt.title(title_string, size = 18)
			plt.xlim(self.Rmin, self.Rmax); plt.ylim(self.Zmin, self.Zmax)
			C = plt.colorbar(cs, pad = 0.01, extend = 'neither', format = '%.2g')
			C.set_label('B$_\\phi$ [T]', rotation = 270, va = 'bottom')
			if spots: 
				for i,j in index: 
					plt.plot(self.R[k,i,j], self.Z[k,i,j], 'ko', mfc = 'none')
		
			if BRBZ:
				plt.rcParams['font.size'] = 18
				plt.rcParams['font.family'] = 'Arial'
				plt.figure(figsize = (7,9)); plt.subplot(111, aspect = 'equal')
				cs = plt.contourf(self.R[k,:,:], self.Z[k,:,:], self.BR[k,:,:], np.linspace(-BRrange,BRrange,128))
				plt.xlabel('R [m]'); plt.ylabel('Z [m]'); plt.title(title_string, size = 18)
				plt.xlim(self.Rmin, self.Rmax); plt.ylim(self.Zmin, self.Zmax)
				C = plt.colorbar(cs, pad = 0.01, extend = 'neither', format = '%.2g')
				C.set_label('B$_R$ [T]', rotation = 270, va = 'bottom')
				if spots: 
					for i,j in index: 
						plt.plot(self.R[k,i,j], self.Z[k,i,j], 'ko', mfc = 'none')

				plt.figure(figsize = (7,9)); plt.subplot(111, aspect = 'equal')
				cs = plt.contourf(self.R[k,:,:], self.Z[k,:,:], self.BZ[k,:,:], np.linspace(-BZrange,BZrange,128))
				plt.xlabel('R [m]'); plt.ylabel('Z [m]'); plt.title(title_string, size = 18)
				plt.xlim(self.Rmin, self.Rmax); plt.ylim(self.Zmin, self.Zmax)
				C = plt.colorbar(cs, pad = 0.01, extend = 'neither', format = '%.2g')
				C.set_label('B$_Z$ [T]', rotation = 270, va = 'bottom')
				if spots: 
					for i,j in index: 
						plt.plot(self.R[k,i,j], self.Z[k,i,j], 'ko', mfc = 'none')
		else:
			plt.figure(figsize = (7,9)); plt.subplot(111, aspect = 'equal')
			cs = plt.contourf(self.R[k,:,:], self.Z[k,:,:], Bmod, np.linspace(0,4,128), cmap = plt.cm.prism)
			plt.xlabel('R [m]'); plt.ylabel('Z [m]'); plt.title(title_string, size = 18)
			plt.xlim(self.Rmin, self.Rmax); plt.ylim(self.Zmin, self.Zmax)
			#C = plt.colorbar(cs, pad = 0.01, extend = 'neither', format = '%.2g')
			#C.set_label('|B| [T]', rotation = 270, va = 'bottom')
			if spots: 
				for i,j in index: 
					plt.plot(self.R[k,i,j], self.Z[k,i,j], 'ko', mfc = 'none')
		
	
	def write(self, filename):
		"save R, phi, Z, BR, Bphi, BZ in file; '.h5' creates hdf5; ascii otherwise"
		# hdf5 file
		if '.h5' in filename:
			# reverse phi direction to make it right handed increasing
			R = self.R
			v = self.phi
			Z = self.Z
			BR = self.BR
			Bphi = self.Bphi
			BZ = self.BZ
			Pressure = self.Pressure

			# if array is not c-contiguous, make a copy
			if not R.flags.contiguous: R = R.copy()
			if not v.flags.contiguous: v = v.copy()
			if not Z.flags.contiguous: Z = Z.copy()
			if not BR.flags.contiguous: BR = BR.copy()
			if not Bphi.flags.contiguous: Bphi = Bphi.copy()
			if not BZ.flags.contiguous: BZ = BZ.copy()
			if not Pressure.flags.contiguous: Pressure = Pressure.copy()
		
			with h5py.File(filename, 'w') as f:
				dset = f.create_dataset('R', R.shape, R.dtype, maxshape = (None,None,None))
				dset.write_direct(R)		
				dset = f.create_dataset('phi', v.shape, v.dtype, maxshape = (None,None,None))
				dset.write_direct(v)		
				dset = f.create_dataset('Z', Z.shape, Z.dtype, maxshape = (None,None,None))
				dset.write_direct(Z)		
				dset = f.create_dataset('BR', BR.shape, BR.dtype, maxshape = (None,None,None))
				dset.write_direct(BR)		
				dset = f.create_dataset('Bphi', Bphi.shape, Bphi.dtype, maxshape = (None,None,None))
				dset.write_direct(Bphi)		
				dset = f.create_dataset('BZ', BZ.shape, BZ.dtype, maxshape = (None,None,None))
				dset.write_direct(BZ)		
				dset = f.create_dataset('Pressure', Pressure.shape, Pressure.dtype, maxshape = (None,None,None))
				dset.write_direct(Pressure)		
		
		# ascii file
		else:
			R = self.R.flatten()
			v = self.phi.flatten()
			Z = self.Z.flatten()
			BR = self.BR.flatten()
			Bphi = self.Bphi.flatten()
			BZ = self.BZ.flatten()
			Pressure = self.Pressure.flatten()
			
			fmt = ' 21.15e'		# output format
			with open(filename,'w') as f:
				f.write('# NR \t Nphi \t NZ\n')
				f.write('# ' + str(self.NR) + '\t' + str(self.Np) + '\t' + str(self.NZ) + '\n')
				f.write('# R[m]               \t phi[rad]             \t Z[m]                 \t BR[T]                \t Bphi[T]              \t BZ[T]              \t Pres[Pa]\n')
				for i in xrange(len(R)): 
					f.write(format(R[i],fmt) + "\t" + format(v[i],fmt) + "\t" + format(Z[i],fmt) + "\t" + format(BR[i],fmt) + "\t" + format(Bphi[i],fmt) + "\t" + format(BZ[i],fmt) + "\t" + format(Pressure[i],fmt) + "\n")
				
	
	def load(self, filename = 'xpand.dat'):
		"load R, phi, Z, BR, Bphi, BZ from file"
		if '.h5' in filename:
			with h5py.File(filename, 'r') as f:
				dset = f['R']; R = np.zeros((dset.shape)); dset.read_direct(R);
				dset = f['phi']; v = np.zeros((dset.shape)); dset.read_direct(v);
				dset = f['Z']; Z = np.zeros((dset.shape)); dset.read_direct(Z);
				dset = f['BR']; BR = np.zeros((dset.shape)); dset.read_direct(BR);
				dset = f['Bphi']; Bphi = np.zeros((dset.shape)); dset.read_direct(Bphi);
				dset = f['BZ']; BZ = np.zeros((dset.shape)); dset.read_direct(BZ);
				dset = f['Pressure']; Pressure = np.zeros((dset.shape)); dset.read_direct(Pressure);
			
			self.Np, self.NZ, self.NR = R.shape
		else:
			data = np.loadtxt(filename)
			NRNZ = len(data[data[:,1] == data[0,1],1])
			NRNP = len(data[data[:,2] == data[0,2],2])
			NPNZ = len(data[data[:,0] == data[0,0],0])
			self.NZ = int(round((NRNZ * NPNZ / float(NRNP))**0.5))
			self.NR = NRNZ/self.NZ
			self.Np = NPNZ/self.NZ
			print 'NR =', self.NR, ', Nphi =', self.Np, ', NZ =', self.NZ

			R = data[:,0].reshape(self.Np, self.NZ, self.NR)
			v = data[:,1].reshape(self.Np, self.NZ, self.NR)
			Z = data[:,2].reshape(self.Np, self.NZ, self.NR)
			BR = data[:,3].reshape(self.Np, self.NZ, self.NR)
			Bphi = data[:,4].reshape(self.Np, self.NZ, self.NR)
			BZ = data[:,5].reshape(self.Np, self.NZ, self.NR)
			Pressure = data[:,6].reshape(self.Np, self.NZ, self.NR)
			
		self.R = R
		self.phi = v
		self.Z = Z
		self.BR = BR
		self.Bphi = Bphi
		self.BZ = BZ
		self.Pressure = Pressure
		
		self.Rmin = self.R.min()
		self.Rmax = self.R.max()
		self.Zmin = self.Z.min()
		self.Zmax = self.Z.max()
		self.pmin = self.phi.min()
		self.pmax = self.phi.max()
		
		self.inside_all = []
		for k in xrange(self.Np):
			r, p, z = self.R[k,:,:].flatten(), self.phi[k,0,0], self.Z[k,:,:].flatten()	
			chk_inside = chk_inside_class(self.W, p)
			inside = chk_inside(r,z)
			self.inside_all.append(inside)


	def div(self, fullGrid = False):
		"""
		computes the divergence of B in cylindrical coordinates
		div B = BR/R + dBR/dR + 1/R dBphi/dphi + dBZ/dZ
		"""
		Rh = 0.5*(self.R[:,1::,1::] + self.R[:,1::,0:-1])
		Zh = 0.5*(self.Z[:,1::,1::] + self.Z[:,0:-1,1::])
		BR_Zh = 0.5*(self.BR[:,1::,:] + self.BR[:,0:-1,:])
		BZ_Rh = 0.5*(self.BZ[:,:,1::] + self.BZ[:,:,0:-1])
		dR = (self.Rmax - self.Rmin)/float(self.NR-1)
		dZ = (self.Zmax - self.Zmin)/float(self.NZ-1)
		
		divh = ((self.R[:,1::,1::]*BR_Zh[:,:,1::] - self.R[:,1::,0:-1]*BR_Zh[:,:,0:-1])/Rh/dR 
				+ (BZ_Rh[:,1::,:] - BZ_Rh[:,0:-1,:])/dZ)
		
		if self.R.shape[0] > 1:		# else: only 2D divergence with dBp = 0
			Bp_half_m1 = np.zeros(divh.shape)
			Bp_half_p1 = np.zeros(divh.shape)
			dphi = (self.phi.max() - self.phi.min())/float(self.Np-1)
			Bp_half_m1[0,:,:] = 0.25*(self.Bphi[-1,1::,0:-1] + self.Bphi[-1,0:-1,0:-1] + self.Bphi[-1,1::,1::] + self.Bphi[-1,0:-1,1::])
			Bp_half_m1[1::,:,:] = 0.25*(self.Bphi[0:-1,1::,0:-1] + self.Bphi[0:-1,0:-1,0:-1] + self.Bphi[0:-1,1::,1::] + self.Bphi[0:-1,0:-1,1::])
			Bp_half_p1[-1,:,:] = 0.25*(self.Bphi[0,1::,0:-1] + self.Bphi[0,0:-1,0:-1] + self.Bphi[0,1::,1::] + self.Bphi[0,0:-1,1::])
			Bp_half_p1[0:-1,:,:] = 0.25*(self.Bphi[1::,1::,0:-1] + self.Bphi[1::,0:-1,0:-1] + self.Bphi[1::,1::,1::] + self.Bphi[1::,0:-1,1::])

			divh += 0.5*(Bp_half_p1 - Bp_half_m1)/dphi/Rh;
						
		if fullGrid:
			from scipy.interpolate import griddata
			divB = np.zeros(self.R.shape)
			for k in xrange(self.Np):
				divB[k,:,:] = griddata((Rh[k,:,:].flatten(), Zh[k,:,:].flatten()), divh[k,:,:].flatten(), (self.R[0,:,:], self.Z[0,:,:]), method='linear', fill_value = 0)
			return divB
		return Rh, Zh, divh
		
		
	def plot_div(self, k = 0):
		title_string = str(k) + ': Angle = ' + str(round(self.phi[k,0,0],3))
		R,Z,divB = self.div()
		plt.figure(figsize = (7,9))
		plt.contourf(R[k,:,:], Z[k,:,:], np.log10(np.abs(divB[k,:,:])), np.linspace(-14,0,128), extend = 'both')
		plt.axes().set_aspect('equal')
		C = plt.colorbar(pad = 0.01, extend = 'both', ticks = np.arange(-14,1))
		C.set_label('$\\log_{10}(|\\nabla\\cdot B|)$', rotation = 270, va = 'bottom')
		plt.xlabel('R [m]')
		plt.ylabel('Z [m]')
		plt.title(title_string, fontsize = 18)
		
	
	def _init_vacuumBfield(self, mgrid = None):
		""" 
		constructs vacuum magnetic field interpolation; reads mgrid file, or 'mgrid_file' from wout 
		returns lists BRvac_func[k], Bphivac_func[k], BZvac_func[k], k = 0,...,kp-1, 
		which have interpolation functions for each toroidal plane k
		"""
		# read mgrid and set grid
		if mgrid == None: mgrid = self.W.data['mgrid_file']
		print 'Read MGRID file:', mgrid
		data = openNETCDF(mgrid)
		R = np.linspace(data['rmin'], data['rmax'], data['ir'])
		Z = np.linspace(data['zmin'], data['zmax'], data['jz'])
		dphi = 2*np.pi/data['nfp'] / data['kp']
		self.mgridNp = data['kp']
	
		# all fields from the coils
		extcur = self.W.data['extcur']
		BRmgrid, Bphimgrid, BZmgrid = 0, 0, 0
		for i in xrange(self.W.data['nextcur']):
			BRmgrid += extcur[i] * data['br_'+format(i+1,'03d')]
			Bphimgrid += extcur[i] * data['bp_'+format(i+1,'03d')]
			BZmgrid += extcur[i] * data['bz_'+format(i+1,'03d')]
			
		# rz depends on u and v; rz & v are flat!
		self.BRvac_func = []
		self.Bphivac_func = []
		self.BZvac_func = []
		
		for i in xrange(data['kp']):
			self.BRvac_func.append(interp.RectBivariateSpline(R, Z, BRmgrid[i,:,:].T))
			self.Bphivac_func.append(interp.RectBivariateSpline(R, Z, Bphimgrid[i,:,:].T))
			self.BZvac_func.append(interp.RectBivariateSpline(R, Z, BZmgrid[i,:,:].T))

	
	def get_vacuumB(self, R, phi, Z):
		""" 
		Computes the vacuum magnetic field BR, Bphi, BZ at R, phi, Z location(s).
		Assumes R,Z are flat arrays or scalar, and phi is scalar
		"""
		if not hasattr(self, 'BRvac_func'): self._init_vacuumBfield()
		dphi = 2*np.pi / self.mgridNp
		k = int(round(phi/dphi)) % self.mgridNp	# nearest neighbor approximation
		BR = self.BRvac_func[k].ev(R,Z)
		Bphi = self.Bphivac_func[k].ev(R,Z)
		BZ = self.BZvac_func[k].ev(R,Z)
		return BR, Bphi, BZ
		
		
	def vacuumB(self):
		"""Computes vacuum B-field for the entire grid"""
		BRvac, Bpvac, BZvac = np.zeros(self.R.shape), np.zeros(self.R.shape), np.zeros(self.R.shape)
		for k in xrange(self.Np):
			br, bp, bz = self.get_vacuumB(self.R[k,:,:].flatten(),self.phi[k,0,0],self.Z[k,:,:].flatten())
			BRvac[k,:,:] = br.reshape(self.NZ, self.NR)
			Bpvac[k,:,:] = bp.reshape(self.NZ, self.NR)
			BZvac[k,:,:] = bz.reshape(self.NZ, self.NR)
		return BRvac, Bpvac, BZvac

	
	
# ----------------------------------------------------------------------------------------
# --- End of Class -----------------------------------------------------------------------

def openNETCDF(Filename):
	"""Load netCDF file into dictionary """
	DATA = {}
	from warnings import filterwarnings
	filterwarnings('ignore')
	data = Dataset(Filename, 'r')
	varnames = data.variables.keys()
	for n in varnames:
		if(data.variables[n].size > 1):										# arrays (all types)
			DATA[n] = np.array(data.variables[n][:])
			if(data.variables[n].dtype == 'S1'):							# character array
				if(data.variables[n].size == data.variables[n].shape[0]):	# 1-D						
					DATA[n] = ''.join(data.variables[n][:]).strip()
				else:														# 2-D
					DATA[n] = []
					for i in xrange(data.variables[n].shape[0]):
						DATA[n].append(''.join(data.variables[n][i,:]).strip())
		else:																# single variable
			if(data.variables[n].dtype == 'float'):							# float
				DATA[n] = np.float64(data.variables[n][:])
			elif(data.variables[n].dtype == 'int32'):						# int
				DATA[n] = int(data.variables[n][:])
			elif(data.variables[n].dtype == 'S1'):							# char
				try: DATA[n] = ''.join(data.variables[n][:])	# at fixed bndy: mgrid_mode exists but is masked -> error
				except: DATA[n] = data.variables[n][:]
			else:															# unknown
				print 'unknown datatype in variable:', n
				DATA[n] = data.variables[n][:]
	return DATA

	
	
	
	
	
	
	
	
	
	
	
	
			

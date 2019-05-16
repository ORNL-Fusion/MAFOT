"""
Code draws histogramms for Connection length and Penetration depth
specifically designed for MAFOT footprint data

You can eiter read in one or two footprint files
then you get both histogramms

useLogScale aktivates a logarithmic scale for better visibility
to avoid negative numbers, log(data + 1) is ploted
the axis tick labes are manipulated so that the values are correct
"""

import numpy as np

def makeFLF(file, file2 = None, useLogScale = False, M = 200, SOL_Range_Lc = 0.075, SOL_Range_psi = 0.997,
             Lc_cutoff = 2.6, Lc_ft = 0.13, Lc_long = 0.4, correct_psi = True, save = False, show = True,
             flfLim = None): 
	"""
	Input
	 file            filename of footprint
	 file2           None, if not wanted
	 useLogScale     0: linear y-Axis   1: log(y+1) y-Axis
	 M               Number of bars
	 SOL_Range_Lc    Minimum connection length; D3D = 0.075; ITER = 0.22; NSTX = 0.046
	 SOL_Range_psi   Maximum penetration depth; D3D = 0.997; ITER = 0.997; NSTX = 0.999
	 Lc_cutoff       Cutoff Lc x-Axis at value [km], if value > max(Lx): no cutoff
	 Lc_ft           max Flux Tube con. Length = 2pol. Turns; D3D = 0.13; ITER = 0.44; NSTX = 0.085
	 Lc_long         another con. Length; D3D = 0.4; ITER = 1; NSTX = 0.2
	 correct_psi     set psi = 1 for Lc < SOL_Range_Lc
	"""

	idx = file[::-1].find('/')
	if(idx == -1): tag = file[5:-4] 
	else: tag = file[-idx+5:-4]
	tag = tag.replace('.','_')
	if useLogScale: tag += '_log'
	if file2 is not None: 
		if tag[0:3] == '_in': tag = tag.replace('_in','_all')
		if tag[0:3] == '_ou': tag = tag.replace('_out','_all')
		if tag[0:3] == '_sh': tag = tag.replace('_shelf','_all')
	
	Lc, psimin = readData(file, file2, SOL_Range_Lc, correct_psi)
	_,_, psi_min, flf98, flf95 = psiHist(psimin, M, SOL_Range_psi, useLogScale, save, tag, show, flfLim)
	_,_, Lc_average, flfft, flflong = LcHist(Lc, M, SOL_Range_Lc, Lc_ft, Lc_long, Lc_cutoff, useLogScale, save, tag, show)
	
	return psi_min, flf98, flf95, Lc_average, flfft, flflong


def readData(file, file2 = None, SOL_Range_Lc = 0.075, correct_psi = True):
	with open(file) as f:	# this is much faster than np.loadtxt(file), but does the same
		data = np.array([line.strip().split() for line in f if '#' not in line], np.float64)
		
	if 'lam' in file:
		idx = data[:,2] < data[:,2].max()
		Lc = data[:,3][idx]
		psimin = data[:,4][idx]
	else:
		Lc = data[:,3]
		psimin = data[:,4]

	if file2 is not None:
		with open(file2) as f:	# this is much faster than np.loadtxt(file), but does the same
			data2 = np.array([line.strip().split() for line in f if '#' not in line], np.float64)

		if 'lam' in file2:
			idx = data2[:,2] < data2[:,2].max()
			Lc = np.append(Lc,data2[:,3][idx])
			psimin = np.append(psimin,data2[:,4][idx])
		else:
			Lc = np.append(Lc,data2[:,3])
			psimin = np.append(psimin,data2[:,4])

	# Correct psi: use Lc to set psi = 0 in SOL
	if correct_psi: psimin[Lc <= SOL_Range_Lc] = 1
	
	return Lc, psimin
	

def psiHist(psimin, M, SOL_Range_psi, useLogScale = False, save = False, tag = None, show = True, flfLim = None):
	zmin = psimin.min() - 0.01
	zmax = SOL_Range_psi
	dz = (zmax - zmin)/M
	N = len(psimin)

	P = np.zeros(M)
	k = np.ceil((psimin - zmin)/dz)
	count = 0
	for i in xrange(N):
		if(k[i] >= M): continue     # Values larger than zmax are ignored
		P[int(k[i])] += 1
		count += 1

	# Normalize
	P *= 100.0/count

	# psi Axis
	Px = zmin + np.arange(M)*dz + dz/2.0
	
	# Statistical data
	psi_min = np.round(psimin.min()*1000)/1000.0
	flf1a = 100 - np.round(10*np.sum(P[Px <= 0.98]))/10.0
	flf2a = 100 - np.round(10*np.sum(P[Px <= 0.95]))/10.0
	if flfLim is not None:
		flf3a = 100 - np.round(1000*np.sum(P[Px <= flfLim]))/1000.0
		print 'Field line fraction outside of', flfLim, 'is:', flf3a
	
	# make plot
	if show: barPlot(Px, P, useLogScale, 'psi', psi_min, flf1a, flf2a, save = save, tag = tag)
	
	return Px, P, psi_min, flf1a, flf2a


def LcHist(Lc, M, SOL_Range_Lc, Lc_ft, Lc_long, Lc_cutoff, useLogScale = False, save = False, tag = None, show = True):
	zmin = SOL_Range_Lc
	zmax = Lc.max()
	dz = (zmax - zmin)/M
	N = len(Lc)

	L = np.zeros(M)
	k = np.ceil((Lc - zmin)/dz)
	count = 0
	for i in xrange(N):
		if(k[i] <= 0): continue      # All values smaller than zmin are ignored
		if(k[i] >= M): k[i] = M-1    # All values greater than zmax are set to zmax
		L[int(k[i])] += 1
		count += 1

	# Normalize
	L *= 100.0/count

	# Lc Axis
	Lx = zmin + np.arange(M)*dz + dz/2.0
	
	# Statistical data
	Lc_average = np.round(np.sum(Lc[Lc > zmin])/len(Lc[Lc > zmin])*100)/100.0
	flf1b = np.round(10*np.sum(L[Lx <= Lc_ft]))/10.0
	flf2b = np.round(10*np.sum(L[Lx <= Lc_long]))/10.0
	
	# Cutoff
	L = L[Lx <= Lc_cutoff]
	Lx = Lx[Lx <= Lc_cutoff]
	
	# make plot
	if show: barPlot(Lx, L, useLogScale, 'Lc', Lc_average, flf1b, flf2b, Lc_ft, Lc_long, save = save, tag = tag)

	return Lx, L, Lc_average, flf1b, flf2b


def barPlot(x,y,useLogScale,type,a,b,c,d = None,e = None, save = False, tag = None):
	import matplotlib.pyplot as plt
	
	fig = plt.figure()

	dx = x[2] - x[1]
	if useLogScale: y = y+1;
	plt.bar(x,y, width = 0.8*dx);

	if type == 'Lc': plt.xlim(x.min() - dx, x.max() + dx)
	else: plt.xlim(x.min() - dx, 1)
	
	if useLogScale:
		plt.yscale('log')
		plt.gca().set_yticks([1,1.2,1.4,1.7,2,3,4,5,6,8,11,21,31,41,51,71])
		plt.gca().set_yticklabels(['0','0.2','0.4','0.7','1','2','3','4','5','7','10','20','30','40','50','70'])
		plt.ylim(1, 1.1*y.max())
		va = 'top'
	else:  
		plt.ylim(0, 1.1*y.max())
		va = 'center'
	plt.grid(True)

	if type == 'Lc': plt.xlabel('L$_c$ [km]')
	else: plt.xlabel('$\\psi_{min}$')
	plt.ylabel('field line fraction [%]')

	# Inlet
	if type == 'Lc':
		bla = '$\\langle L_{c}\\rangle$ = ' + str(a) + ' km\n' + \
		      'flf$_{' + str(int(d*1000)) + 'm}$ = ' + str(b) + '%\n' + \
		      'flf$_{' + str(int(e*1000)) + 'm}$ = ' + str(c) + '%\n'
		xloc = x.max() - 4*dx
		ha = 'right'
	else:
		bla = '$\\psi_{min}$ = ' + str(a) + '\n' + \
		      'flf$_{98}$ = ' + str(b) + '%\n' + \
		      'flf$_{95}$ = ' + str(c) + '%\n'
		xloc = x.min() + 4*dx
		ha = 'left'
		
	plt.text(xloc,0.92*y.max(),bla, backgroundcolor = 'w', verticalalignment = va, horizontalalignment = ha)
	
	if tag is None: tag = ''
	else: tag = '_' + tag
	if save: fig.savefig('flf_' + type + tag + '.eps', dpi = (300), bbox_inches = 'tight')



#!/usr/bin/env python
# Analysis of psimin occurencies from MAFOT footprint
# Input: file, filepath
# Output: p(psimin) = 0		value of psimin (if any) where probability density function is equal to 0.
#	      P95				integral over psi of probability desnity function equal to 0.95
#         Lcmin				connection length that separates SOL from lobes

import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps 
import getopt


def readfile(file):
    with open(file) as f:
        data = np.array([line.strip().split() for line in f if '#' not in line], np.float64)
    return data


def setvariables(data):
	phi	= data[:,0]
	s 	= data[:,1]
	Nt	= data[:,2]
	Lc	= data[:,3]
	psimin 	= data[:,4]
	R	= data[:,5]
	Z 	= data[:,6]
	LcPsimin= data[:,7]

	return phi, s, Nt, Lc, psimin, R, Z, LcPsimin
	
	
def getLcmin(Lc, psimin):
	SOL = np.where(psimin > 1.)
	Lcmin = Lc[SOL]
	Lcmin.sort()
	return Lcmin[-3::].mean()
	

def getmask(Lc, psimin, Lcmin = 0.06):
	mask = np.zeros(len(Lc), dtype = bool)
	pfr = np.where((psimin < 1.) & (Lc < Lcmin))
	Lobes = np.where((psimin < 1.) & (Lc > Lcmin))
	SOL = np.where(psimin > 1.)

	return SOL,Lobes,pfr


def getoccurancies(psiLobes, Nbins):
	minpsimin = min(psiLobes)
	maxpsimin = max(psiLobes)
	dpsi = (maxpsimin - minpsimin)/Nbins				# 0.0004
	psi = np.linspace(np.floor(minpsimin*1./dpsi)/(1./dpsi),1.,Nbins)

	return psi


def main(file, Lcmin = None, fontsize = 18, legend_loc = (0.08, 0.67), Nbins = 800, span = None, excludeLcmax = True):	
	if span is None: 
		if Nbins > 150: span = 2*int(Nbins/150) + 1
		else: span = 3

	data = readfile(file)
	phi, s, Nt, Lc, psimin, R, Z, LcPsimin = setvariables(data) 
	
	# remove field lines that did not finish trace
	if excludeLcmax:
		Ntmax = Nt.max()
		use = np.where(Nt < Ntmax)[0]
		phi, s, Nt, Lc, psimin, R, Z, LcPsimin = phi[use], s[use], Nt[use], Lc[use], psimin[use], R[use], Z[use], LcPsimin[use]
	
	if Lcmin is None: Lcmin = getLcmin(Lc, psimin)

	SOL, Lobes, pfr = getmask(Lc, psimin, Lcmin)
	n_binspsi = getoccurancies(psimin[Lobes], Nbins)
	print('Number of psi bins = ', np.size(n_binspsi), '  PDF Smoothing span = ', span)
	
	P_95 = np.zeros(np.size(n_binspsi))
	kernel = sps.gaussian_kde(psimin[Lobes], bw_method = 0.002)
	for i in range(np.size(n_binspsi)):
		P_95[i] = kernel.integrate_box_1d(min(n_binspsi),n_binspsi[i])

	lcfs_95 = n_binspsi[max(np.where(P_95 < 0.05)[0])]
	counts, bins = np.histogram(psimin,bins = n_binspsi)
	try: lcfs_0 = n_binspsi[max(np.where(counts == 0)[0]) + 1]
	except: lcfs_0 = False
	
	#Plottig Section
	fig,ax = plt.subplots(figsize=(9, 6), tight_layout = True)
	ax.plot(Lc[SOL],psimin[SOL],'b .', ms = 1)
	ax.plot(Lc[Lobes],psimin[Lobes], 'k .', ms = 1)
	ax.plot(Lc[pfr],psimin[pfr], 'r .', ms = 1)
	plt.plot([0,Lc.max()],[lcfs_95,lcfs_95],'k-', lw = 2)
	if lcfs_0: plt.plot([0,Lc.max()],[lcfs_0,lcfs_0],'g-', lw = 2)
	plt.plot([Lcmin,Lcmin],[0.985,1.01],'b-', lw = 2)

	# Custom Axis title
	ax.set_xlabel("$L_c [km]$", fontsize = fontsize)
	ax.set_ylabel("$\\psi_{min}$", fontsize = fontsize)
	ax.tick_params(axis='both', labelsize = fontsize)
	if legend_loc is not None: ax.legend(("SOL: $\\psi_{min} > 1$ ", "Lobes: $\\psi_{min} < 1, \\; L_c > L_{c,min}$", "Pfr $\\psi_{min} < 1, \\; L_c < L_{c,min}$", "$P_{95}(\\psi_{min})$", "$p(\\psi_{min}) = 0 $", "Lcmin"), loc = 1, handlelength=2, fontsize = fontsize-6)

	fig,ax1 = plt.subplots(figsize=(10,6), tight_layout = True)
	bins = np.size(n_binspsi)
	ax1.hist(psimin[Lobes], bins = n_binspsi, density = True, color = 'skyblue', edgecolor = 'skyblue', linewidth = 0.5, rwidth = 0.5)
	ax1.set_xlabel("$\\psi_{min}$", fontsize = fontsize)
	ax1.set_ylabel("field line fraction [%]", fontsize = fontsize)
	ax1.tick_params(axis = 'both', labelsize = fontsize)
	z = np.convolve(kernel.pdf(n_binspsi), np.ones(span)/span, mode = 'same')
	ax1.plot(n_binspsi, z, color = 'orange', lw = 2)

	ax2 = ax1.twinx()
	ax2.plot(n_binspsi,P_95, color = 'red')
	ax2.set_ylabel("$\\it{P(\\psi_{min})}$", fontsize = fontsize, color = 'red')
	ax2.tick_params(axis = 'y', labelcolor = 'red', labelsize = fontsize)
	ax2.axvline(n_binspsi[max(np.where(P_95 < 0.05)[0])], color = 'black', linewidth = 2, label = 'axvline - full height')
	ax2.text(lcfs_95+0.001, 0.9, format(lcfs_95,'.4f'), color = 'black')

	if lcfs_0:
		ax2.axvline(lcfs_0, color = 'green', linewidth = 2, label = 'axvline - full height')
		ax2.text(lcfs_0+0.001, 0.8, format(lcfs_0,'.4f'), color = 'green')
		print('psi(p = 0) = ', np.round(lcfs_0,4))
	else:
		print('No psimin with p(psi_{min}) = 0')
	print('psi_{min,95} = ', np.round(lcfs_95,4))
	print('Lcmin = ', np.round(Lcmin,4))

	if legend_loc is not None: fig.legend(("$\\it{p(\\psi_{min})}$","$Hist \\; \\psi_{min}$", "$P(\\psi_{min}) = \\int_{\\psi_{min}}^{1} p(\\psi) d \\psi)$", "$P_{95}(\\psi_{min})$", "$p(\\psi_{min}) = 0 $"), loc = legend_loc, fontsize = fontsize-6)
	plt.show()


# ----------------------------------------------------------------------------------------
# --- Launch main() ----------------------------------------------------------------------
if __name__ == '__main__':
	# this is the classic way, used here for backwards compatibility with basic python installs
	Lcmin = None
	fontsize = 18
	loc = (0.12, 0.62)
	Nbins = 800

	opts, args = getopt.gnu_getopt(sys.argv[1:], "hL:f:l:N:", ["help", "Lcmin=", "fontsize=", "loc=", "Nbins="])
	for o, a in opts:
		if o in ("-h", "--help"):
			print ("usage: d3dplot.py [-h] [-L LCMIN] pathname")
			print ("")
			print ("Analyse MAFOT footprint for LCFS")
			print ("")
			print ("positional arguments:")
			print ("  pathname              file name or (full or rel.) pathname")
			print ("")
			print ("optional arguments:")
			print ("  -h, --help            show this help message and exit")
			print ("  -f, --fontsize <Arg>  Font Size")
			print ("  -l, --loc <Arg>       legend location")
			print ("  -L, --Lcmin <Arg>     Guess for Lcmin")
			print ("  -N, --Nbins <Arg>     Number of bins for histogram")
			print ("")
			sys.exit()
		elif o in ("-L", "--Lcmin"):
			Lcmin = float(a)
		elif o in ("-f", "--fontsize"):
			fontsize = int(a)
		elif o in ("-l", "--loc"):
			if(',' in a): 
				loc = a.split(',')
				loc = (float(loc[0]), float(loc[1]))
			elif(a.isdigit()):
				loc = int(a)
			else:
				loc = a
		elif o in ("-N", "--Nbins"):
			Nbins = int(a)
		else:
			raise AssertionError("unknown option")
			
	main(args[0], Lcmin = Lcmin, fontsize = fontsize, legend_loc = loc, Nbins = Nbins, span = None, excludeLcmax = True)


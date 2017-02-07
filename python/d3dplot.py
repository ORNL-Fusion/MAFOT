import os, sys, socket
import getopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap

HOME = os.getenv('HOME')
HOST = socket.gethostname()

def d3dplot(pathname, printme = False, coordinates = 'psi', what = 'psimin', machine = 'd3d', 
			tag = None, graphic = 'png', physical = 1, b = None, N = 60, Title = None,
			typeOfPlot = 'contourf', xlimit = None, ylimit = None, figwidth = None, figheight = None, 
			latex = True, cmap = 'jet'):
	"""
	plot MAFOT results
	--- user input ------------------
	pathname    path/filename
	printme     True: save to File, False: no saving
	coordinates 'RZ', 'psi', 'pest', 'phi'
	what        'Lc', 'psimin', 'psimax', 'psiav'
	machine     'd3d', 'iter'
	tag         arbitary string, attached to the file name of the saved figure
	graphic     'png', 'eps', 'pdf'
	physical    0: t native, 1: t in RZ, 2: t in cm
	b           array, user defined color range (e.g. use linspace), None: defaults are used
	N           int, Number of color levels, used only in default color array
	Title       string, Figure title
	typeOfPlot  'contourf', 'imshow'
	x-,ylimit   tuple (min,max) for the respective axis limits, None: defaults are used
	latex		True: create labels with LaTeX, False: use unicode chars (on systems that do not support LaTeX)
	cmap		string, specifying a matplotlib colormap
	---------------------------------
	"""

	# --- set path from pathname ---
	idx = pathname[::-1].find('/')	# returns location of last '/' in pathname or -1 if not found
	if(idx == -1):
		path = './'
		filename = pathname
	else:
		idx *= -1
		path = pathname[0:idx - 1] + '/'	# path with a final '/'
		filename = pathname[idx::]

	if not (path[-1] == '/'): path += '/'
	if(path[0] == '~'): path = HOME + path[1::]

	# --- auto-set footprint defaults from filename key-word ---
	if(filename[0:4] == 'foot'):
		coordinates = 'phi'
		if(filename[5:7] == 'in'): target = 'in'
		elif(filename[5:7] == 'ou'): target = 'out'
		elif(filename[5:7] == 'sh'): target = 'shelf'
		elif(filename[5:7] == 'sa'): target = 'sas'
		else: 
			print 'cannot identify target'
			return
	else: target = None
	
	# --- define reversed cmap ---
	if cmap[-2::] == '_r': cmap_r = cmap[0:-2]
	else: cmap_r = cmap + '_r'
		
	# --- shot and time from file header ---
	with open(path + filename, 'r') as f:
		# Skip first 3 lines
		for i in xrange(3): head = f.readline()
		# read Shot and Time
		head = f.readline().split()
		shot = int(head[-1])
		head = f.readline().split()
		time = int(head[-1])
		
	# --- read data -------------------
	data = readfile(path + filename)
	Nth, Npsi = gridSize(data)
	if(coordinates == 'RZ') | (coordinates == 'phi'): Nth, Npsi = Npsi, Nth
	x = data[:,0].reshape(Nth,Npsi)
	y = data[:,1].reshape(Nth,Npsi)
	Lc = data[:,3].reshape(Nth,Npsi)
	psimin = data[:,4].reshape(Nth,Npsi)
	try:
		psimax = data[:,5].reshape(Nth,Npsi)
		psiav = data[:,6].reshape(Nth,Npsi)
		use_psimaxav = True
	except:
		use_psimaxav = False
	try:
		pitch = data[:,7].reshape(Nth,Npsi)
		yaw = data[:,8].reshape(Nth,Npsi)
		use_pitch_yaw = True
	except:
		use_pitch_yaw = False

	# --- set xy data -----------------
	if(coordinates == 'RZ'):
		xLabel = 'R [m]'
		yLabel = 'Z [m]'

	elif(coordinates == 'psi'):
		if not latex:
			xLabel = u'\u03B8' + ' [rad]'
			yLabel = u'\u03c8'
		else:
			xLabel = '$\\theta$ $\\mathrm{[rad]}$'
			yLabel = '$\\psi$'
	
	elif(coordinates == 'pest'):
		if not latex:
			xLabel = u'\u03B8' + '$_p$ [rad]'
			yLabel = u'\u03c8'
		else:
			xLabel = '$\\theta_p$ $\\mathrm{[rad]}$'
			yLabel = '$\\psi$'

	elif(coordinates == 'phi'):
		if('d3d' in machine):
			if not latex: xLabel = u'\u03C6' + ' [rad]'
			else: xLabel = '$\\varphi$ $\\mathrm{[rad]}$'
			yLabel = 't'
			if(physical > 0):
				x = (360 - x*180.0/np.pi)[:,::-1]	# reverse x
				if not latex: xLabel = u'\u03C6' + ' [deg]'
				else: xLabel = '$\\phi$ $\\mathrm{[deg]}$'
				yLabel = 'R [m]'
				if(target== 'in'):
					if(physical == 1):
						y[y >= 0] *= 0.14
						y[y < 0] *= 0.1959
						y = -1.223 - y
						y = y[::-1,:]	# reverse y
						yLabel = 'Z [m]'
					elif(physical == 2): 
						y *= 19.59
						yLabel = 't [cm]'
				elif(target == 'out'): y = 1.153 + y*0.219
				elif(target == 'shelf'): y = 1.372 + y*0.219
				elif(target == 'sas'):
					physical = 2
					yLabel = 't [cm]'
					y = y*101.693189
		elif(machine == 'iter'):
			if not latex: xLabel = u'\u03C6' + ' [rad]'
			else: xLabel = '$\\varphi$ $\\mathrm{[rad]}$'
			yLabel = 't [cm]'
			if(physical > 0):
				x = x*180.0/np.pi
				if not latex: xLabel = u'\u03C6' + ' [deg]'
				else: xLabel = '$\\phi$ $\\mathrm{[deg]}$'
								
	else:
		print 'coordinates: Unknown input'
		return

	# --- set z data ------------------
	if(what == 'Lc'):
		if(b == None): 
			if('d3d' in machine): b = np.linspace(0.075,0.4,N)
			elif(machine == 'iter'): b = np.linspace(0.22,1.8,N)
		z = Lc
		if not latex: C_label = '$L_{c}$ [km]'
		else: C_label = '$L_{c}$ $\\mathrm{[km]}$'
		#usecolormap = cm.jet	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.get_cmap(cmap)._segmentdata
		
	elif(what == 'psimin'):
		if(b == None): b = np.linspace(0.88,1.02,N)
		z = psimin
		if not latex: C_label = u'\u03c8' + '$_{Min}$'
		else: C_label = '$\\psi_{Min}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.get_cmap(cmap_r)._segmentdata

	elif(what == 'psimax'):
		if not use_psimaxav:
			print "psimax data is not available in this file"
			return
		if(b == None): b = np.linspace(0.88,1.02,N)
		z = psimax
		if not latex: C_label = u'\u03c8' + '$_{Max}$'
		else: C_label = '$\\psi_{Max}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.get_cmap(cmap_r)._segmentdata
		
	elif(what == 'psiav'):
		if not use_psimaxav:
			print "psiav data is not available in this file"
			return
		if(b == None): b = np.linspace(0.88,1.02,N)
		z = psiav
		if not latex: C_label = u'\u03c8' + '$_{av}$'
		else: C_label = '$\\psi_{av}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.get_cmap(cmap_r)._segmentdata

	elif(what == 'pitch'):
		if not use_pitch_yaw:
			print "pitch angle data is not available in this file"
			return
		if(b == None): b = np.linspace(-5,20,N)
		z = pitch/np.pi*180
		if not latex: C_label = u'\u03b1' + '$_{p}$ [deg]'
		else: C_label = '$\\alpha_{p}$ [deg]'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.get_cmap(cmap)._segmentdata

	elif(what == 'yaw'):
		if not use_pitch_yaw:
			print "yaw angle data is not available in this file"
			return
		if(b == None): b = np.linspace(-10,10,N)
		z = yaw/np.pi*180
		if not latex: C_label = u'\u03b1' + '$_{r}$  [deg]'
		else: C_label = '$\\alpha_{r}$  [deg]'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.get_cmap(cmap)._segmentdata
		
	else: 
		print 'what: Unknown input'
		return
		
	usecolormap = LinearSegmentedColormap('my_cmap', cdict, len(b))
	
	# correct for PFR
	if('d3d' in machine): Lcmin = 0.075
	elif(machine == 'iter'): Lcmin = 0.22
	
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		if('psi' in what):
			z[Lc < Lcmin] = 1.01*b.max()
	elif not (coordinates == 'phi'): 
		if(what == 'Lc'): z[(y >= 1) & (z >= 4)] = b.min()
		elif('psi' in what): z[z < 0.5] = b.max()
	
	# reverse axes for footprint
	if(coordinates == 'phi') & (physical > 0):
		z = z[:,::-1]		# reverse x-axis
		if(target== 'in') & (physical == 1):
			z = z[::-1,:]	# reverse y-axis
		
	# correct orientation in RZ plot and footprint
	if(coordinates == 'RZ') | (coordinates == 'phi'): x, y, z = x.T, y.T, z.T
	
	# correct, if last x point is larger than 2pi -> x mod 2pi is very small
	if(x[-1,0] < x[-2,0]): x = x[0:-1,:]; y = y[0:-1,:]; z = z[0:-1,:]

	# --- Pest theta ------------------
	if(coordinates == 'pest'):
		try:
			import scipy.interpolate as inter
			import Misc.sfc_class as sfc_class
			sfc = sfc_class.straight_field_line_coordinates(shot, time)
			xp = sfc.ev(x,y)
			for i in xrange(Npsi):
				spline = inter.UnivariateSpline(xp[:,i], z[:,i], s = 0)
				z[:,i] = spline(x[:,i])
		except:
			raise ImportError('Pest coordinates not available. Choose other coordinates')

	# --- layout ----------------------
	#rcParams['text.latex.preamble'] = [r'\usepackage{times}']#\usepackage{amsmath}
	if 'Anaconda' in sys.version:
		plt.rcParams['font.sans-serif'] = 'Arial'	# anaconda
		plt.rcParams['font.serif'] = 'Times'	# anaconda
	else:
		plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Bitstream Vera Sans']	# enthought
		plt.rcParams['font.serif'] = ['Times', 'Times New Roman', 'Bitstream Vera Serif']	# enthought
	font = {'family' : 'sans-serif',
			'weight' : 'normal', # normal
			'size'   : 22} #18
	plt.matplotlib.rc('font', **font)
	if not latex: latexFontSize = font['size']
	else: latexFontSize = 24

	if(coordinates == 'RZ'): 
		width = 9
		height = 7
		xlabel_size = font['size']
		ylabel_size = font['size']
	else: 
		width = 10
		height = 6
		xlabel_size = latexFontSize
		ylabel_size = latexFontSize
	if(coordinates == 'phi'):
		ylabel_size = font['size']
		
	if figwidth == None: figwidth = width
	if figheight == None: figheight = height
	
	C_label_size = latexFontSize

	# --- make plot -------------------
	F = plt.figure(figsize = (figwidth,figheight))
	if(typeOfPlot == 'contourf'):
		plt.contour(x, y, z, b, cmap = usecolormap, extend = 'both')
		cs = plt.contourf(x, y, z, b, cmap = usecolormap, extend = 'both')
	else:
		cs = plt.imshow(z.T, extent = [x.min(), x.max(), y.min(), y.max()], cmap = usecolormap, origin = 'lower', vmin = b.min(), vmax = b.max(), aspect = 'auto', interpolation = 'bilinear')

	# --- additional plots ------------
	if(coordinates == 'RZ'): # plot wall
		wall = get_wall(machine)
		plt.plot(wall[:,0], wall[:,1], 'k--', linewidth = 2)
		
	if(coordinates == 'phi'):
		if (target == 'in'):	# plot tile limit
			if(physical == 1) & ('d3d' in machine): plt.plot([x.min(), x.max()], [-1.223, -1.223], 'k--', linewidth = 1.5)
			else: plt.plot([x.min(), x.max()], [0, 0], 'k--', linewidth = 1.5)
		elif (target == 'sas') & ('d3d' in machine):	# plot sas edges
			if(physical == 0):
				plt.plot([x.min(), x.max()], [0.17276, 0.17276], 'k--', linewidth = 1.5)
				plt.plot([x.min(), x.max()], [0.41, 0.41], 'k--', linewidth = 1.5)			
				plt.plot([x.min(), x.max()], [0.28363, 0.28363], 'k--', linewidth = 1.5)
				plt.plot([x.min(), x.max()], [0.296, 0.296], 'k--', linewidth = 1.5)			
			else:
				plt.plot([x.min(), x.max()], [17.569, 17.569], 'k--', linewidth = 1.5)
				plt.plot([x.min(), x.max()], [41.736, 41.736], 'k--', linewidth = 1.5)			
				plt.plot([x.min(), x.max()], [28.843, 28.843], 'k--', linewidth = 1.5)
				plt.plot([x.min(), x.max()], [30.11, 30.11], 'k--', linewidth = 1.5)			
		else:	# plot the pump limit / edge of shelf nose location
			if(physical == 1) & ('d3d' in machine): plt.plot([x.min(), x.max()], [1.372, 1.372], 'k--', linewidth = 1.5)
			else: plt.plot([x.min(), x.max()], [1, 1], 'k--', linewidth = 1.5)
			
	if(coordinates == 'RZ') & (what in ['pitch', 'yaw']): # plot plasma boundary
		plt.contour(x, y, Lc.T, [Lcmin], colors = 'k', linewidths = 2)

	# --- Axes  -----------------------
	if(coordinates == 'RZ'): 
		plt.axes().set_aspect('equal')
		plt.xlim(x.min(),x.max())
		plt.ylim(y.min(),y.max())
		try: plt.locator_params(axis = 'x', nbins = 5)
		except: pass
	elif(coordinates == 'phi') | (coordinates == 'psi'): 
		plt.xlim(x.min(),x.max())
		plt.ylim(y.min(),y.max())
	elif(coordinates == 'pest'): 
		plt.ylim(y.min(), min([y.max(), 1.0]))
	if(coordinates == 'phi') & (target== 'in') & (physical == 2): 
		plt.gca().invert_yaxis()
		
	if not (xlimit == None): plt.xlim(xlimit)
	if not (ylimit == None): plt.ylim(ylimit)

	plt.xlabel(xLabel, size = xlabel_size)
	plt.ylabel(yLabel, size = ylabel_size, labelpad = 10)
	
	# adjust axis position on screen; leaves savefile unaffected
	pos = plt.gca().get_position().get_points()	# get bounding-box
	pos[0,0] += 0.025; pos[1,0] += 0.025	# x-position: left; right
	pos[0,1] += 0.05; pos[1,1] += 0.05		# y-position: bottom; top
	plt.gca().set_position(plt.matplotlib.transforms.Bbox(pos))
	
	# --- Colorbar --------------------
	# set ticks
	steps = (b.max() - b.min())/10.0
	digit = int(np.floor(np.log10(steps)))
	factor = 10**digit
	steps = int(np.ceil(steps/factor))*factor
	bmin = int(np.ceil(b.min()/factor))*factor
	myticks = np.arange(bmin, b.max()+steps, steps)
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		if('psi' in what): 
			myticks = myticks[myticks <= b.max()]
			myticks[-1] = b.max()
		else:
			myticks[0] = b.min()
			
	if 'Anaconda' in sys.version:
		if('psi' in what): 
			myticks = myticks[myticks <= b.max()]
	
	myticks[np.abs(myticks) < 1e-10] = 0
	
	# If the "extend" argument is given, contourf sets the data limits to some odd extension of the actual data.
	# Resetting the data limits, after plotting the contours, forces set_over & set_under colors to show,
	# once the colorbar is called after the reset.
	if(typeOfPlot == 'contourf'): cs.set_clim(b.min(),b.max())
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		if('psi' in what): cs.cmap.set_over('w')
		elif(what == 'Lc'): cs.cmap.set_under('w')

	# show colorbar
	C = plt.colorbar(cs, pad = 0.01, extend = 'both', format = '%.3g', ticks = myticks)
	if 'Anaconda' in sys.version: C.set_label(C_label, rotation = 270, size = C_label_size, va = 'bottom')	# anaconda
	elif (HOST == 'head.cluster'): C.set_label(C_label, rotation = 270, size = C_label_size, va = 'bottom') # Drop Cluster
	else: C.set_label(C_label, rotation = 270, size = C_label_size)	# enthought
	
	# add SOL label in RZ plot and footprint
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		myticklabels = [str(item) for item in myticks]
		if('psi' in what): myticklabels[-1] = 'SOL'
		elif(what == 'Lc'): myticklabels[0] = 'SOL'
		C.ax.set_yticklabels(myticklabels)
						
	# --- Title -----------------------
	if(Title == True): plt.title(str(shot) + ',  ' + str(time) + 'ms', size = 18)
	elif not (Title == None): plt.title(Title, size = 18)

	# --- save Figure to file ---------
	if not (tag == None): printme = True
	if printme: 
		if(tag == None): F.savefig(path + filename[0:-4] + '.' + graphic, dpi = (300), bbox_inches = 'tight')
		else: F.savefig(path + filename[0:-4] + '_' + tag + '.' + graphic, dpi = (300), bbox_inches = 'tight')


# --- data = readfile(name) ---
# reads in ascii data file
def readfile(name):
	with open(name, 'r') as f:
		# Skip header, comment char is '#'
		head = f.readline().split()
		while head[0][0] == '#':
			head = f.readline().split()
		
		# last head line is already first data line
		data = [[float(x) for x in head]]	# data is a 2D array <- double brackets
		
		# read data to end of file
		for line in f:
			line = line.split()			
			data.append([float(x) for x in line])
			
	# cast data into numpy array and return
	data = np.array(data)
	return data
	
	
# --- Nx, Ny = get_gridSize(data) ---
# determines 2D grid sizes
# data = output from readfile
# assumes that the first two columns of data are the flattend mashgrid of x and y respectively
def gridSize(data):
	x = data[:,0]
	Ny = len(x[x == x[0]])	# number of dublicates == number of grid points on other axis
	y = data[:,1]
	Nx = len(y[y == y[0]])
	
	if not (len(x) == Nx*Ny):
		raise RuntimeError('Grid size could not be determined')
		
	return Nx, Ny


def get_wall(machine):
	if machine == 'd3d':
		wall = [[1.41903,	1.34783],[1.28008,	1.34773],[1.27999,	1.33104],
				[1.24803,	1.2539],[1.22784,	1.22143],[1.20913,	1.20165],
				[1.19011,	1.18772],[1.16185,	1.1759],[1.11593,	1.16606],
				[1.04209,	1.16223],[1.02923,	1.21657],[1.00088,	1.21721],
				[1.01777,	1.13839],[1.01779,	1.00132],[1.01619,	1.00132],
				[1.01627,	-1.21079],[1.01608,	-1.22884],[1.15285,	-1.3664],
				[1.19734,	-1.3664],[1.24185,	-1.3664],[1.28635,	-1.3664],
				[1.33085,	-1.3664],[1.37535,	-1.3664],[1.41985,	-1.3664],
				[1.41985,	-1.329],[1.372,	-1.329],[1.372,-1.25],
				[1.372,	-1.25],[1.57,	-1.25],[1.76801,	-1.25],
				[1.76801,	-1.21104],[1.78603,	-1.17428],[2.13401,	-0.973024],
				[2.37704,	-0.389023],[2.35342,	-0.400247],[2.35116,	-0.337078],
				[2.354,	0],[2.35082,	0.204946],[2.3523,	0.400178],
				[2.37704,	0.389023],[2.1282,	0.993424],[2.0699,	1.03975],
				[1.78499,	1.07688],[1.647,	1.07675],[1.60799,	1.09525],
				[1.37198,	1.29208],[1.3721,	1.30954],[1.41897,	1.31017]]	
	elif machine == 'd3d_sas':
		wall = [[ 1.016  ,  0.     ],[ 1.016  ,  0.964  ],[ 1.016  ,  0.968  ],[ 1.016  ,  1.001  ],
				[ 1.016  ,  1.019  ],[ 1.016  ,  1.077  ],[ 1.016  ,  1.07   ],[ 1.016  ,  1.096  ],
				[ 1.016  ,  1.113  ],[ 1.016  ,  1.138  ],[ 1.016  ,  1.147  ],[ 1.012  ,  1.165  ],
				[ 1.001  ,  1.217  ],[ 1.029  ,  1.217  ],[ 1.042  ,  1.1624 ],[ 1.046  ,  1.16238],
				[ 1.056  ,  1.1626 ],[ 1.097  ,  1.1645 ],[ 1.108  ,  1.16594],[ 1.116  ,  1.16591],
				[ 1.134  ,  1.16896],[ 1.148  ,  1.17175],[ 1.162  ,  1.17556],[ 1.181  ,  1.183  ],
				[ 1.182  ,  1.1835 ],[ 1.185  ,  1.185  ],[ 1.19   ,  1.188  ],[ 1.195  ,  1.191  ],
				[ 1.201  ,  1.196  ],[ 1.209  ,  1.202  ],[ 1.215  ,  1.208  ],[ 1.222  ,  1.214  ],
				[ 1.228  ,  1.221  ],[ 1.234  ,  1.231  ],[ 1.239  ,  1.238  ],[ 1.242  ,  1.244  ],
				[ 1.248  ,  1.254  ],[ 1.258  ,  1.278  ],[ 1.263  ,  1.29   ],[ 1.28   ,  1.331  ],
				[ 1.28   ,  1.347  ],[ 1.28   ,  1.348  ],[ 1.31   ,  1.348  ],[ 1.328  ,  1.348  ],
				[ 1.361  ,  1.348  ],[ 1.38   ,  1.348  ],[ 1.419  ,  1.348  ],[ 1.419  ,  1.31   ],
				[ 1.372  ,  1.31   ],[ 1.37167,  1.29238],[ 1.37003,  1.28268],[ 1.36688,  1.25644],
				[ 1.36719,  1.22955],[ 1.37178,  1.19576],[ 1.37224,  1.19402],[ 1.38662,  1.16487],
				[ 1.38708,  1.16421],[ 1.40382,  1.15696],[ 1.41127,  1.1573 ],[ 1.41857,  1.16132],
				[ 1.421  ,  1.164  ],[ 1.48663,  1.2405 ],[ 1.4973 ,  1.23458],[ 1.49762,  1.23428],
				[ 1.49745,  1.23174],[ 1.49275,  1.2133 ],[ 1.4926 ,  1.21061],[ 1.49261,  1.20486],
				[ 1.49279,  1.20214],[ 1.4934 ,  1.19642],[ 1.4947 ,  1.18511],[ 1.49622,  1.1607 ],
				[ 1.47981,  1.12426],[ 1.48082,  1.12256],[ 1.48149,  1.12138],[ 1.48646,  1.11692],
				[ 1.49095,  1.11439],[ 1.50305,  1.11244],[ 1.59697,  1.09489],[ 1.6255 ,  1.0853 ],
				[ 1.63752,  1.07988],[ 1.647  ,  1.077  ],[ 1.785  ,  1.077  ],[ 2.07   ,  1.04   ],
				[ 2.128  ,  0.993  ],[ 2.245  ,  0.709  ],[ 2.33956,  0.46143],[ 2.34708,  0.41583],
				[ 2.34913,  0.27218],[ 2.35103,  0.17018],[ 2.35158,  0.07012],[ 2.35125, -0.03179],
				[ 2.35051, -0.14435],[ 2.34965, -0.21483],[ 2.3487 , -0.32669],[ 2.3476 , -0.38677],
				[ 2.3402 , -0.45304],[ 2.32942, -0.47757],[ 2.134  , -0.973  ],[ 1.786  , -1.174  ],
				[ 1.768  , -1.211  ],[ 1.768  , -1.25   ],[ 1.682  , -1.25   ],[ 1.372  , -1.25   ],
				[ 1.372  , -1.329  ],[ 1.42   , -1.329  ],[ 1.42   , -1.363  ],[ 1.273  , -1.363  ],
				[ 1.153  , -1.363  ],[ 1.016  , -1.223  ],[ 1.016  , -1.223  ],[ 1.016  , -0.83   ],
				[ 1.016  , -0.8    ],[ 1.016  , -0.415  ],[ 1.016  , -0.4    ],[  1.016,  -0.001]]
	elif machine == 'iter':
		wall = [[4.05460000,   -2.50630000],[4.05460000,   -1.50000000],[4.05460000,   -0.48360000],
				[4.05460000,   0.53280000],[4.05460000,   1.54920000],[4.05460000,   2.56560000],
				[4.05460000,   3.58200000],[4.32000000,   4.32400000],[4.91280000,   4.71150000],
				[5.76290000,   4.53230000],[6.59610000,   3.89340000],[7.47630000,   3.08330000],
				[7.94290000,   2.40240000],[8.27940000,   1.68140000],[8.40350000,   0.63290000],
				[8.31540000,   -0.42150000],[7.90780000,   -1.34170000],[7.29200000,   -2.25700000],
				[6.27560000,   -3.04610000],[6.16250000,   -3.23970000],[5.97380000,   -3.28690000],
				[5.80690000,   -3.38700000],[5.67620000,   -3.53110000],[5.59300000,   -3.70700000],
				[5.56640000,   -3.89850000],[5.56440000,   -3.99970000],[5.55740000,   -4.09980000],
				[5.55740000,   -4.42530000],[5.55740000,   -4.42810000],[5.55736000,   -4.55894000],
				[5.54940000,   -4.55070000],[5.45480000,   -4.45610000],[5.36020000,   -4.38150000],
				[5.26550000,   -4.28680000],[5.24580000,   -3.98890000],[5.14250000,   -3.84210000],
				[4.99140000,   -3.74540000],[4.81490000,   -3.71300000],[4.63930000,   -3.75000000],
				[4.48570000,   -3.91300000],[4.38470000,   -3.90500000],[4.28380000,   -3.89710000],
				[4.18280000,   -3.88910000],[4.17395000,   -3.88824000],[4.22630000,   -3.77750000],
				[4.29550000,   -3.63990000],[4.37140000,   -3.48420000],[4.40040000,   -3.40880000],
				[4.44610000,   -3.28470000],[4.50950000,   -3.11880000],[4.50020000,   -2.94610000],
				[4.43480000,   -2.78610000],[4.31980000,   -2.65690000],[4.16650000,   -2.57300000],
				[4.05460000,   -2.50630000]]
	else:
		raise AssertionError('unkown machine')
	return np.array(wall)
	


# ----------------------------------------------------------------------------------------
# --- Launch main() ----------------------------------------------------------------------
if __name__ == '__main__':
	# this is the more modern way
	"""
	import argparse
	import textwrap
	parser = argparse.ArgumentParser(description = 'Plot MAFOT output', 
				formatter_class = argparse.RawDescriptionHelpFormatter,
				epilog = textwrap.dedent('''\
                Examples: d3dplot.py foot_in_test.dat
                          d3dplot.py /path/to/gfile/foot_in_test.dat -p -c RZ
                          d3dplot.py lam_psi_nopr.dat -w Lc'''))

	parser.add_argument('pathname', help = 'file name or (full or rel.) pathname', type = str)
	parser.add_argument('-c','--coordinates', help = 'Coordinate sytem for plot, options are: RZ, psi (default), phi, pest', type = str, default = 'psi')
	parser.add_argument('-w','--what', help = 'Data to plot, options are: Lc, psimin (default), psimax, psiav', type = str, default = 'psimin')
	parser.add_argument('-m','--machine', help = 'Machine, options are: d3d (default), iter', type = str, default = 'd3d')
	parser.add_argument('-p','--printme', help = 'Save to File', action = 'store_true', default = False)
	parser.add_argument('-g','--graphic', help = 'Graphic extension, options are: png (default), eps, pdf', type = str, default = 'png')
	parser.add_argument('-t','--tag', help = 'Arbitary string, attached to the file name of the saved figure', type = str, default = None)
	parser.add_argument('-P','--physical', help = 'Type of y-Axis in fooprint. 0: native, 1: RZ (default), 2: t in cm', type = int, default = 1)
	parser.add_argument('-N', help = 'Number of color levels in plot, default = 60', type = int, default = 60)
	parser.add_argument('-T','--Title', help = 'Figure title, default = None', type = str, default = None)	
	parser.add_argument('-i','--imshow', help = 'Use imshow instead of contourf (default)', action = 'store_true', default = False)
	parser.add_argument('-W','--figwidth', help = 'Force width of figure from default to value', type = int, default = None)
	parser.add_argument('-H','--figheight', help = 'Force height of figure from default to value', type = int, default = None)
	
	args = parser.parse_args()
	if args.imshow: toP = 'imshow'
	else: toP = 'contourf'
	d3dplot(args.pathname, printme = args.printme, coordinates = args.coordinates, what = args.what, machine = args.machine, 
			tag = args.tag, graphic = args.graphic, physical = args.physical, b = None, N = args.N, Title = args.Title,
			typeOfPlot = toP, xlimit = None, ylimit = None, figwidth = args.figwidth, figheight = args.figheight)
	"""	

	# this is the classic way, used here for backwards compatibility with basic python installs
	coordinates = 'psi'
	what = 'psimin'
	machine = 'd3d'
	printme = False
	graphic = 'png'
	tag = None
	physical = 1
	N = 60
	cmap = 'jet'
	Title = None
	toP = 'contourf'
	figwidth = None
	figheight = None
	#if(HOST == 'head.cluster'): latex = False	# Drop does not support LaTeX
	#else: latex = True
	latex = True		# Drop now supports LaTeX too
	b = None;		set_color = False
	xlim = None
	ylim = None

	opts, args = getopt.gnu_getopt(sys.argv[1:], "hc:w:m:pg:t:P:N:C:T:iW:H:Ub:x:y:", ["help", "coord=", "what=", 
																 "machine=", "printme", "graphic=", 
																 "tag=", "physical=", "cmap=", "Title=", "imshow", 	 
																 "figwidth=", "figheight=", "unicode",
																 "range", "xlim", "ylim"])
	for o, a in opts:
		if o in ("-h", "--help"):
			print "usage: d3dplot.py [-h] [-c COORDINATES] [-w WHAT] [-m MACHINE] [-p]"
			print "                  [-g GRAPHIC] [-t TAG] [-P PHYSICAL] [-N N] [-C CMAP]"
			print "                  [-T TITLE] [-i] [-W FIGWIDTH] [-H FIGHEIGHT] [-U]"
			print "                  [-b MIN,MAX] [-x MIN,MAX] [-y MIN,MAX]"
			print "                  pathname"
			print ""
			print "Plot MAFOT output"
			print ""
			print "positional arguments:"
			print "  pathname              file name or (full or rel.) pathname"
			print ""
			print "optional arguments:"
			print "  -h, --help            show this help message and exit"
			print "  -c, --coord <Arg>     Coordinate sytem for plot. <Arg> = "
			print "                        RZ, psi(default), phi, pest"
			print "  -w, --what <Arg>      Data to plot. <Arg> = "
			print "                        Lc, psimin (default), psimax, psiav"
			print "  -m, --machine <Arg>   Machine. <Arg> = d3d (default), iter"
			print "  -p, --printme         Save to File"
			print "  -g, --graphic <Arg>   Graphic extension. <Arg> = "
			print "                        png (default), eps, pdf"
			print "  -t, --tag <Arg>       Arbitary string, attached to the file name of"
			print "                        the saved figure"
			print "  -P, --physical <Arg>  Type of y-Axis in fooprint. <Arg> = "
			print "                        0: native, 1: RZ(default), 2: t in cm"
			print "  -N <Arg>              Number of color levels in plot, default = 60"
			print "  -C, --cmap <Arg>      colormap, default = jet"
			print "  -T, --Title <Arg>     Figure title, default = None"
			print "  -i, --imshow          Use imshow instead of contourf (default)"
			print "  -W, --figwidth <Arg>  Force width of figure from default to <Arg>"
			print "  -H, --figheight <Arg> Force height of figure from default to <Arg>"
			print "  -U, --unicode         Use Unicode instead of LaTeX (default) in labels"
			print "                        Note: Drop Cluster uses Unicode by default"
			print "  -b, --range <Arg>     ColorBar range: Min,Max (no spaces)"
			print "  -x, --xlim <Arg>      X-Axis range: Min,Max (no spaces)"
			print "  -y, --ylim <Arg>      Y-Axis range: Min,Max (no spaces)"
			print ""
			print "Examples: d3dplot.py foot_in_test.dat"
			print "          d3dplot.py /path/to/gfile/foot_in_test.dat -p -c RZ"
			print "          d3dplot.py lam_psi_nopr.dat -w Lc"
			sys.exit()
		elif o in ("-c", "--coordinates"):
			coordinates = a
		elif o in ("-w", "--what"):
			what = a
		elif o in ("-m", "--machine"):
			machine = a
		elif o in ("-p", "--printme"):
			printme = True
		elif o in ("-g", "--graphic"):
			graphic = a
		elif o in ("-t", "--tag"):
			tag = a
		elif o in ("-P", "--physical"):
			physical = int(a)
		elif o in ("-N",):
			N = int(a)
		elif o in ("-C", "--cmap"):
			cmap = a
		elif o in ("-T", "--Title"):
			Title = a
		elif o in ("-i", "--imshow"):
			toP = 'imshow'
		elif o in ("-W", "--figwidth"):
			figwidth = float(a)
		elif o in ("-H", "--figheight"):
			figheight = float(a)
		elif o in ("-U", "--unicode"):
			latex = False
		elif o in ("-b", "--range"):
			range = a.split(',')
			cmin, cmax = float(range[0]), float(range[1])
			set_color = True
		elif o in ("-x", "--xlim"):
			range = a.split(',')
			xlim = (float(range[0]), float(range[1]))
		elif o in ("-y", "--ylim"):
			range = a.split(',')
			ylim = (float(range[0]), float(range[1]))
		else:
			raise AssertionError("unknown option")
			
	if set_color: b = np.linspace(cmin,cmax,N)
	
	d3dplot(args[0], printme = printme, coordinates = coordinates, what = what, machine = machine, 
			tag = tag, graphic = graphic, physical = physical, b = b, N = N, Title = Title,
			typeOfPlot = toP, xlimit = xlim, ylimit = ylim, figwidth = figwidth, figheight = figheight,
			latex = latex, cmap = cmap)

	plt.show()


import os
import pickle
from numpy import *
from matplotlib.pyplot import *
import scipy.interpolate as inter
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap

import readfile as rf
from load_gfile_d3d import *
import myColorMaps; reload(myColorMaps)

HOME = os.getenv('HOME')


def d3dplot(filename, path = './', printme = False, tag = None, graphic = 'png', typeOfPlot = 'contourf', 
			coordinates = 'psi', what = 'psimin', target = None, physical = 1, b = None, N = 256, Title = None,
			myticks = None, showSOL = True):
	# --- user input ------------------
	# printme		True: save to File, False: no saving
	# tag			arbitary string, attached to the file name of the saved figure
	# graphic		'png', 'eps', 'pdf'
	# typeOfPlot	'contourf', 'imshow'
	# coordinates 	'RZ', 'psi', 'pest', 'phi'
	# what 			'Lc', 'psimin', 'psimax', 'psiav'
	# target		'in', 'out', 'shelf'
	# physical		0: t native, 1: t in RZ, 2: t in cm
	# b				array, user defined color range (e.g. use linspace), None: defaults are used
	# N				int, Number of color levels, used only in default color array
	# Title			string, Figure title
	# mytick		array of where ticks are drawn on colorbar
	# useSOL		bool, replace SOL tick with word 'SOL' on colorbar

	if not (path[-1] == '/'): path += '/'
	if(path[0] == '~'): path = HOME + path[1::]
	
	if(filename[0:4] == 'foot') & (target == None):
		coordinates = 'phi'
		if(filename[5:7] == 'in'): target = 'in'
		elif(filename[5:7] == 'ou'): target = 'out'
		elif(filename[5:7] == 'sh'): target = 'shelf'
		else: 
			print 'Please specify target'
			return
		
	# --- shot and time from fiel header ---
	with open(path + filename, 'r') as f:
		# Skip first 3 lines
		for i in xrange(3): head = f.readline()
		# read Shot and Time
		head = f.readline().split()
		shot = int(head[-1])
		head = f.readline().split()
		time = int(head[-1])
		
	# --- read data -------------------
	data = rf.readfile(path + filename)
	Nth, Npsi = rf.gridSize(data)
	if(coordinates == 'RZ') | (coordinates == 'phi'): Nth, Npsi = Npsi, Nth
	x = data[:,0].reshape(Nth,Npsi)
	y = data[:,1].reshape(Nth,Npsi)
	Lc = data[:,3].reshape(Nth,Npsi)
	psimin = data[:,4].reshape(Nth,Npsi)
	if (coordinates == 'psi') | (coordinates == 'pest'):
		psimax = data[:,5].reshape(Nth,Npsi)
		psiav = data[:,6].reshape(Nth,Npsi)

	# --- set xy data -----------------
	if(coordinates == 'RZ'):
		xLabel = 'R [m]'
		yLabel = 'Z [m]'

	elif(coordinates == 'psi'):
		xLabel = '$\\theta$ $\\mathrm{[rad]}$'
		yLabel = '$\\psi$'
	
	elif(coordinates == 'pest'):
		xLabel = '$\\theta_p$ $\\mathrm{[rad]}$'
		yLabel = '$\\psi$'

	elif(coordinates == 'phi'):
		xLabel = '$\\varphi$ $\\mathrm{[rad]}$'
		yLabel = 't'
		if(physical > 0):
			x = (360 - x*180.0/pi)[:,::-1]	# reverse x
			xLabel = '$\\phi$ $\\mathrm{[deg]}$'
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
	else:
		print 'coordinates: Unknown input'
		return

	# --- set z data ------------------
	if(what == 'Lc'):
		if(b == None): b = linspace(0.075,0.4,N)
		z = Lc
		C_label = '$L_{c}$ $\\mathrm{[km]}$'
		#usecolormap = cm.jet	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = cm.jet._segmentdata
		
	elif(what == 'psimin'):
		if(b == None): b = linspace(0.88,1.02,N)
		z = psimin
		C_label = '$\\psi_{Min}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = cm.jet_r._segmentdata

	elif(what == 'psimax'):
		if(b == None): b = linspace(0.88,1.02,N)
		z = psimax
		C_label = '$\\psi_{Max}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = cm.jet_r._segmentdata
		
	elif(what == 'psiav'):
		if(b == None): b = linspace(0.88,1.02,N)
		z = psiav
		C_label = '$\\psi_{av}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = cm.jet_r._segmentdata
		
	else: 
		print 'what: Unknown input'
		return
		
	usecolormap = LinearSegmentedColormap('my_cmap', cdict, len(b))
	
	# correct for PFR
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		if not (what == 'Lc'):
			z[Lc < 0.075] = 1.01*b.max()
	elif not (coordinates == 'phi'): 
		if(what == 'Lc'): z[(y >= 1) & (z > 4)] = b.min()
		else: z[z < 0.5] = b.max()
	
	# reverse axes for footprint
	if(coordinates == 'phi') & (physical > 0):
		z = z[:,::-1]		# reverse x-axis
		if(target== 'in') & (physical == 1):
			z = z[::-1,:]	# reverse y-axis
			if(typeOfPlot == 'imshow'):
				x0, y0 = meshgrid(linspace(x.min(), x.max(), Npsi), linspace(y.min(), y.max(), Nth))
				z = inter.griddata((x.flatten(),y.flatten()), z.flatten(), (x0,y0), method = 'linear')
				x, y = x0, y0
		
	# correct orientation in RZ plot and footprint
	if(coordinates == 'RZ') | (coordinates == 'phi'): x, y, z = x.T, y.T, z.T
	
	# correct, if last x point is larger than 2pi -> x mod 2pi is very small
	if(x[-1,0] < x[-2,0]): x = x[0:-1,:]; y = y[0:-1,:]; z = z[0:-1,:]

	# --- Pest theta ------------------
	if(coordinates == 'pest'):
		import sfc_class
		sfc = sfc_class.straight_field_line_coordinates(shot, time)
		xp = sfc.ev(x,y)
		for i in xrange(Npsi):
			spline = inter.UnivariateSpline(xp[:,i], z[:,i], s = 0)
			z[:,i] = spline(x[:,i])

	# --- layout ----------------------
	#rcParams['text.latex.preamble'] = [r'\usepackage{times}']#\usepackage{amsmath}
	rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
	rcParams['font.serif'] = ['Times', 'Times New Roman']
	font = {'family' : 'serif',
			'weight' : 'normal', # normal
			'size'   : 22} #18
	matplotlib.rc('font', **font)
	latexFontSize = 24

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
	
	C_label_size = latexFontSize

	# --- make plot -------------------
	F = figure(figsize = (width,height))
	if(typeOfPlot == 'contourf'):
		cs = contourf(x, y, z, b, cmap = usecolormap, extend = 'both')
	else:
		cs = imshow(z.T, extent = [x.min(), x.max(), y.min(), y.max()], cmap = usecolormap, origin = 'lower', vmin = b.min(), vmax = b.max(), aspect = 'auto', interpolation = 'bilinear')

	# --- Axes  -----------------------
	if(coordinates == 'RZ'): 
		axes().set_aspect('equal')
		xlim(x.min(),x.max())
		ylim(y.min(),y.max())
		locator_params(axis = 'x', nbins = 5)
	elif(coordinates == 'phi'): 
		xlim(x.min(),x.max())
		ylim(y.min(),y.max())
	elif(coordinates == 'pest'): 
		ylim(y.min(), min([y.max(), 1.0]))
	if(coordinates == 'phi') & (target== 'in') & (physical == 2): 
		gca().invert_yaxis()

	xlabel(xLabel, size = xlabel_size)
	ylabel(yLabel, size = ylabel_size, labelpad = 10)
	
	# adjust axis position on screen; leaves savefile unaffected
	pos = gca().get_position().get_points()	# get bounding-box
	pos[0,0] += 0.025; pos[1,0] += 0.025	# x-position: left; right
	pos[0,1] += 0.05; pos[1,1] += 0.05		# y-position: bottom; top
	gca().set_position(matplotlib.transforms.Bbox(pos))
	
	# --- Colorbar --------------------
	# set ticks
	if(myticks == None):
		steps = (b.max() - b.min())/10.0
		digit = int(floor(log10(steps)))
		factor = 10**digit
		steps = int(ceil(steps/factor))*factor
		bmin = int(ceil(b.min()/factor))*factor
		myticks = arange(bmin, b.max()+steps, steps)
		if(coordinates == 'RZ') | (coordinates == 'phi'): 
			if(what == 'psimin'): 
				myticks = myticks[myticks <= b.max()]
				myticks[-1] = b.max()
			else:
				myticks[0] = b.min()
	
		myticks[abs(myticks) < 1e-10] = 0
	
	# If the "extend" argument is given, contourf sets the data limits to some odd extension of the actual data.
	# Resetting the data limits, after plotting the contours, forces set_over & set_under colors to show,
	# once the colorbar is called after the reset.
	if(typeOfPlot == 'contourf'): cs.set_clim(b.min(),b.max())
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		if(what == 'psimin'): cs.cmap.set_over('w')
		else: cs.cmap.set_under('w')

	# show colorbar
	C = colorbar(cs, pad = 0.01, extend = 'both', format = '%.3g', ticks = myticks)
	C.set_label(C_label, rotation = 270, size = C_label_size)
	
	# add SOL label in RZ plot and footprint
	if((coordinates == 'RZ') | (coordinates == 'phi')) & (showSOL): 
		myticklabels = [item.get_text() for item in C.ax.get_yticklabels()]
		if(what == 'psimin'): myticklabels[-1] = 'SOL'
		else: myticklabels[0] = 'SOL'
		C.ax.set_yticklabels(myticklabels)
				
	# --- additional plots ------------
	if(coordinates == 'RZ'): # plot wall
		wall = loadtxt(HOME + '/c++/d3d/wall.dat')
		plot(wall[:,2], wall[:,3], 'k--', linewidth = 2)
		
	if(coordinates == 'phi') & (target== 'in'):	# plot tile limit
		if(physical == 1): plot([x.min(), x.max()], [-1.223, -1.223], 'k--', linewidth = 1.5)
		else: plot([x.min(), x.max()], [0, 0], 'k--', linewidth = 1.5)
		
	# --- Title -----------------------
	if(Title == True): title(str(shot) + ',  ' + str(time) + 'ms', size = 18)
	elif not (Title == None): title(Title, size = 18)

	# --- save Figure to file ---------
	if not (tag == None): printme = True
	if printme: 
		if(tag == None): F.savefig(path + filename[0:-4] + '.' + graphic, dpi = (300), bbox_inches = 'tight')
		else: F.savefig(path + filename[0:-4] + '_' + tag + '.' + graphic, dpi = (300), bbox_inches = 'tight')











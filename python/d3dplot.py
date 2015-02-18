import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap

HOME = os.getenv('HOME')

def d3dplot(pathname, printme = False, coordinates = 'psi', what = 'psimin', machine = 'd3d', 
			tag = None, graphic = 'png', physical = 1, b = None, N = 60, Title = None,
			typeOfPlot = 'contourf', xlimit = None, ylimit = None, figwidth = None, figheight = None):
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
		else: 
			print 'cannot identify target'
			return
	else: target = None
		
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
	if ((coordinates == 'psi') | (coordinates == 'pest')) & (machine == 'd3d'):
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
		if(machine == 'd3d'):
			xLabel = '$\\varphi$ $\\mathrm{[rad]}$'
			yLabel = 't'
			if(physical > 0):
				x = (360 - x*180.0/np.pi)[:,::-1]	# reverse x
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
		elif(machine == 'iter'):
			xLabel = '$\\varphi$ $\\mathrm{[rad]}$'
			yLabel = 't [cm]'
			if(physical > 0):
				x = x*180.0/np.pi
				xLabel = '$\\phi$ $\\mathrm{[deg]}$'			
	else:
		print 'coordinates: Unknown input'
		return

	# --- set z data ------------------
	if(what == 'Lc'):
		if(b == None): 
			if(machine == 'd3d'): b = np.linspace(0.075,0.4,N)
			elif(machine == 'iter'): b = np.linspace(0.22,1.8,N)
		z = Lc
		C_label = '$L_{c}$ $\\mathrm{[km]}$'
		#usecolormap = cm.jet	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.jet._segmentdata
		
	elif(what == 'psimin'):
		if(b == None): b = np.linspace(0.88,1.02,N)
		z = psimin
		C_label = '$\\psi_{Min}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.jet_r._segmentdata

	elif(what == 'psimax'):
		if(b == None): b = np.linspace(0.88,1.02,N)
		z = psimax
		C_label = '$\\psi_{Max}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.jet_r._segmentdata
		
	elif(what == 'psiav'):
		if(b == None): b = np.linspace(0.88,1.02,N)
		z = psiav
		C_label = '$\\psi_{av}$'
		#usecolormap = cm.jet_r	#  'myjet', 'jet' or cm.jet, cm.jet_r
		cdict = plt.cm.jet_r._segmentdata
		
	else: 
		print 'what: Unknown input'
		return
		
	usecolormap = LinearSegmentedColormap('my_cmap', cdict, len(b))
	
	# correct for PFR
	if(machine == 'd3d'): Lcmin = 0.075
	elif(machine == 'iter'): Lcmin = 0.22
	
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		if not (what == 'Lc'):
			z[Lc < Lcmin] = 1.01*b.max()
	elif not (coordinates == 'phi'): 
		if(what == 'Lc'): z[(y >= 1) & (z >= 4)] = b.min()
		else: z[z < 0.5] = b.max()
	
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
	plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Bitstream Vera Sans']
	plt.rcParams['font.serif'] = ['Times', 'Times New Roman', 'Bitstream Vera Serif']
	font = {'family' : 'sans-serif',
			'weight' : 'normal', # normal
			'size'   : 22} #18
	plt.matplotlib.rc('font', **font)
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
		
	if figwidth == None: figwidth = width
	if figheight == None: figheight = height
	
	C_label_size = latexFontSize

	# --- make plot -------------------
	F = plt.figure(figsize = (figwidth,figheight))
	if(typeOfPlot == 'contourf'):
		cs = plt.contourf(x, y, z, b, cmap = usecolormap, extend = 'both')
	else:
		cs = plt.imshow(z.T, extent = [x.min(), x.max(), y.min(), y.max()], cmap = usecolormap, origin = 'lower', vmin = b.min(), vmax = b.max(), aspect = 'auto', interpolation = 'bilinear')

	# --- Axes  -----------------------
	if(coordinates == 'RZ'): 
		plt.axes().set_aspect('equal')
		plt.xlim(x.min(),x.max())
		plt.ylim(y.min(),y.max())
		plt.locator_params(axis = 'x', nbins = 5)
	elif(coordinates == 'phi'): 
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
		if(what == 'psimin'): 
			myticks = myticks[myticks <= b.max()]
			myticks[-1] = b.max()
		else:
			myticks[0] = b.min()
	
	myticks[np.abs(myticks) < 1e-10] = 0
	
	# If the "extend" argument is given, contourf sets the data limits to some odd extension of the actual data.
	# Resetting the data limits, after plotting the contours, forces set_over & set_under colors to show,
	# once the colorbar is called after the reset.
	if(typeOfPlot == 'contourf'): cs.set_clim(b.min(),b.max())
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		if(what == 'psimin'): cs.cmap.set_over('w')
		else: cs.cmap.set_under('w')

	# show colorbar
	C = plt.colorbar(cs, pad = 0.01, extend = 'both', format = '%.3g', ticks = myticks)
	C.set_label(C_label, rotation = 270, size = C_label_size)
	
	# add SOL label in RZ plot and footprint
	if(coordinates == 'RZ') | (coordinates == 'phi'): 
		myticklabels = [item.get_text() for item in C.ax.get_yticklabels()]
		if(what == 'psimin'): myticklabels[-1] = 'SOL'
		else: myticklabels[0] = 'SOL'
		C.ax.set_yticklabels(myticklabels)
				
	# --- additional plots ------------
	if(coordinates == 'RZ'): # plot wall
		if(machine == 'd3d'): 
			wall = np.loadtxt(HOME + '/c++/d3d/wall.dat')
			plt.plot(wall[:,2], wall[:,3], 'k--', linewidth = 2)
		elif(machine == 'iter'): 
			wall = np.loadtxt(HOME + '/c++/iter/wall.dat')
			plt.plot(wall[:,0], wall[:,1], 'k--', linewidth = 2)
		
	if(coordinates == 'phi') & (target== 'in'):	# plot tile limit
		if(physical == 1) & (machine == 'd3d'): plt.plot([x.min(), x.max()], [-1.223, -1.223], 'k--', linewidth = 1.5)
		else: plt.plot([x.min(), x.max()], [0, 0], 'k--', linewidth = 1.5)
		
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
	Title = None
	toP = 'contourf'
	figwidth = None
	figheight = None

	import sys, getopt
	opts, args = getopt.gnu_getopt(sys.argv[1:], "hc:w:m:pg:t:P:N:T:iW:H:", ["help", "coord=", "what=", 
																 "machine=", "printme", "graphic=", 
																 "tag=", "physical=", "Title=", "imshow", 	 
																 "figwidth=", "figheight="])
	for o, a in opts:
		if o in ("-h", "--help"):
			print "usage: d3dplot.py [-h] [-c COORDINATES] [-w WHAT] [-m MACHINE] [-p]"
			print "                  [-g GRAPHIC] [-t TAG] [-P PHYSICAL] [-N N] [-T TITLE] [-i]"
			print "                  [-W FIGWIDTH] [-H FIGHEIGHT]"
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
			print "  -T, --Title <Arg>     Figure title, default = None"
			print "  -i, --imshow          Use imshow instead of contourf (default)"
			print "  -W, --figwidth <Arg>  Force width of figure from default to <Arg>"
			print "  -H, --figheight <Arg> Force height of figure from default to <Arg>"
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
		elif o in ("-T", "--Title"):
			Title = a
		elif o in ("-i", "--imshow"):
			toP = 'imshow'
		elif o in ("-W", "--figwidth"):
			figwidth = float(a)
		elif o in ("-H", "--figheight"):
			figheight = float(a)
		else:
			raise AssertionError("unknown option")
		
	d3dplot(args[0], printme = printme, coordinates = coordinates, what = what, machine = machine, 
			tag = tag, graphic = graphic, physical = physical, b = None, N = N, Title = Title,
			typeOfPlot = toP, xlimit = None, ylimit = None, figwidth = figwidth, figheight = figheight)

	plt.show()


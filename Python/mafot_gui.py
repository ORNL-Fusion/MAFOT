import os
import Tkinter as tk
import ttk
import socket
from subprocess import call

HOME = os.getenv('HOME')
HOST = socket.gethostname()
MPIRUN = 'mpirun'

# -------------------------------------------------------------------------------------------------------------
# --- common input --------------------------------------------------------------------------------------------
class Common_gui:
	def __init__(self, frame):
		self.frame = frame
	
		# --- read parameterfile, if it is there ---
		self.ControlFileFound = True
		if os.path.isfile('_plot.dat'):
			shot, time, gpath, _ = self.readControlFile('_plot.dat')
		elif os.path.isfile('_fix.dat'):
			shot, time, gpath, _ = self.readControlFile('_fix.dat')
		elif os.path.isfile('_inner.dat'):
			shot, time, gpath, _ = self.readControlFile('_inner.dat')
		elif os.path.isfile('_outer.dat'):
			shot, time, gpath, _ = self.readControlFile('_outer.dat')
		elif os.path.isfile('_shelf.dat'):
			shot, time, gpath, _ = self.readControlFile('_shelf.dat')
		elif os.path.isfile('_lam.dat'):
			shot, time, gpath, _ = self.readControlFile('_lam.dat')
		elif os.path.isfile('_lam_psi.dat'):
			shot, time, gpath, _ = self.readControlFile('_lam_psi.dat')
		else:
			shot = ''
			time = ''
			gpath = HOME + '/c++/d3d/gfiles'
			self.ControlFileFound = False
			
		if not (gpath[-1] == '/'): gpath += '/'
		
		okayCommand = frame.register(self.isOkay)
		notOkayCommand = frame.register(self.isNotOkay)
		
		# --- Machine ---
		self.MachFlag = tk.StringVar(); self.MachFlag.set('d3d'); row = 0
		tk.Radiobutton(frame, text = 'DIII-D', variable = self.MachFlag, value = 'd3d', 
			command = self.set_Machine).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'ITER', variable = self.MachFlag, value = 'iter', 
			command = self.set_Machine, state=tk.DISABLED).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'NSTX', variable = self.MachFlag, value = 'nstx', 
			command = self.set_Machine, state=tk.DISABLED).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'MAST', variable = self.MachFlag, value = 'mast', 
			command = self.set_Machine, state=tk.DISABLED).grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Machine").grid(column = 1, row = row, sticky = tk.E )

		# --- Shot number ---
		self.Shot = tk.StringVar(); self.Shot.set(str(shot)); row += 1
		tk.Entry(frame, width = 7, textvariable = self.Shot).grid(column = 2, row = row, sticky = tk.W + tk.E)
		tk.Label(frame, text = "Shot").grid(column = 1, row = row, sticky = tk.E)

		# --- Time ---
		self.Time = tk.StringVar(); self.Time.set(str(time)); #row += 1
		tk.Entry(frame, width = 7, textvariable = self.Time).grid(column = 4, row = row, sticky = tk.W + tk.E)
		tk.Label(frame, text = "Time").grid(column = 3, row = row, sticky = tk.E)

		# --- gPath ---
		self.gPath = tk.StringVar(); self.gPath.set(gpath); row += 1
		tk.Entry(frame, width = 7, textvariable = self.gPath).grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)
		tk.Label(frame, text = "Path to g-file").grid(column = 1, row = row, sticky = tk.E)
		
		# --- Tag ---
		self.tag = tk.StringVar(); self.tag.set(''); row += 1; self.row_write = row
		self.tag_entry = tk.Entry(frame, width = 7, textvariable = self.tag, validate='all',
        	validatecommand = (okayCommand, '%P'), invalidcommand = notOkayCommand)
		self.tag_entry.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)
		self.tag_label = tk.Label(frame, text = "File Tag")
		self.tag_label.grid(column = 1, row = row, sticky = tk.E)
		
		# --- cwd ---
		self.path = tk.StringVar(); self.path.set('./'); row += 1
		self.path_entry = tk.Entry(frame, width = 7, textvariable = self.path)
		self.path_entry.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)
		self.path_label = tk.Label(frame, text = "Working Dir")
		self.path_label.grid(column = 1, row = row, sticky = tk.E)
		
		# --- invisible separator ---
		row += 1
		separator = tk.Label(frame, text = "")
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W)

		# --- adjust style for all ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		
		self.tag_entry.focus()
		
		# --- default setup is for D3D ---
		self.d3d_tabs(frame)
	
	
	# --- validity check of Tag ---
	# returns True if tag is alphanumeric or has _ + - as special chars
	# else returns False
	def isOkay(self, tag):
		if(len(tag) == 0): return True
		else: return tag.translate(None, '_+-').isalnum()
		
	# prints error Message, if isOkay == False
	def isNotOkay(self):
		print 'Warning: Invalid Input Character in File Tag. Only alphanumeric and + - _ are allowed'
				
		
	# --- Tabs for DIII-D ---
	def d3d_tabs(self, frame):
		if not self.ControlFileFound:
			self.gPath.set(HOME + '/c++/d3d/gfiles/')

		# --- start tabs ---
		self.nb = ttk.Notebook(frame)
		self.nb.grid(column = 1, row = 7, columnspan = 10)
		
		# --- Tab 0: set plot tab ---
		self.tabframe1 = tk.Frame(frame)
		self.tabframe1.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe1.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe1, text = 'dtplot')
		dtplot_gui(self.tabframe1, self)

		# --- Tab 1: set fix tab ---
		self.tabframe2 = tk.Frame(frame)
		self.tabframe2.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe2.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe2, text = 'dtfix')
		dtfix_gui(self.tabframe2, self)
	
		# --- Tab 2: set man tab ---
		self.tabframe3 = tk.Frame(frame)
		self.tabframe3.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe3.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe3, text = 'dtman')
		dtman_gui(self.tabframe3, self)
	
		# --- Tab 3: set foot tab ---
		self.tabframe4 = tk.Frame(frame)
		self.tabframe4.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe4.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe4, text = 'dtfoot')
		dtfoot_gui(self.tabframe4, self)

		# --- Tab 4: set laminar tab ---
		self.tabframe5 = tk.Frame(frame)
		self.tabframe5.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe5.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe5, text = 'dtlaminar')	
		dtlam_gui(self.tabframe5, self)
		
		# --- Tab 5: set info tab ---
		self.tabframe6 = tk.Frame(frame)
		self.tabframe6.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe6.rowconfigure(24, weight=1)
		self.nb.add(self.tabframe6, text = 'Info')
		info_gui(self.tabframe6)
		
		
	# --- Tabs for ITER ---
	def iter_tabs(self, frame):
		if not self.ControlFileFound:
			self.gPath.set(HOME + '/c++/iter/gfiles/')
			
		# --- start tabs ---
		self.nb = ttk.Notebook(frame)
		self.nb.grid(column = 1, row = 7, columnspan = 10)
		
		# --- Tab 0: set plot tab ---
		self.tabframe1 = tk.Frame(frame)
		self.tabframe1.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe1.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe1, text = 'iterplot')
		
		# --- Tab 1: set fix tab ---
		self.tabframe2 = tk.Frame(frame)
		self.tabframe2.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe2.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe2, text = 'iterfix')
	
		# --- Tab 2: set man tab ---
		self.tabframe3 = tk.Frame(frame)
		self.tabframe3.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe3.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe3, text = 'iterman')
	
		# --- Tab 3: set foot tab ---
		self.tabframe4 = tk.Frame(frame)
		self.tabframe4.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe4.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe4, text = 'iterfoot')

		# --- Tab 4: set laminar tab ---
		self.tabframe5 = tk.Frame(frame)
		self.tabframe5.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe5.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe5, text = 'iterlaminar')	
		
		# --- Tab 5: set info tab ---
		self.tabframe6 = tk.Frame(frame)
		self.tabframe6.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe6.rowconfigure(24, weight=1)
		self.nb.add(self.tabframe6, text = 'Info')
		info_gui(self.tabframe6)


	# --- Tabs for NSTX ---
	def nstx_tabs(self, frame):
		if not self.ControlFileFound:
			self.gPath.set(HOME + '/c++/nstx/gfiles/')

		# --- start tabs ---
		self.nb = ttk.Notebook(frame)
		self.nb.grid(column = 1, row = 7, columnspan = 10)
		
		# --- Tab 0: set plot tab ---
		self.tabframe1 = tk.Frame(frame)
		self.tabframe1.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe1.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe1, text = 'nstxplot')
		
		# --- Tab 1: set fix tab ---
		self.tabframe2 = tk.Frame(frame)
		self.tabframe2.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe2.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe2, text = 'nstxfix')
	
		# --- Tab 2: set man tab ---
		self.tabframe3 = tk.Frame(frame)
		self.tabframe3.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe3.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe3, text = 'nstxman')
	
		# --- Tab 3: set foot tab ---
		self.tabframe4 = tk.Frame(frame)
		self.tabframe4.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe4.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe4, text = 'nstxfoot')

		# --- Tab 4: set laminar tab ---
		self.tabframe5 = tk.Frame(frame)
		self.tabframe5.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe5.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe5, text = 'nstxlaminar')	
		
		# --- Tab 5: set info tab ---
		self.tabframe6 = tk.Frame(frame)
		self.tabframe6.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe6.rowconfigure(24, weight=1)
		self.nb.add(self.tabframe6, text = 'Info')
		info_gui(self.tabframe6)


	# --- Tabs for MAST ---
	def mast_tabs(self, frame):
		if not self.ControlFileFound:
			self.gPath.set(HOME + '/c++/mast/gfiles/')

		# --- start tabs ---
		self.nb = ttk.Notebook(frame)
		self.nb.grid(column = 1, row = 7, columnspan = 10)
		
		# --- Tab 0: set plot tab ---
		self.tabframe1 = tk.Frame(frame)
		self.tabframe1.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe1.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe1, text = 'mastplot')
		
		# --- Tab 1: set fix tab ---
		self.tabframe2 = tk.Frame(frame)
		self.tabframe2.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe2.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe2, text = 'mastfix')
	
		# --- Tab 2: set man tab ---
		self.tabframe3 = tk.Frame(frame)
		self.tabframe3.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe3.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe3, text = 'mastman')
	
		# --- Tab 3: set foot tab ---
		self.tabframe4 = tk.Frame(frame)
		self.tabframe4.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe4.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe4, text = 'mastfoot')

		# --- Tab 4: set laminar tab ---
		self.tabframe5 = tk.Frame(frame)
		self.tabframe5.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe5.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe5, text = 'mastlaminar')	
		
		# --- Tab 5: set info tab ---
		self.tabframe6 = tk.Frame(frame)
		self.tabframe6.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe6.rowconfigure(24, weight=1)
		self.nb.add(self.tabframe6, text = 'Info')
		info_gui(self.tabframe6)


	# --- set Tabs for different Machines ---
	def set_Machine(self):
		if(self.MachFlag.get() == 'iter'):
			self.nb.forget(self.tabframe1); del self.tabframe1
			self.nb.forget(self.tabframe2); del self.tabframe2
			self.nb.forget(self.tabframe3); del self.tabframe3
			self.nb.forget(self.tabframe4); del self.tabframe4
			self.nb.forget(self.tabframe5); del self.tabframe5
			self.nb.forget(self.tabframe6); del self.tabframe6
			self.iter_tabs(self.frame)
		elif(self.MachFlag.get() == 'nstx'):
			self.nb.forget(self.tabframe1); del self.tabframe1
			self.nb.forget(self.tabframe2); del self.tabframe2
			self.nb.forget(self.tabframe3); del self.tabframe3
			self.nb.forget(self.tabframe4); del self.tabframe4
			self.nb.forget(self.tabframe5); del self.tabframe5
			self.nb.forget(self.tabframe6); del self.tabframe6
			self.nstx_tabs(self.frame)
		elif(self.MachFlag.get() == 'mast'):
			self.nb.forget(self.tabframe1); del self.tabframe1
			self.nb.forget(self.tabframe2); del self.tabframe2
			self.nb.forget(self.tabframe3); del self.tabframe3
			self.nb.forget(self.tabframe4); del self.tabframe4
			self.nb.forget(self.tabframe5); del self.tabframe5
			self.nb.forget(self.tabframe6); del self.tabframe6
			self.mast_tabs(self.frame)
		else:	# d3d
			self.nb.forget(self.tabframe1); del self.tabframe1
			self.nb.forget(self.tabframe2); del self.tabframe2
			self.nb.forget(self.tabframe3); del self.tabframe3
			self.nb.forget(self.tabframe4); del self.tabframe4
			self.nb.forget(self.tabframe5); del self.tabframe5
			self.nb.forget(self.tabframe6); del self.tabframe6
			self.d3d_tabs(self.frame)


	# --- read Control File, and return its values ---
	def readControlFile(self, name):
		with open(name, 'r') as f:
			# Skip first line
			head = f.readline().split()
			# read Shot and Time
			head = f.readline().split()
			shot = int(head[2])
			time = int(head[-1][0:-2])
			# read gPath if there, else this is the first data point
			head = f.readline().split()
			if head[0][0] == '#':
				gpath = head[-1]
				data = []
			else:
				gpath = HOME + '/c++/d3d/gfiles'
				data = [float(head[-1])]
		
			# read data to end of file
			for line in f:
				line = line.split()			
				data.append(float(line[-1]))
			
		return shot, time, gpath, data


# -------------------------------------------------------------------------------------------------------------
# --- dtplot --------------------------------------------------------------------------------------------------
class dtplot_gui:
	def __init__(self, frame, M):
		
		# --- read parameterfile, if it is there ---
		if os.path.isfile('_plot.dat'):
			_,_,_, data = M.readControlFile('_plot.dat')
		else: # defaults
			data = [0, 300, 0.42, 0.57, 0, 0, 40, 0, 1, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 3.141592653589793, 6.283185307179586]
	
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path

		# --- grid  type ---
		# define... 
		self.createFlag = tk.StringVar(); row = 0
		if(data[12] == 0) | (data[12] == 1): self.createFlag.set('polar')
		elif(data[12] == 3) | (data[12] == 4): self.createFlag.set('psi')
		else: self.createFlag.set('RZ')

		# ...and set grid-type RadioButton
		tk.Radiobutton(frame, text = 'RZ', variable = self.createFlag, value = 'RZ', 
			command = self.refresh_grid_labels).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'psi_n', variable = self.createFlag, value = 'psi', 
			command = self.refresh_grid_labels).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Polar', variable = self.createFlag, value = 'polar', 
			command = self.refresh_grid_labels).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Coordinate Type").grid(column = 1, row = row, sticky = tk.E )

		# --- x -> theta or R ---
		row += 1
		self.x_label = tk.Label(frame, text = "")
		self.x_label.grid(column = 1, row = row, sticky = tk.E)
		
		# Min
		self.xmin = tk.StringVar(); 
		if(self.createFlag.get() == 'RZ'): self.xmin.set(str(data[2]))
		else: self.xmin.set(str(data[4]))
		tk.Entry(frame, width = 16, textvariable = self.xmin).grid(column = 2, row = row, columnspan = 2)
		xmin_label = tk.Label(frame, text = "  Min")
		xmin_label.grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.xmax = tk.StringVar();
		if(self.createFlag.get() == 'RZ'): self.xmax.set(str(data[3]))
		else: self.xmax.set(str(data[5]))
		tk.Entry(frame, width = 16, textvariable = self.xmax).grid(column = 4, row = row, columnspan = 2)
		xmax_label = tk.Label(frame, text = "Max")
		xmax_label.grid(column = 4, row = row, sticky = tk.W )

		# Nx -> Nth or NR
		self.Nx = tk.StringVar(); row += 1
		if(self.createFlag.get() == 'polar'): self.Nx.set(str(int(data[0])))
		elif(self.createFlag.get() == 'psi'): self.Nx.set(str(int(data[0])))
		else: self.Nx.set(str(int(data[6])))
		tk.Entry(frame, width = 16, textvariable = self.Nx).grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )
		
		self.pi_text = tk.Text(frame, height= 1, width = 30, bd  = 0, takefocus = 0, bg = frame.cget('bg'), relief = tk.FLAT)
		self.pi_text.grid(column = 4, row = row, columnspan = 2); self.pi_row = row
		self.pi_text.insert(1.0, 'pi = 3.141593 2pi = 6.283185')
		self.pi_text.configure(state = "disabled")
		
		# --- y -> r, psi or Z ---
		row += 1
		self.y_label = tk.Label(frame, text = "")
		self.y_label.grid(column = 1, row = row, sticky = tk.E )
		
		# Min
		self.ymin = tk.StringVar();
		if(self.createFlag.get() == 'RZ'): self.ymin.set(str(data[4]))
		else: self.ymin.set(str(data[2]))
		tk.Entry(frame, width = 16, textvariable = self.ymin).grid(column = 2, row = row, columnspan = 2)
		ymin_label = tk.Label(frame, text = "  Min")
		ymin_label.grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.ymax = tk.StringVar();
		if(self.createFlag.get() == 'RZ'): self.ymax.set(str(data[5]))
		else: self.ymax.set(str(data[3]))
		tk.Entry(frame, width = 16, textvariable = self.ymax).grid(column = 4, row = row, columnspan = 2)
		ymax_label = tk.Label(frame, text = "Max")
		ymax_label.grid(column = 4, row = row, sticky = tk.W )

		# Ny -> Nr, Npsi or NZ
		self.Ny = tk.StringVar(); row += 1
		if(self.createFlag.get() == 'polar'): self.Ny.set(str(int(data[6])))
		elif(self.createFlag.get() == 'psi'): self.Ny.set(str(int(data[6])))
		else: self.Ny.set(str(int(data[0])))
		tk.Entry(frame, width = 16, textvariable = self.Ny).grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )

		self.refresh_grid_labels()
		
		# --- toroidal turns ---
		self.itt = tk.StringVar(); self.itt.set(str(int(data[1]))); row += 1
		tk.Entry(frame, width = 7, textvariable = self.itt).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "tor. Iterations").grid(column = 1, row = row, sticky = tk.E )

		# --- phistart ---
		self.phistart = tk.StringVar(); self.phistart.set(str(data[7])); row += 1
		tk.Entry(frame, width = 7, textvariable = self.phistart).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "tor. Angle [deg]").grid(column = 1, row = row, sticky = tk.E )
		tk.Label(frame, text = "For Machine coord. use negative angles").grid(column = 3, row = row, columnspan = 3, sticky = tk.W )

		# --- MapDirection ---
		self.MapDirection = tk.IntVar(); self.MapDirection.set(int(data[8])); row += 1
		tk.Radiobutton(frame, text = '+1', variable = self.MapDirection, value = 1).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = '-1', variable = self.MapDirection, value = -1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'both', variable = self.MapDirection, value = 0).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Map Direction").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- coils ---
		self.useFcoil = tk.IntVar(); self.useFcoil.set(int(data[13])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useFcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useFcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use F-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useCcoil = tk.IntVar(); self.useCcoil.set(int(data[14])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useCcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useCcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use C-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useIcoil = tk.IntVar(); self.useIcoil.set(int(data[15])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useIcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useIcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use I-coil").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator2 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator2.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- M3DC1 ---
		self.useM3DC1 = tk.IntVar(); self.useM3DC1.set(int(data[10])); row += 1
		tk.Radiobutton(frame, text = 'g-file Vacuum', variable = self.useM3DC1, value = -1,
			command = self.activate_response).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Equilibrium', variable = self.useM3DC1, value = 0,
			command = self.activate_response).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Perturbation', variable = self.useM3DC1, value = 1,
			command = self.activate_response).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Eq + Pert', variable = self.useM3DC1, value = 2,
			command = self.activate_response).grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use M3DC1").grid(column = 1, row = row, sticky = tk.E )

		self.response = tk.IntVar(); self.response.set(int(data[9])); row += 1
		self.response_R1 = tk.Radiobutton(frame, text = 'Off', variable = self.response, value = 0)
		self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.response_R2 = tk.Radiobutton(frame, text = 'On', variable = self.response, value = 1)
		self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Plasma response").grid(column = 1, row = row, sticky = tk.E )
		
		self.activate_response()

		# --- separator ---
		row += 1
		separator3 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator3.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- particles ---
		self.sigma = tk.IntVar(); self.sigma.set(int(data[16])); row += 1
		tk.Radiobutton(frame, text = 'Field lines', variable = self.sigma, value = 0, 
			command = self.show_particle_params).grid(column = 2, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Co pass', variable = self.sigma, value = 1,
			command = self.show_particle_params).grid(column = 3, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Counter pass', variable = self.sigma, value = -1,
			command = self.show_particle_params).grid(column = 4, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Label(frame, text = "Orbits").grid(column = 1, row = row, sticky = tk.E + tk.N)

		self.charge = tk.IntVar(); self.charge.set(int(data[17])); row += 1; self.row_particle = row
		self.charge_R1 = tk.Radiobutton(frame, text = 'Electrons', variable = self.charge, value = -1)
		self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.charge_R2 = tk.Radiobutton(frame, text = 'Ions', variable = self.charge, value = 1)
		self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.charge_label = tk.Label(frame, text = "Species")
		self.charge_label.grid(column = 1, row = row, sticky = tk.E )

		self.Ekin = tk.StringVar(); self.Ekin.set(str(data[18])); row += 1
		self.Ekin_entry = tk.Entry(frame, width = 7, textvariable = self.Ekin)
		self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.Ekin_label = tk.Label(frame, text = "kin Energy [keV]")
		self.Ekin_label.grid(column = 1, row = row, sticky = tk.E )

		self.Lambda = tk.StringVar(); self.Lambda.set(str(data[19]));
		self.Lambda_entry = tk.Entry(frame, width = 7, textvariable = self.Lambda)
		self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.Lambda_label = tk.Label(frame, text = "Energy ratio")	
		self.Lambda_label.grid(column = 3, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator4 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator4.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- Filament ---
		self.useFilament = tk.StringVar(); self.useFilament.set(str(int(data[20]))); row += 1
		tk.Entry(frame, width = 7, textvariable = self.useFilament).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "# current Filaments").grid(column = 1, row = row, sticky = tk.E )

		# --- run Button ---
		row += 1
		runButton = tk.Button(frame, text = "Run dtplot", command = self.run_funct)
		runButton.grid(column = 1, row = row, columnspan = 5, sticky = tk.W + tk.E)

		# --- adjust style for all ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.show_particle_params()


	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		if(self.Shot.get() == '') | (self.Time.get() == ''):
			print 'You must enter a Shot # and a Time'
			return

		cwd = os.getcwd()
		if(self.path.get()[0] == '.'):
			path = cwd + self.path.get()[1::]
		else: 
			path = self.path.get()
		if not (path[-1] == '/'): path += '/'
		chk = not (cwd + '/' == path)
		
		if chk: os.chdir(path)	
		self.writeControlFile('_plot.dat')
		if(HOST == 'head.cluster'):		# Drop Cluster
			self.write_qsub_file(self.tag.get())
			call('qsub run_job', shell = True)
		else:
			call('dtplot _plot.dat ' + self.tag.get() + ' &', shell = True)		
		if chk: os.chdir(cwd)


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, tag):
		with open('run_job', 'w') as f:
			f.write('#$ -N P' + tag + '\n')
			f.write('#$ -cwd \n')
			f.write('#$ -o /home/wingen/work/batch.out \n')
			f.write('#$ -e /home/wingen/work/batch.err \n')
			f.write('#$ -S /bin/bash \n')
			f.write('#$ -V \n')
			f.write('#$ -q all.q \n')
			f.write('dtplot _plot.dat ' + tag + '\n')
			

	# --- Change Labels on grid variables, depending on createFlag ---
	def refresh_grid_labels(self):
		if(self.createFlag.get() == 'polar'):
			self.x_label.configure(text = "theta [rad]")
			self.y_label.configure(text = "r [m]")
			self.pi_text.grid(column = 4, row = self.pi_row, columnspan = 2)
		elif(self.createFlag.get() == 'psi'):
			self.x_label.configure(text = "theta [rad]")
			self.y_label.configure(text = "psi_n")
			self.pi_text.grid(column = 4, row = self.pi_row, columnspan = 2)
		else:
			self.x_label.configure(text = "R [m]")
			self.y_label.configure(text = "Z [m]")
			self.pi_text.grid_forget()
			
			
	# --- turn on/off response radiobutton ---
	def activate_response(self):
		if(self.useM3DC1.get() == -1):
			self.response_R1.configure(state=tk.DISABLED)
			self.response_R2.configure(state=tk.DISABLED)
		else:
			self.response_R1.configure(state=tk.NORMAL)
			self.response_R2.configure(state=tk.NORMAL)

	
	# --- Show or Hide Particle Options, depending on sigma ---
	def show_particle_params(self):
		if not (self.sigma.get() == 0):
			row = self.row_particle
			self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			row = self.row_particle + 1
			self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Ekin_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Lambda_label.grid(column = 3, row = row, sticky = tk.E , padx=5, pady=5)
		else:
			self.charge_R1.grid_forget()
			self.charge_R2.grid_forget()
			self.charge_label.grid_forget()
			self.Ekin_entry.grid_forget()
			self.Ekin_label.grid_forget()
			self.Lambda_entry.grid_forget()
			self.Lambda_label.grid_forget()
	
	
	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		with open(name, 'w') as f:
			f.write('# Parameterfile for DIII-D Programs\n')
			f.write('# Shot: ' + self.Shot.get() + '\tTime: ' + self.Time.get() + 'ms\n')
			f.write('# Path: ' + self.gPath.get() + '\n')
			
			if(self.createFlag.get() == 'polar'): f.write('free_Parameter=\t0\n')
			elif(self.createFlag.get() == 'psi'): f.write('Nth=\t' + self.Nx.get() + '\n')
			else: f.write('NZ=\t' + self.Ny.get() + '\n')
			
			f.write('itt=\t' + self.itt.get() + '\n')
			
			if(self.createFlag.get() == 'polar'):
				f.write('rmin=\t' + self.ymin.get() + '\n')
				f.write('rmax=\t' + self.ymax.get() + '\n')
				f.write('thmin=\t' + self.xmin.get() + '\n')
				f.write('thmax=\t' + self.xmax.get() + '\n')
				f.write('N=\t' + self.Ny.get() + '\n')
			elif(self.createFlag.get() == 'psi'):
				f.write('psimin=\t' + self.ymin.get() + '\n')
				f.write('psimax=\t' + self.ymax.get() + '\n')
				f.write('thmin=\t' + self.xmin.get() + '\n')
				f.write('thmax=\t' + self.xmax.get() + '\n')
				f.write('Npsi=\t' + self.Ny.get() + '\n')
			else:
				f.write('Rmin=\t' + self.xmin.get() + '\n')
				f.write('Rmax=\t' + self.xmax.get() + '\n')
				f.write('Zmin=\t' + self.ymin.get() + '\n')
				f.write('Zmax=\t' + self.ymax.get() + '\n')
				f.write('NR=\t' + self.Nx.get() + '\n')

			f.write('phistart(deg)=\t' + self.phistart.get() + '\n')
			f.write('MapDirection=\t' + str(self.MapDirection.get()) + '\n')
			f.write('PlasmaResponse(0=no,1=yes)=\t' + str(self.response.get()) + '\n')
			f.write('Field(-1=M3D-C1_off,0=Eq,1=I-coil,2=both)=\t' + str(self.useM3DC1.get()) + '\n')
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')
			
			if(self.createFlag.get() == 'polar'): f.write('createPoints(0=setr,3=setpsi,5=setR)=\t0\n')
			elif(self.createFlag.get() == 'psi'): f.write('createPoints(0=setr,3=setpsi,5=setR)=\t3\n')
			else: f.write('createPoints(0=setr,3=setpsi,5=setR)=\t5\n')
			
			f.write('useFcoil(0=no,1=yes)=\t' + str(self.useFcoil.get()) + '\n')
			f.write('useCcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
			f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
			f.write('ParticleDirection(1=pass,-1=co-pass,0=field-lines)=\t' + str(self.sigma.get()) + '\n')
			f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge.get()) + '\n')
			f.write('Ekin[keV]=\t' + self.Ekin.get() + '\n')
			f.write('lambda=\t' + self.Lambda.get() + '\n')
			f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')
		
		
# -------------------------------------------------------------------------------------------------------------
# --- dtfix ---------------------------------------------------------------------------------------------------
class dtfix_gui:
	def __init__(self, frame, M):
		
		# --- read parameterfile, if it is there ---
		if os.path.isfile('_fix.dat'):
			_,_,_, data = M.readControlFile('_fix.dat')
		else: # defaults
			data = [1e-4, 0, 1, 1.3, 4.1, 4.6, 900, 0, 1, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 3.141592653589793, 6.283185307179586]
	
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		
		# --- shift ---
		self.shift = data[0]
		
		# --- MapDirection ---
		self.MapDirection = int(data[8]); 

		# --- Period of Hyperbolic point ---
		# define... 
		self.HypRPt = tk.IntVar(); self.HypRPt.set(2);
		self.HypPt = tk.StringVar(); self.HypPt.set(str(1)); row = 0 
		
		# ..., set RadioButton for Separatrix
		tk.Radiobutton(frame, text = 'lower X', variable = self.HypRPt, value = 1,
			command = self.activate_entrys).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'upper X', variable = self.HypRPt, value = -1,
			command = self.activate_entrys).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'manual', variable = self.HypRPt, value = 2, 
			command = self.activate_entrys).grid(column = 4, row = row, sticky = tk.W + tk.E )
		# ... and entry for others
		self.period_entry = tk.Entry(frame, width = 4, textvariable = self.HypPt)
		self.period_entry.grid(column = 5, row = row, sticky = tk.W)
		tk.Label(frame, text = "Period").grid(column = 1, row = row, sticky = tk.E )

		# --- x -> theta ---
		row += 1
		tk.Label(frame, text = "theta [rad]").grid(column = 1, row = row, sticky = tk.E)
		
		# Min
		self.xmin = tk.StringVar(); self.xmin.set(str(data[4]))
		self.xmin_entry = tk.Entry(frame, width = 16, textvariable = self.xmin)
		self.xmin_entry.grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "  Min").grid(column = 2, row = row, sticky = tk.W )
		
		# Max
		self.xmax = tk.StringVar(); self.xmax.set(str(data[5]))
		self.xmax_entry = tk.Entry(frame, width = 16, textvariable = self.xmax)
		self.xmax_entry.grid(column = 4, row = row, columnspan = 2)
		tk.Label(frame, text = "Max").grid(column = 4, row = row, sticky = tk.W )

		# Nx -> Nth
		self.Nx = tk.StringVar();self.Nx.set(str(int(data[6]**0.5))); row += 1
		self.Nx_entry = tk.Entry(frame, width = 16, textvariable = self.Nx)
		self.Nx_entry.grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )
		
		self.pi_text = tk.Text(frame, height= 1, width = 30, bd  = 0, takefocus = 0, bg = frame.cget('bg'), relief = tk.FLAT)
		self.pi_text.grid(column = 4, row = row, columnspan = 2); self.pi_row = row
		self.pi_text.insert(1.0, 'pi = 3.141593 2pi = 6.283185')
		self.pi_text.configure(state = "disabled")
		
		# --- y -> r ---
		row += 1
		tk.Label(frame, text = "r [m]").grid(column = 1, row = row, sticky = tk.E )
		
		# Min
		self.ymin = tk.StringVar(); self.ymin.set(str(data[2]))
		self.ymin_entry = tk.Entry(frame, width = 16, textvariable = self.ymin)
		self.ymin_entry.grid(column = 2, row = row, columnspan = 2 )
		tk.Label(frame, text = "  Min").grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.ymax = tk.StringVar(); self.ymax.set(str(data[3]))
		self.ymax_entry = tk.Entry(frame, width = 16, textvariable = self.ymax)
		self.ymax_entry.grid(column = 4, row = row, columnspan = 2 )
		tk.Label(frame, text = "Max").grid(column = 4, row = row, sticky = tk.W )

		# Ny -> Nr
		self.Ny = tk.StringVar();self.Ny.set(str(int(data[6]**0.5))); row += 1
		self.Ny_entry = tk.Entry(frame, width = 16, textvariable = self.Ny)
		self.Ny_entry.grid(column = 2, row = row, columnspan = 2 )
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )

		self.activate_entrys
		
		# --- invisible separator ---
		row += 1
		tk.Label(frame, text = "").grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)

		# --- phistart ---
		self.phistart = tk.StringVar(); self.phistart.set(str(data[7])); row += 1
		tk.Entry(frame, width = 7, textvariable = self.phistart).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "tor. Angle [deg]").grid(column = 1, row = row, sticky = tk.E )
		tk.Label(frame, text = "For Machine coord. use negative angles").grid(column = 3, row = row, columnspan = 3, sticky = tk.W )

		# --- invisible separator ---
		row += 1
		tk.Label(frame, text = "").grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)

		# --- separator ---
		row += 1
		separator = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- coils ---
		self.useFcoil = tk.IntVar(); self.useFcoil.set(int(data[13])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useFcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useFcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use F-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useCcoil = tk.IntVar(); self.useCcoil.set(int(data[14])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useCcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useCcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use C-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useIcoil = tk.IntVar(); self.useIcoil.set(int(data[15])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useIcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useIcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use I-coil").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator2 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator2.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- M3DC1 ---
		self.useM3DC1 = tk.IntVar(); self.useM3DC1.set(int(data[10])); row += 1
		tk.Radiobutton(frame, text = 'g-file Vacuum', variable = self.useM3DC1, value = -1,
			command = self.activate_response).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Equilibrium', variable = self.useM3DC1, value = 0,
			command = self.activate_response).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Perturbation', variable = self.useM3DC1, value = 1,
			command = self.activate_response).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Eq + Pert', variable = self.useM3DC1, value = 2,
			command = self.activate_response).grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use M3DC1").grid(column = 1, row = row, sticky = tk.E )

		self.response = tk.IntVar(); self.response.set(int(data[9])); row += 1
		self.response_R1 = tk.Radiobutton(frame, text = 'Off', variable = self.response, value = 0)
		self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.response_R2 = tk.Radiobutton(frame, text = 'On', variable = self.response, value = 1)
		self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Plasma response").grid(column = 1, row = row, sticky = tk.E )
		
		self.activate_response()

		# --- separator ---
		row += 1
		separator3 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator3.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- particles ---
		self.sigma = tk.IntVar(); self.sigma.set(int(data[16])); row += 1
		tk.Radiobutton(frame, text = 'Field lines', variable = self.sigma, value = 0, 
			command = self.show_particle_params).grid(column = 2, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Co pass', variable = self.sigma, value = 1,
			command = self.show_particle_params).grid(column = 3, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Counter pass', variable = self.sigma, value = -1,
			command = self.show_particle_params).grid(column = 4, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Label(frame, text = "Orbits").grid(column = 1, row = row, sticky = tk.E + tk.N)

		self.charge = tk.IntVar(); self.charge.set(int(data[17])); row += 1; self.row_particle = row
		self.charge_R1 = tk.Radiobutton(frame, text = 'Electrons', variable = self.charge, value = -1)
		self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.charge_R2 = tk.Radiobutton(frame, text = 'Ions', variable = self.charge, value = 1)
		self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.charge_label = tk.Label(frame, text = "Species")
		self.charge_label.grid(column = 1, row = row, sticky = tk.E )

		self.Ekin = tk.StringVar(); self.Ekin.set(str(data[18])); row += 1
		self.Ekin_entry = tk.Entry(frame, width = 7, textvariable = self.Ekin)
		self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.Ekin_label = tk.Label(frame, text = "kin Energy [keV]")
		self.Ekin_label.grid(column = 1, row = row, sticky = tk.E )

		self.Lambda = tk.StringVar(); self.Lambda.set(str(data[19]));
		self.Lambda_entry = tk.Entry(frame, width = 7, textvariable = self.Lambda)
		self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.Lambda_label = tk.Label(frame, text = "Energy ratio")	
		self.Lambda_label.grid(column = 3, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator4 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator4.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- Filament ---
		self.useFilament = tk.StringVar(); self.useFilament.set(str(int(data[20]))); row += 1
		tk.Entry(frame, width = 7, textvariable = self.useFilament).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "# current Filaments").grid(column = 1, row = row, sticky = tk.E )

		# --- run Button ---
		row += 1
		runButton = tk.Button(frame, text = "Run dtfix", command = self.run_funct)
		runButton.grid(column = 1, row = row, columnspan = 5, sticky = tk.W + tk.E)

		# --- adjust style for all ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.show_particle_params()


	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		if(self.Shot.get() == '') | (self.Time.get() == ''):
			print 'You must enter a Shot # and a Time'
			return

		cwd = os.getcwd()
		if(self.path.get()[0] == '.'):
			path = cwd + self.path.get()[1::]
		else: 
			path = self.path.get()
		if not (path[-1] == '/'): path += '/'
		chk = not (cwd + '/' == path)
		
		if chk: os.chdir(path)	
		self.writeControlFile('_fix.dat')
		if(HOST == 'head.cluster'):		# Drop Cluster
			self.write_qsub_file(int(self.HypPt.get()), self.tag.get())
			call('qsub run_job', shell = True)
		else:
			call('dtfix _fix.dat ' + str(int(self.HypPt.get())) + ' ' + self.tag.get() + ' &', shell = True)		
		if chk: os.chdir(cwd)
	
			
	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, period, tag):
		with open('run_job', 'w') as f:
			f.write('#$ -N fix' + tag + '\n')
			f.write('#$ -cwd \n')
			f.write('#$ -o /home/wingen/work/batch.out \n')
			f.write('#$ -e /home/wingen/work/batch.err \n')
			f.write('#$ -S /bin/bash \n')
			f.write('#$ -V \n')
			f.write('#$ -q all.q \n')
			f.write('dtfix _fix.dat ' + str(period) + ' ' + tag + '\n')
			

	# --- turn on/off period entry ---
	def activate_entrys(self):
		if(abs(self.HypRPt.get()) == 1):
			self.period_entry.configure(state=tk.DISABLED)
			if(self.HypRPt.get() == 1): # lower X-point
				self.xmin.set(str(4.1))
				self.xmax.set(str(4.6))
			elif(self.HypRPt.get() == -1): # upper X-point
				self.xmin.set(str(1.8))
				self.xmax.set(str(2.2))
			self.Nx.set(str(30))
			self.ymin.set(str(1))
			self.ymax.set(str(1.3))
			self.Ny.set(str(30))
			self.HypPt.set(str(1))
			self.xmin_entry.configure(state=tk.DISABLED)
			self.xmax_entry.configure(state=tk.DISABLED)
			self.Nx_entry.configure(state=tk.DISABLED)
			self.ymin_entry.configure(state=tk.DISABLED)
			self.ymax_entry.configure(state=tk.DISABLED)
			self.Ny_entry.configure(state=tk.DISABLED)
		else:
			self.period_entry.configure(state=tk.NORMAL)
			self.xmin_entry.configure(state=tk.NORMAL)
			self.xmax_entry.configure(state=tk.NORMAL)
			self.Nx_entry.configure(state=tk.NORMAL)
			self.ymin_entry.configure(state=tk.NORMAL)
			self.ymax_entry.configure(state=tk.NORMAL)
			self.Ny_entry.configure(state=tk.NORMAL)
		
			
	# --- turn on/off response radiobutton ---
	def activate_response(self):
		if(self.useM3DC1.get() == -1):
			self.response_R1.configure(state=tk.DISABLED)
			self.response_R2.configure(state=tk.DISABLED)
		else:
			self.response_R1.configure(state=tk.NORMAL)
			self.response_R2.configure(state=tk.NORMAL)

	
	# --- Show or Hide Particle Options, depending on sigma ---
	def show_particle_params(self):
		if not (self.sigma.get() == 0):
			row = self.row_particle
			self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			row = self.row_particle + 1
			self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Ekin_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Lambda_label.grid(column = 3, row = row, sticky = tk.E , padx=5, pady=5)
		else:
			self.charge_R1.grid_forget()
			self.charge_R2.grid_forget()
			self.charge_label.grid_forget()
			self.Ekin_entry.grid_forget()
			self.Ekin_label.grid_forget()
			self.Lambda_entry.grid_forget()
			self.Lambda_label.grid_forget()
	
	
	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		N = int(self.Nx.get()) * int(self.Ny.get())
		with open(name, 'w') as f:
			f.write('# Parameterfile for DIII-D Programs\n')
			f.write('# Shot: ' + self.Shot.get() + '\tTime: ' + self.Time.get() + 'ms\n')
			f.write('# Path: ' + self.gPath.get() + '\n')
			f.write('shift=\t' + str(self.shift) + '\n')
			f.write('itt=\t0\n')			
			f.write('rmin=\t' + self.ymin.get() + '\n')
			f.write('rmax=\t' + self.ymax.get() + '\n')
			f.write('thmin=\t' + self.xmin.get() + '\n')
			f.write('thmax=\t' + self.xmax.get() + '\n')
			f.write('N=\t' + str(N) + '\n')
			f.write('phistart(deg)=\t' + self.phistart.get() + '\n')
			f.write('MapDirection=\t' + str(self.MapDirection) + '\n')
			f.write('PlasmaResponse(0=no,1=yes)=\t' + str(self.response.get()) + '\n')
			f.write('Field(-1=M3D-C1_off,0=Eq,1=I-coil,2=both)=\t' + str(self.useM3DC1.get()) + '\n')
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')			
			f.write('createPoints(0=setr,3=setpsi,5=setR)=\t0\n')			
			f.write('useFcoil(0=no,1=yes)=\t' + str(self.useFcoil.get()) + '\n')
			f.write('useCcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
			f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
			f.write('ParticleDirection(1=pass,-1=co-pass,0=field-lines)=\t' + str(self.sigma.get()) + '\n')
			f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge.get()) + '\n')
			f.write('Ekin[keV]=\t' + self.Ekin.get() + '\n')
			f.write('lambda=\t' + self.Lambda.get() + '\n')
			f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')
		
		
# -------------------------------------------------------------------------------------------------------------
# --- dtman ---------------------------------------------------------------------------------------------------
class dtman_gui:
	def __init__(self, frame, M):
		
		# --- read parameterfile, if it is there ---
		if os.path.isfile('_fix.dat'):
			_,_,_, data = M.readControlFile('_fix.dat')
		else: # defaults
			data = [1e-4, 0, 1, 1.3, 4.1, 4.6, 900, 0, 1, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 3.141592653589793, 6.283185307179586]
	
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		
		self.xmin = data[4]
		self.xmax = data[5]
		self.ymin = data[2]
		self.ymax = data[3]
		self.N = data[6]
		
		okayCommand = frame.register(self.isOkay)
		notOkayCommand = frame.register(self.isNotOkay)

		# --- type ---
		# define... 
		self.Type = tk.IntVar(); row = 0
		if(data[8] >= 0): self.Type.set(1)
		else: self.Type.set(-1)

		# ...and set type RadioButton
		tk.Radiobutton(frame, text = 'unstable', variable = self.Type, value = 1, 
			command = self.set_type).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'stable', variable = self.Type, value = -1, 
			command = self.set_type).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Manifold").grid(column = 1, row = row, sticky = tk.E )

		# --- shift ---
		self.shift = tk.StringVar(); self.shift.set(str(abs(data[0]))); row += 1
		tk.Entry(frame, width = 7, textvariable = self.shift).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "|Shift|").grid(column = 1, row = row, sticky = tk.E )

		self.shift_sign = tk.IntVar()
		if(data[0] >= 0): self.shift_sign.set(1)
		else: self.shift_sign.set(-1)
		tk.Radiobutton(frame, text = 'right', variable = self.shift_sign, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'left', variable = self.shift_sign, value = -1).grid(column = 4, row = row, sticky = tk.W + tk.E )
				
		# --- info Text ---
		row += 1
		tk.Label(frame, text = "Note: For upper X switch 'left' and 'right'").grid(column = 2, row = row, columnspan = 4, sticky = tk.E + tk.W, ipady = 1)

		# --- dtfix file info ---
		row += 1
		tk.Label(frame, text = "dtfix File").grid(column = 1, row = row, sticky = tk.E )
		
		self.fixfile_tag = tk.StringVar(); self.fixfile_tag.set(self.tag.get())
		tk.Entry(frame, width = 22, textvariable = self.fixfile_tag, validate='all',
        	validatecommand = (okayCommand, '%P'), invalidcommand = notOkayCommand).grid(column = 2, row = row, columnspan = 2, sticky = tk.E)
		tk.Label(frame, text = " Tag").grid(column = 2, row = row, sticky = tk.W )
		
		self.cptag = tk.Button(frame, text = "Copy File Tag", command = self.read_file_tag)
		self.cptag.grid(column = 4, row = row, sticky = tk.W)
		
		self.fixfile_period = tk.StringVar(); self.fixfile_period.set(str(1))
		tk.Entry(frame, width = 4, textvariable = self.fixfile_period).grid(column = 5, row = row, sticky = tk.E)
		tk.Label(frame, text = "Period").grid(column = 5, row = row, sticky = tk.W )
		
		# --- invisible separator ---
		row += 1
		tk.Label(frame, text = "").grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)
		
		# --- invisible separator ---
		row += 1
		tk.Label(frame, text = "").grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)

		# --- phistart ---
		self.phistart = tk.StringVar(); self.phistart.set(str(data[7])); row += 1
		tk.Entry(frame, width = 7, textvariable = self.phistart).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "tor. Angle [deg]").grid(column = 1, row = row, sticky = tk.E )
		tk.Label(frame, text = "For Machine coord. use negative angles").grid(column = 3, row = row, columnspan = 3, sticky = tk.W )

		# --- MapDirection ---
		self.MapDirection = tk.IntVar(); self.MapDirection.set(int(data[8])); row += 1
		tk.Radiobutton(frame, text = '+1', variable = self.MapDirection, value = 1, state = tk.DISABLED).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = '-1', variable = self.MapDirection, value = -1, state = tk.DISABLED).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'both', variable = self.MapDirection, value = 0, state = tk.DISABLED).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Map Direction").grid(column = 1, row = row, sticky = tk.E )

		self.set_type()

		# --- separator ---
		row += 1
		separator = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- coils ---
		self.useFcoil = tk.IntVar(); self.useFcoil.set(int(data[13])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useFcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useFcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use F-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useCcoil = tk.IntVar(); self.useCcoil.set(int(data[14])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useCcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useCcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use C-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useIcoil = tk.IntVar(); self.useIcoil.set(int(data[15])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useIcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useIcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use I-coil").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator2 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator2.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- M3DC1 ---
		self.useM3DC1 = tk.IntVar(); self.useM3DC1.set(int(data[10])); row += 1
		tk.Radiobutton(frame, text = 'g-file Vacuum', variable = self.useM3DC1, value = -1,
			command = self.activate_response).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Equilibrium', variable = self.useM3DC1, value = 0,
			command = self.activate_response).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Perturbation', variable = self.useM3DC1, value = 1,
			command = self.activate_response).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Eq + Pert', variable = self.useM3DC1, value = 2,
			command = self.activate_response).grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use M3DC1").grid(column = 1, row = row, sticky = tk.E )

		self.response = tk.IntVar(); self.response.set(int(data[9])); row += 1
		self.response_R1 = tk.Radiobutton(frame, text = 'Off', variable = self.response, value = 0)
		self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.response_R2 = tk.Radiobutton(frame, text = 'On', variable = self.response, value = 1)
		self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Plasma response").grid(column = 1, row = row, sticky = tk.E )
		
		self.activate_response()

		# --- separator ---
		row += 1
		separator3 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator3.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- particles ---
		self.sigma = tk.IntVar(); self.sigma.set(int(data[16])); row += 1
		tk.Radiobutton(frame, text = 'Field lines', variable = self.sigma, value = 0, 
			command = self.show_particle_params).grid(column = 2, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Co pass', variable = self.sigma, value = 1,
			command = self.show_particle_params).grid(column = 3, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Counter pass', variable = self.sigma, value = -1,
			command = self.show_particle_params).grid(column = 4, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Label(frame, text = "Orbits").grid(column = 1, row = row, sticky = tk.E + tk.N)

		self.charge = tk.IntVar(); self.charge.set(int(data[17])); row += 1; self.row_particle = row
		self.charge_R1 = tk.Radiobutton(frame, text = 'Electrons', variable = self.charge, value = -1)
		self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.charge_R2 = tk.Radiobutton(frame, text = 'Ions', variable = self.charge, value = 1)
		self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.charge_label = tk.Label(frame, text = "Species")
		self.charge_label.grid(column = 1, row = row, sticky = tk.E )

		self.Ekin = tk.StringVar(); self.Ekin.set(str(data[18])); row += 1
		self.Ekin_entry = tk.Entry(frame, width = 7, textvariable = self.Ekin)
		self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.Ekin_label = tk.Label(frame, text = "kin Energy [keV]")
		self.Ekin_label.grid(column = 1, row = row, sticky = tk.E )

		self.Lambda = tk.StringVar(); self.Lambda.set(str(data[19]));
		self.Lambda_entry = tk.Entry(frame, width = 7, textvariable = self.Lambda)
		self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.Lambda_label = tk.Label(frame, text = "Energy ratio")	
		self.Lambda_label.grid(column = 3, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator4 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator4.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- Filament ---
		self.useFilament = tk.StringVar(); self.useFilament.set(str(int(data[20]))); row += 1
		tk.Entry(frame, width = 7, textvariable = self.useFilament).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "# current Filaments").grid(column = 1, row = row, sticky = tk.E )

		# --- run Button ---
		row += 1
		runButton = tk.Button(frame, text = "Run dtman", command = self.run_funct)
		runButton.grid(column = 1, row = row, columnspan = 5, sticky = tk.W + tk.E)

		# --- adjust style for all ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.cptag.grid_configure(pady = 1)
		self.show_particle_params()


	# --- validity check of FixFile Tag ---
	# returns True if tag is alphanumeric or has _ + - as special chars
	# else returns False
	def isOkay(self, tag):
		if(len(tag) == 0): return True
		else: return tag.translate(None, '_+-').isalnum()
		
	# prints error Message, if isOkay == False
	def isNotOkay(self):
		print 'Warning: Invalid Input Character in Tag. Only alphanumeric and + - _ are allowed'


	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		if(self.Shot.get() == '') | (self.Time.get() == ''):
			print 'You must enter a Shot # and a Time'
			return

		cwd = os.getcwd()
		if(self.path.get()[0] == '.'):
			path = cwd + self.path.get()[1::]
		else: 
			path = self.path.get()
		if not (path[-1] == '/'): path += '/'
		chk = not (cwd + '/' == path)
		
		if(self.fixfile_tag.get() == ''): fixtag = ''
		else: fixtag = '_' + self.fixfile_tag.get()
		
		if chk: os.chdir(path)	
		self.writeControlFile('_fix.dat')
		fixfile = 'fix_' + str(int(self.fixfile_period.get())) + fixtag + '.dat'
		if(HOST == 'head.cluster'):		# Drop Cluster
			self.write_qsub_file(fixfile, self.tag.get())
			call('qsub run_job', shell = True)
		else:
			call('dtman _fix.dat ' + fixfile + ' ' + self.tag.get() + ' &', shell = True)		
		if chk: os.chdir(cwd)


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, fixfile, tag):
		with open('run_job', 'w') as f:
			f.write('#$ -N M' + tag + '\n')
			f.write('#$ -cwd \n')
			f.write('#$ -o /home/wingen/work/batch.out \n')
			f.write('#$ -e /home/wingen/work/batch.err \n')
			f.write('#$ -S /bin/bash \n')
			f.write('#$ -V \n')
			f.write('#$ -q all.q \n')
			f.write('dtman _fix.dat ' + fixfile + ' ' +  tag + '\n')
			

	# --- read tag and paste it into dtfix tag ---
	def read_file_tag(self):
		self.fixfile_tag.set(self.tag.get())
		

	# --- Change manifold type ---
	def set_type(self):
		if(self.Type.get() == 1): # unstable
			self.MapDirection.set(1)
			self.shift_sign.set(1)
		else:
			self.MapDirection.set(-1)
			self.shift_sign.set(-1)
			
			
	# --- turn on/off response radiobutton ---
	def activate_response(self):
		if(self.useM3DC1.get() == -1):
			self.response_R1.configure(state=tk.DISABLED)
			self.response_R2.configure(state=tk.DISABLED)
		else:
			self.response_R1.configure(state=tk.NORMAL)
			self.response_R2.configure(state=tk.NORMAL)

	
	# --- Show or Hide Particle Options, depending on sigma ---
	def show_particle_params(self):
		if not (self.sigma.get() == 0):
			row = self.row_particle
			self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			row = self.row_particle + 1
			self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Ekin_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Lambda_label.grid(column = 3, row = row, sticky = tk.E , padx=5, pady=5)
		else:
			self.charge_R1.grid_forget()
			self.charge_R2.grid_forget()
			self.charge_label.grid_forget()
			self.Ekin_entry.grid_forget()
			self.Ekin_label.grid_forget()
			self.Lambda_entry.grid_forget()
			self.Lambda_label.grid_forget()
	
	
	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		shift = abs(float(self.shift.get())) * self.shift_sign.get()
		with open(name, 'w') as f:
			f.write('# Parameterfile for DIII-D Programs\n')
			f.write('# Shot: ' + self.Shot.get() + '\tTime: ' + self.Time.get() + 'ms\n')
			f.write('# Path: ' + self.gPath.get() + '\n')
			f.write('shift=\t' + str(shift) + '\n')
			f.write('itt=\t0\n')		
			f.write('rmin=\t' + str(self.ymin) + '\n')
			f.write('rmax=\t' + str(self.ymax) + '\n')
			f.write('thmin=\t' + str(self.xmin) + '\n')
			f.write('thmax=\t' + str(self.xmax) + '\n')
			f.write('N=\t' + str(self.N) + '\n')
			f.write('phistart(deg)=\t' + self.phistart.get() + '\n')
			f.write('MapDirection=\t' + str(self.MapDirection.get()) + '\n')
			f.write('PlasmaResponse(0=no,1=yes)=\t' + str(self.response.get()) + '\n')
			f.write('Field(-1=M3D-C1_off,0=Eq,1=I-coil,2=both)=\t' + str(self.useM3DC1.get()) + '\n')
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')			
			f.write('createPoints(0=setr,3=setpsi,5=setR)=\t0\n')			
			f.write('useFcoil(0=no,1=yes)=\t' + str(self.useFcoil.get()) + '\n')
			f.write('useCcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
			f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
			f.write('ParticleDirection(1=pass,-1=co-pass,0=field-lines)=\t' + str(self.sigma.get()) + '\n')
			f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge.get()) + '\n')
			f.write('Ekin[keV]=\t' + self.Ekin.get() + '\n')
			f.write('lambda=\t' + self.Lambda.get() + '\n')
			f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')
		

# -------------------------------------------------------------------------------------------------------------
# --- dtfoot --------------------------------------------------------------------------------------------------
class dtfoot_gui:
	def __init__(self, frame, M):
		
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.readControlFile = M.readControlFile

		# --- target  type ---
		# define... 
		self.TargetFlag = tk.IntVar(); self.TargetFlag.set(1); row = 0

		# ...and set grid-type RadioButton
		tk.Radiobutton(frame, text = 'Inner', variable = self.TargetFlag, value = 1, 
			command = self.read_defaults).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Outer', variable = self.TargetFlag, value = 2, 
			command = self.read_defaults).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Shelf', variable = self.TargetFlag, value = 3, 
			command = self.read_defaults).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Target").grid(column = 1, row = row, sticky = tk.E )
		
		# --- x -> phi ---
		row += 1
		self.x_label = tk.Label(frame, text = "phi [rad]")
		self.x_label.grid(column = 1, row = row, sticky = tk.E)
		
		# Min
		self.xmin = tk.StringVar()
		tk.Entry(frame, width = 16, textvariable = self.xmin).grid(column = 2, row = row, columnspan = 2)
		xmin_label = tk.Label(frame, text = "  Min")
		xmin_label.grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.xmax = tk.StringVar()
		tk.Entry(frame, width = 16, textvariable = self.xmax).grid(column = 4, row = row, columnspan = 2)
		xmax_label = tk.Label(frame, text = "Max")
		xmax_label.grid(column = 4, row = row, sticky = tk.W )

		# Nx -> Nphi
		self.Nx = tk.StringVar(); row += 1
		tk.Entry(frame, width = 16, textvariable = self.Nx).grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )
		
		self.pi_text = tk.Text(frame, height= 1, width = 30, bd  = 0, takefocus = 0, bg = frame.cget('bg'), relief = tk.FLAT)
		self.pi_text.grid(column = 4, row = row, columnspan = 2); self.pi_row = row
		self.pi_text.insert(1.0, 'pi = 3.141593 2pi = 6.283185')
		self.pi_text.configure(state = "disabled")
		
		# --- y -> t ---
		row += 1
		self.y_label = tk.Label(frame, text = "")
		self.y_label.grid(column = 1, row = row, sticky = tk.E )
		
		# Min
		self.ymin = tk.StringVar()
		tk.Entry(frame, width = 16, textvariable = self.ymin).grid(column = 2, row = row, columnspan = 2)
		ymin_label = tk.Label(frame, text = "  Min")
		ymin_label.grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.ymax = tk.StringVar()
		tk.Entry(frame, width = 16, textvariable = self.ymax).grid(column = 4, row = row, columnspan = 2)
		ymax_label = tk.Label(frame, text = "Max")
		ymax_label.grid(column = 4, row = row, sticky = tk.W )

		# Ny -> Nt
		self.Ny = tk.StringVar(); row += 1
		tk.Entry(frame, width = 16, textvariable = self.Ny).grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )
		
		# --- Info Text ---
		row += 1
		self.Info = tk.Label(frame, text = "")
		self.Info.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)

		# --- toroidal turns ---
		self.itt = tk.StringVar(); row += 1
		tk.Entry(frame, width = 7, textvariable = self.itt).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "tor. Iterations").grid(column = 1, row = row, sticky = tk.E )

		# --- MapDirection ---
		self.MapDirection = tk.IntVar(); row += 1
		tk.Radiobutton(frame, text = '+1', variable = self.MapDirection, value = 1, state = tk.DISABLED).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = '-1', variable = self.MapDirection, value = -1, state = tk.DISABLED).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'both', variable = self.MapDirection, value = 0, state = tk.DISABLED).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Map Direction").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- coils ---
		self.useFcoil = tk.IntVar(); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useFcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useFcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use F-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useCcoil = tk.IntVar(); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useCcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useCcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use C-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useIcoil = tk.IntVar(); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useIcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useIcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use I-coil").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator2 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator2.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- M3DC1 ---
		self.useM3DC1 = tk.IntVar(); row += 1
		tk.Radiobutton(frame, text = 'g-file Vacuum', variable = self.useM3DC1, value = -1,
			command = self.activate_response).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Equilibrium', variable = self.useM3DC1, value = 0,
			command = self.activate_response).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Perturbation', variable = self.useM3DC1, value = 1,
			command = self.activate_response).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Eq + Pert', variable = self.useM3DC1, value = 2,
			command = self.activate_response).grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use M3DC1").grid(column = 1, row = row, sticky = tk.E )

		self.response = tk.IntVar(); row += 1
		self.response_R1 = tk.Radiobutton(frame, text = 'Off', variable = self.response, value = 0)
		self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.response_R2 = tk.Radiobutton(frame, text = 'On', variable = self.response, value = 1)
		self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Plasma response").grid(column = 1, row = row, sticky = tk.E )
		
		# --- separator ---
		row += 1
		separator3 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator3.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- particles ---
		self.sigma = tk.IntVar(); row += 1
		tk.Radiobutton(frame, text = 'Field lines', variable = self.sigma, value = 0, 
			command = self.show_particle_params).grid(column = 2, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Co pass', variable = self.sigma, value = 1,
			command = self.show_particle_params).grid(column = 3, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Counter pass', variable = self.sigma, value = -1,
			command = self.show_particle_params).grid(column = 4, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Label(frame, text = "Orbits").grid(column = 1, row = row, sticky = tk.E + tk.N)

		self.charge = tk.IntVar(); row += 1; self.row_particle = row
		self.charge_R1 = tk.Radiobutton(frame, text = 'Electrons', variable = self.charge, value = -1)
		self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.charge_R2 = tk.Radiobutton(frame, text = 'Ions', variable = self.charge, value = 1)
		self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.charge_label = tk.Label(frame, text = "Species")
		self.charge_label.grid(column = 1, row = row, sticky = tk.E )

		self.Ekin = tk.StringVar(); row += 1
		self.Ekin_entry = tk.Entry(frame, width = 7, textvariable = self.Ekin)
		self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.Ekin_label = tk.Label(frame, text = "kin Energy [keV]")
		self.Ekin_label.grid(column = 1, row = row, sticky = tk.E )

		self.Lambda = tk.StringVar();
		self.Lambda_entry = tk.Entry(frame, width = 7, textvariable = self.Lambda)
		self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.Lambda_label = tk.Label(frame, text = "Energy ratio")	
		self.Lambda_label.grid(column = 3, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator4 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator4.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- Filament ---
		self.useFilament = tk.StringVar(); row += 1
		tk.Entry(frame, width = 7, textvariable = self.useFilament).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "# current Filaments").grid(column = 1, row = row, sticky = tk.E )

		# --- number of processes for mpi ---
		self.nproc = tk.StringVar(); self.nproc.set(str(4)); row += 1
		self.nproc_entry = tk.Entry(frame, width = 4, textvariable = self.nproc)
		self.nproc_entry.grid(column = 1, row = row, sticky = tk.E)
		tk.Label(frame, text = "       # Procs").grid(column = 1, row = row, sticky = tk.W )
		
		# --- run Button ---
		runButton = tk.Button(frame, text = "Run dtfoot", command = self.run_funct)
		runButton.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)

		# --- adjust style for all ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		
		# --- read parameterfile, if it is there ---
		self.read_defaults()
		self.activate_response()
		self.show_particle_params()


	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		if(self.Shot.get() == '') | (self.Time.get() == ''):
			print 'You must enter a Shot # and a Time'
			return
		
		cwd = os.getcwd()
		if(self.path.get()[0] == '.'):
			path = cwd + self.path.get()[1::]
		else: 
			path = self.path.get()
		if not (path[-1] == '/'): path += '/'
		chk = not (cwd + '/' == path)
		
		if chk: os.chdir(path)	
		if(self.TargetFlag.get() == 1):
			self.writeControlFile('_inner.dat')
			if(HOST == 'head.cluster'):		# Drop Cluster
				self.write_qsub_file(int(self.nproc.get()), self.tag.get(), 'in')
				call('qsub run_mpijob', shell = True)
			else:
				call(MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' dtfoot_mpi _inner.dat ' + self.tag.get() + ' &', shell = True)		
		elif(self.TargetFlag.get() == 2):
			self.writeControlFile('_outer.dat')
			if(HOST == 'head.cluster'):		# Drop Cluster
				self.write_qsub_file(int(self.nproc.get()), self.tag.get(), 'out')
				call('qsub run_mpijob', shell = True)
			else:
				call(MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' dtfoot_mpi _outer.dat ' + self.tag.get() + ' &', shell = True)		
		else:
			self.writeControlFile('_shelf.dat')
			if(HOST == 'head.cluster'):		# Drop Cluster
				self.write_qsub_file(int(self.nproc.get()), self.tag.get(), 'shelf')
				call('qsub run_mpijob', shell = True)
			else:
				call(MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' dtfoot_mpi _shelf.dat ' + self.tag.get() + ' &', shell = True)		
		if chk: os.chdir(cwd)
		
		
	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, type):
		with open('run_mpijob', 'w') as f:
			f.write('#$ -N F' + tag + '\n')
			f.write('#$ -cwd \n')
			f.write('#$ -o /home/wingen/work/batch.out \n')
			f.write('#$ -e /home/wingen/work/batch.err \n')
			f.write('#$ -S /bin/bash \n')
			f.write('#$ -V \n')
			f.write('#$ -q all.q \n')
			f.write('#$ -pe mpi ' + str(nproc) + ' \n')
			f.write('module load openmpi-1.6/gcc \n')
			if(type == 'in'):
				f.write('mpirun -n ${NSLOTS} dtfoot_mpi _inner.dat ' + tag + '\n')
			elif(type == 'out'):
				f.write('mpirun -n ${NSLOTS} dtfoot_mpi _outer.dat ' + tag + '\n')
			elif(type == 'shelf'):
				f.write('mpirun -n ${NSLOTS} dtfoot_mpi _shelf.dat ' + tag + '\n')
			
	
	# --- read initial settings from file or set defaults ---
	def read_defaults(self):
		if(self.TargetFlag.get() == 1):
			self.y_label.configure(text = "t [-1 <--> 1]")
			self.Info.configure(text = "t < 0: Centerpost upwards, t > 0: 45deg Tile downwards, t = 0: Connection point")
			if os.path.isfile('_inner.dat'):
				_,_,_, data = self.readControlFile('_inner.dat')
			else: # defaults
				data = [500, 300, -0.1, 0.3, 0, 6.283185307179586, 400, 0, -1, 0, -1, 1, 2, 1, 1, 1, 0, 1, 
						100, 0.1, 0, 3.141592653589793, 6.283185307179586]
		elif(self.TargetFlag.get() == 2):
			self.y_label.configure(text = "t [0 <--> 1]")
			self.Info.configure(text = "t > 0: Divertor Floor outwards, t = 0: Connection 45deg Tile, t = 1: Pump Entry")
			if os.path.isfile('_outer.dat'):
				_,_,_, data = self.readControlFile('_outer.dat')
			else: # defaults
				data = [500, 300, 0.9, 1.0, 0, 6.283185307179586, 100, 0, 1, 0, -1, 2, 2, 1, 1, 1, 0, 1, 
						100, 0.1, 0, 3.141592653589793, 6.283185307179586]
		else:
			self.y_label.configure(text = "t [0 <--> 1]")
			self.Info.configure(text = "t > 0: Shelf outwards, t = 0: Nose edge")
			if os.path.isfile('_shelf.dat'):
				_,_,_, data = self.readControlFile('_shelf.dat')
			else: # defaults
				data = [500, 300, 0.0, 0.1, 0, 6.283185307179586, 100, 0, 1, 0, -1, 3, 2, 1, 1, 1, 0, 1, 
						100, 0.1, 0, 3.141592653589793, 6.283185307179586]

		self.xmin.set(str(data[4]))
		self.xmax.set(str(data[5]))
		self.Nx.set(str(int(data[0])))
		self.ymin.set(str(data[2]))
		self.ymax.set(str(data[3]))
		self.Ny.set(str(int(data[6])))
		self.itt.set(str(int(data[1])))
		self.MapDirection.set(int(data[8]))
		self.useFcoil.set(int(data[13]))
		self.useCcoil.set(int(data[14]))
		self.useIcoil.set(int(data[15]))
		self.useM3DC1.set(int(data[10]))
		self.response.set(int(data[9]))
		self.sigma.set(int(data[16]))
		self.charge.set(int(data[17]))
		self.Ekin.set(str(data[18]))
		self.Lambda.set(str(data[19]))
		self.useFilament.set(str(int(data[20])))			
		self.data = data

			
	# --- turn on/off response radiobutton ---
	def activate_response(self):
		if(self.useM3DC1.get() == -1):
			self.response_R1.configure(state=tk.DISABLED)
			self.response_R2.configure(state=tk.DISABLED)
		else:
			self.response_R1.configure(state=tk.NORMAL)
			self.response_R2.configure(state=tk.NORMAL)

	
	# --- Show or Hide Particle Options, depending on sigma ---
	def show_particle_params(self):
		if not (self.sigma.get() == 0):
			row = self.row_particle
			self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			row = self.row_particle + 1
			self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Ekin_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Lambda_label.grid(column = 3, row = row, sticky = tk.E , padx=5, pady=5)
		else:
			self.charge_R1.grid_forget()
			self.charge_R2.grid_forget()
			self.charge_label.grid_forget()
			self.Ekin_entry.grid_forget()
			self.Ekin_label.grid_forget()
			self.Lambda_entry.grid_forget()
			self.Lambda_label.grid_forget()
	
	
	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		with open(name, 'w') as f:
			f.write('# Parameterfile for DIII-D Programs\n')
			f.write('# Shot: ' + self.Shot.get() + '\tTime: ' + self.Time.get() + 'ms\n')
			f.write('# Path: ' + self.gPath.get() + '\n')
			f.write('Nphi=\t' + self.Nx.get() + '\n')
			f.write('itt=\t' + self.itt.get() + '\n')	
			f.write('tmin=\t' + self.ymin.get() + '\n')
			f.write('tmax=\t' + self.ymax.get() + '\n')
			f.write('phimin=\t' + self.xmin.get() + '\n')
			f.write('phimax=\t' + self.xmax.get() + '\n')
			f.write('Nt=\t' + self.Ny.get() + '\n')
			f.write('phistart(deg)=\t0\n')
			f.write('MapDirection=\t' + str(self.MapDirection.get()) + '\n')
			f.write('PlasmaResponse(0=no,1=yes)=\t' + str(self.response.get()) + '\n')
			f.write('Field(-1=M3D-C1_off,0=Eq,1=I-coil,2=both)=\t' + str(self.useM3DC1.get()) + '\n')
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t' + str(self.TargetFlag.get()) + '\n')			
			f.write('createPoints(2=target)=\t2\n')	
			f.write('useFcoil(0=no,1=yes)=\t' + str(self.useFcoil.get()) + '\n')
			f.write('useCcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
			f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
			f.write('ParticleDirection(1=pass,-1=co-pass,0=field-lines)=\t' + str(self.sigma.get()) + '\n')
			f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge.get()) + '\n')
			f.write('Ekin[keV]=\t' + self.Ekin.get() + '\n')
			f.write('lambda=\t' + self.Lambda.get() + '\n')
			f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')
		

# -------------------------------------------------------------------------------------------------------------
# --- dtlam ---------------------------------------------------------------------------------------------------
class dtlam_gui:
	def __init__(self, frame, M):
		
		# --- read parameterfile, if it is there ---
		if os.path.isfile('_lam_psi.dat'):
			_,_,_, data = M.readControlFile('_lam_psi.dat')
		else: # defaults
			data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 0, -1, 1, 3, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 3.141592653589793, 6.283185307179586]
		if os.path.isfile('_lam.dat'):
			_,_,_, data = M.readControlFile('_lam.dat')
		else: # defaults
			data = [930, 200, 1.0, 1.45, -1.367, -0.902, 900, 0, 0, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 3.141592653589793, 6.283185307179586]

		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.readControlFile = M.readControlFile

		# --- grid  type ---
		# define... 
		self.createFlag = tk.StringVar(); row = 0
		if(data[12] == 3) | (data[12] == 4): self.createFlag.set('psi')
		else: self.createFlag.set('RZ')

		# ...and set grid-type RadioButton
		tk.Radiobutton(frame, text = 'RZ', variable = self.createFlag, value = 'RZ', 
			command = self.refresh_grid_labels).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'psi_n', variable = self.createFlag, value = 'psi', 
			command = self.refresh_grid_labels).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Coordinate Type").grid(column = 1, row = row, sticky = tk.E )

		# --- x -> theta or R ---
		row += 1
		self.x_label = tk.Label(frame, text = "")
		self.x_label.grid(column = 1, row = row, sticky = tk.E)
		
		# Min
		self.xmin = tk.StringVar(); 
		if(self.createFlag.get() == 'RZ'): self.xmin.set(str(data[2]))
		else: self.xmin.set(str(data[4]))
		tk.Entry(frame, width = 16, textvariable = self.xmin).grid(column = 2, row = row, columnspan = 2)
		xmin_label = tk.Label(frame, text = "  Min")
		xmin_label.grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.xmax = tk.StringVar();
		if(self.createFlag.get() == 'RZ'): self.xmax.set(str(data[3]))
		else: self.xmax.set(str(data[5]))
		tk.Entry(frame, width = 16, textvariable = self.xmax).grid(column = 4, row = row, columnspan = 2)
		xmax_label = tk.Label(frame, text = "Max")
		xmax_label.grid(column = 4, row = row, sticky = tk.W )

		# Nx -> Nth or NR
		self.Nx = tk.StringVar(); row += 1
		if(self.createFlag.get() == 'psi'): self.Nx.set(str(int(data[0])))
		else: self.Nx.set(str(int(data[6])))
		tk.Entry(frame, width = 16, textvariable = self.Nx).grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )
		
		self.pi_text = tk.Text(frame, height= 1, width = 30, bd  = 0, takefocus = 0, bg = frame.cget('bg'), relief = tk.FLAT)
		self.pi_text.grid(column = 4, row = row, columnspan = 2); 
		self.pi_row = row
		self.pi_text.insert(1.0, 'pi = 3.141593 2pi = 6.283185')
		self.pi_text.configure(state = "disabled")
		
		# --- y -> r, psi or Z ---
		row += 1
		self.y_label = tk.Label(frame, text = "")
		self.y_label.grid(column = 1, row = row, sticky = tk.E )
		
		# Min
		self.ymin = tk.StringVar();
		if(self.createFlag.get() == 'RZ'): self.ymin.set(str(data[4]))
		else: self.ymin.set(str(data[2]))
		tk.Entry(frame, width = 16, textvariable = self.ymin).grid(column = 2, row = row, columnspan = 2)
		ymin_label = tk.Label(frame, text = "  Min")
		ymin_label.grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.ymax = tk.StringVar();
		if(self.createFlag.get() == 'RZ'): self.ymax.set(str(data[5]))
		else: self.ymax.set(str(data[3]))
		tk.Entry(frame, width = 16, textvariable = self.ymax).grid(column = 4, row = row, columnspan = 2)
		ymax_label = tk.Label(frame, text = "Max")
		ymax_label.grid(column = 4, row = row, sticky = tk.W )

		# Ny -> Nr, Npsi or NZ
		self.Ny = tk.StringVar(); row += 1
		if(self.createFlag.get() == 'psi'): self.Ny.set(str(int(data[6])))
		else: self.Ny.set(str(int(data[0])))
		tk.Entry(frame, width = 16, textvariable = self.Ny).grid(column = 2, row = row, columnspan = 2)
		tk.Label(frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )

		# --- toroidal turns ---
		self.itt = tk.StringVar(); self.itt.set(str(int(data[1]))); row += 1
		tk.Entry(frame, width = 7, textvariable = self.itt).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "tor. Iterations").grid(column = 1, row = row, sticky = tk.E )

		# --- phistart ---
		self.phistart = tk.StringVar(); self.phistart.set(str(data[7])); row += 1
		tk.Entry(frame, width = 7, textvariable = self.phistart).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "tor. Angle [deg]").grid(column = 1, row = row, sticky = tk.E )
		tk.Label(frame, text = "For Machine coord. use negative angles").grid(column = 3, row = row, columnspan = 3, sticky = tk.W )

		# --- MapDirection ---
		self.MapDirection = tk.IntVar(); self.MapDirection.set(int(data[8])); row += 1
		tk.Radiobutton(frame, text = '+1', variable = self.MapDirection, value = 1).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = '-1', variable = self.MapDirection, value = -1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'both', variable = self.MapDirection, value = 0).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Map Direction").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- coils ---
		self.useFcoil = tk.IntVar(); self.useFcoil.set(int(data[13])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useFcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useFcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use F-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useCcoil = tk.IntVar(); self.useCcoil.set(int(data[14])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useCcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useCcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use C-coil").grid(column = 1, row = row, sticky = tk.E )

		self.useIcoil = tk.IntVar(); self.useIcoil.set(int(data[15])); row += 1
		tk.Radiobutton(frame, text = 'Off', variable = self.useIcoil, value = 0).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'On', variable = self.useIcoil, value = 1).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use I-coil").grid(column = 1, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator2 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator2.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- M3DC1 ---
		self.useM3DC1 = tk.IntVar(); self.useM3DC1.set(int(data[10])); row += 1
		tk.Radiobutton(frame, text = 'g-file Vacuum', variable = self.useM3DC1, value = -1,
			command = self.activate_response).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Equilibrium', variable = self.useM3DC1, value = 0,
			command = self.activate_response).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Perturbation', variable = self.useM3DC1, value = 1,
			command = self.activate_response).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'Eq + Pert', variable = self.useM3DC1, value = 2,
			command = self.activate_response).grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "use M3DC1").grid(column = 1, row = row, sticky = tk.E )

		self.response = tk.IntVar(); self.response.set(int(data[9])); row += 1
		self.response_R1 = tk.Radiobutton(frame, text = 'Off', variable = self.response, value = 0)
		self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.response_R2 = tk.Radiobutton(frame, text = 'On', variable = self.response, value = 1)
		self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Plasma response").grid(column = 1, row = row, sticky = tk.E )
		
		self.activate_response()

		# --- separator ---
		row += 1
		separator3 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator3.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- particles ---
		self.sigma = tk.IntVar(); self.sigma.set(int(data[16])); row += 1
		tk.Radiobutton(frame, text = 'Field lines', variable = self.sigma, value = 0, 
			command = self.show_particle_params).grid(column = 2, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Co pass', variable = self.sigma, value = 1,
			command = self.show_particle_params).grid(column = 3, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(frame, text = 'Counter pass', variable = self.sigma, value = -1,
			command = self.show_particle_params).grid(column = 4, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Label(frame, text = "Orbits").grid(column = 1, row = row, sticky = tk.E + tk.N)

		self.charge = tk.IntVar(); self.charge.set(int(data[17])); row += 1; self.row_particle = row
		self.charge_R1 = tk.Radiobutton(frame, text = 'Electrons', variable = self.charge, value = -1)
		self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.charge_R2 = tk.Radiobutton(frame, text = 'Ions', variable = self.charge, value = 1)
		self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.charge_label = tk.Label(frame, text = "Species")
		self.charge_label.grid(column = 1, row = row, sticky = tk.E )

		self.Ekin = tk.StringVar(); self.Ekin.set(str(data[18])); row += 1
		self.Ekin_entry = tk.Entry(frame, width = 7, textvariable = self.Ekin)
		self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.Ekin_label = tk.Label(frame, text = "kin Energy [keV]")
		self.Ekin_label.grid(column = 1, row = row, sticky = tk.E )

		self.Lambda = tk.StringVar(); self.Lambda.set(str(data[19]));
		self.Lambda_entry = tk.Entry(frame, width = 7, textvariable = self.Lambda)
		self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.Lambda_label = tk.Label(frame, text = "Energy ratio")	
		self.Lambda_label.grid(column = 3, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator4 = tk.Label(frame, text = "---------------------------------------------------------------------------------------")
		separator4.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- Filament ---
		self.useFilament = tk.StringVar(); self.useFilament.set(str(int(data[20]))); row += 1
		tk.Entry(frame, width = 7, textvariable = self.useFilament).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "# current Filaments").grid(column = 1, row = row, sticky = tk.E )

		# --- number of processes for mpi ---
		self.nproc = tk.StringVar(); self.nproc.set(str(4)); row += 1
		self.nproc_entry = tk.Entry(frame, width = 4, textvariable = self.nproc)
		self.nproc_entry.grid(column = 1, row = row, sticky = tk.E)
		tk.Label(frame, text = "       # Procs").grid(column = 1, row = row, sticky = tk.W )
		
		# --- run Button ---
		runButton = tk.Button(frame, text = "Run dtlaminar", command = self.run_funct)
		runButton.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)

		# --- adjust style for all ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.show_particle_params()
		self.refresh_grid_labels()
		

	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		if(self.Shot.get() == '') | (self.Time.get() == ''):
			print 'You must enter a Shot # and a Time'
			return

		cwd = os.getcwd()
		if(self.path.get()[0] == '.'):
			path = cwd + self.path.get()[1::]
		else: 
			path = self.path.get()
		if not (path[-1] == '/'): path += '/'
		chk = not (cwd + '/' == path)
		
		if chk: os.chdir(path)	
		if(self.createFlag.get() == 'psi'):
			self.writeControlFile('_lam_psi.dat')
			if(HOST == 'head.cluster'):		# Drop Cluster
				self.write_qsub_file(int(self.nproc.get()), self.tag.get(), 'psi')
				call('qsub run_mpijob', shell = True)
			else:
				call(MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' dtlaminar_mpi _lam_psi.dat ' + self.tag.get() + ' &', shell = True)
		else:
			self.writeControlFile('_lam.dat')
			if(HOST == 'head.cluster'):		# Drop Cluster
				self.write_qsub_file(int(self.nproc.get()), self.tag.get())
				call('qsub run_mpijob', shell = True)
			else:
				call(MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' dtlaminar_mpi _lam.dat ' + self.tag.get() + ' &', shell = True)
		if chk: os.chdir(cwd)
		
	
	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, type = 'RZ'):
		with open('run_mpijob', 'w') as f:
			f.write('#$ -N L' + tag + '\n')
			f.write('#$ -cwd \n')
			f.write('#$ -o /home/wingen/work/batch.out \n')
			f.write('#$ -e /home/wingen/work/batch.err \n')
			f.write('#$ -S /bin/bash \n')
			f.write('#$ -V \n')
			f.write('#$ -q all.q \n')
			f.write('#$ -pe mpi ' + str(nproc) + ' \n')
			f.write('module load openmpi-1.6/gcc \n')
			if(type == 'psi'):
				f.write('mpirun -n ${NSLOTS} dtlaminar_mpi _lam_psi.dat ' + tag + '\n')
			else:
				f.write('mpirun -n ${NSLOTS} dtlaminar_mpi _lam.dat ' + tag + '\n')
			

	# --- Change Labels on grid variables, depending on createFlag ---
	def refresh_grid_labels(self):
		if(self.createFlag.get() == 'psi'):
			self.x_label.configure(text = "theta [rad]")
			self.y_label.configure(text = "psi_n")
			self.pi_text.grid(column = 4, row = self.pi_row, columnspan = 2)
			if os.path.isfile('_lam_psi.dat'):
				_,_,_, data = self.readControlFile('_lam_psi.dat')
			else: # defaults
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 0, -1, 1, 3, 1, 1, 1, 0, 1, 
						100, 0.1, 0, 3.141592653589793, 6.283185307179586]
			self.xmin.set(str(data[4]))
			self.xmax.set(str(data[5]))
			self.Nx.set(str(int(data[0])))
			self.ymin.set(str(data[2]))
			self.ymax.set(str(data[3]))
			self.Ny.set(str(int(data[6])))

		else:
			self.x_label.configure(text = "R [m]")
			self.y_label.configure(text = "Z [m]")
			self.pi_text.grid_forget()
			if os.path.isfile('_lam.dat'):
				_,_,_, data = self.readControlFile('_lam.dat')
			else: # defaults
				data = [930, 200, 1.0, 1.45, -1.367, -0.902, 900, 0, 0, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
						100, 0.1, 0, 3.141592653589793, 6.283185307179586]
			self.xmin.set(str(data[2]))
			self.xmax.set(str(data[3]))
			self.Nx.set(str(int(data[6])))
			self.ymin.set(str(data[4]))
			self.ymax.set(str(data[5]))
			self.Ny.set(str(int(data[0])))
			
			
	# --- turn on/off response radiobutton ---
	def activate_response(self):
		if(self.useM3DC1.get() == -1):
			self.response_R1.configure(state=tk.DISABLED)
			self.response_R2.configure(state=tk.DISABLED)
		else:
			self.response_R1.configure(state=tk.NORMAL)
			self.response_R2.configure(state=tk.NORMAL)

	
	# --- Show or Hide Particle Options, depending on sigma ---
	def show_particle_params(self):
		if not (self.sigma.get() == 0):
			row = self.row_particle
			self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			row = self.row_particle + 1
			self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Ekin_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Lambda_label.grid(column = 3, row = row, sticky = tk.E , padx=5, pady=5)
		else:
			self.charge_R1.grid_forget()
			self.charge_R2.grid_forget()
			self.charge_label.grid_forget()
			self.Ekin_entry.grid_forget()
			self.Ekin_label.grid_forget()
			self.Lambda_entry.grid_forget()
			self.Lambda_label.grid_forget()
	
	
	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		with open(name, 'w') as f:
			f.write('# Parameterfile for DIII-D Programs\n')
			f.write('# Shot: ' + self.Shot.get() + '\tTime: ' + self.Time.get() + 'ms\n')
			f.write('# Path: ' + self.gPath.get() + '\n')
			
			if(self.createFlag.get() == 'psi'): f.write('Nth=\t' + self.Nx.get() + '\n')
			else: f.write('NZ=\t' + self.Ny.get() + '\n')
			
			f.write('itt=\t' + self.itt.get() + '\n')
			
			if(self.createFlag.get() == 'psi'):
				f.write('psimin=\t' + self.ymin.get() + '\n')
				f.write('psimax=\t' + self.ymax.get() + '\n')
				f.write('thmin=\t' + self.xmin.get() + '\n')
				f.write('thmax=\t' + self.xmax.get() + '\n')
				f.write('Npsi=\t' + self.Ny.get() + '\n')
			else:
				f.write('Rmin=\t' + self.xmin.get() + '\n')
				f.write('Rmax=\t' + self.xmax.get() + '\n')
				f.write('Zmin=\t' + self.ymin.get() + '\n')
				f.write('Zmax=\t' + self.ymax.get() + '\n')
				f.write('NR=\t' + self.Nx.get() + '\n')

			f.write('phistart(deg)=\t' + self.phistart.get() + '\n')
			f.write('MapDirection=\t' + str(self.MapDirection.get()) + '\n')
			f.write('PlasmaResponse(0=no,1=yes)=\t' + str(self.response.get()) + '\n')
			f.write('Field(-1=M3D-C1_off,0=Eq,1=I-coil,2=both)=\t' + str(self.useM3DC1.get()) + '\n')
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')
			
			if(self.createFlag.get() == 'psi'): f.write('createPoints(0=setR,3=setpsi)=\t3\n')
			else: f.write('createPoints(0=setR,3=setpsi)=\t0\n')
			
			f.write('useFcoil(0=no,1=yes)=\t' + str(self.useFcoil.get()) + '\n')
			f.write('useCcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
			f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
			f.write('ParticleDirection(1=pass,-1=co-pass,0=field-lines)=\t' + str(self.sigma.get()) + '\n')
			f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge.get()) + '\n')
			f.write('Ekin[keV]=\t' + self.Ekin.get() + '\n')
			f.write('lambda=\t' + self.Lambda.get() + '\n')
			f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')
		
		
# -------------------------------------------------------------------------------------------------------------
# --- info tab ------------------------------------------------------------------------------------------------
class info_gui:
	def __init__(self, frame):
	
		row = 0
		self.info_text = tk.Text(frame, height= 20, width = 80, bd  = 0, takefocus = 0, 
								 bg = frame.cget('bg'), relief = tk.FLAT, wrap=tk.WORD)
		self.info_text.grid(column = 1, row = row, columnspan = 5, padx=10, pady=10); 
		self.info_text.insert(1.0, 
		'MAFOT Control GUI for DIII-D, ITER, NSTX & MAST \n\n'
		'MAFOT Version 3.1 \n'
		'GUI Version 1.0 \n'
		'Author: Andreas Wingen \n\n'
		'The GUI creates/reads/modifies the respective MAFOT control files in the working '
		'directory and launches the respective MAFOT tool binary. \n'
		'You can at any time skip the GUI, edit the control file by hand and launch the code manually. \n'
		'MAFOT output is saved as plain ascii files. The output data file includes a '
		'header (comment character is #), which shows all the parameters used in the run and defines the data columns. '
		'Output also includes a log file which can savely be deleted after a successfull run. \n'
		'If errors occur during the run, the log file might include information about it. \n')
		self.info_text.configure(state = "disabled")
		self.info_text.tag_add('line1', 1.0, 2.0)
		self.info_text.tag_config('line1', justify = tk.CENTER, font = ('Times', '16', 'bold'))
		self.info_text.tag_add('versions', 3.0, 6.0)
		self.info_text.tag_config('versions', justify = tk.CENTER, font = ('Times', '12'))

	
# -------------------------------------------------------------------------------------------------------------
# --- main ----------------------------------------------------------------------------------------------------
def main():
	# --- set proper shell environment ---
	global MPIRUN
	if(HOST == 'r2d2.gat.com'):
		LD_LIBRARY_PATH = os.environ['LD_LIBRARY_PATH']
		LD_LIBRARY_PATH += ':/home/wingen/lib/64/blitz/lib:/home/wingen/lib/64:/usr/local/openmpi/lib:/usr/local/m3dc1/lib:/usr/local/hdf5/lib'
		os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
		MPIRUN = '/usr/local/openmpi/bin/mpirun'
	elif(HOST == 'head.cluster'):
		LD_LIBRARY_PATH = os.environ['LD_LIBRARY_PATH']
		LD_LIBRARY_PATH += ':/opt/intel/composerxe-2011.0.084/compiler/lib/intel64:/opt/intel/composerxe-2011.0.084/mpirt/lib/intel64'
		LD_LIBRARY_PATH += ':/opt/intel/composerxe-2011.0.084/mkl/lib/intel64:/home/wingen/lib/64/blitz/lib:/home/wingen/lib/64:/home/wingen/lib/64/hdf5/lib'
		os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH

	# --- set main window ---
	root = tk.Tk()
	root.title("MAFOT Control")
	
	# --- set main frame ---
	mainframe = tk.Frame(root)
	mainframe.grid(column=0, row=0, sticky = tk.N + tk.W + tk.E + tk.S)
	mainframe.columnconfigure(0, weight=1)
	mainframe.rowconfigure(0, weight=1)

	# --- make gui ---
	Common_gui(mainframe)
	
	# --- run gui ---
	root.mainloop()

if __name__ == '__main__':
	main()



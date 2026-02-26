#!/usr/bin/env python
import sys, os, re
if sys.version_info[0] < 3:
	import Tkinter as tk
	import ttk
else:
	import tkinter as tk
	import tkinter.ttk as ttk
import socket
import numpy as np
from subprocess import call

#from Misc.autocompleteEntry import AutocompleteEntry 

HOME = os.getenv('HOME')
HOST = socket.gethostname()
MPIRUN = 'mpirun'

# -------------------------------------------------------------------------------------------------------------
# --- common input --------------------------------------------------------------------------------------------
class Common_gui:
	def __init__(self, frame):
		self.frame = frame
		self.MachFlag = tk.StringVar(); self.MachFlag.set('dt')
		self.Shot = tk.StringVar()
		self.Time = tk.StringVar()
		self.gPath = tk.StringVar()
		self.tag = tk.StringVar()
		self.path = tk.StringVar()
		self.previousPath = None
		self.previousgPath = None
		self.translateMachFlag = {'dt':'DIII-D', 'iter':'ITER', 'nstx':'NSTX', 
									'mast':'MAST', 'cmod':'C-Mod', 'any':'Any'}
		self.labelFont = None
		self.textFont = None
	
		okayCommand = self.frame.register(self.isOkay)
		notOkayCommand = self.frame.register(self.isNotOkay)
		makecwd = self.frame.register(self.makeCWD)
		makeshottime = self.frame.register(self.makeShotTime)
		
		# --- Machine ---
		row = 0
		tk.Radiobutton(self.frame, text = 'DIII-D', variable = self.MachFlag, value = 'dt', 
			command = self.set_Machine).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'ITER', variable = self.MachFlag, value = 'iter', 
			command = self.set_Machine).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'NSTX', variable = self.MachFlag, value = 'nstx', 
			command = self.set_Machine).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'MAST', variable = self.MachFlag, value = 'mast', 
			command = self.set_Machine).grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'ANY', variable = self.MachFlag, value = 'any', 
			command = self.set_Machine).grid(column = 6, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "Machine").grid(column = 1, row = row, sticky = tk.E )

		# --- gPath ---
		row += 1
		#tk.Entry(self.frame, width = 7, textvariable = self.gPath).grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)
		self.gFile_entry = AutocompleteEntry(os.listdir('.'), self.frame, listboxLength = 6, 
				width = 50, textvariable = self.gPath, font = self.textFont, validate = 'focusout', validatecommand = makeshottime)
		self.gFile_entry.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)
		tk.Label(self.frame, text = "g-file Pathname").grid(column = 1, row = row, sticky = tk.E)

		# --- Shot number ---
		row += 1
		tk.Entry(self.frame, width = 7, textvariable = self.Shot).grid(column = 2, row = row, sticky = tk.W + tk.E)
		tk.Label(self.frame, text = "Shot").grid(column = 1, row = row, sticky = tk.E)

		# --- Time ---
		tk.Entry(self.frame, width = 7, textvariable = self.Time).grid(column = 4, row = row, sticky = tk.W + tk.E)
		tk.Label(self.frame, text = "Time").grid(column = 3, row = row, sticky = tk.E)

		# --- Tag ---
		row += 1; self.row_write = row
		self.tag_entry = tk.Entry(self.frame, width = 7, textvariable = self.tag, validate = 'all',
        	validatecommand = (okayCommand, '%P'), invalidcommand = notOkayCommand)
		self.tag_entry.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)
		self.tag_label = tk.Label(self.frame, text = "File Tag")
		self.tag_label.grid(column = 1, row = row, sticky = tk.E)
		
		# --- cwd ---
		row += 1
		self.path_entry = tk.Entry(self.frame, width = 7, textvariable = self.path, validate = 'focusout',
        	validatecommand = makecwd)
		self.path_entry.grid(column = 2, row = row, columnspan = 3, sticky = tk.W + tk.E)
		self.path_label = tk.Label(self.frame, text = "Working Dir")
		self.path_label.grid(column = 1, row = row, sticky = tk.E)
		self.path_entry.bind('<Return>', self.makeCWD)
		
		self.reload_default = tk.IntVar(); self.reload_default.set(0)
		tk.Checkbutton(self.frame, text = 'reload\nsettings', variable = self.reload_default).grid(column = 5, row = row, sticky = tk.W)
		
		# --- invisible separator ---
		row += 1
		separator = tk.Label(self.frame, text = "")
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W)

		# --- adjust style for all ---
		for child in self.frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		
		self.tag_entry.focus()
		
		# --- default setup is for D3D ---
		self.path.set('./')
		self.set_defaults()
		self.set_tabs()
		self.set_Machine()
		self.tag.set('')
		self.reload_default.set(1) # this has to be set last at startup, since it affects the makeCWD callback
	
	
	# --- set default values ---
	def set_defaults(self):
		# set path if it already exists:
		try:
			path = self.path.get()
		except:
			path = './'
			
		# read parameterfile, if it is there
		default_gpath = path
		self.ControlFileFound = True
		if os.path.isfile(path + '_plot.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_plot.dat',default_gpath)
		elif os.path.isfile(path + '_fix.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_fix.dat',default_gpath)
		elif os.path.isfile(path + '_inner.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_inner.dat',default_gpath)
		elif os.path.isfile(path + '_outer.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_outer.dat',default_gpath)
		elif os.path.isfile(path + '_lam.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_lam.dat',default_gpath)
		elif os.path.isfile(path + '_lam_psi.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_lam_psi.dat',default_gpath)
		elif os.path.isfile(path + '_shelf.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_shelf.dat',default_gpath)
		elif os.path.isfile(path + '_innerup.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_innerup.dat',default_gpath)
		elif os.path.isfile(path + '_innerdwn.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_innerdwn.dat',default_gpath)
		elif os.path.isfile(path + '_outerup.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_outerup.dat',default_gpath)
		elif os.path.isfile(path + '_outerdwn.dat'):
			shot, time, gpath, _,_ = readControlFile(path + '_outerdwn.dat',default_gpath)
		else:
			shot = ''
			time = ''
			gpath = path
			self.ControlFileFound = False
			
		#if not (gpath[-1] == '/'): gpath += '/'
		#self.Shot.set(str(shot))
		#self.Time.set(str(time))
		#self.gPath.set(gpath)
		self.gFile_entry.update_baseSearchPath(gpath)
		self.makeShotTime(defShot = shot, defTime = time)


	
	# --- validity check of Tag ---
	# returns True if tag is alphanumeric or has _ + - as special chars
	# else returns False
	def isOkay(self, tag):
		tag = tag.translate(dict.fromkeys('_+- '))
		if(len(tag) == 0): return True
		else: return tag.isalnum()
		
	# prints error Message, if isOkay == False
	def isNotOkay(self):
		print ('Warning: Invalid Input Character in File Tag. Only alphanumeric and + - _ are allowed')
		
	
	# --- callback, executed whenever Working Dir looses focus --- 
	def makeCWD(self, event = None):
		if self.path.get() == self.previousPath:
			self.path_entry.configure(validate = 'focusout')
			return True			
		if len(self.path.get()) < 1: self.path.set('./')
		if not (self.path.get()[-1] == '/'): self.path.set(self.path.get() + '/')
		self.previousPath = self.path.get()
		if(self.reload_default.get() == 1):
			self.set_defaults()
			self.plot_tab.set_defaults()
			self.fix_tab.set_defaults()
			self.man_tab.set_defaults()
			self.foot_tab.set_defaults()
			self.lam_tab.set_defaults()
		
		self.path_entry.configure(validate = 'focusout')
		return True
		
		
	# --- callback, executed whenever gPath looses focus --- 
	def makeShotTime(self, event = None, defShot = None, defTime = None):
		if self.gPath.get() == self.previousgPath:
			self.gFile_entry.configure(validate = 'focusout')
			return True			

		gfileNam = self.gPath.get()
		idx = gfileNam[::-1].find('/')
		if(idx == -1):
			gpath = '.'
			gfile = gfileNam
		else:
			idx *= -1
			gpath = gfileNam[0:idx - 1]  # path without a final '/'
			gfile = gfileNam[idx::]
		
		try:
			idx = gfile.find('.')
			fmtstr = '0' + str(idx - 1) + 'd'
			shot, time = int(gfile[1:idx]), gfile[idx + 1::]
			if '.' in time:
				idx = time.find('.')
				time = time[0:idx]
			if '_' in time:
				idx = time.find('_')
				time = time[0:idx]
			time = int(time)
			self.Shot.set(str(shot))
			self.Time.set(str(time))
		except:
			if defShot is not None: self.Shot.set(str(defShot))
			if defTime is not None: self.Time.set(str(defTime))

		self.previousgPath = self.gPath.get()		
		self.gFile_entry.configure(validate = 'focusout')
		return True


	# --- Tabs ---
	def set_tabs(self):
		# --- start tabs ---
		self.nb = ttk.Notebook(self.frame)
		self.nb.grid(column = 1, row = 7, columnspan = 10)
		
		# --- Tab 0: set plot tab ---
		self.tabframe1 = tk.Frame(self.frame)
		self.tabframe1.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe1.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe1, text = 'dtplot')
		self.plot_tab = set_plot_tab(self.tabframe1, self)

		# --- Tab 1: set fix tab ---
		self.tabframe2 = tk.Frame(self.frame)
		self.tabframe2.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe2.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe2, text = 'dtfix')
		self.fix_tab = set_fix_tab(self.tabframe2, self)
	
		# --- Tab 2: set man tab ---
		self.tabframe3 = tk.Frame(self.frame)
		self.tabframe3.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe3.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe3, text = 'dtman')
		self.man_tab = set_man_tab(self.tabframe3, self)
	
		# --- Tab 3: set foot tab ---
		self.tabframe4 = tk.Frame(self.frame)
		self.tabframe4.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe4.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe4, text = 'dtfoot')
		self.foot_tab = set_foot_tab(self.tabframe4, self)

		# --- Tab 4: set laminar tab ---
		self.tabframe5 = tk.Frame(self.frame)
		self.tabframe5.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe5.rowconfigure(16, weight=1)
		self.nb.add(self.tabframe5, text = 'dtlaminar')	
		self.lam_tab = set_lam_tab(self.tabframe5, self)
		
		# --- Tab 5: set info tab ---
		self.tabframe6 = tk.Frame(self.frame)
		self.tabframe6.grid(column = 1, row = 7, sticky = tk.N + tk.W + tk.E + tk.S)
		self.tabframe6.rowconfigure(24, weight=1)
		self.nb.add(self.tabframe6, text = 'Info')
		info_gui(self.tabframe6)


	# --- set different Machines ---
	def set_machine_elements(self):
		self.plot_tab.set_machine_elements(self.MachFlag.get())
		self.plot_tab.set_defaults()
		self.fix_tab.set_machine_elements(self.MachFlag.get())
		self.fix_tab.set_defaults()
		self.man_tab.set_machine_elements(self.MachFlag.get())
		self.foot_tab.set_machine_elements(self.MachFlag.get())
		self.foot_tab.set_target_elements()
		self.foot_tab.set_defaults()
		self.lam_tab.set_machine_elements(self.MachFlag.get())
		self.lam_tab.set_defaults()
	
	
	# --- set Tabs for different Machines ---
	def set_Machine(self):
		if(self.MachFlag.get() == 'any'):
			#if not self.ControlFileFound:
			#	self.gPath.set(HOME + '/work/gfiles/')
			self.nb.tab(0, text='anyplot')
			self.nb.tab(1, text='anyfix')
			self.nb.tab(2, text='anyman')
			self.nb.tab(3, text='anyfoot')
			self.nb.tab(4, text='anylaminar')
		elif(self.MachFlag.get() == 'iter'):
			#if not self.ControlFileFound:
			#	self.gPath.set(HOME + '/c++/iter/gfiles/')
			self.nb.tab(0, text='iterplot')
			self.nb.tab(1, text='iterfix')
			self.nb.tab(2, text='iterman')
			self.nb.tab(3, text='iterfoot')
			self.nb.tab(4, text='iterlaminar')
		elif(self.MachFlag.get() == 'nstx'):
			#if not self.ControlFileFound:
			#	self.gPath.set(HOME + '/c++/nstx/gfiles/')
			self.nb.tab(0, text='nstxplot')
			self.nb.tab(1, text='nstxfix')
			self.nb.tab(2, text='nstxman')
			self.nb.tab(3, text='nstxfoot')
			self.nb.tab(4, text='nstxlaminar')
		elif(self.MachFlag.get() == 'mast'):
			#if not self.ControlFileFound:
			#	self.gPath.set(HOME + '/c++/mast/gfiles/')
			self.nb.tab(0, text='mastplot')
			self.nb.tab(1, text='mastfix')
			self.nb.tab(2, text='mastman')
			self.nb.tab(3, text='mastfoot')
			self.nb.tab(4, text='mastlaminar')
		else:	# d3d
			#if not self.ControlFileFound:
			#	self.gPath.set(HOME + '/c++/d3d/gfiles/')
			self.nb.tab(0, text='dtplot')
			self.nb.tab(1, text='dtfix')
			self.nb.tab(2, text='dtman')
			self.nb.tab(3, text='dtfoot')
			self.nb.tab(4, text='dtlaminar')
		self.set_machine_elements()


# -------------------------------------------------------------------------------------------------------------
# --- machine specific elements -------------------------------------------------------------------------------
class set_machine():
	def __init__(self, frame, M):
		self.frame = frame
		self.MachFlag = M.MachFlag
		self.gPath = M.gPath
		self.path = M.path
		self.translateMachFlag = M.translateMachFlag
		self.labelFont = None
		self.textFont = None


	def define_machine_elements(self):
		self.useIcoil = tk.IntVar()
		self.useCcoil = tk.IntVar()
		self.useFcoil = tk.IntVar()
		self.useBus = tk.IntVar()
		self.useBcoil = tk.IntVar()

		# --- 3D Fields Label ---
		self.machine_3D_Label = tk.Label(self.frame, font = self.labelFont)
		# --- Icoil checkbox ---
		self.machine_3D_chk1 = tk.Checkbutton(self.frame, variable = self.useIcoil, font = self.labelFont)
		# --- Ccoil checkbox ---
		self.machine_3D_chk2 = tk.Checkbutton(self.frame, variable = self.useCcoil, font = self.labelFont)
		# --- Error Fields Label ---
		self.machine_ef_Label = tk.Label(self.frame, font = self.labelFont)
		# --- useEF checkbox ---
		self.machine_ef_chk1 = tk.Checkbutton(self.frame, variable = self.useFcoil, font = self.labelFont)
		# --- useBus checkbox ---
		self.machine_ef_chk2 = tk.Checkbutton(self.frame, variable = self.useBus, font = self.labelFont)
		# --- Bcoil checkbox ---
		self.machine_ef_chk3 = tk.Checkbutton(self.frame, variable = self.useBcoil, font = self.labelFont)


	def set_machine_elements(self, MachFlag, row = 10):
		try: self.forget_elements()
		except: pass
		if(MachFlag == 'any'):
			self.set_any_elements(row)
		elif(MachFlag == 'iter'):
			self.set_iter_elements(row)
		elif(MachFlag == 'nstx'):
			self.set_nstx_elements(row)
		elif(MachFlag == 'mast'):
			self.set_mast_elements(row)
		else:	# d3d
			self.set_d3d_elements(row)


	def forget_elements(self):
		self.machine_3D_Label.grid_forget()
		self.machine_3D_chk1.grid_forget()
		self.machine_3D_chk2.grid_forget()
		self.machine_ef_Label.grid_forget()
		self.machine_ef_chk1.grid_forget()
		self.machine_ef_chk2.grid_forget()
		self.machine_ef_chk3.grid_forget()


	def set_machine_defaults(self, data, MachFlag, flag = None):
		if data is None: data = self.tool_defaults(flag)
		
		data_dic = {'0-8':data[0:9],'itt':str(int(data[1])),'phistart':repr(data[7]),'MapDirection':int(data[8]),
				'useM3DC1':int(data[10]),'response':int(data[9]),'selectField':int(data[10]),
				'target':int(data[11]),'create':int(data[12]),
				'sigma':int(data[16]),'charge':int(data[17]),'Ekin':repr(data[18]),'Lambda':repr(data[19]),'Mass':repr(data[20])}
		
		if(MachFlag == 'any'):
			data_dic = self.set_any_defaults(data, data_dic)
		elif(MachFlag == 'iter'):
			data_dic = self.set_iter_defaults(data, data_dic)
		elif(MachFlag == 'nstx'):
			data_dic = self.set_nstx_defaults(data, data_dic)
		elif(MachFlag == 'mast'):
			data_dic = self.set_mast_defaults(data, data_dic)
		else:	# d3d
			data_dic = self.set_d3d_defaults(data, data_dic)
		
		self.set_common_defaults(data_dic)
		return data_dic


	def write_coils(self, MachFlag, f):
		if(MachFlag == 'any'):
			self.write_any_coils(f)
		elif(MachFlag == 'iter'):
			self.write_iter_coils(f)
		elif(MachFlag == 'nstx'):
			self.write_nstx_coils(f)
		elif(MachFlag == 'mast'):
			self.write_mast_coils(f)
		else:	# d3d
			self.write_d3d_coils(f)


	# --- Overloaded Functions -----------------------------------------------------------
	def tool_defaults(self, flag):	# gets overloaded by each tool
		return np.zeros(100)

	def set_common_defaults(self, data):	# gets overloaded by common tab
		return	
			
			
	# --- DIII-D -------------------------------------------------------------------------
	def set_d3d_elements(self, row):
		self.machine_3D_Label.configure(text = "3-D Fields: ")
		self.machine_3D_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		self.machine_3D_chk1.configure(text = ' I-coils')
		self.machine_3D_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)
		self.machine_3D_chk2.configure(text = ' C-coils')
		self.machine_3D_chk2.grid(column = 3, row = row, sticky = tk.E + tk.W, padx=5, pady=5)

		row += 1
		self.machine_ef_Label.configure(text = "Error Fields: ")
		self.machine_ef_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		self.machine_ef_chk1.configure(text = ' F-coil')
		self.machine_ef_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)
		self.machine_ef_chk2.configure(text = ' Bus')
		self.machine_ef_chk2.grid(column = 3, row = row, sticky = tk.E + tk.W, padx=5, pady=5)
		self.machine_ef_chk3.configure(text = ' B-coil')
		self.machine_ef_chk3.grid(column = 4, row = row, sticky = tk.E + tk.W, padx=5, pady=5)

	
	def set_d3d_defaults(self, data, data_dic):
		if len(data) < 26:
			if len(data) == 23: data = np.append(data[0:21],[0, 0, 0, 3.141592653589793, 6.283185307179586])
			elif len(data) == 24: data = np.append(data[0:22],[0, 0, 3.141592653589793, 6.283185307179586])
			elif len(data) == 25: data = np.append(data[0:23],[0, 3.141592653589793, 6.283185307179586])

		dic = {'useFcoil':int(data[13]),'useCcoil':int(data[14]),'useIcoil':int(data[15]),
			'useBus':int(data[22]),'useBcoil':int(data[23]),
			'useFilament':str(int(data[21]))}
		for key in dic: data_dic[key] = dic[key]
		
		self.useFcoil.set(data_dic['useFcoil'])
		self.useCcoil.set(data_dic['useCcoil'])
		self.useIcoil.set(data_dic['useIcoil'])
		self.useBus.set(data_dic['useBus'])
		self.useBcoil.set(data_dic['useBcoil'])
		return data_dic
		
	
	def write_d3d_coils(self,f):
		f.write('useFcoil(0=no,1=yes)=\t' + str(self.useFcoil.get()) + '\n')
		f.write('useCcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
		f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
	
	
	def write_d3d_errorFileds(self,f):
		f.write('useBusError(0=no,1=yes)=\t' + str(self.useBus.get()) + '\n')
		f.write('useBcoilError(0=no,1=yes)=\t' + str(self.useBcoil.get()) + '\n')
		

	# --- ANY ---------------------------------------------------------------------------
	def set_any_elements(self, row):
		self.machine_3D_Label.configure(text = "3-D Fields: ")
		self.machine_3D_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		#self.machine_3D_chk1.configure(text = ' ITER-coils')
		#self.machine_3D_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)


	def set_any_defaults(self, data, data_dic):
		data_dic['useIcoil'] = 0
		data_dic['useFilament'] = str(int(data[21]))
		self.useIcoil.set(data_dic['useIcoil'])
		return data_dic


	def write_any_coils(self,f):
		f.write('useFcoil(0=no,1=yes)=\t0\n')
		f.write('useCcoil(0=no,1=yes)=\t0\n')
		f.write('useIcoil(0=no,1=yes)=\t0\n')


	# --- ITER ---------------------------------------------------------------------------
	def set_iter_elements(self, row):
		self.machine_3D_Label.configure(text = "3-D Fields: ")
		self.machine_3D_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		self.machine_3D_chk1.configure(text = ' ITER-coils')
		self.machine_3D_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)


	def set_iter_defaults(self, data, data_dic):
		data_dic['useIcoil'] = int(data[13])
		data_dic['useFilament'] = str(int(data[14]))
		self.useIcoil.set(data_dic['useIcoil'])
		return data_dic


	def write_iter_coils(self,f):
		f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
	

	# --- NSTX ---------------------------------------------------------------------------
	def set_nstx_elements(self, row):
		self.machine_3D_Label.configure(text = "3-D Fields: ")
		self.machine_3D_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		self.machine_3D_chk1.configure(text = ' EC-coils')
		self.machine_3D_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)


	def set_nstx_defaults(self, data, data_dic):
		data_dic['useIcoil'] = int(data[13])
		data_dic['useFilament'] = str(int(data[14]))		
		self.useIcoil.set(data_dic['useIcoil'])
		return data_dic


	def write_nstx_coils(self,f):
		f.write('useECcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
	

	# --- MAST ---------------------------------------------------------------------------
	def set_mast_elements(self, row):
		self.machine_3D_Label.configure(text = "3-D Fields: ")
		self.machine_3D_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		self.machine_3D_chk1.configure(text = ' I-coils')
		self.machine_3D_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)
		self.machine_3D_chk2.configure(text = ' EC-coils')
		self.machine_3D_chk2.grid(column = 3, row = row, sticky = tk.E + tk.W, padx=5, pady=5)


	def set_mast_defaults(self, data, data_dic):
		data_dic['useCcoil'] = int(data[13])
		data_dic['useIcoil'] = int(data[14])
		data_dic['useFilament'] = str(int(data[15]))
		self.useCcoil.set(data_dic['useCcoil'])
		self.useIcoil.set(data_dic['useIcoil'])
		return data_dic


	def write_mast_coils(self,f):
		f.write('useECcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
		f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')


	# --- CMOD ---------------------------------------------------------------------------
	def write_cmod_headerLines(self,f):
		f.write('# Parameterfile for ' + self.translateMachFlag['cmod'] + ' Programs\n')
		f.write('# Shot: ' + format(int(self.Shot.get()),'010d') + '\tTime: ' + format(int(self.Time.get()),'04d') + 'ms\n')
		f.write('# Path: ' + self.gPath.get() + '\n')



# -------------------------------------------------------------------------------------------------------------
# --- common tab elements -------------------------------------------------------------------------------------
class common_tab(set_machine):	# inherit set_machine class
	def __init__(self, frame, M):
		self.frame = frame
		self.MachFlag = M.MachFlag
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.translateMachFlag = M.translateMachFlag
		self.gFile_entry = M.gFile_entry
		self.labelFont = None
		self.textFont = None


	def set_MinMax_elements(self, row):
		self.xmin = tk.StringVar()
		self.xmax = tk.StringVar()
		self.Nx = tk.StringVar()
		self.ymin = tk.StringVar();
		self.ymax = tk.StringVar();
		self.Ny = tk.StringVar(); 
		
		# Min
		self.xmin_entry = tk.Entry(self.frame, width = 17, textvariable = self.xmin, font = self.textFont)
		self.xmin_entry.grid(column = 2, row = row, columnspan = 2)
		tk.Label(self.frame, text = "  Min", font = self.labelFont).grid(column = 2, row = row, sticky = tk.W )
		
		# Max
		self.xmax_entry = tk.Entry(self.frame, width = 17, textvariable = self.xmax, font = self.textFont)
		self.xmax_entry.grid(column = 4, row = row, columnspan = 2)
		tk.Label(self.frame, text = "Max", font = self.labelFont).grid(column = 4, row = row, sticky = tk.W )

		# Nx -> Nth
		row += 1
		self.Nx_entry = tk.Entry(self.frame, width = 17, textvariable = self.Nx, font = self.textFont)
		self.Nx_entry.grid(column = 2, row = row, columnspan = 2)
		tk.Label(self.frame, text = "     #", font = self.labelFont).grid(column = 2, row = row, sticky = tk.W )
		
		self.pi_text = tk.Text(self.frame, height= 1, width = 30, bd  = 0, takefocus = 0, bg = self.frame.cget('bg'), relief = tk.FLAT)
		self.pi_text.grid(column = 4, row = row, columnspan = 2); self.pi_row = row
		self.pi_text.insert(1.0, 'pi = 3.141593 2pi = 6.283185')
		self.pi_text.configure(state = "disabled")
		
		row += 1
		
		# Min
		self.ymin_entry = tk.Entry(self.frame, width = 17, textvariable = self.ymin, font = self.textFont)
		self.ymin_entry.grid(column = 2, row = row, columnspan = 2 )
		tk.Label(self.frame, text = "  Min", font = self.labelFont).grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.ymax_entry = tk.Entry(self.frame, width = 17, textvariable = self.ymax, font = self.textFont)
		self.ymax_entry.grid(column = 4, row = row, columnspan = 2 )
		tk.Label(self.frame, text = "Max", font = self.labelFont).grid(column = 4, row = row, sticky = tk.W )

		# Ny -> Nr
		row += 1
		self.Ny_entry = tk.Entry(self.frame, width = 17, textvariable = self.Ny, font = self.textFont)
		self.Ny_entry.grid(column = 2, row = row, columnspan = 2 )
		tk.Label(self.frame, text = "     #", font = self.labelFont).grid(column = 2, row = row, sticky = tk.W )


	def set_toroidatTurn_element(self, row):
		# --- toroidal turns ---
		self.itt = tk.StringVar()
		tk.Entry(self.frame, width = 7, textvariable = self.itt, font = self.textFont).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "tor. Iterations", font = self.labelFont).grid(column = 1, row = row, sticky = tk.E )
		return row + 1
	

	def set_phistart_element(self, row):
		# --- phistart ---
		self.phistart = tk.StringVar()
		tk.Entry(self.frame, width = 7, textvariable = self.phistart, font = self.textFont).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "tor. Angle [deg]", font = self.labelFont).grid(column = 1, row = row, sticky = tk.E )
		tk.Label(self.frame, text = "For Machine coord. use negative angles", font = self.labelFont).grid(column = 3, row = row, columnspan = 3, sticky = tk.W )
		return row + 1


	def set_MapDirection_element(self, row, state = tk.NORMAL):
		self.MapDirection = tk.IntVar()
		tk.Radiobutton(self.frame, text = '+1', variable = self.MapDirection, value = 1, state = state, font = self.labelFont).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = '-1', variable = self.MapDirection, value = -1, state = state, font = self.labelFont).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'both', variable = self.MapDirection, value = 0, state = state, font = self.labelFont).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "Map Direction", font = self.labelFont).grid(column = 1, row = row, sticky = tk.E )
		return row + 1
		
		
	def set_runMode_elements(self, row, SIESTA = True):
		okayCommand = self.frame.register(self.okay_response)
		notOkayCommand = self.frame.register(self.notokay_response)
		fileFound = self.frame.register(self.file_Found)

		self.selectField = tk.IntVar()
		self.useM3DC1 = tk.IntVar()
		self.resp = tk.IntVar()
		self.response = tk.IntVar()
		self.nresp = tk.StringVar()
		self.SIESTA = SIESTA
		self.specialFieldFile1 = tk.StringVar()
		self.specialFieldFile2 = tk.StringVar()
		self.specialFieldFile3 = tk.StringVar()
		self.specialFieldFile4 = tk.StringVar()
		self.GPECscale = tk.StringVar()

		# --- separator ---
		separator2 = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------", font = self.labelFont)
		separator2.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- magnetic field ---
		row += 1
		tk.Radiobutton(self.frame, text = 'g-file Vacuum', variable = self.selectField, value = -1,
			command = self.activate_response, font = self.labelFont).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, width = 8, text = 'M3D-C1', variable = self.selectField, value = 0,
			command = self.activate_response, font = self.labelFont).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'VMEC', variable = self.selectField, value = -3,
			command = self.activate_response, font = self.labelFont).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, width = 9, text = 'GPEC', variable = self.selectField, value = -4,
			command = self.activate_response, font = self.labelFont).grid(column = 5, row = row, sticky = tk.W + tk.E )
		#rb = tk.Radiobutton(self.frame, width = 9, text = 'SIESTA', variable = self.selectField, value = -2,
		#	command = self.activate_response, font = self.labelFont)
		#rb.grid(column = 5, row = row, sticky = tk.W + tk.E )	# SIESTA button removed and replaced by GPEC
		tk.Label(self.frame, text = "Field", font = self.labelFont).grid(column = 1, row = row, sticky = tk.E )
		#if not self.SIESTA: rb.configure(state = tk.DISABLED)

		# --- M3DC1 ---
		row += 1; self.row_M3DC1 = row
		self.useM3DC1_R1 = tk.Radiobutton(self.frame, text = 'Equilibrium', variable = self.useM3DC1, value = 0, font = self.labelFont)
		self.useM3DC1_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.useM3DC1_R2 = tk.Radiobutton(self.frame, text = 'Pert.', variable = self.useM3DC1, value = 1, font = self.labelFont)
		self.useM3DC1_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.useM3DC1_R3 = tk.Radiobutton(self.frame, text = 'Eq + Pert', variable = self.useM3DC1, value = 2, font = self.labelFont)
		self.useM3DC1_R3.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.useM3DC1_label = tk.Label(self.frame, text = "use M3DC1", font = self.labelFont)
		self.useM3DC1_label.grid(column = 1, row = row, sticky = tk.E )

		row += 1
		self.response_R1 = tk.Radiobutton(self.frame, text = 'Off', variable = self.resp, value = 0, command = self.make_response, font = self.labelFont)
		self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.response_R2 = tk.Radiobutton(self.frame, text = 'On', variable = self.resp, value = 1, command = self.make_response, font = self.labelFont)
		self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.response_label = tk.Label(self.frame, text = "Plasma response", font = self.labelFont)
		self.response_label.grid(column = 1, row = row, sticky = tk.E )

		self.nresp.set('0');
		self.nresp_entry = tk.Entry(self.frame, width = 4, textvariable = self.nresp, validate = 'all',
        	validatecommand = (okayCommand, '%P'), invalidcommand = notOkayCommand, font = self.textFont)
		self.nresp_entry.grid(column = 4, row = row, sticky = tk.E)
		self.nresp_label = tk.Label(self.frame, text = "       Time", font = self.labelFont)
		self.nresp_label.grid(column = 4, row = row, sticky = tk.W )

		# --- VMEC ---
		row = self.row_M3DC1
		#self.specialFieldFile1_entry = tk.Entry(self.frame, width = 10, textvariable = self.specialFieldFile1, 
		#	validate = 'focusout', validatecommand = (fileFound, '%d','%P'))
		self.specialFieldFile1_entry = AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
				width = 50, textvariable = self.specialFieldFile1, font = self.textFont)
		self.specialFieldFile1_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
		self.specialFieldFile1_label = tk.Label(self.frame, text = "VMEC wout: ", font = self.labelFont)
		self.specialFieldFile1_label.grid(column = 1, row = row, sticky = tk.E )
		#self.specialFieldFile1_entry.bind('<Return>', lambda event: self.file_Found(event,self.specialFieldFile1.get()))
		
		row += 1

		#self.specialFieldFile2_entry = tk.Entry(self.frame, width = 10, textvariable = self.specialFieldFile2, 
		#	validate = 'focusout', validatecommand = (fileFound, '%d','%P'))
		self.specialFieldFile2_entry = AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
				width = 50, textvariable = self.specialFieldFile2, font = self.textFont)
		self.specialFieldFile2_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
		self.specialFieldFile2_label = tk.Label(self.frame, text = "Xpand data: ", font = self.labelFont)
		self.specialFieldFile2_label.grid(column = 1, row = row, sticky = tk.E )
		#self.specialFieldFile2_entry.bind('<Return>', lambda event: self.file_Found(event,self.specialFieldFile2.get()))

		# --- SIESTA ---
		# This code is basically dead and will probably never be used again
		row = self.row_M3DC1
		#self.specialFieldFile3_entry = tk.Entry(self.frame, width = 10, textvariable = self.specialFieldFile3, 
		#	validate = 'focusout', validatecommand = (fileFound, '%d','%P'))
		#self.specialFieldFile3_entry = AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
		#		width = 50, textvariable = self.specialFieldFile3, font = self.textFont)
		#self.specialFieldFile3_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
		#self.specialFieldFile3_label = tk.Label(self.frame, text = "SIESTA data: ", font = self.labelFont)
		#self.specialFieldFile3_label.grid(column = 1, row = row, sticky = tk.E )
		#self.specialFieldFile3_entry.bind('<Return>', lambda event: self.file_Found(event,self.specialFieldFile3.get()))

		# --- GPEC ---
		# This replaces the old SIESTA selection
		row = self.row_M3DC1
		self.GPECscale.set('1');
		self.specialFieldFile4_entry = AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
				width = 50, textvariable = self.specialFieldFile4, font = self.textFont)
		self.specialFieldFile4_entry.grid(column = 2, row = row, columnspan = 3,sticky = tk.E+tk.W)
		self.specialFieldFile4_label = tk.Label(self.frame, text = "GPEC data: ", font = self.labelFont)
		self.specialFieldFile4_label.grid(column = 1, row = row, sticky = tk.E )
		self.GPECscale_entry = tk.Entry(self.frame, width = 4, textvariable = self.GPECscale, font = self.textFont)
		self.GPECscale_entry.grid(column = 5, row = row, sticky = tk.E, padx=15)
		self.GPECscale_label = tk.Label(self.frame, text = "  Scale", font = self.labelFont)
		self.GPECscale_label.grid(column = 5, row = row, sticky = tk.W)

		return self.row_M3DC1 + 2


	def set_particle_elements(self, row):
		self.sigma = tk.IntVar()
		self.charge = tk.IntVar()
		self.Ekin = tk.StringVar()
		self.Lambda = tk.StringVar()
		self.Mass = tk.StringVar()
		self.useFilament = tk.StringVar()
		self.nproc = tk.StringVar()
		self.Ekin_array = None
		self.Ekin_idx = 0
		self.ErProfileFile = tk.StringVar()

		# --- separator ---
		separator3 = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------", font = self.labelFont)
		separator3.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- particles ---
		row += 1
		tk.Radiobutton(self.frame, text = 'Field lines', variable = self.sigma, value = 0, 
			command = self.show_particle_params, font = self.labelFont).grid(column = 2, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(self.frame, text = 'Co pass', variable = self.sigma, value = 1,
			command = self.show_particle_params, font = self.labelFont).grid(column = 3, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(self.frame, text = 'Counter pass', variable = self.sigma, value = -1,
			command = self.show_particle_params, font = self.labelFont).grid(column = 4, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Label(self.frame, text = "Orbits", font = self.labelFont).grid(column = 1, row = row, sticky = tk.E + tk.N)

		row += 1; self.row_particle = row
		self.charge_R1 = tk.Radiobutton(self.frame, text = 'Electrons', variable = self.charge, 
			value = -1, command = self.setMass, font = self.labelFont)
		self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.charge_R2 = tk.Radiobutton(self.frame, text = 'Ions', variable = self.charge, 
			value = 1, command = self.setMass, font = self.labelFont)
		self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.charge_label = tk.Label(self.frame, text = "Species", font = self.labelFont)
		self.charge_label.grid(column = 1, row = row, sticky = tk.E )

		self.Mass_entry = tk.Entry(self.frame, width = 7, textvariable = self.Mass, font = self.textFont)
		self.Mass_entry.grid(column = 4, row = row, sticky = tk.E )
		self.Mass_label = tk.Label(self.frame, text = "   Mass", font = self.labelFont)
		self.Mass_label.grid(column = 4, row = row, sticky = tk.W )

		row += 1
		self.Ekin_entry = tk.Entry(self.frame, width = 7, textvariable = self.Ekin, font = self.textFont)
		self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.Ekin_label = tk.Label(self.frame, text = "kin Energy [keV]", font = self.labelFont)
		self.Ekin_label.grid(column = 1, row = row, sticky = tk.E )

		self.Lambda_entry = tk.Entry(self.frame, width = 7, textvariable = self.Lambda, font = self.textFont)
		self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.Lambda_label = tk.Label(self.frame, text = "Energy ratio", font = self.labelFont)	
		self.Lambda_label.grid(column = 3, row = row, sticky = tk.E )
		
		row += 1
		self.ErProfileFile_entry = AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
				width = 50, textvariable = self.ErProfileFile, font = self.textFont)
		self.ErProfileFile_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
		self.ErProfileFile_label = tk.Label(self.frame, text = "Er Profile: ", font = self.labelFont)
		self.ErProfileFile_label.grid(column = 1, row = row, sticky = tk.E )
		
		# --- separator ---
		row += 1
		separator4 = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------", font = self.labelFont)
		separator4.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- Filament ---
		row += 1
		tk.Entry(self.frame, width = 7, textvariable = self.useFilament, font = self.textFont).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "# current Filaments", font = self.labelFont).grid(column = 1, row = row, sticky = tk.E )

		return row + 1
		
		
	def set_run_mpi_elements(self, row):
		# --- number of processes for mpi ---
		self.nproc.set(str(4))
		self.nproc_entry = tk.Entry(self.frame, width = 4, textvariable = self.nproc, font = self.textFont)
		self.nproc_entry.grid(column = 1, row = row, sticky = tk.E)
		tk.Label(self.frame, text = "       # Procs", font = self.labelFont).grid(column = 1, row = row, sticky = tk.W )


	def set_bottom_elements(self, row, SIESTA = True):
		# --- separator ---
		separator = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------", font = self.labelFont)
		separator.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		row += 4
		
		# --- Coils ---
		self.define_machine_elements()
				
		# --- Modes ---
		row = self.set_runMode_elements(row, SIESTA = SIESTA)
		
		# --- Particles ---
		row = self.set_particle_elements(row)
		
		return row


	# --- turn on/off MinMax entry ---
	def activate_MinMax_entrys(self, flag):
		if flag:
			self.xmin_entry.configure(state=tk.DISABLED)
			self.xmax_entry.configure(state=tk.DISABLED)
			self.Nx_entry.configure(state=tk.DISABLED)
			self.ymin_entry.configure(state=tk.DISABLED)
			self.ymax_entry.configure(state=tk.DISABLED)
			self.Ny_entry.configure(state=tk.DISABLED)
		else:
			self.xmin_entry.configure(state=tk.NORMAL)
			self.xmax_entry.configure(state=tk.NORMAL)
			self.Nx_entry.configure(state=tk.NORMAL)
			self.ymin_entry.configure(state=tk.NORMAL)
			self.ymax_entry.configure(state=tk.NORMAL)
			self.Ny_entry.configure(state=tk.NORMAL)


	# --- turn on/off response radiobutton ---
	def activate_response(self):
		self.activate_M3DC1_response()
		self.activate_VMEC_response()
		self.activate_GPEC_response()


	def activate_M3DC1_response(self):
		if(self.selectField.get() < 0):
			self.useM3DC1.set(0)
			self.useM3DC1_R1.grid_forget()
			self.useM3DC1_R2.grid_forget()
			self.useM3DC1_R3.grid_forget()
			self.useM3DC1_label.grid_forget()
			self.response_R1.grid_forget()
			self.response_R2.grid_forget()
			self.response_label.grid_forget()
			self.nresp_entry.grid_forget()
			self.nresp_label.grid_forget()
		else:
			row = self.row_M3DC1
			self.useM3DC1_R1.grid(column = 2, row = row, sticky = tk.W + tk.E, padx=5, pady=5)
			self.useM3DC1_R2.grid(column = 3, row = row, sticky = tk.W + tk.E, padx=5, pady=5)
			self.useM3DC1_R3.grid(column = 4, row = row, sticky = tk.W + tk.E, padx=5, pady=5)
			self.useM3DC1_label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
			row = self.row_M3DC1 + 1
			self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E, padx=5, pady=5)
			self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E, padx=5, pady=5)
			self.response_label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
			self.nresp_entry.grid(column = 4, row = row, sticky = tk.E, padx=5, pady=5)
			self.nresp_label.grid(column = 4, row = row, sticky = tk.W, padx=5, pady=5)


	def activate_VMEC_response(self):
		if(self.selectField.get() == -3):
			self.specialFieldFile1_label.grid(column = 1, row = self.row_M3DC1, sticky = tk.E, padx=5, pady=5)
			self.specialFieldFile1_entry.grid(column = 2, row = self.row_M3DC1, columnspan = 4, sticky = tk.E+tk.W, padx=5, pady=5)
			self.specialFieldFile2_label.grid(column = 1, row = self.row_M3DC1+1, sticky = tk.E, padx=5, pady=5)
			self.specialFieldFile2_entry.grid(column = 2, row = self.row_M3DC1+1, columnspan = 4, sticky = tk.E+tk.W, padx=5, pady=5)
		else:
			self.specialFieldFile1_label.grid_forget()
			self.specialFieldFile1_entry.grid_forget()
			self.specialFieldFile2_label.grid_forget()
			self.specialFieldFile2_entry.grid_forget()
			self.specialFieldFile1.set('')
			self.specialFieldFile2.set('')


	def activate_GPEC_response(self):
		if(self.selectField.get() == -4):
			self.specialFieldFile4_label.grid(column = 1, row = self.row_M3DC1, sticky = tk.E, padx=5, pady=5)
			self.specialFieldFile4_entry.grid(column = 2, row = self.row_M3DC1, columnspan = 3, sticky = tk.E+tk.W, padx=5, pady=5)
			self.GPECscale_entry.grid(column = 5, row = self.row_M3DC1, sticky = tk.E, padx=15)
			self.GPECscale_label.grid(column = 5, row = self.row_M3DC1, sticky = tk.W )
			
			path = os.path.abspath(self.path.get())
			if not (path[-1] == '/'): path += '/'
			if os.path.isfile(path + 'gpecsup.in'):
				with open(path + 'gpecsup.in') as f:
					input = f.readlines()
				data = input[0].strip().split()
				self.specialFieldFile4.set(data[0])
				self.GPECscale.set(data[1])
				if self.specialFieldFile4_entry.listboxUp:
					self.specialFieldFile4_entry.listbox.destroy()
					self.specialFieldFile4_entry.tw.destroy()
					self.specialFieldFile4_entry.listboxUp = False

		else:
			self.specialFieldFile4_label.grid_forget()
			self.specialFieldFile4_entry.grid_forget()
			self.GPECscale_label.grid_forget()
			self.GPECscale_entry.grid_forget()
			self.specialFieldFile4.set('')
			self.GPECscale.set('1.0')


	def activate_SIESTA_response(self):
		pass
# 		if(self.selectField.get() == -2):
# 			self.specialFieldFile3_label.grid(column = 1, row = self.row_M3DC1, sticky = tk.E, padx=5, pady=5)
# 			self.specialFieldFile3_entry.grid(column = 2, row = self.row_M3DC1, columnspan = 4, sticky = tk.E+tk.W, padx=5, pady=5)
# 			self.createFlag.set('psi')
# 			self.refresh_grid_labels()
# 			self.create_R1.configure(state=tk.DISABLED)
# 			self.create_R2.configure(text = 's,u')
# 			try: self.create_R3.configure(state=tk.DISABLED)
# 			except: pass
# 		else:
# 			self.specialFieldFile3_label.grid_forget()
# 			self.specialFieldFile3_entry.grid_forget()
# 			self.specialFieldFile3.set('')
# 			if self.SIESTA:
# 				self.refresh_grid_labels()
# 				self.create_R1.configure(state=tk.NORMAL)
# 				self.create_R2.configure(text = 'psi_n')
# 				try: self.create_R3.configure(state=tk.NORMAL)
# 				except: pass


	def make_VMEC_SIESTA_shell_flags(self):
		flag = ''
		if(self.selectField.get() == -3):
			vmecfile = self.specialFieldFile1.get()
			xpandfile = self.specialFieldFile2.get()
			if self.file_Found(file = vmecfile):
				flag += ' -V ' + vmecfile
			if self.file_Found(file = xpandfile):
				flag += ' -X ' + xpandfile
		elif(self.selectField.get() == -2):
			siestafile = self.specialFieldFile3.get()
			if self.file_Found(file = siestafile):
				flag += ' -S ' + siestafile	
		return flag
		
	
	def write_sup_dot_in_file(self):
		path = os.path.abspath(self.path.get())
		if not (path[-1] == '/'): path += '/'
		if(self.selectField.get() == -4):
			with open(path + 'gpecsup.in', 'w') as f:
				f.write(self.specialFieldFile4.get() + '\t' + self.GPECscale.get() + '\n')
				
			
	# --- set response variable ---
	def make_response(self):
		if(self.resp.get() == 0):
			self.nresp_entry.configure(state=tk.DISABLED)
			self.nresp_label.configure(state=tk.DISABLED)
			self.response.set(0)
			self.nresp.set('0')
		else:
			self.nresp_entry.configure(state=tk.NORMAL)
			self.nresp_label.configure(state=tk.NORMAL)
			
			path = self.path_to_m3dc1()
			#print('Searching for M3D-C1 files in ', path)
			nresp = 1
			while(os.path.isfile(path + 'time_' + format(nresp,'03d') + '.h5')):
				nresp += 1
			nresp -= 1
				
			if (not os.path.isfile(path + 'time_000.h5')): print('Warning, no time_xxx.h5 file found in', path )
			#else: print('Highest time_xxx.h5 file available is: time_' + format(nresp,'03d') + '.h5')
			self.nresp.set(str(nresp))
			self.response.set(int(self.nresp.get()))


	# --- validity check of response ---
	def path_to_m3dc1(self):
		path = os.path.abspath(self.path.get())
		if not (path[-1] == '/'): path += '/'
		if os.path.isfile(path + 'm3dc1sup.in'):
			with open(path + 'm3dc1sup.in') as f:
				input = f.readlines()
			c1 = input[0].strip().split()[0]
			idx = c1[::-1].find('/')	# returns location of last '/' in c1 or -1 if not found
			if(idx == -1): 
				c1path = path
			else:
				idx *= -1
				c1path = c1[0:idx]	# path with a final '/'
				if not (c1path[0] == '/'): c1path = path + c1path	# c1path is a relative path, make it absolute using location of m3dc1sup.in
		else: c1path = path
		c1path = os.path.abspath(c1path)
		if not (c1path[-1] == '/'): c1path += '/'
		return c1path
		
		
	# returns True if time_xxx.h5 file exists
	# else returns False
	def okay_response(self, nresp):
		path = self.path_to_m3dc1()
		try: 
			nresp = int(nresp)
		except: 
			if not nresp:
				return True	# empty string
			else: 
				return False
		
		if os.path.isfile(path + 'time_' + format(nresp,'03d') + '.h5'):
			if nresp == 0: self.resp.set(0)
			else: self.resp.set(1)
			self.response.set(nresp)
			return True
		elif nresp == 0: 
			return True	# not even time_000.h5 exists
		else: 
			return False


	# sets response to highest, if okay_response == False
	def notokay_response(self):
		print ('invalid Entry: requested time_xxx.h5 file not found in', self.path_to_m3dc1())


	# --- Validate the file entries ---
	def file_Found(self, event = None, file = ''):
		if len(file) == 0: return False
		if not '/' in file:
			file = os.path.abspath(self.path.get()) + '/' + file
		if not os.path.isfile(file):
			print ('Requested file not found:', file)
			return False	# file not found
		else:
			#self.specialFieldFile1_entry.configure(validate = 'focusout')
			#self.specialFieldFile2_entry.configure(validate = 'focusout')
			return True		# file found
		
	
	# --- Show or Hide Particle Options, depending on sigma ---
	def show_particle_params(self):
		if not (self.sigma.get() == 0):
			row = self.row_particle
			self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.charge_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			self.Mass_entry.grid(column = 4, row = row, sticky = tk.E , padx=5, pady=5)
			self.Mass_label.grid(column = 4, row = row, sticky = tk.W , padx=5, pady=5)
			row = self.row_particle + 1
			self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Ekin_label.grid(column = 1, row = row, sticky = tk.E , padx=5, pady=5)
			self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E , padx=5, pady=5)
			self.Lambda_label.grid(column = 3, row = row, sticky = tk.E , padx=5, pady=5)
			row = self.row_particle + 2
			self.ErProfileFile_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W, padx=5, pady=5)
			self.ErProfileFile_label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		else:
			self.charge_R1.grid_forget()
			self.charge_R2.grid_forget()
			self.charge_label.grid_forget()
			self.Ekin_entry.grid_forget()
			self.Ekin_label.grid_forget()
			self.Lambda_entry.grid_forget()
			self.Lambda_label.grid_forget()
			self.Mass_entry.grid_forget()
			self.Mass_label.grid_forget()
			self.ErProfileFile_entry.grid_forget()
			self.ErProfileFile_label.grid_forget()
			self.ErProfileFile.set('')


	def make_ErProfile_shell_flag(self):
		flag = ''
		ErFile = self.ErProfileFile.get()
		if self.file_Found(file = ErFile):
			flag = ' -E ' + ErFile	
		return flag


	def set_Ekin(self):
		Ekin = self.Ekin.get()
		if ',' in Ekin: 
			Ekin = Ekin.split(',')
			self.Ekin_array = np.array([np.float64(item) for item in Ekin])
		elif ':' in Ekin: 
			Ekin = Ekin.split(':')
			if len(Ekin) > 2: self.Ekin_array = np.arange(np.float64(Ekin[0]), np.float64(Ekin[1]) + 0.5*np.float64(Ekin[2]), np.float64(Ekin[2]))
			else: self.Ekin_array = np.arange(np.float64(Ekin[0]), np.float64(Ekin[1]) + 0.1, 0.25)
		else: 
			self.Ekin_array = np.array([np.float64(Ekin)])
			
			
	def setMass(self):
		if self.charge.get() == -1:
			self.Mass_entry.configure(state=tk.DISABLED)
			self.Mass.set('1')
		else:
			self.Mass_entry.configure(state=tk.NORMAL)
			self.Mass.set('2')


	# --- Default Functions --------------------------------------------------------------

	# --- read parameter file ---
	def read_par_file(self, names):
		# --- read parameterfile, if it is there ---
		for name in names:
			if os.path.isfile(self.path.get() + name):
				_,_,_, data, head = readControlFile(self.path.get() + name, self.gPath.get())
				if self.translateMachFlag[self.MachFlag.get()] in head: break
				else: data = None
			else: data = None
		return data


	# --- set default values ---
	def set_common_defaults(self, data):
		self.itt.set(data['itt'])
		self.phistart.set(data['phistart'])
		self.MapDirection.set(data['MapDirection'])
		
		if (data['useM3DC1'] >= 0):
			self.selectField.set(0)
			self.useM3DC1.set(data['useM3DC1'])
		else:
			self.selectField.set(data['selectField'])
			self.useM3DC1.set(0)
			
		self.response.set(data['response'])
		if self.response.get() == 0: 
			self.resp.set(0)
			self.nresp.set('0')
			self.nresp_entry.configure(state=tk.DISABLED)
			self.nresp_label.configure(state=tk.DISABLED)
		else:
			self.resp.set(1)
			self.nresp.set(str(self.response.get()))
			self.nresp_entry.configure(state=tk.NORMAL)
			self.nresp_label.configure(state=tk.NORMAL)
		self.sigma.set(data['sigma'])
		self.charge.set(data['charge'])
		self.Ekin.set(data['Ekin'])
		self.Lambda.set(data['Lambda'])
		self.Mass.set(data['Mass'])
		self.useFilament.set(data['useFilament'])
		
		self.activate_response()
		self.show_particle_params()
		
		self.specialFieldFile1_entry.update_baseSearchPath(self.path.get())
		self.specialFieldFile2_entry.update_baseSearchPath(self.path.get())
		#self.specialFieldFile3_entry.update_baseSearchPath(self.path.get())
		self.specialFieldFile4_entry.update_baseSearchPath(self.path.get())
		self.ErProfileFile_entry.update_baseSearchPath(self.path.get())
		

	# --- Function, executed when Button is pressed ---
	def run_common_funct(self, name, shellCall, data, shellFlags = ''):
		if(self.Shot.get() == '') | (self.Time.get() == ''):
			print ('You must enter a Shot # and a Time')
			return

		# make gPath absolute & check if gPath ends with a /  if so, append the conventional gfile name 
		self.gPath.set(os.path.abspath(self.gPath.get()))
		self.gFile_entry.selection(None)
		#if (self.gPath.get()[-1] == '/'): self.gPath.set(self.gPath + 'g' + format(int(self.Shot.get()),'06d') + '.' + format(int(self.Time.get()),'05d'))
	
		# convert relative path to absolute
		cwd = os.getcwd()
		path = os.path.abspath(self.path.get())
		if not (path[-1] == '/'): path += '/'
		chk = not (cwd + '/' == path)
				
		# set flags
		# do this before the os.chdir, because the os.path.abspath recognizes the chdir, but self.path stays the same
		shellFlags += self.make_VMEC_SIESTA_shell_flags()
		shellFlags += self.make_ErProfile_shell_flag()
			
		# change to working dir, write contol file(s), launch code, and return to original dir
		if chk: os.chdir(path)
		
		self.set_Ekin()
		NE = len(self.Ekin_array)
		if NE > 1: 
			sc = shellCall.split(name)
			shellCall_i = []
			shellCall_slurm = sc[0] + name[0:-4] + '_${SLURM_ARRAY_TASK_ID}' + '.dat' + sc[1] + '_${SLURM_ARRAY_TASK_ID}'
		else:
			shellCall_slurm = shellCall

		for i in range(NE):
			if NE > 1: 
				name_i = name[0:-4] + '_' + str(i) + '.dat'
				shellCall_i.append(sc[0] + name_i + sc[1] + '_' + str(i))
			else: 
				name_i = name
				shellCall_i = [shellCall]
			self.Ekin_idx = i
			self.writeControlFile(name_i)
		
		self.write_sup_dot_in_file()
		
		if(('iris' in HOST) | ('pppl.gov' in HOST) | ('omega' in HOST)):			# GA Iris or Omega Cluster & PPPL Portal cluster
			self.write_qsub_file(data, self.tag.get(), shellCall_slurm + shellFlags)
			#call('sbatch ./mafot.sbatch', shell = True)
			print ('To run, type: sbatch ./mafot.sbatch')
		else:
			#for i in range(NE): call(shellCall_i[i] + shellFlags + ' &', shell = True)
			#print 'running in dir:', os.getcwd()
			for i in range(NE): print ('To run, type: ' + shellCall_i[i] + shellFlags + ' &')
			
			
		if chk: os.chdir(cwd)


	# --- Write qsub File on Clusters ---
	# here: shellFlags already embedded in shellCall
	def write_common_qsub_file(self, tooltag, nproc, tag, shellCall, mpi = True):
		import getpass
		user = getpass.getuser()
		if('iris' in HOST):			# Iris Cluster
			with open('mafot.sbatch', 'w') as f:
				f.write('#!/bin/bash' + '\n')
				if len(self.Ekin_array) > 1:
					f.write('#SBATCH --array=0-' + str(len(self.Ekin_array)-1) + '\n')
				#f.write('#SBATCH --job-name=mafot' + tooltag + tag.translate(None, '_+- ') + '\n') # python 2.7 syntax
				f.write('#SBATCH --job-name=mafot' + tooltag + tag.translate(dict.fromkeys('_+- ')) + '\n')
				f.write('#SBATCH -p preemptable' + '\n')
				f.write('#SBATCH -o batch_mafot.out' + '\n')
				if mpi:
					f.write('#SBATCH -n ' + str(nproc) + ' \n')
				f.write('#SBATCH -t 120' + '\n')
				f.write('#SBATCH --mem-per-cpu=1G' + '\n')
				if (user == 'wingen') & (len(self.Ekin_array) > 1): f.write('##SBATCH --mail-user=' + user + '@fusion.gat.com' + '\n')
				else: f.write('#SBATCH --mail-user=' + user + '@fusion.gat.com' + '\n')
				if len(self.Ekin_array) > 1:
					if user == 'wingen': f.write('##SBATCH --mail-type=FAIL,ARRAY_TASKS  ## or BEGIN,END,FAIL,ALL,ARRAY_TASKS ' + '\n')
					else: f.write('#SBATCH --mail-type=FAIL,ARRAY_TASKS  ## or BEGIN,END,FAIL,ALL,ARRAY_TASKS ' + '\n')
				else:
					f.write('#SBATCH --mail-type=FAIL  ## or BEGIN,END,FAIL,ALL ' + '\n')
				f.write('#SBATCH --export=ALL' + '\n')
				f.write('##module purge' + '\n')
				f.write('module load default-paths' + '\n')
				f.write('module load mafot' + '\n')
				f.write('module list' + '\n')
				f.write(shellCall + '\n')
		if('omega' in HOST):			# Omega Cluster
			with open('mafot.sbatch', 'w') as f:
				f.write('#!/bin/bash' + '\n')
				if len(self.Ekin_array) > 1:
					f.write('#SBATCH --array=0-' + str(len(self.Ekin_array)-1) + '\n')
				#f.write('#SBATCH --job-name=mafot' + tooltag + tag.translate(None, '_+- ') + '\n') # python 2.7 syntax
				f.write('#SBATCH --job-name=mafot' + tooltag + tag.translate(dict.fromkeys('_+- ')) + '\n')
				f.write('#SBATCH -p preemptable' + '\n')
				f.write('#SBATCH -o batch_mafot.out' + '\n')
				if mpi:
					f.write('#SBATCH -n ' + str(nproc) + ' \n')
				f.write('#SBATCH -t 120' + '\n')
				#f.write('#SBATCH --mem-per-cpu=1G' + '\n')
				if (user == 'wingen') & (len(self.Ekin_array) > 1): f.write('##SBATCH --mail-user=' + user + '@fusion.gat.com' + '\n')
				else: f.write('#SBATCH --mail-user=' + user + '@fusion.gat.com' + '\n')
				if len(self.Ekin_array) > 1:
					if user == 'wingen': f.write('##SBATCH --mail-type=FAIL,ARRAY_TASKS  ## or BEGIN,END,FAIL,ALL,ARRAY_TASKS ' + '\n')
					else: f.write('#SBATCH --mail-type=FAIL,ARRAY_TASKS  ## or BEGIN,END,FAIL,ALL,ARRAY_TASKS ' + '\n')
				else:
					f.write('#SBATCH --mail-type=FAIL  ## or BEGIN,END,FAIL,ALL ' + '\n')
				f.write('#SBATCH --export=ALL' + '\n')
				#f.write('export FI_PROVIDER=tcp' + '\n')
				f.write('module purge' + '\n')
				f.write('module load mafot' + '\n')
				f.write('module list' + '\n')
				#sc = shellCall.split()
				#shellCallSrun = 'srun -n ' + str(nproc) + ' --mpi=pmi2 --label ' + ' '.join(sc[3::]) # replace mpirun with srun
				#f.write(shellCallSrun + '\n')  # Work around due to bug in slurm with gcc11 on Omega that fails mpirun command
				f.write(shellCall + '\n')  # This works now again using gcc8 on Omega
		elif('pppl.gov' in HOST):			# PPPL Portal cluster
			with open('mafot.sbatch', 'w') as f:
				f.write('#!/bin/bash -vx' + '\n')
				if len(self.Ekin_array) > 1:
					f.write('#SBATCH --array=0-' + str(len(self.Ekin_array)-1) + '\n')
				f.write('#SBATCH --job-name=mafot' + tooltag + tag.translate(dict.fromkeys('_+- ')) + '\n')
				f.write('#SBATCH -p general' + '\n')
				#f.write('#SBATCH -o batch_mafot.out' + '\n')
				if mpi:
					f.write('#SBATCH -n ' + str(nproc) + ' \n')
				f.write('#SBATCH -t 02:00:00' + '\n')
				f.write('#SBATCH --mem-per-cpu=10000' + '\n')
				f.write(shellCall + '\n')
		else:
			print ('Unknown Host. No Queue submission file created')


	# --- write common Control File parts ------------------------------------------------
	def write_headerLines(self, MachFlag, f):
		if(MachFlag == 'cmod'):
			self.write_cmod_headerLines(f)	# special, because the shot number has 10 digits, instead of 6, like all others
		else:	# all others
			self.write_common_headerLines(f, self.translateMachFlag[MachFlag])


	def write_common_headerLines(self, f, name):
		f.write('# Parameterfile for ' + name + ' Programs\n')
		f.write('# Shot: ' + format(int(self.Shot.get()),'06d') + '\tTime: ' + format(int(self.Time.get()),'04d') + 'ms\n')
		f.write('# Path: ' + self.gPath.get() + '\n')


	def write_Ctrl_resonse(self, f):
		f.write('PlasmaResponse(0=no,>1=yes)=\t' + str(self.response.get()) + '\n')
		f.write('Field(-4=GPEC,-3=VMEC,-2=SIESTA,-1=gfile,M3DC1:0=Eq,1=I-coil,2=both)=\t' + str(self.selectField.get() + self.useM3DC1.get()) + '\n')


	def write_Ctrl_particles(self, f):
		f.write('ParticleDirection(1=co-pass,-1=ctr-pass,0=field-lines)=\t' + str(self.sigma.get()) + '\n')
		f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge.get()) + '\n')
		f.write('Ekin[keV]=\t' + str(self.Ekin_array[self.Ekin_idx]) + '\n')
		f.write('lambda=\t' + self.Lambda.get() + '\n')
		f.write('Mass=\t' + self.Mass.get() + '\n')


	def write_Ctrl_center(self, f):
		f.write('phistart(deg)=\t' + self.phistart.get() + '\n')
		f.write('MapDirection=\t' + str(self.MapDirection.get()) + '\n')
		self.write_Ctrl_resonse(f)
	
	
	def write_Ctrl_bottom(self, f):
			self.write_coils(self.MachFlag.get(), f)
			if self.MachFlag.get() in ['iter', 'nstx', 'mast']:
				f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			if self.MachFlag.get() in ['iter', 'nstx']:
				f.write('useTe_profile(0=no)=	0\n')
			self.write_Ctrl_particles(f)
			if self.MachFlag.get() in ['any']:
				f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			if self.MachFlag.get() in ['dt']:
				f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
				self.write_d3d_errorFileds(f)
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')


	# --- Overloaded Functions -----------------------------------------------------------
	
	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		return

	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, shellCall):
		return


# -------------------------------------------------------------------------------------------------------------
# --- plot ----------------------------------------------------------------------------------------------------
class set_plot_tab(common_tab):		# inherit common tab class
	def __init__(self, frame, M):
		self.frame = frame
		self.MachFlag = M.MachFlag
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.translateMachFlag = M.translateMachFlag
		self.gFile_entry = M.gFile_entry
		self.labelFont = None
		self.textFont = None

		# --- grid  type ---
		# define... 
		self.createFlag = tk.StringVar(); row = 0

		# ...and set grid-type RadioButton
		self.create_R1 = tk.Radiobutton(frame, text = 'RZ', variable = self.createFlag, value = 'RZ', 
			command = self.refresh_grid_labels)
		self.create_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.create_R2 = tk.Radiobutton(frame, text = 'psi_n', variable = self.createFlag, value = 'psi', 
			command = self.refresh_grid_labels)
		self.create_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.create_R3 = tk.Radiobutton(frame, text = 'Polar', variable = self.createFlag, value = 'polar', 
			command = self.refresh_grid_labels)
		self.create_R3.grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Coordinate Type").grid(column = 1, row = row, sticky = tk.E )

		# --- x -> theta or R ---
		row += 1
		self.x_label = tk.Label(frame, text = "")
		self.x_label.grid(column = 1, row = row, sticky = tk.E)
		
		self.set_MinMax_elements(row)
		
		# --- y -> r, psi or Z ---
		row += 2
		self.y_label = tk.Label(frame, text = "")
		self.y_label.grid(column = 1, row = row, sticky = tk.E )

		row += 2
		
		# --- toroidal turns ---
		row = self.set_toroidatTurn_element(row)

		# --- phistart ---
		row = self.set_phistart_element(row)
		
		# --- MapDirection ---
		row = self.set_MapDirection_element(row)
		
		# --- Coils & Mode & Particles ---
		row = self.set_bottom_elements(row)		
		
		# --- Run ---
		self.set_run_mpi_elements(row)
		runButton = tk.Button(frame, text = "Run", command = self.run_funct)
		runButton.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)

		# --- adjust style ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.set_defaults()


	# --- Change Labels on grid variables, depending on createFlag ---
	def refresh_grid_labels(self):
		if(self.createFlag.get() == 'polar'):
			self.x_label.configure(text = "theta [rad]")
			self.y_label.configure(text = "r [m]")
			self.pi_text.grid(column = 4, row = self.pi_row, columnspan = 2)
		elif(self.createFlag.get() == 'psi'):
			if(self.selectField.get() == -2):
				self.x_label.configure(text = "u [rad]")
				self.y_label.configure(text = "s")
			else:
				self.x_label.configure(text = "theta [rad]")
				self.y_label.configure(text = "psi_n")
			self.pi_text.grid(column = 4, row = self.pi_row, columnspan = 2)
		else:
			self.x_label.configure(text = "R [m]")
			self.y_label.configure(text = "Z [m]")
			self.pi_text.grid_forget()
	

	# --- define machine specific default values ---
	def tool_defaults(self, flag):	
		if self.MachFlag.get() == 'any':
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 1, 100, 0.1, 2, 0,
					3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'iter':
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 2,
					3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0, 1, 
					100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 0, -1, 1, 3, 1, 1, 0, 0, 1, 
					100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		else:
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 0, -1, 1, 3, 1, 1, 1, 0, 1, 
					100, 0.1, 2, 0, 0, 0, 3.141592653589793, 6.283185307179586]
		return data


	# --- set default values ---
	def set_defaults(self):
		data = self.read_par_file(['_plot.dat'])
		data = self.set_machine_defaults(data, self.MachFlag.get())
										
		if(data['create'] == 0) | (data['create'] == 1): self.createFlag.set('polar')
		elif(data['create'] == 3) | (data['create'] == 4): self.createFlag.set('psi')
		else: self.createFlag.set('RZ')
		
		if(self.createFlag.get() == 'RZ'): 
			self.xmin.set(repr(data['0-8'][2]))
			self.xmax.set(repr(data['0-8'][3]))
			self.ymin.set(repr(data['0-8'][4]))
			self.ymax.set(repr(data['0-8'][5]))
			self.Nx.set(str(int(data['0-8'][6])))
			self.Ny.set(str(int(data['0-8'][0])))
		else: 
			self.xmin.set(repr(data['0-8'][4]))
			self.xmax.set(repr(data['0-8'][5]))
			self.ymin.set(repr(data['0-8'][2]))
			self.ymax.set(repr(data['0-8'][3]))
			self.Nx.set(str(int(data['0-8'][0])))
			self.Ny.set(str(int(data['0-8'][6])))
			
		self.refresh_grid_labels()
		

	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		tool = self.MachFlag.get() + 'plot_mpi'
		name = '_plot.dat'
		shellCall = MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' ' + tool + ' ' + name + ' ' + self.tag.get()
		self.run_common_funct(name, shellCall, int(self.nproc.get()))


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, shellCall):
		tooltag = 'P'
		self.write_common_qsub_file(tooltag, nproc, tag, shellCall, mpi = True)


	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		with open(name, 'w') as f:
			self.write_headerLines(self.MachFlag.get(), f)
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
			self.write_Ctrl_center(f)
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')
			if(self.createFlag.get() == 'polar'): f.write('createPoints(0=setr,3=setpsi,5=setR)=\t0\n')
			elif(self.createFlag.get() == 'psi'): f.write('createPoints(0=setr,3=setpsi,5=setR)=\t3\n')
			else: f.write('createPoints(0=setr,3=setpsi,5=setR)=\t5\n')
			self.write_Ctrl_bottom(f)


# -------------------------------------------------------------------------------------------------------------
# --- fix -----------------------------------------------------------------------------------------------------
class set_fix_tab(common_tab):		# inherit common tab class
	def __init__(self, frame, M):
		self.frame = frame
		self.MachFlag = M.MachFlag
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.translateMachFlag = M.translateMachFlag
		self.gFile_entry = M.gFile_entry
		self.labelFont = None
		self.textFont = None
		
		# --- unused Variables ---
		self.itt = tk.StringVar(); self.itt.set('0')
		self.MapDirection = tk.IntVar(); self.MapDirection.set(0)

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
		self.x_label = tk.Label(frame, text = "R [m]")
		self.x_label.grid(column = 1, row = row, sticky = tk.E)
		
		self.set_MinMax_elements(row)
		
		# --- y -> r ---
		row += 2
		self.y_label = tk.Label(frame, text = "Z [m]")
		self.y_label.grid(column = 1, row = row, sticky = tk.E )
		 
		row += 2
				
		# --- grid  type ---
		self.createFlag = tk.StringVar(); self.createFlag.set('RZ')
		self.create_R1 = tk.Radiobutton(frame, text = 'RZ', variable = self.createFlag, value = 'RZ', 
			command = self.refresh_grid_labels)
		self.create_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.create_R2 = tk.Radiobutton(frame, text = 'psi_n', variable = self.createFlag, value = 'psi', 
			command = self.refresh_grid_labels)
		self.create_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Coordinate Type").grid(column = 1, row = row, sticky = tk.E )
		row += 1

		# --- phistart ---
		row = self.set_phistart_element(row)

		# --- invisible separator ---
		tk.Label(frame, text = "").grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)
		row += 1
		
		# --- Coils & Mode & Particles ---
		row = self.set_bottom_elements(row, SIESTA = False)		
		
		# --- Run ---
		runButton = tk.Button(frame, text = "Run", command = self.run_funct)
		runButton.grid(column = 1, row = row, columnspan = 5, sticky = tk.W + tk.E)
		row += 1

		# --- adjust style ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.pi_text.grid_forget()
		self.set_defaults()


	# --- Change Labels on grid variables, depending on createFlag ---
	def refresh_grid_labels(self):
		if(self.createFlag.get() == 'psi'):
			self.x_label.configure(text = "theta [rad]")
			self.y_label.configure(text = "psi_n")
			self.pi_text.grid(column = 4, row = self.pi_row, columnspan = 2)
		else:
			self.x_label.configure(text = "R [m]")
			self.y_label.configure(text = "Z [m]")
			self.pi_text.grid_forget()


	# --- turn on/off period entry ---
	def activate_entrys(self):
		data = self.read_par_file(['_fix.dat'])
		if(abs(self.HypRPt.get()) == 1):
			data = self.tool_defaults(None)	# always overwrite whats in any existing file
			self.period_entry.configure(state=tk.DISABLED)
			self.createFlag.set('RZ')
			self.refresh_grid_labels()
			self.create_R1.configure(state=tk.DISABLED)
			self.create_R2.configure(state=tk.DISABLED)
			if(self.HypRPt.get() == 1): # lower X-point
				self.ymin.set(repr(data[4]))
				self.ymax.set(repr(data[5]))
			elif(self.HypRPt.get() == -1): # upper X-point
				self.ymin.set(repr(-data[5]))
				self.ymax.set(repr(-data[4]))
			self.HypPt.set(str(1))
		else:
			if data is None: data = self.tool_defaults(None)
			self.ymin.set(repr(data[4]))
			self.ymax.set(repr(data[5]))
			self.period_entry.configure(state=tk.NORMAL)
			self.create_R1.configure(state=tk.NORMAL)
			self.create_R2.configure(state=tk.NORMAL)
		self.Nx.set(str(int(data[6]**0.5)))
		self.xmin.set(repr(data[2]))
		self.xmax.set(repr(data[3]))
		self.Ny.set(str(int(data[6]**0.5)))
		self.activate_MinMax_entrys(abs(self.HypRPt.get()) == 1)
		

	# --- define machine specific default values ---
	def tool_defaults(self, flag):	
		if self.MachFlag.get() == 'any':
			data = [1e-4, 0, 0, 0, 0, 0, 900, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 1, 100, 0.1, 2, 0,
					3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'iter':
			data = [1e-4, 0, 4, 5.85, -4.4, -3.3, 900, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 2,
					3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			data = [1e-4, 0, 0.25, 0.55, -1.5, -1, 900, 0, 1, 0, -1, 1, 5, 1, 0, 0, 0, 1, 
					100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			data = [1e-4, 0, 0.4, 0.9, -1.5, -1.0, 900, 0, 1, 0, -1, 1, 5, 1, 1, 0, 0, 1, 
					100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		else:
			data = [1e-4, 0, 1.1, 1.6, -1.45, -0.8, 900, 0, 1, 0, -1, 1, 5, 1, 1, 1, 0, 1, 
					100, 0.1, 2, 0, 0, 0, 3.141592653589793, 6.283185307179586]
		return data


	# --- set default values ---
	def set_defaults(self):
		data = self.read_par_file(['_fix.dat'])
		data = self.set_machine_defaults(data, self.MachFlag.get())						
		self.shift = data['0-8'][0]
		if(data['create'] == 3) | (data['create'] == 4): self.createFlag.set('psi')
		else: self.createFlag.set('RZ')	
		if(self.createFlag.get() == 'RZ'): 
			self.xmin.set(repr(data['0-8'][2]))
			self.xmax.set(repr(data['0-8'][3]))
			self.ymin.set(repr(data['0-8'][4]))
			self.ymax.set(repr(data['0-8'][5]))
			self.Nx.set(str(int(data['0-8'][6])))
			self.Ny.set(str(int(data['0-8'][0])))
		else: 
			self.xmin.set(repr(data['0-8'][4]))
			self.xmax.set(repr(data['0-8'][5]))
			self.ymin.set(repr(data['0-8'][2]))
			self.ymax.set(repr(data['0-8'][3]))
			self.Nx.set(str(int(data['0-8'][0])))
			self.Ny.set(str(int(data['0-8'][6])))
		self.activate_entrys()
		self.refresh_grid_labels()
		

	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		tool = self.MachFlag.get() + 'fix'
		name = '_fix.dat'
		shellCall = tool + ' ' + name + ' ' + str(int(self.HypPt.get())) + ' ' + self.tag.get()
		self.run_common_funct(name, shellCall, int(self.HypPt.get()))


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, period, tag, shellCall):
		tooltag = 'fix'
		self.write_common_qsub_file(tooltag, 1, tag, shellCall, mpi = False)


	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		N = int(self.Nx.get()) * int(self.Ny.get())
		with open(name, 'w') as f:
			self.write_headerLines(self.MachFlag.get(), f)
			f.write('shift=\t' + repr(self.shift) + '\n')
			f.write('itt=\t0\n')
			if(self.createFlag.get() == 'RZ'): 	
				f.write('Rmin=\t' + self.xmin.get() + '\n')
				f.write('Rmax=\t' + self.xmax.get() + '\n')
				f.write('Zmin=\t' + self.ymin.get() + '\n')
				f.write('Zmax=\t' + self.ymax.get() + '\n')
			else:
				f.write('psimin=\t' + self.ymin.get() + '\n')
				f.write('psimax=\t' + self.ymax.get() + '\n')
				f.write('thmin=\t' + self.xmin.get() + '\n')
				f.write('thmax=\t' + self.xmax.get() + '\n')
			f.write('N=\t' + str(N) + '\n')
			self.write_Ctrl_center(f)
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')
			if(self.createFlag.get() == 'RZ'): 
				f.write('createPoints(0=setr,4=randompsi,5=setR)=\t5\n')
			else:
				f.write('createPoints(0=setr,4=randompsi,5=setR)=\t4\n')	
			self.write_Ctrl_bottom(f)


# -------------------------------------------------------------------------------------------------------------
# --- dtman ---------------------------------------------------------------------------------------------------
class set_man_tab(common_tab):		# inherit common tab class
	def __init__(self, frame, M):
		self.frame = frame
		self.MachFlag = M.MachFlag
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.translateMachFlag = M.translateMachFlag
		self.gFile_entry = M.gFile_entry
		self.labelFont = None
		self.textFont = None

		okayCommand = self.frame.register(M.isOkay)
		notOkayCommand = self.frame.register(M.isNotOkay)

		# --- unused Variables ---
		self.itt = tk.StringVar(); self.itt.set('0')

		# --- type ---
		# define... 
		self.Type = tk.IntVar(); row = 0

		# ...and set type RadioButton
		tk.Radiobutton(frame, text = 'unstable', variable = self.Type, value = 1, 
			command = self.set_type).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(frame, text = 'stable', variable = self.Type, value = -1, 
			command = self.set_type).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Manifold").grid(column = 1, row = row, sticky = tk.E )

		# --- shift ---
		self.shift = tk.StringVar(); row += 1
		tk.Entry(frame, width = 7, textvariable = self.shift).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "|Shift|").grid(column = 1, row = row, sticky = tk.E )

		self.shift_sign = tk.IntVar()
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
		row += 1

		# --- phistart ---
		row = self.set_phistart_element(row)
		
		# --- MapDirection ---
		row = self.set_MapDirection_element(row, state = tk.DISABLED)

		# --- Coils & Mode & Particles ---
		row = self.set_bottom_elements(row, SIESTA = False)		
		
		# --- Run ---
		runButton = tk.Button(frame, text = "Run", command = self.run_funct)
		runButton.grid(column = 1, row = row, columnspan = 5, sticky = tk.W + tk.E)
		row += 1

		# --- adjust style ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.set_defaults()


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
			

	# --- define machine specific default values ---
	def tool_defaults(self, flag):	
		if self.MachFlag.get() == 'any':
			data = [1e-4, 0, 0, 0, 0, 0, 900, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 1, 100, 0.1, 2, 0,
					3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'iter':
			data = [1e-4, 0, 4, 5.85, -4.4, -3.3, 900, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 2,
					3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			data = [1e-4, 0, 0.25, 0.55, -1.5, -1, 900, 0, 1, 0, -1, 1, 0, 1, 0, 0, 0, 1, 
					100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			data = [1e-4, 0, 0.4, 0.9, -1.5, -1.0, 900, 0, 1, 0, -1, 1, 5, 1, 1, 0, 0, 1, 
					100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		else:
			data = [1e-4, 0, 1.1, 1.6, -1.45, -0.8, 900, 0, 1, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
					100, 0.1, 2, 0, 0, 0, 3.141592653589793, 6.283185307179586]
		return data


	# --- set default values ---
	def set_defaults(self):
		data = self.read_par_file(['_fix.dat'])
		data = self.set_machine_defaults(data, self.MachFlag.get())
										
		self.xmin = data['0-8'][4]
		self.xmax = data['0-8'][5]
		self.ymin = data['0-8'][2]
		self.ymax = data['0-8'][3]
		self.N = int(data['0-8'][6])
		
		if(data['0-8'][8] >= 0): self.Type.set(1)
		else: self.Type.set(-1)
			
		self.shift.set(repr(abs(data['0-8'][0])))
		
		if(data['0-8'][0] >= 0): self.shift_sign.set(1)
		else: self.shift_sign.set(-1)
		
		self.set_type()


	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		name = '_fix.dat'
		tool = self.MachFlag.get() + 'man'
		if(self.fixfile_tag.get() == ''): fixtag = ''
		else: fixtag = '_' + self.fixfile_tag.get()			
		fixfile = 'fix_' + str(int(self.fixfile_period.get())) + fixtag + '.dat'
		shellCall = tool + ' ' + name + ' ' + fixfile + ' ' + self.tag.get()
		self.run_common_funct(name, shellCall, fixfile)


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, fixfile, tag, shellCall):
		tooltag = 'M'
		self.write_common_qsub_file(tooltag, 1, tag, shellCall, mpi = False)


	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		shift = abs(float(self.shift.get())) * self.shift_sign.get()
		with open(name, 'w') as f:
			self.write_headerLines(self.MachFlag.get(), f)
			f.write('shift=\t' + repr(shift) + '\n')
			f.write('itt=\t0\n')		
			f.write('rmin=\t' + repr(self.ymin) + '\n')
			f.write('rmax=\t' + repr(self.ymax) + '\n')
			f.write('thmin=\t' + repr(self.xmin) + '\n')
			f.write('thmax=\t' + repr(self.xmax) + '\n')
			f.write('N=\t' + str(self.N) + '\n')
			self.write_Ctrl_center(f)
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')			
			f.write('createPoints(0=setr,3=setpsi,5=setR)=\t0\n')			
			self.write_Ctrl_bottom(f)


# -------------------------------------------------------------------------------------------------------------
# --- foot ----------------------------------------------------------------------------------------------------
class set_foot_tab(common_tab):		# inherit common tab class
	def __init__(self, frame, M):
		self.frame = frame
		self.MachFlag = M.MachFlag
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.translateMachFlag = M.translateMachFlag
		self.gFile_entry = M.gFile_entry
		self.labelFont = None
		self.textFont = None

		# --- unused Variables ---
		self.phistart = tk.StringVar(); self.phistart.set('0')
		
		# --- target  type ---
		row = 0
		# define... 
		self.TargetFlag = tk.IntVar()

		# ...and set grid-type RadioButton
		self.target_RB1 = tk.Radiobutton(frame, variable = self.TargetFlag, value = 1, command = self.refresh_grid_labels)
		self.target_RB2 = tk.Radiobutton(frame, variable = self.TargetFlag, value = 2, command = self.refresh_grid_labels)
		self.target_RB3 = tk.Radiobutton(frame, variable = self.TargetFlag, value = 3, command = self.refresh_grid_labels)
		self.target_RB4 = tk.Radiobutton(frame, variable = self.TargetFlag, value = 4, command = self.refresh_grid_labels)
		self.target_RB0 = tk.Radiobutton(frame, variable = self.TargetFlag, value = 0, command = self.refresh_grid_labels)
		tk.Label(frame, text = "Target").grid(column = 1, row = row, sticky = tk.W + tk.E )
		
		self.UpgradeFlag = tk.IntVar()
		self.U_chkBtn = tk.Checkbutton(frame, text = 'U', variable = self.UpgradeFlag, command = self.refresh_grid_labels)
				
		# --- x -> phi ---
		row += 1
		self.x_label = tk.Label(frame, text = "phi [rad]")
		self.x_label.grid(column = 1, row = row, sticky = tk.E)
		
		self.set_MinMax_elements(row)
				
		# --- y -> t ---
		row += 2
		self.y_label = tk.Label(frame, text = "")
		self.y_label.grid(column = 1, row = row, sticky = tk.E )
		
		row += 2
			
		# --- Info Text ---
		self.Info = tk.Label(frame, text = "")
		self.Info.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)
		row += 1
		
		# --- toroidal turns ---
		row = self.set_toroidatTurn_element(row)
		
		# --- MapDirection ---
		row = self.set_MapDirection_element(row, state = tk.NORMAL)

		# --- Coils & Mode & Particles ---
		row = self.set_bottom_elements(row, SIESTA = False)		
		
		# --- Run ---
		self.set_run_mpi_elements(row)
		runButton = tk.Button(frame, text = "Run", command = self.run_funct)
		runButton.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)

		# --- adjust style ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.set_defaults()


	# --- Change first row (row = 0), depending on Machine ---
	def set_target_elements(self):
		row = 0
		if self.MachFlag.get() == 'any':
			self.target_RB0.configure(text = 'Full Wall')		
			self.target_RB0.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB1.grid_forget()
			self.target_RB2.grid_forget()
			self.target_RB3.grid_forget()
			self.target_RB4.grid_forget()
			self.U_chkBtn.grid_forget()
		elif self.MachFlag.get() == 'iter':
			self.target_RB0.configure(text = 'Full Wall')		
			self.target_RB1.configure(text = 'Inner')
			self.target_RB2.configure(text = 'Outer')
			self.target_RB0.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB1.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 4, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid_forget()
			self.target_RB4.grid_forget()
			self.U_chkBtn.grid_forget()
		elif self.MachFlag.get() == 'nstx':
			self.target_RB0.configure(text = 'Full Wall')
			self.target_RB1.configure(text = 'Inn-up')
			self.target_RB2.configure(text = 'Out-up')
			self.target_RB3.configure(text = 'Inn-dwn')
			self.target_RB4.configure(text = 'Out-dwn')
			self.target_RB0.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB1.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 4, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid_forget()
			self.target_RB4.grid(column = 5, row = row, sticky = tk.W + tk.E )
			self.U_chkBtn.grid(column = 1, row = row, sticky = tk.E)
		elif self.MachFlag.get() == 'mast':
			self.target_RB0.configure(text = 'Full Wall')		
			self.target_RB1.configure(text = 'Inner')
			self.target_RB2.configure(text = 'Outer')
			self.target_RB0.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB1.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 4, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid_forget()
			self.target_RB4.grid_forget()
			self.U_chkBtn.grid_forget()
		else:
			self.target_RB0.configure(text = 'Full Wall')		
			self.target_RB1.configure(text = 'Inner')
			self.target_RB2.configure(text = 'Outer')
			self.target_RB3.configure(text = 'Shelf')
			self.target_RB0.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB1.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 4, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid(column = 5, row = row, sticky = tk.W + tk.E )
			self.target_RB4.grid_forget()
			self.U_chkBtn.grid_forget()
		

	# --- Change Labels on grid variables, depending on TargetFlag ---
	def refresh_grid_labels(self):
		if self.MachFlag.get() == 'any':
			self.y_label.configure(text = "Swall [m]")
			self.Info.configure(text = "s = 0: Centerpost midplane, s > 0: counter clock wise")
		elif self.MachFlag.get() == 'iter':
			if(self.TargetFlag.get() == 1):
				self.y_label.configure(text = "t [cm]")
				self.Info.configure(text = "t in [-166 <--> 31] is length along the wall; t < 0: curve upwards, t > 0: linear outwards")
			elif(self.TargetFlag.get() == 2):
				self.y_label.configure(text = "t [cm]")
				self.Info.configure(text = "t in [-39 <--> 163] is length along the wall; t < 0: linear inwards, t > 0: curve upwards")
			else:
				self.y_label.configure(text = "Swall [m]")
				self.Info.configure(text = "s = 0: Centerpost midplane, s > 0: counter clock wise")
		elif self.MachFlag.get() == 'nstx':
			if(self.TargetFlag.get() == 1):
				self.y_label.configure(text = "Z [m]")
				if self.UpgradeFlag.get() == 1:
					self.Info.configure(text = "NSTX-U inner target, upper divertor, Z in [1.05 <--> 1.578]")
				else:
					self.Info.configure(text = "inner target, upper divertor, Z in [1.1714 <--> 1.578]")
			elif(self.TargetFlag.get() == 2):
				self.y_label.configure(text = "R [m]")
				if self.UpgradeFlag.get() == 1:
					self.Info.configure(text = "NSTX-U outer target, upper divertor, R in [0.435 <--> 1.0433]")
				else:
					self.Info.configure(text = "outer target, upper divertor, R in [0.2979 <--> 1.0433]")
			elif(self.TargetFlag.get() == 3):
				self.y_label.configure(text = "Z [m]")
				if self.UpgradeFlag.get() == 1:
					self.Info.configure(text = "NSTX-U inner target, lower divertor, Z in [-1.578 <--> -1.05]")
				else:
					self.Info.configure(text = "inner target, lower divertor, Z in [-1.578 <--> -1.1714]")
			elif(self.TargetFlag.get() == 4):
				self.y_label.configure(text = "R [m]")
				if self.UpgradeFlag.get() == 1:
					self.Info.configure(text = "NSTX-U outer target, lower divertor, R in [0.435 <--> 1.0433]")
				else:
					self.Info.configure(text = "outer target, lower divertor, R in [0.2979 <--> 1.0433]")
			else:
				self.y_label.configure(text = "Swall [m]")
				self.Info.configure(text = "s = 0: Centerpost midplane, s > 0: counter clock wise")
		elif self.MachFlag.get() == 'mast':
			if(self.TargetFlag.get() == 1):
				self.y_label.configure(text = "Z [m]")
				self.Info.configure(text = "inner target, Z in [-1.6835 <--> -1.22909]")
			elif(self.TargetFlag.get() == 2):
				self.y_label.configure(text = "R [m]")
				self.Info.configure(text = "outer target, R in [0.7835 <--> 1.9]")
			else:
				self.y_label.configure(text = "Swall [m]")
				self.Info.configure(text = "s = 0: Centerpost midplane, s > 0: counter clock wise")
		else:	# d3d
			if(self.TargetFlag.get() == 3):
				self.y_label.configure(text = "t [0 <--> 1]")
				self.Info.configure(text = "t > 0: Shelf outwards, t = 0: Nose edge")
			elif(self.TargetFlag.get() == 2):
				self.y_label.configure(text = "t [0 <--> 1]")
				self.Info.configure(text = "t > 0: Divertor Floor outwards, t = 0: Connection 45deg Tile, t = 1: Pump Entry")
			elif(self.TargetFlag.get() == 1):
				self.y_label.configure(text = "t [-1 <--> 1]")
				self.Info.configure(text = "t < 0: Centerpost upwards, t > 0: 45deg Tile downwards, t = 0: Connection point")
			else:
				self.y_label.configure(text = "Swall [m]")
				self.Info.configure(text = "s = 0: Centerpost midplane, s > 0: counter clock wise")

		name = self.set_name()
		data = self.read_par_file([name])
		if data is None: data = self.tool_defaults(self.TargetFlag.get())
		self.xmin.set(repr(data[4]))
		self.xmax.set(repr(data[5]))
		self.Nx.set(str(int(data[0])))
		self.ymin.set(repr(data[2]))
		self.ymax.set(repr(data[3]))
		self.Ny.set(str(int(data[6])))
		self.MapDirection.set(int(data[8]))		# same index for all machines !!!


	# --- define machine specific default values ---
	def tool_defaults(self, flag):	
		if self.MachFlag.get() == 'any':
			data = [50, 300, 10, 80, 0, 6.283185307179586, 90, 0, -1, 0, -1, 0, 2, 0, 0, 0, 0, 1, 100, 0.1, 2, 0,
						3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'iter':
			if(flag == 1):
				data = [500, 300, -60, 10, 0, 6.283185307179586, 700, 0, 1, 0, -1, 1, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
						3.141592653589793, 6.283185307179586]
			elif(flag == 2):
				data = [500, 300, -10, 80, 0, 6.283185307179586, 900, 0, -1, 0, -1, 2, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
						3.141592653589793, 6.283185307179586]
			else:
				data = [500, 300, 10, 80, 0, 6.283185307179586, 900, 0, -1, 0, -1, 2, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
						3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			if(flag == 1):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, 1.05, 1.578, 0, 6.283185307179586, 400, 0, 1, 0, -1, 1, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, 1.1714, 1.578, 0, 6.283185307179586, 400, 0, 1, 0, -1, 1, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
			elif(flag == 2):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, 0.435, 1.0433, 0, 6.283185307179586, 400, 0, -1, 0, -1, 2, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, 0.2979, 1.0433, 0, 6.283185307179586, 400, 0, -1, 0, -1, 2, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
			elif(flag == 3):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, -1.578, -1.05, 0, 6.283185307179586, 400, 0, -1, 0, -1, 3, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, -1.578, -1.1714, 0, 6.283185307179586, 400, 0, -1, 0, -1, 3, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
			elif(flag == 4):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, 0.435, 1.0433, 0, 6.283185307179586, 400, 0, 1, 0, -1, 4, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, 0.2979, 1.0433, 0, 6.283185307179586, 400, 0, 1, 0, -1, 4, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
							3.141592653589793, 6.283185307179586]
			else:
				data = [500, 500, 0.3, 1.0, 0, 6.283185307179586, 400, 0, 1, 0, -1, 4, 2, 1, 0, 0, 0, 1, 100, 0.1, 2,
						3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			if(flag == 1):
				data = [500, 500, -1.6835, -1.229, 0, 6.283185307179586, 400, 0, -1, 0, -1, 1, 2, 1, 1, 0, 0, 1, 
							100, 0.1, 2, 3.141592653589793, 6.283185307179586]
			elif(flag == 2):
				data = [500, 500, 0.7835, 1.9, 0, 6.283185307179586, 400, 0, 1, 0, -1, 2, 2, 1, 1, 0, 0, 1, 
							100, 0.1, 2, 3.141592653589793, 6.283185307179586]
			else:
				data = [500, 500, 0.7835, 1.9, 0, 6.283185307179586, 400, 0, 1, 0, -1, 2, 2, 1, 1, 0, 0, 1, 
							100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		else:	# d3d
			if(flag == 3):
				data = [500, 300, 0.0, 0.1, 0, 6.283185307179586, 100, 0, 1, 0, -1, 3, 2, 1, 1, 1, 0, 1, 
							100, 0.1, 2, 0, 0, 0, 3.141592653589793, 6.283185307179586]		
			elif(flag == 2):
				data = [500, 300, 0.9, 1.0, 0, 6.283185307179586, 100, 0, 1, 0, -1, 2, 2, 1, 1, 1, 0, 1, 
							100, 0.1, 2, 0, 0, 0, 3.141592653589793, 6.283185307179586]
			elif(flag == 1):
				data = [500, 300, -0.1, 0.3, 0, 6.283185307179586, 400, 0, -1, 0, -1, 1, 2, 1, 1, 1, 0, 1, 
							100, 0.1, 2, 0 ,0 ,0, 3.141592653589793, 6.283185307179586]
			else:
				data = [500, 300, 0.0, 0.1, 0, 6.283185307179586, 400, 0, -1, 0, -1, 1, 2, 1, 1, 1, 0, 1, 
							100, 0.1, 2, 0 ,0 ,0, 3.141592653589793, 6.283185307179586]
		return data


	# --- set default values ---
	def set_defaults(self):
		if self.MachFlag.get() == 'any':
			data = self.read_par_file(['_foot.dat'])
		elif self.MachFlag.get() == 'iter':
			data = self.read_par_file(['_foot.dat','_inner.dat','_outer.dat'])
		elif self.MachFlag.get() == 'nstx':
			data = self.read_par_file(['_foot.dat','_innerup.dat','_outerup.dat','_innerdwn.dat','_outerdwn.dat'])
		elif self.MachFlag.get() == 'mast':
			data = self.read_par_file(['_foot.dat','_inner.dat','_outer.dat'])
		else:
			data = self.read_par_file(['_foot.dat','_inner.dat','_outer.dat', '_shelf.dat'])
		data = self.set_machine_defaults(data, self.MachFlag.get(), flag = 1)
		if self.MachFlag.get() == 'any': self.TargetFlag.set(0)
		else: self.TargetFlag.set(int(data['target']))
		self.UpgradeFlag.set(0)
		self.refresh_grid_labels()


	# --- set control file name ---
	def set_name(self):
		if(self.TargetFlag.get() == 0): name = '_foot.dat'
		if self.MachFlag.get() == 'iter':
			if(self.TargetFlag.get() == 2): name = '_outer.dat'
			elif(self.TargetFlag.get() == 1): name = '_inner.dat'
		elif self.MachFlag.get() == 'nstx':
			if(self.TargetFlag.get() == 4): name = '_outerdwn.dat'
			elif(self.TargetFlag.get() == 3): name = '_innerdwn.dat'
			elif(self.TargetFlag.get() == 2): name = '_outerup.dat'
			elif(self.TargetFlag.get() == 1): name = '_innerup.dat'
		elif self.MachFlag.get() == 'mast':
			if(self.TargetFlag.get() == 2): name = '_outer.dat'
			elif(self.TargetFlag.get() == 1): name = '_inner.dat'
		else:
			if(self.TargetFlag.get() == 3): name = '_shelf.dat'
			elif(self.TargetFlag.get() == 2): name = '_outer.dat'
			elif(self.TargetFlag.get() == 1): name = '_inner.dat'
		return name
		

	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		name = self.set_name()
		tool = self.MachFlag.get() + 'foot_mpi'
		shellCall = MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' ' + tool + ' ' + name + ' ' + self.tag.get()
		self.run_common_funct(name, shellCall, int(self.nproc.get()))


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, shellCall):
		tooltag = 'F'
		self.write_common_qsub_file(tooltag, nproc, tag, shellCall, mpi = True)


	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		with open(name, 'w') as f:
			self.write_headerLines(self.MachFlag.get(), f)
			f.write('Nphi=\t' + self.Nx.get() + '\n')
			f.write('itt=\t' + self.itt.get() + '\n')	
			f.write('tmin=\t' + self.ymin.get() + '\n')
			f.write('tmax=\t' + self.ymax.get() + '\n')
			f.write('phimin=\t' + self.xmin.get() + '\n')
			f.write('phimax=\t' + self.xmax.get() + '\n')
			f.write('Nt=\t' + self.Ny.get() + '\n')
			self.write_Ctrl_center(f)
			f.write('target(0=wall,1=inner,2=outer,3=shelf)=\t' + str(self.TargetFlag.get()) + '\n')			
			f.write('createPoints(2=target)=\t2\n')	
			self.write_Ctrl_bottom(f)
		

# -------------------------------------------------------------------------------------------------------------
# --- lam -----------------------------------------------------------------------------------------------------
class set_lam_tab(common_tab):		# inherit common tab class
	def __init__(self, frame, M):
		self.frame = frame
		self.MachFlag = M.MachFlag
		self.Shot = M.Shot
		self.Time = M.Time
		self.gPath = M.gPath
		self.tag = M.tag
		self.path = M.path
		self.translateMachFlag = M.translateMachFlag
		self.gFile_entry = M.gFile_entry
		self.labelFont = None
		self.textFont = None
		
		fileFound = self.frame.register(self.file_Found)
		
		row = 0
		# --- grid  type ---
		# define... 
		self.createFlag = tk.StringVar()

		# ...and set grid-type RadioButton
		self.create_R1 = tk.Radiobutton(frame, text = 'RZ', variable = self.createFlag, value = 'RZ', 
			command = self.refresh_grid_labels)
		self.create_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.create_R2 = tk.Radiobutton(frame, text = 'psi_n', variable = self.createFlag, value = 'psi', 
			command = self.refresh_grid_labels)
		self.create_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.create_R3 = tk.Radiobutton(frame, text = 'manual', variable = self.createFlag, value = 'manual', 
			command = self.refresh_grid_labels)
		self.create_R3.grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(frame, text = "Coordinate Type").grid(column = 1, row = row, sticky = tk.E )

		# --- x -> theta or R ---
		row += 1
		self.x_label = tk.Label(frame, text = "")
		self.x_label.grid(column = 1, row = row, sticky = tk.E)

		self.set_MinMax_elements(row)
		
		self.pointsFile = tk.StringVar()
		self.pointsFile_entry = AutocompleteEntry(os.listdir(self.path.get()), frame, listboxLength = 6, 
				width = 50, textvariable = self.pointsFile)
		self.pointsFile_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
				
		# --- y -> r, psi or Z ---
		row += 2
		self.y_label = tk.Label(frame, text = "")
		self.y_label.grid(column = 1, row = row, sticky = tk.E )
		
		row += 1

		# --- spare interior ---
		self.spareInterior = tk.IntVar()
		tk.Checkbutton(frame, text = ' only for psi > ', variable = self.spareInterior,
			command = self.activate_spareInterior).grid(column = 4, row = row, columnspan = 2, sticky = tk.W)

		self.spareInterior_value = tk.StringVar()
		self.spareInterior_entry = tk.Entry(frame, width = 6, textvariable = self.spareInterior_value)
		self.spareInterior_entry.grid(column = 5, row = row, sticky = tk.W)
		
		row += 1

		# --- toroidal turns ---
		row = self.set_toroidatTurn_element(row)
		
		# --- phistart ---
		row = self.set_phistart_element(row)
		
		# --- MapDirection ---
		row = self.set_MapDirection_element(row)

		# --- Coils & Mode & Particles ---
		row = self.set_bottom_elements(row)		
		
		# --- Run ---
		self.set_run_mpi_elements(row)
		runButton = tk.Button(frame, text = "Run", command = self.run_funct)
		runButton.grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)

		# --- adjust style ---
		for child in frame.winfo_children(): child.grid_configure(padx=5, pady=5)
		self.set_defaults()


	# --- Change Labels on grid variables, depending on createFlag ---
	def refresh_grid_labels(self):
		name = self.set_name()
		data = self.read_par_file([name])
		if data is None: 
			data = self.tool_defaults(self.createFlag.get())
			use_defaults = True
		else:
			use_defaults = False
		if(self.createFlag.get() == 'psi'):
			if(self.selectField.get() == -2):
				self.x_label.configure(text = "u [rad]")
				self.y_label.configure(text = "s")
				if use_defaults: data[2] = 0.1; data[3] = 1.0
			else:
				self.x_label.configure(text = "theta [rad]")
				self.y_label.configure(text = "psi_n")
			self.pi_text.grid(column = 4, row = self.pi_row, columnspan = 2)
			self.xmin.set(repr(data[4]))
			self.xmax.set(repr(data[5]))
			self.ymin.set(repr(data[2]))
			self.ymax.set(repr(data[3]))
			self.Nx.set(str(int(data[0])))
			self.Ny.set(str(int(data[6])))
		else:
			self.x_label.configure(text = "R [m]")
			self.y_label.configure(text = "Z [m]")
			self.pi_text.grid_forget()
			self.xmin.set(repr(data[2]))
			self.xmax.set(repr(data[3]))
			self.ymin.set(repr(data[4]))
			self.ymax.set(repr(data[5]))
			self.Nx.set(str(int(data[6])))
			self.Ny.set(str(int(data[0])))
		self.activate_MinMax_entrys(self.createFlag.get() == 'manual')
		if(self.createFlag.get() == 'manual'):
			self.x_label.configure(text = "Points data: ")
			self.pointsFile_entry.grid(column = 2, row = 1, columnspan = 4,sticky = tk.E+tk.W)
		else:
			self.pointsFile_entry.grid_forget()


	def activate_spareInterior(self):
		if self.spareInterior.get() == 1:
			self.spareInterior_entry.configure(state = tk.NORMAL)
		else:
			self.spareInterior_entry.configure(state = tk.DISABLED)
		

	# --- define machine specific default values ---
	def tool_defaults(self, flag):	
		if self.MachFlag.get() == 'any':
			if(flag == 'psi'):
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 0, -1, 1, 3, 0, 0, 0, 0, 1, 100, 0.1, 2, 0,
						3.141592653589793, 6.283185307179586]
			else:
				data = [520, 200, 0, 0, 0, 0, 500, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, 100, 0.1, 2, 0,
						3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'iter':
			if(flag == 'psi'):
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 0, -1, 1, 3, 1, 0, 0, 0, 1, 100, 0.1, 2, 
						3.141592653589793, 6.283185307179586]
			else:
				data = [520, 200, 4.0, 6.5, -4.6, -2.0, 500, 0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 2,
						3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			if(flag == 'psi'):
				data = [300, 500, 0.88, 1.02, 0, 6.283185307179586, 150, 0, 0, 0, -1, 1, 3, 1, 0, 0, 0, 1, 
						100, 0.1, 2, 3.141592653589793, 6.283185307179586]
			else:
				data = [100, 500, 0.17, 0.9, -1.65, -1.0, 100, 0, 0, 0, -1, 1, 0, 1, 0, 0, 0, 1, 
						100, 0.1, 2, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			if(flag == 'psi'):
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 0, -1, 1, 3, 1, 1, 0, 0, 1, 
						100, 0.1, 2, 3.141592653589793, 6.283185307179586]
			else:
				data = [610, 200, 0.195, 1.9, -1.825, 0.6, 850, 0, 0, 0, -1, 1, 0, 1, 1, 0, 0, 1, 
						100, 0.1, 3.141592653589793, 6.283185307179586]
		else:
			if(flag == 'psi'):
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 0, -1, 1, 3, 1, 1, 1, 0, 1, 
						100, 0.1, 2, 0, 0, 0, 3.141592653589793, 6.283185307179586]
			else:
				data = [930, 200, 1.0, 1.45, -1.367, -0.902, 900, 0, 0, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
						100, 0.1, 2, 0, 0, 0, 3.141592653589793, 6.283185307179586]
		return data


	# --- set default values ---
	def set_defaults(self):
		data = self.read_par_file(['_lam_psi.dat','_lam.dat'])
		data = self.set_machine_defaults(data, self.MachFlag.get(), flag = 'RZ')
		if(data['create'] == 3) | (data['create'] == 4): self.createFlag.set('psi')
		else: self.createFlag.set('RZ')
		self.spareInterior.set(0)
		self.spareInterior_value.set('0.85')
		self.activate_spareInterior()
		self.refresh_grid_labels()
		self.pointsFile_entry.update_baseSearchPath(self.path.get())

		
	# --- set control file name ---
	def set_name(self):
		if(self.createFlag.get() == 'psi'):
			name = '_lam_psi.dat'
		else:
			name = '_lam.dat'
		return name


	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		name = self.set_name()
		tool = self.MachFlag.get() + 'laminar_mpi'
		shellCall = MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' ' + tool + ' ' + name + ' ' + self.tag.get()
		shellFlags = ''
		if(self.createFlag.get() == 'manual'):
			pointsfile = self.pointsFile.get()
			if self.file_Found(file = pointsfile):
				shellFlags += ' -P ' + pointsfile
		if self.spareInterior.get() == 1:
			try: spareValue = np.float64(self.spareInterior_value.get())
			except: spareValue = -1
			if (spareValue > 0) & (spareValue < 1):
				shellFlags += ' -s -l ' + str(spareValue)
			else: print ('invalid entry for psi limit')
		self.run_common_funct(name, shellCall, int(self.nproc.get()), shellFlags)


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, shellCall):
		tooltag = 'L'
		self.write_common_qsub_file(tooltag, nproc, tag, shellCall, mpi = True)


	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		with open(name, 'w') as f:
			self.write_headerLines(self.MachFlag.get(), f)			
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
			self.write_Ctrl_center(f)
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')		
			if(self.createFlag.get() == 'psi'): f.write('createPoints(0=setR,3=setpsi)=\t3\n')
			else: f.write('createPoints(0=setR,3=setpsi)=\t0\n')
			self.write_Ctrl_bottom(f)


# -------------------------------------------------------------------------------------------------------------
# --- AutocompleteEntry ---------------------------------------------------------------------------------------
class AutocompleteEntry(tk.Entry):
	def __init__(self, autocompleteList, *args, **kwargs):
		# set working dir
		self.init_path = kwargs.pop('cwd', os.getcwd() + '/')
		self.path = self.init_path
		self.include_path = False

		# Listbox length
		self.listboxLength = kwargs.pop('listboxLength', 8)
		
		# Entry variable
		self.var = kwargs.pop('textvariable', tk.StringVar())

		# Custom matches function
		if 'matchesFunction' in kwargs:
			self.matchesFunction = kwargs['matchesFunction']
			del kwargs['matchesFunction']
		else:
			def matches(fieldValue, acListEntry):
				pattern = re.compile('.*' + re.escape(fieldValue) + '.*', re.IGNORECASE)
				return re.match(pattern, acListEntry)
				
			self.matchesFunction = matches
		
		tk.Entry.__init__(self, *args, textvariable = self.var, **kwargs)
		self.focus()

		self.autocompleteList = autocompleteList
		
		self.var.trace('w', self.changed)
		self.bind("<Right>", self.selection)
		self.bind("<Return>", self.selection)
		self.bind("<Up>", self.moveUp)
		self.bind("<Down>", self.moveDown)
		
		self.listboxUp = False
		
		
	def separate_dir(self):
		var = self.var.get()
		idx = var[::-1].find('/')	# returns location of last '/' in var or -1 if not found
		if(idx == -1):
			self.path = '.'
		else:
			idx *= -1
			self.path = var[0:idx - 1] + '/'	# path with a final '/', works also if var[0:idx - 1] == '', which means var is like '/blabla'
			var = var[idx::]
		if not self.path[0] == '/':	# path is not absolute
			self.path = os.path.abspath(self.init_path + self.path) + '/'
		return var


	def update_baseSearchPath(self, path):
		if os.path.isfile(path):
			self.init_path = os.path.dirname(os.path.abspath(path)) + '/'	# path is an existing filepath
		elif os.path.isdir(path):
			self.init_path = os.path.abspath(path) + '/'	# path is an existing path
		else:
			self.init_path = os.path.abspath('./') + '/'	# path does not exist
		self.path = self.init_path
		self.update_list()
		
		
	def update_list(self, newList = None):
		if newList is None:
			self.autocompleteList = os.listdir(self.path)
		else:
			self.autocompleteList = newList


	def changed(self, name, index, mode):
		if self.var.get() == '':
			if self.include_path:
				self.path = self.init_path
				self.update_list()
				self.include_path = False
			if self.listboxUp:
				self.listbox.destroy()
				self.tw.destroy()
				self.listboxUp = False
		else:
			words = self.comparison()
			if words:
				if not self.listboxUp:
					# creates a toplevel window
					self.tw = tk.Toplevel(self)
					x = y = 0
					x, y, cx, cy = self.bbox("insert")
					x += self.winfo_rootx()
					y += self.winfo_rooty() + self.winfo_height()
					# Leaves only the label and removes the app window
					self.tw.wm_overrideredirect(True)
					self.tw.wm_geometry("+%d+%d" % (x, y))
					self.listbox = tk.Listbox(self.tw, width=self["width"], height=self.listboxLength)
					self.listbox.bind("<Button-1>", self.selection)
					self.listbox.bind("<Right>", self.selection)
					self.listbox.bind("<Return>", self.selection)
					self.listbox.pack(ipadx = 1)
					self.listboxUp = True
				
				self.listbox.delete(0, tk.END)
				for w in words:
					self.listbox.insert(tk.END,w)
				self.listbox.selection_set(first=0)
				self.listbox.activate(0)
			else:
				if self.listboxUp:
					self.listbox.destroy()
					self.tw.destroy()
					self.listboxUp = False

		
	def selection(self, event):
		if self.listboxUp:
			if event is None: pass
			elif event.type == '4':		# button press event
				index = self.listbox.index("@%s,%s" % (event.x, event.y))
				self.listbox.selection_set(first=index)
				self.listbox.activate(index)
			if self.include_path:
				self.var.set(self.path + self.listbox.get(tk.ACTIVE))
			else:
				self.var.set(self.listbox.get(tk.ACTIVE))
			self.listbox.destroy()
			self.tw.destroy()
			self.listboxUp = False
			self.icursor(tk.END)


	def moveUp(self, event):
		if self.listboxUp:
			if self.listbox.curselection() == ():
				index = '0'
			else:
				index = self.listbox.curselection()[0]
				
			if index != '0':				
				self.listbox.selection_clear(first=index)
				index = str(int(index) - 1)
				
				self.listbox.see(index) # Scroll!
				self.listbox.selection_set(first=index)
				self.listbox.activate(index)


	def moveDown(self, event):
		if self.listboxUp:
			if self.listbox.curselection() == ():
				index = '0'
			else:
				index = self.listbox.curselection()[0]
				
			if index != tk.END:						
				self.listbox.selection_clear(first=index)
				index = str(int(index) + 1)
				
				self.listbox.see(index) # Scroll!
				self.listbox.selection_set(first=index)
				self.listbox.activate(index) 


	def comparison(self):
		var = self.var.get()
		if '/' in var:
			var = self.separate_dir()
			self.update_list()
			self.include_path = True
		if self.var.get()[-1] == '/':
			return [ w for w in self.autocompleteList]
		else:
			return [ w for w in self.autocompleteList if self.matchesFunction(var, w) ]


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
		'MAFOT Version 5.7 \n'
		'GUI Version 3.21 \n'
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
# --- single functions ----------------------------------------------------------------------------------------

# --- read Control File, and return its values ---
def readControlFile(name, default_gpath):
	with open(name, 'r') as f:
		# Skip first line
		head1 = f.readline().split()
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
			gpath = default_gpath
			data = [float(head[-1])]
	
		# read data to end of file
		for line in f:
			line = line.split()			
			data.append(float(line[-1]))
		
	return shot, time, gpath, data, head1

	
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
		LD_LIBRARY_PATH += ':/opt/intel/composerxe-2011.0.084/mkl/lib/intel64'
		LD_LIBRARY_PATH += ':/home/wingen/lib/64/blitz/lib:/home/wingen/lib/64:/home/wingen/lib/64/hdf5/lib:/home/wingen/lib/64/m3dc1/lib'
		os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
		if not os.path.exists(HOME + '/work'):
			os.makedirs(HOME + '/work')

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


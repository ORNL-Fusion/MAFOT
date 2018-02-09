import os
import Tkinter as tk
import ttk
import socket
import numpy as np
from subprocess import call

import autocompleteEntry as acE 

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
		self.translateMachFlag = {'dt':'DIII-D', 'iter':'ITER', 'nstx':'NSTX', 
									'mast':'MAST', 'cmod':'C-Mod'}
	
		okayCommand = self.frame.register(self.isOkay)
		notOkayCommand = self.frame.register(self.isNotOkay)
		makecwd = self.frame.register(self.makeCWD)
		
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
		tk.Radiobutton(self.frame, text = 'C-Mod', variable = self.MachFlag, value = 'cmod', 
			command = self.set_Machine, state = tk.DISABLED).grid(column = 6, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "Machine").grid(column = 1, row = row, sticky = tk.E )

		# --- Shot number ---
		row += 1
		tk.Entry(self.frame, width = 7, textvariable = self.Shot).grid(column = 2, row = row, sticky = tk.W + tk.E)
		tk.Label(self.frame, text = "Shot").grid(column = 1, row = row, sticky = tk.E)

		# --- Time ---
		tk.Entry(self.frame, width = 7, textvariable = self.Time).grid(column = 4, row = row, sticky = tk.W + tk.E)
		tk.Label(self.frame, text = "Time").grid(column = 3, row = row, sticky = tk.E)

		# --- gPath ---
		row += 1
		tk.Entry(self.frame, width = 7, textvariable = self.gPath).grid(column = 2, row = row, columnspan = 4, sticky = tk.W + tk.E)
		tk.Label(self.frame, text = "Path to g-file").grid(column = 1, row = row, sticky = tk.E)
		
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
		default_gpath = HOME + '/c++/' + self.MachFlag.get() + '/gfiles'
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
			if self.MachFlag.get() == 'dt':
				gpath = HOME + '/c++/d3d/gfiles'
			else:
				gpath = HOME + '/c++/' + self.MachFlag.get() + '/gfiles'
			self.ControlFileFound = False
			
		if not (gpath[-1] == '/'): gpath += '/'
		self.Shot.set(str(shot)); 
		self.Time.set(str(time))
		self.gPath.set(gpath)

	
	# --- validity check of Tag ---
	# returns True if tag is alphanumeric or has _ + - as special chars
	# else returns False
	def isOkay(self, tag):
		tag = tag.translate(None, '_+- ')
		if(len(tag) == 0): return True
		else: return tag.isalnum()
		
	# prints error Message, if isOkay == False
	def isNotOkay(self):
		print 'Warning: Invalid Input Character in File Tag. Only alphanumeric and + - _ are allowed'
		
	
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
		if(self.MachFlag.get() == 'iter'):
			if not self.ControlFileFound:
				self.gPath.set(HOME + '/c++/iter/gfiles/')
			self.nb.tab(0, text='iterplot')
			self.nb.tab(1, text='iterfix')
			self.nb.tab(2, text='iterman')
			self.nb.tab(3, text='iterfoot')
			self.nb.tab(4, text='iterlaminar')
		elif(self.MachFlag.get() == 'nstx'):
			if not self.ControlFileFound:
				self.gPath.set(HOME + '/c++/nstx/gfiles/')
			self.nb.tab(0, text='nstxplot')
			self.nb.tab(1, text='nstxfix')
			self.nb.tab(2, text='nstxman')
			self.nb.tab(3, text='nstxfoot')
			self.nb.tab(4, text='nstxlaminar')
		elif(self.MachFlag.get() == 'mast'):
			if not self.ControlFileFound:
				self.gPath.set(HOME + '/c++/mast/gfiles/')
			self.nb.tab(0, text='mastplot')
			self.nb.tab(1, text='mastfix')
			self.nb.tab(2, text='mastman')
			self.nb.tab(3, text='mastfoot')
			self.nb.tab(4, text='mastlaminar')
		else:	# d3d
			if not self.ControlFileFound:
				self.gPath.set(HOME + '/c++/d3d/gfiles/')
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


	def define_machine_elements(self):
		self.useIcoil = tk.IntVar()
		self.useCcoil = tk.IntVar()
		self.useFcoil = tk.IntVar()
		self.useBus = tk.IntVar()
		self.useBcoil = tk.IntVar()

		# --- 3D Fields Label ---
		self.machine_3D_Label = tk.Label(self.frame)
		# --- Icoil checkbox ---
		self.machine_3D_chk1 = tk.Checkbutton(self.frame, variable = self.useIcoil)
		# --- Ccoil checkbox ---
		self.machine_3D_chk2 = tk.Checkbutton(self.frame, variable = self.useCcoil)
		# --- Error Fields Label ---
		self.machine_ef_Label = tk.Label(self.frame)
		# --- useEF checkbox ---
		self.machine_ef_chk1 = tk.Checkbutton(self.frame, variable = self.useFcoil)
		# --- useBus checkbox ---
		self.machine_ef_chk2 = tk.Checkbutton(self.frame, variable = self.useBus)
		# --- Bcoil checkbox ---
		self.machine_ef_chk3 = tk.Checkbutton(self.frame, variable = self.useBcoil)


	def set_machine_elements(self, MachFlag, row = 10):
		try: self.forget_elements()
		except: pass
		if(MachFlag == 'iter'):
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
		if(MachFlag == 'iter'):
			data = self.set_iter_defaults(data, flag)
		elif(MachFlag == 'nstx'):
			data = self.set_nstx_defaults(data, flag)
		elif(MachFlag == 'mast'):
			data = self.set_mast_defaults(data, flag)
		else:	# d3d
			data = self.set_d3d_defaults(data, flag)
		return data


	def write_coils(self, MachFlag, f):
		if(MachFlag == 'iter'):
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

	
	def set_d3d_defaults(self, data, flag = None):
		if data is None:
			data = self.tool_defaults(flag)
		elif len(data) < 26:
			if len(data) == 23: data = np.append(data[0:21],[0, 0, 0, 3.141592653589793, 6.283185307179586])
			elif len(data) == 24: data = np.append(data[0:22],[0, 0, 3.141592653589793, 6.283185307179586])
			elif len(data) == 25: data = np.append(data[0:23],[0, 3.141592653589793, 6.283185307179586])

		data_dic = {'0-8':data[0:9],'itt':str(int(data[1])),'phistart':repr(data[7]),'MapDirection':int(data[8]),
		'useFcoil':int(data[13]),'useCcoil':int(data[14]),'useIcoil':int(data[15]),'useBus':int(data[22]),'useBcoil':int(data[23]),
		'useM3DC1':int(data[10]),'response':int(data[9]),'selectField':int(data[10]),
		'sigma':int(data[16]),'charge':int(data[17]),'Ekin':repr(data[18]),'Lambda':repr(data[19]),
		'useFilament':str(int(data[20]))}
		
		data_dic['target'] = data[11]
		data_dic['create'] = data[12]
		
		self.useFcoil.set(data_dic['useFcoil'])
		self.useCcoil.set(data_dic['useCcoil'])
		self.useIcoil.set(data_dic['useIcoil'])
		self.useBus.set(data_dic['useBus'])
		self.useBcoil.set(data_dic['useBcoil'])
		
		self.set_common_defaults(data_dic)
		return data_dic
		
	
	def write_d3d_coils(self,f):
		f.write('useFcoil(0=no,1=yes)=\t' + str(self.useFcoil.get()) + '\n')
		f.write('useCcoil(0=no,1=yes)=\t' + str(self.useCcoil.get()) + '\n')
		f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
	
	
	def write_d3d_errorFileds(self,f):
		f.write('useBusError(0=no,1=yes)=\t' + str(self.useBus.get()) + '\n')
		f.write('useBcoilError(0=no,1=yes)=\t' + str(self.useBcoil.get()) + '\n')
		

	# --- ITER ---------------------------------------------------------------------------
	def set_iter_elements(self, row):
		self.machine_3D_Label.configure(text = "3-D Fields: ")
		self.machine_3D_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		self.machine_3D_chk1.configure(text = ' ITER-coils')
		self.machine_3D_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)


	def set_iter_defaults(self, data, flag = None):
		if data is None: data = self.tool_defaults(flag)

		data_dic = {'0-8':data[0:9],'itt':str(int(data[1])),'phistart':repr(data[7]),'MapDirection':int(data[8]),
		'useIcoil':int(data[11]),
		'useM3DC1':int(data[19]),'response':int(data[18]),'selectField':int(data[19]),
		'sigma':int(data[14]),'charge':int(data[15]),'Ekin':repr(data[16]),'Lambda':repr(data[17]),
		'useFilament':str(int(data[12]))}
		
		data_dic['target'] = data[9]
		data_dic['create'] = data[10]
		
		self.useIcoil.set(data_dic['useIcoil'])
		
		self.set_common_defaults(data_dic)
		return data_dic


	def write_iter_coils(self,f):
		f.write('useIcoil(0=no,1=yes)=\t' + str(self.useIcoil.get()) + '\n')
	

	# --- NSTX ---------------------------------------------------------------------------
	def set_nstx_elements(self, row):
		self.machine_3D_Label.configure(text = "3-D Fields: ")
		self.machine_3D_Label.grid(column = 1, row = row, sticky = tk.E, padx=5, pady=5)
		self.machine_3D_chk1.configure(text = ' EC-coils')
		self.machine_3D_chk1.grid(column = 2, row = row, sticky = tk.E + tk.W, padx=5, pady=5)


	def set_nstx_defaults(self, data, flag = None):
		if data is None: data = self.tool_defaults(flag)
		
		data_dic = {'0-8':data[0:9],'itt':str(int(data[1])),'phistart':repr(data[7]),'MapDirection':int(data[8]),
		'useIcoil':int(data[11]),
		'useM3DC1':int(data[19]),'response':int(data[18]),'selectField':int(data[19]),
		'sigma':int(data[14]),'charge':int(data[15]),'Ekin':repr(data[16]),'Lambda':repr(data[17]),
		'useFilament':str(int(data[12]))}
		
		data_dic['target'] = data[9]
		data_dic['create'] = data[10]
		
		self.useIcoil.set(data_dic['useIcoil'])
		
		self.set_common_defaults(data_dic)
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


	def set_mast_defaults(self, data, flag = None):
		if data is None: data = self.tool_defaults(flag)

		data_dic = {'0-8':data[0:9],'itt':str(int(data[1])),'phistart':repr(data[7]),'MapDirection':int(data[8]),
		'useCcoil':int(data[13]),'useIcoil':int(data[14]),
		'useM3DC1':int(data[10]),'response':int(data[9]),'selectField':int(data[10]),
		'sigma':int(data[16]),'charge':int(data[17]),'Ekin':repr(data[18]),'Lambda':repr(data[19]),
		'useFilament':str(int(data[15]))}
		
		data_dic['target'] = data[11]
		data_dic['create'] = data[12]
		
		self.useCcoil.set(data_dic['useCcoil'])
		self.useIcoil.set(data_dic['useIcoil'])
		
		self.set_common_defaults(data_dic)
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


	def set_MinMax_elements(self, row):
		self.xmin = tk.StringVar()
		self.xmax = tk.StringVar()
		self.Nx = tk.StringVar()
		self.ymin = tk.StringVar();
		self.ymax = tk.StringVar();
		self.Ny = tk.StringVar(); 
		
		# Min
		self.xmin_entry = tk.Entry(self.frame, width = 17, textvariable = self.xmin)
		self.xmin_entry.grid(column = 2, row = row, columnspan = 2)
		tk.Label(self.frame, text = "  Min").grid(column = 2, row = row, sticky = tk.W )
		
		# Max
		self.xmax_entry = tk.Entry(self.frame, width = 17, textvariable = self.xmax)
		self.xmax_entry.grid(column = 4, row = row, columnspan = 2)
		tk.Label(self.frame, text = "Max").grid(column = 4, row = row, sticky = tk.W )

		# Nx -> Nth
		row += 1
		self.Nx_entry = tk.Entry(self.frame, width = 17, textvariable = self.Nx)
		self.Nx_entry.grid(column = 2, row = row, columnspan = 2)
		tk.Label(self.frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )
		
		self.pi_text = tk.Text(self.frame, height= 1, width = 30, bd  = 0, takefocus = 0, bg = self.frame.cget('bg'), relief = tk.FLAT)
		self.pi_text.grid(column = 4, row = row, columnspan = 2); self.pi_row = row
		self.pi_text.insert(1.0, 'pi = 3.141593 2pi = 6.283185')
		self.pi_text.configure(state = "disabled")
		
		row += 1
		
		# Min
		self.ymin_entry = tk.Entry(self.frame, width = 17, textvariable = self.ymin)
		self.ymin_entry.grid(column = 2, row = row, columnspan = 2 )
		tk.Label(self.frame, text = "  Min").grid(column = 2, row = row, sticky = tk.W )

		# Max
		self.ymax = tk.StringVar(); 
		self.ymax_entry = tk.Entry(self.frame, width = 17, textvariable = self.ymax)
		self.ymax_entry.grid(column = 4, row = row, columnspan = 2 )
		tk.Label(self.frame, text = "Max").grid(column = 4, row = row, sticky = tk.W )

		# Ny -> Nr
		self.Ny = tk.StringVar(); row += 1
		self.Ny_entry = tk.Entry(self.frame, width = 17, textvariable = self.Ny)
		self.Ny_entry.grid(column = 2, row = row, columnspan = 2 )
		tk.Label(self.frame, text = "     #").grid(column = 2, row = row, sticky = tk.W )


	def set_toroidatTurn_element(self, row):
		# --- toroidal turns ---
		self.itt = tk.StringVar()
		tk.Entry(self.frame, width = 7, textvariable = self.itt).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "tor. Iterations").grid(column = 1, row = row, sticky = tk.E )
		return row + 1
	

	def set_phistart_element(self, row):
		# --- phistart ---
		self.phistart = tk.StringVar()
		tk.Entry(self.frame, width = 7, textvariable = self.phistart).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "tor. Angle [deg]").grid(column = 1, row = row, sticky = tk.E )
		tk.Label(self.frame, text = "For Machine coord. use negative angles").grid(column = 3, row = row, columnspan = 3, sticky = tk.W )
		return row + 1


	def set_MapDirection_element(self, row, state = tk.NORMAL):
		self.MapDirection = tk.IntVar()
		tk.Radiobutton(self.frame, text = '+1', variable = self.MapDirection, value = 1, state = state).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = '-1', variable = self.MapDirection, value = -1, state = state).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'both', variable = self.MapDirection, value = 0, state = state).grid(column = 4, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "Map Direction").grid(column = 1, row = row, sticky = tk.E )
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

		# --- separator ---
		separator2 = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------")
		separator2.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- magnetic field ---
		row += 1
		tk.Radiobutton(self.frame, text = 'g-file Vacuum', variable = self.selectField, value = -1,
			command = self.activate_response).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, width = 8, text = 'M3D-C1', variable = self.selectField, value = 0,
			command = self.activate_response).grid(column = 3, row = row, sticky = tk.W + tk.E )
		tk.Radiobutton(self.frame, text = 'VMEC', variable = self.selectField, value = -3,
			command = self.activate_response).grid(column = 4, row = row, sticky = tk.W + tk.E )
		rb = tk.Radiobutton(self.frame, width = 9, text = 'SIESTA', variable = self.selectField, value = -2,
			command = self.activate_response)
		rb.grid(column = 5, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "Field").grid(column = 1, row = row, sticky = tk.E )
		if not self.SIESTA: rb.configure(state = tk.DISABLED)

		# --- M3DC1 ---
		row += 1; self.row_M3DC1 = row
		self.useM3DC1_R1 = tk.Radiobutton(self.frame, text = 'Equilibrium', variable = self.useM3DC1, value = 0)
		self.useM3DC1_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.useM3DC1_R2 = tk.Radiobutton(self.frame, text = 'Pert.', variable = self.useM3DC1, value = 1)
		self.useM3DC1_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.useM3DC1_R3 = tk.Radiobutton(self.frame, text = 'Eq + Pert', variable = self.useM3DC1, value = 2)
		self.useM3DC1_R3.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.useM3DC1_label = tk.Label(self.frame, text = "use M3DC1")
		self.useM3DC1_label.grid(column = 1, row = row, sticky = tk.E )

		row += 1
		self.response_R1 = tk.Radiobutton(self.frame, text = 'Off', variable = self.resp, value = 0, command = self.make_response)
		self.response_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.response_R2 = tk.Radiobutton(self.frame, text = 'On', variable = self.resp, value = 1, command = self.make_response)
		self.response_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.response_label = tk.Label(self.frame, text = "Plasma response")
		self.response_label.grid(column = 1, row = row, sticky = tk.E )

		self.nresp.set('0');
		self.nresp_entry = tk.Entry(self.frame, width = 4, textvariable = self.nresp, validate = 'all',
        	validatecommand = (okayCommand, '%P'), invalidcommand = notOkayCommand)
		self.nresp_entry.grid(column = 4, row = row, sticky = tk.E)
		self.nresp_label = tk.Label(self.frame, text = "       Time")
		self.nresp_label.grid(column = 4, row = row, sticky = tk.W )

		# --- VMEC ---
		row = self.row_M3DC1
		#self.specialFieldFile1_entry = tk.Entry(self.frame, width = 10, textvariable = self.specialFieldFile1, 
		#	validate = 'focusout', validatecommand = (fileFound, '%d','%P'))
		self.specialFieldFile1_entry = acE.AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
				width = 59, textvariable = self.specialFieldFile1)
		self.specialFieldFile1_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
		self.specialFieldFile1_label = tk.Label(self.frame, text = "VMEC wout: ")
		self.specialFieldFile1_label.grid(column = 1, row = row, sticky = tk.E )
		#self.specialFieldFile1_entry.bind('<Return>', lambda event: self.file_Found(event,self.specialFieldFile1.get()))
		
		row += 1

		#self.specialFieldFile2_entry = tk.Entry(self.frame, width = 10, textvariable = self.specialFieldFile2, 
		#	validate = 'focusout', validatecommand = (fileFound, '%d','%P'))
		self.specialFieldFile2_entry = acE.AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
				width = 59, textvariable = self.specialFieldFile2)
		self.specialFieldFile2_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
		self.specialFieldFile2_label = tk.Label(self.frame, text = "Xpand data: ")
		self.specialFieldFile2_label.grid(column = 1, row = row, sticky = tk.E )
		#self.specialFieldFile2_entry.bind('<Return>', lambda event: self.file_Found(event,self.specialFieldFile2.get()))

		# --- SIESTA ---
		row = self.row_M3DC1
		#self.specialFieldFile3_entry = tk.Entry(self.frame, width = 10, textvariable = self.specialFieldFile3, 
		#	validate = 'focusout', validatecommand = (fileFound, '%d','%P'))
		self.specialFieldFile3_entry = acE.AutocompleteEntry(os.listdir(self.path.get()), self.frame, listboxLength = 6, 
				width = 59, textvariable = self.specialFieldFile3)
		self.specialFieldFile3_entry.grid(column = 2, row = row, columnspan = 4,sticky = tk.E+tk.W)
		self.specialFieldFile3_label = tk.Label(self.frame, text = "SIESTA data: ")
		self.specialFieldFile3_label.grid(column = 1, row = row, sticky = tk.E )
		#self.specialFieldFile3_entry.bind('<Return>', lambda event: self.file_Found(event,self.specialFieldFile3.get()))

		return self.row_M3DC1 + 2


	def set_particle_elements(self, row):
		self.sigma = tk.IntVar()
		self.charge = tk.IntVar()
		self.Ekin = tk.StringVar()
		self.Lambda = tk.StringVar()
		self.useFilament = tk.StringVar()
		self.nproc = tk.StringVar()

		# --- separator ---
		separator3 = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------")
		separator3.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- particles ---
		row += 1
		tk.Radiobutton(self.frame, text = 'Field lines', variable = self.sigma, value = 0, 
			command = self.show_particle_params).grid(column = 2, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(self.frame, text = 'Co pass', variable = self.sigma, value = 1,
			command = self.show_particle_params).grid(column = 3, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Radiobutton(self.frame, text = 'Counter pass', variable = self.sigma, value = -1,
			command = self.show_particle_params).grid(column = 4, row = row, sticky = tk.W + tk.E + tk.N)
		tk.Label(self.frame, text = "Orbits").grid(column = 1, row = row, sticky = tk.E + tk.N)

		row += 1; self.row_particle = row
		self.charge_R1 = tk.Radiobutton(self.frame, text = 'Electrons', variable = self.charge, value = -1)
		self.charge_R1.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.charge_R2 = tk.Radiobutton(self.frame, text = 'Ions', variable = self.charge, value = 1)
		self.charge_R2.grid(column = 3, row = row, sticky = tk.W + tk.E )
		self.charge_label = tk.Label(self.frame, text = "Species")
		self.charge_label.grid(column = 1, row = row, sticky = tk.E )

		row += 1
		self.Ekin_entry = tk.Entry(self.frame, width = 7, textvariable = self.Ekin)
		self.Ekin_entry.grid(column = 2, row = row, sticky = tk.W + tk.E )
		self.Ekin_label = tk.Label(self.frame, text = "kin Energy [keV]")
		self.Ekin_label.grid(column = 1, row = row, sticky = tk.E )

		self.Lambda_entry = tk.Entry(self.frame, width = 7, textvariable = self.Lambda)
		self.Lambda_entry.grid(column = 4, row = row, sticky = tk.W + tk.E )
		self.Lambda_label = tk.Label(self.frame, text = "Energy ratio")	
		self.Lambda_label.grid(column = 3, row = row, sticky = tk.E )

		# --- separator ---
		row += 1
		separator4 = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------")
		separator4.grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W )

		# --- Filament ---
		row += 1
		tk.Entry(self.frame, width = 7, textvariable = self.useFilament).grid(column = 2, row = row, sticky = tk.W + tk.E )
		tk.Label(self.frame, text = "# current Filaments").grid(column = 1, row = row, sticky = tk.E )

		return row + 1
		
		
	def set_run_mpi_elements(self, row):
		# --- number of processes for mpi ---
		self.nproc.set(str(4))
		self.nproc_entry = tk.Entry(self.frame, width = 4, textvariable = self.nproc)
		self.nproc_entry.grid(column = 1, row = row, sticky = tk.E)
		tk.Label(self.frame, text = "       # Procs").grid(column = 1, row = row, sticky = tk.W )


	def set_bottom_elements(self, row, SIESTA = True):
		# --- separator ---
		separator = tk.Label(self.frame, text = "---------------------------------------------------------------------------------------")
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
		self.activate_SIESTA_response()


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


	def activate_SIESTA_response(self):
		if(self.selectField.get() == -2):
			self.specialFieldFile3_label.grid(column = 1, row = self.row_M3DC1, sticky = tk.E, padx=5, pady=5)
			self.specialFieldFile3_entry.grid(column = 2, row = self.row_M3DC1, columnspan = 4, sticky = tk.E+tk.W, padx=5, pady=5)
			self.createFlag.set('psi')
			self.refresh_grid_labels()
			self.create_R1.configure(state=tk.DISABLED)
			self.create_R2.configure(text = 's,u')
			try: self.create_R3.configure(state=tk.DISABLED)
			except: pass
		else:
			self.specialFieldFile3_label.grid_forget()
			self.specialFieldFile3_entry.grid_forget()
			if self.SIESTA:
				self.refresh_grid_labels()
				self.create_R1.configure(state=tk.NORMAL)
				self.create_R2.configure(text = 'psi_n')
				try: self.create_R3.configure(state=tk.NORMAL)
				except: pass


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
			nresp = 1
			while(os.path.isfile(path + 'time_' + format(nresp,'03d') + '.h5')):
				nresp += 1
				
			self.nresp.set(str(nresp - 1))
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
			return c1path
		else: return path
		
		
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
		print 'invalid Entry: requested time_xxx.h5 file not found in', self.path_to_m3dc1()


	# --- Validate the file entries ---
	def file_Found(self, event = None, file = ''):
		if len(file) == 0: return False
		if not '/' in file:
			file = os.path.abspath(self.path.get()) + '/' + file
		if not os.path.isfile(file):
			print 'Requested file not found:', file
			return False	# file not found
		else:
			self.specialFieldFile1_entry.configure(validate = 'focusout')
			self.specialFieldFile2_entry.configure(validate = 'focusout')
			return True		# file found
		
	
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
		self.useFilament.set(data['useFilament'])
		
		self.activate_response()
		self.show_particle_params()
		
		self.specialFieldFile1_entry.update_baseSearchPath(self.path.get())
		self.specialFieldFile2_entry.update_baseSearchPath(self.path.get())
		self.specialFieldFile3_entry.update_baseSearchPath(self.path.get())
		

	# --- Function, executed when Button is pressed ---
	def run_common_funct(self, name, shellCall, data, shellFlags = ''):
		if(self.Shot.get() == '') | (self.Time.get() == ''):
			print 'You must enter a Shot # and a Time'
			return

		# make gPath absolute & check that gPath ends with a /
		self.gPath.set(os.path.abspath(self.gPath.get()))
		if not (self.gPath.get()[-1] == '/'): self.gPath.set(self.gPath.get() + '/')
	
		# convert relative path to absolute
		cwd = os.getcwd()
		path = os.path.abspath(self.path.get())
		if not (path[-1] == '/'): path += '/'
		chk = not (cwd + '/' == path)
		
		# set flags
		# do this before the os.chdir, because the os.path.abspath recognizes the chdir, but self.path stays the same
		shellFlags += self.make_VMEC_SIESTA_shell_flags()
			
		# change to working dir, write contol file(s), launch code, and return to original dir
		if chk: os.chdir(path)
		self.writeControlFile(name)
		if(HOST == 'head.cluster'):		# Drop Cluster
			self.write_qsub_file(data, self.tag.get(), shellFlags)
			#call('qsub run_MAFOTjob', shell = True)
			print 'qsub run_MAFOTjob'
		else:
			#call(shellCall + shellFlags + ' &', shell = True)
			print 'running in dir:', os.getcwd()
			print shellCall + shellFlags + ' &'
		if chk: os.chdir(cwd)


	# --- Write qsub File on Drop Cluster ---
	# here: shellFlags already embedded in shellCall
	def write_common_qsub_file(self, tooltag, nproc, tag, shellCall, mpi = True):
		with open('run_MAFOTjob', 'w') as f:
			f.write('#$ -N ' + tooltag + tag.translate(None, '_+- ') + '\n')
			f.write('#$ -cwd \n')
			f.write('#$ -o ' + HOME + '/work/batch.out \n')
			f.write('#$ -e ' + HOME + '/work/batch.err \n')
			f.write('#$ -S /bin/bash \n')
			f.write('#$ -V \n')
			f.write('#$ -q all.q \n')
			if mpi:
				f.write('#$ -pe mpi ' + str(nproc) + ' \n')
				f.write('source /etc/profile.d/modules.sh\n')
				f.write('module load openmpi-1.6/gcc \n')
			f.write(shellCall + '\n')


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
		f.write('Field(-3=VMEC,-2=SIESTA,-1=gfile,M3DC1:0=Eq,1=I-coil,2=both)=\t' + str(self.selectField.get() + self.useM3DC1.get()) + '\n')


	def write_Ctrl_particles(self, f):
		f.write('ParticleDirection(1=pass,-1=co-pass,0=field-lines)=\t' + str(self.sigma.get()) + '\n')
		f.write('PartileCharge(-1=electrons,>=1=ions)=\t' + str(self.charge.get()) + '\n')
		f.write('Ekin[keV]=\t' + self.Ekin.get() + '\n')
		f.write('lambda=\t' + self.Lambda.get() + '\n')


	def write_Ctrl_center(self, f):
		f.write('phistart(deg)=\t' + self.phistart.get() + '\n')
		f.write('MapDirection=\t' + str(self.MapDirection.get()) + '\n')
		if self.MachFlag.get() in ['dt', 'mast']:
			self.write_Ctrl_resonse(f)
	
	
	def write_Ctrl_bottom(self, f):
			self.write_coils(self.MachFlag.get(), f)
			if self.MachFlag.get() in ['iter', 'nstx', 'mast']:
				f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
			if self.MachFlag.get() in ['iter', 'nstx']:
				f.write('useTe_profile(0=no)=	0\n')
			self.write_Ctrl_particles(f)
			if self.MachFlag.get() in ['iter', 'nstx']:
				self.write_Ctrl_resonse(f)
			if self.MachFlag.get() in ['dt']:
				f.write('useFilament(0=no)=\t' + self.useFilament.get() + '\n')
				f.write('useTe_profile(0=no)=	0\n')
				self.write_d3d_errorFileds(f)
			f.write('pi=\t3.141592653589793\n')
			f.write('2*pi=\t6.283185307179586\n')


	# --- Overloaded Functions -----------------------------------------------------------
	
	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		return

	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, shellFlags):
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
		if self.MachFlag.get() == 'iter':
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 
					0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 1, 0, 1, 0, 0, 0, 1, 
					100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 1, 0, 1, 1, 0, 0, 1, 
					100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		else:
			data = [0, 300, 0.6, 0.95, 0, 0, 40, 0, 1, 0, -1, 1, 3, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 0, 0, 0, 3.141592653589793, 6.283185307179586]
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
	def write_qsub_file(self, nproc, tag, shellFlags):
		tool = self.MachFlag.get() + 'plot_mpi'
		tooltag = 'P'
		shellCall = 'mpirun -n ' + str(nproc) + ' ' + tool + ' _plot.dat ' + tag + shellFlags
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
		tk.Label(frame, text = "R [m]").grid(column = 1, row = row, sticky = tk.E)
		
		self.set_MinMax_elements(row)
		
		# --- y -> r ---
		row += 2
		tk.Label(frame, text = "Z [m]").grid(column = 1, row = row, sticky = tk.E )
		
		row += 2
				
		# --- invisible separator ---
		tk.Label(frame, text = "").grid(column = 1, row = row, columnspan = 5, sticky = tk.E + tk.W, ipady = 1)
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


	# --- turn on/off period entry ---
	def activate_entrys(self):
		data = self.read_par_file(['_fix.dat'])
		if(abs(self.HypRPt.get()) == 1):
			data = self.tool_defaults(None)	# always overwrite whats in any existing file
			self.period_entry.configure(state=tk.DISABLED)
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
		self.Nx.set(str(int(data[6]**0.5)))
		self.xmin.set(repr(data[2]))
		self.xmax.set(repr(data[3]))
		self.Ny.set(str(int(data[6]**0.5)))
		self.activate_MinMax_entrys(abs(self.HypRPt.get()) == 1)
		

	# --- define machine specific default values ---
	def tool_defaults(self, flag):	
		if self.MachFlag.get() == 'iter':
			data = [1e-4, 0, 4, 5.85, -4.4, -3.3, 900, 0, 1, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 
					0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			data = [1e-4, 0, 0.25, 0.55, -1.5, -1, 900, 0, 1, 1, 5, 1, 0, 0, 0, 1, 
					100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			data = [1e-4, 0, 0.25, 0.55, -1.5, -1, 900, 0, 1, 1, 5, 1, 1, 0, 0, 1, 
					100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		else:
			data = [1e-4, 0, 1.1, 1.6, -1.45, -0.8, 900, 0, 1, 0, -1, 1, 5, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 0, 0, 0, 3.141592653589793, 6.283185307179586]
		return data


	# --- set default values ---
	def set_defaults(self):
		data = self.read_par_file(['_fix.dat'])
		data = self.set_machine_defaults(data, self.MachFlag.get())						
		self.shift = data['0-8'][0]
		self.activate_entrys()
		

	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		tool = self.MachFlag.get() + 'fix'
		name = '_fix.dat'
		shellCall = tool + ' ' + name + ' ' + str(int(self.HypPt.get())) + ' ' + self.tag.get()
		self.run_common_funct(name, shellCall, int(self.HypPt.get()))


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, period, tag, shellFlags):
		tool = self.MachFlag.get() + 'fix'
		tooltag = 'fix'
		shellCall = tool + ' _fix.dat ' + str(period) + ' ' + tag + shellFlags
		self.write_common_qsub_file(tooltag, 1, tag, shellCall, mpi = False)


	# --- write Control File to current working dir ---
	def writeControlFile(self, name):
		N = int(self.Nx.get()) * int(self.Ny.get())
		with open(name, 'w') as f:
			self.write_headerLines(self.MachFlag.get(), f)
			f.write('shift=\t' + repr(self.shift) + '\n')
			f.write('itt=\t0\n')			
			f.write('Rmin=\t' + self.xmin.get() + '\n')
			f.write('Rmax=\t' + self.xmax.get() + '\n')
			f.write('Zmin=\t' + self.ymin.get() + '\n')
			f.write('Zmax=\t' + self.ymax.get() + '\n')
			f.write('N=\t' + str(N) + '\n')
			self.write_Ctrl_center(f)
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t0\n')			
			f.write('createPoints(0=setr,3=setpsi,5=setR)=\t5\n')			
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
		if self.MachFlag.get() == 'iter':
			data = [1e-4, 0, 4, 5.85, -4.4, -3.3, 900, 0, 1, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 
					0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			data = [1e-4, 0, 1, 1.3, 4.1, 4.6, 900, 0, 1, 1, 0, 1, 0, 0, 0, 1, 
					100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			data = [1e-4, 0, 1, 1.3, 4.1, 4.6, 900, 0, 1, 1, 0, 1, 0, 0, 0, 1, 
					100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		else:
			data = [1e-4, 0, 1, 1.3, 4.1, 4.6, 900, 0, 1, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
					100, 0.1, 0, 0, 0, 0, 3.141592653589793, 6.283185307179586]
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
	def write_qsub_file(self, fixfile, tag, shellFlags):
		tool = self.MachFlag.get() + 'man'
		tooltag = 'M'
		shellCall = tool + ' _fix.dat ' + fixfile + ' ' + tag + shellFlags
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
		row = self.set_MapDirection_element(row, state = tk.DISABLED)

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
		if self.MachFlag.get() == 'iter':
			self.target_RB1.configure(text = 'Inner')
			self.target_RB2.configure(text = 'Outer')
			self.target_RB1.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid_forget()
			self.target_RB4.grid_forget()
			self.U_chkBtn.grid_forget()
		elif self.MachFlag.get() == 'nstx':
			self.target_RB1.configure(text = 'Inn-up')
			self.target_RB2.configure(text = 'Out-up')
			self.target_RB3.configure(text = 'In-dwn')
			self.target_RB4.configure(text = 'Out-dwn')
			self.target_RB1.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid(column = 4, row = row, sticky = tk.W + tk.E )
			self.target_RB4.grid(column = 5, row = row, sticky = tk.W + tk.E )
			self.U_chkBtn.grid(column = 1, row = row, sticky = tk.E)
		elif self.MachFlag.get() == 'mast':
			self.target_RB1.configure(text = 'Inner')
			self.target_RB2.configure(text = 'Outer')
			self.target_RB1.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid_forget()
			self.target_RB4.grid_forget()
			self.U_chkBtn.grid_forget()
		else:
			self.target_RB1.configure(text = 'Inner')
			self.target_RB2.configure(text = 'Outer')
			self.target_RB3.configure(text = 'Shelf')		
			self.target_RB1.grid(column = 2, row = row, sticky = tk.W + tk.E )
			self.target_RB2.grid(column = 3, row = row, sticky = tk.W + tk.E )
			self.target_RB3.grid(column = 4, row = row, sticky = tk.W + tk.E )
			self.target_RB4.grid_forget()
			self.U_chkBtn.grid_forget()
		

	# --- Change Labels on grid variables, depending on TargetFlag ---
	def refresh_grid_labels(self):
		if self.MachFlag.get() == 'iter':
			if(self.TargetFlag.get() == 1):
				self.y_label.configure(text = "t [cm]")
				self.Info.configure(text = "t in [-166 <--> 31] is length along the wall; t < 0: curve upwards, t > 0: linear outwards")
			elif(self.TargetFlag.get() == 2):
				self.y_label.configure(text = "t [cm]")
				self.Info.configure(text = "t in [-39 <--> 163] is length along the wall; t < 0: linear inwards, t > 0: curve upwards")
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
		elif self.MachFlag.get() == 'mast':
			if(self.TargetFlag.get() == 1):
				self.y_label.configure(text = "Z [m]")
				self.Info.configure(text = "inner target, Z in [-1.6835 <--> -1.22909]")
			elif(self.TargetFlag.get() == 2):
				self.y_label.configure(text = "R [m]")
				self.Info.configure(text = "outer target, R in [0.7835 <--> 1.9]")
		else:	# d3d
			if(self.TargetFlag.get() == 3):
				self.y_label.configure(text = "t [0 <--> 1]")
				self.Info.configure(text = "t > 0: Shelf outwards, t = 0: Nose edge")
			elif(self.TargetFlag.get() == 2):
				self.y_label.configure(text = "t [0 <--> 1]")
				self.Info.configure(text = "t > 0: Divertor Floor outwards, t = 0: Connection 45deg Tile, t = 1: Pump Entry")
			else:
				self.y_label.configure(text = "t [-1 <--> 1]")
				self.Info.configure(text = "t < 0: Centerpost upwards, t > 0: 45deg Tile downwards, t = 0: Connection point")

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
		if self.MachFlag.get() == 'iter':
			if(flag == 1):
				data = [500, 300, -60, 10, 0, 6.283185307179586, 700, 0, 1, 1, 2, 1, 0, 0, 0, 1, 100, 0.1, 
						0, -1, 3.141592653589793, 6.283185307179586]
			elif(flag == 2):
				data = [500, 300, -10, 80, 0, 6.283185307179586, 900, 0, -1, 2, 2, 1, 0, 0, 0, 1, 100, 0.1, 
						0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			if(flag == 1):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, 1.05, 1.578, 0, 6.283185307179586, 400, 0, 1, 1, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, 1.1714, 1.578, 0, 6.283185307179586, 400, 0, 1, 1, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
			elif(flag == 2):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, 0.435, 1.0433, 0, 6.283185307179586, 400, 0, -1, 2, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, 0.2979, 1.0433, 0, 6.283185307179586, 400, 0, -1, 2, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
			elif(flag == 3):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, -1.578, -1.05, 0, 6.283185307179586, 400, 0, -1, 3, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, -1.578, -1.1714, 0, 6.283185307179586, 400, 0, -1, 3, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
			elif(flag == 4):
				if self.UpgradeFlag.get() == 1:
					data = [500, 500, 0.435, 1.0433, 0, 6.283185307179586, 400, 0, 1, 4, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
				else:
					data = [500, 500, 0.2979, 1.0433, 0, 6.283185307179586, 400, 0, 1, 4, 2, 1, 0, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			if(flag == 1):
				data = [500, 500, -1.6835, -1.229, 0, 6.283185307179586, 400, 0, -1, 3, 2, 1, 1, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
			elif(flag == 2):
				data = [500, 500, 0.7835, 1.9, 0, 6.283185307179586, 400, 0, 1, 4, 2, 1, 1, 0, 0, 1, 100, 0.1, 
							0, -1, 3.141592653589793, 6.283185307179586]
		else:	# d3d
			if(flag == 3):
				data = [500, 300, 0.0, 0.1, 0, 6.283185307179586, 100, 0, 1, 0, -1, 3, 2, 1, 1, 1, 0, 1, 
							100, 0.1, 0, 0, 0, 0, 3.141592653589793, 6.283185307179586]		
			elif(flag == 2):
				data = [500, 300, 0.9, 1.0, 0, 6.283185307179586, 100, 0, 1, 0, -1, 2, 2, 1, 1, 1, 0, 1, 
							100, 0.1, 0, 0, 0, 0, 3.141592653589793, 6.283185307179586]
			else:
				data = [500, 300, -0.1, 0.3, 0, 6.283185307179586, 400, 0, -1, 0, -1, 1, 2, 1, 1, 1, 0, 1, 
							100, 0.1, 0, 0 ,0 ,0, 3.141592653589793, 6.283185307179586]
		return data


	# --- set default values ---
	def set_defaults(self):
		if self.MachFlag.get() == 'iter':
			data = self.read_par_file(['_inner.dat','_outer.dat'])
		elif self.MachFlag.get() == 'nstx':
			data = self.read_par_file(['_innerup.dat','_outerup.dat','_innerdwn.dat','_outerdwn.dat'])
		elif self.MachFlag.get() == 'mast':
			data = self.read_par_file(['_inner.dat','_outer.dat'])
		else:
			data = self.read_par_file(['_inner.dat','_outer.dat', '_shelf.dat'])
		data = self.set_machine_defaults(data, self.MachFlag.get(), flag = 1)
		if(data['target'] == 0): self.TargetFlag.set(1)
		else: self.TargetFlag.set(int(data['target']))
		self.UpgradeFlag.set(0)
		self.refresh_grid_labels()


	# --- set control file name ---
	def set_name(self):
		if self.MachFlag.get() == 'iter':
			if(self.TargetFlag.get() == 2): name = '_outer.dat'
			else: name = '_inner.dat'
		elif self.MachFlag.get() == 'nstx':
			if(self.TargetFlag.get() == 4): name = '_outerdwn.dat'
			elif(self.TargetFlag.get() == 3): name = '_innerdwn.dat'
			elif(self.TargetFlag.get() == 2): name = '_outerup.dat'
			else: name = '_innerup.dat'
		elif self.MachFlag.get() == 'mast':
			if(self.TargetFlag.get() == 2): name = '_outer.dat'
			else: name = '_inner.dat'
		else:
			if(self.TargetFlag.get() == 3): name = '_shelf.dat'
			elif(self.TargetFlag.get() == 2): name = '_outer.dat'
			else: name = '_inner.dat'
		return name
		

	# --- Function, executed when Button is pressed ---
	def run_funct(self):
		name = self.set_name()
		tool = self.MachFlag.get() + 'foot_mpi'
		shellCall = MPIRUN + ' -n ' + str(int(self.nproc.get())) + ' ' + tool + ' ' + name + ' ' + self.tag.get()
		self.run_common_funct(name, shellCall, int(self.nproc.get()))


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, shellFlags):
		name = self.set_name()
		tool = self.MachFlag.get() + 'foot_mpi'
		tooltag = 'F'
		shellCall = 'mpirun -n ' + str(nproc) + ' ' + tool + ' ' + name + ' ' + tag + shellFlags
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
			f.write('target(0=cp,1=inner,2=outer,3=shelf)=\t' + str(self.TargetFlag.get()) + '\n')			
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
		self.pointsFile_entry = acE.AutocompleteEntry(os.listdir(self.path.get()), frame, listboxLength = 6, 
				width = 59, textvariable = self.pointsFile)
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
		if self.MachFlag.get() == 'iter':
			if(flag == 'psi'):
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 1, 3, 1, 0, 0, 0, 1, 100, 0.1, 
						0, -1, 3.141592653589793, 6.283185307179586]
			else:
				data = [520, 200, 4.0, 6.5, -4.6, -2.0, 500, 0, 0, 1, 0, 1, 0, 0, 0, 1, 100, 0.1, 
						0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'nstx':
			if(flag == 'psi'):
				data = [300, 500, 0.88, 1.02, 0, 6.283185307179586, 150, 0, 0, 1, 3, 1, 0, 0, 0, 1, 
						100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
			else:
				data = [100, 500, 0.17, 0.9, -1.65, -1.0, 100, 0, 0, 1, 0, 1, 0, 0, 0, 1, 
						100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		elif self.MachFlag.get() == 'mast':
			if(flag == 'psi'):
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 1, 3, 1, 0, 0, 0, 1, 100, 0.1, 
						0, -1, 3.141592653589793, 6.283185307179586]
			else:
				data = [100, 200, 0.7, 1.9, -1.65, -1.0, 100, 0, 0, 1, 0, 1, 0, 0, 0, 1, 
						100, 0.1, 0, -1, 3.141592653589793, 6.283185307179586]
		else:
			if(flag == 'psi'):
				data = [1200, 200, 0.88, 1.02, 0, 6.283185307179586, 700, 0, 0, 0, -1, 1, 3, 1, 1, 1, 0, 1, 
						100, 0.1, 0, 3.141592653589793, 6.283185307179586]
			else:
				data = [930, 200, 1.0, 1.45, -1.367, -0.902, 900, 0, 0, 0, -1, 1, 0, 1, 1, 1, 0, 1, 
						100, 0.1, 0, 0, 0, 0, 3.141592653589793, 6.283185307179586]
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
			else: print 'invalid entry for psi limit'
		self.run_common_funct(name, shellCall, int(self.nproc.get()), shellFlags)


	# --- Write qsub File on Drop Cluster ---
	def write_qsub_file(self, nproc, tag, shellFlags):
		name = self.set_name()
		tool = self.MachFlag.get() + 'laminar_mpi'
		tooltag = 'L'
		shellCall = 'mpirun -n ' + str(nproc) + ' ' + tool + ' ' + name + ' ' + tag + shellFlags
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
# --- info tab ------------------------------------------------------------------------------------------------
class info_gui:
	def __init__(self, frame):
	
		row = 0
		self.info_text = tk.Text(frame, height= 20, width = 80, bd  = 0, takefocus = 0, 
								 bg = frame.cget('bg'), relief = tk.FLAT, wrap=tk.WORD)
		self.info_text.grid(column = 1, row = row, columnspan = 5, padx=10, pady=10); 
		self.info_text.insert(1.0, 
		'MAFOT Control GUI for DIII-D, ITER, NSTX & MAST \n\n'
		'MAFOT Version 4.02 \n'
		'GUI Version 2.0 \n'
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


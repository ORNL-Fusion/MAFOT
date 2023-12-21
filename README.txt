System requirements:
1) gnu c++ compiler 4.1.2 or newer
2) intel or gnu fortran compiler
3) openmp (comes with gnu c++ compiler 4.2 or newer)

Dependencies:
1) blitz++
2) openmpi
3) HDF5 and netcdf
3) (optional) m3dc1

Dependencies GUI (optional):
1) python 2.7 or 3
2) Tk/Tcl 8.5 (Tkinter with ttk)

Build:
1) Edit make.inc and adjust to local platform. Examples are given in the install/ dir
2) make		 options are: all, d3d, iter, nstx, mast, clean or any of the individual tools

Build GUI (optional):
1) make gui

Set environment
1) check/add path to mpirun to PATH variable
2) add path to all required libraries to LD_LIBRARY_PATH
3) check/add path to created binary to PATH variable

How to run a case:
1) check for g-file
2) check for diiidsup.in file
3) launch mafot_gui.py

or

3) check control file (filename starts with _)
	set shot and time right
	set path to g-file right
4) run code
	$/> mpirun -np <N> -host <hostname> dtlaminar_mpi _lam.dat <arbitrary_single_string> &
	( N: 			   Number of processes 
	  hostname: 		   Names of Hosts, separated by comma
	  arbitrary_single_string: will be added to the output filename, no blanks!
	)

------------------------------------------------------------------------------------
Current Version:
----------------
MAFOT 5.51
GUI 3.1

Version Notes: 
----------------
MAFOT 5.51 -- Dec 2023
- laminar now returns B-field components instead of previously obsolete variables (pitch, yaw amd psimax)

MAFOT 5.5 -- May 2023
- Code can now read arbitrarily named gFiles
- bug fix in structure, where toroidal angle was wrong for dpinit != 1

GUI 3.1 -- Jan 2023
- gFile is now a full pathname, selectable from drop down menu, and sets shot & time, if applicable

MAFOT 5.4.1 -- Feb 2023
- Enable any machine without specific 3D coils

GUI 3.0 -- Jan 2023
- converted to Python3
- PPPL cluster support added

MAFOT 5.4 -- Oct 2021
- Collision class (SULI project) merged into main

MAFOT 5.3 -- Oct 2021
- HEAT inluded in main file structure

MAFOT 5.2 -- Aug 2021
- TCABR tokamak support added, also to d3dplot.py

MAFOT 5.11 -- Nov 2019
- sheath model added

MAFOT 5.1 -- July 2019
- new tools trace and lcfs added
- bug fixes

MAFOT 5.0 -- Jan 2019
- Zeff added to all input files. This required some minor restructuring, which makes old ITER and NSTX input files incompatible.
  For D3D, only the position of useFilaments shifted by one. MAST and C-Mod are unchanged
- electric field for ExB particle drifts added

GUI 2.1 -- Jan 2019
- Zeff added
- electric field input added

MAFOT 4.2 -- Sep 2018
- restart option added to xpand, foot and laminar
- VMEC field line tracing now posible inside s = 1 without xpand
- d3dplot uses now pcolormesh instead of imshow

GUI 2.1 -- Jul 2018
- support for full wall outline in foot, using Swall
- Mapdirection now selectable in foot, because of reversed Bt cases

MAFOT 4.11 -- Jul 2018
- special 2D bisection method added to support s,u calculation from R,Z in xpand

MAFOT 4.1 -- May 2018
- laminar now set up to work with Sheft toolkit on NSTX. Several fixes for NSTX
- 3D wall with any number of toroidal slices
- check for wall crossing within first integration step

MAFOT 4.03 -- Mar 2018
- full wall now accessible for foot_mpi using the swall coordinate and target = 0

GUI 2.0 -- Feb 2018
- complete restructuring of ther code to unify common code segments
- several new features, like use of shell Flags for VMEC/SIESTA, or dtlaminar
- new autocomplete Entry widget for file names, that searches through the directory structure

MAFOT 4.02 -- Feb 2018
- B-coil shift&tilt error added for DIII-D

MAFOT 4.01 -- Jan 2018
- Buswork error field added for DIII-D
- For backward compatibility: c++11 standard no longer required, but recommended
- structure can now use the 3D wall
- compile targets for "each tool in all machines" added to makefile

MAFOT 4.0 -- Dec 2017
- unification of M3DC1 and regular machine file -> *_m3dc1.hxx machine files no longer supported
- simplification of machine specific files by moving larg pieces of common code to mafot.hxx
- bug fix for xpand
- CMOD now fully available (except target plate specification for footprints)
- improved command line options for all tools, including file name support for SIESTA, VMEC and XPAND
- new read-in routine for parameterfile
- 3D field line tracing with structure tool now fully supported
- c++ 2011 standard now required -> compiler flag -std=c++1 added in general

MAFOT 3.9 -- July 2017
- laminar can now compute Lc, psimin, etc. along an EFIT flux surface of constant psi
- d3dplot.py can plot the new laminar plots

MAFOT 3.8 -- May 2017
- C-Mod support added (only plot_mpi for now)
- 3D wall added
- Python "== None" warning fixed

MAFOT 3.75 -- Jul 2016
- divB correction in xpand does not work properly - DO NOT USE
- use of multiple M3DC1 files now allows arbitrary relative phase
- fake island perturbation added (D3D only) 

MAFOT 3.74 -- Apr 2016
- splines separated from EFIT_Class into own header
- function added to read & store header lines from files
- bug fix in VMEC class
- structure tool updated
- divB correction added to xpand
- minor fixes in python tools

GUI 1.43 -- Feb 2016
- bug fix in all plot calls for Drop cluster only where it calls the wrong job file
- replaced all ${NSLOTS} with str(nproc) in Drop cluster job files
- add „source /etc/profile.d/modules.sh“ to the mpirun job batch files for Drop as a workaround for the module load bug in bash

MAFOT 3.73 -- Oct 2015
- shell command line flags enabled in all tools, use -h for details

GUI 1.42 -- Oct 2015
- shell command flags can be entered through the "File Tag" input now

MAFOT 3.72 -- Aug 2015
- Bug fix in EFIT_class
- fix now searches on RZ grid as well; for period = 1: only first result within grid is returned

GUI 1.41 -- Aug 2015
- changed fix to use RZ grid in D3D and NSTX

MAFOT 3.71 -- Aug 2015
- in xpand: enable force axisymmetriy in VMECthrough n0only
- in xpand: use full VMEC B-field for virtual casing; interpolation error in vacuum field much larger 
            than any noise from including the vacuum field

MAFOT 3.7 -- Aug 2015
- NSTX-U divertor targets and NSTX outer inclined wall added

GUI 1.4 -- Aug 2015
- support for new NSTX & NSTX-U walls added

MAFOT 3.61 -- July 2015
- readfile now counts rows and allocates total array at once
	-> XFIELD read in works now with very large files
- M3DC1 response > 1 now possible -> reads proper time_xxx.h5

GUI 1.3 -- July 2015
- new entry field added for M3DC1 response_time > 1

MAFOT 3.6 -- July 2015
- M3DC1 interface is now a separate class

MAFOT 3.52 -- July 2015
- Bug fix for M3DC1: psi_n is now calculated correctly from M3DC1 equilibrium

GUI 1.2 -- May 2015
- All NSTX tools added

MAFOT 3.51 -- May 2015
- Bug fix for NSTX: M3DC1 settings now enabled

MAFOT 3.5 -- Feb 2015
- Include the support for external fields from xpand_mpi. This replaces Diagno.
  code looks for file "xpand.dat" ("diagno.dat" no longer recognized, 
  "xpand.dat" can come from the DIAGNO code)

MAFOT 3.41 -- Feb 2015
- BUG fixes in M3DC1 multimode capability for all machines:
    Equilibrium only used from first modefile, not from all
    <machine>sup.in file no longer required, if m3dc1sup.in is used
    field eval fails are now properly handled

GUI 1.15 -- Oct 2014
- SIESTA support added; plot und laminar only

MAFOT 3.4 -- Oct 2014
- SIESTA support added; code looks for file "siesta.dat"

MAFOT 3.3 -- Oct 2014
- VMEC/DIAGNO support added; code looks for files "diagno.dat" and "wout.nc"
- pragma OMP barrier removed to comply with latest compiler rules

GUI 1.14 -- Oct 2014
- VMEC/DIAGNO support added

GUI 1.131 -- Sep 2014
- g-file path can now be a relative path too

MAFOT 3.2 -- Jul 2014
- M3D-C1 interface updated to use Nate's libfusionio -> resistive wall and multiple M3D-C1 files now enabled
- Tools now unified across machines
- plot is now parallel -> plot_mpi (old, serial plot version still exists)

GUI 1.13 -- Jul 2014
- GUI now uses plot_mpi instead of plot

GUI 1.12 -- Feb 2014
- NSTX foot_mpi added

GUI 1.11 -- Nov 2013
- Bug fixes for ITER

GUI 1.1 -- Nov 2013
- ITER added

GUI 1.0 -- Oct 2013
- first release for D3D MAFOT control

MAFOT 3.1 -- Apr 2012
- M3D-C1 interface added; code looks for file "C1.h5"

MAFOT 3.0 -- Jan 2011
- EFIT, IO and Particle Classes introduced -> replace old includes
- includes unified across machines; the parts that are unique to each machine are in a machine include file

MAFOT 2.0 -- March 2009
- current filaments enabled
- support for multiple machines enabled
- ITER, MAST & NSTX added

MAFOT 1.0 -- May 2008
- first release of D3D MAFOT code

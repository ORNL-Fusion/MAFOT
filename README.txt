System requirements:
1) gnu c++ compiler 4.1.2 or newer
2) intel or gnu fortran compiler
3) openmp (comes with gnu c++ compiler 4.2 or newer)

Dependencies:
1) blitz++
2) openmpi 1.6.5
3) (optional) m3dc1

Dependecies GUI (optional):
1) python 2.7
2) Tk/Tcl 8.5 (Tkinter with ttk)
3) pyinstaller 2.1

Build:
1) Edit make.inc and adjust to local platform. Examples are given in the install/ dir
2) make		 options are: all, d3d, iter, nstx, mast, clean or any of the individual tools

Build GUI (optional):
1) Edit make.inc and set PYINSTALLER to the full pathname of the pyinstaller binary
2) make gui

Set environment
1) check/add path to mpirun to PATH variable
2) add path to all required libraries to LD_LIBRARY_PATH
3) check/add path to created binary to PATH variable

How to run a case:
1) check for g-file
2) check for diiidsup.in file
3) launch mafot_gui

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
MAFOT 3.75
GUI 1.43

Version Notes: 
----------------
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

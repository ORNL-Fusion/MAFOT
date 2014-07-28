System requirements:
1) gnu c++ compiler 4.1.2 or larger
2) intel or gnu fortran compiler
3) openmp (comes with gnu c++ compiler 4.2 or larger)

Dependencies:
1) blitz++
2) openmpi
3) (optional) m3dc1

Build:
1) Edit make.inc and adjust to local platform. Examples are given in the install/ dir
2) make		 options are: all, d3d, iter, nstx, mast, clean or any of the individual tools

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

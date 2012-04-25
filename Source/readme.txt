System requirements:

1) gnu c++ compiler 4.1.2 or larger
2) intel fortran compiler
3) openmp (comes with gnu c++ compiler 4.2 or larger)


How to compile:

1) install blitz
2) install openmpi, add path to mpirun to PATH variable
3) set LD_LIBRARY_PATH correct
4) check paths in compile script
5) run script:
	$/> cxx_ifort_mpi dtlaminar
6) add path to created binary to PATH variable 


How to run a case:

1) check for g-file
2) check for diiidsup.in file
3) check control file (filename starts with _)
	set shot and time right
	set path to g-file right
4) run code
	$/> mpirun -np <N> -host <hostname> dtlaminar_mpi _lam.dat <arbitrary_single_string> &
	( N: 			   Number of processes 
	  hostname: 		   Names of Hosts, separated by comma
	  arbitrary_single_string: will be added to the output filename, no blanks!
	)
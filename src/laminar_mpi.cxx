// Program calculates connection length and penetration depth inside the plasma volume
// with particle-drift and time dependent perturbations
// Fortran subroutines are used for perturbations
// A.Wingen						10.6.11

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	2d connection length data for colored contour plot
//			log-file


// uses Parallel computation by open-mpi and openmp
// IMPORTANT: Node 0 starts two threads using openmp. The Master-thread contols communications with MPI-Slaves. The Slave-thread on Node 0 does the same calculations as the MPI-Slaves
// Launch Syntax:				mpirun -np [NoOfNodes] -hostfile [FileName] dtlaminar_mpi [Parameterfile] [Praefix(optional)] &
// alternative Launch Syntax:	mpirun -np [NoOfNodes] -host [HostNames separated by Commata] dtlaminar_mpi [Parameterfile] [Praefix(optional)] &

// Define
//--------
#define USE_MPI
#if defined(ITER)
	#define program_name "iterlaminar_mpi"
#elif defined(NSTX)
	#define program_name "nstxlaminar_mpi"
#elif defined(MAST)
	#define program_name "mastlaminar_mpi"
#else
	#define program_name "dtlaminar_mpi"
#endif

// Include
//--------
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <mafot.hxx>
#include <omp.h>
#include <unistd.h>

// Prototypes  
//-----------
int pitch_angles(double R, double Z, double phi, double& pitch, double& yaw, EFIT& EQD, IO& PAR);

// Switches
//----------

// Golbal Parameters 
//------------------

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// MPI initialize
MPI::Init(argc, argv);
int mpi_rank = MPI::COMM_WORLD.Get_rank();
int mpi_size = MPI::COMM_WORLD.Get_size();

// Variables
EFIT EQD;
int i,j;
int chk,skip_connect;
double ntor,length,psimin,psimax,psiav;
double pitch,yaw;
double xtmp,ytmp;
Range all = Range::all();

int tag,sender;
double Zmin_slave,Zmax_slave,dz;
Array<double,1> send_Z_limits(Range(1,2));
MPI::Status status;

// Use system time as seed(=idum) for random numbers
double now = zeit();
//long idum = long(now);

// defaults
int spare_interior = 0;					// 0: all points are calculated		1: inside psi = psi_interior_limit results are set to fixed values (code runs faster)
double psi_interior_limit = 0.85;		// psi limit for spare_interior == 1
bool use_3Dwall = false;
LA_STRING wall_file = "none";
double Raxis = 0, Zaxis = 0;
LA_STRING new_axis_loc;
int new_axis_loc_idx;
bool use_inputPointsFile = false;
LA_STRING inputPoints_file = "none";
LA_STRING woutfile = "wout.nc";
LA_STRING xpandfile = "xpand.dat";
LA_STRING siestafile = "siesta.dat";
LA_STRING islandfile = "fakeIslands.in";


// Command line input parsing
int c;
opterr = 0;
while ((c = getopt(argc, argv, "hsl:W:A:P:X:V:S:I:")) != -1)
switch (c)
{
case 'h':
	if(mpi_rank < 1)
	{
		cout << "usage: mpirun -n <cores> dtlaminar_mpi [-h] [-l limit] [-s] [-A Raxis,Zaxis] [-I island] " << endl;
		cout << "                                       [-P points] [-S siesta] [-V wout] [-W wall] [-X xpand] file [tag]" << endl << endl;
		cout << "Calculate field line connection length and penetration depth in a poloidal cross-section." << endl << endl;
		cout << "positional arguments:" << endl;
		cout << "  file          Contol file (starts with '_')" << endl;
		cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
		cout << endl << "optional arguments:" << endl;
		cout << "  -h            show this help message and exit" << endl;
		cout << "  -l            flux limit for spare interior, default = 0.85" << endl;
		cout << "  -s            spare calculation of interior, default = No" << endl;
		cout << "  -A            force magnetic axis location, use: -A Raxis,Zaxis , default: use g-file" << endl;
		cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
		cout << "  -P            use separate input file for initial conditions; argument is the file name; default is None" << endl;
		cout << "                File format of columns: R [m], phi [deg, right-handed coord.], Z [m]" << endl;
		cout << "                Header lines start with '#'; no comment lines between/after data possible" << endl;
		cout << "  -S            filename for SIESTA; default, see below" << endl;
		cout << "  -V            filename for VMEC; default, see below" << endl;
		cout << "  -W            use separate 3D Wall-File; default is 2D wall from EFIT file" << endl;
		cout << "  -X            filename for XPAND; default, see below" << endl;
		cout << endl << "Examples:" << endl;
		cout << "  mpirun -n 4 dtlaminar_mpi _lam.dat blabla" << endl;
		cout << "  mpirun -n 12 dtlaminar_mpi -s -l 0.7 _lam.dat skip_inside0.7" << endl;
		cout << endl << "Infos:" << endl;
		cout << "  To use B-field from M3DC1, set response_field >= 0, and provide file in cwd:" << endl;
		cout << "    m3dc1sup.in    ->  location and scale factor for M3DC1 output C1.h5" << endl;
		cout << "  To use B-field from XPAND, set response_field = -3, and provide files in cwd:" << endl;
		cout << "    xpand.dat      ->  B-field on 3D grid from XPAND; use option -X to specify other filename" << endl;
		cout << "    wout.nc        ->  VMEC output; use option -V to specify other filename" << endl;
		cout << "  To use B-field from SIESTA, set response_field = -2, and provide file in cwd:" << endl;
		cout << "    siesta.dat     ->  B-field on 3D grid; use option -S to specify other filename" << endl;
		cout << "  To use B-field for mock-up islands, set response_field = -10, and provide file in cwd:" << endl;
		cout << "    fakeIslands.in ->  each line gives: Amplitude, pol. mode m, tor. mode n, phase [rad]" << endl;
		cout << "                       use option -I to specify other filename" << endl;
	}
	MPI::Finalize();
	return 0;
case 's':
	spare_interior = 1;
	break;
case 'l':
	psi_interior_limit = atof(optarg);
	break;
case 'W':
	use_3Dwall = true;
	wall_file = optarg;
	break;
case 'A':
	new_axis_loc = LA_STRING(optarg);
	new_axis_loc_idx = new_axis_loc.indexOf(",");
	Raxis = atof(new_axis_loc.left(new_axis_loc_idx-1));
	Zaxis = atof(new_axis_loc.mid(new_axis_loc_idx+1));
	break;
case 'P':
	use_inputPointsFile = true;
	inputPoints_file = optarg;
	break;
case 'I':
	islandfile = optarg;
	break;
case 'S':
	siestafile = optarg;
	break;
case 'V':
	woutfile = optarg;
	break;
case 'X':
	xpandfile = optarg;
	break;
case '?':
	if(mpi_rank < 1)
	{
		if (optopt == 'c')
			fprintf (stderr, "Option -%c requires an argument.\n", optopt);
		else if (isprint (optopt))
			fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		else
			fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
	}
	EXIT;
default:
	EXIT;
}
// Input file names
LA_STRING basename;
LA_STRING praefix = "";
if(argc==optind+2) praefix = "_" + LA_STRING(argv[optind+1]);
if(argc>=optind+1) basename = LA_STRING(argv[optind]);
else {if(mpi_rank < 1) cout << "No Input files -> Abort!" << endl; EXIT;}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// check if enough Nodes have been selected
if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes and restart." << endl; EXIT;}

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + "_Node" + LA_STRING(mpi_rank) + ".dat");
ofs2.precision(16);
ofstream ofs3;

// Output
LA_STRING filenameout = "lam" + praefix + ".dat";
if(mpi_rank < 1) outputtest(filenameout);
MPI::COMM_WORLD.Barrier();	// All Nodes wait for Master

// Read parameter file
if(mpi_rank < 1) cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,11,mpi_rank);

// Read EFIT-data
#ifdef USE_XFIELD
if(PAR.response_field == -3)
{
	VMEC vmec;
	if(mpi_rank < 1) cout << "Read VMEC file" << endl;
	ofs2 << "Read VMEC file" << endl;
	vmec.read(woutfile);
	vmec.get_axis(PAR.phistart/rTOd, Raxis, Zaxis);
}
#endif

#ifdef m3dc1
if(PAR.response_field == 0 || PAR.response_field == 2)
{
	M3D.read_m3dc1sup();
	M3D.open_source(PAR.response, PAR.response_field, -1);
	Raxis = M3D.RmAxis;
	Zaxis = M3D.ZmAxis;
	//if(mpi_rank < 1) cout << Raxis << "\t" << Zaxis << endl;
	M3D.unload();
}
#endif

if(Raxis > 0)
{
	if(mpi_rank < 1) cout << "Shift Equilibrium to new axis: Raxis = " << Raxis << "     Zaxis = " << Zaxis << endl;
	ofs2 << "Shift Equilibrium to new axis: Raxis = " << Raxis << "     Zaxis = " << Zaxis << endl;
}
EQD.ReadData(EQD.Shot,EQD.Time,Raxis,Zaxis);
if(mpi_rank < 1) cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Read input points file. Use format for columns: R [m], phi [deg, right-handed coord.], Z [m]
Array<double,2> inputPoints;
int inputPoints_columns;
if(use_inputPointsFile)
{
	if(mpi_rank < 1) cout << "Using Points from file: " << inputPoints_file << endl;
	ofs2 << "Using Points from file: " << inputPoints_file << endl;
	inputPoints_columns = count_column(inputPoints_file);
	readfile(inputPoints_file, inputPoints_columns, inputPoints);
	PAR.create_flag = 0;	// set to "setRZ" mode by default
}

// Read 3D wall file and add to EQD
if(use_3Dwall)
{
	if(mpi_rank < 1) cout << "Using 3D wall from file: " << wall_file << endl;
	ofs2 << "Using 3D wall from file: " << wall_file << endl;
	EQD.set3Dwall(wall_file);
}

// Set starting parameters
int N_variables = 11;
int NZ_slave = 1;
int N = PAR.NR*PAR.NZ;
int N_slave = PAR.NR*NZ_slave;
int NoOfPackages = int(PAR.NZ/NZ_slave);
int N_slave_last = N_slave;

if(use_inputPointsFile)
{
	N = inputPoints.rows();
	if(N < 10*mpi_size) N_slave = 1;
	else if(N < 100*mpi_size) N_slave = 10;
	else N_slave = 100;
	NoOfPackages = int(N/N_slave);
	if(N_slave*NoOfPackages < N) // if not exact, work on one more package
	{
		N_slave_last = N - N_slave*NoOfPackages;
		NoOfPackages += 1;
	}
	else N_slave_last = N_slave;
	PAR.NR = int(N_slave/NZ_slave);
	PAR.NZ = NZ_slave*NoOfPackages;
}

int N_slave_use = N_slave;	// set to private in omp section, so needs to be reinitialized in there.

// Needed for storage of data while calculating
Array<int,1> write_memory(Range(1,NoOfPackages)); //store which data is written to file (0: not calculated yet, 1: calculated, 2: written to file)
int write_max = 0;
int write_last = 0; //tag of last written data
write_memory = 0;

if(PAR.NZ<=1) {dz=0;}
else dz=(PAR.Zmax-PAR.Zmin)/(PAR.NZ-1);
if(dz == 0) PAR.Zmax = PAR.Zmin;

// Prepare particles
PARTICLE FLT(EQD,PAR,mpi_rank);

MPI::COMM_WORLD.Barrier();	// Syncronize all Nodes

// Master only (Node 0)
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank < 1)
{
	// log file
	ofs3.open("log_" + LA_STRING(program_name) + praefix + "_Master" + ".dat");
	ofs3.precision(16);

	ofs3 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
	ofs3 << "No. of Packages = " << NoOfPackages << " Points per Package = " << N_slave << endl << endl;

	// additional parameters for IO
	PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
	PAR.pv[1].name = "R-grid";			PAR.pv[1].wert = PAR.NR;
	PAR.pv[2].name = "Z-grid";			PAR.pv[2].wert = PAR.NZ;
	PAR.pv[3].name = "Rmin";			PAR.pv[3].wert = PAR.Rmin;
	PAR.pv[4].name = "Rmax";			PAR.pv[4].wert = PAR.Rmax;
	PAR.pv[5].name = "Zmin";			PAR.pv[5].wert = PAR.Zmin;
	PAR.pv[6].name = "Zmax";			PAR.pv[6].wert = PAR.Zmax;
	PAR.pv[7].name = "phistart";		PAR.pv[7].wert = PAR.phistart;
	PAR.pv[8].name = "MapDirection";	PAR.pv[8].wert = PAR.MapDirection;
	PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;
	PAR.pv[10].name = "Ekin";			PAR.pv[10].wert = PAR.Ekin;

	if(PAR.create_flag == 3) {PAR.pv[1].name = "psi-grid"; PAR.pv[2].name = "theta-grid"; PAR.pv[3].name = "psimin"; PAR.pv[4].name = "psimax"; PAR.pv[5].name = "thetamin"; PAR.pv[6].name = "thetamax";}
	if(PAR.create_flag == 6) {PAR.pv[1].name = "theta-grid"; PAR.pv[2].name = "phi-grid"; PAR.pv[3].name = "thetamin"; PAR.pv[4].name = "thetamax"; PAR.pv[5].name = "phimin"; PAR.pv[6].name = "phimax"; PAR.pv[7].name = "psi surface";}

	// Output
	ofstream out(filenameout);
	out.precision(16);
	vector<LA_STRING> var(N_variables);
	var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";  var[5] = "psimax";  var[6] = "psiav";  var[7] = "FL pitch angle";  var[8] = "FL yaw angle";
	var[9] = "theta";  var[10] = "psi";
	if(PAR.create_flag == 3) {var[0] = "theta";  var[1] = "psi";  var[9] = "R[m]";  var[10] = "Z[m]";}
	if(PAR.create_flag == 6) {var[0] = "phi";  var[1] = "theta";  var[9] = "R[m]";  var[10] = "Z[m]";}
	if(PAR.response_field == -2) {var[0] = "u"; var[3] = "s";}
	if(use_inputPointsFile) var[9] = "phi";
	PAR.writeiodata(out,bndy,var);

	// Result array:					Package ID,  Column Number,  Values
	Array<double,3> results_all(Range(1,NoOfPackages),Range(1,N_variables),Range(1,N_slave));
	Array<double,2> recieve(Range(1,N_variables),Range(1,N_slave));
	Array<double,2> slice;
	tag = 1;	// first Package

	// Z array
	// if(use_inputPointsFile) Z_values is not used. Slave works on package 'tag', which already determines the indices of the points to work on in array "inputPoints"
	Array<double,1> Z_values(Range(1,PAR.NZ));
	for(i=1;i<=PAR.NZ;i++) Z_values(i) = PAR.Zmin + (i-1)*dz;

	int sent_packages = 0;
	int recieve_packages = 0;

	#pragma omp parallel shared(results_all,Z_values,sent_packages,recieve_packages,write_memory) private(i,tag,N_slave_use) num_threads(2)
	{
		#pragma omp sections nowait
		{
		#pragma omp section	//-------- Master Thread: controlles comunication ----------------------------------------------------------------------------------------------------------------------
		{
			//#pragma omp barrier	// Syncronize with Slave Thread
			MPI::COMM_WORLD.Barrier();	// Master waits for Slaves

			//ofs3 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;
			ofs3 << "Start Tracer for " << N << " points ... " << endl;

			// Send initial Package to Slaves 
			int workingNodes = 0;
			for(i=1;i<mpi_size;i++)
			{
				if(sent_packages < NoOfPackages)	
				{
					#pragma omp critical
					{
						sent_packages += 1;
						tag = sent_packages;
					}
					send_Z_limits(1) = Z_values((tag-1)*NZ_slave+1);	// Zmin_slave
					send_Z_limits(2) = Z_values(tag*NZ_slave);	// Zmax_slave

					MPI::COMM_WORLD.Send(send_Z_limits.dataFirst(),2,MPI::DOUBLE,i,tag);
					workingNodes += 1;

					ofs3 << "Send Package No.: " << tag << endl;
				}
				else	// more Nodes than Packages -> Send termination signal: tag = 0
				{
					MPI::COMM_WORLD.Send(send_Z_limits.dataFirst(),0,MPI::DOUBLE,i,0);
					ofs3 << "Send termination signal to Node: " << i << endl;
				}
			} // end for(i=1;i<mpi_size;i++)

			// Recieve results and send more
			while(workingNodes > 0)	// workingNodes > 0: Slave still working -> MPI:Revc needed		workingNodes == 0: all Slaves recieved termination signal
			{
				// Recieve Result
				MPI::COMM_WORLD.Recv(recieve.dataFirst(),N_variables*N_slave,MPI::DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,status);
				sender = status.Get_source();
				tag = status.Get_tag();
				ofs3 << "Recieve from Node: " << sender << " Package: " << tag << endl;

				#pragma omp critical
				{
					recieve_packages += 1;
					write_memory(tag) = 1;
					ofs3 << "------------------------------------ Progress: " << recieve_packages << " of " << NoOfPackages << " completed" <<  endl;
				}

				slice.reference(results_all(tag,all,all));
				slice = recieve;

				// Send new Package
				if(sent_packages < NoOfPackages)	
				{
					#pragma omp critical
					{
						sent_packages += 1;
						tag = sent_packages;
					}

					send_Z_limits(1) = Z_values((tag-1)*NZ_slave+1);	// Zmin_slave
					send_Z_limits(2) = Z_values(tag*NZ_slave);	// Zmax_slave

					MPI::COMM_WORLD.Send(send_Z_limits.dataFirst(),2,MPI::DOUBLE,sender,tag);
					ofs3 << "Send again to Node: " << sender << " Package No.: " << tag << endl;
				}
				else	// No Packages left -> Send termination signal: tag = 0
				{
					MPI::COMM_WORLD.Send(send_Z_limits.dataFirst(),0,MPI::DOUBLE,sender,0);
					workingNodes -= 1;
					ofs3 << "Send termination signal to Node: " << sender << endl;
				}

				// Set write_max
				while((write_memory(write_max+1) == 1) && (write_max < NoOfPackages)) write_max++;

				// Write Output to file
				N_slave_use = N_slave;
				for(i=write_last+1; i<=write_max; i++)
				{
					if(i == NoOfPackages) N_slave_use = N_slave_last;	// size of last package may be smaller
					for(j=1;j<=N_slave_use;j++)
					{
						out << results_all(i,1,j);
						for(int k=2;k<=N_variables;k++) out << "\t" << results_all(i,k,j);
						out << endl;
					}
					write_memory(i) = 2;
				}
				write_last = write_max;

			} // end while(recieve_packages <= NoOfPackages)
		} // end omp section: Master thread

		#pragma omp section	//------- Slave Thread on Master Node: does same calculations as Salve Nodes ----------------------------------------------------------------------------------
		{
			// Prepare Perturbation
			prepare_common_perturbations(EQD,PAR,mpi_rank,siestafile,xpandfile,islandfile);
			prep_perturbation(EQD,PAR,mpi_rank);

			//#pragma omp barrier	// Syncronize with Master Thread

			cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
			cout << "Start Tracer for " << N << " points ... " << endl;
			ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;

			// Start working...
			while(1)	
			{
				// Recieve initial conditions to calculate
				#pragma omp critical
				{
					sent_packages += 1;
					tag = sent_packages;
				}
				if(tag > NoOfPackages) break;	// Still work to do?  ->  tag > NoOfPackages: NO, stop working
				if(tag == NoOfPackages) N_slave_use = N_slave_last;	// size of last package may be smaller
				else N_slave_use = N_slave;		// this is a private variable so it is different from the one outside the omp section and needs initialization
				ofs2 << endl;
				ofs2 << "Node: " << mpi_rank << " works on Package: " << tag << endl;
				ofs3 << "Node: " << mpi_rank << " works on Package: " << tag << endl;

				Zmin_slave = Z_values((tag-1)*NZ_slave+1);
				Zmax_slave = Z_values(tag*NZ_slave);

				ofs2 << "Start Tracer for " << N_slave_use << " points ... " << endl;
				for(i=1;i<=N_slave_use;i++)
				{
					// Set and store initial condition
					if(PAR.create_flag == 6)	// creates regular grid from theta and phi at psi = const.
					{
						FLT.set_surface(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave);	// here Rmin = thetamin and Zmin = phimin; max respectively
						results_all(tag,1,i) = FLT.phi;
						results_all(tag,2,i) = FLT.theta;
						xtmp = FLT.R;
						ytmp = FLT.Z;
					}
					else if(PAR.create_flag == 3)	// creates regular grid from theta and psi
					{
						FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave,2);	// here Rmin = psimin and Zmin = thetamin; max respectively
#ifdef USE_SIESTA
						if(PAR.response_field == -2)
						{
							double s,u;
							SIES.get_su(FLT.R, FLT.phi/rTOd, FLT.Z, s, u);
							results_all(tag,1,i) = u;
							results_all(tag,2,i) = s;
							xtmp = FLT.R;
							ytmp = FLT.Z;
						}
						else
						{
							results_all(tag,1,i) = FLT.theta;
							results_all(tag,2,i) = FLT.psi;
							xtmp = FLT.R;
							ytmp = FLT.Z;
						}
#else
						results_all(tag,1,i) = FLT.theta;
						results_all(tag,2,i) = FLT.psi;
						xtmp = FLT.R;
						ytmp = FLT.Z;
#endif
					}
					else						// creates regular grid from R and Z
					{
						FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave);	// matlab requires R to vary first
						results_all(tag,1,i) = FLT.R;
						results_all(tag,2,i) = FLT.Z;
						xtmp = FLT.theta;
						ytmp = FLT.psi;
					}

					// Overwrite previous settings in case a file with initial points is used
					if(use_inputPointsFile)
					{
						FLT.R = inputPoints((tag-1)*N_slave + i, 1);
						FLT.phi = inputPoints((tag-1)*N_slave + i, 2);
						FLT.Z = inputPoints((tag-1)*N_slave + i, 3);
						FLT.get_psi(FLT.R,FLT.Z,FLT.psi);
						FLT.theta = polar_phi(FLT.R-EQD.RmAxis, FLT.Z-EQD.ZmAxis);
						results_all(tag,1,i) = FLT.R;
						results_all(tag,2,i) = FLT.Z;
						xtmp = FLT.phi;
						ytmp = FLT.psi;
					}

					// Integration terminates outside of boundary box
					skip_connect = 0;
					if(outofBndy(FLT.phi,FLT.R,FLT.Z,EQD) == true)
					{
						ntor = 0;
						length = 0;
						psimin = 10;
						psimax = 0;
						psiav = 0;
						skip_connect = 1;
					}

					// Spare the calculation of the interior
					if(spare_interior == 1 && FLT.psi <= psi_interior_limit && FLT.Z > -1.25)
					{
						ntor = 2*PAR.itt;
						length = 4000.0;
						psimin = FLT.psi;
						psimax = FLT.psi;
						psiav = FLT.psi;
						skip_connect = 1;
					}

					chk = pitch_angles(FLT.R, FLT.Z, FLT.phi/rTOd, pitch, yaw, EQD, PAR);
					if(chk == -1) {pitch = 0; yaw = 0;}

					// Follow fieldline to walls
					if(skip_connect == 0) chk = FLT.connect(ntor,length,psimin,psimax,psiav,PAR.itt,PAR.MapDirection);

					// Store results
					results_all(tag,3,i) = ntor;
					results_all(tag,4,i) = length/1000.0;
					results_all(tag,5,i) = psimin;
					results_all(tag,6,i) = psimax;
					results_all(tag,7,i) = psiav;
					results_all(tag,8,i) = pitch;
					results_all(tag,9,i) = yaw;
					results_all(tag,10,i) = xtmp;
					results_all(tag,11,i) = ytmp;

					if(i%100==0) ofs2 << "Trax: " << i << endl;
				} // end for

				#pragma omp critical
				{
					recieve_packages += 1;
					write_memory(tag) = 1;
					ofs3 << "------------------------------------ Progress: " << recieve_packages << " of " << NoOfPackages << " completed" <<  endl;
				}

			}// end while
		}// end omp section: Slave Thread
		} // end omp sections
	} // end omp parallel

	// Write remaining Output to file
	N_slave_use = N_slave;
	for(i=write_last+1;i<=NoOfPackages;i++)
	{
		if(i == NoOfPackages) N_slave_use = N_slave_last;	// size of last package may be smaller
		for(j=1;j<=N_slave_use;j++)
		{
			out << results_all(i,1,j);
			for(int k=2;k<=N_variables;k++) out << "\t" << results_all(i,k,j);
			out << endl;
		}
	}
} // end Master
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	// Prepare Perturbation
	prepare_common_perturbations(EQD,PAR,mpi_rank,siestafile,xpandfile,islandfile);
	prep_perturbation(EQD,PAR,mpi_rank);

	MPI::COMM_WORLD.Barrier();	// Syncronize with Master

	ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;

	// Result array for Slave
	Array<double,2> results(Range(1,N_variables),Range(1,N_slave));

	// Start working...
	while(1)	
	{
		// Recieve initial conditions to calculate
		MPI::COMM_WORLD.Recv(send_Z_limits.dataFirst(),2,MPI::DOUBLE,0,MPI_ANY_TAG,status);

		tag = status.Get_tag();
		ofs2 << endl;
		ofs2 << "Node: " << mpi_rank << " works on Package: " << tag << endl;
		if(tag == 0) break;	// Still work to do?  ->  tag = 0: NO, stop working
		if(tag == NoOfPackages) N_slave_use = N_slave_last;	// size of last package may be smaller

		Zmin_slave = send_Z_limits(1);
		Zmax_slave = send_Z_limits(2);

		ofs2 << "Start Tracer for " << N_slave_use << " points ... " << endl;
		for(i=1;i<=N_slave_use;i++)
		{
			// Set and store initial condition
			if(PAR.create_flag == 6)	// creates regular grid from theta and phi at psi = const.
			{
				FLT.set_surface(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave);	// here Rmin = thetamin and Zmin = phimin; max respectively
				results(1,i) = FLT.phi;
				results(2,i) = FLT.theta;
				xtmp = FLT.R;
				ytmp = FLT.Z;
			}
			else if(PAR.create_flag == 3)	// creates regular grid from theta and psi
			{
				FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave,2);	// here Rmin = psimin and Zmin = thetamin; max respectively
#ifdef USE_SIESTA
				if(PAR.response_field == -2)
				{
					double s,u;
					SIES.get_su(FLT.R, FLT.phi/rTOd, FLT.Z, s, u);
					results(1,i) = u;
					results(2,i) = s;
					xtmp = FLT.R;
					ytmp = FLT.Z;
				}
				else
				{
					results(1,i) = FLT.theta;
					results(2,i) = FLT.psi;
					xtmp = FLT.R;
					ytmp = FLT.Z;
				}
#else
				results(1,i) = FLT.theta;
				results(2,i) = FLT.psi;
				xtmp = FLT.R;
				ytmp = FLT.Z;
#endif
			}
			else						// creates regular grid from R and Z
			{
				FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave);	// matlab requires R to vary first
				results(1,i) = FLT.R;
				results(2,i) = FLT.Z;
				xtmp = FLT.theta;
				ytmp = FLT.psi;
			}

			// Overwrite previous settings in case a file with initial points is used
			if(use_inputPointsFile)
			{
				FLT.R = inputPoints((tag-1)*N_slave + i, 1);
				FLT.phi = inputPoints((tag-1)*N_slave + i, 2);
				FLT.Z = inputPoints((tag-1)*N_slave + i, 3);
				FLT.get_psi(FLT.R,FLT.Z,FLT.psi);
				FLT.theta = polar_phi(FLT.R-EQD.RmAxis, FLT.Z-EQD.ZmAxis);
				results(1,i) = FLT.R;
				results(2,i) = FLT.Z;
				xtmp = FLT.phi;
				ytmp = FLT.psi;
			}

			// Integration terminates outside of boundary box
			skip_connect = 0;
			if(outofBndy(FLT.phi,FLT.R,FLT.Z,EQD) == true)
			{
				ntor = 0;
				length = 0;
				psimin = 10;
				psimax = 0;
				psiav = 0;
				skip_connect = 1;
			}

			// Spare the calculation of the interior
			if(spare_interior == 1 && FLT.psi <= psi_interior_limit && FLT.Z > -1.25)
			{
				ntor = 2*PAR.itt;
				length = 4000.0;
				psimin = FLT.psi;
				psimax = FLT.psi;
				psiav = FLT.psi;
				skip_connect = 1;
			}

			chk = pitch_angles(FLT.R, FLT.Z, FLT.phi/rTOd, pitch, yaw, EQD, PAR);
			if(chk == -1) {pitch = 0; yaw = 0;}

			// Follow fieldline to walls
			if(skip_connect == 0) chk = FLT.connect(ntor,length,psimin,psimax,psiav,PAR.itt,PAR.MapDirection);

			// Store results
			results(3,i) = ntor;
			results(4,i) = length/1000.0;
			results(5,i) = psimin;
			results(6,i) = psimax;
			results(7,i) = psiav;
			results(8,i) = pitch;
			results(9,i) = yaw;
			results(10,i) = xtmp;
			results(11,i) = ytmp;

			if(i%100==0) ofs2 << "Trax: " << i << endl;
		} // end for

		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),N_variables*N_slave,MPI::DOUBLE,0,tag);

	}// end while
} // end Slaves
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double now2 = zeit();
if(mpi_rank < 1) 
{
	cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
	ofs3 << "Program terminates normally, Time: " << now2-now  << " s" << endl;
}
ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

// MPI finalize
MPI::Finalize();

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

int pitch_angles(double R, double Z, double phi, double& pitch, double& yaw, EFIT& EQD, IO& PAR)
{
int chk;
double B_R,B_Z,B_phi,B_pol,B_r;
const double theta = polar_phi(R - EQD.RmAxis, Z - EQD.ZmAxis);
const double sint = sin(theta);
const double cost = cos(theta);

chk = getBfield(R, Z, phi, B_R, B_Z, B_phi, EQD, PAR);
if (chk == -1) return -1;

B_r = B_R*cost + B_Z*sint;
B_pol = -B_R*sint + B_Z*cost;

pitch = atan(B_pol/B_phi);
yaw = atan(B_r/B_phi);

return 0;
}

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

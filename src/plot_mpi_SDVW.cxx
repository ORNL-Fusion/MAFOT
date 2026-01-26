// Program calculates Poincare Plot for particle-drift with time dependent perturbations
// Fortran subroutines for Perturbation are used
// A.Wingen						7.6.11

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	Poincare particle drift data file
//			log-file

// uses Parallel computation by open-mpi and openmp
// IMPORTANT: Node 0 starts two threads using openmp. The Master-thread contols communications with MPI-Slaves. The Slave-thread on Node 0 does the same calculations as the MPI-Slaves
// Launch Syntax:				mpirun -np [NoOfNodes] -hostfile [FileName] dtlaminar_mpi [Parameterfile] [Praefix(optional)] &
// alternative Launch Syntax:	mpirun -np [NoOfNodes] -host [HostNames separated by Commata] dtlaminar_mpi [Parameterfile] [Praefix(optional)] &

// Define
//--------
#define USE_MPI
#if defined(ITER)
	#define program_name "iterplot_mpi"
#elif defined(NSTX)
	#define program_name "nstxplot_mpi"
#elif defined(MAST)
	#define program_name "mastplot_mpi"
#elif defined(CMOD)
	#define program_name "cmodplot_mpi"
#elif defined(TCABR)
	#define program_name "tcabrplot_mpi"
#elif defined(ANYM)
	#define program_name "anyplot_mpi"
#else
	#define program_name "dtplot_mpi"
#endif

// Include
//--------
#ifdef USE_MPICH
	#include <mpi.h>
#else
	#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#endif
#include <mafot.hxx>
#include <omp.h>
#include <unistd.h>
#include <vector>
#include <sstream>

// namespaces
//-----------
using namespace blitz;

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
int i,j,n,chk;
EFIT EQD;
Range all = Range::all();
double s,u;

int tag,sender;
int Nmin_slave,Nmax_slave;
Array<int,1> send_N_limits(Range(1,2));
MPI::Status status;

// Use system time as seed(=idum) for random numbers
double now=zeit();

// defaults
LA_STRING woutfile = "wout.nc";
LA_STRING xpandfile = "None";
LA_STRING siestafile = "siesta.dat";
LA_STRING islandfile = "fakeIslands.in";
LA_STRING ErProfileFile = "None";
bool use_ErProfile = false;
bool use_Tprofile = false;
LA_STRING TprofileFile = "prof_t.dat";
LA_STRING NprofileFile = "prof_n.dat";
double f = 0;  // ratio of impurity to hydrogen ions
double zbar = 2;  // average over impurity ion charge states
double rc = 1;
double mc = 1;
bool use_collision = false;
// start positions from file
LA_STRING startposfile = "None";
bool use_startpos = false;
std::vector<double> startR, startZ, startPhi;

// Command line input parsing
int c;
opterr = 0;
bool force_overwrite = false;
bool print_points = false;
while ((c = getopt(argc, argv, "hX:V:S:I:E:T:i:sP:yp")) != -1)
switch (c)
{
case 'h':
	if(mpi_rank < 1)
	{
		cout << "usage: mpirun -n <cores> dtplot_mpi [-h] [-i dpinit] [-s] [-E ErProfile] [-I island] [-S siesta] [-T Tprofile] [-V wout] [-X xpand] file [tag]" << endl << endl;
		cout << "Calculate a Poincare plot." << endl << endl;
		cout << "positional arguments:" << endl;
		cout << "  file          Control file (starts with '_')" << endl;
		cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
		cout << endl << "optional arguments:" << endl;
		cout << "  -h            show this help message and exit" << endl;
		cout << "  -i            change integrator step size; default is 1.0" << endl;
		cout << "  -s            use a simple boundary box instead of real wall for field line termination" << endl;
		cout << "  -E            use electric field with particle drifts. Filename of Er(psi) profile." << endl;
		cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
		cout << "  -S            filename for SIESTA; default, see below" << endl;
		cout << "  -T            use temperature profile with particle drifts. Filename of T(psi) profile." << endl;
		cout << "  -V            filename for VMEC; default, see below" << endl;
		 cout << "  -X            filename for XPAND; default is None" << endl;
		 cout << "  -p            print points passed to solver (from runtime arrays)" << endl;
		 cout << "  -y            force overwrite of output file without prompting" << endl;
		cout << endl << "Examples:" << endl;
		cout << "  mpirun -n 4 dtplot_mpi _plot.dat blabla" << endl;
		cout << endl << "Infos:" << endl;
		cout << "  To use B-field from M3DC1, set response_field >= 0, and provide file in cwd:" << endl;
		cout << "    m3dc1sup.in    ->  location and scale factor for M3DC1 output C1.h5" << endl;
		cout << "  To use B-field from XPAND, set response_field = -3, and provide files in cwd:" << endl;
		cout << "    xpand.dat      ->  B-field on 3D grid from XPAND; use option -X to specify a filename (default is None -> inside VMEC only)" << endl;
		cout << "    wout.nc        ->  VMEC output; use option -V to specify other filename" << endl;
		cout << "  To use B-field from VMEC, inside only (no xpand file given), set response_field = -3, and provide file in cwd:" << endl;
		cout << "    wout.nc        ->  VMEC output; use option -V to specify other filename" << endl;
		cout << "  To use B-field from SIESTA, set response_field = -2, and provide file in cwd:" << endl;
		cout << "    siesta.dat     ->  B-field on 3D grid; use option -S to specify other filename" << endl;
		cout << "  To use B-field for mock-up islands, set response_field = -10, and provide file in cwd:" << endl;
		cout << "    fakeIslands.in ->  each line gives: Amplitude, pol. mode m, tor. mode n, phase [rad]" << endl;
		cout << "                       use option -I to specify other filename" << endl;
		cout << endl << "Current MAFOT version is: " << MAFOT_VERSION << endl;
	}
	MPI::Finalize();
	return 0;
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
case 'E':
	ErProfileFile = optarg;
	use_ErProfile = true;
	break;
case 'T':
	TprofileFile = optarg;
	use_Tprofile = true;
	break;
case 's':
	simpleBndy = 1;
	break;
case 'i':
	dpinit = atof(optarg);
	ilt = int(360.0/dpinit + 0.5);
	break;
case 'P':
	startposfile = optarg;
	use_startpos = true;
	break;
case 'y':
	force_overwrite = true;
	break;
case 'p':
	print_points = true;
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
else	// No Input: Abort
{
	if(mpi_rank < 1) cout << "No Input files -> Abort!" << endl;
	EXIT;
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";
//if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes and restart." << endl; EXIT;}

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + "_Node" + LA_STRING(mpi_rank) + ".dat");
ofs2.precision(16);
ofstream ofs3;

// Output
LA_STRING filenameout = "plot" + praefix + ".dat";
if(mpi_rank < 1 && !force_overwrite) outputtest(filenameout);
MPI::COMM_WORLD.Barrier();	// All Nodes wait for Master

// Read parameter file
if(mpi_rank < 1) cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10,mpi_rank);
if (PAR.NZ > 0) PAR.N = PAR.NR*PAR.NZ;	// fix grid issue for NZ > 1
if (PAR.MapDirection == 0) PAR.MapDirection = 1;

// If requested, read start-position file on master and broadcast to all ranks
// Synchronize startpos option and filename across all ranks so every rank
// participates in the same MPI collectives even when argument parsing
// or mpirun argument ordering differs.
{
	int use_sp = use_startpos ? 1 : 0;
	MPI_Bcast(&use_sp, 1, MPI_INT, 0, MPI_COMM_WORLD);
	use_startpos = (use_sp != 0);
	if(use_startpos) {
		// broadcast filename length and content
		int len = 0;
		if(mpi_rank == 0) len = (int)strlen((const char*)startposfile);
		MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(len > 0) {
			char *buf = new char[len+1];
			if(mpi_rank == 0) memcpy(buf, (const char*)startposfile, len);
			MPI_Bcast(buf, len, MPI_CHAR, 0, MPI_COMM_WORLD);
			buf[len] = '\0';
			if(mpi_rank != 0) startposfile = LA_STRING(buf);
			delete [] buf;
		}
	}
}
if(use_startpos)
{
	int npos = 0;
	if(mpi_rank < 1)
	{
		std::ifstream ifs((const char*)startposfile);
		if (!ifs) {
			std::cerr << "Cannot open start-position file: " << startposfile << " -> Abort!\n";
			MPI_Abort(MPI_COMM_WORLD, 1);   // abort all ranks
		}

		std::string line;
		int lineno = 0;
		while(std::getline(ifs, line)) {
			++lineno;
			size_t first = line.find_first_not_of(" \t\r\n");
			if(first == std::string::npos) continue; // blank line
			if(line[first] == '#') continue; // comment/header line
			std::istringstream iss(line);
			double a, b;
			if(!(iss >> a >> b)) {
				std::cerr << "Warning: could not parse R Z on line " << lineno << " of " << startposfile << " -> skipped\n";
				continue;
			}
			startR.push_back(a);
			startZ.push_back(b);
			startPhi.push_back(PAR.phistart);
		}
		npos = static_cast<int>(startR.size());
	}
	// Broadcast number of positions to all ranks
	{
		int npos_b = npos;
		MPI_Bcast(&npos_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
		npos = npos_b;
	}

	// Do not broadcast full arrays to all ranks here. We'll send package-sized coordinate arrays per-package to slaves.

	// Update number of points
	PAR.N = npos;
	// NOTE: PAR.N is overridden by the start-position file and is broadcast to all ranks above.
	// Make sure we read at least one point and update the PV metadata used for logging.
	if(mpi_rank < 1)
	{
		if(startR.empty())
		{
			ofs2 << "Error: start-position file contains no points -> aborting" << endl;
			cerr << "Error: start-position file contains no points -> aborting" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	// Keep the parameter metadata consistent
	PAR.pv[1].name = "Points"; PAR.pv[1].wert = PAR.N;
}
	// One-time confirmation print on master that start positions were read
	if(use_startpos && mpi_rank == 0 && print_points)
	{
		cout << "Startpos file read: N=" << PAR.N << ", first idx=1 R=" << startR.front() << " Z=" << startZ.front() << " last idx=" << PAR.N << " R=" << startR.back() << " Z=" << startZ.back() << endl;
	}

// Read EFIT-data
double Raxis = 0, Zaxis = 0;
#ifdef USE_XFIELD
if(PAR.response_field == -3)
{
	if(mpi_rank < 1) cout << "Read VMEC file" << endl;
	ofs2 << "Read VMEC file" << endl;
	vmec.read(woutfile);
	vmec.n0only = true;
	vmec.get_axis(PAR.phistart/rTOd,Raxis,Zaxis);	// need axisymmetric axis only
	vmec.n0only = false;
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

EQD.ReadData(EQD.gFile,Raxis,Zaxis);
if(mpi_rank < 1) cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.gFile << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.gFile << endl;

// Read E-field data
if(use_ErProfile)
{
	EQD.ReadEfield(ErProfileFile);
	PAR.useErProfile = 1;
	if(mpi_rank < 1) cout << "Read Er profile: " << ErProfileFile << endl;
	ofs2 << "Read Er profile: " << ErProfileFile << endl;
}

// Read T-profile data
if(use_Tprofile)
{
	EQD.ReadTprofile(TprofileFile);
	PAR.useTprofile = 1;
	if(mpi_rank < 1) cout << "Ignoring Ekin, read T profile: " << TprofileFile << endl;
	ofs2 << "Ignoring Ekin, read T profile: " << TprofileFile << endl;
}

// Set starting parameters
int N = PAR.N;
int N_slave = 1;    // Number of field lines per package (manual control)
int NoOfPackages = int(PAR.N / N_slave);
int N_slave_last = 0; // size of last package if PAR.N not divisible by N_slave
// If not exact, work on one more package and compute last-package size
if(N_slave * NoOfPackages < N) {
	N_slave_last = N - N_slave * NoOfPackages;
	NoOfPackages += 1;
}

// Prepare collisions
COLLISION COL;
if (use_collision) COL.init(TprofileFile, NprofileFile, f, zbar, PAR.Zq, PAR.Mass, rc, mc, mpi_rank);
// Prepare particles
PARTICLE FLT(EQD,PAR,COL,mpi_rank);

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
	PAR.pv[1].name = "Points";			PAR.pv[1].wert = PAR.N;
	PAR.pv[6].name = "phistart";		PAR.pv[6].wert = PAR.phistart;
	PAR.pv[7].name = "MapDirection";	PAR.pv[7].wert = PAR.MapDirection;
	PAR.pv[8].name = "Ekin";			PAR.pv[8].wert = PAR.Ekin;
	PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;

	switch(PAR.create_flag)
	{
	case 2:
		PAR.pv[2].name = "tmin";			PAR.pv[2].wert = PAR.tmin;
		PAR.pv[3].name = "tmax";			PAR.pv[3].wert = PAR.tmax;
		PAR.pv[4].name = "phimin";			PAR.pv[4].wert = PAR.phistart;
		PAR.pv[5].name = "phimax";			PAR.pv[5].wert = PAR.phistart;
		break;
	case 3: case 4:
		PAR.pv[2].name = "psimin";			PAR.pv[2].wert = PAR.rmin;
		PAR.pv[3].name = "psimax";			PAR.pv[3].wert = PAR.rmax;
		PAR.pv[4].name = "thmin";			PAR.pv[4].wert = PAR.thmin;
		PAR.pv[5].name = "thmax";			PAR.pv[5].wert = PAR.thmax;
		break;
	case 5:
		PAR.pv[2].name = "Rmin";			PAR.pv[2].wert = PAR.Rmin;
		PAR.pv[3].name = "Rmax";			PAR.pv[3].wert = PAR.Rmax;
		PAR.pv[4].name = "Zmin";			PAR.pv[4].wert = PAR.Zmin;
		PAR.pv[5].name = "Zmax";			PAR.pv[5].wert = PAR.Zmax;
		break;
	default:
		PAR.pv[2].name = "rmin";			PAR.pv[2].wert = PAR.rmin;
		PAR.pv[3].name = "rmax";			PAR.pv[3].wert = PAR.rmax;
		PAR.pv[4].name = "thmin";			PAR.pv[4].wert = PAR.thmin;
		PAR.pv[5].name = "thmax";			PAR.pv[5].wert = PAR.thmax;
		break;
	}

	// Output
	ofstream out(filenameout);
	out.precision(16);
	vector<LA_STRING> var(6);
	var[0] = "theta[rad]"; var[1] = "r[m]"; var[2] = "phi[deg]"; var[3] = "psi"; var[4] = "R[m]"; var[5] = "Z[m]";
	if(PAR.response_field == -2) {var[0] = "u"; var[3] = "s";}
	PAR.writeiodata(out,bndy,var);

	cout << "Create (0=r-grid, 1=r-random, 2=target, 3=psi-grid, 4=psi-random, 5=R-grid): " << PAR.create_flag << endl;
	if(PAR.create_flag==2) cout << "Target: " << PAR.which_target_plate << endl;
	ofs3 << "Create (0=r-grid, 1=r-random, 2=target, 3=psi-grid, 4=psi-random, 5=R-grid): " << PAR.create_flag << endl;
	if(PAR.create_flag==2) ofs3 << "Target: " << PAR.which_target_plate << endl;

	cout << "Helicity = " << EQD.helicity << endl;
	ofs3 << "Helicity = " << EQD.helicity << endl;

	// Result array:	 Column Number,  Values
	Array<double,2> results(Range(1,6),Range(1,N_slave*PAR.itt));
	Array<double,2> recieve(Range(1,6),Range(1,N_slave*PAR.itt));
	Array<double,2> slice;
	tag = 1;	// first Package

	// number of field lines array
	Array<int,1> N_values(Range(1,PAR.N));
	for(i=1;i<=PAR.N;i++) N_values(i) = i;

	int sent_packages = 0;
	int recieve_packages = 0;

	#pragma omp parallel shared(N_values,sent_packages,recieve_packages) private(i,tag) num_threads(2)
	{
		#pragma omp sections nowait
		{
		#pragma omp section	//-------- Master Thread: controlles comunication ----------------------------------------------------------------------------------------------------------------------
		{
			//#pragma omp barrier	// Syncronize with Slave Thread
			MPI::COMM_WORLD.Barrier();	// Master waits for Slaves

			ofs3 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
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
					send_N_limits(1) = N_values((tag-1)*N_slave+1);	// Nmin_slave
					// compute safe upper index for this package
					{
						int n2_idx = tag * N_slave;
						if(n2_idx > N) n2_idx = N;
						send_N_limits(2) = N_values(n2_idx); // Nmax_slave
					}

					MPI::COMM_WORLD.Send(send_N_limits.dataFirst(),2,MPI::INTEGER,i,tag);
					// If using start-position file, send the small coordinate chunk for this package
					if(use_startpos)
					{
						int n1 = send_N_limits(1);
						int n2 = send_N_limits(2);
						int cnt = n2 - n1 + 1;
						// master (rank 0) holds full arrays startR/startZ
						if(print_points) {
							// concise summary to master log and stdout
							ofs3 << "Master sending package " << tag << " indices [" << n1 << "," << n2 << "] count=" << cnt;
							if(cnt <= 10) {
								ofs3 << " points:";
								for(int pi = n1; pi <= n2; ++pi) ofs3 << " idx=" << pi << ":R=" << startR[pi-1] << ",Z=" << startZ[pi-1] << ";";
							} else {
								ofs3 << " sample: idx=" << n1 << ":R=" << startR[n1-1] << ",Z=" << startZ[n1-1] << " ... idx=" << n2 << ":R=" << startR[n2-1] << ",Z=" << startZ[n2-1];
							}
							ofs3 << "\n";
							if(mpi_rank < 1) cout << "Master sending package " << tag << " indices [" << n1 << "," << n2 << "] count=" << cnt << endl;
						}
						MPI::COMM_WORLD.Send(&startR[n1-1], cnt, MPI::DOUBLE, i, tag+10000);
						MPI::COMM_WORLD.Send(&startZ[n1-1], cnt, MPI::DOUBLE, i, tag+20000);
					}
					workingNodes += 1;

					ofs3 << "Send Package No.: " << tag << endl;
				}
				else	// more Nodes than Packages -> Send termination signal: tag = 0
				{
					MPI::COMM_WORLD.Send(send_N_limits.dataFirst(),0,MPI::INTEGER,i,0);
					ofs3 << "Send termination signal to Node: " << i << endl;
				}
			} // end for(i=1;i<mpi_size;i++)

			// Recieve results and send more
			while(workingNodes > 0)	// workingNodes > 0: Slave still working -> MPI:Revc needed		workingNodes == 0: all Slaves recieved termination signal
			{
				// Recieve Result
				MPI::COMM_WORLD.Recv(recieve.dataFirst(),6*N_slave*PAR.itt,MPI::DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,status);
				sender = status.Get_source();
				tag = status.Get_tag();
				ofs3 << "Recieve from Node: " << sender << " Package: " << tag << endl;

				// write results to file
				#pragma omp critical
				{
					for(j=1;j<=N_slave*PAR.itt;j++)	if(recieve(4,j) > 0) out << recieve(1,j) << "\t" << recieve(2,j) << "\t" << recieve(3,j) << "\t" << recieve(4,j) << "\t" << recieve(5,j) << "\t" << recieve(6,j) << endl;
					recieve_packages += 1;
					ofs3 << "------------------------------------ Progress: " << recieve_packages << " of " << NoOfPackages << " completed" <<  endl;
				}

				// Send new Package
				if(sent_packages < NoOfPackages)
				{
					#pragma omp critical
					{
						sent_packages += 1;
						tag = sent_packages;
					}

					send_N_limits(1) = N_values((tag-1)*N_slave+1);	// Nmin_slave
					// compute safe upper index for this package
					{
						int n2_idx = tag * N_slave;
						if(n2_idx > N) n2_idx = N;
						send_N_limits(2) = N_values(n2_idx); // Nmax_slave
					}

					MPI::COMM_WORLD.Send(send_N_limits.dataFirst(),2,MPI::INTEGER,sender,tag);
					// send coordinates for this package
					if(use_startpos)
					{
						int n1 = send_N_limits(1);
						int n2 = send_N_limits(2);
						int cnt = n2 - n1 + 1;
						MPI::COMM_WORLD.Send(&startR[n1-1], cnt, MPI::DOUBLE, sender, tag+10000);
						MPI::COMM_WORLD.Send(&startZ[n1-1], cnt, MPI::DOUBLE, sender, tag+20000);
					}
					ofs3 << "Send again to Node: " << sender << " Package No.: " << tag << endl;
				}
				else	// No Packages left -> Send termination signal: tag = 0
				{
					MPI::COMM_WORLD.Send(send_N_limits.dataFirst(),0,MPI::INTEGER,sender,0);
					workingNodes -= 1;
					ofs3 << "Send termination signal to Node: " << sender << endl;
				}
			} // end while(recieve_packages <= NoOfPackages)
		} // end omp section: Master thread

		#pragma omp section	//------- Slave Thread on Master Node: does same calculations as Salve Nodes ----------------------------------------------------------------------------------
		{
			// Prepare Perturbation
			prepare_common_perturbations(EQD,PAR,mpi_rank,siestafile,xpandfile,islandfile);
			prep_perturbation(EQD,PAR,mpi_rank);

			// each process gets a different seed
			long idum=long(now) + 1000*mpi_rank;

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
				ofs2 << endl;
				ofs2 << "Node: " << mpi_rank << " works on Package: " << tag << endl;
				#pragma omp critical
				{
					ofs3 << "Node: " << mpi_rank << " works on Package: " << tag << endl;
				}

				Nmin_slave = N_values((tag-1)*N_slave+1);
				// compute safe upper index for this package
				{
					int n2_idx = tag * N_slave;
					if(n2_idx > N) n2_idx = N;
					Nmax_slave = N_values(n2_idx);
				}

				// Get Poincar� section
				ofs2 << "Start Tracer for " << N_slave << " points ... " << endl;
				// Print actual values that will be passed to the solver for this package (once, concise)
				if(print_points && use_startpos) {
					int p_n1 = Nmin_slave;
					int p_n2 = Nmax_slave;
					int p_cnt = p_n2 - p_n1 + 1;
					cout << "Master will pass to solver package " << tag << " indices [" << p_n1 << "," << p_n2 << "] count=" << p_cnt;
					if(p_cnt <= 10) {
						cout << " points:";
						for(int pi = p_n1; pi <= p_n2; ++pi) cout << " idx=" << pi << ":R=" << startR[pi-1] << ",Z=" << startZ[pi-1] << ";";
					} else {
						cout << " sample: idx=" << p_n1 << ":R=" << startR[p_n1-1] << " ... idx=" << p_n2 << ":R=" << startR[p_n2-1];
					}
					cout << endl;
					ofs2 << "Master will pass to solver package " << tag << " indices [" << p_n1 << "," << p_n2 << "] count=" << p_cnt << "\n";
				}
				for(n=Nmin_slave;n<=Nmax_slave;n++)
				{
					// Set initial values in FLT
					if(use_startpos)
					{
						int idx = n-1; // n is 1-based
						FLT.R = startR[idx];
						FLT.Z = startZ[idx];
						// phi is set from parameters, not the input file
						FLT.phi = PAR.phistart;
						FLT.get_psi(FLT.R, FLT.Z, FLT.psi, FLT.phi/rTOd);
						FLT.theta = FLT.get_theta();
						// reset trajectory counters
						FLT.Lc = 0; FLT.psimin = 10; FLT.psimax = 0; FLT.psiav = 0; FLT.steps = 0; FLT.dpsidLcav = 0;
						if(FLT.sigma != 0 && PAR.useTprofile == 1) { FLT.set_Energy(); }
					}
					else
					{
						switch(PAR.create_flag)
						{
							case 5: 		// grid from R, Z
								FLT.set(n,PAR.N,PAR.Rmin,PAR.Rmax,PAR.Zmin,PAR.Zmax,PAR.NZ,0);
								break;
							case 4: 		// random from psi, theta
								FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,2);
								break;
							case 3: 		// grid from psi, theta
								FLT.set(n,PAR.N,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,PAR.NZ,2);
								break;
							case 2: 		// grid on target from t, phi
								start_on_target(n,PAR.Nt,1,PAR.tmin,PAR.tmax,PAR.phistart,PAR.phistart,EQD,PAR,FLT);
								break;
							case 1: 		// random from r, theta
								FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,1);
								break;
							default: 	// grid from r, theta
								FLT.set(n,PAR.N,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,1,1);
								break;
						}
					}

					// use psi as control parameter and set to -1
					results(4,all) = -1;

					// Integrate
					for(i=1;i<=PAR.itt;i++)
					{
						chk = FLT.mapstep(PAR.MapDirection);
						if(chk<0) {ofs2 << "mapit: wall hit" << endl; break;}	// particle has left system

						if(fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) << endl;
						FLT.phi=PAR.MapDirection*i*dpinit*ilt + PAR.phistart;

						// Store results
#ifdef USE_SIESTA
						if(PAR.response_field == -2)
						{
							SIES.get_su(FLT.R, FLT.phi/rTOd, FLT.Z, s, u);
							results(1,(n-Nmin_slave)*PAR.itt + i) = u;
							results(4,(n-Nmin_slave)*PAR.itt + i) = s;
						}
						else
						{
							results(1,(n-Nmin_slave)*PAR.itt + i) = FLT.get_theta();
							results(4,(n-Nmin_slave)*PAR.itt + i) = FLT.psi;
						}
#else
						results(1,(n-Nmin_slave)*PAR.itt + i) = FLT.get_theta();
						results(4,(n-Nmin_slave)*PAR.itt + i) = FLT.psi;
#endif
						results(2,(n-Nmin_slave)*PAR.itt + i) = FLT.get_r();
						results(3,(n-Nmin_slave)*PAR.itt + i) = FLT.phi;
						results(5,(n-Nmin_slave)*PAR.itt + i) = FLT.R;
						results(6,(n-Nmin_slave)*PAR.itt + i) = FLT.Z;
					} // end for i
					ofs2 << "Trax: " << n << "\t" << "Steps: " << i-1 << endl;
				} // end for n

				#pragma omp critical
				{
					for(j=1;j<=N_slave*PAR.itt;j++)	if(results(4,j) > 0) out << results(1,j) << "\t" << results(2,j) << "\t" << results(3,j) << "\t" << results(4,j) << "\t" << results(5,j) << "\t" << results(6,j) << endl;
					recieve_packages += 1;
					ofs3 << "------------------------------------ Progress: " << recieve_packages << " of " << NoOfPackages << " completed" <<  endl;
				}
			} // end while
		} // end omp section: Slave Thread
		} // end omp sections
	} // end omp parallel
} // end Master
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	bool printed_confirmation = false;
	// Prepare Perturbation
	prepare_common_perturbations(EQD,PAR,mpi_rank,siestafile,xpandfile,islandfile);
	prep_perturbation(EQD,PAR,mpi_rank);

	// each process gets a different seed
	long idum=long(now) + 1000*mpi_rank;

	MPI::COMM_WORLD.Barrier();	// Syncronize with Master

	ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;

	// Result array for Slave
	Array<double,2> results(Range(1,6),Range(1,N_slave*PAR.itt));

	// Start working...
	while(1)
	{
		// Recieve initial conditions to calculate
		MPI::COMM_WORLD.Recv(send_N_limits.dataFirst(),2,MPI::INTEGER,0,MPI_ANY_TAG,status);

		tag = status.Get_tag();
		ofs2 << endl;
		ofs2 << "Node: " << mpi_rank << " works on Package: " << tag << endl;
		if(tag == 0) break;	// Still work to do?  ->  tag = 0: NO, stop working

		Nmin_slave = send_N_limits(1);
		Nmax_slave = send_N_limits(2);

		// If start positions come from file, receive the coordinate chunk for this package
		std::vector<double> pkgR, pkgZ;
		if(use_startpos && tag > 0)
		{
			int cnt = Nmax_slave - Nmin_slave + 1;
			pkgR.resize(cnt);
			pkgZ.resize(cnt);
			MPI::COMM_WORLD.Recv(pkgR.data(), cnt, MPI::DOUBLE, 0, tag+10000, status);
			MPI::COMM_WORLD.Recv(pkgZ.data(), cnt, MPI::DOUBLE, 0, tag+20000, status);
			if(print_points && !printed_confirmation) {
				// concise one-line confirmation only once
				cout << "Node " << mpi_rank << " received startpos package " << tag << " indices [" << Nmin_slave << "," << Nmax_slave << "] count=" << cnt << endl;
				ofs2 << "Node " << mpi_rank << " received startpos package " << tag << " indices [" << Nmin_slave << "," << Nmax_slave << "] count=" << cnt << "\n";
				printed_confirmation = true;
			}
		}

		ofs2 << "Start Tracer for " << N_slave << " points ... " << endl;

		// Get Poincar� section
		for(n=Nmin_slave;n<=Nmax_slave;n++)
		{
			// Set initial values in FLT
			if(use_startpos)
			{
				int pkgIdx = n - Nmin_slave; // index into pkgR/pkgZ
				FLT.R = pkgR[pkgIdx];
				FLT.Z = pkgZ[pkgIdx];
				// phi is set from parameters, not the input file
				FLT.phi = PAR.phistart;
				FLT.get_psi(FLT.R, FLT.Z, FLT.psi, FLT.phi/rTOd);
				FLT.theta = FLT.get_theta();
				// reset trajectory counters
				FLT.Lc = 0; FLT.psimin = 10; FLT.psimax = 0; FLT.psiav = 0; FLT.steps = 0; FLT.dpsidLcav = 0;
				if(FLT.sigma != 0 && PAR.useTprofile == 1) { FLT.set_Energy(); }
			}
			// Print actual package contents once (concise)
			if(print_points && use_startpos && n==Nmin_slave) {
				int cnt = Nmax_slave - Nmin_slave + 1;
				cout << "Node " << mpi_rank << " will pass to solver package " << tag << " indices [" << Nmin_slave << "," << Nmax_slave << "] count=" << cnt;
				if(cnt <= 10) {
					cout << " points:";
					for(int pi = 0; pi < cnt; ++pi) cout << " idx=" << (Nmin_slave + pi) << ":R=" << pkgR[pi] << ",Z=" << pkgZ[pi] << ";";
				} else {
					cout << " sample: idx=" << Nmin_slave << ":R=" << pkgR.front() << " ... idx=" << Nmax_slave << ":R=" << pkgR.back();
				}
				cout << endl;
			}
			else
			{
				switch(PAR.create_flag)
				{
					case 5: 		// grid from R, Z
						FLT.set(n,PAR.N,PAR.Rmin,PAR.Rmax,PAR.Zmin,PAR.Zmax,PAR.NZ,0);
						break;
					case 4: 		// random from psi, theta
						FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,2);
						break;
					case 3: 		// grid from psi, theta
						FLT.set(n,PAR.N,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,PAR.NZ,2);
						break;
					case 2: 		// grid on target from t, phi
						start_on_target(n,PAR.Nt,1,PAR.tmin,PAR.tmax,PAR.phistart,PAR.phistart,EQD,PAR,FLT);
						break;
					case 1: 		// random from r, theta
						FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,1);
						break;
					default: 	// grid from r, theta
						FLT.set(n,PAR.N,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,1,1);
						break;
				}
			}

			// use psi as control parameter and set to -1
			results(4,all) = -1;

			// Integrate
			for(i=1;i<=PAR.itt;i++)
			{
				chk = FLT.mapstep(PAR.MapDirection);
				if(chk<0) {ofs2 << "mapit: wall hit" << endl; break;}	// particle has left system

				if(fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) << endl;
				FLT.phi=PAR.MapDirection*i*dpinit*ilt + PAR.phistart;

				// Store results
#ifdef USE_SIESTA
				if(PAR.response_field == -2)
				{
					SIES.get_su(FLT.R, FLT.phi/rTOd, FLT.Z, s, u);
					results(1,(n-Nmin_slave)*PAR.itt + i) = u;
					results(4,(n-Nmin_slave)*PAR.itt + i) = s;
				}
				else
				{
					results(1,(n-Nmin_slave)*PAR.itt + i) = FLT.get_theta();
					results(4,(n-Nmin_slave)*PAR.itt + i) = FLT.psi;
				}
#else
				results(1,(n-Nmin_slave)*PAR.itt + i) = FLT.get_theta();
				results(4,(n-Nmin_slave)*PAR.itt + i) = FLT.psi;

#endif
				results(2,(n-Nmin_slave)*PAR.itt + i) = FLT.get_r();
				results(3,(n-Nmin_slave)*PAR.itt + i) = FLT.phi;
				results(5,(n-Nmin_slave)*PAR.itt + i) = FLT.R;
				results(6,(n-Nmin_slave)*PAR.itt + i) = FLT.Z;
			} // end for i
			ofs2 << "Trax: " << n << "\t" << "Steps: " << i-1 << endl;
		} // end for n
		
		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),6*N_slave*PAR.itt,MPI::DOUBLE,0,tag);

	}// end while
} // end Slaves
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double now2=zeit();
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
} //end main

//----------------------- End of Main -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

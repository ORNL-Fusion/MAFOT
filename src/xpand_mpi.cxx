// Calculate magnetic field outside of VMEC boundary

// Define
//--------
//#define BZ_DEBUG
#define USE_MPI
#define program_name "xpand_mpi"

// Include
//--------
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <omp.h>
#include <unistd.h>
#include <andi.hxx>
#include <vmec_class.hxx>
#include <xpand_class.hxx>
#include <adapt_gauss_kronrod_class.hxx>

// namespaces
//-----------
using namespace blitz;

// Prototypes
//-----------
Array<double,1> evGK(double Rin, double phiin, double Zin, AdaptiveGK& GKx, AdaptiveGK& GKy, BFIELDVC& bvc, double epsabs, double epsrel);
Array<double,1> B_fkt(double u, double v, BFIELDVC& bvc);
Array<double,1> B_u(double v, AdaptiveGK& GK, BFIELDVC& bvc, double epsabs, double epsrel);

class FUNCTIONx : public FtnBase
{
private:
	typedef Array<double,1> defn_t(double, double, BFIELDVC&);
	defn_t& function_;

	BFIELDVC& bvc_;
	double a;
	Array<double,1> out;

public:
	FUNCTIONx(defn_t& function, double A, BFIELDVC& bvc) : function_(function), bvc_(bvc)
	{
		a = A;
		out.resize(3);
	}
	virtual Array<double,1> operator() (double x) { out = function_(x, a, bvc_); return out;};
};

class FUNCTIONy : public FtnBase
{
private:
	typedef Array<double,1> defn_t(double, AdaptiveGK&, BFIELDVC&, double, double);
	defn_t& function_;

	AdaptiveGK& GK_;
	BFIELDVC& bvc_;
	double epsabs_;
	double epsrel_;
	Array<double,1> out;

public:
	FUNCTIONy(defn_t& function, AdaptiveGK& GK, BFIELDVC& bvc, double epsabs, double epsrel) : function_(function), GK_(GK), bvc_(bvc)
	{
		epsabs_ = epsabs;
		epsrel_ = epsrel;
		out.resize(3);
	}
	virtual Array<double,1> operator() (double x) { out = function_(x, GK_, bvc_, epsabs_, epsrel_); return out;};
};


// Switches
//----------
const bool ALLOW_INSIDE = true;	// compute field inside of VMEC LCFS

// Golbal Parameters
//------------------
ofstream logfile;

// Main Program
//--------------
int main(int argc, char *argv[])
{
// MPI initialize
MPI::Init(argc, argv);
int mpi_rank = MPI::COMM_WORLD.Get_rank();
int mpi_size = MPI::COMM_WORLD.Get_size();

// Variables
int i,j, N, idx, k;
double R, phi, Z, s, u, Pres;
Array<double,1> B(3), Bvac(3), Bvci(3);
double now=zeit();
Range all = Range::all();
int c;
ofstream ofs3;

int tag,sender;
int Nmin_slave,Nmax_slave;
Array<int,1> send_N_limits(Range(1,2));
MPI::Status status;

// adaptive integrator defaults
double epsabs = 1e-6;
double epsrel = 1e-4;
int limit = 1000;
int degree = 10;
bool use_GK = true;			// true: use adaptive Gauss-Kronrod integration;   false: use adaptive Simpson integration
bool VMEC_n0only = false;	// true: force VMEC to be axisymmetric, false: use as it is
LA_STRING mgrid_file = "None";
LA_STRING points_file = "points.dat";

bool VC_INSIDE = false;	// also run virtual casing inside of VMEC LCFS

// Command line input parsing
opterr = 0;
while ((c = getopt(argc, argv, "hP:a:r:i:m:SM:AI")) != -1)
switch (c)
{
case 'h':
	if(mpi_rank < 1)
	{
		cout << "usage: mpirun -n <cores> xpand_mpi [-h] [-P pts] [-a epsabs] [-r epsrel] [-i limit] [-m degree] [-S] [-M mgrid] [-A] wout [tag]" << endl << endl;
		cout << "Calculate magnetic field outside of VMEC boundary." << endl << endl;
		cout << "positional arguments:" << endl;
		cout << "  wout          VMEC wout-file name" << endl;
		cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
		cout << endl << "optional arguments:" << endl;
		cout << "  -h            show this help message and exit" << endl;
		cout << "  -P            use Points-File pts; default: 'points.dat' in the current working dir" << endl;
		cout << "  -a            set absolute tolerance; default 1e-6" << endl;
		cout << "  -r            set relative tolerance; default 1e-4" << endl;
		cout << "  -i            set maximum refinement limit, GK only; default 1000" << endl;
		cout << "  -m            set Gauss-Legendre degree (2*m+1), GK only; default m=10" << endl;
		cout << "  -S            use adaptive Simpson instead of Gauss-Kronrod" << endl;
		cout << "  -M            use Mgrid-File mgrid; default: 'mgrid_file' from wout" << endl;
		cout << "  -A            force axisymmetry -> use n=0 only in VMEC & toroidally average vacuum field" << endl;
		cout << "  -I            also run virtual casing inside of VMEC LCFS" << endl;
		cout << endl << "Examples:" << endl;
		cout << "  mpirun -n 4 xpand_mpi wout.nc" << endl;
		cout << "  mpirun -n 12 xpand_mpi -i 1500 -P ./here/my_points.dat wout.nc test" << endl;
		cout << "  mpirun -n 12 xpand_mpi -r 1e-8 -S -M /home/shared/mgrid_d3d.nc wout.nc test2" << endl;
	}
	MPI::Finalize();
	return 0;
case 'P':
	points_file = optarg;
	break;
case 'a':
	epsabs = atof(optarg);
	break;
case 'r':
	epsrel = atof(optarg);
	break;
case 'i':
	limit = atof(optarg);
	break;
case 'm':
	degree = atof(optarg);
	break;
case 'S':
	use_GK = false;
	break;
case 'M':
	mgrid_file = optarg;
	break;
case 'A':
	VMEC_n0only = true;
	break;
case 'I':
	VC_INSIDE = true;
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
LA_STRING wout_name;
LA_STRING praefix = "";
if(argc==optind+2) praefix = "_" + LA_STRING(argv[optind+1]);
if(argc>=optind+1) wout_name = LA_STRING(argv[optind]);
else {if(mpi_rank < 1) cout << "No Input files -> Abort!" << endl; EXIT;}

// check if enough Nodes have been selected
if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes and restart." << endl; EXIT;}

// Init classes
if(mpi_rank < 1)
{
	cout << "Initialize..." << endl;
	if(VC_INSIDE) cout << "Running VC inside s = 1 and correct VMEC field..." << endl;
	if(not ALLOW_INSIDE) cout << "Forcing inside s = 1 as outside..." << endl;
}
AdaptiveGK GKx(limit, degree, 3);
AdaptiveGK GKy(limit, degree, 3);
VMEC wout(wout_name, VMEC_n0only);
BFIELDVC bvc(wout, epsabs, epsrel, 14, mgrid_file);
INSIDE_VMEC inside(wout);

// read input
Array<double,2> points;
vector<LA_STRING> points_header;
if(mpi_rank < 1) cout << "Read points..." << endl;
int points_header_lines = readFileHeader(points_file, points_header);
readfile(points_file, 3, points);
N = points.rows();
double phi_old = points(1,2);

// Output
LA_STRING filenameout = "xpand" + praefix + ".dat";
if(mpi_rank < 1) outputtest(filenameout);

// Set starting parameters
int N_slave = 10;	// Number of points per package
int NoOfPackages = int(N/N_slave);
if(NoOfPackages < mpi_size)
{
	N_slave = int(N/mpi_size);	// this is already an integer division !!!
	if(N_slave < 1) N_slave = 1;
	NoOfPackages = int(N/N_slave);
}
int N_rest = N - NoOfPackages*N_slave;
const int N_transmit = 10;

// Needed for storage of data while calculating
Array<int,1> write_memory(Range(1,NoOfPackages)); //store which data is written to file (0: not calculated yet, 1: calculated, 2: written to file)
int write_max = 0;
int write_last = 0; //tag of last written data
write_memory = 0;

MPI::COMM_WORLD.Barrier();	// Syncronize all Nodes

// Master only (Node 0)
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank < 1)
{
	ofstream out(filenameout);
	out.precision(16);
	if(VMEC_n0only) out << "# Xpand results from axisymmetric VMEC file " << wout_name << endl;
	else out << "# Xpand results from VMEC file " << wout_name << endl;
	for(i=0;i<points_header_lines;i++) out << points_header[i] << endl;
	out << "# R[m]      \t phi[rad]   \t Z[m]       \t BR[T]      \t Bphi[T]    \t BZ[T]      \t Pressure[Pa] \t BRvac[T]   \t Bphivac[T] \t BZvac[T]   \t BRvci[T]   \t Bphivci[T] \t BZvci[T]" << endl;

	// log file
	if(use_GK) logfile.open("log_" + LA_STRING(program_name) + praefix + "_GK" + ".dat");
	ofs3.open("log_" + LA_STRING(program_name) + praefix + "_Master" + ".dat");
	ofs3.precision(16);
	if(use_GK) ofs3 << "Use Gauss-Kronrod integrator with limit = " << limit << " and degree = " << degree << endl;
	else ofs3 << "Use adaptive Simpson integrator" << endl;
	if(VMEC_n0only) ofs3 << "Using VMEC in axisymmetric mode" << endl;
	ofs3 << "Tolerances are: epsabs = " << epsabs << "\t epsrel = " << epsrel << endl;
	ofs3 << "Calculate B-field for " << N << " points" << endl;
	ofs3 << "No. of Packages = " << NoOfPackages << " Points per Package = " << N_slave << endl << endl;

	// Result array:	 Column Number,  Values
	Array<double,3> results_all(Range(1,NoOfPackages),Range(1,N_transmit),Range(1,N_slave));
	Array<double,2> recieve(Range(1,N_transmit),Range(1,N_slave));
	Array<double,2> slice;
	tag = 1;	// first Package

	// number of field lines array
	Array<int,1> N_values(Range(1,N));
	for(i=1;i<=N;i++) N_values(i) = i;

	int sent_packages = 0;
	int recieve_packages = 0;
	int count = N;

	#pragma omp parallel shared(results_all,N_values,sent_packages,recieve_packages,count) private(i,tag) num_threads(2)
	{
		#pragma omp sections nowait
		{
		#pragma omp section	//-------- Master Thread: controlles comunication ----------------------------------------------------------------------------------------------------------------------
		{
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
					send_N_limits(2) = N_values(tag*N_slave);	// Nmax_slave

					MPI::COMM_WORLD.Send(send_N_limits.dataFirst(),2,MPI::INTEGER,i,tag);
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
				MPI::COMM_WORLD.Recv(recieve.dataFirst(),N_transmit*N_slave,MPI::DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,status);
				sender = status.Get_source();
				tag = status.Get_tag();
				ofs3 << "Recieve from Node: " << sender << " Package: " << tag << endl;

				#pragma omp critical
				{
					count -= N_slave;
					recieve_packages += 1;
					write_memory(tag) = 1;
					ofs3 << "------------------------------------ Progress: " << recieve_packages << " of " << NoOfPackages << " completed" <<  endl;
				}
				slice.reference(results_all(tag,all,all));
				slice = recieve;

				// show progress
				cout << "\rDone: " << int(10000*double(N-count)/double(N))/100.0 << "%   " << flush;

				// Send new Package
				if(sent_packages < NoOfPackages)
				{
					#pragma omp critical
					{
						sent_packages += 1;
						tag = sent_packages;
					}

					send_N_limits(1) = N_values((tag-1)*N_slave+1);	// Nmin_slave
					send_N_limits(2) = N_values(tag*N_slave);	// Nmax_slave

					MPI::COMM_WORLD.Send(send_N_limits.dataFirst(),2,MPI::INTEGER,sender,tag);
					ofs3 << "Send again to Node: " << sender << " Package No.: " << tag << endl;
				}
				else	// No Packages left -> Send termination signal: tag = 0
				{
					MPI::COMM_WORLD.Send(send_N_limits.dataFirst(),0,MPI::INTEGER,sender,0);
					workingNodes -= 1;
					ofs3 << "Send termination signal to Node: " << sender << endl;
				}

				// Set write_max
				while((write_memory(write_max+1) == 1) && (write_max < NoOfPackages)) write_max++;

				// Write Output to file
				for(i=write_last+1; i<=write_max; i++)
				{
					for(j=1;j<=N_slave;j++)
					{
						idx = N_values((i-1)*N_slave+j);
						out << points(idx,1) << "\t" << points(idx,2) << "\t" << points(idx,3);
						if(VC_INSIDE) {for(k=1;k<=N_transmit;k++) out << "\t" << results_all(i,k,j);}
						else {for(k=1;k<=7;k++) out << "\t" << results_all(i,k,j);}
						out << endl;
					}
					write_memory(i) = 2;
				}
				write_last = write_max;

			} // end while(recieve_packages <= NoOfPackages)
		} // end omp section: Master thread

		#pragma omp section	//------- Slave Thread on Master Node: does same calculations as Salve Nodes ----------------------------------------------------------------------------------
		{
			cout << "Calculate B-field for " << N << " points..." << endl;

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
				ofs3 << "Node: " << mpi_rank << " works on Package: " << tag << endl;

				// Prepare inside
				inside.init(phi_old);
				for(i=1;i<=N_slave;i++)
				{
					idx = N_values((tag-1)*N_slave+i);
					R = points(idx,1);
					phi = points(idx,2);
					Z = points(idx,3);
					if(phi != phi_old)
					{
						inside.init(phi);
						phi_old = phi;
					}

					Bvac = bvc.get_vacuumB(R, phi, Z, VMEC_n0only);
					Bvci = 0;

					if(inside.check(R,Z) && ALLOW_INSIDE)
					{
						wout.get_su(R, phi, Z, s, u);
						Pres = wout.presf.ev(s);
						wout.get_B2D(s, u, phi, B(0), B(1), B(2));
						if(VC_INSIDE)
						{
							if(use_GK) Bvci = evGK(R, phi, Z, GKx, GKy, bvc, epsabs, epsrel);
							else Bvci = bvc.ev(R, phi, Z);
							B += Bvac + Bvci;	// Bvci ~ -Bvac
						}
					}
					else
					{
						if(use_GK) B = evGK(R, phi, Z, GKx, GKy, bvc, epsabs, epsrel);
						else B = bvc.ev(R, phi, Z);
						B += Bvac;
						Pres = wout.presf.y(wout.ns); // just the value of presf at LCFS; same as: wout.presf.ev(1.0)
					}

					// Store results
					results_all(tag,1,i) = B(0);
					results_all(tag,2,i) = B(1);
					results_all(tag,3,i) = B(2);
					results_all(tag,4,i) = Pres;
					results_all(tag,5,i) = Bvac(0);
					results_all(tag,6,i) = Bvac(1);
					results_all(tag,7,i) = Bvac(2);
					results_all(tag,8,i) = Bvci(0);
					results_all(tag,9,i) = Bvci(1);
					results_all(tag,10,i) = Bvci(2);
				} // end for

				#pragma omp critical
				{
					count -= N_slave;
					recieve_packages += 1;
					write_memory(tag) = 1;
					ofs3 << "------------------------------------ Progress: " << recieve_packages << " of " << NoOfPackages << " completed" <<  endl;
				}

			}// end while
		}// end omp section: Slave Thread
		} // end omp sections
	} // end omp parallel

	// Write remaining Output to file
	for(i=write_last+1;i<=NoOfPackages;i++)
	{
		for(j=1;j<=N_slave;j++)
		{
			idx = N_values((i-1)*N_slave+j);
			out << points(idx,1) << "\t" << points(idx,2) << "\t" << points(idx,3);
			if(VC_INSIDE) {for(k=1;k<=N_transmit;k++) out << "\t" << results_all(i,k,j);}
			else {for(k=1;k<=7;k++) out << "\t" << results_all(i,k,j);}
			out << endl;
		}
	}

	// calculate the N_rest remaining points
	for(i=N-N_rest+1;i<=N;i++)
	{
		R = points(i,1);
		phi = points(i,2);
		Z = points(i,3);
		if(phi != phi_old)
		{
			inside.init(phi);
			phi_old = phi;
		}

		Bvac = bvc.get_vacuumB(R, phi, Z, VMEC_n0only);
		Bvci = 0;

		if(inside.check(R,Z) && ALLOW_INSIDE)
		{
			wout.get_su(R, phi, Z, s, u);
			Pres = wout.presf.ev(s);
			wout.get_B2D(s, u, phi, B(0), B(1), B(2));
			if(VC_INSIDE)
			{
				if(use_GK) Bvci = evGK(R, phi, Z, GKx, GKy, bvc, epsabs, epsrel);
				else Bvci = bvc.ev(R, phi, Z);
				B += Bvac + Bvci;
			}
		}
		else
		{
			if(use_GK) B = evGK(R, phi, Z, GKx, GKy, bvc, epsabs, epsrel);
			else B = bvc.ev(R, phi, Z);
			B += Bvac;
			Pres = wout.presf.y(wout.ns);
		}
		// show progress
		count -= 1;
		cout << "\rDone: " << int(10000*double(N-count)/double(N))/100.0 << "%   " << flush;
		// Output
		if(VC_INSIDE) out << R << "\t" << phi << "\t" << Z << "\t" << B(0) << "\t" << B(1) << "\t" << B(2) << "\t" << Pres << "\t" << Bvac(0) << "\t" << Bvac(1) << "\t" << Bvac(2) << "\t" << Bvci(0) << "\t" << Bvci(1) << "\t" << Bvci(2) << endl;
		else out << R << "\t" << phi << "\t" << Z << "\t" << B(0) << "\t" << B(1) << "\t" << B(2) << "\t" << Pres << "\t" << Bvac(0) << "\t" << Bvac(1) << "\t" << Bvac(2) << endl;
	}
} // end Master
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	// Result array for Slave
	Array<double,2> results(Range(1,N_transmit),Range(1,N_slave));

	// Start working...
	while(1)
	{
		// Recieve initial conditions to calculate
		MPI::COMM_WORLD.Recv(send_N_limits.dataFirst(),2,MPI::INTEGER,0,MPI_ANY_TAG,status);

		tag = status.Get_tag();
		if(tag == 0) break;	// Still work to do?  ->  tag = 0: NO, stop working

		Nmin_slave = send_N_limits(1);
		Nmax_slave = send_N_limits(2);

		// Prepare inside
		inside.init(phi_old);
		for(i=1;i<=N_slave;i++)
		{
			idx = Nmin_slave + i - 1;
			R = points(idx,1);
			phi = points(idx,2);
			Z = points(idx,3);
			if(phi != phi_old)
			{
				inside.init(phi);
				phi_old = phi;
			}

			Bvac = bvc.get_vacuumB(R, phi, Z, VMEC_n0only);
			Bvci = 0;

			if(inside.check(R,Z) && ALLOW_INSIDE)
			{
				wout.get_su(R, phi, Z, s, u);
				Pres = wout.presf.ev(s);
				wout.get_B2D(s, u, phi, B(0), B(1), B(2));
				if(VC_INSIDE)
				{
					if(use_GK) Bvci = evGK(R, phi, Z, GKx, GKy, bvc, epsabs, epsrel);
					else Bvci = bvc.ev(R, phi, Z);
					B += Bvac + Bvci;
				}
			}
			else
			{
				if(use_GK) B = evGK(R, phi, Z, GKx, GKy, bvc, epsabs, epsrel);
				else B = bvc.ev(R, phi, Z);
				B += Bvac;
				Pres = wout.presf.y(wout.ns);
			}

			// Store results
			results(1,i) = B(0);
			results(2,i) = B(1);
			results(3,i) = B(2);
			results(4,i) = Pres;
			results(5,i) = Bvac(0);
			results(6,i) = Bvac(1);
			results(7,i) = Bvac(2);
			results(8,i) = Bvci(0);
			results(9,i) = Bvci(1);
			results(10,i) = Bvci(2);
		} // end for

		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),N_transmit*N_slave,MPI::DOUBLE,0,tag);

	}// end while
} // end Slaves
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double now2 = zeit();
if(mpi_rank < 1)
{
	cout << endl << "Program terminates normally, Time: " << now2-now  << " s" << endl;
	ofs3 << "Program terminates normally, Time: " << now2-now  << " s" << endl;
}

// MPI finalize
MPI::Finalize();

return 0;
}  //end of main


//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- ev -----------------------------------------------------------------------------------------------------------------
Array<double,1> evGK(double Rin, double phiin, double Zin, AdaptiveGK& GKx, AdaptiveGK& GKy, BFIELDVC& bvc, double epsabs, double epsrel)
{
Array<double,1> integ(3);
double error;
bvc.R = Rin; bvc.phi = phiin; bvc.Z = Zin; 		// load into member variables
FUNCTIONy fBR(B_u, GKx, bvc, epsabs, epsrel);
GKy.qag(fBR, 0, pi2, epsabs, epsrel, integ, error);
integ /= -4*pi;	// mu0/4pi and mu0*H = B; because K = n x H originally; integ is the B-field in carthesian coordinates!!!

// transform back to cylindrical
double Bx, By;
Bx = integ(0); By = integ(1);
integ(0) =  Bx*cos(phiin) + By*sin(phiin);
integ(1) = -Bx*sin(phiin) + By*cos(phiin);
return integ;
}


// --- B_fkt -----------------------------------------
Array<double,1> B_fkt(double u, double v, BFIELDVC& bvc)
{
Array<double,1> out(3);
out = bvc.integs(u, v);
return out;
}

// --- BR_u -----------------------------------------
Array<double,1> B_u(double v, AdaptiveGK& GK, BFIELDVC& bvc, double epsabs, double epsrel)
{
double error;
Array<double,1> result(3);
FUNCTIONx f(B_fkt, v, bvc);
GK.qag(f, 0, pi2, epsabs, epsrel, result, error);	// B(v = const)
return result;
}


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------




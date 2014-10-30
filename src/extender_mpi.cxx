// Calculate magnetic field outside of VMEC boundary

// Define
//--------
//#define BZ_DEBUG
#define USE_MPI
#define program_name "extender_mpi"

// Include
//--------
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <andi.hxx>
#include <vmec_class.hxx>
#include <extender_class.hxx>
#include <omp.h>
#include <unistd.h>

// namespaces
//-----------
using namespace blitz;

// Prototypes
//-----------

// Switches
//----------

// Golbal Parameters
//------------------

// Main Program
//--------------
int main(int argc, char *argv[])
{
// MPI initialize
MPI::Init(argc, argv);
int mpi_rank = MPI::COMM_WORLD.Get_rank();
int mpi_size = MPI::COMM_WORLD.Get_size();
if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes and restart." << endl; EXIT;}

// Variables
int i,j, N, idx;
double R, phi, Z, s, u;
Array<double,1> B(3), Bvac(3);
double now=zeit();
Range all = Range::all();
int c;

int tag,sender;
int Nmin_slave,Nmax_slave;
Array<int,1> send_N_limits(Range(1,2));
MPI::Status status;

// adaptive Simpson accuracy defaults
double epsabs = 1e-6;
double epsrel = 1e-4;

// Command line input parsing
opterr = 0;
while ((c = getopt(argc, argv, "ha:r:")) != -1)
switch (c)
{
case 'h':
	if(mpi_rank < 1)
	{
		cout << "usage: mpirun -n <cores> extender_mpi [-h] [-a epsabs] [-r epsrel] wout [tag]" << endl << endl;
		cout << "Calculate magnetic field outside of VMEC boundary." << endl << endl;
		cout << "positional arguments:" << endl;
		cout << "  wout          VMEC wout-file name" << endl;
		cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
		cout << endl << "optional arguments:" << endl;
		cout << "  -h            show this help message and exit" << endl;
		cout << "  -a            set absolute tolerance for Simpson; default 1e-6" << endl;
		cout << "  -r            set relative tolerance for Simpson; default 1e-4" << endl;
		cout << endl << "Examples:" << endl;
		cout << "  mpirun -n 4 extender_mpi wout.nc" << endl;
		cout << "  mpirun -n 12 extender_mpi -a 1e-8 wout.nc test" << endl;
		cout << "  mpirun -n 12 extender_mpi -a 1e-8 -r 1e-6 wout.nc test2" << endl;
	}
	MPI::Finalize();
	return 0;
case 'a':
	epsabs = atof(optarg);
	break;
case 'r':
	epsrel = atof(optarg);
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
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}

// log file
//ofs2.open("log_" + LA_STRING(program_name) + praefix + ".dat");
ofstream ofs3;

// Init classes
if(mpi_rank < 1) cout << "Initialize..." << endl;
VMEC wout(wout_name);
POTENTIAL Pot(wout, epsabs, epsrel, 14);
INSIDE_VMEC inside(wout);

// read input
Array<double,2> points;
if(mpi_rank < 1) cout << "Read points..." << endl;
readfile("points.dat", 3, points);
N = points.rows();
double phi_old = points(0,2);

// Output
LA_STRING filenameout = "ext" + praefix + ".dat";
if(mpi_rank < 1) outputtest(filenameout);

// Set starting parameters
int N_slave = 10;	// Number of points per package
int NoOfPackages = int(N/N_slave);
int N_rest = N - NoOfPackages*N_slave;

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
	out << "# Extender results from VMEC file " << wout_name << endl;
	out << "# " << N << " points" << endl;
	out << "# R[m]      \t phi[rad]   \t Z[m]       \t BR[T]      \t Bphi[T]    \t BZ[T]      \t BRvac[T]   \t Bphivac[T] \t BZvac[T]" << endl;

	// log file
	ofs3.open("log_" + LA_STRING(program_name) + praefix + "_Master" + ".dat");
	ofs3.precision(16);
	ofs3 << "Calculate B-field for " << N << " points" << endl;
	ofs3 << "No. of Packages = " << NoOfPackages << " Points per Package = " << N_slave << endl << endl;

	// Result array:	 Column Number,  Values
	Array<double,3> results_all(Range(1,NoOfPackages),Range(1,6),Range(1,N_slave));
	Array<double,2> recieve(Range(1,6),Range(1,N_slave));
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
				MPI::COMM_WORLD.Recv(recieve.dataFirst(),6*N_slave,MPI::DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,status);
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
						out << points(idx,1) << "\t" << points(idx,2) << "\t" << points(idx,3) << "\t" << results_all(i,1,j) << "\t" << results_all(i,2,j) << "\t" << results_all(i,3,j) << "\t" << results_all(i,4,j) << "\t" << results_all(i,5,j) << "\t" << results_all(i,6,j) << endl;
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

					if(inside.check(R,Z))
					{
						wout.get_su(R, phi, Z, s, u);
						wout.get_B2D(s, u, phi, B(0), B(1), B(2));
						Bvac = 0;
					}
					else
					{
						B = Pot.grad(R, phi, Z);
						Bvac = Pot.get_vacuumB(R, phi, Z);
						B += Bvac;
					}

					// Store results
					results_all(tag,1,i) = B(0);
					results_all(tag,2,i) = B(1);
					results_all(tag,3,i) = B(2);
					results_all(tag,4,i) = Bvac(0);
					results_all(tag,5,i) = Bvac(1);
					results_all(tag,6,i) = Bvac(2);
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
			out << points(idx,1) << "\t" << points(idx,2) << "\t" << points(idx,3) << "\t" << results_all(i,1,j) << "\t" << results_all(i,2,j) << "\t" << results_all(i,3,j) << "\t" << results_all(i,4,j) << "\t" << results_all(i,5,j) << "\t" << results_all(i,6,j) << endl;
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

		if(inside.check(R,Z))
		{
			wout.get_su(R, phi, Z, s, u);
			wout.get_B2D(s, u, phi, B(0), B(1), B(2));
			Bvac = 0;
		}
		else
		{
			B = Pot.grad(R, phi, Z);
			Bvac = Pot.get_vacuumB(R, phi, Z);
			B += Bvac;
		}
		// show progress
		count -= 1;
		cout << "\rDone: " << int(10000*double(N-count)/double(N))/100.0 << "%   " << flush;
		// Output
		out << R << "\t" << phi << "\t" << Z << "\t" << B(0) << "\t" << B(1) << "\t" << B(2) << "\t" << Bvac(0) << "\t" << Bvac(1) << "\t" << Bvac(2) << endl;
	}
} // end Master
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	// Result array for Slave
	Array<double,2> results(Range(1,6),Range(1,N_slave));

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

			if(inside.check(R,Z))
			{
				wout.get_su(R, phi, Z, s, u);
				wout.get_B2D(s, u, phi, B(0), B(1), B(2));
				Bvac = 0;
			}
			else
			{
				B = Pot.grad(R, phi, Z);
				Bvac = Pot.get_vacuumB(R, phi, Z);
				B += Bvac;
			}

			// Store results
			results(1,i) = B(0);
			results(2,i) = B(1);
			results(3,i) = B(2);
			results(4,i) = Bvac(0);
			results(5,i) = Bvac(1);
			results(6,i) = Bvac(2);
		} // end for

		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),6*N_slave,MPI::DOUBLE,0,tag);

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



//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------




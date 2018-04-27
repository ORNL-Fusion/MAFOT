// Program calculates connection length and penetration depth footprints on the target plates
// include drifts and time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						20.06.11

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	2d footprint data for colored contour plot
//			log-file


// uses Parallel computation by open-mpi and openmp
// IMPORTANT: Node 0 starts two threads using openmp. The Master-thread contols communications with MPI-Slaves. The Slave-thread on Node 0 does the same calculations as the MPI-Slaves
// Launch Syntax:				mpirun -np [NoOfNodes] -hostfile [FileName] dtfoot_mpi [Parameterfile] [Praefix(optional)] &
// alternative Launch Syntax:	mpirun -np [NoOfNodes] -host [HostNames separated by Commata] dtfoot_mpi [Parameterfile] [Praefix(optional)] &

// Define
//--------
#define USE_MPI
#if defined(ITER)
	#define program_name "iterfoot_mpi"
#elif defined(NSTX)
	#define program_name "nstxfoot_mpi"
#elif defined(MAST)
	#define program_name "mastfoot_mpi"
#else
	#define program_name "dtfoot_mpi"
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

// Prototypes  
//-----------

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
int i,j;
int chk;
double t,ntor,length,psimin;
Range all = Range::all();
EFIT EQD;

int tag,sender;
double tmin_slave,tmax_slave,dt;
Array<double,1> send_t_limits(Range(1,2));
MPI::Status status;

// Use system time as seed(=idum) for random numbers
double now = zeit();

// defaults
LA_STRING woutfile = "wout.nc";
LA_STRING xpandfile = "xpand.dat";
LA_STRING siestafile = "siesta.dat";
LA_STRING islandfile = "fakeIslands.in";

// Command line input parsing
int c;
opterr = 0;
while ((c = getopt(argc, argv, "hX:V:S:I:")) != -1)
switch (c)
{
case 'h':
	if(mpi_rank < 1)
	{
		cout << "usage: mpirun -n <cores> dtfoot_mpi [-h] [-I island] [-S siesta] [-V wout] [-X xpand] file [tag]" << endl << endl;
		cout << "Calculate field line connection length and penetration depth on the vessel wall." << endl << endl;
		cout << "positional arguments:" << endl;
		cout << "  file          Contol file (starts with '_')" << endl;
		cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
		cout << endl << "optional arguments:" << endl;
		cout << "  -h            show this help message and exit" << endl;
		cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
		cout << "  -S            filename for SIESTA; default, see below" << endl;
		cout << "  -V            filename for VMEC; default, see below" << endl;
		cout << "  -X            filename for XPAND; default, see below" << endl;
		cout << endl << "Examples:" << endl;
		cout << "  mpirun -n 4 dtfoot_mpi _inner.dat blabla" << endl;
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
else	// No Input: Abort
{
	if(mpi_rank < 1) cout << "No Input files -> Abort!" << endl;
	EXIT;
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";
if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes and restart." << endl; EXIT;}

// Read Parameterfile
if(mpi_rank < 1) cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,11,mpi_rank);

// Set target type for output-filename
LA_STRING type;
switch(PAR.which_target_plate)
{
case 0:
	type = "_wall";
	break;
case 1:
	#if defined(NSTX)
		type = "_inup";
	#else
		type = "_in";
	#endif
	break;
case 2:
	#if defined(NSTX)
		type = "_outup";
	#else
		type = "_out";
	#endif
	break;
case 3:
	#if defined(NSTX)
		type = "_indwn";
	#else
		type = "_shelf";
	#endif
	break;
case 4:
	#if defined(NSTX)
		type = "_outdwn";
	#else
		type = "_sas";
	#endif
	break;
case 5:
	type = "_bsas";
	break;
default:
	type = "";
	break;
}

// log file
ofs2.open("log_" + LA_STRING(program_name) + type + praefix + "_Node" + LA_STRING(mpi_rank) + ".dat");
ofs2.precision(16);
ofstream ofs3;

// Output
LA_STRING filenameout = "foot" + type + praefix + ".dat";
if(mpi_rank < 1) outputtest(filenameout);
MPI::COMM_WORLD.Barrier();	// All Nodes wait for Master

// Read EFIT-data
double Raxis = 0, Zaxis = 0;
#ifdef USE_XFIELD
if(PAR.response_field == -3)
{
	VMEC vmec;
	double Raxisv,Zaxisv,v;
	if(mpi_rank < 1) cout << "Read VMEC file" << endl;
	ofs2 << "Read VMEC file" << endl;
	vmec.read(woutfile);
	for(i=0;i<40;i++)
	{
		v = i*pi2/40.0;
		vmec.get_axis(v,Raxisv,Zaxisv);
		Raxis += Raxisv;
		Zaxis += Zaxisv;
	}
	Raxis /= 40.0;
	Zaxis /= 40.0;
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

EQD.ReadData(EQD.Shot,EQD.Time,Raxis,Zaxis);
if(mpi_rank < 1) cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
int N_variables = 7;
int Nt_slave = 1;
int N = PAR.Nt*PAR.Nphi;
int N_slave = PAR.Nphi*Nt_slave;
int NoOfPackages = int(PAR.Nt/Nt_slave);

// Needed for storage of data while calculating
Array<int,1> write_memory(Range(1,NoOfPackages)); //store which data is written to file (0: not calculated yet, 1: calculated, 2: written to file)
int write_max = 0;
int write_last = 0; //tag of last written data
write_memory = 0;

if(PAR.Nt<=1) {dt=0;}
else dt=(PAR.tmax-PAR.tmin)/(PAR.Nt-1);
if(dt == 0) PAR.tmax = PAR.tmin;

// Use box as boundary
simpleBndy = 1;

// Prepare particles
PARTICLE FLT(EQD,PAR,mpi_rank);

MPI::COMM_WORLD.Barrier();	// Syncronize all Nodes

// Master only (Node 0)
//--------------------------------------------------------------------------------------------------
if(mpi_rank < 1)
{
	// log file
	ofs3.open("log_" + LA_STRING(program_name) + type + praefix + "_Master" + ".dat");
	ofs3.precision(16);

	ofs3 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
	ofs3 << "No. of Packages = " << NoOfPackages << " Points per Package = " << N_slave << endl << endl;

	// additional parameters for IO
	PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
	PAR.pv[1].name = "t-grid";			PAR.pv[1].wert = PAR.Nt;
	PAR.pv[2].name = "phi-grid";		PAR.pv[2].wert = PAR.Nphi;
	PAR.pv[3].name = "tmin";			PAR.pv[3].wert = PAR.tmin;
	PAR.pv[4].name = "tmax";			PAR.pv[4].wert = PAR.tmax;
	PAR.pv[5].name = "phimin";			PAR.pv[5].wert = PAR.phimin;
	PAR.pv[6].name = "phimax";			PAR.pv[6].wert = PAR.phimax;
	PAR.pv[7].name = "MapDirection";	PAR.pv[7].wert = PAR.MapDirection;
	PAR.pv[8].name = "Ekin";			PAR.pv[8].wert = PAR.Ekin;
	PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;
	PAR.pv[10].name = "total length of wall";	PAR.pv[10].wert = EQD.Swall_max;

	// Output
	ofstream out(filenameout);
	out.precision(16);
	vector<LA_STRING> var(7);
	var[0] = "phi[rad]";  var[1] = "length t";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";  var[5] = "R [m]";  var[6] = "Z [m]";
	PAR.writeiodata(out,bndy,var);

	// Result array:					Package ID,  Column Number,  Values
	Array<double,3> results_all(Range(1,NoOfPackages),Range(1,N_variables),Range(1,N_slave));
	Array<double,2> recieve(Range(1,N_variables),Range(1,N_slave));
	Array<double,2> slice;
	tag = 1;	// first Package

	// t array
	Array<double,1> t_values(Range(1,PAR.Nt));
	for(i=1;i<=PAR.Nt;i++) t_values(i) = PAR.tmin + (i-1)*dt;

	int sent_packages = 0;
	int recieve_packages = 0;

	#pragma omp parallel shared(results_all,t_values,sent_packages,recieve_packages,write_memory) private(i,tag) num_threads(2)
	{
		#pragma omp sections nowait
		{
		#pragma omp section	//-------- Master Thread: controlles comunication ----------------------------------------------------------------------------------------------------------------------
		{
			//#pragma omp barrier	// Syncronize with Slave Thread
			MPI::COMM_WORLD.Barrier();	// Master waits for Slaves

			//ofs2 << "Target (0=vertical, 1=45�, 2=horizontal, 3=shelf): " << which_target_plate << endl;
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
					send_t_limits(1) = t_values((tag-1)*Nt_slave+1);	// tmin_slave
					send_t_limits(2) = t_values(tag*Nt_slave);	// tmax_slave

					MPI::COMM_WORLD.Send(send_t_limits.dataFirst(),2,MPI::DOUBLE,i,tag);
					workingNodes += 1;

					ofs3 << "Send Package No.: " << tag << endl;
				}
				else	// more Nodes than Packages -> Send termination signal: tag = 0
				{
					MPI::COMM_WORLD.Send(send_t_limits.dataFirst(),0,MPI::DOUBLE,i,0);
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

					send_t_limits(1) = t_values((tag-1)*Nt_slave+1);	// tmin_slave
					send_t_limits(2) = t_values(tag*Nt_slave);	// tmax_slave

					MPI::COMM_WORLD.Send(send_t_limits.dataFirst(),2,MPI::DOUBLE,sender,tag);
					ofs3 << "Send again to Node: " << sender << " Package No.: " << tag << endl;
				}
				else	// No Packages left -> Send termination signal: tag = 0
				{
					MPI::COMM_WORLD.Send(send_t_limits.dataFirst(),0,MPI::DOUBLE,sender,0);
					workingNodes -= 1;
					ofs3 << "Send termination signal to Node: " << sender << endl;
				}

				// Set write_max
				while((write_memory(write_max+1) == 1) && (write_max < NoOfPackages)) write_max++;

				// Write Output to file
				for(i=write_last+1; i<=write_max; i++)
				{
					//for(j=1;j<=N_slave;j++)	out << results_all(i,1,j) << "\t" << results_all(i,2,j) << "\t" << results_all(i,3,j) << "\t" << results_all(i,4,j) << "\t" << results_all(i,5,j) << endl;
					for(j=1;j<=N_slave;j++)
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

			cout << "Target: " << PAR.which_target_plate << endl;
			cout << "Start Tracer for " << N << " points ... " << endl;
			ofs2 << "Target: " << PAR.which_target_plate << endl;

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
				ofs3 << "Node: " << mpi_rank << " works on Package: " << tag << endl;

				tmin_slave = t_values((tag-1)*Nt_slave+1);
				tmax_slave = t_values(tag*Nt_slave);

				ofs2 << "Start Tracer for " << N_slave << " points ... " << endl;
				for(i=1;i<=N_slave;i++)
				{
					// Set and store initial condition
					t = start_on_target(i,Nt_slave,PAR.Nphi,tmin_slave,tmax_slave,PAR.phimin,PAR.phimax,EQD,PAR,FLT);
					results_all(tag,1,i) = FLT.phi/rTOd;	// phi in rad
					results_all(tag,2,i) = t;
					results_all(tag,6,i) = FLT.R;
					results_all(tag,7,i) = FLT.Z;

					chk = FLT.connect(ntor,length,psimin,PAR.itt,PAR.MapDirection);

					//Store results
					results_all(tag,3,i) = ntor;
					results_all(tag,4,i) = length/1000.0;
					results_all(tag,5,i) = psimin;

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
	for(i=write_last+1;i<=NoOfPackages;i++)
	{
		//for(j=1;j<=N_slave;j++)	out << results_all(i,1,j) << "\t" << results_all(i,2,j) << "\t" << results_all(i,3,j) << "\t" << results_all(i,4,j) << "\t" << results_all(i,5,j) << endl;
		for(j=1;j<=N_slave;j++)
		{
			out << results_all(i,1,j);
			for(int k=2;k<=N_variables;k++) out << "\t" << results_all(i,k,j);
			out << endl;
		}
	}
} // end Master
//---------------------------------------------------------------------------------------------------


// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	// Prepare Perturbation
	prepare_common_perturbations(EQD,PAR,mpi_rank,siestafile,xpandfile,islandfile);
	prep_perturbation(EQD,PAR,mpi_rank);

	MPI::COMM_WORLD.Barrier();	// Syncronize with Master

	ofs2 << "Target: " << PAR.which_target_plate << endl;

	// Result array for Slave
	Array<double,2> results(Range(1,N_variables),Range(1,N_slave));

	// Start working...
	while(1)	
	{
		// Recieve initial conditions to calculate
		MPI::COMM_WORLD.Recv(send_t_limits.dataFirst(),2,MPI::DOUBLE,0,MPI_ANY_TAG,status);

		tag = status.Get_tag();
		ofs2 << endl;
		ofs2 << "Node: " << mpi_rank << " works on Package: " << tag << endl;
		if(tag == 0) break;	// Still work to do?  ->  tag = 0: NO, stop working

		tmin_slave = send_t_limits(1);
		tmax_slave = send_t_limits(2);

		ofs2 << "Start Tracer for " << N_slave << " points ... " << endl;
		for(i=1;i<=N_slave;i++)
		{
			// Set and store initial condition
			t = start_on_target(i,Nt_slave,PAR.Nphi,tmin_slave,tmax_slave,PAR.phimin,PAR.phimax,EQD,PAR,FLT);
			results(1,i) = FLT.phi/rTOd;	// phi in rad
			results(2,i) = t;
			results(6,i) = FLT.R;
			results(7,i) = FLT.Z;

			chk = FLT.connect(ntor,length,psimin,PAR.itt,PAR.MapDirection);

			//Store results
			results(3,i) = ntor;
			results(4,i) = length/1000.0;
			results(5,i) = psimin;

			if(i%100==0) ofs2 << "Trax: " << i << endl;
		}

		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),N_variables*N_slave,MPI::DOUBLE,0,tag);

	}// end while
} // end Slaves
//---------------------------------------------------------------------------------------------------

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



//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

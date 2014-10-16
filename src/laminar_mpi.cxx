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
#if defined(ITER)
	#if defined(m3dc1)
		#include <iter_m3dc1.hxx>
	#else
		#include <iter.hxx>
	#endif
#elif defined(NSTX)
	#if defined(m3dc1)
		#include <nstx_m3dc1.hxx>
	#else
		#include <nstx.hxx>
	#endif
#elif defined(MAST)
	#if defined(m3dc1)
		#include <mast_m3dc1.hxx>
	#else
		#include <mast.hxx>
	#endif
#else
	#if defined(m3dc1)
		#include <d3d_m3dc1.hxx>
	#else
		#include <d3d.hxx>
	#endif
#endif
#include <omp.h>

// Prototypes  
//-----------

// Switches
//----------
const int spare_interior = 1;	// 0: all points are calculated		1: inside psi = psi_interior_limit results are set to fixed values (code runs faster)

// Golbal Parameters 
//------------------
const double psi_interior_limit = 0.85;		// psi limit for spare_interior == 1

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// MPI initialize
MPI::Init(argc, argv);
int mpi_rank = MPI::COMM_WORLD.Get_rank();
int mpi_size = MPI::COMM_WORLD.Get_size();
if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes and restart." << endl; EXIT;}

// Variables
EFIT EQD;
int i,j;
int chk,skip_connect;
double ntor,length,psimin,psimax,psiav;
Range all = Range::all();

int tag,sender;
double Zmin_slave,Zmax_slave,dz;
Array<double,1> send_Z_limits(Range(1,2));
MPI::Status status;

// Use system time as seed(=idum) for random numbers
double now = zeit();
//long idum = long(now);

// Input file names
LA_STRING basename;
LA_STRING praefix = "";
if(argc==3) praefix = "_" + LA_STRING(argv[2]);
if(argc>=2) basename = LA_STRING(argv[1]);
else	// No Input: Abort
{
	if(mpi_rank < 1) cout << "No Input files -> Abort!" << endl;
	EXIT;
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

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

ofs2 << "ok" << endl;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
if(mpi_rank < 1) cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
int NZ_slave = 1;
int N = PAR.NR*PAR.NZ;
int N_slave = PAR.NR*NZ_slave;
int NoOfPackages = int(PAR.NZ/NZ_slave);

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

	// Output
	ofstream out(filenameout);
	out.precision(16);
	vector<LA_STRING> var(7);
	var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";  var[5] = "psimax";  var[6] = "psiav";
	if(PAR.create_flag == 3) {var[0] = "theta";  var[1] = "psi";}
	PAR.writeiodata(out,bndy,var);

	// Result array:					Package ID,  Column Number,  Values
	Array<double,3> results_all(Range(1,NoOfPackages),Range(1,7),Range(1,N_slave));
	Array<double,2> recieve(Range(1,7),Range(1,N_slave));
	Array<double,2> slice;
	tag = 1;	// first Package

	// Z array
	Array<double,1> Z_values(Range(1,PAR.NZ));
	for(i=1;i<=PAR.NZ;i++) Z_values(i) = PAR.Zmin + (i-1)*dz;

	int sent_packages = 0;
	int recieve_packages = 0;

	#pragma omp parallel shared(results_all,Z_values,sent_packages,recieve_packages,write_memory) private(i,tag) num_threads(2)
	{
		#pragma omp sections nowait
		{
		#pragma omp section	//-------- Master Thread: controlles comunication ----------------------------------------------------------------------------------------------------------------------
		{
			//#pragma omp barrier	// Syncronize with Slave Thread
			MPI::COMM_WORLD.Barrier();	// Master waits for Slaves

			cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
			cout << "Start Tracer for " << N << " points ... " << endl;
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
				MPI::COMM_WORLD.Recv(recieve.dataFirst(),7*N_slave,MPI::DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,status);
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
				for(i=write_last+1; i<=write_max; i++)
				{
					for(j=1;j<=N_slave;j++)	out << results_all(i,1,j) << "\t" << results_all(i,2,j) << "\t" << results_all(i,3,j) << "\t" << results_all(i,4,j) << "\t" << results_all(i,5,j) << "\t" << results_all(i,6,j) << "\t" << results_all(i,7,j) << endl;
					write_memory(i) = 2;
				}
				write_last = write_max;

			} // end while(recieve_packages <= NoOfPackages)
		} // end omp section: Master thread

		#pragma omp section	//------- Slave Thread on Master Node: does same calculations as Salve Nodes ----------------------------------------------------------------------------------
		{
			// Prepare Perturbation
			prep_perturbation(EQD,PAR,mpi_rank);

			//#pragma omp barrier	// Syncronize with Master Thread

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
				ofs3 << "Node: " << mpi_rank << " works on Package: " << tag << endl;

				Zmin_slave = Z_values((tag-1)*NZ_slave+1);
				Zmax_slave = Z_values(tag*NZ_slave);

				ofs2 << "Start Tracer for " << N_slave << " points ... " << endl;
				for(i=1;i<=N_slave;i++)
				{
					// Set and store initial condition
					if(PAR.create_flag == 3)	// creates regular grid from theta and psi
					{
						FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave,2);	// here Rmin = psimin and Zmin = thetamin; max respectively
						results_all(tag,1,i) = FLT.get_theta();
						results_all(tag,2,i) = FLT.psi;
					}
					else						// creates regular grid from R and Z
					{
						FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave);	// matlab requires R to vary first
						results_all(tag,1,i) = FLT.R;
						results_all(tag,2,i) = FLT.Z;
					}

					// Integration terminates outside of boundary box
					skip_connect = 0;
					if(outofBndy(FLT.R,FLT.Z,EQD) == true)
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

					// Follow fieldline to walls
					if(skip_connect == 0) chk = FLT.connect(ntor,length,psimin,psimax,psiav,PAR.itt,PAR.MapDirection);

					// Store results
					results_all(tag,3,i) = ntor;
					results_all(tag,4,i) = length/1000.0;
					results_all(tag,5,i) = psimin;
					results_all(tag,6,i) = psimax;
					results_all(tag,7,i) = psiav;

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
		for(j=1;j<=N_slave;j++)	out << results_all(i,1,j) << "\t" << results_all(i,2,j) << "\t" << results_all(i,3,j) << "\t" << results_all(i,4,j) << "\t" << results_all(i,5,j) << "\t" << results_all(i,6,j) << "\t" << results_all(i,7,j) << endl;
	}
} // end Master
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	// Prepare Perturbation
	prep_perturbation(EQD,PAR,mpi_rank);

	MPI::COMM_WORLD.Barrier();	// Syncronize with Master

	ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;

	// Result array for Slave
	Array<double,2> results(Range(1,7),Range(1,N_slave));

	// Start working...
	while(1)	
	{
		// Recieve initial conditions to calculate
		MPI::COMM_WORLD.Recv(send_Z_limits.dataFirst(),2,MPI::DOUBLE,0,MPI_ANY_TAG,status);

		tag = status.Get_tag();
		ofs2 << endl;
		ofs2 << "Node: " << mpi_rank << " works on Package: " << tag << endl;
		if(tag == 0) break;	// Still work to do?  ->  tag = 0: NO, stop working

		Zmin_slave = send_Z_limits(1);
		Zmax_slave = send_Z_limits(2);

		ofs2 << "Start Tracer for " << N_slave << " points ... " << endl;
		for(i=1;i<=N_slave;i++)
		{
			// Set and store initial condition
			if(PAR.create_flag == 3)	// creates regular grid from theta and psi
			{
				FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave,2);	// here Rmin = psimin and Zmin = thetamin; max respectively
				results(1,i) = FLT.get_theta();
				results(2,i) = FLT.psi;
			}
			else						// creates regular grid from R and Z
			{
				FLT.set(i,N_slave,PAR.Rmin,PAR.Rmax,Zmin_slave,Zmax_slave,NZ_slave);	// matlab requires R to vary first
				results(1,i) = FLT.R;
				results(2,i) = FLT.Z;
			}

			// Integration terminates outside of boundary box
			skip_connect = 0;
			if(outofBndy(FLT.R,FLT.Z,EQD) == true)
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

			// Follow fieldline to walls
			if(skip_connect == 0) chk = FLT.connect(ntor,length,psimin,psimax,psiav,PAR.itt,PAR.MapDirection);

			// Store results
			results(3,i) = ntor;
			results(4,i) = length/1000.0;
			results(5,i) = psimin;
			results(6,i) = psimax;
			results(7,i) = psiav;

			if(i%100==0) ofs2 << "Trax: " << i << endl;
		} // end for

		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),7*N_slave,MPI::DOUBLE,0,tag);

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
if(PAR.response_field >= 0) m3dc1_unload_file_();
#endif

// MPI finalize
MPI::Finalize();

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

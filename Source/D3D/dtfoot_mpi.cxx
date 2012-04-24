// Program calculates connection length and penetration depth for D3D
// include drifts and time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						18.11.09

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	2d footprint data for colored contour plot
//			log-file


// uses Parallel computation by open-mpi and openmp
// IMPORTANT: Node 0 starts two threads using openmp. The Master-thread contols communications with MPI-Slaves. The Slave-thread on Node 0 does the same calculations as the MPI-Slaves
// Launch Syntax:				mpirun -np [NoOfNodes] -hostfile [FileName] dtfoot_mpi [Parameterfile] [Praefix(optional)] &
// alternative Launch Syntax:	mpirun -np [NoOfNodes] -host [HostNames separated by Commata] dtfoot_mpi [Parameterfile] [Praefix(optional)] &

// Define
//--------
//#define BZ_DEBUG
#define program_name "dtfoot_mpi"
#define USE_MPI

// Include
//--------
#include <andi.hxx>
#include <efit_class.hxx>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <d3d-drift.hxx>
#include <omp.h>

// Prototypes  

// Switches
const int useparfile = 1;	// 0: additional parameters are set in the code		1: All parameters are read from file

// Golbal Parameters 
EFIT EQD;
double GAMMA;
double eps0;
double Ix;

// Main Program
//--------------
int main(int argc, char *argv[])
{
// MPI initialize
MPI::Init(argc, argv);
mpi_rank = MPI::COMM_WORLD.Get_rank();
mpi_size = MPI::COMM_WORLD.Get_size();
if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes or start dtlaminar for a non-parallel run." << endl; EXIT;}

// Variables
int i,j;
int chk;
int MapDirection;	// set within the Code, depending on which_target_plate; 1: positive phi-direction	-1: negative phi-direction
double phi,phistart,psimin;
double t,ntor,length;
double omegac;
Array<double,1> xa(Range(1,2));
Range all = Range::all();

int tag,sender;
double tmin_slave,tmax_slave,dt;
Array<double,1> send_t_limits(Range(1,2));
MPI::Status status;

// Use system time as seed(=idum) for random numbers
double now = zeit();
long idum = long(now);

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

// Read Parameterfile
vector<double> startvec;
if(mpi_rank < 1) cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);

// Set target type for output-filename
LA_STRING type;
switch(which_target_plate)
{
case 0:
	type = "_cp";
	break;
case 1:
	type = "_in";
	break;
case 2:
	type = "_out";
	break;
case 3:
	type = "_shelf";
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
EQD.ReadData(EQD.Shot,EQD.Time);
if(mpi_rank < 1) cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
int itt = 15;
double tmin = 0.76;
double tmax = 0.85;
double phimin = 0;
double phimax = pi2;
int Nt = 50;
int Nphi = 20;

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

int Nt_slave = 1;

if(useparfile==1)
{
	if(mpi_rank < 1) cout << "All parameters are read from file" << endl;
	ofs2 << "All parameters are read from file" << endl;
	itt = int(startvec[1]);
	tmin = startvec[2];
	tmax = startvec[3];
	phimin = startvec[4];
	phimax = startvec[5];
	Nt = int(startvec[6]);
	Nphi = int(startvec[0]);
	MapDirection = int(startvec[8]);
	Ekin = startvec[18];
	lambda = startvec[19];
}
int N = Nt*Nphi;
//if(which_target_plate==2 || which_target_plate==3) MapDirection = 1;	// Only in lower single null discharges
//else  MapDirection = -1;												// Only in lower single null discharges

//if(Nt%Nt_slave != 0) Nt_slave = 1;	// if NZ is not a multiple of NZ_slave then only one NR row is send at each time -> more communication 
int N_slave = Nphi*Nt_slave;
int NoOfPackages = int(Nt/Nt_slave);

if(Nt<=1) {dt=0;}
else dt=(tmax-tmin)/(Nt-1);
if(dt == 0) tmax = tmin;

// Extend lower boundary to prevent horizontal plate to be outside of boundary
bndy[2] = -1.4;	 // originally  bndy[2] = -1.367 and plate at -1.3664-Z0

// Prepare particles
if(sigma==0)
{
	if(mpi_rank < 1) cout << "Field lines are calculated" << endl;
	ofs2 << "Field lines are calculated" << endl;
	GAMMA = 1;	omegac = 1;	 eps0 = 1;	Ix = 1;
}
else
{
	if(Zq>=1) // Ions
	{
		GAMMA = 1 + Ekin/(E0p*Massnumber);	// relativistic gamma factor 1/sqrt(1-v^2/c^2)
		omegac = e*EQD.Bt0/(mp*Massnumber);	// normalized gyro frequency (SI-System)
		if(mpi_rank < 1) cout << "Ions are calculated" << endl;
		ofs2 << "Ions are calculated" << endl;
	}
	else // Electrons
	{
		Zq = -1;	// default!
		GAMMA = 1 + Ekin/(E0e*Massnumber);	// relativistic gamma factor 1/sqrt(1-v^2/c^2)
		omegac = e*EQD.Bt0/(me*Massnumber);	// normalized gyro frequency (SI-System)
		if(mpi_rank < 1) cout << "Electrons are calculated" << endl;
		ofs2 << "Electrons are calculated" << endl;
	}
	eps0 = c*c/omegac/omegac/EQD.R0/EQD.R0;	// normalized rest energy
	Ix = -0.5/double(Zq)*eps0*((lambda*(GAMMA-1)+1)*(lambda*(GAMMA-1)+1)-1);
	if(mpi_rank < 1) cout << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
	ofs2 << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
}

MPI::COMM_WORLD.Barrier();	// Syncronize all Nodes

// Master only (Node 0)
//--------------------------------------------------------------------------------------------------
if(mpi_rank < 1)
{
	// log file
	ofs3.open("log_" + LA_STRING(program_name) + type + praefix + "_Master" + ".dat");
	ofs3.precision(16);

	ofs3 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
	ofs3 << "F-coil: " << useFcoil << "\t" << "C-coil: " << useCcoil << "\t" << "I-coil: " << useIcoil << endl << endl;
	ofs3 << "No. of Packages = " << NoOfPackages << " Points per Package = " << N_slave << endl << endl;

	// additional parameters for IO
	const int psize=10;
	parstruct * parvec = new parstruct[psize];
	parvec[0].name = "Max. Iterations";	parvec[0].wert = itt;
	parvec[1].name = "t-grid";			parvec[1].wert = Nt;
	parvec[2].name = "phi-grid";		parvec[2].wert = Nphi;
	parvec[3].name = "tmin";			parvec[3].wert = tmin;
	parvec[4].name = "tmax";			parvec[4].wert = tmax;
	parvec[5].name = "phimin";			parvec[5].wert = phimin;
	parvec[6].name = "phimax";			parvec[6].wert = phimax;
	parvec[7].name = "MapDirection";	parvec[7].wert = MapDirection;
	parvec[8].name = "Ekin";			parvec[8].wert = Ekin;
	parvec[9].name = "energy ratio lambda";	parvec[9].wert = lambda;

	// Output
	ofstream out(filenameout);
	out.precision(16);
	vector<LA_STRING> var(5);
	var[0] = "phi[rad]";  var[1] = "length t";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";
	writeiodata(out,var,parvec,psize,parfilename);

	// Result array:					Package ID,  Column Number,  Values
	Array<double,3> results_all(Range(1,NoOfPackages),Range(1,5),Range(1,N_slave)); 
	Array<double,2> recieve(Range(1,5),Range(1,N_slave));
	Array<double,2> slice;
	tag = 1;	// first Package

	// t array
	Array<double,1> t_values(Range(1,Nt));
	for(i=1;i<=Nt;i++) t_values(i) = tmin + (i-1)*dt;

	int sent_packages = 0;
	int recieve_packages = 0;

	#pragma omp parallel shared(results_all,t_values,sent_packages,recieve_packages) private(i,tag) num_threads(2)
	{
		#pragma omp sections nowait
		{
		#pragma omp section	//-------- Master Thread: controlles comunication ----------------------------------------------------------------------------------------------------------------------
		{
			#pragma omp barrier	// Syncronize with Slave Thread
			MPI::COMM_WORLD.Barrier();	// Master waits for Slaves

			cout << "Target (0=vertical, 1=45°, 2=horizontal, 3=shelf): " << which_target_plate << endl;
			cout << "Start Tracer for " << N << " points ... " << endl;
			//ofs2 << "Target (0=vertical, 1=45°, 2=horizontal, 3=shelf): " << which_target_plate << endl;
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
					workingNodes -= 1;
					ofs3 << "Send termination signal to Node: " << i << endl;
				}
			} // end for(i=1;i<mpi_size;i++)

			// Recieve results and send more
			while(workingNodes > 0)	// workingNodes > 0: Slave still working -> MPI:Revc needed		workingNodes == 0: all Slaves recieved termination signal
			{
				// Recieve Result
				MPI::COMM_WORLD.Recv(recieve.dataFirst(),5*N_slave,MPI::DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,status);
				sender = status.Get_source();
				tag = status.Get_tag();
				ofs3 << "Recieve from Node: " << sender << " Package: " << tag << endl;

				#pragma omp critical
				{
					recieve_packages += 1;
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

			} // end while(recieve_packages <= NoOfPackages)
		} // end omp section: Master thread

		#pragma omp section	//------- Slave Thread on Master Node: does same calculations as Salve Nodes ----------------------------------------------------------------------------------
		{
			// Prepare Perturbation
			prep_perturbation();

			#pragma omp barrier	// Syncronize with Master Thread

			ofs2 << "Target (0=vertical, 1=45°, 2=horizontal, 3=shelf): " << which_target_plate << endl;

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
					t = start_on_target(i,Nt_slave,Nphi,tmin_slave,tmax_slave,phimin,phimax,xa(2),xa(1),phistart);
					phi = phistart*rTOd;	//phi in deg;

					chk = connect(xa,phi,itt,MapDirection,ntor,length,psimin);

					//Store results
					results_all(tag,1,i) = phistart;
					results_all(tag,2,i) = t;
					results_all(tag,3,i) = ntor;
					results_all(tag,4,i) = length/1000.0;
					results_all(tag,5,i) = psimin;

					if(i%100==0) ofs2 << "Trax: " << i << endl;
				} // end for

				#pragma omp critical
				{
					recieve_packages += 1;
					ofs3 << "------------------------------------ Progress: " << recieve_packages << " of " << NoOfPackages << " completed" <<  endl;
				}

			}// end while
		}// end omp section: Slave Thread
		} // end omp sections
	} // end omp parallel

	// Write Output to file
	for(i=1;i<=NoOfPackages;i++)
	{
		for(j=1;j<=N_slave;j++)	out << results_all(i,1,j) << "\t" << results_all(i,2,j) << "\t" << results_all(i,3,j) << "\t" << results_all(i,4,j) << "\t" << results_all(i,5,j) << endl;
	}
} // end Master
//---------------------------------------------------------------------------------------------------


// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	// Prepare Perturbation
	prep_perturbation();

	MPI::COMM_WORLD.Barrier();	// Syncronize with Master

	ofs2 << "Target (0=vertical, 1=45°, 2=horizontal, 3=shelf): " << which_target_plate << endl;

	// Result array for Slave
	Array<double,2> results(Range(1,5),Range(1,N_slave));

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
			t = start_on_target(i,Nt_slave,Nphi,tmin_slave,tmax_slave,phimin,phimax,xa(2),xa(1),phistart);
			phi = phistart*rTOd;	//phi in deg;

			chk = connect(xa,phi,itt,MapDirection,ntor,length,psimin);

			//Store results
			results(1,i) = phistart;
			results(2,i) = t;
			results(3,i) = ntor;
			results(4,i) = length/1000.0;
			results(5,i) = psimin;

			if(i%100==0) ofs2 << "Trax: " << i << endl;
		}

		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),5*N_slave,MPI::DOUBLE,0,tag);			

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

// MPI finalize
MPI::Finalize();

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

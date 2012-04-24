// Program calculates connection length and penetration depth for NSTX inside the plasma volume
// for ITER-Drift with Time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						16.3.11

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	2d connection length data for colored contour plot
//			log-file


// uses Parallel computation by open-mpi and openmp
// IMPORTANT: Node 0 starts two threads using openmp. The Master-thread contols communications with MPI-Slaves. The Slave-thread on Node 0 does the same calculations as the MPI-Slaves
// Launch Syntax:				mpirun -np [NoOfNodes] -hostfile [FileName] dtlaminar_mpi [Parameterfile] [Praefix(optional)] &
// alternative Launch Syntax:	mpirun -np [NoOfNodes] -host [HostNames separated by Commata] dtlaminar_mpi [Parameterfile] [Praefix(optional)] &

// Define
//--------
//#define BZ_DEBUG
#define program_name "nstxlaminar_mpi"
#define USE_MPI

// Include
//--------
#include <andi.hxx>
#include <efit_class_nstx.hxx>
#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#include <nstx-drift.hxx>
#include <omp.h>

// Prototypes  
int follow(double& R, double& Z, double& phi, double& psi, int MapDirection);

// Switches
const int useparfile = 1;	// 0: additional parameters are set in the code		1: All parameters are read from file
const int spare_interior = 1;	// 0: all points are calculated		1: inside psi=0.9 results are set to fixed values (code runs faster)

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
mpi_rank = MPI::COMM_WORLD.Get_rank();		// declared in d3d-drift.hxx
mpi_size = MPI::COMM_WORLD.Get_size();		// declared in d3d-drift.hxx
if(mpi_size < 2 && mpi_rank < 1) {cout << "Too few Nodes selected. Please use more Nodes or start dtlaminar for a non-parallel run." << endl; EXIT;}

// Variables
int i,j;
int chk,skip_connect;
double R,Z,phi,psi,psimin;
double ntor,length,dummy;
double omegac;
Array<double,1> xa(Range(1,2));
Range all = Range::all();

int tag,sender;
double Zmin_slave,Zmax_slave,dz;
Array<double,1> send_Z_limits(Range(1,2));
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

// Read parameter file
vector<double> startvec;
if(mpi_rank < 1) cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);

// Set area type for output-filename
LA_STRING type;
if(startvec[5] > 2 && startvec[4] > -2) type = "_up";
else type = "";

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + "_Node" + LA_STRING(mpi_rank) + ".dat");
ofs2.precision(16);
ofstream ofs3;

// Output
LA_STRING filenameout = "lam" + type + praefix + ".dat";
if(mpi_rank < 1) outputtest(filenameout);
MPI::COMM_WORLD.Barrier();	// All Nodes wait for Master

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
if(mpi_rank < 1) cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
int itt = 40;
double Rmin = 4 - EQD.RmAxis;
double Rmax = 6 - EQD.RmAxis;
double Zmin = -4.6 - EQD.ZmAxis;
double Zmax = -0.5 - EQD.ZmAxis;
int NR = 200;
int NZ = 300;
double phistart = 0;
int MapDirection = 0;	// 1: positive phi-direction	-1: negative phi-direction	0: both directions

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

int NZ_slave = 1;

if(useparfile==1)
{
	if(mpi_rank < 1) cout << "All parameters are read from file" << endl;
	ofs2 << "All parameters are read from file" << endl;
	itt = int(startvec[1]);
	Rmin = startvec[2] - EQD.RmAxis;
	Rmax = startvec[3] - EQD.RmAxis;
	Zmin = startvec[4] - EQD.ZmAxis;
	Zmax = startvec[5] - EQD.ZmAxis;
	NR = int(startvec[6]);
	NZ = int(startvec[0]);
	phistart = startvec[7];
	MapDirection = int(startvec[8]);
	Ekin = startvec[16];
	lambda = startvec[17];
}
int N = NR*NZ;
//if(NZ%NZ_slave != 0) NZ_slave = 1;	// if NZ is not a multiple of NZ_slave then only one NR row is send at each time -> more communication 
int N_slave = NR*NZ_slave;
int NoOfPackages = int(NZ/NZ_slave);

if(NZ<=1) {dz=0;}
else dz=(Zmax-Zmin)/(NZ-1);
if(dz == 0) Zmax = Zmin;

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
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank < 1)
{
	// log file
	ofs3.open("log_" + LA_STRING(program_name) + praefix + "_Master" + ".dat");
	ofs3.precision(16);

	ofs3 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
	ofs3 << "NSTX-coil (0 = off, 1 = on): " << useIcoil << endl << endl;
	ofs3 << "No. of Packages = " << NoOfPackages << " Points per Package = " << N_slave << endl << endl;

	// additional parameters for IO
	const int psize = 11;
	parstruct * parvec = new parstruct[psize];
	parvec[0].name = "Max. Iterations";	parvec[0].wert = itt;
	parvec[1].name = "R-grid";			parvec[1].wert = NR;
	parvec[2].name = "Z-grid";			parvec[2].wert = NZ;
	parvec[3].name = "Rmin";			parvec[3].wert = Rmin + EQD.RmAxis;
	parvec[4].name = "Rmax";			parvec[4].wert = Rmax + EQD.RmAxis;
	parvec[5].name = "Zmin";			parvec[5].wert = Zmin + EQD.ZmAxis;
	parvec[6].name = "Zmax";			parvec[6].wert = Zmax + EQD.ZmAxis;
	parvec[7].name = "phistart";		parvec[7].wert = phistart;
	parvec[8].name = "MapDirection";	parvec[8].wert = MapDirection;
	parvec[9].name = "energy ratio lambda";	parvec[9].wert = lambda;
	parvec[10].name = "Ekin";			parvec[10].wert = Ekin;

	// Output
	ofstream out(filenameout);
	out.precision(16);
	vector<LA_STRING> var(5);
	var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";
	writeiodata(out,var,parvec,psize,parfilename);

	// Result array:					Package ID,  Column Number,  Values
	Array<double,3> results_all(Range(1,NoOfPackages),Range(1,5),Range(1,N_slave)); 
	Array<double,2> recieve(Range(1,5),Range(1,N_slave));
	Array<double,2> slice;
	tag = 1;	// first Package

	// Z array
	Array<double,1> Z_values(Range(1,NZ));
	for(i=1;i<=NZ;i++) Z_values(i) = Zmin + (i-1)*dz;

	int sent_packages = 0;
	int recieve_packages = 0;

	#pragma omp parallel shared(results_all,Z_values,sent_packages,recieve_packages) private(i,tag) num_threads(2)
	{
		#pragma omp sections nowait
		{
		#pragma omp section	//-------- Master Thread: controlles comunication ----------------------------------------------------------------------------------------------------------------------
		{
			#pragma omp barrier	// Syncronize with Slave Thread
			MPI::COMM_WORLD.Barrier();	// Master waits for Slaves

			cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;
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

			} // end while(recieve_packages <= NoOfPackages)
		} // end omp section: Master thread

		#pragma omp section	//------- Slave Thread on Master Node: does same calculations as Salve Nodes ----------------------------------------------------------------------------------
		{
			// Prepare Perturbation
			prep_perturbation();

			#pragma omp barrier	// Syncronize with Master Thread

			ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;

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
					set(i,N_slave,Zmin_slave,Zmax_slave,Rmin,Rmax,Z,R,NZ_slave);	// swap x and y, because matlab requires R to vary first
					
					// Integration terminates outside of boundary box
					skip_connect = 0;
					if(outofBndy(R+EQD.RmAxis,Z+EQD.ZmAxis,simpleBndy) == true)
					{
						ntor = 0;
						length = 0;
						psimin = 10;
						skip_connect = 1;
					}

					// Toroidal coordinates
					xa(1) = atan(Z/R);	// theta
					if(R<0) xa(1) += pi;
					if(R>=0 && Z<0) xa(1) += pi2;
					xa(2) = sqrt(R*R+Z*Z);	// r
					phi = phistart;

					// Spare the calculation of the interior
					EQD.get_psi(R+EQD.RmAxis,Z+EQD.ZmAxis,psi,dummy,dummy); // get psi
					if(spare_interior == 1 && psi <= 0.9) 
					{
						ntor = 2*itt;
						length = 4000.0;
						psimin = 0.9;
						skip_connect = 1;
					}

					// Follow fieldline to walls
					if(skip_connect == 0) chk = connect(xa,phi,itt,MapDirection,ntor,length,psimin);

					//Store results
					results_all(tag,1,i) = R + EQD.RmAxis;
					results_all(tag,2,i) = Z + EQD.ZmAxis;
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
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Slaves only (Exeption: only one Node is used <=> no MPI or mpi_size = 1)
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(mpi_rank > 0)
{
	// Prepare Perturbation
	prep_perturbation();

	MPI::COMM_WORLD.Barrier();	// Syncronize with Master

	ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;

	// Result array for Slave
	Array<double,2> results(Range(1,5),Range(1,N_slave));

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
			set(i,N_slave,Zmin_slave,Zmax_slave,Rmin,Rmax,Z,R,NZ_slave);	// swap x and y, because matlab requires R to vary first
			
			// Integration terminates outside of boundary box
			skip_connect = 0;
			if(outofBndy(R+EQD.RmAxis,Z+EQD.ZmAxis,simpleBndy) == true) 
			{
				ntor = 0;
				length = 0;
				psimin = 10;
				skip_connect = 1;
			}

			// Toroidal coordinates
			xa(1) = atan(Z/R);	// theta
			if(R<0) xa(1) += pi;
			if(R>=0 && Z<0) xa(1) += pi2;
			xa(2) = sqrt(R*R+Z*Z);	// r
			phi = phistart;

			// Spare the calculation of the interior
			EQD.get_psi(R+EQD.RmAxis,Z+EQD.ZmAxis,psi,dummy,dummy); // get psi
			if(spare_interior == 1 && psi <= 0.9) 
			{
				ntor = 2*itt;
				length = 4000.0;
				psimin = 0.9;
				skip_connect = 1;
			}
			
			// Follow fieldline to walls
			if(skip_connect == 0) chk = connect(xa,phi,itt,MapDirection,ntor,length,psimin);

			//Store results
			results(1,i) = R + EQD.RmAxis;
			results(2,i) = Z + EQD.ZmAxis;
			results(3,i) = ntor;
			results(4,i) = length/1000.0;
			results(5,i) = psimin;

			if(i%100==0) ofs2 << "Trax: " << i << endl;
		} // end for

		// Send results to Master
		MPI::COMM_WORLD.Send(results.dataFirst(),5*N_slave,MPI::DOUBLE,0,tag);			

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

// MPI finalize
MPI::Finalize();

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------- follow ----------------------------
int follow(double& R, double& Z, double& phi, double& psi, int MapDirection)
{
int chk;
double dummy;
double dphi = MapDirection*dpinit/rTOd;
double phi_rad = phi/rTOd;	// phi in rad

Array<double,1> y(2);
y(0) = R;	// R
y(1) = Z;	// Z

chk = rkint(nvar,10,dphi,y,phi_rad);
if(chk<0) return -1;	// particle has left system

// return coordinates
R = y(0); // R 
Z = y(1); // Z 

// phi back in deg
phi = phi_rad*rTOd;

// Get normalized flux
EQD.get_psi(y(0),y(1),psi,dummy,dummy);

return 0;
}


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

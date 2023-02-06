// Header-File for any unspecified machine (and CMOD)
// Only Machine specific subroutines
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						61.2.23

// Define
//--------
#ifndef ANYM_INCLUDED
#define ANYM_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD);

// -------------- global Parameters ---------------------------------------------------------------------------------------
//double bndy[4] = {0.44, 0.91, -0.45, 0.45};	// Boundary
double bndy[4] = {0.001, 13.0, -10.0, 10.0};	// Boundary Box: Rmin, Rmax, Zmin, Zmax for all tokamaks smaller than DEMO

// extern
#ifdef m3dc1
	extern M3DC1 M3D;
#endif

extern Array<double,4> field;
extern fakeIsland FISLD;

extern ofstream ofs2;

//---------------- getBfield ----------------------------------------------------------------------------------------------
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR)
{
int chk;
chk = getBfield_general(R,Z,phi,B_R,B_Z,B_phi,EQD,PAR);
return chk;
}

//---------- prep_perturbation --------------------------------------------------------------------------------------------
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank, LA_STRING supPath)
{
int chk;

#ifdef m3dc1
	// Prepare loading M3D-C1
	if(PAR.response_field >= 0) chk = M3D.read_m3dc1sup(supPath);
	else chk = 0;
#else
	chk = 0;
#endif

// Load machine sup file here to get perturbation coil currents
// None to load here

#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{
		//if(chk == -1) M3D.scale_from_coils(d3icoil_.curntIc, nIloops, nIloops);	// no m3dc1sup.in file found -> scale from diiidsup.in file
		if(chk == -1)
		{
			if(mpi_rank < 1) cout << "m3dc1sup.in not found --> Abort!" << endl;
			EXIT;
		}
		M3D.load(PAR, mpi_rank);
	}
	else
	{
		if(mpi_rank < 1) cout << "Using g-file!" << endl;
		ofs2 << "Using g-file!" << endl;
	}
#else
	if(mpi_rank < 1) cout << "Using g-file!" << endl;
	ofs2 << "Using g-file!" << endl;
#endif

return;
}

//---------------- point_along_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t=[0 1]
// end points of target are explicitly defined here
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD)
{
Array<double,1> p1(Range(1,2)),p2(Range(1,2)),d(Range(1,2));
double R1 = 0,Z1 = 0;	// upper or left Point
double R2 = 0,Z2 = 0;	// lower or right Point

switch(target)
{
//case 1:
//	break;
//case 2:
//	break;
default:
	ofs2 << "No target specified" << endl;	// Here this function does nothing, and will abort the run if called
	EXIT;
	break;
}

p1(1) = R1;	 p1(2) = Z1;
p2(1) = R2;	 p2(2) = Z2;
d = p2 - p1;

// Coordinates
p = p1 + t*d;
}

#endif // ANYM_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

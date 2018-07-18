// Header-File for the C-Mod Programs
// Only Machine specific subroutines
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						21.3.17

// Define
//--------
#ifndef CMOD_INCLUDED
#define CMOD_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD);

// -------------- global Parameters ---------------------------------------------------------------------------------------
double bndy[4] = {0.44, 0.91, -0.45, 0.45};	// Boundary

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
	ofs2 << "No target specified" << endl;
	EXIT;
	break;
}

p1(1) = R1;	 p1(2) = Z1;
p2(1) = R2;	 p2(2) = Z2;
d = p2 - p1;

// Coordinates
p = p1 + t*d;
}

#endif // CMOD_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

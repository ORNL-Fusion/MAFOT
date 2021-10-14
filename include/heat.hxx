// Header-File for the HEAT Drift Programs
// Only Machine specific subroutines
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						22.6.11
// T.Looby						20200421

/*
Note from TL: 20201116
These are the things I changed for HEAT
  1)  added heat header file heat.hxx
  2)  added heat to io_class
  3)  commented cout line in structure main loop (averaged statistics)
  4)  changed size in structure loop from itt*360/dphi to itt
  5)  added heat stuff to makefile
  6)  added HEAT program_name in structure if statement
  7)  added 0 padded shot number to io_class
  8)  changed simpleBndy to 1 in mafot.hxx
  9)  added include <heat.hxx> to mafot.hxx in preprocessor directives conditional
  10) changed point_along_target() to be a dummy function in this file
  11) fixed Fpol exterpolation in mafot.hxx getBfield_general()

for all these edits I have included my initials "TL" in a comment, so you can
search for TL to find the line I am talking about

this enables HEAT to run with a machine agnostic code, while simultaneously
not messing with Andreas' code.  If a user really needs to include machine
specific coils (as may be true for M3DC1 runs?) this may need to be adapted.
Here this is for 2D plasmas (ie EFIT)
*/


// Define
//--------
#ifndef HEAT_INCLUDED
#define HEAT_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");

// ------------ Set Parameters for fortran --------------------------------------------------------------------------------
const int mxbands = 1;
const int mxloops = 6;
const int mxsegs = 26;
const int nturns = 2;

// Output of iterigeom_, set in: prep_perturbation()
int kuse[mxbands][mxloops];
int nbands;
int nloops[mxbands];
int nsegs[mxbands][mxloops];
double xs[mxbands][mxloops][mxsegs][3];
double dvs[mxbands][mxloops][mxsegs][4];
double curntw[mxbands][mxloops];

// ------------------- Fortran Common Blocks ------------------------------------------------------------------------------
extern "C"
{
	extern struct{double pi,twopi,cir,rtd,dtr;} consts_;
	extern struct{double ECadj[mxbands][3];
				  double ECcur[mxbands][mxloops];} currents_;
}

// ----------------- Fortran Routines -------------------------------------------------------------------------------------
extern "C"
{
	void nstxecgeom_(int *kuse, const int *nbands, int nloops[], int *nsegs, double *xs, double *dvs, double *curntw);
	void polygonb_(const int *loopsdim, const int *segsdim, const int *nloops, int nsegs[], int kuse[],
					double *xs, double *dvs, double curnt[],
					double *x, double *y, double *z, double *bx, double *by, double *bz);
}

// -------------- global Parameters ---------------------------------------------------------------------------------------
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
double X,Y,bx,by,bz;
double B_X,B_Y;
double sinp,cosp;

B_R = 0; B_phi = 0; B_Z = 0;

sinp = sin(phi);
cosp = cos(phi);

X = R*cosp;
Y = R*sinp;

chk = getBfield_general(R,Z,phi,B_R,B_Z,B_phi,EQD,PAR);

//Added 20210901 by TL for testing
//cout << "===TESTING GFILE BFIELD INTERPOLATION===" << endl;
//cout << R << endl;
//cout << Z << endl;
//cout << phi << endl;
//cout << X << endl;
//cout << Y << endl;
//cout << B_R << endl;
//cout << B_Z << endl;
//cout << B_phi << endl;

if(chk==-1) {return -1;}

B_X = 0;	B_Y = 0;
// I-coil perturbation field
if(PAR.useIcoil==1)
{
	for(int i=0;i<nbands;i++)
	{
		bx = 0;	by = 0;	bz = 0;

		polygonb_(&mxloops, &mxsegs, &nloops[i], &nsegs[i][0], &kuse[i][0],
					&xs[i][0][0][0], &dvs[i][0][0][0], &curntw[i][0],
					&X, &Y, &Z, &bx, &by, &bz);
		B_X += bx;
		B_Y += by;
		B_Z += bz;
	}
}

// Transform B_perturbation = (B_X, B_Y, B_Z) to cylindrical coordinates and add
B_R += B_X*cosp + B_Y*sinp;
B_phi += -B_X*sinp + B_Y*cosp;

return 0;
}

//---------- prep_perturbation --------------------------------------------------------------------------------------------
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank, LA_STRING supPath)
{
cout << "Preparing HEAT perturbation" << endl;
int i,j;
int chk;
LA_STRING line;	// entire line is read by ifstream
ifstream in;

// Set common blocks parameters
consts_.pi = pi;
consts_.twopi = pi2;
consts_.cir = 360.0;
consts_.rtd = 360.0/pi2;
consts_.dtr = 1.0/consts_.rtd;

#ifdef m3dc1
	// Prepare loading M3D-C1
	if(PAR.response_field >= 0) chk = M3D.read_m3dc1sup(supPath);
	else chk = 0;
#else
	chk = 0;
#endif

// Read heatsub.in file, if coils or M3D-C1 are on
if(PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	in.open(supPath + "heatsup.in");
	if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open heatsup.in file " << endl; EXIT;}

	for(i=1;i<=5;i++) {in >> line;} // Skip 5 lines
	for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) in >> currents_.ECcur[i][j];}		// Read coil currents

	in >> line;	// Skip line
	for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) in >> currents_.ECadj[i][j];}		// Read coil adjustments

	in.close();	// close file
	in.clear();	// reset ifstream for next use
}

//Modified 20191118 by TL
//HEAT does this thousands of times so we don't want to print every time
#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{
		if(chk == -1) M3D.scale_from_coils(currents_.ECcur[0], mxloops, mxloops*mxbands, currents_.ECadj[0][2]);	// no m3dc1sup.in file found -> scale from heatsup.in file, there is only one band
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

if(mpi_rank < 1) cout << "HEAT-coil (0 = off, 1 = on): " << PAR.useIcoil << endl << endl;
ofs2 << "HEAT-coil (0 = off, 1 = on): " << PAR.useIcoil << endl;

// Write Currents to log files (Check if corretly read in)
ofs2 << "Currents:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) ofs2 << currents_.ECcur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "Adjustments:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) ofs2 << currents_.ECadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << endl;

// Set ECoil geometry
//Will probably need to adapt in future for HEAT machine specific coils?
//for now this if is commented
// if(PAR.useIcoil==1) nstxecgeom_(&kuse[0][0],&nbands,&nloops[0],&nsegs[0][0],&xs[0][0][0][0],&dvs[0][0][0][0],&curntw[0][0]);
}

//---------------- point_along_target ----------------------------------------------------------------------------------------
// We use HEAT to do geometry, not MAFOT, so this is just a dummy function now
// Returns garbage
// Modofied by TL 20210114
//
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t in m
// inner target: t = Z
// outer target: t = R
// points of targets are set to fixed values here
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD)
{
target=1;
switch(target)
{
default:
	ofs2 << "start_on_target: No target specified, because using HEAT" << endl;
	EXIT;
	break;
}
}


#endif // HEAT_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

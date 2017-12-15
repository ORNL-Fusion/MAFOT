// Header-File for the MAST Programs
// Only Machine specific subroutines
// uses Nate Ferraro's M3D-C1 plasma response code output, fixed filename: C1.h5
// Plasma response can be for Equilibrium, or I-coils, or both
// ECC-coils are not yet included in Plasma response
// ++++++ IMPORTANT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// uses separate Libraries for the M3D-C1 routines by Nate Ferraro
// use -Dm3dc1 when compiling -> this define activates this part of the code
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						12.2.13

// Define
//--------
#ifndef MAST_INCLUDED
#define MAST_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
//void IO::readiodata(char* name, int mpi_rank);								// declared in IO class, defined here
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT);

// ------------ Set Parameters for fortran --------------------------------------------------------------------------------
const int mxIbands = 2;
const int mxIloops = 12;
const int mxIsegs = 100;
const int niturns = 4;

const int mxECbands = 1;
const int mxECloops = 4;
const int mxECsegs = 26;
const int ncturns = 3;

// Output of mastigeom_, set in: prep_perturbation()
int kuseI[mxIbands][mxIloops];
int nIbands;
int nIloops[mxIbands];
int nIsegs[mxIbands][mxIloops];
double xsI[mxIbands][mxIloops][mxIsegs][3];
double dvsI[mxIbands][mxIloops][mxIsegs][4];
double curntwI[mxIbands][mxIloops];

// Output of mastecgeom_, set in: prep_perturbation()
int kuseEC[mxECbands][mxECloops];
int nECbands;
int nECloops[mxECbands];
int nECsegs[mxECbands][mxECloops];
double xsEC[mxECbands][mxECloops][mxECsegs][3];
double dvsEC[mxECbands][mxECloops][mxECsegs][4];
double curntwEC[mxECbands][mxECloops];

// ------------------- Fortran Common Blocks ------------------------------------------------------------------------------
extern "C" 
{
	extern struct{double pi,twopi,cir,rtd,dtr;} consts_;
	extern struct{double Iadj[mxIbands][4];
				  double Icur[mxIbands][mxIloops];} icurrents_;
	extern struct{double ECadj[mxECbands][3];
				  double ECcur[mxECbands][mxECloops];} eccurrents_;
}

// ----------------- Fortran Routines -------------------------------------------------------------------------------------
extern "C" 
{
	void mastigeom_(int *kuse, const int *nbands, int nloops[], int *nsegs, double *xs, double *dvs, double *curntw);
	void mastecgeom_(int *kuse, const int *nbands, int nloops[], int *nsegs, double *xs, double *dvs, double *curntw);
	void polygonb_(const int *loopsdim, const int *segsdim, const int *nloops, int nsegs[], int kuse[],
					double *xs, double *dvs, double curnt[], 
					double *x, double *y, double *z, double *bx, double *by, double *bz);
}

// -------------- global Parameters ---------------------------------------------------------------------------------------
double bndy[4] = {0.195, 1.9, -1.8251, 1.8251};	// Boundary; EFIT boundary = {0.06, 2, -2, 2}

// extern
#ifdef USE_SIESTA
	extern SIESTA SIES;
#endif
#ifdef USE_XFIELD
	extern XFIELD XPND;
#endif
#ifdef m3dc1
	extern M3DC1 M3D;
#endif

extern Array<double,4> field;
extern fakeIsland FISLD;

extern ofstream ofs2;

// ---------------------- IO Member functions -----------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------------------------------

// ----------------- readiodata -------------------------------------------------------------------------------------------
void IO::readiodata(char* name, int mpi_rank)
{
LA_STRING input;	// !!! LA_STRING reads entire line, string reads only one word !!!

// Get ShotNr and ShotTime from Parameterfile
ifstream in;
in.open(name);
if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open file " << name << endl; EXIT;}
in >> input;	// Skip first line
in >> input;	// second line gives shot number and time
EQDr.Shot = input.mid(9,6);	// 6 characters starting at index 9 of input string
EQDr.Time = input.mid(22,4); // 4 characters starting at index 22 of input string

// Get Path to g-File from Parameterfile (optional), Linux only!!!
in >> input;	
if(input[1] == '#') 
{
	input = input.mid(input.indexOf('/'));	// all characters of input string starting from index of first '/' in string
	if(input.right(1) != '/') input = input.left(input.length()-1);			// last char in string can be '\r' (Carriage return) --> Error!;  this removes last char if necessary 
	if(input.indexOf(' ') > 1) input = input.left(input.indexOf(' ')-1);		// blanks or comments that follow path are removed from string 
	EQDr.Path = input;
}
in.close();

// Get Parameters
vector<double> vec;
readparfile(name,vec);

// private Variables
filename = name;

// Map Parameters
itt = int(vec[1]);
phistart = vec[7];
MapDirection = int(vec[8]);

// t grid (footprints only)
Nt = int(vec[6]); 
tmin = vec[2]; 
tmax = vec[3];			

// phi grid (footprints only)
Nphi = int(vec[0]); 
phimin = vec[4]; 
phimax = vec[5];		

// R grid (laminar only)
NR = int(vec[6]); 
Rmin = vec[2]; 
Rmax = vec[3];		

// Z grid (laminar only)
NZ = int(vec[0]); 
Zmin = vec[4]; 
Zmax = vec[5];		

// r grid
N = int(vec[6]); 
Nr = int(sqrt(N)); 
rmin = vec[2]; 
rmax = vec[3];		

// theta grid
Nth = int(sqrt(N)); 
thmin = vec[4]; 
thmax = vec[5];		

// Particle Parameters
Ekin = vec[18];
lambda = vec[19];
verschieb = vec[0];

// Set switches
which_target_plate = int(vec[11]);
create_flag = int(vec[12]);
useCcoil = int(vec[13]);
useIcoil = int(vec[14]);
sigma = int(vec[16]);
Zq = int(vec[17]);
useFilament = int(vec[15]);

// M3D-C1 parameter
response = int(vec[9]);
response_field = int(vec[10]);

// Set unused Parameters to defaults
useFcoil = 0;
useTprofile = 0;
}

//------------ End of IO Member functions ---------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

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
if(chk==-1) {return -1;}

B_X = 0;	B_Y = 0;
// ECC-coil perturbation field
if(PAR.useCcoil==1)
{
	for(int i=0;i<nECbands;i++)
	{
		bx = 0;	by = 0;	bz = 0;

		polygonb_(&mxECloops, &mxECsegs, &nECloops[i], &nECsegs[i][0], &kuseEC[i][0],
					&xsEC[i][0][0][0], &dvsEC[i][0][0][0], &curntwEC[i][0],
					&X, &Y, &Z, &bx, &by, &bz);
		B_X += bx;
		B_Y += by;
		B_Z += bz;
	}
}

// I-coil perturbation field
if(PAR.useIcoil==1)
{
	for(int i=0;i<nIbands;i++)
	{
		bx = 0;	by = 0;	bz = 0;

		polygonb_(&mxIloops, &mxIsegs, &nIloops[i], &nIsegs[i][0], &kuseI[i][0],
					&xsI[i][0][0][0], &dvsI[i][0][0][0], &curntwI[i][0],
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

// Read mastsub.in file, if coils or M3D-C1 are on
if(PAR.useCcoil == 1 || PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	in.open(supPath + "mastsup.in");
	if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open mastsup.in file " << endl; EXIT;}

	for(i=1;i<=5;i++) in >> line;	// Skip 5 lines
	for(i=0;i<mxECbands;i++) {for(j=0;j<mxECloops;j++) in >> eccurrents_.ECcur[i][j];}	// Read EC-coil currents

	in >> line;	// Skip line
	for(i=0;i<mxECbands;i++) {for(j=0;j<3;j++) in >> eccurrents_.ECadj[i][j];}			// Read EC-coil adjustments

	in >> line;	// Skip line
	for(i=0;i<mxIbands;i++) {for(j=0;j<mxIloops;j++) in >> icurrents_.Icur[i][j];}		// Read I-coil currents

	in >> line;	// Skip line
	for(i=0;i<mxIbands;i++) {for(j=0;j<4;j++) in >> icurrents_.Iadj[i][j];}				// Read I-coil adjustments

	in.close();	// close file
	in.clear();	// reset ifstream for next use
}

#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{
		if(chk == -1) M3D.scale_from_coils(icurrents_.Icur[1], mxIloops, mxIloops*mxIbands, icurrents_.Iadj[1][3]);	// no m3dc1sup.in file found -> scale from mastsup.in file, use only the lower band for scaling
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

if(mpi_rank < 1) cout << "ECC-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << endl << endl;
ofs2 << "ECC-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << endl << endl;

// Write Currents to log files (Check if corretly read in)
ofs2 << "EC-Currents:" << endl;
for(i=0;i<mxECbands;i++) {for(j=0;j<mxECloops;j++) ofs2 << eccurrents_.ECcur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "EC-Adjustments:" << endl;
for(i=0;i<mxECbands;i++) {for(j=0;j<3;j++) ofs2 << eccurrents_.ECadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << "I-Currents:" << endl;
for(i=0;i<mxIbands;i++) {for(j=0;j<mxIloops;j++) ofs2 << icurrents_.Icur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "I-Adjustments:" << endl;
for(i=0;i<mxIbands;i++) {for(j=0;j<4;j++) ofs2 << icurrents_.Iadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << endl;

// Set EC-coil geometry
if(PAR.useCcoil==1) mastecgeom_(&kuseEC[0][0],&nECbands,&nECloops[0],&nECsegs[0][0],&xsEC[0][0][0][0],&dvsEC[0][0][0][0],&curntwEC[0][0]);

// Set I-coil geometry
if(PAR.useIcoil==1) mastigeom_(&kuseI[0][0],&nIbands,&nIloops[0],&nIsegs[0][0],&xsI[0][0][0][0],&dvsI[0][0][0][0],&curntwI[0][0]);
}

//---------------- start_on_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t=[0 1]
// end points of target are explicitly defined here
// Position (R0,Z0) of magnetic axis is required
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT)
{
int i_p = 0;
int i_phi = 0;
int N = Np*Nphi;
int target;
double dp,dphi,t;
Array<double,1> p1(Range(1,2)),p2(Range(1,2)),p(Range(1,2)),d(Range(1,2));

// Magnetic Axis
//const double R0=EQD.RmAxis;
//const double Z0=EQD.ZmAxis;

// Grid stepsizes and t
if(Np == 1) tmax = tmin;
if(Nphi == 1) phimax = phimin;
if(N <= 1) {dp = 0; dphi = 0;}
else
{
	dp = (tmax - tmin)/(N-1);
	dphi = (phimax - phimin)/(N-1);
}
if(dp == 0) i_phi = i - 1;
if(dphi == 0) i_p = i - 1;
if(dp != 0 && dphi != 0)
{
	dp = (tmax - tmin)/double(Np-1);
	dphi = (phimax - phimin)/double(Nphi-1);
	i_phi = (i - 1)%Nphi;
	i_p = int(double(i-1)/double(Nphi));
}
t = tmin + i_p*dp;

// Postion of Target-Plate
double R1 = 0, Z1 = 0;	// upper or left Point
double R2 = 0, Z2 = 0;	// lower or right Point
switch(target)
{
case 1:	// inner target plate
	R1 = 0.28;	Z1 = -1.6835;
	R2 = 0.28;	Z2 = -1.22909;
	break;
case 2:	// outer target plate
	R1 = 0.7835;	Z1 = -1.825;
	R2 = 1.9;		Z2 = -1.825;
	break;
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

FLT.R = p(1);
FLT.Z = p(2);
FLT.phi = (phimin + dphi*i_phi)*rTOd;	// phi in deg
FLT.get_psi(p(1),p(2),FLT.psi);

FLT.Lc = 0;
FLT.psimin = 10;

if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}
return t;
}

#endif // MAST_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

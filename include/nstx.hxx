// Header-File for the NSTX Drift Programs 
// Only Machine specific subroutines
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						22.6.11

// Define
//--------
#ifndef NSTX_INCLUDED
#define NSTX_INCLUDED

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
double bndy[4] = {0.185, 1.57, -1.63, 1.63};	// Boundary Box: Rmin, Rmax, Zmin, Zmax	

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
//if(vec.size()<19) {if(mpi_rank < 1) cout << "Fail to read all parameters, File incomplete" << endl; EXIT;}

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
Ekin = vec[16];
lambda = vec[17];
verschieb = vec[0];

// Set switches
which_target_plate = int(vec[9]);
create_flag = int(vec[10]);
useFcoil = int(vec[11]);
useCcoil = int(vec[11]);
useIcoil = int(vec[11]);
useFilament = int(vec[12]);
useTprofile = int(vec[13]);
sigma = int(vec[14]);
Zq = int(vec[15]);

// M3D-C1 parameter
response = int(vec[18]);
response_field = int(vec[19]);
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

// Read nstxsub.in file, if coils or M3D-C1 are on
if(PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	in.open(supPath + "nstxsup.in");
	if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open nstxsup.in file " << endl; EXIT;}

	for(i=1;i<=5;i++) {in >> line;} // Skip 5 lines
	for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) in >> currents_.ECcur[i][j];}		// Read coil currents

	in >> line;	// Skip line
	for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) in >> currents_.ECadj[i][j];}		// Read coil adjustments

	in.close();	// close file
	in.clear();	// reset ifstream for next use
}

#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{
		if(chk == -1) M3D.scale_from_coils(currents_.ECcur[0], mxloops, mxloops*mxbands, currents_.ECadj[0][2]);	// no m3dc1sup.in file found -> scale from nstxsup.in file, there is only one band
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

if(mpi_rank < 1) cout << "NSTX-coil (0 = off, 1 = on): " << PAR.useIcoil << endl << endl;
ofs2 << "NSTX-coil (0 = off, 1 = on): " << PAR.useIcoil << endl;

// Write Currents to log files (Check if corretly read in)
ofs2 << "Currents:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) ofs2 << currents_.ECcur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "Adjustments:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) ofs2 << currents_.ECadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << endl;

// Set ECoil geometry
if(PAR.useIcoil==1) nstxecgeom_(&kuse[0][0],&nbands,&nloops[0],&nsegs[0][0],&xs[0][0][0][0],&dvs[0][0][0][0],&curntw[0][0]);
}

//---------------- start_on_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t in m
// inner target: t = Z
// outer target: t = R
// points of targets are set to fixed values here
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

// Grid stepsizes and t
if(Np == 1) tmax = tmin;
if(Nphi == 1) phimax = phimin;
if(N<=1) {dp = 0; dphi = 0;}
else
{
	dp = (tmax-tmin)/(N-1);
	dphi = (phimax-phimin)/(N-1);
}
if(dp==0) i_phi = i-1;
if(dphi==0) i_p = i-1;
if(dp!=0 && dphi!=0) 
{
	dp = (tmax-tmin)/double(Np-1);
	dphi = (phimax-phimin)/double(Nphi-1);
	i_phi = (i-1)%Nphi;
	i_p = int(double(i-1)/double(Nphi));
}
t = tmin + i_p*dp;	// t in m

// Postion of Target-Plate
double R1,Z1;	// (inner,lower) or (outer, left) target Corner
double R2,Z2;	// (inner,upper) or (outer, right) target Corner

target = PAR.which_target_plate;
if(PAR.which_target_plate == 4 && t > 0.5712) target = 5;
if(PAR.which_target_plate == 2 && t > 0.5712) target = 6;
if(PAR.which_target_plate == 40 && t > 0.5712) target = 5;
if(PAR.which_target_plate == 20 && t > 0.5712) target = 6;
if(PAR.which_target_plate == 10 && t < 1.27) target = 80;
if(PAR.which_target_plate == 30 && t > -1.27) target = 70;

switch(target)
{
case 1:	// upper inner target plate, length 40.66 cm
	R1 = 0.2794;	Z1 = 1.1714;
	R2 = 0.2794;	Z2 = 1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 2:	// upper outer target plate, length 27.33 cm
	R1 = 0.2979;	Z1 = 1.6034;
	R2 = 0.5712;	Z2 = 1.6034;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 3:	// lower inner target plate, length 40.66 cm
	R1 = 0.2794;	Z2 = -1.1714;
	R2 = 0.2794;	Z1 = -1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 4:	// lower outer target plate, length 27.33 cm
	R1 = 0.2979;	Z1 = -1.6034;
	R2 = 0.5712;	Z2 = -1.6034;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 5:	// lower outer inclined wall at R > 0.6 and small horizontal piece before that
	R1 = 0.617;		Z1 = -1.628;
	R2 = 1.0433;	Z2 = -1.4603;
	p(1) = t;
	if(t < R1) p(2) = Z1;
	else p(2) = (Z2-Z1)/(R2-R1)*t + (Z1*R2-R1*Z2)/(R2-R1);
	if(t < 0.5712 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
case 6:	// upper outer declined wall at R > 0.6
	R1 = 0.617;		Z1 = 1.628;
	R2 = 1.0433;	Z2 = 1.4603;
	p(1) = t;
	if(t < R1) p(2) = Z1;
	else p(2) = (Z2-Z1)/(R2-R1)*t + (Z1*R2-R1*Z2)/(R2-R1);
	if(t < 0.5712 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
case 10:	// NSTX-U upper inner target plate
	R1 = 0.4150;	Z1 = 1.27;
	R2 = 0.4150;	Z2 = 1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 20:	// NSTX-U upper outer target plate
	R1 = 0.4350;	Z1 = 1.6234;
	R2 = 0.5712;	Z2 = 1.6234;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 30:	// NSTX-U lower inner target plate
	R1 = 0.4150;	Z2 = -1.27;
	R2 = 0.4150;	Z1 = -1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 40:	// NSTX-U lower outer target plate
	R1 = 0.4350;	Z1 = -1.6234;
	R2 = 0.5712;	Z2 = -1.6234;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 70:	// NSTX-U lower inner inclined wall at R < 0.4
	R1 = 0.415;		Z1 = -1.27;
	R2 = 0.315;		Z2 = -1.05;
	p(2) = t;
	p(1) = (R2-R1)/(Z2-Z1)*t + (R1*Z2-Z1*R2)/(Z2-Z1);
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
case 80:	// NSTX-U upper inner inclined wall at R < 0.4
	R1 = 0.415;		Z1 = 1.27;
	R2 = 0.315;		Z2 = 1.05;
	p(2) = t;
	p(1) = (R2-R1)/(Z2-Z1)*t + (R1*Z2-Z1*R2)/(Z2-Z1);
	if(t < Z2 || t > Z1){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
default:
	ofs2 << "start_on_target: No target specified" << endl;
	EXIT;
	break;
}

// Coordinates
FLT.R = p(1);
FLT.Z = p(2);
FLT.phi = (phimin + dphi*i_phi)*rTOd;	// phi in deg
FLT.get_psi(p(1),p(2),FLT.psi);

FLT.Lc = 0;
FLT.psimin = 10;

if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}
return t;
}

#endif // NSTX_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

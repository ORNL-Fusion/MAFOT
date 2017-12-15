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
//void IO::readiodata(char* name, int mpi_rank);								// declared in IO class, defined here
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT);

// -------------- global Parameters ---------------------------------------------------------------------------------------
double bndy[4] = {0.44, 0.91, -0.45, 0.45};	// Boundary

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
EQDr.Shot = input.mid(9,10);	// 10 characters starting at index 9 of input string
EQDr.Time = input.mid(26,4); 	// 4 characters starting at index 26 of input string

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
if(vec.size()<22) {if(mpi_rank < 1) cout << "Fail to read all parameters, File incomplete" << endl; EXIT;}

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
useFcoil = int(vec[13]);
useCcoil = int(vec[14]);
useIcoil = int(vec[15]);
sigma = int(vec[16]);
Zq = int(vec[17]);
useFilament = int(vec[20]);

// External fields parameter
response = int(vec[9]);
response_field = int(vec[10]);

if(vec[21]>1) useTprofile = 0;
else useTprofile = int(vec[21]);
}

//------------ End of IO Member functions ---------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

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

//---------------- start_on_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t=[0 1]
// end points of target are explicitly defined here
// Position (R0,Z0) of magnetic axis is required
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT)
{
int i_p=0;
int i_phi=0;
int N=Np*Nphi;
int target;
double dp,dphi,t;
Array<double,1> p1(Range(1,2)),p2(Range(1,2)),p(Range(1,2)),d(Range(1,2));

// Magnetic Axis
//const double R0=EQD.RmAxis;
//const double Z0=EQD.ZmAxis;

// Grid stepsizes and t
if(Np == 1) tmax = tmin;
if(Nphi == 1) phimax = phimin;
if(N<=1) {dp=0; dphi=0;}
else
{
	dp=(tmax-tmin)/(N-1);
	dphi=(phimax-phimin)/(N-1);
}
if(dp==0) i_phi=i-1;
if(dphi==0) i_p=i-1;
if(dp!=0 && dphi!=0) 
{
	dp=(tmax-tmin)/double(Np-1);
	dphi=(phimax-phimin)/double(Nphi-1);
	i_phi=(i-1)%Nphi;
	i_p=int(double(i-1)/double(Nphi));
}
t=tmin+i_p*dp;
if(PAR.which_target_plate==1 && t<0) target=0;
else target=PAR.which_target_plate;

// Postion of Target-Plate
double R1 = 0,Z1 = 0;	// upper or left Point
double R2 = 0,Z2 = 0;	// lower or right Point
switch(target)
{
default:
	ofs2 << "No target specified" << endl;
	EXIT;
	break;
}
p1(1) = R1;	 p1(2) = Z1;
p2(1) = R2;	 p2(2) = Z2;
d = p2 - p1;
if(target == 0) d(2) *= -1;

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

#endif // CMOD_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

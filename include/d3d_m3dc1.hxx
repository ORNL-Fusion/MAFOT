// Header-File for the DIII-D Programs 
// Only Machine specific subroutines
// uses Nate Ferraro's M3D-C1 plasma response code output, fixed filename: C1.h5
// Plasma response can be for Equilibrium, or I-coils, or both
// C-coils and F-coils are not yet included in Plasma response
// ++++++ IMPORTANT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Libraries for the M3D-C1 routines only exist in Nate's u-drive account at GA
// use -Dm3dc1 when compiling -> this define activates this part of the code
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						3.5.12

// Define
//--------
#ifndef D3D_M3DC1_INCLUDED
#define D3D_M3DC1_INCLUDED

// Include
//--------
#include <m3dc1_class.hxx>
#include <fakeIsland_class.hxx>

// --------------- Prototypes ---------------------------------------------------------------------------------------------
//void IO::readiodata(char* name, int mpi_rank);								// declared in IO class, defined here
//void IO::writeiodata(ofstream& out, double bndy[], vector<LA_STRING>& var);	// declared in IO class, defined here

bool outofBndy(double x, double y, EFIT& EQD);
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT);

// ------------ Set Parameters for fortran --------------------------------------------------------------------------------
const int nFc = 18;
const int nFlps = 2*nFc;
const int nIloops = 12;
const int nIsegs = 14;
const int nCloops = 6;
const int nCsegs = 10;

// Global Variables: have to be known during integration for perturbations, set in: prep_perturbation()
int kuseF[nFlps];
int kuseC[nCloops];
int kuseI[nIloops];
int nccsegs[nCloops];
int nicsegs[nIloops];

// ------------------- Fortran Common Blocks ------------------------------------------------------------------------------
extern "C" 
{
	extern struct{double pi,twopi,cir,rtd,dtr;} consts_;
	extern struct{double fcshft[nFc],fashft[nFc],fctilt[nFc],fatilt[nFc],fcur[nFc];} d3pfer_;
	extern struct{double amat[nFlps][3][3];
				  double origin[nFlps][3]; 
				  double rcoil[nFlps], curlps[nFlps];} d3pflps_;
	extern struct{double dsbp,alfsbp,dthbp,alfthbp;
				  int iplasbp,ipdir,lbpol;} eqpol_;

	extern struct{double curntIc[nIloops]; 
				  double addanglIU, addanglIL, scaleIU, scaleIL;} d3icoil_;
	extern struct{double xis[nIloops][nIsegs][3]; 
				  double divs[nIloops][nIsegs][4];} d3iloops_;

	extern struct{double curntC[nCloops], curntw[nCloops]; 
				  double addanglC, scaleC;} d3ccoil_;
	extern struct{double xcs[nCloops][nCsegs][3];
				  double dcvs[nCloops][nCsegs][4];} d3cloops_;
}

// ----------------- Fortran Routines -------------------------------------------------------------------------------------
extern "C" 
{
	void d3pfgeom_(int kuse[]);
	void d3igeom_(int kuse[]);
	void d3cgeom_(int kuse[]);
	void d3pferrb_(int kuse[], double *x, double *y, double *z, double *bxf, double *byf, double *bzf);
	void polygonb_(const int *loopsdim, const int *segsdim, const int *nloops, int nsegs[], int kuse[],
					double *xs, double *dvs, double curnt[], 
					double *x, double *y, double *z, double *bx, double *by, double *bz);
}

// -------------- global Parameters ---------------------------------------------------------------------------------------
int simpleBndy = 0;		// 0: use real wall as boundaries, 1: use simple boundary box
double bndy[4] = {1.0, 2.4, -1.367, 1.36};	// Boundary

Array<double,4> field;	// default constructed

#ifdef USE_SIESTA
	SIESTA SIES;
#endif
#ifdef USE_XFIELD
	XFIELD XPND;
#endif

M3DC1 M3D;
fakeIsland FISLD;

// ------------------ log file --------------------------------------------------------------------------------------------
ofstream ofs2;

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

// M3D-C1 parameter
response = int(vec[9]);
response_field = int(vec[10]);

//if(vec[9]>1) {cout << "Plasma Response Parameters are missing in file -> Abort!" << endl; EXIT;}

if(vec[21]>1) useTprofile = 0;
else useTprofile = int(vec[21]);
}

// ------------------- writeiodata ----------------------------------------------------------------------------------------
void IO::writeiodata(ofstream& out, double bndy[], vector<LA_STRING>& var)
{
int i;
out << "# " << program_name << endl;
out << "#-------------------------------------------------" << endl;
out << "### Parameterfile: " << filename << endl;
out << "# Shot: " << EQDr.Shot << endl;
out << "# Time: " << EQDr.Time << endl;
out << "#-------------------------------------------------" << endl;
out << "### M3D-C1:" << endl;
out << "# Plasma response (0=no, >1=yes): " << response << endl;
out << "# Field (-1=M3D-C1 off, 0=Eq, 1=I-coil, 2=both): " << response_field << endl;
out << "#-------------------------------------------------" << endl;
out << "### Switches:" << endl;
out << "# F-coil active (0=no, 1=yes): " << useFcoil << endl;
out << "# C-coil active (0=no, 1=yes): " << useCcoil << endl;
out << "# I-coil active (0=no, 1=yes): " << useIcoil << endl;
out << "# No. of current filaments (0=none): " << useFilament << endl;
out << "# Use Temperature Profile (0=off, 1=on): " << useTprofile << endl;
out << "# Target (0=cp, 1=inner, 2=outer, 3=shelf): " << which_target_plate << endl;
out << "# Create Points (0=r-grid, 1=r-random, 2=target, 3=psi-grid, 4=psi-random, 5=RZ-grid): " << create_flag << endl;
out << "# Direction of particles (1=co-pass, -1=count-pass, 0=field lines): " << sigma << endl;
out << "# Charge number of particles (=-1:electrons, >=1:ions): " << Zq << endl;
out << "# Boundary (0=Wall, 1=Box): " << simpleBndy << endl;
out << "#-------------------------------------------------" << endl;
out << "### Global Parameters:" << endl;
out << "# Steps till Output (ilt): " << ilt << endl;
out << "# Step size (dpinit): " << dpinit << endl;
out << "# Boundary Rmin: " << bndy[0] << endl;
out << "# Boundary Rmax: " << bndy[1] << endl;
out << "# Boundary Zmin: " << bndy[2] << endl;
out << "# Boundary Zmax: " << bndy[3] << endl;
out << "# Magnetic Axis: R0: " << EQDr.RmAxis << endl;
out << "# Magnetic Axis: Z0: " << EQDr.ZmAxis << endl;
out << "#-------------------------------------------------" << endl;
out << "### additional Parameters:" << endl;
for(i=0;i<psize;++i)
{
	out << "# " << pv[i].name << ": " << pv[i].wert << endl;
}
out << "#-------------------------------------------------" << endl;
out << "### Data:" << endl;
out << "# ";
for(i=0;i<int(var.size());i++) out << var[i] << "     ";
out << endl;
out << "#" << endl;
}

//------------ End of IO Member functions ---------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------- outofBndy ---------------------------------------------------------------------------------------------------
// Check if (x,y) is out of the torus. Returns 0 if (x,y) 
// is in boundary an 1 if (x,y) is out of boundary. 
// simpleBndy = 0; use real wall as boundaries
// simpleBndy = 1: use simple boundary box
bool outofBndy(double x, double y, EFIT& EQD)
{
if(isnan(x) || isnan(y) || isinf(x) || isinf(y)) return false;
switch(simpleBndy)
{
case 0:
	return outofRealBndy(x,y,EQD);
	break;
case 1:
	if(x<bndy[0] || x>bndy[1] || y<bndy[2] || y>bndy[3]) return true;	//  bndy[4] = {1.0, 2.4, -1.367, 1.36};
	break;
default:
    cout << "simpleBndy switch has a wrong value!" << endl;
}
return false;
}

//---------------- getBfield ----------------------------------------------------------------------------------------------
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR)
{
int i, chk, chk2;
double psi,dpsidr,dpsidz;
double F;
double X,Y,bx,by,bz;
double B_X,B_Y;
double sinp,cosp;
double coord[3], b_field[3];

coord[0] = R; coord[1] = phi; coord[2] = Z;
B_R = 0; B_phi = 0; B_Z = 0;

sinp = sin(phi);
cosp = cos(phi);

X = R*cosp;
Y = R*sinp;

// Equilibrium field
switch(PAR.response_field)
{
#ifdef USE_XFIELD
case -3:
	XPND.get_B(R, phi, Z, B_R, B_phi, B_Z);
	break;
#endif
#ifdef USE_SIESTA
case -2:
	SIES.get_B(R, phi, Z, B_R, B_phi, B_Z);
	break;
#endif
case -1: case 1: case -10:	// Vacuum equilibrium field from g file
	// get normalized poloidal Flux psi (should be chi in formulas!)
	chk = EQD.get_psi(R,Z,psi,dpsidr,dpsidz);
	if(chk==-1) {ofs2 << "Point is outside of EFIT grid" << endl; B_R=0; B_Z=0; B_phi=1; return -1;}	// integration of this point terminates

	// Equilibrium field
	F = EQD.get_Fpol(psi);
	B_R = dpsidz/R;
	B_phi = F/R;	//B_phi = EQD.Bt0*EQD.R0/R;
	B_Z = -dpsidr/R;
	break;

case 0: case 2: 	// M3D-C1: equilibrium field or total field
	for(i=0;i<M3D.nfiles;i++)
	{
		chk = fio_eval_field(M3D.imag[i], coord, b_field);
		if(chk != 0) // field eval failed, probably outside of M3DC1 domain -> fall back to g-file equilibrium
		{
			chk2 = EQD.get_psi(R,Z,psi,dpsidr,dpsidz);
			if(chk2 == -1) {ofs2 << "Point is outside of EFIT grid" << endl; B_R=0; B_Z=0; B_phi=1; return -1;}	// integration of this point terminates

			// Equilibrium field
			F = EQD.get_Fpol(psi);
			B_R = dpsidz/R;
			B_phi = F/R;	//B_phi = EQD.Bt0*EQD.R0/R;
			B_Z = -dpsidr/R;
			break;	// break the for loop
		}
		else
		{
			B_R += b_field[0];
			B_phi += b_field[1];
			B_Z += b_field[2];
		}
	}
	break;
}

// M3D-C1: I-coil perturbation field only, coils are turned off in prep_perturbation
if(PAR.response_field == 1)
{
	for(i=0;i<M3D.nfiles;i++)
	{
		coord[1] = phi + M3D.phase[i];
		chk = fio_eval_field(M3D.imag[i], coord, b_field);
		if(chk != 0) {b_field[0] = 0; b_field[1] = 0; b_field[2] = 0; break;}
		B_R += b_field[0];
		B_phi += b_field[1];
		B_Z += b_field[2];
		//coord[1] = phi;
	}
}

B_X = 0;	B_Y = 0;
// F-coil perturbation field
bx = 0;	by = 0;	bz = 0;
if(PAR.useFcoil==1) d3pferrb_(&kuseF[0], &X, &Y, &Z, &bx, &by, &bz);
B_X += bx;
B_Y += by;
B_Z += bz;

// C-coil perturbation field
bx = 0;	by = 0;	bz = 0;
if(PAR.useCcoil==1) polygonb_(&nCloops, &nCsegs, &nCloops, &nccsegs[0], &kuseC[0],
						  &d3cloops_.xcs[0][0][0], &d3cloops_.dcvs[0][0][0], &d3ccoil_.curntw[0], 
						  &X, &Y, &Z, &bx, &by, &bz);
B_X += bx;
B_Y += by;
B_Z += bz;

// I-coil perturbation field
bx = 0;	by = 0;	bz = 0;
if(PAR.useIcoil==1) polygonb_(&nIloops, &nIsegs, &nIloops, &nicsegs[0], &kuseI[0],
						  &d3iloops_.xis[0][0][0], &d3iloops_.divs[0][0][0], &d3icoil_.curntIc[0], 
						  &X, &Y, &Z, &bx, &by, &bz);
B_X += bx;
B_Y += by;
B_Z += bz;

// Field of any current filament
bx = 0;	by = 0;	bz = 0;
if(PAR.useFilament>0) get_filament_field(R,phi,Z,field,bx,by,bz,EQD);

B_X += bx;
B_Y += by;
B_Z += bz;

// Transform B_perturbation = (B_X, B_Y, B_Z) to cylindrical coordinates and add
B_R += B_X*cosp + B_Y*sinp;
B_phi += -B_X*sinp + B_Y*cosp;

if(PAR.response_field == -10)
{
	bx = 0;	bz = 0;
	FISLD.get_B(R,phi,Z,bx,bz,EQD);
	B_R += bx;
	B_Z += bz;
}

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

eqpol_.dsbp = 0.0;		// no shift !!!
eqpol_.dthbp = 0.0;		// no tilt !!!
eqpol_.alfsbp = 85.0*consts_.dtr;
eqpol_.alfthbp = 110.0*consts_.dtr;
eqpol_.ipdir = 1;		// positive Ip direction

d3ccoil_.scaleC = 1.0;
d3ccoil_.addanglC = 0.0;

d3icoil_.addanglIU = 0.0;
d3icoil_.addanglIL = 0.0;
d3icoil_.scaleIU = 1.0;
d3icoil_.scaleIL = 1.0;

// Prepare loading M3D-C1
if(PAR.response_field >= 0) chk = M3D.read_m3dc1sup(supPath);
else chk = 0;

// Read diiidsub.in file, if coils or M3D-C1 are on
if(PAR.useFcoil == 1 || PAR.useCcoil == 1 || PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	in.open(supPath + "diiidsup.in");
	if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open diiidsup.in file " << endl; EXIT;}

	for(i=1;i<=4;i++) in >> line;	// Skip 4 lines
	for(i=0;i<nFc;i++) in >> d3pfer_.fcur[i];		// Read F-coil currents

	in >> line;	// Skip line
	for(i=0;i<nCloops;i++) in >> d3ccoil_.curntC[i];		// Read C-coil currents

	in >> line;	// Skip line
	for(i=0;i<nIloops;i++) in >> d3icoil_.curntIc[i];		// Read I-coil currents

	in.close();	// close file
	in.clear();	// reset ifstream for next use
}

// Read C1.h5 file
if(PAR.response_field >= 0)
{
	if(chk == -1) M3D.scale_from_coils(d3icoil_.curntIc, nIloops, nIloops);	// no m3dc1sup.in file found -> scale from diiidsup.in file
	M3D.load(PAR, mpi_rank);
}
else
{
	if(mpi_rank < 1) cout << "Using g-file!" << endl;
	ofs2 << "Using g-file!" << endl;
}

// Prepare SIESTA
#ifdef USE_SIESTA
	if(PAR.response_field == -2)
	{
		if(mpi_rank < 1) cout << "Read SIESTA file" << endl;
		ofs2 << "Read SIESTA file" << endl;
		SIES.read("siesta.dat");
	}
#endif

// Prepare XFIELD
#ifdef USE_XFIELD
	if(PAR.response_field == -3)
	{
		if(mpi_rank < 1) cout << "Read XFIELD file" << endl;
		ofs2 << "Read XFIELD file" << endl;
		XPND.read("xpand.dat");
		if(mpi_rank < 1) cout << "NR = " << XPND.NR << "\t Nphi = " << XPND.Np-1 << "\t NZ = " << XPND.NZ << endl;
	}
#endif

if(mpi_rank < 1) cout << "F-coil: " << PAR.useFcoil << "\t" << "C-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << endl << endl;
ofs2 << "F-coil: " << PAR.useFcoil << "\t" << "C-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << endl << endl;

// Write I-coil currents to log files (Check if corretly read in)
ofs2 << "I-coil currents:" << endl;
for(i=0;i<nIloops/2;i++) ofs2 << d3icoil_.curntIc[i] << "\t";
ofs2 << endl;
for(i=nIloops/2;i<nIloops;i++) ofs2 << d3icoil_.curntIc[i] << "\t";
ofs2 << endl;

// Set F-coil geometry
if(PAR.useFcoil==1)
{
	for(i=0;i<nFlps;i++) kuseF[i]=0;	// kuseF is set inside the subroutine !?!
	d3pfgeom_(&kuseF[0]);	
}

// Set C-coil geometry
if(PAR.useCcoil==1)
{
	for(i=0;i<nCloops;i++) {if(d3ccoil_.curntC[i] != 0) kuseC[i]=1; else kuseC[i]=0;}
	for(i=0;i<nCloops;i++) nccsegs[i] = nCsegs;
	d3cgeom_(&kuseC[0]);
}

// Set I-coil geometry
if(PAR.useIcoil==1)
{
	for(i=0;i<nIloops;i++) {if(d3icoil_.curntIc[i] != 0) kuseI[i]=1; else kuseI[i]=0;}
	for(i=0;i<nIloops;i++) nicsegs[i] = nIsegs;
	d3igeom_(&kuseI[0]);
}

// Prepare filaments
if(PAR.useFilament>0)
{
	if(mpi_rank < 1) cout << "Interpolated filament field is used" << endl;
	ofs2 << "Interpolated filament field is used" << endl;
	in.open("filament_all.in");
	if(in.fail()==1)
	{
		if(mpi_rank == 1) cout << "Unable to open filament_all.in file. Please run fi_prepare." << endl; 
		EXIT;
	}
	else	// Read field on grid from file
	{
		// Set field size
		field.resize(Range(1,3),Range(0,359),Range(0,EQD.NR+1),Range(0,EQD.NZ+1));

		// Skip 3 lines
		in >> line;	
		if(mpi_rank < 1) cout << line.mid(3) << endl;
		ofs2 << line.mid(3) << endl;
		in >> line;	
		if(mpi_rank < 1) cout << line.mid(3) << endl;
		ofs2 << line.mid(3) << endl;
		in >> line;	

		// Read data
		for(int k=0;k<360;k++)
		{
			for(i=0;i<=EQD.NR+1;i++)
			{
				for(int j=0;j<=EQD.NZ+1;j++)
				{
					in >> field(1,k,i,j);
					in >> field(2,k,i,j);
					in >> field(3,k,i,j);
				}
			}
		}
		in.close();
	}
	in.clear();
	if(mpi_rank < 1) cout << endl;
	ofs2 << endl;
}

// Prepare fake Islands
if(PAR.response_field == -10)
{
	if(mpi_rank < 1) cout << "Read Fake Islands file" << endl;
	ofs2 << "Read Fake Islands file" << endl;
	FISLD.read("fakeIslands.in");
	FISLD.get_surfaces(EQD);
	if(mpi_rank < 1)
	{
		cout << "Amplitude: " << FISLD.A << endl;
		cout << "m: " << FISLD.m << endl;
		cout << "n: " << FISLD.n << endl;
		cout << "Phase: " << FISLD.delta << endl;
		cout << "Location: " << FISLD.psi0 << endl;
	}
}
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
Array<double,1> R,Z,S;	// Curve
int idx;
double x,Smax;

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
case 0:	// 19.59cm (same length as inner target) vertical wall above inner target, t = 0 -> -1, t=0 <=> P1 at inner target
	R1=1.016;	Z1=-1.223;
	R2=1.016;	Z2=-1.0271;
	break;
case 1:	// inner target plate
	R1=1.016;	Z1=-1.223;
	R2=1.153;	Z2=-1.363;
	break;
case 2:	// outer target plate
	R1=1.153;	Z1=-1.363;
	R2=1.372;	Z2=-1.363;	// normally R2=1.42, but that is inside the pump
	break;
case 3:	// 21.9cm (same length as outer target) horizontal shelf above pump to outer target
	R1=1.372;	Z1=-1.25;
	R2=1.591;	Z2=-1.25;
	break;
case 4:	// SAS divertor at upper outer divertor;  here t is dimensionless length along the wall;  t = 1 is same as Smax = 1.01693189 m; t = 0 is at the upper pump exit
	N = 36;
	R.resize(Range(1,N));
	Z.resize(Range(1,N));
	S.resize(Range(1,N));

	R = 1.372  ,  1.37167,  1.37003,  1.36688,  1.36719,  1.37178,
	        1.37224,  1.38662,  1.38708,  1.40382,  1.41127,  1.41857,
	        1.421  ,  1.48663,  1.4973 ,  1.49762,  1.49745,  1.49275,
	        1.4926 ,  1.49261,  1.49279,  1.4934 ,  1.4947 ,  1.49622,
	        1.47981,  1.48082,  1.48149,  1.48646,  1.49095,  1.50305,
	        1.59697,  1.6255 ,  1.63752,  1.647  ,  1.785  ,  2.07;
	Z = 1.31   ,  1.29238,  1.28268,  1.25644,  1.22955,  1.19576,
	        1.19402,  1.16487,  1.16421,  1.15696,  1.1573 ,  1.16132,
	        1.164  ,  1.2405 ,  1.23458,  1.23428,  1.23174,  1.2133 ,
	        1.21061,  1.20486,  1.20214,  1.19642,  1.18511,  1.1607 ,
	        1.12426,  1.12256,  1.12138,  1.11692,  1.11439,  1.11244,
	        1.09489,  1.0853 ,  1.07988,  1.077  ,  1.077  ,  1.04;
	Smax = 1.01693189;

	if(tmin < 0 || tmax > 1) ofs2 << "start_on_target: Warning, t out of range" << endl;
	S(1) = 0;
	idx = 1;
	for(int i=2;i<=N;i++)
	{
		S(i) = S(i-1) + sqrt((R(i)-R(i-1))*(R(i)-R(i-1)) + (Z(i)-Z(i-1))*(Z(i)-Z(i-1)));	//length of curve in m
		if(S(i) < Smax*t) idx = i;
		else break;
	}
	p1(1) = R(idx);		p1(2) = Z(idx);
	p2(1) = R(idx+1);		p2(2) = Z(idx+1);
	d = p2 - p1;
	x = (Smax*t - S(idx))/sqrt(d(1)*d(1)+d(2)*d(2));	// rescale t in m (like S); x is dimensionless in [0,1]
	p = p1 + x*d;
	break;
default:
	ofs2 << "No target specified" << endl;
	EXIT;
	break;
}

if(target <= 3)
{
	p1(1) = R1;	 p1(2) = Z1;
	p2(1) = R2;	 p2(2) = Z2;
	d = p2 - p1;
	if(target == 0) d(2) *= -1;

	// Coordinates
	p = p1 + t*d;
}

FLT.R = p(1);
FLT.Z = p(2);
FLT.phi = (phimin + dphi*i_phi)*rTOd;	// phi in deg
FLT.get_psi(p(1),p(2),FLT.psi);

FLT.Lc = 0;
FLT.psimin = 10;

if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}
return t;
}

#endif // D3D_M3DC1_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

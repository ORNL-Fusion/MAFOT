// Header-File for the DIII-D Drift Programs modified for JET
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						20.2.08
// M.Rack						18.02.2011

// uses interpolated filament fields 

// Include
//--------
#include <mystructs.hxx>

// -------------------- typedef -------------------------------------------------------------
typedef struct {string name; double wert;} parstruct;

// --------------- Prototypes ----------------------------------------------------------------
void setMessageCoils(ostream &out);
void setMessageTarget(ostream &out);
LA_STRING setTargetType(void);

void readiodata(char* name, vector<double>& vec);
void writeiodata(ofstream& out, vector<LA_STRING>& var, parstruct * pv, int psize, char* name=0);

int mapit(double& theta, double& r, double& phistart, double& psi, int itt, int MapDirection=1, 
		  int get_parameter=0, double* parameter=0);
int mapit(double& theta, double& r, double& phistart, double& psi, int itt, int MapDirection, 
		  double eps, double& good);

int mapstep(Array<double,1>& y, double& theta, double& r, double& phi, double& psi, int MapDirection=1, 
		  int get_parameter=0, double* parameter=0);
int mapstep(Array<double,1>& y, double& theta, double& r, double& phi, double& psi, int MapDirection, 
			double eps, double& good);

int connect(Array<double,1>& xa, double phistart, int itt, int MapDirection, double& ntor, double& length, double& psimin);

void getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi);
void dgls(double x, Array<double,1> y, Array<double,1>& dydx);

int rkint(int nvar, int nstep, double dx, Array<double,1>& y, double& x,
		  int get_parameter=0, double* parameter=0);

void prep_perturbation(void);

void get_filament_field(double R, double phi, double Z, Array<double,4>& field, double& bx, double& by, double& bz);
void prepare_filament_field(Array<double,4>& field);
double filament_fields_ex(double x, double y, double z, Array<double,2>& data, int N, double I, double dIdt,
					   double& B_x, double& B_y, double& B_z, double& dAdt, double cosp, double sinp);

void set(int i, int N, double xmin, double xmax, double ymin, double ymax, double& x, double& y, int N_x=1);
void create(long& idum, double xmin, double xmax, double ymin, double ymax, double& x, double& y);
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 double& r, double& theta, double& phi);

int setMapDirection(void);

void rungekutta4(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout);

int odeint(int nvar, Array<double,1>& ystart, double xstart, double& xend, 
			double eps, double h1, double hmin, int& nok, int& nbad);
void rkqs(int nvar, Array<double,1>& y, Array<double,1>& dydx, double& x, 
		  double htry, double eps, Array<double,1>& yscal, double& hdid, double& hnext);
void rungekutta5(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout,
				 Array<double,1>& yerr);
//void rungekutta5(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout);

void bcuderiv_square(Array<double,2>& y, int j, int k, double d1, double d2, 
					 Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12);
void bcuint_square(Array<double,1>& Ra, Array<double,1>& Za, double dR, double dZ, Array<double,2>& field,
			double R, double Z, double& y, double& y1, double& y2);

void Bfield_EFCC(double &Bx, double &By, double &Bz, double x, double y, double z);
void Bfield_EF(double &Bx, double &By, double &Bz, double x, double y, double z);

bool outofBndy(double x, double y, int simpleBndy);
bool outofRealBndy(double x, double y);

// --------------------Global Variables------------------------------------------------------
//have to be known during integration for perturbations, set in: prep_perturbation()
const int nEFCC = 4;	// number of Coils
Array<double,1> curntEFCC(nEFCC);	// Array for EFCC current
Array<double,1> curntEF(nEFCC);		// Array for Error Field current
Array<vektor,1> EFCC1(14);
Array<vektor,1> EFCC3(16);
Array<vektor,1> EFCC5(14);
Array<vektor,1> EFCC7(14);
Array<double,3> filament_data(Range(1,1),Range(0,1),Range(1,5));	// default size
Array<double,4> field;	// default constructed

// ------------------- extern Parameters ----------------------------------------------------
extern EFIT EQD;
extern double GAMMA;
extern double eps0;
extern double Ix;

// -------------------- Switches ------------------------------------------------------------
int simpleBndy = 0;		// 0: real boundary 	1: rectangle boundary

int useEF = 0;			// 0: no	1: yes
int useEFCC = 0;		// 0: no	1: yes
int useFilament = 0;		// 0: no	>= 1: Number of Filaments to be included 

int which_target_plate;	// 0: inner target	1: outer target
int create_flag;	// 0: fixed grid	1: random numbers	2: Start on target

int sigma;	// 1: co-passing particles		-1: count-passing particles		0: field lines only
int Zq;		// Charge number: 1: ions are calculated	-1: electrons are calculated

// -------------- global Parameters --------------------------------------------------------
const double c=299792458;			// speed of light in m/s
const double e=1.602176462*1e-19;	// elementary charge in C
const double me=9.10938188*1e-31;	// Electron mass in kg
const double mp=1.67262158*1e-27;	// Proton mass in kg
const double E0e=me*c*c/e/1000.0;	// Rest energy of Electrons in [keV] 
const double E0p=mp*c*c/e/1000.0;	// Rest energy of Protons in [keV] 

const double pi = LA_PI;
const double pi2 = 2*pi;
const double rTOd = 180.0/pi;	// rad to deg

const int Massnumber=1;				// Mass number

const int nvar = 2;	// Number of Variables
const int ilt = 360;	// Steps till Output
const double dpinit = 1.0;	// step size of phi in [deg]
double bndy[4] = {1.8060, 3.8910, -1.7460, 2.0180};	// simple Boundary (x_min, x_max, y_min, y_max)

// -------------- MPI Parameters -----------------------------------------------------------
int mpi_rank = 0;
int mpi_size = 0;

#ifdef USE_MPI
#define EXIT MPI::COMM_WORLD.Abort(0)
#else
#define EXIT exit(0)
#endif

// ------------------ log file -----------------------------------------------
#ifndef program_name		// checks if program_name is defined
#define program_name "default"	// if not, set program_name to default
#endif						// end

ofstream ofs2;

// ---------------------- IO -------------------------------------------------------
// ---------------------------------------------------------------------------------

// ----------------- some machine specific output (messages) ---------------------------------
void setMessageCoils(ostream &out){ out << "EF: " << useEF << "\t" << "EFCC: " << useEFCC << endl << endl; }
void setMessageTarget(ostream &out){ out << "Target (0=inner, 1=outer): " << which_target_plate << endl; }
LA_STRING setTargetType(void)
{
  switch(which_target_plate)
  {
    case 0:
      return "_inn";
      break;
    case 1:
      return "_out";
      break;
    default:
      return "";
      break;
  }
}

// ----------------- readiodata ------------------------------------------------------------
void readiodata(char* name, vector<double>& vec)
{
LA_STRING input;	// !!! LA_STRING reads entire line, string reads only one word !!!

// Get ShotNr and ShotTime from Parameterfile
ifstream in;
in.open(name);
if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open file " << name << endl; EXIT;}
in >> input;	// Skip first line
in >> input;	// second line gives shot number and time
EQD.Shot = input.mid(9,6);	// 6 characters starting at index 9 of input string
EQD.Time = input.mid(22,5); 	// 5 characters starting at index 22 of input string

// Get Path to g-File from Parameterfile (optional), Linux only!!!
in >> input;	
if(input[1] == '#') 
{
	input = input.mid(input.indexOf('/'));	// all characters of input string starting from index of first '/' in string
	if(input.right(1) != '/') input = input.left(input.length()-1);			// last char in string can be '\r' (Carriage return) --> Error!;  this removes last char if necessary 
	if(input.indexOf(' ') > 1) input = input.left(input.indexOf(' ')-1);		// blanks or comments that follow path are removed from string 
	EQD.Path = input;
}
in.close();

// Get Parameters
readparfile(name,vec);

// Set switches
which_target_plate = int(vec[11]);
create_flag = int(vec[12]);

if(vec[13]>1) 
{
	if(mpi_rank < 1) cout << "Coil flags are undefined in file " << name << endl;
	ofs2 << "Coil flags are undefined in file " << name << endl;
}
else
{
	useEF = int(vec[13]);
	useEFCC = int(vec[14]);
	//useIcoil = int(vec[15]);	//Not used in JET
}

if(vec.size()>=22)
{
	sigma = int(vec[16]);
	Zq = int(vec[17]);
}
else {if(mpi_rank < 1) cout << "Fail to read particle parameters" << endl; EXIT;}

if(vec.size()>=23)
{
	useFilament = int(vec[20]);
}
else {if(mpi_rank < 1) cout << "No current filaments included" << endl;}
}

// ------------------- writeiodata ----------------------------------------------------------
void writeiodata(ofstream& out, vector<LA_STRING>& var, parstruct * pv, int psize, char* name)
{
int i;
out << "# " << program_name << endl;
out << "#-------------------------------------------------" << endl;
out << "### Parameterfile: " << name << endl;
out << "# Shot: " << EQD.Shot << endl;
out << "# Time: " << EQD.Time << endl;
out << "#-------------------------------------------------" << endl;
out << "### Switches:" << endl;
out << "# EF active (0=no, 1=yes): " << useEF << endl;
out << "# EFCC active (0=no, 1=yes): " << useEFCC << endl;
out << "# No. of current filaments (0=none): " << useFilament << endl;
out << "# Target (0=inner, 1=outer): " << which_target_plate << endl;
out << "# Create Points (0=grid, 1=random, 2=target): " << create_flag << endl;
out << "# Direction of particles (1=co-pass, -1=count-pass, 0=field lines): " << sigma << endl;
out << "# Charge number of particles (=-1:electrons, >=1:ions): " << Zq << endl;
out << "#-------------------------------------------------" << endl;
out << "### Global Parameters:" << endl;
out << "# Steps till Output (ilt): " << ilt << endl;
out << "# Step size (dpinit): " << dpinit << endl;
if (simpleBndy)
{
  out << "# Boundary Rmin: " << bndy[0] << endl;
  out << "# Boundary Rmax: " << bndy[1] << endl;
  out << "# Boundary Zmin: " << bndy[2] << endl;
  out << "# Boundary Zmax: " << bndy[3] << endl;
} else {
  out << "# real Boundary is used" << endl;
}
out << "# Magnetic Axis: R0: " << EQD.RmAxis << endl;
out << "# Magnetic Axis: Z0: " << EQD.ZmAxis << endl;
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

// -----------------------------------------------------------------------------------------------
// ----------------- ende IO ---------------------------------------------------------------------

//---------------- mapit ----------------------------------------------------------------------------------------------------
int mapit(double& theta, double& r, double& phistart, double& psi, int itt, int MapDirection, 
		  int get_parameter, double parameter[])
{
int i,chk;
double phi;
double L = 0;
Array<double,1> y(nvar); // Array to set initial conditions

y(0) = r*cos(theta) + EQD.RmAxis;	// R
y(1) = r*sin(theta) + EQD.ZmAxis;	// Z
phi = phistart;

// Integrate
for(i=1;i<=itt;i++)
{
	chk = mapstep(y,theta,r,phi,psi,MapDirection,get_parameter,parameter);
	if(chk<0) {ofs2 << "mapit: wall hit" << endl; break;}	// particle has left system

	if(fabs(phi-MapDirection*i*dpinit*ilt-phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(phi-MapDirection*i*dpinit*ilt-phistart) << endl;
	phi=MapDirection*i*dpinit*ilt+phistart;

	if(get_parameter!=0) L += parameter[0];
}
phistart = phi;
if(get_parameter!=0) parameter[0] = L;
return chk;
}

//---------------- mapit ----------------------------------------------------------------------------------------------------
// same as mapit above, but uses adaptive step size integrator, see mapstep
// eps is the desired accuracy
// goodsum is number of successfull steps relative to total number of steps, summed over all iterations
int mapit(double& theta, double& r, double& phistart, double& psi, int itt, int MapDirection, 
		  double eps, double& goodsum)
{
int i,chk;
double phi;
double good;
Array<double,1> y(nvar); // Array to set initial conditions

y(0) = r*cos(theta) + EQD.RmAxis;	// R
y(1) = r*sin(theta) + EQD.ZmAxis;	// Z
phi = phistart;
goodsum = 0; 

// Integrate
for(i=1;i<=itt;i++)
{
	chk = mapstep(y,theta,r,phi,psi,MapDirection,eps,good);
	if(chk<0) {ofs2 << "mapit: wall hit" << endl; break;}	// particle has left system
	goodsum += good;

	if(fabs(phi-MapDirection*i*dpinit*ilt-phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(phi-MapDirection*i*dpinit*ilt-phistart) << endl;
	phi=MapDirection*i*dpinit*ilt+phistart;
}
goodsum /= double(i-1);
phistart = phi;
return chk;
}

//---------------- mapstep --------------------------------------------------------------------------------------------------
int mapstep(Array<double,1>& y, double& theta, double& r, double& phi, double& psi, int MapDirection, 
		  int get_parameter, double parameter[])
{
int chk;
double dummy;
double Rm,Zm;
double dphi = MapDirection*dpinit/rTOd;
double phi_rad = phi/rTOd;	// phi in rad

// integrate one full toroidal turn
chk = rkint(nvar,ilt,dphi,y,phi_rad,get_parameter,parameter);	
if(chk<0) return -1;	// particle has left system

// Get poloidal coordinates
Rm = y(0) - EQD.RmAxis; // R relative to magnetic axis
Zm = y(1) - EQD.ZmAxis; // Z relative to magnetic axis
r = sqrt(Rm*Rm + Zm*Zm);
theta = atan(Zm/Rm);
if(Rm<0) theta += pi;
if(Rm>=0 && Zm<0) theta += pi2;

// phi back in deg
phi = phi_rad*rTOd;

// Get normalized flux
EQD.get_psi(y(0),y(1),psi,dummy,dummy);

if(get_parameter!=0) 
{
	if(psi<parameter[1]) parameter[1] = psi;
}

return 0;
}

//---------------- mapstep --------------------------------------------------------------------------------------------------
// same as mapstep above, but uses adaptive step size integrator: odeint
// eps is desired accuracy
// good is number of successfull steps relative to total number of steps
int mapstep(Array<double,1>& y, double& theta, double& r, double& phi, double& psi, int MapDirection, 
			double eps, double& good)
{
int chk;
double dummy;
double Rm,Zm;
double dphi = MapDirection*dpinit/rTOd/10.0;	// Minimum step size
double phistart = phi/rTOd;	// phi in rad
double phiend = phistart + pi2*MapDirection;

int nok,nbad;

// integrate one full toroidal turn
chk = odeint(nvar,y,phistart,phiend,eps,100*dphi,dphi,nok,nbad);	
if(chk<0) return -1;	// particle has left system

// Get poloidal coordinates
Rm = y(0) - EQD.RmAxis; // R relative to magnetic axis
Zm = y(1) - EQD.ZmAxis; // Z relative to magnetic axis
r = sqrt(Rm*Rm + Zm*Zm);
theta = atan(Zm/Rm);
if(Rm<0) theta += pi;
if(Rm>=0 && Zm<0) theta += pi2;

// phi back in deg
phi = phiend*rTOd;

// Get normalized flux
EQD.get_psi(y(0),y(1),psi,dummy,dummy);

//ofs2 << "Good steps: " << nok << "\t" << "Bad steps: " << nbad << endl;
good = nok/double(nok+nbad);
return 0;
}

//------------------ connect --------------------------------------------------------------------------
// Integration goes in respecive direction, depending on MapDirection
// MapDirection=0 means both directions are calculated and results are added
// phistart has to be in deg!
int connect(Array<double,1>& xa, double phistart, int itt, int MapDirection, double& ntor, double& length, double& psimin)
{
int chk;
double r,theta,phi,psi;
double parameter[2] = {0, 10};

theta = modulo2pi(xa(1));	
r = xa(2);
phi = phistart;
length = 0;

// positive phi direction
if(MapDirection >= 0)
{
	chk = mapit(theta,r,phi,psi,itt,1,1,parameter);
	ntor = fabs(phi-phistart)/360.0;
	length += parameter[0];
	psimin = parameter[1];
}
else
{
	ntor = 0;
	length = 0;
	psimin = 10;
}

// negative phi direction
theta = modulo2pi(xa(1));	
r = xa(2);
phi = phistart;
parameter[0] = 0;  parameter[1] = 10;

if(MapDirection <= 0)
{
	chk = mapit(theta,r,phi,psi,itt,-1,1,parameter);
	ntor += fabs(phi-phistart)/360.0;
	length += parameter[0];
	if(parameter[1] < psimin) psimin = parameter[1];
}

return 0;
}

//---------------- getBfield ----------------------------------------------------------------------------------------------
void getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi)
{
int chk;
double psi,dpsidr,dpsidz;
double F;
double X,Y,bx,by,bz;
double B_X,B_Y;
double sinp,cosp;

sinp = sin(phi);
cosp = cos(phi);

X = R*cosp;
Y = R*sinp;

// get normalized poloidal Flux psi (should be chi in formulas!)
chk = EQD.get_psi(R,Z,psi,dpsidr,dpsidz);
if(chk==-1) {ofs2 << "Point is outside of EFIT grid" << endl; B_R=0; B_Z=0; B_phi=1; return;}	// integration of this point terminates 

// Equilibrium field
F = EQD.get_Fpol(psi);
B_R = dpsidz/R;
B_phi = F/R;
//B_phi = EQD.Bt0*EQD.R0/R;
B_Z = -dpsidr/R;

B_X = 0;	B_Y = 0;

// EF perturbation field
bx = 0;	by = 0;	bz = 0;
if(useEF==1) Bfield_EF(bx, by, bz, X, Y, Z);
B_X += bx;
B_Y += by;
B_Z += bz;

// EFCC-coil perturbation field
bx = 0;	by = 0;	bz = 0;
if(useEFCC==1) Bfield_EFCC(bx, by, bz, X, Y, Z);
B_X += bx;
B_Y += by;
B_Z += bz;

// Field of any current filament
bx = 0;	by = 0;	bz = 0;
if(useFilament>0) get_filament_field(R,phi,Z,field,bx,by,bz);
B_X += bx;
B_Y += by;
B_Z += bz;

//double dummy,I;
//Range all = Range::all();
//Array<double,2> slice;
//for(int i=1;i<=useFilament;i++)	// no filament fields, if useFilament == 0
//{
//	bx = 0;	by = 0;	bz = 0;
//	//I = get_current(i,0,dummy);	// time is fixed at t = 0; dummy is used for dIdt
//	I = filament_data(i,0,2);	// current is stored in here, see prep_perturbation
//	// filament_data(i,0,1) gives number of rows of data; index 0 is not used in filament_fields_ex
//	slice.reference(filament_data(i,all,all));
//	dummy = filament_fields_ex(X,Y,Z,slice,int(filament_data(i,0,1)),I,0,	// dIdt is set to 0 and fixed <- constant current
//							   bx,by,bz,dummy,cosp,sinp);			// dummy is used instead of A_phi and dAdt
//	B_X += bx;
//	B_Y += by;
//	B_Z += bz;
//}

// Transform B_perturbation = (B_X, B_Y, B_Z) to cylindrical coordinates and add
B_R += B_X*cosp + B_Y*sinp;
B_phi += -B_X*sinp + B_Y*cosp;
}

//---------------- dgls ---------------------------------------------------------------------------------------------------
// Type the differential equations (dgls) as they are written, while x is the independent variable
// Here: x = phi, y(0) = R, y(1) = Z, dydx(0) = dR/dphi, dydx(1) = dZ/dphi
void dgls(double x, Array<double,1> y, Array<double,1>& dydx)
{
double B_R,B_Z,B_phi;
double S;
getBfield(y(0),y(1),x,B_R,B_Z,B_phi);

dydx(0) = y(0)*B_R/B_phi;
dydx(1) = y(0)*B_Z/B_phi;

if(sigma !=0)
{
	S = eps0*(GAMMA*GAMMA-1)-2*EQD.R0*Ix/y(0);
	if(S<0) {ofs2 << "dgls: Sqrt argument negative => Abort" << endl; EXIT;}	// Error -> Abort program
	S = sqrt(S);
	dydx(1) += -sigma/double(Zq)*(y(0)*S + EQD.R0*Ix/S);	// sign corrected: sigma = +1 is indeed co-passing
}
}

//----------------- rkint -------------------------------------------------------------------------------------------------
//Starting from initial values y[0..nvar-1] known at x=x1 use Runge-Kutta
//to advance nstep equal increments to x2=x1+nstep*dx. The user-supplied routine dgls(x,v,dvdx)
//evaluates derivatives. Results after nstep Steps are stored in y[0..nvar-1]
//if flag 'get_parameter' is set: additional parameters are returned
//parameter[0]: length of trajectory, parameter[1]: minimal psi during integration
int rkint(int nvar, int nstep, double dx, Array<double,1>& y, double& x, 
		  int get_parameter, double parameter[])
{
int k;
double x1 = x;	//Store first value (helps reduce Error in x)
double psi;
double L = 0;
double dummy;
Array<double,1> yout(nvar),dydx(nvar);

//Take nstep steps
for (k=1;k<=nstep;k++) 
{ 
	dgls(x,y,dydx);
	rungekutta4(y,dydx,nvar,x,dx,yout);
	x = x1 + k*dx; // Better than x+=dx

	// Integration terminates outside of boundary box
	if(outofBndy(yout(0),yout(1),simpleBndy)) return -1;

	if(get_parameter!=0) 
	{
		L += sqrt((yout(0)-y(0))*(yout(0)-y(0)) + (yout(1)-y(1))*(yout(1)-y(1)) + 0.25*(yout(0)+y(0))*(yout(0)+y(0))*dx*dx);
		EQD.get_psi(y(0),y(1),psi,dummy,dummy);
		if(psi<parameter[1]) parameter[1] = psi;
	}

	y = yout;
} 
if(get_parameter!=0) parameter[0] = L;
return 0;
}

//---------- prep_perturbation --------------------------------------
void prep_perturbation(void)
{
int i;
LA_STRING line;	// entire line is read by ifstream


if(mpi_rank < 1) setMessageCoils(cout);
setMessageCoils(ofs2);

// Read jetsup.in file
ifstream in;
in.open("jetsup.in");
if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open jetsup.in file " << endl; EXIT;}

for(i=1;i<=4;i++) in >> line;	// Skip 4 lines
for(i=0;i<nEFCC;i++) in >> curntEF(i);		// Read EF currents 

in >> line;	// Skip line
for(i=0;i<nEFCC;i++) in >> curntEFCC(i);	// Read EFCC-coil currents

in.close();	// close file
in.clear();	// reset ifstream for next use

// Set EFCC geometry
if((useEF==1)|(useEFCC==1))
{
  //set EFCC1 geometry
  EFCC1( 0).set(-3.324, 6.170, 2.131);
  EFCC1( 1).set(-2.987, 5.560, 2.454);
  EFCC1( 2).set(-2.987, 5.560, 3.090);
  EFCC1( 3).set( 2.987, 5.560, 3.090);
  EFCC1( 4).set( 2.987, 5.560, 2.498);
  EFCC1( 5).set( 3.324, 6.170, 1.805);
  EFCC1( 6).set( 3.324, 6.170,-2.990);
  EFCC1( 7).set( 3.455, 5.313,-2.990);
  EFCC1( 8).set( 3.338, 4.424,-2.990);
  EFCC1( 9).set( 1.125, 5.341,-2.990);
  EFCC1(10).set(-1.125, 5.341,-2.990);
  EFCC1(11).set(-3.338, 4.424,-2.990);
  EFCC1(12).set(-3.455, 5.313,-2.990);
  EFCC1(13).set(-3.324, 6.170,-2.990);
  //set EFCC3 geometry
  EFCC3( 0).set(-6.170,-3.324, 2.131);
  EFCC3( 1).set(-5.561,-2.987, 2.453);
  EFCC3( 2).set(-5.561,-2.987, 2.835);
  EFCC3( 3).set(-5.660,-2.735, 3.087);
  EFCC3( 4).set(-5.660, 2.735, 3.087);
  EFCC3( 5).set(-5.561, 2.987, 2.835);
  EFCC3( 6).set(-5.561, 2.987, 2.453);
  EFCC3( 7).set(-6.170, 3.324, 2.131);
  EFCC3( 8).set(-6.170, 3.324,-2.992);
  EFCC3( 9).set(-5.120, 3.493,-2.992);
  EFCC3(10).set(-4.400, 3.398,-2.992);
  EFCC3(11).set(-5.342, 1.125,-2.992);
  EFCC3(12).set(-5.342,-1.125,-2.992);
  EFCC3(13).set(-4.400,-3.398,-2.992);
  EFCC3(14).set(-5.120,-3.493,-2.992);
  EFCC3(15).set(-6.170,-3.324,-2.992);
  //set EFCC5 geometry
  EFCC5( 0).set(-3.324,-6.170, 2.131);
  EFCC5( 1).set(-2.987,-5.560, 2.454);
  EFCC5( 2).set(-2.987,-5.560, 3.090);
  EFCC5( 3).set( 2.987,-5.560, 3.090);
  EFCC5( 4).set( 2.987,-5.560, 2.498);
  EFCC5( 5).set( 3.324,-6.170, 1.805);
  EFCC5( 6).set( 3.324,-6.170,-2.990);
  EFCC5( 7).set( 3.455,-5.313,-2.990);
  EFCC5( 8).set( 3.338,-4.424,-2.990);
  EFCC5( 9).set( 1.125,-5.341,-2.990);
  EFCC5(10).set(-1.125,-5.341,-2.990);
  EFCC5(11).set(-3.338,-4.424,-2.990);
  EFCC5(12).set(-3.455,-5.313,-2.990);
  EFCC5(13).set(-3.324,-6.170,-2.990);
  //set EFCC7 geometry
  EFCC7( 0).set( 6.170,-3.324, 2.131);
  EFCC7( 1).set( 5.560,-2.987, 2.454);
  EFCC7( 2).set( 5.560,-2.987, 3.090);
  EFCC7( 3).set( 5.660, 2.987, 3.090);
  EFCC7( 4).set( 5.660, 2.987, 2.498);
  EFCC7( 5).set( 6.170, 3.324, 1.805);
  EFCC7( 6).set( 6.170, 3.324,-2.990);
  EFCC7( 7).set( 5.313, 3.455,-2.990);
  EFCC7( 8).set( 4.424, 3.338,-2.990);
  EFCC7( 9).set( 5.341, 1.125,-2.990);
  EFCC7(10).set( 5.341,-1.125,-2.990);
  EFCC7(11).set( 4.424,-3.338,-2.990);
  EFCC7(12).set( 5.313,-3.455,-2.990);
  EFCC7(13).set( 6.170,-3.324,-2.990);
}

// Prepare filaments
if(useFilament>0)
{
	if(mpi_rank < 1) cout << "Interpolated filament field is used" << endl;
	ofs2 << "Interpolated filament field is used" << endl;
	in.open("filament_all.in");
	if(in.fail()==1)
	{
		if(mpi_size == 0) prepare_filament_field(field);	// Read filament files and get field on grid; writes filament_all.in file
		else {if(mpi_rank == 1) cout << "Unable to open filament_all.in file. Please run fi_prepare." << endl; EXIT;}
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
}

//------------------ get_filament_field --------------------------------------------------------
// determines sum of magnetic fields of all current filaments at point (R,phi,Z) by interpolation
// field contains total field on EFIT grid
void get_filament_field(double R, double phi, double Z, Array<double,4>& field, double& bx, double& by, double& bz)
{
int xyz;
double dummy;
Range all = Range::all();
Array<double,2> slice;
Array<double,1> y(Range(1,3));

// Get phi slice: phi in deg and rounded to integer
int k = int(phi*rTOd+0.5*sign(phi));	// has to be set between 0 and 359
while(k<0) k += 360;	// move k to a positive value
k = k%360;	// k now between 0 and 359

// Get field in carthesian coordinates
for(xyz=1;xyz<=3;xyz++)
{
	slice.reference(field(xyz,k,all,all));
	bcuint_square(EQD.R,EQD.Z,EQD.dR,EQD.dZ,slice,R,Z,y(xyz),dummy,dummy);
}
bx = y(1);
by = y(2);
bz = y(3);
}

//------------------ prepare_filament_field ----------------------------------------------------
// Reads filament files and calculates total field on EFIT grid
// phi is discretized in 1 deg. steps, according to integrator step size  -->  no interpolation in phi!!!!!
// Bicubic interpolation on R and Z.
void prepare_filament_field(Array<double,4>& field)
{
// Variables
int i,j,k;
ofstream out;

// firstIndex: 1 = B_X, 2 = B_Y, 3 = B_Z;   secondIndex: phi(0,...,359);   thirdIndex: R(0,1,...,NR,NR+1);   forthIndex: Z(0,1,...,NZ,NZ+1); 
// first and last index in R and Z are over the boundary to calculate derivatives at boundary
field.resize(Range(1,3),Range(0,359),Range(0,EQD.NR+1),Range(0,EQD.NZ+1));

// Read current-filament paths
//-----------------------------
ifstream in;
LA_STRING filamentfile,line;
string word;

Array<double,2> data;

if(useFilament>0) 
{
	//write field to file
	out.open("filament_all.in");
	out.precision(16);

	cout << "Number of current filaments used: " << useFilament << endl; 
	ofs2 << "Number of current filaments used: " << useFilament << endl; 
	out << "# Number of current filaments used: " << useFilament << endl; 
	cout << "Filament currents: ";
	ofs2 << "Filament currents: ";
	out << "# Filament currents: ";
}
for(i=1;i<=useFilament;i++)
{
	filamentfile = "filament" + LA_STRING(i) + ".in";
	readfile(filamentfile,5,data);

	// Store data in filament_data -> Adjust size if necessary 
	if(data.rows()>=filament_data.cols()) filament_data.resizeAndPreserve(useFilament,data.rows()+1,5);
	for(int j=1;j<=data.rows();j++) for(int k=1;k<=5;k++) filament_data(i,j,k) = data(j,k);
	filament_data(i,0,1) = data.rows();	// number of rows in data is stored at col-index 0 and depth-index 1
	data.free();	// data no longer used; prepare data for next use

	// read current from filament file -> store it in filament_data(i,0,2)
	in.open(filamentfile);	// location in the file needed !!!!!!!!!
	in >> line;
	in >> line;
	in >> word;	in >> word;
	in >> filament_data(i,0,2); 
	in.close();
	cout << i << ": " << filament_data(i,0,2) << "A" << "\t";
	ofs2 << i << ": " << filament_data(i,0,2) << "A" << "\t";
	out << i << ": " << filament_data(i,0,2) << "A" << "\t";
}
if(useFilament>0) 
{
	cout << endl;
	ofs2 << endl;
	out << endl << "#" << endl;
}

// Calculate total field of filaments on grid
//--------------------------------------------
double R,phi,X,Y,Z;
double bx,by,bz;
double B_X,B_Y,B_Z;
double sinp,cosp;
double dummy,I;
Range all = Range::all();
Array<double,2> slice;
Array<double,2> ddR(Range(1,EQD.NR),Range(1,EQD.NZ)),ddZ(Range(1,EQD.NR),Range(1,EQD.NZ)),d2dRdZ(Range(1,EQD.NR),Range(1,EQD.NZ));

// set grid
cout << "Calculate total field..." << endl;
ofs2 << "Calculate total field..." << endl;
for(k=0;k<360;k++)
{
	phi = double(k)/rTOd;
	sinp = sin(phi);
	cosp = cos(phi);

	for(i=0;i<=EQD.NR+1;i++)
	{
		if(i==0) R = EQD.R(1) - EQD.dR;
		if(i==EQD.NR+1) R = EQD.R(EQD.NR) + EQD.dR;
		if(i>0 && i<=EQD.NR) R = EQD.R(i);
		X = R*cosp;
		Y = R*sinp;

		for(j=0;j<=EQD.NZ+1;j++) 
		{
			if(j==0) Z = EQD.Z(1) - EQD.dZ;
			if(j==EQD.NZ+1) Z = EQD.Z(EQD.NZ) + EQD.dZ;
			if(j>0 && j<=EQD.NZ) Z = EQD.Z(j);
			B_X = 0;	B_Y = 0;	B_Z = 0;	

			// Sum up fields of any current filament
			for(int l=1;l<=useFilament;l++)	// no filament fields, if useFilament == 0
			{
				bx = 0;	by = 0;	bz = 0;
				//I = get_current(l,0,dummy);	// time is fixed at t = 0; dummy is used for dIdt
				I = filament_data(l,0,2);	// current is stored in here, see prep_perturbation
				// filament_data(l,0,1) gives number of rows of data; index 0 is not used in filament_fields_ex
				slice.reference(filament_data(l,all,all));
				dummy = filament_fields_ex(X,Y,Z,slice,int(filament_data(l,0,1)),I,0,	// dIdt is set to 0 and fixed <- constant current
										   bx,by,bz,dummy,cosp,sinp);			// dummy is used instead of A_phi and dAdt
				B_X += bx;
				B_Y += by;
				B_Z += bz;
			}

			// Store total field
			field(1,k,i,j) = B_X;
			field(2,k,i,j) = B_Y;
			field(3,k,i,j) = B_Z;

			out << B_X << "\t" << B_Y << "\t" << B_Z << endl;

		} // end for j
	} // end for i
	//cout << k << "\t" << flush;
} // end for k
//cout << endl;
cout << "done" << endl;
ofs2 << "done" << endl;
}

//------------- filament_fields_ex ----------------------------------------
// returns toroidal component A_phi of vector potential at position (x, y, z) at time t.
// data = data(1..N, 1..5) contains the pathway of the filament with columns X, Y, Z, R, phi
// Current I(t) in [A] and its time derivative dIdt(t) has to be provided. t is not needed!
// data usually has the form (R, Z, phi, psi) after reading from file -> has to be converted in main program
// before calling this function
// B_x, B_y and B_z are the components of the magnetic field.
// dAdt is the time derivative of A_phi
// Uses exact solution for vector potential, NO Taylor-series
// sinp = sin(phi) and cosp = cos(phi); phi is the toroidal angle of position, phi is not needed
double filament_fields_ex(double x, double y, double z, Array<double,2>& data, int N, double I, double dIdt,
					   double& B_x, double& B_y, double& B_z, double& dAdt, double cosp, double sinp)
{
int i,vor;
double dist,dist2;
double jx,jy,jz;
double dx1,dy1,dz1;
double dx2,dy2,dz2;
double d1,d2;
double dAx,dAy,dAz,A_part;
double A_x,A_y,A_phi;
double f;

B_x = 0;	B_y = 0;	B_z = 0;	
A_x = 0;	A_y = 0;

// calculate fields
for(i=1;i<N;i++)
{
	// Direction of current, without I
	jx = data(i+1,1) - data(i,1);
	jy = data(i+1,2) - data(i,2);
	jz = data(i+1,3) - data(i,3);
	dist = sqrt(jx*jx + jy*jy + jz*jz);
	jx /= dist;
	jy /= dist;
	jz /= dist;

	// Distance between position and filament
	dx1 = x - data(i,1);
	dy1 = y - data(i,2);
	dz1 = z - data(i,3);
	dx2 = x - data(i+1,1);
	dy2 = y - data(i+1,2);
	dz2 = z - data(i+1,3);

	dist = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);	
	dist2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

	// Set help variables
	d1 = dist - (jx*dx1 + jy*dy1 + jz*dz1);
	d2 = dist2 - (jx*dx2 + jy*dy2 + jz*dz2);
	vor = 1;

	// Help variables can get 0 <- not allowed --> Use another antiderivative (differs only on a constant -> same A)
	if(d1 <= 0 || d2 <= 0)
	{
		d1 = dist + (jx*dx1 + jy*dy1 + jz*dz1);
		d2 = dist2 + (jx*dx2 + jy*dy2 + jz*dz2);
		vor = -1;
	}

	// Vector potential, still without I !!!!
	A_part = vor*log(d2 / d1);
	A_x += A_part*jx;
	A_y += A_part*jy;
	//A_z += A_part*jz;

	// magnetic field
	dAx = vor*((dx2/dist2 - jx)/d2 - (dx1/dist - jx)/d1);
	dAy = vor*((dy2/dist2 - jy)/d2 - (dy1/dist - jy)/d1);
	dAz = vor*((dz2/dist2 - jz)/d2 - (dz1/dist - jz)/d1);

	B_x += jz*dAy - jy*dAz;
	B_y += jx*dAz - jz*dAx;
	B_z += jy*dAx - jx*dAy;
}
// transform A_x and A_y to A_phi
A_phi = A_y*cosp - A_x*sinp;

// Multiply with constants: mu0/4pi*I/L, (mu0/4pi = 1e-7)
f = 1e-7;
dAdt = f*dIdt*A_phi;

f *= I;
A_phi *= f;

B_x *= f;
B_y *= f;
B_z *= f;

return A_phi;
}

//-------------- set --------------------------------------------
// creates initial condition i: (x,y) on a grid with N points
// values taken from the intervals [xmin,xmax] and [ymin,ymax] respectively
// x and y are output variables
// if N_x is spezified, N=N_x*N_y has to be used -> Grid size is then N_x X N_y
// otherwise Grid size is about sqrt(N) X sqrt(N)
// y is varied first, x second
void set(int i, int N, double xmin, double xmax, double ymin, double ymax, double& x, double& y, int N_x)
{
double dx,dy;
int i_x=0,i_y=0;
int N_y;

if(N<=1) {dx=0; dy=0;}
else
{
	dx=(xmax-xmin)/(N-1);
	dy=(ymax-ymin)/(N-1);
}

if(dx==0) i_y=i-1;
if(dy==0) i_x=i-1;
if(dx!=0 && dy!=0) 
{
	if(N_x<=1) N_x=int(sqrt(double(N))+0.5);
	N_y=int(double(N)/double(N_x)+0.5);
	i_y=(i-1)%N_y;
	i_x=int(double(i-1)/double(N_y));
	dx=(xmax-xmin)/(N_x-1);
	dy=(ymax-ymin)/(N_y-1);
}

x=xmin+i_x*dx;
y=ymin+i_y*dy;
}

//-------------- create --------------------------------------------
// creates an initial condition (x,y) using random numbers
// values taken from the intervals [xmin,xmax] and [ymin,ymax] respectively
// x and y are output variables
void create(long& idum, double xmin, double xmax, double ymin, double ymax, double& x, double& y)
{
const double dx=xmax-xmin;
const double dy=ymax-ymin;
double z;

z=ran0(idum);
x=xmin+z*dx;

z=ran0(idum);
y=ymin+z*dy;
}

//---------------- start_on_target --------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t=[0 81.1966] for the inner target, t=[0 525.9103] for the outer target
// edges (koordinates and length) of target are explicitly defined here
// Position (R0,Z0) of magnetic axis is required
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 double& r, double& theta, double& phi)
{
int i_p=0;
int i_phi=0;
int N=Np*Nphi;
int Nedge,Nedge_tmp;	// number of edges (that describe the target)
int k=1;		// the k-th edge is selected
int target;
int IDtarget_start,IDtarget_end,IDdirection=1;	// IDs range (and direction) of target in EQD.wall
int IDtarget_max = EQD.wall.numElements()/2;
double dp,dphi,t,ts,suml;
Array<double,1> p(Range(1,2));

// Magnetic Axis
const double R0=EQD.RmAxis;
const double Z0=EQD.ZmAxis;

// Grid stepsizes and t
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
target=which_target_plate; 
ts = t;			// use ts because t needed for return value

if(t<0)
{
  ofs2 << "t have to be positive" << endl;
  EXIT;
}

// Coordinates of Target-Plate, IDs has to be modified for different machines
IDtarget_start = 97;
if (target==0) IDtarget_end = 80;	//  87.6966 cm inner target, R,Z in m
else IDtarget_end = 159;		// 525.9103 cm outer target, R,Z in m
Nedge = IDtarget_end-IDtarget_start;
if (Nedge < 0)
{
  Nedge *= -1;
  IDdirection = -1;
}

Array<double,1> R(Range(1,Nedge+1)),Z(Range(1,Nedge+1)),l(Range(1,Nedge));

if ((target>=0) && (target<=1))
{
  if (IDtarget_end > IDtarget_max)
  {
    Nedge_tmp = IDtarget_max - IDtarget_start;
    R(Range(1,Nedge_tmp+1)) = EQD.wall(Range(IDtarget_start*2-1,IDtarget_max*2-1,2))-R0;
    R(Range(Nedge_tmp+2,Nedge+1)) = EQD.wall(Range(3,(IDtarget_end-IDtarget_max+1)*2-1,2))-R0;
    Z(Range(1,Nedge_tmp+1)) = EQD.wall(Range(IDtarget_start*2,IDtarget_max*2,2))-Z0;
    Z(Range(Nedge_tmp+2,Nedge+1)) = EQD.wall(Range(4,(IDtarget_end-IDtarget_max+1)*2,2))-Z0;
  } else {
    R = EQD.wall(Range(IDtarget_start*2-1,IDtarget_end*2-1,IDdirection*2))-R0;
    Z = EQD.wall(Range(IDtarget_start*2,IDtarget_end*2,IDdirection*2))-Z0;
  }
} else {
  ofs2 << "No target specified" << endl;
  EXIT;
}  
    
// calculate legth of edges; l,suml in cm 
for(int k=1; k<=Nedge; k++)
  l(k) = sqrt((R(k+1)-R(k))*(R(k+1)-R(k)) + (Z(k+1)-Z(k))*(Z(k+1)-Z(k)))*100;
suml = sum(l);

while (ts > 0)
{
  if (k > Nedge)
  {
    ofs2 << "t is greater then the length of the selected target (" << suml << " [cm])" << endl;
    EXIT;
  }
  ts -= l(k);
  k++;
}

// Coordinates
if (ts==0)
{
  p(1) = R(k);
  p(2) = Z(k);
} else {
  ts /= l(k-1);
  p(1) = (R(k)-R(k-1))*ts + R(k);
  p(2) = (Z(k)-Z(k-1))*ts + Z(k);
}
// Get poloidal coordinates
theta=atan(p(2)/p(1));
if(p(1)<0) theta += pi;
if(p(1)>=0 && p(2)<0) theta += pi2;
r=sqrt(p(1)*p(1)+p(2)*p(2));
phi=phimin+dphi*i_phi;
return t;
}

//--------------- setMapDirection ---------------------------------------------------------------------------------------------
//used in dtfoot, sets the MapDirection depending on which target plate the fieldline is startet.
int setMapDirection(void)
{
  if(which_target_plate==1) return 1;
  else  return -1;
}

//--------------- rungekutta4 ---------------------------------------------------------------------------------------------
//Given values for the variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
//fourth-order Runge-Kutta method to advance the solution over an interval h and return the
//incremented variables as yout[0..n-1], which need not be a distinct array from y. The user
//supplies the routine dgl(x,y,dydx), which returns derivatives dydx at x.
void rungekutta4(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout)
{
int i;
double xh,hh,h6;
Array<double,1> dym(n),dyt(n),yt(n);

hh=h*0.5;
h6=h/6.0;
xh=x+hh;

//First step
for (i=0;i<n;i++) yt(i)=y(i)+hh*dydx(i); 

//Second step
dgls(xh,yt,dyt); 
for (i=0;i<n;i++) yt(i)=y(i)+hh*dyt(i);

//Third step
dgls(xh,yt,dym); 
for (i=0;i<n;i++) 
{
	yt(i)=y(i)+h*dym(i);
	dym(i) += dyt(i);
}

//Fourth step
dgls(x+h,yt,dyt); 
for (i=0;i<n;i++) yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0*dym(i)); //Accumulate increments with proper weights
 
}

//----------------- odeint --------------------------------------------------------------------------
//User storage for intermediate results. Preset kmax and dxsav in the calling program. If kmax !=
//0 results are stored at approximate intervals dxsav in the arrays xp[1..kount], yp[1..nvar]
//[1..kount], where kount is output by odeint. Defining declarations for these variables, with
//memory allocations xp[1..kmax] and yp[1..nvar][1..kmax] for the arrays, should be in
//the calling program.
//extern int kmax,kount;
//extern float *xp,**yp,dxsav;
//
//Runge-Kutta driver with adaptive stepsize control. Integrate starting values ystart[1..nvar]
//from xstart to xend with accuracy eps, storing intermediate results in global variables. h1 should
//be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can be zero). On
//output nok and nbad are the number of good and bad (but retried and fixed) steps taken, and
//ystart is replaced by values at the end of the integration interval. 
//xend is replaced with actual achived x at the end of the integration. dgls is the user-supplied
//routine for calculating the right-hand side derivative, while rkqs is the name of the stepper
//routine to be used.
int odeint(int nvar, Array<double,1>& ystart, double xstart, double& xend, 
			double eps, double h1, double hmin, int& nok, int& nbad)
{
int i;
double x,hnext,hdid,h;
static const double tiny = 1.0e-30;
static const double MAXSTP = 10000;
//double xsav;

Array<double,1> yscal(nvar),y(nvar),dydx(nvar);

x = xstart;
h = fabs(h1)*sign(xend-xstart);
nok = nbad = 0;
//kount = 0;

y = ystart;
//if (kmax > 0) xsav = x - dxsav*2.0; //Assures storage of first step.

//Take at most MAXSTP steps.
for (i=1;i<=MAXSTP;i++) 
{ 
	// Integration terminates outside of boundary box
	if(outofBndy(y(0),y(1),simpleBndy)) return -1;

	dgls(x,y,dydx);

	//Scaling used to monitor accuracy. This general-purpose choice can be modified if need be.
	//for (i=0;i<nvar;i++) yscal(i) = fabs(y(i)) + fabs(dydx(i)*h) + tiny;
	yscal = fabs(y) + fabs(dydx*h) + tiny;

	// Store intermediate results.
	//if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) 
	//{
	//	xp[++kount]=x; 
	//	for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
	//	xsav=x;
	//}

	if ((x+h-xend)*(x+h-xstart) > 0.0) h = xend - x; // If stepsize can overshoot, decrease.
	rkqs(nvar,y,dydx,x,h,eps,yscal,hdid,hnext);

	if (hdid == h) ++nok; 
	else ++nbad;

	//Are we done?
	if ((x-xend)*(xend-xstart) >= 0.0) 
	{ 
		ystart = y;
		//if (kmax > 0) 
		//{
		//	xp[++kount]=x; //Save final step.
		//	for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
		//}
		xend = x;
		return 0; //Normal exit.
	}
	if (fabs(hnext) <= hmin) 
	{
		ofs2 << "Step size too small in odeint" << endl;
		h = hmin;
	}
	else h = hnext;
}
ofs2 << "Too many steps in routine odeint" << endl;
return -1;
}

//------------- rkqs -----------------------------------------------------------------------------------
//Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and
//adjust stepsize. Input are the dependent variable vector y[1..nvar] and its derivative dydx[1..nvar]
//at the starting value of the independent variable x. Also input are the stepsize to be attempted
//htry, the required accuracy eps, and the vector yscal[1..n] against which the error is
//scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was
//actually accomplished, and hnext is the estimated next stepsize. dgls is the user-supplied
//routine that computes the right-hand side derivatives.
//max(a,b) Maximum of two values.
void rkqs(int nvar, Array<double,1>& y, Array<double,1>& dydx, double& x, 
		  double htry, double eps, Array<double,1>& yscal, double& hdid, double& hnext)
{
int i;
double errmax,h,htemp,xnew;
static const double SAFETY = 0.9;
static const double PGROW = -0.2;
static const double PSHRNK = -0.25;
static const double ERRCON = 1.89e-4;	//The value ERRCON equals (5/SAFETY) raised to the power (1/PGROW), see use below.

Array<double,1> yerr(nvar),ytemp(nvar);

h = htry; //Set stepsize to the initial trial value.

for (;;) 
{
	rungekutta5(y,dydx,nvar,x,h,ytemp,yerr); //Take a step.
	errmax = 0.0; //Evaluate accuracy.

	for (i=0;i<nvar;i++) errmax = max(errmax,fabs(yerr(i)/yscal(i)));
	errmax /= eps; //Scale relative to required tolerance.
	if (errmax <= 1.0) break; //Step succeeded. Compute size of next step.

	htemp = SAFETY*h*pow(errmax,PSHRNK);

	//Truncation error too large, reduce stepsize.
	h = (h >= 0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));	//No more than a factor of 10.
	xnew = x + h;
	if (xnew == x) {ofs2 << "stepsize underflow in rkqs" << endl; return;}
}

if (errmax > ERRCON) hnext = SAFETY*h*pow(errmax,PGROW);
else hnext = 5.0*h; // No more than a factor of 5 increase.
hdid = h;
x += hdid;
y = ytemp;
}

//--------------- rungekutta5 ---------------------------------------------------------------------------------------------
//Given values for the variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
//fifth-order Runge-Kutta method to advance the solution over an interval h and return the
//incremented variables as yout[0..n-1], which need not be a distinct array from y. The user
//supplies the routine dgls(x,y,dydx), which returns derivatives dydx at x.
//Returns an estimate of the error in yerr
void rungekutta5(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout,
				 Array<double,1>& yerr)
{
int i;
double xh;

static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;
static double c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0;
static double b21 = 0.2; 
static double b31 = 3.0/40.0, b32 = 9.0/40.0; 
static double b41 = 0.3, b42 = -0.9, b43 = 1.2;  
static double b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0, b54 = 35.0/27.0; 
static double b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0;
static double dc5 = -277.00/14336.0;
double dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0, dc4 = c4-13525.0/55296.0, dc6 = c6-0.25;

Array<double,1> yt(n),k2(n),k3(n),k4(n),k5(n),k6(n);

//First step
for (i=0;i<n;i++) yt(i)=y(i)+b21*h*dydx(i);

//Second step
xh=x+h*a2;
dgls(xh,yt,k2);
for (i=0;i<n;i++) yt(i)=y(i)+b31*h*dydx(i)+b32*h*k2(i);

//Third step
xh=x+h*a3;
dgls(xh,yt,k3);
for (i=0;i<n;i++) yt(i)=y(i)+b41*h*dydx(i)+b42*h*k2(i)+b43*h*k3(i);

//Forth step
xh=x+h*a4;
dgls(xh,yt,k4);
for (i=0;i<n;i++) yt(i)=y(i)+b51*h*dydx(i)+b52*h*k2(i)+b53*h*k3(i)+b54*h*k4(i);  

//Fifth step
xh=x+h*a5;
dgls(xh,yt,k5);
for (i=0;i<n;i++) yt(i)=y(i)+b61*h*dydx(i)+b62*h*k2(i)+b63*h*k3(i)+b64*h*k4(i)+b65*h*k5(i); 

//Last step
xh=x+h*a6;
dgls(xh,yt,k6);
for (i=0;i<n;i++) yout(i)=y(i)+c1*h*dydx(i)+c3*h*k3(i)+c4*h*k4(i)+c6*h*k6(i);

//Estimate error as difference between fourth and fifth order methods.
for (i=0;i<n;i++)
yerr(i)=h*(dc1*dydx(i)+dc3*k3(i)+dc4*k4(i)+dc5*k5(i)+dc6*k6(i));
}

//------------------ bcuderiv_square --------------------------------------------------------------------------------------
//Calculates the gradients y1 and y2 and the cross-derivative y12 numerically
//second order accuracy
//needs function y on the entire grid and the indices (j,k) of the lower left corner of the grid square
//and grid-distances d1 and d2 (equidistant grid required)
//efit_class version modified: use only for grid points at a time
//y(0...NR+1,0...NZ+1), y1(1...4), y2(1...4), y12(1...4)
void bcuderiv_square(Array<double,2>& y, int j, int k, double d1, double d2, 
					 Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12)
{
const double d1d2 = d1*d2;

y1(1) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(1) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(1) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;

j += 1;
y1(2) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(2) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(2) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;

k += 1;
y1(3) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(3) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(3) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;

j -= 1;
y1(4) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(4) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(4) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;
}

//------------------ bcuint_square ------------------------------------------------------------------------------------------------
//Bicubic interpolation. Input quantities are Ra and Za containing the grid coordinates,
//field contains the function values on the entire grid
//R and Z are the coordinates of the desired point for
//the interpolation. The interpolated function value is returned as y, and the interpolated
//gradient values as y1 and y2. 
//Modified from efit_class version: This routine calls bcucof and bcuderiv_square.
void bcuint_square(Array<double,1>& Ra, Array<double,1>& Za, double dR, double dZ, Array<double,2>& field,
			double R, double Z, double& y, double& y1, double& y2)
{
int i,j,k;
double t,u;
Array<double,2> c(Range(1,4),Range(1,4));
Range all = Range::all();

// Determine grid square where (R,Z) is in
j = int((R-Ra(1))/dR) + 1;
k = int((Z-Za(1))/dZ) + 1;
if(j>Ra.rows() || j<1 || k>Za.rows() || k<1)	{ofs2 << "Point outside of grid" << endl; EXIT;}

// Get derivatives
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
bcuderiv_square(field,j,k,dR,dZ,y1_sq,y2_sq,y12_sq);

// Get the c�s.
y_sq(1) = field(j,k); y_sq(2) = field(j+1,k); y_sq(3) = field(j+1,k+1); y_sq(4) = field(j,k+1);
bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,c);

// Interpolate
t=(R-Ra(j))/dR;
u=(Z-Za(k))/dZ;

y = y1 = y2 = 0.0;
for (i=4;i>=1;i--) 
{ 
	y = t*y + ((c(i,4)*u + c(i,3))*u + c(i,2))*u + c(i,1);
	y2 = t*y2 + (3.0*c(i,4)*u + 2.0*c(i,3))*u + c(i,2);
	y1 = u*y1 + (3.0*c(4,i)*t + 2.0*c(3,i))*t + c(2,i);
}

y1 /= dR;
y2 /= dZ;
}

//--------------- rungekutta5 ---------------------------------------------------------------------------------------------
// OLD VERSION
//Given values for the variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
//fifth-order Runge-Kutta method to advance the solution over an interval h and return the
//incremented variables as yout[0..n-1], which need not be a distinct array from y. The user
//supplies the routine dgls(x,y,dydx), which returns derivatives dydx at x.
//void rungekutta5(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout)
//{
//int i;
//double xh;
//double a2,a3,a4,a5,a6;
//double c1,c3,c4,c6;
//double b21,b31,b32,b41,b42,b43,b51,b52,b53,b54;
//double b61,b62,b63,b64,b65;
//
//a2=0.2;	a3=0.3;	a4=0.6;	a5=1.0;	a6=0.875;
//c1 = 37.0/378.0;  c3 = 250.0/621.0; c4 = 125.0/594.0; c6 = 512.0/1771.0;
//b21 = 0.2; 
//b31 = 3.0/40.0; b32 = 9.0/40.0; 
//b41 = 0.3; b42 = -0.9; b43 = 1.2;  
//b51 = -11.0/54.0; b52 = 2.5; b53 = -70.0/27.0; b54 = 35.0/27.0; 
//b61 = 1631.0/55296.0; b62 = 175.0/512.0; b63 = 575.0/13824.0; b64 = 44275.0/110592.0; b65 = 253.0/4096.0;
//
//Array<double,1> yt(n),k2(n),k3(n),k4(n),k5(n),k6(n);
//
////First step
//for (i=0;i<n;i++) yt(i)=y(i)+b21*h*dydx(i);
//
////Second step
//xh=x+h*a2;
//dgls(xh,yt,k2);
//for (i=0;i<n;i++) yt(i)=y(i)+b31*h*dydx(i)+b32*h*k2(i);
//
////Third step
//xh=x+h*a3;
//dgls(xh,yt,k3);
//for (i=0;i<n;i++) yt(i)=y(i)+b41*h*dydx(i)+b42*h*k2(i)+b43*h*k3(i);
//
////Forth step
//xh=x+h*a4;
//dgls(xh,yt,k4);
//for (i=0;i<n;i++) yt(i)=y(i)+b51*h*dydx(i)+b52*h*k2(i)+b53*h*k3(i)+b54*h*k4(i);  
//
////Fifth step
//xh=x+h*a5;
//dgls(xh,yt,k5);
//for (i=0;i<n;i++) yt(i)=y(i)+b61*h*dydx(i)+b62*h*k2(i)+b63*h*k3(i)+b64*h*k4(i)+b65*h*k5(i); 
//
////Last step
//xh=x+h*a6;
//dgls(xh,yt,k6);
//for (i=0;i<n;i++) yout(i)=y(i)+c1*h*dydx(i)+c3*h*k3(i)+c4*h*k4(i)+c6*h*k6(i);
//}

//--------------- Bfield_linf ---------------------------------------------------------------------------------------------
/* fast Biot-Savart integrator for linear conductor
taken from Paper Hanson & Hirshman, Phys. Plasmas 9, 4410 (2002)
parameters:
B: magnetic field
r: position at which the magnetic field should be calculated
rss: start position of the linear conductor
rse: end position of the linear conductor
I: current */
void Bfield_linf(vektor &B, vektor r, vektor rss, vektor rse, double I)
{
  double mu0d4pi,L,nrmrss,nrmrse,nsum,fak,nenner;
  vektor e,rmrss,rmrse;
  
  mu0d4pi = 1e-7;
  L = (rse-rss).norm();
  e = (rse-rss)/L;
  rmrss = r-rss;
  rmrse = r-rse;
  nrmrss = rmrss.norm();
  nrmrse = rmrse.norm();
  nsum = nrmrss+nrmrse;
  nenner = nsum*nsum-L*L;
  if (nenner == 0)
  {
    B.x = 0;
    B.y = 0;
    B.z = 0;
  } else {
    fak = 2.*mu0d4pi*I*L*nsum/nrmrss/nrmrse/nenner;
    B = e.cross(rmrss)*fak;
  }
}

//--------------- Bfield_rahmen ---------------------------------------------------------------------------------------------
/* calculates the magnetic field of a polygon with N edges
therefor it sum over N linear conductor by using Bfield_linf
parameters:
B: magnetic field
r: position at which the magnetic field should be calculated
re: array of all corners of the polygon
N: number of edges
I: current */
void Bfield_rahmen(vektor &B, vektor r, blitz::Array<vektor,1> re, int N, double I)
{
  vektor Btmp;
  
  B.x = 0;
  B.y = 0;
  B.z = 0;
  for (int i=0; i<(N-1); i++)
  {
    Bfield_linf(Btmp,r,re(i),re(i+1),I);
    B += Btmp;
  }
  Bfield_linf(Btmp,r,re(N-1),re(0),I);
  B += Btmp;
}

//--------------- Bfield_EFCC ---------------------------------------------------------------------------------------------
/* magnetic field of the EFCC of JET. Geometry of EFCC given by Dr. Liang.
parameters:
B: magnetic field
r: position at which the magnetic field should be calculated
Current and coil geometry is stored in global variables*/
//IMPORTANT for n=1: I = [+ + + +]; n=2: I = [- + + -];
void Bfield_EFCC(double &Bx, double &By, double &Bz, double x, double y, double z)
{
  vektor B,Btmp,r;
  
  B.set(0,0,0);
  r.set(x,y,z);
  
  // EFCC1:
  Bfield_rahmen(Btmp,r,EFCC1,14,curntEFCC(0));
  B += Btmp;
  
  // EFCC3:
  Bfield_rahmen(Btmp,r,EFCC3,16,curntEFCC(1));
  B += Btmp;
  
  // EFCC5:
  Bfield_rahmen(Btmp,r,EFCC5,14,curntEFCC(2));
  B += Btmp;
  
  // EFCC7:
  Bfield_rahmen(Btmp,r,EFCC7,14,curntEFCC(3));
  B += Btmp;
  
  Bx = B.x;
  By = B.y;
  Bz = B.z;
}

//--------------- Bfield_EF ---------------------------------------------------------------------------------------------
/* magnetic field of the Error Fields of JET, simulated by EFCC.. Geometry of EFCC given by Dr. Liang.
parameters:
B: magnetic field
r: position at which the magnetic field should be calculated
Current and coil geometry is stored in global variables*/
//IMPORTANT for n=1: I = [+ + + +]; n=2: I = [- + + -]; Coils are implementated in strange directions ;-)
void Bfield_EF(double &Bx, double &By, double &Bz, double x, double y, double z)
{
  vektor B,Btmp,r;
  
  B.set(0,0,0);
  r.set(x,y,z);
  
  // EFCC1:
  Bfield_rahmen(Btmp,r,EFCC1,14,curntEF(0));
  B += Btmp;
  
  // EFCC3:
  Bfield_rahmen(Btmp,r,EFCC3,16,curntEF(1));
  B += Btmp;
  
  // EFCC5:
  Bfield_rahmen(Btmp,r,EFCC5,14,curntEF(2));
  B += Btmp;
  
  // EFCC7:
  Bfield_rahmen(Btmp,r,EFCC7,14,curntEF(3));
  B += Btmp;
  
  Bx = B.x;
  By = B.y;
  Bz = B.z;
}

//----------- outofBndy ----------------------------------
// Check if (x,y) is out of the torus. Returns 0 if (x,y) 
// is in boundary an 1 if (x,y) is out of boundary. 
// simpleBndy = 0; use real boundaries
// simpleBndy = 1: use rectangle boundaries
bool outofBndy(double x, double y, int simpleBndy)
{
  switch(simpleBndy)
  {
    case 0:
      return outofRealBndy(x,y);
      break;
    case 1:
      return (x<bndy[0] || x>bndy[1] || y<bndy[2] || y>bndy[3]);	// bndy[4] = {x_min, x_max, y_min, y_max}
      break;
    default:
      cout << "simpleBndy switch has a wrong value!" << endl;
  }
}

//------------ outofRealBndy -------------------------------
// Check if (x,y) is out of the torus. It uses the jordan curve theorem.
bool outofRealBndy(double x, double y)
{
  int wn = 0;

  double x1,y1;
  double x2 = EQD.wall(1);	//R1
  double y2 = EQD.wall(2);	//Z1

  bool startUeber = (y2 >= y) ? 1 : 0;
  for(int i=3; i<(2*EQD.Nwall); i=i+2)
  {
    x1 = x2;
    y1 = y2;
    x2 = EQD.wall(i);
    y2 = EQD.wall(i+1);
    
    /*if((y1==y2) && (y==y1))
    {
      if(((x1<=x)&&(x<=x2)) || ((x2<=x)&&(x<=x1))) return 0;
    } 
    else if((x1==x2) && (x==x1))
    {
      if(((y1<=y)&&(y<=y2)) || ((y2<=y)&&(y<=y1))) return 0;
    } else {
      a = (x-x1)*(y2-y1)-(y-y1)*(x2-x1);
      if((a <= 1e-15) && (a >= -1e-15)) return 0;	//necessary for use in dtfoot
    }*/
    bool endUeber = (y2 >= y) ? 1 : 0;
    if(startUeber != endUeber) 
    {
      if((y2 - y)*(x2 - x1) <= (y2 - y1)*(x2 - x))
      {
        if(endUeber) { wn++; }
      } else {
        if(!endUeber) { wn--; }
      }
    }
    startUeber = endUeber;
  }
  return wn == 0;
}
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

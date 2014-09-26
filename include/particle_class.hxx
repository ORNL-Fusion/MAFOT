// PARTICLE Class contains all parameters regarding position and trajectory
// needs input variables from EFIT class and IO class
// includes member functions to integrate particle trajectory
// needs extern subroutine 'getBfield' as declared in the Prototypes
// needs extern subroutine 'outofBndy' as declared in the Prototypes
// used by all D3D, ITER and NSTX drift programs
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						1.06.11


// Define
//--------
#ifndef PARTICLE_CLASS_INCLUDED
#define PARTICLE_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;

// extern Prototypes, defined in Machine-specific header
#ifdef USE_SIESTA
	void getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR, SIESTA& SIES);
#else
	void getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
#endif
bool outofBndy(double x, double y, EFIT& EQD);	

// Prototypes  
void get_Energy(double psi, double& Enorm, double& dEnorm);
double get_Lmfp(double Ekin);
void getRZ(double x, double y, double& R, double &Z, EFIT& EQD);
double bisec(double psi, double theta, double a, double b, EFIT& EQD, int flag);

// Integrator Parameters
const int nvar = 2;				// Number of Variables
const int ilt = 360;			// Steps till Output
const double dpinit = 1.0;		// step size of phi in [deg]

// Golbal Parameters 
extern ofstream ofs2;

//--------- Begin Class PARTICLE ----------------------------------------------------------------------------------------------
class PARTICLE
{
private:
// Member Variables
	EFIT& EQDr;		// Only a Reference, not a copy
	IO& PARr;		// Only a Reference, not a copy
#ifdef USE_SIESTA
	SIESTA& SIESr;	// Only a Reference, not a copy
#endif

	double GAMMA;
	double eps0;
	double Ix;
	double mc2;
	int steps;		// total number of integration steps along trajectory

// Member-Functions
	void dgls(double x, Array<double,1> y, Array<double,1>& dydx);
	int rkint(int nvar, int nstep, double dx, Array<double,1>& y, double& x);
	void rungekutta4(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout);

public: 
// Member Variables
	double R;		// cylindrical major radius [m]
	double Z;		// cylindrical vertical coordinate [m]
	double phi;		// toroidal angle [deg]
	double psi;		// normalized flux

	double Lc;		// connection length [m]
	double psimin;	// minimum of normalized flux reached by trajectory
	double psimax;	// maximum of normalized flux reached by trajectory
	double psiav;	// average normalized flux along trajectory

	double Ekin;	// kinitic particle energy [keV]
	int sigma;		// 1: co-passing particles		-1: count-passing particles		0: field lines only
	int Zq;			// Charge number: 1: ions are calculated	-1: electrons are calculated

	double Lmfp_total;	// Sum of all mean free paths along the trajectory (PARr.useTprofile == 1 only)

// Constructors
#ifdef USE_SIESTA
	PARTICLE(EFIT& EQD, IO& PAR, SIESTA& SIES, int mpi_rank=0);		// Default Constructor
#else
	PARTICLE(EFIT& EQD, IO& PAR, int mpi_rank=0);		// Default Constructor
#endif

// Member-Operators
	PARTICLE& operator =(const PARTICLE& FLT);								// Operator =
	friend std::ostream& operator <<(std::ostream& out, PARTICLE& FLT);		// Operator <<

// Member-Functions
	double get_r();
	double get_theta();
	void set_Energy();
	int mapit(const int itt,int MapDirection=1);
	int mapstep(int MapDirection=1, int nstep=ilt);
	int connect(double& ntor, double& length, double& psimintotal, double& psimaxtotal, double& psiavtotal, const int itt, int MapDirection=0);	// with psimax
	int connect(double& ntor, double& length, double& psimintotal, const int itt, int MapDirection=0);						// without psimax (old version, kept for compatibility reasons)

	void set(int i, int N, double Rmin, double Rmax, double Zmin, double Zmax, int N_Z=1, int flag=0);
	void create(long& idum, double Rmin, double Rmax, double Zmin, double Zmax, int flag=0);
	void convertRZ(double theta, double r);

}; //end of class

//------------------------ End of Class -----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Default Constructor -------------------------------------------------------------------------------------------
// Constructor list necessary to set EQDr and PARr since references cannot be empty
#ifdef USE_SIESTA
	PARTICLE::PARTICLE(EFIT& EQD, IO& PAR, SIESTA& SIES, int mpi_rank): EQDr(EQD), PARr(PAR), SIESr(SIES)
#else
	PARTICLE::PARTICLE(EFIT& EQD, IO& PAR, int mpi_rank): EQDr(EQD), PARr(PAR)
#endif
{
// Variables
double omegac;

const double c=299792458;			// speed of light in m/s
const double e=1.602176462*1e-19;	// elementary charge in C
const double me=9.10938188*1e-31;	// Electron mass in kg
const double mp=1.67262158*1e-27;	// Proton mass in kg
const double E0e=me*c*c/e/1000.0;	// Rest energy of Electrons in [keV] 
const double E0p=mp*c*c/e/1000.0;	// Rest energy of Protons in [keV] 
const int Massnumber=1;				// Mass number

// Public Member Variables

R = 0;				// cylindrical major radius [m]
Z = 0;				// cylindrical vertical coordinate [m]
phi = 0;			// toroidal angle [deg]
psi = 0;			// normalized flux

Lc = 0;				// connection length [m]
psimin = 10;		// minimum of normalized flux reached by trajectory
psimax = 0;			// maximum of normalized flux reached by trajectory
psiav = 0;
steps = 0;

Ekin = PAR.Ekin;	// kinitic particle energy
sigma = PAR.sigma;	// 1: co-passing particles		-1: count-passing particles		0: field lines only
Zq = PAR.Zq;		// Charge number: 1: ions are calculated	-1: electrons are calculated

Lmfp_total = 0;

// dependent Private Member Variables
if(sigma == 0)
{
	if(mpi_rank < 1) cout << "Field lines are calculated" << endl;
	ofs2 << "Field lines are calculated" << endl;
	GAMMA = 1;	omegac = 1;	 eps0 = 1;	Ix = 1;  mc2 = 1;
}
else
{
	if(Zq >= 1) // Ions
	{
		mc2 = E0p*Massnumber;						// Rest Energy
		omegac = e*EQD.Bt0/(mp*Massnumber);		// normalized gyro frequency (SI-System)
		if(mpi_rank < 1) cout << "Ions are calculated" << endl;
		ofs2 << "Ions are calculated" << endl;
	}
	else // Electrons
	{
		Zq = -1;									// default!
		mc2 = E0e*Massnumber;						// Rest Energy
		omegac = e*EQD.Bt0/(me*Massnumber);		// normalized gyro frequency (SI-System)
		if(mpi_rank < 1) cout << "Electrons are calculated" << endl;
		ofs2 << "Electrons are calculated" << endl;
	}
	GAMMA = 1 + Ekin/mc2;							// relativistic gamma factor 1/sqrt(1-v^2/c^2)
	eps0 = c*c/omegac/omegac/EQD.R0/EQD.R0;		// normalized rest energy
	Ix = -0.5/double(Zq)*eps0*((PAR.lambda*(GAMMA-1)+1)*(PAR.lambda*(GAMMA-1)+1)-1);
	if(mpi_rank < 1) cout << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
	ofs2 << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
}
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
PARTICLE& PARTICLE::operator =(const PARTICLE& FLT)
{
if (this == &FLT) return(*this);	    // if: x=x
// Private Member Variables
// References to EQD and PAR remain unchanged! There is only one EQD and one PAR
GAMMA = FLT.GAMMA;
eps0 = FLT.eps0;
Ix = FLT.Ix;
mc2 = FLT.mc2;
steps = FLT.steps;

// Public Member Variables
R = FLT.R;	
Z = FLT.Z;	
phi = FLT.phi;	
psi = FLT.psi;	

Lc = FLT.Lc;	
psimin = FLT.psimin;
psimax = FLT.psimax;
psiav = FLT.psiav;

Ekin = FLT.Ekin;
sigma = FLT.sigma;
Zq = FLT.Zq;

Lmfp_total = FLT.Lmfp_total;
return(*this);
}

//--------- Operator << ---------------------------------------------------------------------------------------------------
ostream& operator <<(ostream& out, PARTICLE& FLT)
{
out << "--- Position ---" << endl;
out << "R = " << FLT.R << endl;	
out << "phi = " << FLT.phi << endl;	
out << "Z = " << FLT.Z << endl;	
out << "r = " << FLT.get_r() << endl;
out << "theta = " << FLT.get_theta() << endl;
out << "psi = " << FLT.psi << endl;	

out << "--- Trajectory ---" << endl;
out << "Lc = " << FLT.Lc << endl;	
out << "psimin = " << FLT.psimin << endl;
out << "psimax = " << FLT.psimax << endl;
out << "psiav = " << FLT.psiav/FLT.steps << endl;

out << "--- Properties ---" << endl;
out << "Ekin = " << FLT.Ekin << endl;
out << "Type (0=field line) = " << FLT.sigma << endl;
out << "Charge = " << FLT.Zq << endl;

return out;
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------- get_r --------------------------------------------------------------------------------------------------
double PARTICLE::get_r()
{
double Rm,Zm;
Rm = R - EQDr.RmAxis; // R relative to magnetic axis
Zm = Z - EQDr.ZmAxis; // Z relative to magnetic axis
return sqrt(Rm*Rm + Zm*Zm);
}

//---------------- get_theta ----------------------------------------------------------------------------------------------
double PARTICLE::get_theta()
{
double Rm,Zm,theta;
Rm = R - EQDr.RmAxis; // R relative to magnetic axis
Zm = Z - EQDr.ZmAxis; // Z relative to magnetic axis

theta = atan(Zm/Rm);
if(Rm<0) theta += pi;
if(Rm>=0 && Zm<0) theta += pi2;

return theta;
}

//---------------- set_Energy ---------------------------------------------------------------------------------------------
void PARTICLE::set_Energy()
{
double Enorm,dummy;

get_Energy(psi,Enorm,dummy);
Ekin = Enorm*PARr.Ekin;
GAMMA = 1 + Ekin/mc2;
Ix = -0.5/double(Zq)*eps0*((PARr.lambda*(GAMMA-1)+1)*(PARr.lambda*(GAMMA-1)+1)-1);
}

//---------------- mapit --------------------------------------------------------------------------------------------------
int PARTICLE::mapit(const int itt, int MapDirection)
{
int chk = 0;
const double phistart = phi;

// Integrate
for(int i=1;i<=itt;i++)
{
	chk = mapstep(MapDirection);
	if(chk<0) {break;}	// particle has left system

	if(fabs(phi-MapDirection*i*dpinit*ilt-phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(phi-MapDirection*i*dpinit*ilt-phistart) << endl;
	phi=MapDirection*i*dpinit*ilt+phistart;
}
return chk;
}

//---------------- mapstep ------------------------------------------------------------------------------------------------
int PARTICLE::mapstep(int MapDirection, int nstep)
{
int chk;
double phi_rad = phi/rTOd;	// phi in rad
const double dphi = MapDirection*dpinit/rTOd;
Array<double,1> y(nvar); // Array to set initial conditions

y(0) = R;
y(1) = Z;

// integrate one full toroidal turn
chk = rkint(nvar,nstep,dphi,y,phi_rad);	
if(chk<0) return -1;	// particle has left system

R = y(0);
Z = y(1);
phi = phi_rad*rTOd;		// phi back in deg

return 0;
}

//------------------ connect ----------------------------------------------------------------------------------------------
// Integration goes in respecive direction, depending on MapDirection
// MapDirection=0 means both directions are calculated and results are added
// phistart has to be in deg!
int PARTICLE::connect(double& ntor, double& length, double& psimintotal, double& psimaxtotal, double& psiavtotal, const int itt, int MapDirection)
{
int chk;
const double Rstart = R;
const double Zstart = Z;
const double phistart = phi;
const double psistart = psi;

const double Estart = Ekin;
const double GAMMAstart = GAMMA;
const double Ixstart = Ix;
const double Lmfpstart = Lmfp_total;

// positive phi direction
if(MapDirection >= 0)
{
	chk = mapit(itt,1);
	ntor = fabs(phi-phistart)/360.0;
	length = Lc;
	psimintotal = psimin;
	psimaxtotal = psimax;
	psiavtotal = psiav;
}
else
{
	ntor = 0;
	length = 0;
	psimintotal = 10;
	psimaxtotal = 0;
	psiavtotal = 0;
}

// negative phi direction
R = Rstart;
Z = Zstart;
phi = phistart;
psi = psistart;
Lc = 0;  
psimin = 10;
psimax = 0;
psiav = 0;

Ekin = Estart;
GAMMA = GAMMAstart;
Ix = Ixstart;
Lmfp_total = Lmfpstart;

if(MapDirection <= 0)
{
	chk = mapit(itt,-1);
	ntor += fabs(phi-phistart)/360.0;
	length += Lc;
	if(psimin < psimintotal) psimintotal = psimin;
	if(psimax > psimaxtotal) psimaxtotal = psimax;
	psiavtotal += psiav;
}
psiavtotal /= steps;

return 0;
}

//------------------ connect ----------------------------------------------------------------------------------------------
// Integration goes in respecive direction, depending on MapDirection
// MapDirection=0 means both directions are calculated and results are added
// phistart has to be in deg!
int PARTICLE::connect(double& ntor, double& length, double& psimintotal, const int itt, int MapDirection)
{
int chk;
const double Rstart = R;
const double Zstart = Z;
const double phistart = phi;
const double psistart = psi;

const double Estart = Ekin;
const double GAMMAstart = GAMMA;
const double Ixstart = Ix;
const double Lmfpstart = Lmfp_total;

// positive phi direction
if(MapDirection >= 0)
{
	chk = mapit(itt,1);
	ntor = fabs(phi-phistart)/360.0;
	length = Lc;
	psimintotal = psimin;
}
else
{
	ntor = 0;
	length = 0;
	psimintotal = 10;
}

// negative phi direction
R = Rstart;
Z = Zstart;
phi = phistart;
psi = psistart;
Lc = 0;
psimin = 10;

Ekin = Estart;
GAMMA = GAMMAstart;
Ix = Ixstart;
Lmfp_total = Lmfpstart;

if(MapDirection <= 0)
{
	chk = mapit(itt,-1);
	ntor += fabs(phi-phistart)/360.0;
	length += Lc;
	if(psimin < psimintotal) psimintotal = psimin;
}

return 0;
}

//-------------- set ------------------------------------------------------------------------------------------------------
// creates initial condition i: (R,Z) on a grid with N points
// values taken from the intervals [Rmin,Rmax] and [Zmin,Zmax] respectively
// if N_Z is spezified, N=N_Z*N_R has to be used -> Grid size is then N_Z X N_R
// otherwise Grid size is about sqrt(N) X sqrt(N)
// R is varied first, Z second
void PARTICLE::set(int i, int N, double Rmin, double Rmax, double Zmin, double Zmax, int N_Z, int flag)
{
double dZ,dR,dummy;
double x,y;
int i_Z=0,i_R=0;
int N_R;

if(N<=1) {dR=0; dZ=0;}
else
{
	dR=(Rmax-Rmin)/(N-1);
	dZ=(Zmax-Zmin)/(N-1);
}

if(dR==0) i_Z=i-1;
if(dZ==0) i_R=i-1;
if(dZ!=0 && dR!=0) 
{
	if(N_Z<=1) N_Z=int(sqrt(double(N))+0.5);
	N_R=int(double(N)/double(N_Z)+0.5);
	i_R=(i-1)%N_R;
	i_Z=int(double(i-1)/double(N_R));
	dR=(Rmax-Rmin)/(N_R-1);
	dZ=(Zmax-Zmin)/(N_Z-1);
}
x = Rmin + i_R*dR;
y = Zmin + i_Z*dZ;

phi = PARr.phistart;

switch(flag)
{
case 2:		// get R, Z from x = psi and y = theta
#ifdef USE_SIESTA
	if(PARr.response_field == -2) SIESr.get_RZ(x, y, phi, R, Z);
	else getRZ(x, y, R, Z, EQDr);
#else
	getRZ(x, y, R, Z, EQDr);
#endif
	psi = x;
	break;
case 1:		// get R, Z from r and theta
	R = x*cos(y) + EQDr.RmAxis;
	Z = x*sin(y) + EQDr.ZmAxis;
	EQDr.get_psi(R,Z,psi,dummy,dummy);
	break;
default:
	R = x;
	Z = y;
	EQDr.get_psi(R,Z,psi,dummy,dummy);
	break;
}

Lc = 0;
psimin = 10;
psimax = 0;
psiav = 0;
steps = 0;

if(sigma != 0 && PARr.useTprofile == 1) {set_Energy(); Lmfp_total = get_Lmfp(Ekin);}
}

//-------------- create ---------------------------------------------------------------------------------------------------
// creates an initial condition (R,Z) using random numbers
// values taken from the intervals [Rmin,Rmax] and [Zmin,Zmax] respectively
// R and Z are output variables
void PARTICLE::create(long& idum, double Rmin, double Rmax, double Zmin, double Zmax, int flag)
{
const double dR=Rmax-Rmin;
const double dZ=Zmax-Zmin;
double v,dummy;
double x,y;

v=ran0(idum);
x=Rmin+v*dR;

v=ran0(idum);
y=Zmin+v*dZ;

phi = PARr.phistart;

switch(flag)
{
case 2:		// get R, Z from psi and theta
#ifdef USE_SIESTA
	if(PARr.response_field == -2) SIESr.get_RZ(x, y, phi, R, Z);
	else getRZ(x, y, R, Z, EQDr);
#else
	getRZ(x, y, R, Z, EQDr);
#endif
	psi = x;
	break;
case 1:		// get R, Z from r and theta
	R = x*cos(y) + EQDr.RmAxis;
	Z = x*sin(y) + EQDr.ZmAxis;
	EQDr.get_psi(R,Z,psi,dummy,dummy);
	break;
default:
	R = x;
	Z = y;
	EQDr.get_psi(R,Z,psi,dummy,dummy);
	break;
}

Lc = 0;
psimin = 10;
psimax = 0;
psiav = 0;
steps = 0;

if(sigma != 0 && PARr.useTprofile == 1) {set_Energy(); Lmfp_total = get_Lmfp(Ekin);}
}

//-------------- convertRZ ---------------------------------------------------------------------------------------------------
// Sets R and Z member variables from polar coordinates theta and r
void PARTICLE::convertRZ(double theta, double r)
{
	R = r*cos(theta) + EQDr.RmAxis;
	Z = r*sin(theta) + EQDr.ZmAxis;
}

//----------------------- End of Public Member Functions ------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- Private Member Functions ------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------- dgls ---------------------------------------------------------------------------------------------------
// Type the differential equations (dgls) as they are written, while x is the independent variable
// Here: x = phi, y(0) = R, y(1) = Z, dydx(0) = dR/dphi, dydx(1) = dZ/dphi
void PARTICLE::dgls(double x, Array<double,1> y, Array<double,1>& dydx)
{
double B_R,B_Z,B_phi;
double S;
#ifdef USE_SIESTA
	getBfield(y(0),y(1),x,B_R,B_Z,B_phi,EQDr,PARr,SIESr);
#else
	getBfield(y(0),y(1),x,B_R,B_Z,B_phi,EQDr,PARr);
#endif

dydx(0) = y(0)*B_R/B_phi;
dydx(1) = y(0)*B_Z/B_phi;

if(sigma != 0)
{
	S = eps0*(GAMMA*GAMMA-1)-2*EQDr.R0*Ix/y(0);
	if(S<0) {ofs2 << "dgls: Sqrt argument negative => Abort" << endl; EXIT;}	// Error -> Abort program
	S = sqrt(S);
	dydx(1) += -sigma/double(Zq)*(y(0)*S + EQDr.R0*Ix/S);	// sign corrected: sigma = +1 is indeed co-passing
}
}

//----------------- rkint -------------------------------------------------------------------------------------------------
//Starting from initial values y[0..nvar-1] known at x=x1 use Runge-Kutta
//to advance nstep equal increments to x2=x1+nstep*dx. The user-supplied routine dgls(x,v,dvdx)
//evaluates derivatives. Results after nstep Steps are stored in y[0..nvar-1]
int PARTICLE::rkint(int nvar, int nstep, double dx, Array<double,1>& y, double& x)
{
int k;
double x1 = x;	//Store first value (helps reduce Error in x)
double dummy;
double Lmfp;
Array<double,1> yout(nvar),dydx(nvar);

//Take nstep steps
for (k=1;k<=nstep;k++) 
{ 
	dgls(x,y,dydx);
	rungekutta4(y,dydx,nvar,x,dx,yout);
	x = x1 + k*dx; // Better than x+=dx

	// Integration terminates outside of boundary box
	if(outofBndy(yout(0),yout(1),EQDr) == true) return -1;

	// Get additional Parameter
	Lc += sqrt((yout(0)-y(0))*(yout(0)-y(0)) + (yout(1)-y(1))*(yout(1)-y(1)) + 0.25*(yout(0)+y(0))*(yout(0)+y(0))*dx*dx);
	EQDr.get_psi(yout(0),yout(1),psi,dummy,dummy);
	if(psi < psimin) psimin = psi;
	if(psi > psimax) psimax = psi;
	psiav += psi;
	steps += 1;

	y = yout;

	// Set Particle Energy
	if(sigma != 0 && PARr.useTprofile == 1 && Lc > Lmfp_total)			// else: already set in Constructor
	{
		set_Energy();
		Lmfp = get_Lmfp(Ekin);
		Lmfp_total += Lmfp;
	}
} 
return 0;
}

//--------------- rungekutta4 ---------------------------------------------------------------------------------------------
//Given values for the variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
//fourth-order Runge-Kutta method to advance the solution over an interval h and return the
//incremented variables as yout[0..n-1], which need not be a distinct array from y. The user
//supplies the routine dgl(x,y,dydx), which returns derivatives dydx at x.
void PARTICLE::rungekutta4(Array<double,1> y, Array<double,1> dydx, int n, double x, double h, Array<double,1>& yout)
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

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------- get_Energy ----------------------------------------------
void get_Energy(double psi, double& Enorm, double& dEnorm)
{
double xs = 0.975;	// Symmetry point in Pedestal
double dw = 0.04;	// half of Pedestal width

Enorm = 0.5*tanh(2*(xs-psi)/dw) + 0.5;
dEnorm = -1/dw*(1-(2*Enorm-1)*(2*Enorm-1));
}

//--------- get_Lmfp -------------------------------------------------
double get_Lmfp(double Ekin)
{
const double L_Debye = sqrt(Ekin*1000)*7.43e-7;			// Debye length [m], density 1e14 cm^-3, Ekin in [keV], needed in [eV]
const double Lambda_Coulomb = 4*pi*1e20*pow(L_Debye,3);	// density set to constant 1e20 m^-3 = 1e14 cm^-3
return 64*pi*L_Debye*Lambda_Coulomb/log(Lambda_Coulomb);
}

//--------- getRZ ----------------------------------------------------
// calculate (R,Z) from (theta,psi), x = psi, y = theta
// exclude the private flux region <--> two points (R,Z) with the same psi exist
// optimized for lower single null discharges
// wall(Range(1,toEnd,2)) are all R coordinates of wall, same with lcfs
// wall(Range(2,toEnd,2)) are all Z coordiantes of wall, same with lcfs
void getRZ(double x, double y, double& R, double &Z, EFIT& EQD)
{
if(y < pi/4 || y > 7*pi/4)
{
	// Z = Z(R,theta)
	if(x <= 1) R = bisec(x, y, EQD.RmAxis, max(EQD.lcfs(Range(1,toEnd,2)))+EQD.dR, EQD, 0);
	else R = bisec(x, y, EQD.RmAxis, max(EQD.wall(Range(1,toEnd,2))), EQD, 0);
	Z = (R - EQD.RmAxis)*tan(y) + EQD.ZmAxis;
}
if(y > 3*pi/4 && y < 5*pi/4)
{
	// Z = Z(R,theta)
	if(x <= 1) R = bisec(x, y, min(EQD.lcfs(Range(1,toEnd,2)))-EQD.dR, EQD.RmAxis, EQD, 0);
	else R = bisec(x, y, min(EQD.wall(Range(1,toEnd,2))), EQD.RmAxis, EQD, 0);
	Z = (R - EQD.RmAxis)*tan(y) + EQD.ZmAxis;
}
if(y >= pi/4 && y <= 3*pi/4)
{
	// R = R(Z,theta)
	if(x <= 1) Z = bisec(x, y, EQD.ZmAxis, max(EQD.lcfs(Range(2,toEnd,2)))+EQD.dZ, EQD, 1);
	else Z = bisec(x, y, EQD.ZmAxis, max(EQD.wall(Range(2,toEnd,2))), EQD, 1);
	R = (Z - EQD.ZmAxis)/tan(y) + EQD.RmAxis;
}
if(y >= 5*pi/4 && y <= 7*pi/4)
{
	// R = R(Z,theta)
	if(x <= 1) Z = bisec(x, y, min(EQD.lcfs(Range(2,toEnd,2))), EQD.ZmAxis, EQD, 1);
	else Z = bisec(x, y, min(EQD.wall(Range(2,toEnd,2))), EQD.ZmAxis, EQD, 1);
	R = (Z - EQD.ZmAxis)/tan(y) + EQD.RmAxis;
}
}

//--------- bisec ----------------------------------------------------
double bisec(double psi, double theta, double a, double b, EFIT& EQD, int flag)
{
double xo,xu,x;
double R,Z;
double f,dummy;
const double eps = 1e-14;

x = a;
if (flag == 0)	// Z = Z(R,theta)
{
	Z = (x - EQD.RmAxis)*tan(theta) + EQD.ZmAxis;
	if(outofBndy(x, Z, EQD)) f = 1.2;
	else EQD.get_psi(x, Z, f, dummy, dummy);
	f -= psi;
}
else			// R = R(Z,theta)
{
	R = (x - EQD.ZmAxis)/tan(theta) + EQD.RmAxis;
	if(outofBndy(R, x, EQD)) f = 1.2;
	else EQD.get_psi(R, x, f, dummy, dummy);
	f -= psi;
}

if(f > 0)
{
	xo = a;
	xu = b;
}
else
{
	xo = b;
	xu = a;
}

while(fabs(xo-xu) > eps)
{
	x = (xo + xu)/2;
	if (flag == 0)	// Z = Z(R,theta)
	{
		Z = (x - EQD.RmAxis)*tan(theta) + EQD.ZmAxis;
		if(outofBndy(x, Z, EQD)) f = 1.2;
		else EQD.get_psi(x, Z, f, dummy, dummy);
		f -= psi;
	}
	else			// R = R(Z,theta)
	{
		R = (x - EQD.ZmAxis)/tan(theta) + EQD.RmAxis;
		if(outofBndy(R, x, EQD)) f = 1.2;
		else EQD.get_psi(R, x, f, dummy, dummy);
		f -= psi;
	}
	if(f > 0) xo = x;
	else xu = x;
}

return x;
}

#endif //  PARTICLE_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

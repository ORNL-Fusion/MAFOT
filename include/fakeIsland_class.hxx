// Fake Islands perturbation
// A.Wingen						17.5.16


// Define
//--------
#ifndef fakeIsland_CLASS_INCLUDED
#define fakeIsland_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <andi.hxx>
#include <efit_class.hxx>
using namespace blitz;

// Prototypes
//-----------

// Golbal Parameters
//------------------

//--------- Begin Class fakeIsland ----------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// creates a fake perturbation field that opens islands on rational surfaces
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class fakeIsland
{
private:
	// Member Variables
	int N;					// number of modes

	// Member-Functions

public:
	// Member Variables
	Array<double,1> m;			// poloidal mode number
	Array<double,1> n;			// toroidal mode number
	Array<double,1> A;		// mode amplitude
	Array<double,1> delta;	// mode phase
	Array<double,1> psi0;	// psi location of q(psi0) = m/n

	// Constructors
	fakeIsland();								// Default Constructor

	// Member-Operators
	fakeIsland& operator =(const fakeIsland& d);	// Operator =

	// Member-Functions
	void read(LA_STRING filename);
	void get_surfaces(EFIT& EQD);
	void get_B(double R, double phi, double Z, double& BR, double &BZ, EFIT& EQD);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
fakeIsland::fakeIsland()
{
TinyVector <int,1> index(1);
N = 1;
A.resize(N);		A.reindexSelf(index);
m.resize(N);		m.reindexSelf(index);
n.resize(N);		n.reindexSelf(index);
delta.resize(N);	delta.reindexSelf(index);
psi0.resize(N);		psi0.reindexSelf(index);
}


//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
fakeIsland& fakeIsland::operator = (const fakeIsland& d)
{
if (this == &d) return(*this);	    // if: x=x
N = d.N;
A.reference(d.A);
m.reference(d.m);
n.reference(d.n);
delta.reference(d.delta);
psi0.reference(d.psi0);
return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// read config file
void fakeIsland::read(LA_STRING filename)
{
Range all = Range::all();
Array<double,2> data;
readfile(filename, 4, data);
N = data.rows();
A.resize(N); 		A = data(all,1);
m.resize(N); 		m = data(all,2);
n.resize(N); 		n = data(all,3);
delta.resize(N); 	delta = data(all,4);
}

//-------------------------------------------------------------------------------------------------------------------------
// find location of m/n rational surfaces in psi space
void fakeIsland::get_surfaces(EFIT& EQD)
{
int i;
double psi,a,b,f;
const double eps = 1e-12;
psi0.resize(N);

for(i=1;i<=N;i++)
{
	a = 0; b = 1;
	while(fabs(b-a) > eps)
	{
		psi = 0.5*(a+b);
		f = EQD.get_q(psi) - fabs(m(i)/double(n(i)));
		if(f < 0) a = psi;
		else b = psi;
	}
	psi0(i) = psi;
}
}

//-------------------------------------------------------------------------------------------------------------------------
// return magnetic field at any point (R,phi,Z)
void fakeIsland::get_B(double R, double phi, double Z, double& BR, double &BZ, EFIT& EQD)
{
int i, chk;
double psi,theta,dpsidr,dpsidz,f,ph,dthdr,dthdz,r2;
double Ap = 0;			// vector potential, toroidal component only
double dAdp = 0;		// dAp/dpsi
double dAdth = 0;		// dAp/dtheta
const double width = 0.1;	// localization (= width of Gaussian) of perturbation in psi space

theta = polar_phi(R-EQD.RmAxis, Z-EQD.ZmAxis);
r2 = (R-EQD.RmAxis)*(R-EQD.RmAxis) + (Z-EQD.ZmAxis)*(Z-EQD.ZmAxis);
dthdr = -(Z-EQD.ZmAxis)/r2;
dthdz = (R-EQD.RmAxis)/r2;
chk = EQD.get_psi(R,Z,psi,dpsidr,dpsidz);

for(i=1;i<=N;i++)
{
	f = A(i)*exp(-(psi-psi0(i))*(psi-psi0(i))/width/width);
	ph = m(i)*theta - n(i)*phi + delta(i);
	Ap += f*cos(ph);
	dAdp += -2*(psi-psi0(i))/width/width*f*cos(ph);
	dAdth += -f*m(i)*sin(ph);

}
BR = -dAdp*dpsidz - dAdth*dthdz;		// -dAp/dZ
BZ = Ap/R + dAdp*dpsidr + dAdth*dthdr;	// 1/R d(R*Ap)/dR = Ap/R + dAp/dR
}

//------------------------ End of Class fakeIsland ------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


#endif //  fakeIsland_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

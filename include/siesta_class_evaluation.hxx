// SIESTA Class reads and interpolates B-fields from SIESTA
// used by all D3D, ITER and NSTX drift programs
// uses arrays and multiple-arrays from blitz-Library
// all arrays start with index 1
// A.Wingen						30.07.14


// Define
//--------
#ifndef SIESTA_CLASS_INCLUDED
#define SIESTA_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <netcdf.h>
#include <andi.hxx>
#include <efit_class.hxx>
#include <vmec_class.hxx>
using namespace blitz;

// Prototypes
//-----------
void spline(Array<double,1>& x, Array<double,1>& y, int n, double yp1, double ypn, Array<double,1>& y2);
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, int n, double x, double& y, double& yx);

// Golbal Parameters
//------------------
extern ofstream ofs2;

//--------- Begin Class SIESTA_SPECTRAL -----------------------------------------------------------------------------------
class SIESTA_SPECTRAL
{
private:
	// Parameter
	static const double sqrt2 = 1.4142135623730951;

	// Member Variables
	double ds;		// s array step size
	Array<double,3> d2ymns;	// d^2/ds^2 ymns(s), output of spline and input of splint
	Array<double,3> d2ymnc;	// d^2/ds^2 ymnc(s), output of spline and input of splint

	// Member-Functions
	double orthonorm(int n, int m);

public:
	// Member Variables
	LA_STRING id; 	// name of spectral variable

	int parity;		// parity = 1 <-> cosine,	parity = -1 <-> sine, 	parity = 0 <-> both
	int ns;			// number of s points
	int ntor;		// number of toroidal modes n	-ntor -> ntor
	int mpol;		// number of poloidal modes m	0 -> mpol

	Array<double,1> S;	// s array

	Array<double,3> ymns; 	// sine spectral data
	Array<double,3> ymnc; 	// cosine spectral data

	// Constructors
	SIESTA_SPECTRAL();								// Default Constructor

	// Member-Operators
	SIESTA_SPECTRAL& operator =(const SIESTA_SPECTRAL& spec);	// Operator =

	// Member-Functions
	void set(LA_STRING id0, int parity0, int ns0, int ntor, int mpol, Array<double,1>& S0);	// initialize all public members but spectral data
	void Vspline();											// prepare interpolation in s
	double Vsplint(double s, int n, int m, int par);				// evaluate spline at s for mode m, n, and parity par
	double Vsplint(double s, int n, int m, int par, double& dyds);	// same as above, but with first derivative
	double ev(double s, double u, double v);				// evaluate Fourier series at (s,u,v)
	double ev(double s, double u, double v, double& dyds, double& dydu, double& dydv);	// same as above but with derivatives

}; //end of class

//------------------------ End of Class SIESTA_SPECTRAL--------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
SIESTA_SPECTRAL::SIESTA_SPECTRAL()
{
ns = 0;
ntor = 0;
mpol = 0;
id = "None";
parity = 0;
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
SIESTA_SPECTRAL& SIESTA_SPECTRAL::operator =(const SIESTA_SPECTRAL& spec)
{
if (this == &spec) return(*this);	    // if: x=x
ds = spec.ds;
d2ymns.reference(spec.d2ymns.copy());
d2ymnc.reference(spec.d2ymnc.copy());

id = spec.id;
parity = spec.parity;
ns = spec.ns;
ntor = spec.ntor;
mpol = spec.mpol;

S.reference(spec.S);

ymns.reference(spec.ymns.copy());
ymnc.reference(spec.ymnc.copy());

return(*this);
}

//---------------------------- orthonorm ----------------------------------------------------------------------------------
double SIESTA_SPECTRAL::orthonorm(int n, int m)
{
if((n == 0) && (m == 0)) return 1;
else return sqrt2;
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------------------- set ----------------------------------------------------------------------------------------
// initalize all public & private member variables, except the spectral data and its derivative
void SIESTA_SPECTRAL::set(LA_STRING id0, int parity0, int ns0, int ntor0, int mpol0, Array<double,1>& S0)
{
id = id0;
parity = parity0;
ns = ns0;
ntor = ntor0;
mpol = mpol0;
S.reference(S0);
ds = 1.0 / (ns - 1);
}

//---------------------------- spline -------------------------------------------------------------------------------------
// prepare 1D spline interpolation of ymn(s)  ->  constructs d2ymn
void SIESTA_SPECTRAL::Vspline()
{
int n,m;
double d1,dn;
Array<double,1> slice, d2slice;
Range all = Range::all();

TinyVector <int,3> index3(1,-ntor,0);	// Array ranges

// cosine series or both
if(parity >= 0)
{
	d2ymnc.resize(ns, 2*ntor+1, mpol+1);	d2ymnc.reindexSelf(index3);
	for(n=-ntor;n<=ntor;n++)
	{
		for(m=0;m<=mpol;m++)
		{
			d1 = (ymnc(2,n,m) - ymnc(1,n,m)) / ds;
			dn = (ymnc(ns,n,m) - ymnc(ns-1,n,m)) / ds;
			slice.reference(ymnc(all,n,m));
			d2slice.reference(d2ymnc(all,n,m));
			spline(S, slice, ns, d1, dn, d2slice);
		}
	}
}

// sine series or both
if(parity <= 0)
{
	d2ymns.resize(ns, 2*ntor+1, mpol+1);	d2ymns.reindexSelf(index3);
	for(n=-ntor;n<=ntor;n++)
	{
		for(m=0;m<=mpol;m++)
		{
			d1 = (ymns(2,n,m) - ymns(1,n,m)) / ds;
			dn = (ymns(ns,n,m) - ymns(ns-1,n,m)) / ds;
			slice.reference(ymns(all,n,m));
			d2slice.reference(d2ymns(all,n,m));
			spline(S, slice, ns, d1, dn, d2slice);
		}
	}
}
}

//---------------------------- splint -------------------------------------------------------------------------------------
// evaluates spline at s
double SIESTA_SPECTRAL::Vsplint(double s, int n, int m, int par)
{
double dyds;
return Vsplint(s, n, m, par, dyds);
}

//-----------------------------------------------
// ... with first derivative
double SIESTA_SPECTRAL::Vsplint(double s, int n, int m, int par, double& dyds)
{
double y;
double orth = orthonorm(n,m);
Array<double,1> slice, d2slice;
Range all = Range::all();

// cosine series
if(par == 1)
{
	slice.reference(ymnc(all,n,m));
	d2slice.reference(d2ymnc(all,n,m));
	splint(S, slice, d2slice, ns, s, y, dyds);
}
else // sine series
{
	slice.reference(ymns(all,n,m));
	d2slice.reference(d2ymns(all,n,m));
	splint(S, slice, d2slice, ns, s, y, dyds);
}
dyds *= orth;
return y * orth;
}

//---------------------------- ev -----------------------------------------------------------------------------------------
// evaluates Fourier series at location (s,u,v)
// orthonorm is part of splint!
double SIESTA_SPECTRAL::ev(double s, double u, double v)
{
int n,m;
double sinuv, cosuv;
double spl;
double y = 0;

// cosine series
if(parity >= 0)
{
	for(n=-ntor;n<=ntor;n++)
	{
		for(m=0;m<=mpol;m++)
		{
			cosuv = cos(m*u + n*v);
			spl = Vsplint(s, n, m, 1);
			y += spl * cosuv;
		}
	}
}

// sine series
if(parity <= 0)
{
	for(n=-ntor;n<=ntor;n++)
	{
		for(m=0;m<=mpol;m++)
		{
			sinuv = sin(m*u + n*v);
			spl = Vsplint(s, n, m, -1);
			y += spl * sinuv;
		}
	}
}

return y;
}
//-----------------------------------------------
// ... and return all derivatives too
double SIESTA_SPECTRAL::ev(double s, double u, double v, double& dyds, double& dydu, double& dydv)
{
int m,n;
double sinuv, cosuv;
double spl, dsplds;
double y = 0;

dyds = 0;
dydu = 0;
dydv = 0;

// cosine series
if(parity == 1)
{
	for(n=-ntor;n<=ntor;n++)
	{
		for(m=0;m<=mpol;m++)
		{
			sinuv = sin(m*u + n*v);
			cosuv = cos(m*u + n*v);
			spl = Vsplint(s, n, m, 1, dsplds);
			y += spl * cosuv;
			dyds += dsplds * cosuv;
			dydu += -m * spl * sinuv;
			dydv += -n * spl * sinuv;
		}
	}
}

// sine series
if(parity == -1)
{
	for(n=-ntor;n<=ntor;n++)
	{
		for(m=0;m<=mpol;m++)
		{
			sinuv = sin(m*u + n*v);
			cosuv = cos(m*u + n*v);
			spl = Vsplint(s, n, m, -1, dsplds);
			y += spl * sinuv;
			dyds += dsplds * sinuv;
			dydu += m * spl * cosuv;
			dydv += n * spl * cosuv;
		}
	}
}

// both series
if(parity == 0)
{
	for(n=-ntor;n<=ntor;n++)
	{
		for(m=0;m<=mpol;m++)
		{
			sinuv = sin(m*u + n*v);
			cosuv = cos(m*u + n*v);

			spl = Vsplint(s, n, m, 1, dsplds);
			y += spl * cosuv;
			dyds += dsplds * cosuv;
			dydu += -m * spl * sinuv;
			dydv += -n * spl * sinuv;

			spl = Vsplint(s, n, m, -1, dsplds);
			y += spl * sinuv;
			dyds += dsplds * sinuv;
			dydu += m * spl * cosuv;
			dydv += n * spl * cosuv;
		}
	}
}
return y;
}
//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------------------
//--------- Begin Class SIESTA --------------------------------------------------------------------------------------------
class SIESTA
{
private:
// Parameter
	static const double smin = 0;
	static const double smax = 1;

// Member Variables
	int nrad;
	double b_factor;
	LA_STRING wout_file;
	bool lasym;

// Member-Functions
	void prepare_splines(void);		// prepare all 1D splines in s by calling spline from efit_class
	void get_RZ(double s, double u, double v, double& r, double& z, double& drds, double& drdu, double& drdv, double& dzds, double& dzdu, double& dzdv); // same as public, but with derivatives
	//int newton2D(double r, double phi, double z, double& s, double& u);	// 2D Newton procedure to find s,u from a given R,Z
	//double bisec(double r0, double phi, double z0, double u, double Raxis, double Zaxis, double a = 0, double b = 1);	// preconditioner to give a rough estimate of s to use as initial condition in newton2D

public: 
// Member Variables
	int ntor;
	int mpol;
	int ns;
	int nshalf;

	Array<double,1> S;
	Array<double,1> Shalf;

	SIESTA_SPECTRAL bsupsmn;
	SIESTA_SPECTRAL bsupumn;
	SIESTA_SPECTRAL bsupvmn;
	SIESTA_SPECTRAL presmn;

	VMEC wout;

// Constructors
	SIESTA();								// Default Constructor

// Member-Operators
	SIESTA& operator =(const SIESTA& SIES);	// Operator =

// Member-Functions
	void read(LA_STRING filename);			// read in B-field tracer file from SIESTA
	void get_B(double r, double phi, double z, double& br, double& bphi, double& bz);	// evaluate B at any location (r,phi,z) inside the SIESTA boundary
	void get_RZ(double s, double u, double v, double& r, double& z);	// evaluate R,Z at any location (s,u,v)
	//void get_su(double r, double phi, double z, double& s, double& u, double sstart = -1, double ustart = -1);	// find s,u at any location (r,phi,z) inside the SIESTA boundary
	void get_su(double r, double phi, double z, double& s, double& u);
}; //end of class

//------------------------ End of Class -----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
SIESTA::SIESTA()
{
nrad = 0;
b_factor = 0;
wout_file = " ";
lasym = false;
ntor = 0;
mpol = 0;
nrad = 1;
ns = nrad;
nshalf = ns-1;
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
SIESTA& SIESTA::operator =(const SIESTA& SIES)
{
if (this == &SIES) return(*this);	    // if: x=x

nrad = SIES.nrad;
b_factor = SIES.b_factor;
wout_file = SIES.wout_file;
lasym = SIES.lasym;

ntor = SIES.ntor;
mpol = SIES.mpol;
ns = SIES.ns;
nshalf = SIES.nshalf;

S.reference(SIES.S.copy());
Shalf.reference(SIES.Shalf.copy());

bsupsmn = SIES.bsupsmn;
bsupumn = SIES.bsupumn;
bsupvmn = SIES.bsupvmn;
presmn = SIES.presmn;

wout = SIES.wout;

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------------------- read ---------------------------------------------------------------------------------------
// read SIESTA .nc file and set member variables
void SIESTA::read(LA_STRING filename)
{
// Variables
int i;
int chk, ncid, varid;
Range all = Range::all();

// Input
i=filename.length();
if(filename.right(4) == ".dat") filename = filename(1,i-4) + ".nc";
chk = nc_open(filename, NC_NOWRITE, &ncid);
if(chk != 0) {cout << "Unable to open file " << filename << endl; EXIT;}

// Read the dimensions
chk = nc_inq_varid(ncid, "ntor", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &ntor);	// read
chk = nc_inq_varid(ncid, "mpol", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &mpol);	// read
chk = nc_inq_varid(ncid, "nrad", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &nrad);	// read
ns = nrad;
nshalf = nrad - 1;

// read and reindex data: (1 -> nshalf, -ntor -> ntor, 0 -> mpol)
TinyVector <int,3> index3(1,-ntor,0);	// Array range
Array<double,3> input;
input.resize(ns, 2*ntor+1, mpol+1);

chk = nc_inq_varid(ncid, "Bsups_r_m_n_", &varid);		// get variable id
chk = nc_get_var_double(ncid, varid, input.data());		// read
bsupsmn.ymns.resize(nshalf, 2*ntor+1, mpol+1);				// set size
bsupsmn.ymns.reindexSelf(index3);							// set indices
bsupsmn.ymns = input(Range(1,nshalf),all,all).copy();		// drop the first (index = 0) s value

chk = nc_inq_varid(ncid, "Bsupu_r_m_n_", &varid);
chk = nc_get_var_double(ncid, varid, input.data());
bsupumn.ymnc.resize(nshalf, 2*ntor+1, mpol+1);
bsupumn.ymnc.reindexSelf(index3);
bsupumn.ymnc = input(Range(1,nshalf),all,all).copy();

chk = nc_inq_varid(ncid, "Bsupv_r_m_n_", &varid);
chk = nc_get_var_double(ncid, varid, input.data());
bsupvmn.ymnc.resize(nshalf, 2*ntor+1, mpol+1);
bsupvmn.ymnc.reindexSelf(index3);
bsupvmn.ymnc = input(Range(1,nshalf),all,all).copy();

chk = nc_inq_varid(ncid, "pres_r_m_n_", &varid);
chk = nc_get_var_double(ncid, varid, input.data());
presmn.ymnc.resize(nshalf, 2*ntor+1, mpol+1);
presmn.ymnc.reindexSelf(index3);
presmn.ymnc = input(Range(1,nshalf),all,all).copy();

// wout file name
char a[500];
chk = nc_inq_varid(ncid, "wout_file", &varid);		// get variable id
chk = nc_get_var_text(ncid, varid, &a[0]);			// read
wout_file = a;
int idx = wout_file.indexOf(".nc");
wout_file = wout_file.left(idx+2);

// close file
chk = nc_close(ncid);

//---- prepare stuff ----------------------

// set S arrays
TinyVector <int,1> index(1);	// Array range
S.resize(ns); 			S.reindexSelf(index);
Shalf.resize(nshalf);	Shalf.reindexSelf(index);
double ds = (smax - smin) / (ns - 1);
for(i=1;i<=ns;i++) S(i) = smin + (i-1)*ds;
for(i=1;i<=nshalf;i++) Shalf(i) = S(i) + ds/2;

// read VMEC file
wout.read(wout_file);

// get normalization factor for B-field
b_factor = sqrt(1.0/fabs(wout.wb)) / pi2;

// up-down symmetry? lasym = False <-> yes;	lasym = True <-> no
lasym = wout.lasym;

prepare_splines();
}

//-------------------------------- get_B ----------------------------------------------------------------------------------
// evaluate br, bphi, bz at any arbitrary location (r, phi, z) inside the SIESTA boundary; r is the major radius here!
void SIESTA::get_B(double r, double phi, double z, double& br, double& bphi, double& bz)
{
double s, u, v;
double bsups;
double bsupu,bsupv;
double R,Z,g,gb;
double dRds, dRdu, dRdv, dZds, dZdu, dZdv;

phi = modulo2pi(phi);	// make sure phi in [0, 2pi]
v = phi;

get_su(r, phi, z, s, u);
bsups = bsupsmn.ev(s, u, v);
bsupu = bsupumn.ev(s, u, v);
bsupv = bsupvmn.ev(s, u, v);

get_RZ(s, u, v, R, Z, dRds, dRdu, dRdv, dZds, dZdu, dZdv);
g = wout.gmn.ev(s*s, u, v) * 2*s;

gb = g*b_factor;
bsups /= gb;
bsupu /= gb;
bsupv /= gb;

br = dRdu*bsupu + dRdv*bsupv + dRds*bsups;
bz = dZdu*bsupu + dZdv*bsupv + dZds*bsups;
bphi = R * bsupv;
}

//------------------------- get_RZ ---------------------------------------------------------------------------------------
// evaluates r,z at loaction (s,u,v)
void SIESTA::get_RZ(double s, double u, double v, double& r, double& z)
{
r = wout.rmn.ev(s*s, u, v);
z = wout.zmn.ev(s*s, u, v);
}

////------------------------- get_su ---------------------------------------------------------------------------------------
//// find s,u at any location (r,phi,z) inside the SIESTA boundary
//// sstart and ustart are optional initial guesses of s and u
//void SIESTA::get_su(double r, double phi, double z, double& s, double& u, double sstart, double ustart)
//{
//double err, Raxis, Zaxis;
//
//if((ustart == -1) || (sstart == -1)) wout.get_axis(phi, Raxis, Zaxis);
//
//// set initial guesses
//if(ustart == -1) u = polar_phi(r - Raxis, z - Zaxis);
//else u = ustart;
//if(sstart == -1) s = bisec(r, phi, z, u, Raxis, Zaxis);
//else s = sstart;
//
//// run Newton
//u = modulo2pi(u);	// make sure u in [0, 2pi]
//err = newton2D(r,phi,z,s,u);
//if(err > 0) cout << "SIESTA Newton2D: no convergence; remaining error: " << err << endl;
//
//}

//------------------------- get_su ---------------------------------------------------------------------------------------
// find s,u at any location (r,phi,z) inside the SIESTA boundary
// sstart and ustart are optional initial guesses of s and u
void SIESTA::get_su(double r, double phi, double z, double& s, double& u)
{
wout.get_su(r, phi, z, s, u);
int imax = 12;
while(s < 0 || s > 1)
{
	wout.get_su(r, phi, z, s, u, -1, -1, imax);
	imax += 2;
}
s = sqrt(s);
}

//----------------------- End of Public Member Functions ------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- Private Member Functions ------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//------------------------- get_RZ ----------------------------------------------------------------------------------------
// evaluates r,z at loaction (s,u,v) with derivatives
void SIESTA::get_RZ(double s, double u, double v, double& r, double& z, double& drds, double& drdu, double& drdv, double& dzds, double& dzdu, double& dzdv)
{
r = wout.rmn.ev(s*s, u, v, drds, drdu, drdv);
z = wout.zmn.ev(s*s, u, v, dzds, dzdu, dzdv);
drds *= 2*s;
dzds *= 2*s;
}

//--------------------- prepare_splines -----------------------------------------------------------------------------------
// get the 2. derivative in s-direction, which is input to splint for all modes n, m
void SIESTA::prepare_splines(void)
{
if(lasym)
{
	bsupsmn.set("bsupsmn", 0, nshalf, ntor, mpol, Shalf);
	bsupumn.set("bsupumn", 0, nshalf, ntor, mpol, Shalf);
	bsupvmn.set("bsupvmn", 0, nshalf, ntor, mpol, Shalf);
	presmn.set("presmn", 0, nshalf, ntor, mpol, Shalf);
}
else
{
	bsupsmn.set("bsupsmns", -1, nshalf, ntor, mpol, Shalf);
	bsupumn.set("bsupumnc", 1, nshalf, ntor, mpol, Shalf);
	bsupvmn.set("bsupvmnc", 1, nshalf, ntor, mpol, Shalf);
	presmn.set("presmnc", 1, nshalf, ntor, mpol, Shalf);
}
bsupsmn.Vspline();
bsupumn.Vspline();
bsupvmn.Vspline();
presmn.Vspline();
}

////------------------------ newton2D ----------------------------------------------------------------------------------------
//// 2D Newton procedure
//int SIESTA::newton2D(double r0, double phi0, double z0, double& s, double& u)
//{
//const int imax = 20;
//const double delta = 1e-12;
//
//int i;
//double r,z,drds,drdu,drdv,dzds,dzdu,dzdv;
//double det,delta_s,delta_u,fr,fz,err;
//
//for(i=0;i<imax;i++)
//{
//	get_RZ(s, u, phi0, r, z, drds, drdu, drdv, dzds, dzdu, dzdv);
//
//	fr = r - r0;
//	fz = z - z0;
//
//	det = drds * dzdu - drdu * dzds;
//	delta_s = (dzdu * fr - drdu * fz) / det;
//	delta_u = (drds * fz - dzds * fr) / det;
//
//	err = sqrt(delta_s*delta_s + delta_u*delta_u);
//	if(err < delta) return 0;	// convergence
//
//	s -= delta_s;
//	u -= delta_u;
//	u = modulo2pi(u);	// make sure u in [0, 2pi]
//}
//return err;	// no convergence
//}
//
////------------------------ bisec ----------------------------------------------------------------------------------------
//// bisection with only imax steps to find a crude estimate of s; use as preconditioner for newton2D
//// set u fixed as the geometric angle
//// assume f(a) < 0 & f(b) > 0 with f = (r(s,u,v) - Raxis)^2 + (z(s,u,v) - Zaxis)^2 - (r0 - Raxis)^2 - (z0 - Zaxis)^2
//// and v = phi0
//double SIESTA::bisec(double r0, double phi0, double z0, double u, double Raxis, double Zaxis, double a, double b)
//{
//int i;
//double s,rminor,r,z,f;
//
//// Iterations
//const int imax = 10;
//
//// reference point
//const double rminor0 = (r0 - Raxis)*(r0 - Raxis) + (z0 - Zaxis)*(z0 - Zaxis);
//
//for(i=1;i<=imax;i++)
//{
//	s = (a + b)/2.0;
//	get_RZ(s, u, phi0, r, z);
//	r -= Raxis;
//	z -= Zaxis;
//	rminor = r*r + z*z;
//	f = rminor - rminor0;
//	if(f > 0) b = s;
//	else a = s;
//}
//return s;
//}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  SIESTA_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

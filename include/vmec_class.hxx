// VMEC Class reads some VMEC data
// used by SIESTA class for now
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						4.08.14


// Define
//--------
#ifndef VMEC_CLASS_INCLUDED
#define VMEC_CLASS_INCLUDED

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
using namespace blitz;

// Prototypes
//-----------
void spline(Array<double,1>& x, Array<double,1>& y, int n, double yp1, double ypn, Array<double,1>& y2);
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, int n, double x, double& y, double& yx);

// Golbal Parameters
//------------------
extern ofstream ofs2;


//--------- Begin Class VMEC_SPECTRAL -------------------------------------------------------------------------------------
class VMEC_SPECTRAL
{
private:
	// Member Variables
	double ds;		// s array step size
	Array<double,2> d2ymns;	// d^2/ds^2 ymns(s), output of spline and input of splint
	Array<double,2> d2ymnc;	// d^2/ds^2 ymnc(s), output of spline and input of splint

public:
	// Member Variables
	LA_STRING id; 	// name of spectral variable

	int parity;		// parity = 1 <-> cosine,	parity = -1 <-> sine, 	parity = 0 <-> both
	int ns;			// number of s points
	int mnmax;		// number of m,n modes

	Array<int,1> xn;	// all n modes
	Array<int,1> xm;	// all m modes
	Array<double,1> S;	// s array

	Array<double,2> ymns; 	// sine spectral data
	Array<double,2> ymnc; 	// cosine spectral data

	// Constructors
	VMEC_SPECTRAL();								// Default Constructor

	// Member-Operators
	VMEC_SPECTRAL& operator =(const VMEC_SPECTRAL& spec);	// Operator =

	// Member-Functions
	void set(LA_STRING id0, int parity0, int ns0, int mnmax0, Array<int,1>& xn0, Array<int,1>& xm0, Array<double,1>& S0);	// initialize all public members but spectral data
	void Vspline();											// prepare interpolation in s
	double Vsplint(double s, int i, int par);				// evaluate spline at s for m,n mode i, and parity par
	double Vsplint(double s, int i, int par, double& dyds);	// same as above, but with first derivative
	double ev(double s, double u, double v);				// evaluate Fourier series at (s,u,v)
	double ev(double s, double u, double v, double& dyds, double& dydu, double& dydv);	// same as above but with derivatives

}; //end of class

//------------------------ End of Class VMEC_SPECTRAL----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
VMEC_SPECTRAL::VMEC_SPECTRAL()
{
ns = 0;
id = "None";
parity = 0;
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
VMEC_SPECTRAL& VMEC_SPECTRAL::operator =(const VMEC_SPECTRAL& spec)
{
if (this == &spec) return(*this);	    // if: x=x
ds = spec.ds;
d2ymns.reference(spec.d2ymns.copy());
d2ymnc.reference(spec.d2ymnc.copy());

id = spec.id;
parity = spec.parity;
ns = spec.ns;
mnmax = spec.mnmax;

xn.reference(spec.xn);
xm.reference(spec.xm);
S.reference(spec.S);

ymns.reference(spec.ymns.copy());
ymnc.reference(spec.ymnc.copy());

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------------------- set ----------------------------------------------------------------------------------------
// initalize all public & private member variables, except the spectral data and its derivative
void VMEC_SPECTRAL::set(LA_STRING id0, int parity0, int ns0, int mnmax0, Array<int,1>& xn0, Array<int,1>& xm0, Array<double,1>& S0)
{
id = id0;
parity = parity0;
ns = ns0;
mnmax = mnmax0;
xn.reference(xn0);
xm.reference(xm0);
S.reference(S0);
ds = 1.0 / (ns - 1);
}

//---------------------------- spline -------------------------------------------------------------------------------------
// prepare 1D spline interpolation of ymn(s)  ->  constructs d2ymn
void VMEC_SPECTRAL::Vspline()
{
int i;
double d1,dn;
Array<double,1> slice, d2slice;
Range all = Range::all();

TinyVector <int,2> index2(1,0);	// Array ranges

// cosine series or both
if(parity >= 0)
{
	d2ymnc.resize(ns,mnmax);	d2ymnc.reindexSelf(index2);
	for(i=0;i<mnmax;i++)
	{
		d1 = (ymnc(2,i) - ymnc(1,i)) / ds;
		dn = (ymnc(ns,i) - ymnc(ns-1,i)) / ds;
		slice.reference(ymnc(all,i));
		d2slice.reference(d2ymnc(all,i));
		spline(S, slice, ns, d1, dn, d2slice);
	}
}

// sine series or both
if(parity <= 0)
{
	d2ymns.resize(ns,mnmax);	d2ymns.reindexSelf(index2);
	for(i=0;i<mnmax;i++)
	{
		d1 = (ymns(2,i) - ymns(1,i)) / ds;
		dn = (ymns(ns,i) - ymns(ns-1,i)) / ds;
		slice.reference(ymns(all,i));
		d2slice.reference(d2ymns(all,i));
		spline(S, slice, ns, d1, dn, d2slice);
	}
}
}

//---------------------------- splint -------------------------------------------------------------------------------------
// evaluates spline at s
double VMEC_SPECTRAL::Vsplint(double s, int i, int par)
{
double dyds;
return Vsplint(s, i, par, dyds);
}

//-----------------------------------------------
// ... with first derivative
double VMEC_SPECTRAL::Vsplint(double s, int i, int par, double& dyds)
{
double y;
Array<double,1> slice, d2slice;
Range all = Range::all();

// cosine series
if(par == 1)
{
	slice.reference(ymnc(all,i));
	d2slice.reference(d2ymnc(all,i));
	splint(S, slice, d2slice, ns, s, y, dyds);
}
else // sine series
{
	slice.reference(ymns(all,i));
	d2slice.reference(d2ymns(all,i));
	splint(S, slice, d2slice, ns, s, y, dyds);
}
return y;
}

//---------------------------- ev -----------------------------------------------------------------------------------------
// evaluates Fourier series at location (s,u,v)
double VMEC_SPECTRAL::ev(double s, double u, double v)
{
int i;
double sinuv, cosuv;
double spl;
double y = 0;

// cosine series
if(parity >= 0)
{
	for(i=0;i<mnmax;i++)
	{
		cosuv = cos(xm(i)*u - xn(i)*v);
		spl = Vsplint(s, i, 1);
		y += spl * cosuv;
	}
}

// sine series
if(parity <= 0)
{
	for(i=0;i<mnmax;i++)
	{
		sinuv = sin(xm(i)*u - xn(i)*v);
		spl = Vsplint(s, i, -1);
		y += spl * sinuv;
	}
}

return y;
}
//-----------------------------------------------
// ... and return all derivatives too
double VMEC_SPECTRAL::ev(double s, double u, double v, double& dyds, double& dydu, double& dydv)
{
int i,m,n;
double sinuv, cosuv;
double spl, dsplds;
double y = 0;

dyds = 0;
dydu = 0;
dydv = 0;

// cosine series
if(parity == 1)
{
	for(i=0;i<mnmax;i++)
	{
		m = xm(i);
		n = xn(i);
		sinuv = sin(m*u - n*v);
		cosuv = cos(m*u - n*v);
		spl = Vsplint(s, i, 1, dsplds);
		y += spl * cosuv;
		dyds += dsplds * cosuv;
		dydu += -m * spl * sinuv;
		dydv += n * spl * sinuv;
	}
}

// sine series
if(parity == -1)
{
	for(i=0;i<mnmax;i++)
	{
		m = xm(i);
		n = xn(i);
		sinuv = sin(m*u - n*v);
		cosuv = cos(m*u - n*v);
		spl = Vsplint(s, i, -1, dsplds);
		y += spl * sinuv;
		dyds += dsplds * sinuv;
		dydu += m * spl * cosuv;
		dydv += -n * spl * cosuv;
	}
}

// both series
if(parity == 0)
{
	for(i=0;i<mnmax;i++)
	{
		m = xm(i);
		n = xn(i);
		sinuv = sin(m*u - n*v);
		cosuv = cos(m*u - n*v);

		spl = Vsplint(s, i, 1, dsplds);
		y += spl * cosuv;
		dyds += dsplds * cosuv;
		dydu += -m * spl * sinuv;
		dydv += n * spl * sinuv;

		spl = Vsplint(s, i, -1, dsplds);
		y += spl * sinuv;
		dyds += dsplds * sinuv;
		dydu += m * spl * cosuv;
		dydv += -n * spl * cosuv;
	}
}
return y;
}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------------------
//--------- Begin Class VMEC ----------------------------------------------------------------------------------------------
class VMEC
{
private:
// Member Variables
	int ntor;

	Array<int,1> xn;
	Array<int,1> xm;

// Member-Functions
	void prepare_splines(void);		// prepare all 1D splines in s by calling spline from efit_class

public: 
// Member Variables
	int ns;
	int nshalf;
	int mnmax;
	double wb;
	bool lasym;

	Array<double,1> S;
	Array<double,1> Shalf;

	VMEC_SPECTRAL rmn;
	VMEC_SPECTRAL zmn;
	VMEC_SPECTRAL gmn;

	Array<double,1> raxis_cc;
	Array<double,1> raxis_cs;
	Array<double,1> zaxis_cc;
	Array<double,1> zaxis_cs;

// Constructors
	VMEC();								// Default Constructor

// Member-Operators
	VMEC& operator =(const VMEC& V);	// Operator =

// Member-Functions
	void read(LA_STRING filename);			// read in wout file from VMEC
	void get_axis(double v, double& Raxis, double& Zaxis);	// get magnetic axis

}; //end of class

//------------------------ End of Class VMEC-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
VMEC::VMEC()
{
ntor = 0;
ns = 0;
nshalf = 0;
mnmax = 0;
wb = 0;
lasym = false;
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
VMEC& VMEC::operator =(const VMEC& V)
{
if (this == &V) return(*this);	    // if: x=x
ntor = V.ntor;

xn.reference(V.xn.copy());
xm.reference(V.xm.copy());

ns = V.ns;
mnmax = V.mnmax;
wb = V.wb;
lasym = V.lasym;

S.reference(V.S.copy());
Shalf.reference(V.Shalf.copy());

rmn = V.rmn;
zmn = V.zmn;
gmn = V.gmn;

raxis_cc.reference(V.raxis_cc.copy());
raxis_cs.reference(V.raxis_cs.copy());
zaxis_cc.reference(V.zaxis_cc.copy());
zaxis_cs.reference(V.zaxis_cs.copy());

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------------------- read ---------------------------------------------------------------------------------------
// read VMEC wout.nc file and set member variables
void VMEC::read(LA_STRING filename)
{
// Variables
int i;
int chk, ncid, varid;
Range all = Range::all();

//---- Read VMEC file ----------------------
// Input
chk = nc_open(filename, NC_NOWRITE, &ncid);
if(chk != 0) {cout << "Unable to open file " << filename << endl; EXIT;}

// Read the dimensions
chk = nc_inq_varid(ncid, "ns", &varid);		// get variable id
chk = nc_get_var_int(ncid, varid, &ns);		// read
chk = nc_inq_varid(ncid, "ntor", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &ntor);	// read
chk = nc_inq_varid(ncid, "mnmax", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &mnmax);	// read
nshalf = ns - 1;

// read lasym
int lasym_in;
chk = nc_inq_varid(ncid, "lasym__logical__", &varid);		// get variable id
chk = nc_get_var_int(ncid, varid, &lasym_in);	// read
lasym = bool(lasym_in);

// read xn & xm
xn.resize(mnmax);
chk = nc_inq_varid(ncid, "xn", &varid);
chk = nc_get_var_int(ncid, varid, xn.data());
xm.resize(mnmax);
chk = nc_inq_varid(ncid, "xm", &varid);
chk = nc_get_var_int(ncid, varid, xm.data());

// read and reindex spectral data: full-mesh:(1 -> ns, 0 -> mnmax-1),  half-mesh:(1 -> nshalf, 0 -> mnmax-1)
TinyVector <int,2> index2(1,0);	// Array range
Array<double,2> input;
input.resize(ns, mnmax);

if(lasym)
{
	chk = nc_inq_varid(ncid, "rmns", &varid);			// get variable id
	chk = nc_get_var_double(ncid, varid, input.data());	// read
	rmn.ymns.resize(ns, mnmax); 						// set size
	rmn.ymns.reindexSelf(index2);						// set indices
	rmn.ymns = input.copy();							// move into place

	chk = nc_inq_varid(ncid, "zmnc", &varid);			// get variable id
	chk = nc_get_var_double(ncid, varid, input.data());	// read
	zmn.ymnc.resize(ns, mnmax); 						// set size
	zmn.ymnc.reindexSelf(index2);						// set indices
	zmn.ymnc = input.copy();							// move into place

	chk = nc_inq_varid(ncid, "gmns", &varid);			// get variable id
	chk = nc_get_var_double(ncid, varid, input.data());	// read
	gmn.ymns.resize(nshalf, mnmax); 						// set size
	gmn.ymns.reindexSelf(index2);							// set indices
	gmn.ymns = input(Range(1,nshalf),all).copy();			// drop the first (index = 0) s value
}

chk = nc_inq_varid(ncid, "rmnc", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, input.data());	// read
rmn.ymnc.resize(ns, mnmax); 						// set size
rmn.ymnc.reindexSelf(index2);							// set indices
rmn.ymnc = input.copy();								// move into place

chk = nc_inq_varid(ncid, "zmns", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, input.data());	// read
zmn.ymns.resize(ns, mnmax); 						// set size
zmn.ymns.reindexSelf(index2);							// set indices
zmn.ymns = input.copy();								// move into place

chk = nc_inq_varid(ncid, "gmnc", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, input.data());	// read
gmn.ymnc.resize(nshalf, mnmax); 						// set size
gmn.ymnc.reindexSelf(index2);							// set indices
gmn.ymnc = input(Range(1,nshalf),all).copy();			// drop the first (index = 0) s value

// read parameter
chk = nc_inq_varid(ncid, "wb", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &wb);		// read

// read axis
if(lasym)
{
	raxis_cs.resize(ntor+1);
	chk = nc_inq_varid(ncid, "raxis_cs", &varid);
	chk = nc_get_var_double(ncid, varid, raxis_cs.data());
	zaxis_cc.resize(ntor+1);
	chk = nc_inq_varid(ncid, "zaxis_cc", &varid);
	chk = nc_get_var_double(ncid, varid, zaxis_cc.data());
	zaxis_cs.resize(ntor+1);
	chk = nc_inq_varid(ncid, "zaxis_cs", &varid);
	chk = nc_get_var_double(ncid, varid, zaxis_cs.data());
}
raxis_cc.resize(ntor+1);
chk = nc_inq_varid(ncid, "raxis_cc", &varid);
chk = nc_get_var_double(ncid, varid, raxis_cc.data());

// close file
chk = nc_close(ncid);

//---- prepare stuff ----------------------

// set S arrays
TinyVector <int,1> index(1);	// Array range
S.resize(ns); 			S.reindexSelf(index);
Shalf.resize(nshalf);	Shalf.reindexSelf(index);
double ds = 1.0 / (ns - 1);
for(i=1;i<=ns;i++) S(i) = (i-1)*ds;
for(i=1;i<=nshalf;i++) Shalf(i) = S(i) + ds/2;

prepare_splines();
}

//------------------------ get_axis ----------------------------------------------------------------------------------------
// evaluate the Fourier series of the magnetic axis at location v
void VMEC::get_axis(double v, double& Raxis, double& Zaxis)
{
int n;
double sinnv,cosnv;

Raxis = 0;
Zaxis = 0;

for(n=0;n<=ntor;n++)
{
	cosnv = cos(n*v);
	Raxis += raxis_cc(n) * cosnv;
	if(lasym)
	{
		sinnv = sin(n*v);
		Raxis += raxis_cs(n) * sinnv;
		Zaxis += zaxis_cs(n) * sinnv + zaxis_cc(n) * cosnv;
	}
}
}

//----------------------- End of Public Member Functions ------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- Private Member Functions ------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------------------- prepare_splines -----------------------------------------------------------------------------------
// get the 2. derivative in s-direction, which is input to splint for all modes n, m
void VMEC::prepare_splines(void)
{
if(lasym)
{
	rmn.set("rmn", 0, ns, mnmax, xn, xm, S);
	zmn.set("zmn", 0, ns, mnmax, xn, xm, S);
	gmn.set("gmn", 0, nshalf, mnmax, xn, xm, Shalf);
}
else
{
	rmn.set("rmnc", 1, ns, mnmax, xn, xm, S);
	zmn.set("zmns", -1, ns, mnmax, xn, xm, S);
	gmn.set("gmnc", 1, nshalf, mnmax, xn, xm, Shalf);
}

rmn.Vspline();
zmn.Vspline();
gmn.Vspline();
}



//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  VMEC_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

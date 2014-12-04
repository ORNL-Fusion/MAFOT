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
	double ev(double s, double u, double v, double& dyds, double& dydu, double& dydv);	// same as above but with 1st derivatives
	double ev(double s, double u, double v, double& dyds, double& dydu, double& dydv, double& dydudv); // same as above but with mixed derivative dudv

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
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
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
VMEC_SPECTRAL& VMEC_SPECTRAL::operator =(const VMEC_SPECTRAL& spec)
{
if (this == &spec) return(*this);	    // if: x=x
ds = spec.ds;
d2ymns.reference(spec.d2ymns);
d2ymnc.reference(spec.d2ymnc);

id = spec.id;
parity = spec.parity;
ns = spec.ns;
mnmax = spec.mnmax;

xn.reference(spec.xn);
xm.reference(spec.xm);
S.reference(spec.S);

ymns.reference(spec.ymns);
ymnc.reference(spec.ymnc);

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
ds = (S(ns) - S(1)) / (ns - 1);
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
		d1 = (-49/20.0*ymnc(1,i) + 6*ymnc(2,i) - 15/2.0*ymnc(3,i) + 20/3.0*ymnc(4,i) - 15/4.0*ymnc(5,i) + 6/5.0*ymnc(6,i) - 1/6.0*ymnc(7,i)) / ds;
		dn = (49/20.0*ymnc(ns,i) - 6*ymnc(ns-1,i) + 15/2.0*ymnc(ns-2,i) - 20/3.0*ymnc(ns-3,i) + 15/4.0*ymnc(ns-4,i) - 6/5.0*ymnc(ns-5,i) + 1/6.0*ymnc(ns-6,i)) / ds;
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
		d1 = (-49/20.0*ymns(1,i) + 6*ymns(2,i) - 15/2.0*ymns(3,i) + 20/3.0*ymns(4,i) - 15/4.0*ymns(5,i) + 6/5.0*ymns(6,i) - 1/6.0*ymns(7,i)) / ds;
		dn = (49/20.0*ymns(ns,i) - 6*ymns(ns-1,i) + 15/2.0*ymns(ns-2,i) - 20/3.0*ymns(ns-3,i) + 15/4.0*ymns(ns-4,i) - 6/5.0*ymns(ns-5,i) + 1/6.0*ymns(ns-6,i)) / ds;
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
		if(s == 0 && xm(i) > 0) continue;	// no m modes on magnetic axis
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
		if(s == 0 && xm(i) > 0) continue;	// no m modes on magnetic axis
		sinuv = sin(xm(i)*u - xn(i)*v);
		spl = Vsplint(s, i, -1);
		y += spl * sinuv;
	}
}

return y;
}
//-----------------------------------------------
// ... and return all 1st derivatives
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
		if(s == 0 && m > 0) continue;	// no m modes on magnetic axis
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
		if(s == 0 && m > 0) continue;	// no m modes on magnetic axis
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
		if(s == 0 && m > 0) continue;	// no m modes on magnetic axis
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

//-----------------------------------------------
// ... ,return all 1st derivatives and the mixed derivative dudv
double VMEC_SPECTRAL::ev(double s, double u, double v, double& dyds, double& dydu, double& dydv, double& dydudv)
{
int i,m,n;
double sinuv, cosuv;
double spl, dsplds;
double y = 0;

dyds = 0;
dydu = 0;
dydv = 0;
dydudv = 0;

// cosine series
if(parity == 1)
{
	for(i=0;i<mnmax;i++)
	{
		m = xm(i);
		n = xn(i);
		if(s == 0 && m > 0) continue;	// no m modes on magnetic axis
		sinuv = sin(m*u - n*v);
		cosuv = cos(m*u - n*v);
		spl = Vsplint(s, i, 1, dsplds);
		y += spl * cosuv;
		dyds += dsplds * cosuv;
		dydu += -m * spl * sinuv;
		dydv += n * spl * sinuv;
		dydudv += n*m * spl * cosuv;
	}
}

// sine series
if(parity == -1)
{
	for(i=0;i<mnmax;i++)
	{
		m = xm(i);
		n = xn(i);
		if(s == 0 && m > 0) continue;	// no m modes on magnetic axis
		sinuv = sin(m*u - n*v);
		cosuv = cos(m*u - n*v);
		spl = Vsplint(s, i, -1, dsplds);
		y += spl * sinuv;
		dyds += dsplds * sinuv;
		dydu += m * spl * cosuv;
		dydv += -n * spl * cosuv;
		dydudv += n*m * spl * sinuv;
	}
}

// both series
if(parity == 0)
{
	for(i=0;i<mnmax;i++)
	{
		m = xm(i);
		n = xn(i);
		if(s == 0 && m > 0) continue;	// no m modes on magnetic axis
		sinuv = sin(m*u - n*v);
		cosuv = cos(m*u - n*v);

		spl = Vsplint(s, i, 1, dsplds);
		y += spl * cosuv;
		dyds += dsplds * cosuv;
		dydu += -m * spl * sinuv;
		dydv += n * spl * sinuv;
		dydudv += n*m * spl * cosuv;

		spl = Vsplint(s, i, -1, dsplds);
		y += spl * sinuv;
		dyds += dsplds * sinuv;
		dydu += m * spl * cosuv;
		dydv += -n * spl * cosuv;
		dydudv += n*m * spl * sinuv;
	}
}
return y;
}

//------------------------ End of Class VMEC_SPECTRAL----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------------------
//--------- Begin Class VMEC ----------------------------------------------------------------------------------------------
class VMEC
{
private:
// Member Variables
	int ntor;
	int ns;
	int nshalf;
	int mnmax;
	int mnmaxpot;

	Array<int,1> xnpot;
	Array<int,1> xmpot;
	Array<int,1> xn;
	Array<int,1> xm;

	Array<double,1> raxis_cc;
	Array<double,1> raxis_cs;
	Array<double,1> zaxis_cc;
	Array<double,1> zaxis_cs;

	Array<double,1> potsin;
	Array<double,1> potcos;

// Member-Functions
	void prepare_splines(void);		// prepare all 1D splines in s by calling spline from efit_class
	int newton2D(double r, double phi, double z, double& s, double& u);	// 2D Newton procedure to find s,u from a given R,Z
	double bisec(double r0, double phi, double z0, double u, double Raxis, double Zaxis, double a = 0, double b = 1, int imax = 10);	// preconditioner to give a rough estimate of s to use as initial condition in newton2D

public: 
// Member Variables
	int nextcur;

	double wb;
	double ctor;

	bool lasym;
	bool lpot;		// potential spectral data available or not

	LA_STRING input_extension;
	LA_STRING mgrid_file;

	Array<double,1> S;
	Array<double,1> Shalf;
	Array<double,1> extcur;

	VMEC_SPECTRAL rmn;
	VMEC_SPECTRAL zmn;
	VMEC_SPECTRAL gmn;
	VMEC_SPECTRAL bsupumn;
	VMEC_SPECTRAL bsupvmn;

// Constructors
	VMEC();								// Default Constructor
	VMEC(const VMEC& wout);				// Copy Constructor
	VMEC(LA_STRING filename);			// Standard Constructor

// Member-Operators
	VMEC& operator =(const VMEC& V);	// Operator =

// Member-Functions
	void read(LA_STRING filename);			// read in wout file from VMEC
	void get_axis(double v, double& Raxis, double& Zaxis);	// get magnetic axis
	void get_B2D(double s, double u, double v, double& BR, double& Bphi, double& BZ);	// get B-field at (s,u,v)
	double pot(double u, double v);			// magnetic scalar potental on s = 1 surface
	double pot(double u, double v, double& dpdu, double& dpdv, double& dpdudv);	// ...with 1st and mixed derivatives
	double get_r(double R, double Z, double Raxis, double Zaxis);
	double get_theta(double R, double Z, double Raxis, double Zaxis);
	void get_su(double R, double phi, double Z, double& s, double& u, double sstart = -1, double ustart = -1, int imax = 10);	// find s,u at any location (R,phi,Z) inside the VMEC boundary

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
VMEC::VMEC()
{
ntor = 0;
ns = 0;
nshalf = 0;
mnmax = 0;
mnmaxpot = 0;
nextcur = 0;

wb = 0;
ctor = 0;

lasym = false;
lpot = false;

input_extension = "None";
mgrid_file = "None";
}

// Copy Constructor
VMEC::VMEC(const VMEC& wout)
{
*this = wout;
}

// Standard Constructor
VMEC::VMEC(LA_STRING filename)
{
read(filename);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
VMEC& VMEC::operator =(const VMEC& V)
{
if (this == &V) return(*this);	    // if: x=x
ntor = V.ntor;
ns = V.ns;
nshalf = V.nshalf;
mnmax = V.mnmax;
mnmaxpot = V.mnmaxpot;
nextcur = V.nextcur;

wb = V.wb;
ctor = V.ctor;

lasym = V.lasym;
lpot = V.lpot;

input_extension = V.input_extension;
mgrid_file = V.mgrid_file;

xn.reference(V.xn);
xm.reference(V.xm);
xnpot.reference(V.xnpot);
xmpot.reference(V.xmpot);

S.reference(V.S);
Shalf.reference(V.Shalf);
extcur.reference(V.extcur);

rmn = V.rmn;
zmn = V.zmn;
gmn = V.gmn;
bsupumn = V.bsupumn;
bsupvmn = V.bsupvmn;

raxis_cc.reference(V.raxis_cc);
raxis_cs.reference(V.raxis_cs);
zaxis_cc.reference(V.zaxis_cc);
zaxis_cs.reference(V.zaxis_cs);

potsin.reference(V.potsin);
potcos.reference(V.potcos);

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
nshalf = ns;

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

// set sign array for s < 0 entry in half-grid spectral variables
Array<int,1> signhalf(mnmax);
for(i=0;i<mnmax;i++) 	// signhalf = (-1)^xm
{
	if(xm(i)%2 == 0) signhalf(i) = 1;
	else signhalf(i) = -1;
}

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
	gmn.ymns = input.copy();							// move into place
	gmn.ymns(1,all) = gmn.ymns(2,all) * signhalf;	// expand beyond magnetic axis: Shalf(1) = -Shalf(2) => Amn(-s)*sin(mu-nv) = Amn(s)*sin(m(u+pi)-nv) = (-1)^m * Amn(s)*sin(mu-nv); same for cos

	chk = nc_inq_varid(ncid, "bsupumns", &varid);			// get variable id
	chk = nc_get_var_double(ncid, varid, input.data());	// read
	bsupumn.ymns.resize(nshalf, mnmax); 						// set size
	bsupumn.ymns.reindexSelf(index2);							// set indices^n
	bsupumn.ymns = input.copy();							// move into place
	bsupumn.ymns(1,all) = bsupumn.ymns(2,all) * signhalf;

	chk = nc_inq_varid(ncid, "bsupvmns", &varid);			// get variable id
	chk = nc_get_var_double(ncid, varid, input.data());	// read
	bsupvmn.ymns.resize(nshalf, mnmax); 						// set size
	bsupvmn.ymns.reindexSelf(index2);							// set indices
	bsupvmn.ymns = input.copy();							// move into place
	bsupvmn.ymns(1,all) = bsupvmn.ymns(2,all) * signhalf;
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
gmn.ymnc = input.copy();							// move into place
gmn.ymnc(1,all) = gmn.ymnc(2,all) * signhalf;

chk = nc_inq_varid(ncid, "bsupumnc", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, input.data());	// read
bsupumn.ymnc.resize(nshalf, mnmax); 						// set size
bsupumn.ymnc.reindexSelf(index2);							// set indices
bsupumn.ymnc = input.copy();							// move into place
bsupumn.ymnc(1,all) = bsupumn.ymnc(2,all) * signhalf;

chk = nc_inq_varid(ncid, "bsupvmnc", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, input.data());	// read
bsupvmn.ymnc.resize(nshalf, mnmax); 						// set size
bsupvmn.ymnc.reindexSelf(index2);							// set indices
bsupvmn.ymnc = input.copy();							// move into place
bsupvmn.ymnc(1,all) = bsupvmn.ymnc(2,all) * signhalf;

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

// read parameter
chk = nc_inq_varid(ncid, "wb", &varid);
chk = nc_get_var_double(ncid, varid, &wb);
chk = nc_inq_varid(ncid, "ctor", &varid);
chk = nc_get_var_double(ncid, varid, &ctor);
chk = nc_inq_varid(ncid, "nextcur", &varid);
chk = nc_get_var_int(ncid, varid, &nextcur);
chk = nc_inq_varid(ncid, "mnmaxpot", &varid);
if(chk == 0) lpot = true;
else lpot = false;
if(lpot) chk = nc_get_var_int(ncid, varid, &mnmaxpot);

// read names
char text[500];
chk = nc_inq_varid(ncid, "input_extension", &varid);
chk = nc_get_var_text(ncid, varid, &text[0]);
input_extension = text; input_extension = input_extension.strip();
text[0] = 0;
chk = nc_inq_varid(ncid, "mgrid_file", &varid);
chk = nc_get_var_text(ncid, varid, &text[0]);
mgrid_file = text; mgrid_file = mgrid_file.left(mgrid_file.indexOf(".nc") + 2);

// read other 1D-arrays
TinyVector <int,1> index1(1);
extcur.resize(nextcur); extcur.reindexSelf(index1);
chk = nc_inq_varid(ncid, "extcur", &varid);
chk = nc_get_var_double(ncid, varid, extcur.data());

if(lpot)
{
	xnpot.resize(mnmaxpot);
	chk = nc_inq_varid(ncid, "xnpot", &varid);
	chk = nc_get_var_int(ncid, varid, xnpot.data());
	xmpot.resize(mnmaxpot);
	chk = nc_inq_varid(ncid, "xmpot", &varid);
	chk = nc_get_var_int(ncid, varid, xmpot.data());
	potsin.resize(mnmaxpot);
	chk = nc_inq_varid(ncid, "potsin", &varid);
	chk = nc_get_var_double(ncid, varid, potsin.data());
	potcos.resize(mnmaxpot);
	chk = nc_inq_varid(ncid, "potcos", &varid);
	chk = nc_get_var_double(ncid, varid, potcos.data());
}

// close file
chk = nc_close(ncid);

//---- prepare stuff ----------------------

// set S arrays
TinyVector <int,1> index(1);	// Array range
S.resize(ns); 			S.reindexSelf(index);
Shalf.resize(nshalf);	Shalf.reindexSelf(index);
double ds = 1.0 / (ns - 1);
for(i=1;i<=ns;i++) S(i) = (i-1)*ds;
for(i=1;i<=nshalf;i++) Shalf(i) = S(i) - ds/2;	// Shalf(1) < 0; Shalf(1) = -Shalf(2)

prepare_splines();
}

//------------------------ get_axis ---------------------------------------------------------------------------------------
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

//------------------------ get_B2D ----------------------------------------------------------------------------------------
void VMEC::get_B2D(double s, double u, double v, double& BR, double& Bphi, double& BZ)
{
double dummy,R,Z,dRdu,dRdv,dZdu,dZdv,bsupu,bsupv;
//if(s > Shalf(nshalf)) // field between s = Shalf(nshalf)  and s = 1
//{
//	BR = 0;
//	BZ = 0;
//	Bphi = 0;
//}
//else	// bsupu & bsupv well defined inside Shalf(nshalf)
//{
	R = rmn.ev(s, u, v, dummy, dRdu, dRdv);
	Z = zmn.ev(s, u, v, dummy, dZdu, dZdv);

	bsupu = bsupumn.ev(s,u,v);
	bsupv = bsupvmn.ev(s,u,v);

	BR = dRdu*bsupu + dRdv*bsupv;
	BZ = dZdu*bsupu + dZdv*bsupv;
	Bphi = R * bsupv;
//}
}

//------------------------ pot --------------------------------------------------------------------------------------------
// evaluates the scalar potential on the s = 1 surface
double VMEC::pot(double u, double v)
{
int i;
double mu_nv;
double out = 0;
if(not lpot) return out;

for(i=0;i<mnmaxpot;i++)
{
	mu_nv = xmpot(i)*u - xnpot(i)*v;
	out += potsin(i) * sin(mu_nv) + potcos(i) * cos(mu_nv);
}
return out;
}

//----------------------------------
// ...with 1st and mixed derivatives
double VMEC::pot(double u, double v, double& dpdu, double& dpdv, double& dpdudv)
{
int i,m,n;
double mu_nv, sinuv, cosuv;
double out = 0;
dpdu = 0; dpdv = 0; dpdudv = 0;
if(not lpot) return out;

for(i=0;i<mnmaxpot;i++)
{
	m = xmpot(i);
	n = xnpot(i);
	mu_nv = m*u - n*v;
	sinuv = sin(mu_nv);
	cosuv = cos(mu_nv);
	out += potsin(i) * sinuv + potcos(i) * cosuv;
	dpdu += m*(potsin(i) * cosuv - potcos(i) * sinuv);
	dpdv += n*(-potsin(i) * cosuv + potcos(i) * sinuv);
	dpdu += m*n*(potsin(i) * sinuv + potcos(i) * cosuv);
}
return out;
}

//------------------------ get_r ------------------------------------------------------------------------------------------
double VMEC::get_r(double R, double Z, double Raxis, double Zaxis)
{
return sqrt((R-Raxis)*(R-Raxis) + (Z-Zaxis)*(Z-Zaxis));
}

//------------------------ get_theta --------------------------------------------------------------------------------------
double VMEC::get_theta(double R, double Z, double Raxis, double Zaxis)
{
double Rm = R - Raxis;
double Zm = Z - Zaxis;
double theta = atan(Zm/Rm);

if(Rm < 0) theta += pi;
if((Rm >= 0) && (Zm < 0)) theta += pi2;

return theta;
}

//------------------------- get_su ---------------------------------------------------------------------------------------
// find s,u at any location (r,phi,z) inside the SIESTA boundary
// sstart and ustart are optional initial guesses of s and u
void VMEC::get_su(double R, double phi, double Z, double& s, double& u, double sstart, double ustart, int imax)
{
double err, Raxis, Zaxis;

if((ustart == -1) || (sstart == -1)) get_axis(phi, Raxis, Zaxis);

// set initial guesses
if(ustart == -1) u = get_theta(R, Z, Raxis, Zaxis);
else u = ustart;
if(sstart == -1) s = bisec(R, phi, Z, u, Raxis, Zaxis, 0, 1, imax);
else s = sstart;

// run Newton
u = modulo2pi(u);	// make sure u in [0, 2pi]
err = newton2D(R,phi,Z,s,u);
if(err > 0) cout << "VMEC Newton2D: no convergence; remaining error: " << err << endl;
}

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
	bsupumn.set("bsupumn", 0, nshalf, mnmax, xn, xm, Shalf);
	bsupvmn.set("bsupvmn", 0, nshalf, mnmax, xn, xm, Shalf);
}
else
{
	rmn.set("rmnc", 1, ns, mnmax, xn, xm, S);
	zmn.set("zmns", -1, ns, mnmax, xn, xm, S);
	gmn.set("gmnc", 1, nshalf, mnmax, xn, xm, Shalf);
	bsupumn.set("bsupumnc", 1, nshalf, mnmax, xn, xm, Shalf);
	bsupvmn.set("bsupvmnc", 1, nshalf, mnmax, xn, xm, Shalf);
}

rmn.Vspline();
zmn.Vspline();
gmn.Vspline();
bsupumn.Vspline();
bsupvmn.Vspline();
}

//------------------------ newton2D ----------------------------------------------------------------------------------------
// 2D Newton procedure
int VMEC::newton2D(double r0, double phi0, double z0, double& s, double& u)
{
const int imax = 20;
const double delta = 1e-12;

int i;
double r,z,drds,drdu,drdv,dzds,dzdu,dzdv;
double det,delta_s,delta_u,fr,fz,err;

for(i=0;i<imax;i++)
{
	r = rmn.ev(s, u, phi0, drds, drdu, drdv);
	z = zmn.ev(s, u, phi0, dzds, dzdu, dzdv);

	fr = r - r0;
	fz = z - z0;

	det = drds * dzdu - drdu * dzds;
	delta_s = (dzdu * fr - drdu * fz) / det;
	delta_u = (drds * fz - dzds * fr) / det;

	err = sqrt(delta_s*delta_s + delta_u*delta_u);
	if(err < delta) return 0;	// convergence

	s -= delta_s;
	u -= delta_u;
	u = modulo2pi(u);	// make sure u in [0, 2pi]
}
return err;	// no convergence
}

//------------------------ bisec ----------------------------------------------------------------------------------------
// bisection with only imax steps to find a crude estimate of s; use as preconditioner for newton2D
// set u fixed as the geometric angle
// assume f(a) < 0 & f(b) > 0 with f = (r(s,u,v) - Raxis)^2 + (z(s,u,v) - Zaxis)^2 - (r0 - Raxis)^2 - (z0 - Zaxis)^2
// and v = phi0
double VMEC::bisec(double r0, double phi0, double z0, double u, double Raxis, double Zaxis, double a, double b, int imax)
{
int i;
double s,rminor,r,z,f;

// Iterations
//const int imax = 10;

// reference point
const double rminor0 = (r0 - Raxis)*(r0 - Raxis) + (z0 - Zaxis)*(z0 - Zaxis);

for(i=1;i<=imax;i++)
{
	s = (a + b)/2.0;
	r = rmn.ev(s, u, phi0) - Raxis;
	z = zmn.ev(s, u, phi0) - Zaxis;
	rminor = r*r + z*z;
	f = rminor - rminor0;
	if(f > 0) b = s;
	else a = s;
}
return s;
}
//------------------------ End of Class VMEC-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Begin Class INSIDE_VMEC -------------------------------------------------------------------------------------
class INSIDE_VMEC
{
private:
	// Member Variables
	int N;
	double phi;
	double Raxis;
	double Zaxis;
	VMEC wout;

	Array<double,1> tha;
	Array<double,1> ra;
	Array<double,1> d2r;

public:
	// Member Variables

	// Constructors
	INSIDE_VMEC();								// Default Constructor
	INSIDE_VMEC(const INSIDE_VMEC& inside);		// Copy Constructor
	INSIDE_VMEC(VMEC woutin);					// set just VMEC Constructor; default: N = 200
	INSIDE_VMEC(VMEC woutin, int Nin);			// set all Constructor

	// Member-Operators
	INSIDE_VMEC& operator =(const INSIDE_VMEC& inside);	// Operator =

	// Member-Functions
	void set_VMEC(VMEC woutin);
	void set_N(int Nin = 200);
	void init(double phiin);
	bool check(double R, double Z);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
INSIDE_VMEC::INSIDE_VMEC()
{
TinyVector <int,1> index(1);	// Array range
N = 0;
ra.reindexSelf(index);
tha.reindexSelf(index);
d2r.reindexSelf(index);
}

// Copy Constructor
INSIDE_VMEC::INSIDE_VMEC(const INSIDE_VMEC& inside)
{
*this = inside;
}

// set just VMEC Constructor
INSIDE_VMEC::INSIDE_VMEC(VMEC woutin)
{
TinyVector <int,1> index(1);	// Array range
set_N();						// default: N = 200
ra.reindexSelf(index);
tha.reindexSelf(index);
d2r.reindexSelf(index);
set_VMEC(woutin);
}

// set all Constructor
INSIDE_VMEC::INSIDE_VMEC(VMEC woutin, int Nin)
{
TinyVector <int,1> index(1);	// Array range
set_N(Nin);
ra.reindexSelf(index);
tha.reindexSelf(index);
d2r.reindexSelf(index);
set_VMEC(woutin);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
INSIDE_VMEC& INSIDE_VMEC::operator =(const INSIDE_VMEC& inside)
{
if (this == &inside) return(*this);	    // if: x=x
N = inside.N;
phi = inside.phi;
Raxis = inside.Raxis;
Zaxis = inside.Zaxis;
wout = inside.wout;
tha.reference(inside.tha);
ra.reference(inside.ra);
d2r.reference(inside.d2r);
return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- set_VMEC ------------------------------------------------------------------------------------------------------
void INSIDE_VMEC::set_VMEC(VMEC woutin)
{
wout = woutin;
}

//--------- set_N ---------------------------------------------------------------------------------------------------------
void INSIDE_VMEC::set_N(int Nin)
{
N = Nin;
tha.resize(2*N);
ra.resize(2*N);
d2r.resize(2*N);
}

//--------- init ----------------------------------------------------------------------------------------------------------
void INSIDE_VMEC::init(double phiin)
{
int i;
double u, Rs, Zs;
phi = phiin;	// v = phi

// magnetic axis
wout.get_axis(phi, Raxis, Zaxis);

// s = 1 surface; ra(N+1) = ra(1); tha(N+1) = tha(1) + pi2; tha & ra have length 2N
double thaold = 0;
for(i=1;i<=N;i++)
{
	u = (i-1)*pi2/N;
	Rs = wout.rmn.ev(1.0, u, phi);
	Zs = wout.zmn.ev(1.0, u, phi);
	ra(i) = wout.get_r(Rs, Zs, Raxis, Zaxis);
	ra(N+i) = ra(i);
	tha(i) = wout.get_theta(Rs, Zs, Raxis, Zaxis);

	// make theta monotonic
	if((tha(i) - thaold) > 5) tha(i) -= pi2;
	if((tha(i) - thaold) < 0) tha(i) += pi2;
	thaold = tha(i);
	tha(N+i) = tha(i);
}
// shift one half of tha into the proper range so that always tha(1) < 0 and tha(2*N) > pi2
if(tha(1) <= 0) tha(Range(N+1,2*N)) += pi2;
else tha(Range(1,N)) -= pi2;

// prepare spline
double dy1 = (ra(N+2) - ra(N))/(tha(N+2) - tha(N));	// dy(1) = dy(N+1)
double dyn = (ra(N+1) - ra(N-1))/(tha(N+1) - tha(N-1)); // dy(2N) = dy(N)
spline(tha, ra, 2*N, dy1, dyn, d2r);
}

//--------- check ---------------------------------------------------------------------------------------------------------
bool INSIDE_VMEC::check(double R, double Z)
{
double rs, dummy;
double theta = wout.get_theta(R, Z, Raxis, Zaxis);
double r = wout.get_r(R, Z, Raxis, Zaxis);

splint(tha, ra, d2r, 2*N, theta, rs, dummy);

return r <= rs;
}

//------------------------ End of Class INSIDE_VMEC------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  VMEC_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

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
#include <andi.hxx>
#include <efit_class.hxx>
using namespace blitz;

// Prototypes
//-----------

// Golbal Parameters
//------------------
extern ofstream ofs2;

//--------- Begin Class PARTICLE ----------------------------------------------------------------------------------------------
class SIESTA
{
private:
// Parameter
	static const double smin = 0;
	static const double smax = 1;

// Member Variables
	double ds, du, dv;
	Array<double,1> Raxis, Zaxis;
	Array<double,1> S,U;
	Array<double,3> dBRds, dBPHIds, dBZds;
	Array<double,3> dBRdu, dBPHIdu, dBZdu;
	Array<double,3> d2BR, d2BPHI, d2BZ;
	Array<double,3> dRds, dZds;
	Array<double,3> dRdu, dZdu;
	Array<double,3> d2R, d2Z;
	Array<double,5> CaBR, CaBPHI, CaBZ, CaR, CaZ;

// Member-Functions
	void prep_interpolation(void);	// initiate the interpolations
	void interpolate_B(double s, double u, int k, double& br, double& bphi, double& bz);	// interpolate all 3 B-field components on the (s,u,v) grid
	void interpolate_RZ(double s, double u, int k, double& r, double& z, double& drds, double& drdu, double& dzds, double& dzdu);	// interpolate R,Z on the (s,u,v) grid
	void find_su(double r, int k, double z, double& s, double& u, double sstart = -1, double ustart = -1);	// phi = v !,  k = index of v-plane -> v = (k-1)*dv; k = 1,...,Nv
	int newton2D(double r, int k, double z, double& s, double& u);	// 2D Newton procedure to find s,u from a given R,Z
	double bisec(double r0, int k, double z0, double u, double a = 0, double b = 1);	// preconditioner to give a rough estimate of s to use as initial condition in newton2D

	//--- old versions, it works, but incredibly slow ---
	void prep_interpolation_2Dsplines(void);	// initiate the interpolations
	void interpolate_B_2Dsplines(double s, double u, int k, double& br, double& bphi, double& bz);	// interpolate all 3 B-field components on the (s,u,v) grid
	void interpolate_RZ_2Dsplines(double s, double u, int k, double& r, double& z, double& drds, double& drdu, double& dzds, double& dzdu);	// interpolate R,Z on the (s,u,v) grid
	//---------------------------------------------------

public: 
// Member Variables
	int Ns, Nu, Nv;

	Array<double,3> R;
	Array<double,3> PHI;
	Array<double,3> Z;
	Array<double,3> BR;
	Array<double,3> BPHI;
	Array<double,3> BZ;

// Constructors
	SIESTA();								// Default Constructor

// Member-Operators
	SIESTA& operator =(const SIESTA& SIES);	// Operator =

// Member-Functions
	void read(LA_STRING filename);			// read in B-field tracer file from SIESTA
	void get_B(double r, double phi, double z, double& br, double& bphi, double& bz);	// evaluate B at any location (r,phi,z) inside the SIESTA boundary
	void get_RZ(double s, double u, double v, double& r, double& z);	// evaluate R,Z at any location (s,u,v)
	void get_su(double r, double phi, double z, double& s, double& u);	// find s,u at any location (r,phi,z) inside the SIESTA boundary
}; //end of class

//------------------------ End of Class -----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
SIESTA::SIESTA()
{
TinyVector <int,1> index(1);
TinyVector <int,3> index3(1,1,1);
TinyVector <int,5> index5(1,1,1,1,1);

Ns = 100;
Nu = 150;
Nv = 80;

ds = (smax - smin) / (Ns-1);	// s goes: smin -> smax
du = pi2/(Nu-1);				// 0 -> 2*pi
dv = pi2/(Nv-1);				// 0 -> 2*pi

S.resize(Ns);				S.reindexSelf(index);
U.resize(Nu);				U.reindexSelf(index);

R.resize(Ns, Nu, Nv);		R.reindexSelf(index3);
PHI.resize(Ns, Nu, Nv);		PHI.reindexSelf(index3);
Z.resize(Ns, Nu, Nv);		Z.reindexSelf(index3);
BR.resize(Ns, Nu, Nv);		BR.reindexSelf(index3);
BPHI.resize(Ns, Nu, Nv);	BPHI.reindexSelf(index3);
BZ.resize(Ns, Nu, Nv);		BZ.reindexSelf(index3);

dBRds.resize(Ns,Nu,Nv);		dBRds.reindexSelf(index3);
dBPHIds.resize(Ns,Nu,Nv);	dBPHIds.reindexSelf(index3);
dBZds.resize(Ns,Nu,Nv);		dBZds.reindexSelf(index3);

dRds.resize(Ns,Nu,Nv);		dRds.reindexSelf(index3);
dZds.resize(Ns,Nu,Nv);		dZds.reindexSelf(index3);

dBRdu.resize(Ns,Nu,Nv);		dBRdu.reindexSelf(index3);
dBPHIdu.resize(Ns,Nu,Nv);	dBPHIdu.reindexSelf(index3);
dBZdu.resize(Ns,Nu,Nv);		dBZdu.reindexSelf(index3);

dRdu.resize(Ns,Nu,Nv);		dRdu.reindexSelf(index3);
dZdu.resize(Ns,Nu,Nv);		dZdu.reindexSelf(index3);

d2BR.resize(Ns,Nu,Nv);		d2BR.reindexSelf(index3);
d2BPHI.resize(Ns,Nu,Nv);	d2BPHI.reindexSelf(index3);
d2BZ.resize(Ns,Nu,Nv);		d2BZ.reindexSelf(index3);

d2R.resize(Ns,Nu,Nv);		d2R.reindexSelf(index3);
d2Z.resize(Ns,Nu,Nv);		d2Z.reindexSelf(index3);

CaBR.resize(Ns-1,Nu-1,Nv-1,4,4);	CaBR.reindexSelf(index5);
CaBPHI.resize(Ns-1,Nu-1,Nv-1,4,4);	CaBPHI.reindexSelf(index5);
CaBZ.resize(Ns-1,Nu-1,Nv-1,4,4);	CaBZ.reindexSelf(index5);
CaR.resize(Ns-1,Nu-1,Nv-1,4,4);		CaR.reindexSelf(index5);
CaZ.resize(Ns-1,Nu-1,Nv-1,4,4);		CaZ.reindexSelf(index5);

Raxis.resize(Nv);			Raxis.reindexSelf(index);
Zaxis.resize(Nv);			Zaxis.reindexSelf(index);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
SIESTA& SIESTA::operator =(const SIESTA& SIES)
{
if (this == &SIES) return(*this);	    // if: x=x

Ns = SIES.Ns;
Nu = SIES.Nu;
Nv = SIES.Nv;

ds = SIES.ds;
du = SIES.du;
dv = SIES.dv;

S.reference(SIES.S.copy());
U.reference(SIES.U.copy());

R.reference(SIES.R.copy());
PHI.reference(SIES.PHI.copy());
Z.reference(SIES.Z.copy());
BR.reference(SIES.BR.copy());
BPHI.reference(SIES.BPHI.copy());
BZ.reference(SIES.BZ.copy());

dBRds.reference(SIES.dBRds.copy());
dBPHIds.reference(SIES.dBPHIds.copy());
dBZds.reference(SIES.dBZds.copy());

dRds.reference(SIES.dRds.copy());
dZds.reference(SIES.dZds.copy());

dBRdu.reference(SIES.dBRdu.copy());
dBPHIdu.reference(SIES.dBPHIdu.copy());
dBZdu.reference(SIES.dBZdu.copy());

dRdu.reference(SIES.dRdu.copy());
dZdu.reference(SIES.dZdu.copy());

d2BR.reference(SIES.d2BR.copy());
d2BPHI.reference(SIES.d2BPHI.copy());
d2BZ.reference(SIES.d2BZ.copy());

d2R.reference(SIES.d2R.copy());
d2Z.reference(SIES.d2Z.copy());

CaBR.reference(SIES.CaBR.copy());
CaBPHI.reference(SIES.CaBPHI.copy());
CaBZ.reference(SIES.CaBZ.copy());
CaR.reference(SIES.CaR.copy());
CaZ.reference(SIES.CaZ.copy());

Raxis.reference(SIES.Raxis.copy());
Zaxis.reference(SIES.Zaxis.copy());

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------------------- read ---------------------------------------------------------------------------------------
// read SIESTA bfield_tracing file and set member variables
void SIESTA::read(LA_STRING filename)
{
// Variables
int i,j,k;
Range all = Range::all();

// Input
ifstream in;
in.open(filename);
if(in.fail()==1) {cout << "Unable to open file " << filename << endl; EXIT;}

// Read the first line, 3 integers, to get dimensions
in >> Nv;
in >> Ns;
in >> Nu;

// set grid; increase Nu & Nv by 1 to make periodic
Nu += 1;
Nv += 1;

ds = (smax - smin) / (Ns-1);	// s goes: smin -> smax
du = pi2/(Nu-1);				// 0 -> 2*pi
dv = pi2/(Nv-1);				// 0 -> 2*pi

S.resize(Ns);
for(i=1;i<=Ns;i++) S(i) = smin + (i-1)*ds;
U.resize(Nu);
for(i=1;i<=Nu;i++) U(i) = (i-1)*du;

// resize data: (0 -> Ns-1, 0 -> Nu-1, 0 -> Nv-1)
R.resize(Ns, Nu, Nv);
PHI.resize(Ns, Nu, Nv);
Z.resize(Ns, Nu, Nv);
BR.resize(Ns, Nu, Nv);
BPHI.resize(Ns, Nu, Nv);
BZ.resize(Ns, Nu, Nv);

// Read data
for(i=1;i<=Nu-1;i++)	// go one less, since Nu is one larger
{
	for(j=1;j<=Nv-1;j++)	// go one less, since Nv is one larger
	{
		for(k=1;k<=Ns;k++)
		{
			in >> R(k,i,j);
			in >> Z(k,i,j);
			in >> PHI(k,i,j);
			in >> BR(k,i,j);
			in >> BZ(k,i,j);
			in >> BPHI(k,i,j);
		}
	}
}
in.close();

// make periodic
R(all,Nu,all) = R(all,1,all);
R(all,all,Nv) = R(all,all,1);
PHI(all,Nu,all) = PHI(all,1,all);
PHI(all,all,Nv) = pi2;
Z(all,Nu,all) = Z(all,1,all);
Z(all,all,Nv) = Z(all,all,1);
BR(all,Nu,all) = BR(all,1,all);
BR(all,all,Nv) = BR(all,all,1);
BPHI(all,Nu,all) = BPHI(all,1,all);
BPHI(all,all,Nv) = BPHI(all,all,1);
BZ(all,Nu,all) = BZ(all,1,all);
BZ(all,all,Nv) = BZ(all,all,1);

// get a guess of the Axis location
Raxis.resize(Nv);
Zaxis.resize(Nv);
for(k=1;k<=Nv;k++)
{
	Raxis(k) = sum(R(1,all,k)) / Nu;
	Zaxis(k) = sum(Z(1,all,k)) / Nu;
}

// prepare the bicubic splines
prep_interpolation();
}

//-------------------------------- get_B ----------------------------------------------------------------------------------
// evaluate br, bphi, bz at any arbitrary location (r, phi, z) inside the SIESTA boundary; r is the major radius here!
void SIESTA::get_B(double r, double phi, double z, double& br, double& bphi, double& bz)
{
int k;
double t;
double sl, su, ul, uu;
double bru, bphiu, bzu;

// locate neighboring phi planes: lower (or matching) phi = (k-1)*dv; upper phi = k*dv
phi = modulo2pi(phi);	// make sure phi in [0, 2pi]
k = int(phi/dv) + 1;
t = phi/dv - k + 1;

// transform to flux coordinates & bicubic spline interpolation of B-field in lower phi plane
find_su(r, k, z, sl, ul);
interpolate_B(sl, ul, k, br, bphi, bz);

if(t > 0)	// not exactly in the v-plane (k-1)*dv
{
	// transform to flux coordinates & bicubic spline interpolation of B-field in upper phi plane
	find_su(r, k+1, z, su, uu, sl, ul);	// use sl and ul as initial guess for su and uu
	interpolate_B(su, uu, k+1, bru, bphiu, bzu);

	// linear interpolation of B-field in phi
	br += t*(bru - br);
	bphi += t*(bphiu - bphi);
	bz += t*(bzu - bz);
}
}

//------------------------- get_RZ ---------------------------------------------------------------------------------
// gets r(s,u), z(s,u) in the v-plane (k-1)*dv through bicubic spline interpolation
// calls the private member function interpolate_RZ and omits the derivatives
void SIESTA::get_RZ(double s, double u, double v, double& r, double& z)
{
int k;
double dummy, t;
double ru,zu;

v = modulo2pi(v);	// make sure v in [0, 2pi]
k = int(v/dv) + 1;
t = v/dv - k + 1;

interpolate_RZ(s, u, k, r, z, dummy, dummy, dummy, dummy);

if(t > 0)	// not exactly in the v-plane (k-1)*dv
{
	interpolate_RZ(s, u, k+1, ru, zu, dummy, dummy, dummy, dummy);

	// linear interpolation of RZ in v
	r += t*(ru - r);
	z += t*(zu - z);
}
}

//------------------------- get_su ---------------------------------------------------------------------------------
// find s,u at any location (r,phi,z) inside the SIESTA boundary
void SIESTA::get_su(double r, double phi, double z, double& s, double& u)
{
int k;
double t,su,uu;

// locate neighboring phi planes: lower (or matching) phi = (k-1)*dv; upper phi = k*dv
phi = modulo2pi(phi);	// make sure phi in [0, 2pi]
k = int(phi/dv) + 1;
t = phi/dv - k + 1;

find_su(r, k, z, s, u);

if(t > 0)	// not exactly in the v-plane (k-1)*dv
{
	find_su(r, k+1, z, su, uu);

	// linear interpolation of s,u in v
	s += t*(su - s);
	u += t*(uu - u);
}
}

//----------------------- End of Public Member Functions ------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- Private Member Functions ------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------------------- prep_interpolation --------------------------------------------------------------------------------
// get the 2. derivative in s-direction, which is input to splint_2D, for all three fields, all u's in all v-planes
void SIESTA::prep_interpolation(void)
{
int i,j,k;
Array<double,2> slice, sliceds, slicedu, sliced2;
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
Range all = Range::all();

dBRds.resize(Ns,Nu,Nv);
dBPHIds.resize(Ns,Nu,Nv);
dBZds.resize(Ns,Nu,Nv);

dRds.resize(Ns,Nu,Nv);
dZds.resize(Ns,Nu,Nv);

dBRdu.resize(Ns,Nu,Nv);
dBPHIdu.resize(Ns,Nu,Nv);
dBZdu.resize(Ns,Nu,Nv);

dRdu.resize(Ns,Nu,Nv);
dZdu.resize(Ns,Nu,Nv);

d2BR.resize(Ns,Nu,Nv);
d2BPHI.resize(Ns,Nu,Nv);
d2BZ.resize(Ns,Nu,Nv);

d2R.resize(Ns,Nu,Nv);
d2Z.resize(Ns,Nu,Nv);

CaBR.resize(Ns-1,Nu-1,Nv-1,4,4);
CaBPHI.resize(Ns-1,Nu-1,Nv-1,4,4);
CaBZ.resize(Ns-1,Nu-1,Nv-1,4,4);
CaR.resize(Ns-1,Nu-1,Nv-1,4,4);
CaZ.resize(Ns-1,Nu-1,Nv-1,4,4);

// Get the derivatives
for(i=1;i<=Nv;i++)
{
	slice.reference(BR(all,all,i));
	sliceds.reference(dBRds(all,all,i));
	slicedu.reference(dBRdu(all,all,i));
	sliced2.reference(d2BR(all,all,i));
	bcuderiv(slice,ds,du,sliceds,slicedu,sliced2);

	slice.reference(BPHI(all,all,i));
	sliceds.reference(dBPHIds(all,all,i));
	slicedu.reference(dBPHIdu(all,all,i));
	sliced2.reference(d2BPHI(all,all,i));
	bcuderiv(slice,ds,du,sliceds,slicedu,sliced2);

	slice.reference(BZ(all,all,i));
	sliceds.reference(dBZds(all,all,i));
	slicedu.reference(dBZdu(all,all,i));
	sliced2.reference(d2BZ(all,all,i));
	bcuderiv(slice,ds,du,sliceds,slicedu,sliced2);

	slice.reference(R(all,all,i));
	sliceds.reference(dRds(all,all,i));
	slicedu.reference(dRdu(all,all,i));
	sliced2.reference(d2R(all,all,i));
	bcuderiv(slice,ds,du,sliceds,slicedu,sliced2);

	slice.reference(Z(all,all,i));
	sliceds.reference(dZds(all,all,i));
	slicedu.reference(dZdu(all,all,i));
	sliced2.reference(d2Z(all,all,i));
	bcuderiv(slice,ds,du,sliceds,slicedu,sliced2);
}

// Make derivatives periodic in u through symmetry
// df1 = (f2 - f1)/du;	dfN = (fN - fN-1)/du;
// periodicity: df1 = dfN = (df1 + dfN)/2
// df1 = dfN = (f2 - f1 + fN - fN-1)/2du = (f2 - fN-1)/2du
//dBRds(all,1,all) = (dBRds(all,1,all) + dBRds(all,Nu,all))/2.0;	dBRds(all,Nu,all) = dBRds(all,1,all);
dBRdu(all,1,all) = (dBRdu(all,1,all) + dBRdu(all,Nu,all))/2.0;	dBRdu(all,Nu,all) = dBRdu(all,1,all);
d2BR(all,1,all) = (d2BR(all,1,all) + d2BR(all,Nu,all))/2.0;		d2BR(all,Nu,all) = d2BR(all,1,all);

//dBPHIds(all,1,all) = (dBPHIds(all,1,all) + dBPHIds(all,Nu,all))/2.0;	dBPHIds(all,Nu,all) = dBPHIds(all,1,all);
dBPHIdu(all,1,all) = (dBPHIdu(all,1,all) + dBPHIdu(all,Nu,all))/2.0;	dBPHIdu(all,Nu,all) = dBPHIdu(all,1,all);
d2BPHI(all,1,all) = (d2BPHI(all,1,all) + d2BPHI(all,Nu,all))/2.0;		d2BPHI(all,Nu,all) = d2BPHI(all,1,all);

//dBZds(all,1,all) = (dBZds(all,1,all) + dBZds(all,Nu,all))/2.0;	dBZds(all,Nu,all) = dBZds(all,1,all);
dBZdu(all,1,all) = (dBZdu(all,1,all) + dBZdu(all,Nu,all))/2.0;	dBZdu(all,Nu,all) = dBZdu(all,1,all);
d2BZ(all,1,all) = (d2BZ(all,1,all) + d2BZ(all,Nu,all))/2.0;		d2BZ(all,Nu,all) = d2BZ(all,1,all);

//dRds(all,1,all) = (dRds(all,1,all) + dRds(all,Nu,all))/2.0;		dRds(all,Nu,all) = dRds(all,1,all);
dRdu(all,1,all) = (dRdu(all,1,all) + dRdu(all,Nu,all))/2.0;		dRdu(all,Nu,all) = dRdu(all,1,all);
d2R(all,1,all) = (d2R(all,1,all) + d2R(all,Nu,all))/2.0;		d2R(all,Nu,all) = d2R(all,1,all);

//dZds(all,1,all) = (dZds(all,1,all) + dZds(all,Nu,all))/2.0;		dZds(all,Nu,all) = dZds(all,1,all);
dZdu(all,1,all) = (dZdu(all,1,all) + dZdu(all,Nu,all))/2.0;		dZdu(all,Nu,all) = dZdu(all,1,all);
d2Z(all,1,all) = (d2Z(all,1,all) + d2Z(all,Nu,all))/2.0;		d2Z(all,Nu,all) = d2Z(all,1,all);

// Get the c's for bcuint, as done by bcucof
for(i=1;i<Ns;i++)
{
	for(j=1;j<Nu;j++)
	{
		for(k=1;k<Nv;k++)
		{
			y_sq(1) = BR(i,j,k); y_sq(2) = BR(i+1,j,k); y_sq(3) = BR(i+1,j+1,k); y_sq(4) = BR(i,j+1,k);
			y1_sq(1) = dBRds(i,j,k); y1_sq(2) = dBRds(i+1,j,k); y1_sq(3) = dBRds(i+1,j+1,k); y1_sq(4) = dBRds(i,j+1,k);
			y2_sq(1) = dBRdu(i,j,k); y2_sq(2) = dBRdu(i+1,j,k); y2_sq(3) = dBRdu(i+1,j+1,k); y2_sq(4) = dBRdu(i,j+1,k);
			y12_sq(1) = d2BR(i,j,k); y12_sq(2) = d2BR(i+1,j,k); y12_sq(3) = d2BR(i+1,j+1,k); y12_sq(4) = d2BR(i,j+1,k);
			slice.reference(CaBR(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,ds,du,slice);

			y_sq(1) = BPHI(i,j,k); y_sq(2) = BPHI(i+1,j,k); y_sq(3) = BPHI(i+1,j+1,k); y_sq(4) = BPHI(i,j+1,k);
			y1_sq(1) = dBPHIds(i,j,k); y1_sq(2) = dBPHIds(i+1,j,k); y1_sq(3) = dBPHIds(i+1,j+1,k); y1_sq(4) = dBPHIds(i,j+1,k);
			y2_sq(1) = dBPHIdu(i,j,k); y2_sq(2) = dBPHIdu(i+1,j,k); y2_sq(3) = dBPHIdu(i+1,j+1,k); y2_sq(4) = dBPHIdu(i,j+1,k);
			y12_sq(1) = d2BPHI(i,j,k); y12_sq(2) = d2BPHI(i+1,j,k); y12_sq(3) = d2BPHI(i+1,j+1,k); y12_sq(4) = d2BPHI(i,j+1,k);
			slice.reference(CaBPHI(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,ds,du,slice);

			y_sq(1) = BZ(i,j,k); y_sq(2) = BZ(i+1,j,k); y_sq(3) = BZ(i+1,j+1,k); y_sq(4) = BZ(i,j+1,k);
			y1_sq(1) = dBZds(i,j,k); y1_sq(2) = dBZds(i+1,j,k); y1_sq(3) = dBZds(i+1,j+1,k); y1_sq(4) = dBZds(i,j+1,k);
			y2_sq(1) = dBZdu(i,j,k); y2_sq(2) = dBZdu(i+1,j,k); y2_sq(3) = dBZdu(i+1,j+1,k); y2_sq(4) = dBZdu(i,j+1,k);
			y12_sq(1) = d2BZ(i,j,k); y12_sq(2) = d2BZ(i+1,j,k); y12_sq(3) = d2BZ(i+1,j+1,k); y12_sq(4) = d2BZ(i,j+1,k);
			slice.reference(CaBZ(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,ds,du,slice);

			y_sq(1) = R(i,j,k); y_sq(2) = R(i+1,j,k); y_sq(3) = R(i+1,j+1,k); y_sq(4) = R(i,j+1,k);
			y1_sq(1) = dRds(i,j,k); y1_sq(2) = dRds(i+1,j,k); y1_sq(3) = dRds(i+1,j+1,k); y1_sq(4) = dRds(i,j+1,k);
			y2_sq(1) = dRdu(i,j,k); y2_sq(2) = dRdu(i+1,j,k); y2_sq(3) = dRdu(i+1,j+1,k); y2_sq(4) = dRdu(i,j+1,k);
			y12_sq(1) = d2R(i,j,k); y12_sq(2) = d2R(i+1,j,k); y12_sq(3) = d2R(i+1,j+1,k); y12_sq(4) = d2R(i,j+1,k);
			slice.reference(CaR(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,ds,du,slice);

			y_sq(1) = Z(i,j,k); y_sq(2) = Z(i+1,j,k); y_sq(3) = Z(i+1,j+1,k); y_sq(4) = Z(i,j+1,k);
			y1_sq(1) = dZds(i,j,k); y1_sq(2) = dZds(i+1,j,k); y1_sq(3) = dZds(i+1,j+1,k); y1_sq(4) = dZds(i,j+1,k);
			y2_sq(1) = dZdu(i,j,k); y2_sq(2) = dZdu(i+1,j,k); y2_sq(3) = dZdu(i+1,j+1,k); y2_sq(4) = dZdu(i,j+1,k);
			y12_sq(1) = d2Z(i,j,k); y12_sq(2) = d2Z(i+1,j,k); y12_sq(3) = d2Z(i+1,j+1,k); y12_sq(4) = d2Z(i,j+1,k);
			slice.reference(CaZ(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,ds,du,slice);
		}
	}
}
}

void SIESTA::prep_interpolation_2Dsplines(void)
{
int i,j;
double d1,dn;
Array<double,1> slice_1,slice_2;
Range all = Range::all();

d2BR.resize(Ns,Nu,Nv);
d2BPHI.resize(Ns,Nu,Nv);
d2BZ.resize(Ns,Nu,Nv);

d2R.resize(Ns,Nu,Nv);
d2Z.resize(Ns,Nu,Nv);

for(j=1;j<=Nv;j++)
{
	for(i=1;i<=Nu;i++)
	{
		d1 = (BR(2,i,j) - BR(1,i,j)) / ds;
		dn = (BR(Ns,i,j) - BR(Ns-1,i,j)) / ds;
		slice_1.reference(BR(all,i,j));
		slice_2.reference(d2BR(all,i,j));
		spline(S,slice_1,Ns,d1,dn,slice_2);

		d1 = (BPHI(2,i,j) - BPHI(1,i,j)) / ds;
		dn = (BPHI(Ns,i,j) - BPHI(Ns-1,i,j)) / ds;
		slice_1.reference(BPHI(all,i,j));
		slice_2.reference(d2BPHI(all,i,j));
		spline(S,slice_1,Ns,d1,dn,slice_2);

		d1 = (BZ(2,i,j) - BZ(1,i,j)) / ds;
		dn = (BZ(Ns,i,j) - BZ(Ns-1,i,j)) / ds;
		slice_1.reference(BZ(all,i,j));
		slice_2.reference(d2BZ(all,i,j));
		spline(S,slice_1,Ns,d1,dn,slice_2);

		d1 = (R(2,i,j) - R(1,i,j)) / ds;
		dn = (R(Ns,i,j) - R(Ns-1,i,j)) / ds;
		slice_1.reference(R(all,i,j));
		slice_2.reference(d2R(all,i,j));
		spline(S,slice_1,Ns,d1,dn,slice_2);

		d1 = (Z(2,i,j) - Z(1,i,j)) / ds;
		dn = (Z(Ns,i,j) - Z(Ns-1,i,j)) / ds;
		slice_1.reference(Z(all,i,j));
		slice_2.reference(d2Z(all,i,j));
		spline(S,slice_1,Ns,d1,dn,slice_2);
	}
}
}

//------------------------- interpolate_B ---------------------------------------------------------------------------------
// gets br(s,u), bphi(s,u), bz(s,u) in the v-plane (k-1)*dv through bicubic spline interpolation
void SIESTA::interpolate_B(double s, double u, int k, double& br, double& bphi, double& bz)
{
double dummy;
Array<double,4> slice;
Range all = Range::all();

slice.reference(CaBR(all,all,k,all,all));
bcuint(S,U,slice,ds,du,s,u,br,dummy,dummy);

slice.reference(CaBPHI(all,all,k,all,all));
bcuint(S,U,slice,ds,du,s,u,bphi,dummy,dummy);

slice.reference(CaBZ(all,all,k,all,all));
bcuint(S,U,slice,ds,du,s,u,bz,dummy,dummy);
}

void SIESTA::interpolate_B_2Dsplines(double s, double u, int k, double& br, double& bphi, double& bz)
{
double dummy;
Array<double,2> slice_1,slice_2;
Range all = Range::all();

slice_1.reference(BR(all,all,k));
slice_2.reference(d2BR(all,all,k));
splint_2D(S, U, slice_1, slice_2, Ns, Nu, s, u, br, dummy, dummy);

slice_1.reference(BPHI(all,all,k));
slice_2.reference(d2BPHI(all,all,k));
splint_2D(S, U, slice_1, slice_2, Ns, Nu, s, u, bphi, dummy, dummy);

slice_1.reference(BZ(all,all,k));
slice_2.reference(d2BZ(all,all,k));
splint_2D(S, U, slice_1, slice_2, Ns, Nu, s, u, bz, dummy, dummy);
}

//------------------------- interpolate_RZ ---------------------------------------------------------------------------------
// gets r(s,u), z(s,u) in the v-plane (k-1)*dv through bicubic spline interpolation
void SIESTA::interpolate_RZ(double s, double u, int k, double& r, double& z, double& drds, double& drdu, double& dzds, double& dzdu)
{
Array<double,4> slice;
Range all = Range::all();

slice.reference(CaR(all,all,k,all,all));
bcuint(S,U,slice,ds,du,s,u,r,drds,drdu);

slice.reference(CaZ(all,all,k,all,all));
bcuint(S,U,slice,ds,du,s,u,z,dzds,dzdu);
}

void SIESTA::interpolate_RZ_2Dsplines(double s, double u, int k, double& r, double& z, double& drds, double& drdu, double& dzds, double& dzdu)
{
Array<double,2> slice_1,slice_2;
Range all = Range::all();

slice_1.reference(R(all,all,k));
slice_2.reference(d2R(all,all,k));
splint_2D(S, U, slice_1, slice_2, Ns, Nu, s, u, r, drds, drdu);

slice_1.reference(Z(all,all,k));
slice_2.reference(d2Z(all,all,k));
splint_2D(S, U, slice_1, slice_2, Ns, Nu, s, u, z, dzds, dzdu);
}

//------------------------- find_su ----------------------------------------------------------------------------------------
// use a 2D Newton procedure to get s,u from R, Z in the (k-1)*dv plane
void SIESTA::find_su(double r, int k, double z, double& s, double& u, double sstart, double ustart)
{
double err;

// set initial guesses
if(ustart == -1) u = polar_phi(r - Raxis(k), z - Zaxis(k));
else u = ustart;
if(sstart == -1) s = bisec(r, k, z, u);
else s = sstart;

// run Newton
u = modulo2pi(u);	// make sure u in [0, 2pi]
err = newton2D(r,k,z,s,u);
if(err > 0) cout << "SIESTA Newton2D: no convergence; remaining error: " << err << endl;
}

//------------------------ newton2D ----------------------------------------------------------------------------------------
// 2D Newton procedure
int SIESTA::newton2D(double r0, int k, double z0, double& s, double& u)
{
const int imax = 20;
const double delta = 1e-12;

int i;
double r,z,drds,drdu,dzds,dzdu;
double det,delta_s,delta_u,fr,fz,err;

for(i=0;i<imax;i++)
{
	interpolate_RZ(s, u, k, r, z, drds, drdu, dzds, dzdu);

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
// assume f(a) < 0 & f(b) > 0 with f = (r(s,u,v_k) - Raxis)^2 + (z(s,u,v_k) - Zaxis)^2 - (r0 - Raxis)^2 - (z0 - Zaxis)^2
// and v_k = (k-1)*dv
double SIESTA::bisec(double r0, int k, double z0, double u, double a, double b)
{
int i;
double s,dummy,rminor,r,z,f;

// Iterations
const int imax = 4;

// reference point
const double rminor0 = (r0 - Raxis(k))*(r0 - Raxis(k)) + (z0 - Zaxis(k))*(z0 - Zaxis(k));

for(i=1;i<=imax;i++)
{
	s = (a + b)/2.0;
	interpolate_RZ(s, u, k, r, z, dummy, dummy, dummy, dummy);
	r -= Raxis(k);
	z -= Zaxis(k);
	rminor = r*r + z*z;
	f = rminor - rminor0;
	if(f > 0) b = s;
	else a = s;
}
return s;
}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  SIESTA_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// Class calculates B-field outside of VMEC boundary
// A.Wingen						9.10.14


// Define
//--------
#ifndef EXTENDER_CLASS_INCLUDED
#define EXTENDER_CLASS_INCLUDED

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

// Golbal Parameters
//------------------


//--------- Begin Class MGRID ---------------------------------------------------------------------------------------------
class MGRID
{
private:
	// Member Variables
	LA_STRING mgrid_file;		// pathname of MGRID file
	Array<double,1> extcur;		// coil currents 1->Ncoils
	int Ncoils;					// number of coils in mgrid
	double dR, dZ, dp;
	Array<double,3> dBRdR, dBPHIdR, dBZdR;
	Array<double,3> dBRdZ, dBPHIdZ, dBZdZ;
	Array<double,3> d2BR, d2BPHI, d2BZ;
	Array<double,5> CaBR, CaBPHI, CaBZ;

	// Member-Functions

public:
	// Member Variables
	int NR, NZ, Np;
	Array<double,1> R;
	Array<double,1> Z;
	Array<double,3> BR;
	Array<double,3> BPHI;
	Array<double,3> BZ;

	// Constructors
	MGRID();								// Default Constructor

	// Member-Operators
	MGRID& operator =(const MGRID& mgrid);	// Operator =

	// Member-Functions
	void read(LA_STRING file, int N, Array<double,1> cur);
	void prep_interpolation(void);	// initiate the interpolations
	void interpolate_B(double r, double z, int k, double& br, double& bphi, double& bz);	// interpolate all 3 vacuum B-field components on the mgrid grid
}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
MGRID::MGRID()
{
TinyVector <int,1> index(1);
TinyVector <int,3> index3(1,1,1);
TinyVector <int,5> index5(1,1,1,1,1);

mgrid_file = "None";
Ncoils = 36;
extcur.resize(Ncoils); extcur.reindexSelf(index);

NR = 129;
NZ = 129;
Np = 48;

dR = 0.001;
dZ = 0.001;
dp = pi2/(Np-1);				// 0 -> 2*pi

R.resize(NR);				R.reindexSelf(index);
Z.resize(NZ);				Z.reindexSelf(index);

// these get transposed in read()
BR.resize(Np, NZ, NR);		BR.reindexSelf(index3);
BPHI.resize(Np, NZ, NR);	BPHI.reindexSelf(index3);
BZ.resize(Np, NZ, NR);		BZ.reindexSelf(index3);

// these are correct
dBRdR.resize(NR, NZ, Np);	dBRdR.reindexSelf(index3);
dBPHIdR.resize(NR, NZ, Np);	dBPHIdR.reindexSelf(index3);
dBZdR.resize(NR, NZ, Np);	dBZdR.reindexSelf(index3);

dBRdZ.resize(NR, NZ, Np);	dBRdZ.reindexSelf(index3);
dBPHIdZ.resize(NR, NZ, Np);	dBPHIdZ.reindexSelf(index3);
dBZdZ.resize(NR, NZ, Np);	dBZdZ.reindexSelf(index3);

d2BR.resize(NR, NZ, Np);	d2BR.reindexSelf(index3);
d2BPHI.resize(NR, NZ, Np);	d2BPHI.reindexSelf(index3);
d2BZ.resize(NR, NZ, Np);	d2BZ.reindexSelf(index3);

CaBR.resize(NR-1,NZ-1,Np,4,4);	CaBR.reindexSelf(index5);
CaBPHI.resize(NR-1,NZ-1,Np,4,4);CaBPHI.reindexSelf(index5);
CaBZ.resize(NR-1,NZ-1,Np,4,4);	CaBZ.reindexSelf(index5);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
MGRID& MGRID::operator =(const MGRID& mgrid)
{
if (this == &mgrid) return(*this);	    // if: x=x

mgrid_file = mgrid.mgrid_file;
Ncoils = mgrid.Ncoils;
extcur.reference(mgrid.extcur);

NR = mgrid.NR;
NZ = mgrid.NZ;
Np = mgrid.Np;

dR = mgrid.dR;
dZ = mgrid.dZ;
dp = mgrid.dp;

R.reference(mgrid.R);
Z.reference(mgrid.Z);

BR.reference(mgrid.BR);
BPHI.reference(mgrid.BPHI);
BZ.reference(mgrid.BZ);

dBRdR.reference(mgrid.dBRdR);
dBPHIdR.reference(mgrid.dBPHIdR);
dBZdR.reference(mgrid.dBZdR);

dBRdZ.reference(mgrid.dBRdZ);
dBPHIdZ.reference(mgrid.dBPHIdZ);
dBZdZ.reference(mgrid.dBZdZ);

d2BR.reference(mgrid.d2BR);
d2BPHI.reference(mgrid.d2BPHI);
d2BZ.reference(mgrid.d2BZ);

CaBR.reference(mgrid.CaBR);
CaBPHI.reference(mgrid.CaBPHI);
CaBZ.reference(mgrid.CaBZ);

return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- read ---------------------------------------------------------------------------------------------------------------
void MGRID::read(LA_STRING file, int N, Array<double,1> cur)
{
// Variables
int i,chk, ncid, varid;
firstIndex k;
stringstream ss;
string s;
LA_STRING varname;

// assign private members
mgrid_file = file;
Ncoils = N;
extcur.resize(Ncoils);
extcur.reference(cur);

// open file
chk = nc_open(mgrid_file, NC_NOWRITE, &ncid);
if(chk != 0) {cout << "Unable to open file " << mgrid_file << endl; EXIT;}

// read grid data
chk = nc_inq_varid(ncid, "ir", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &NR);	// read variable
chk = nc_inq_varid(ncid, "jz", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &NZ);	// read variable
chk = nc_inq_varid(ncid, "kp", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &Np);	// read variable

double rmin, rmax, zmin, zmax;
chk = nc_inq_varid(ncid, "rmin", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &rmin);		// read
chk = nc_inq_varid(ncid, "rmax", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &rmax);		// read
chk = nc_inq_varid(ncid, "zmin", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &zmin);		// read
chk = nc_inq_varid(ncid, "zmax", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &zmax);		// read

R.resize(NR);
Z.resize(NZ);
dR = (rmax - rmin)/double(NR);
dZ = (zmax - zmin)/double(NZ);
dp = pi2/(Np-1);
R = rmin + k*dR;
Z = zmin + k*dZ;

// read B-field data
Array<double,3> input;
input.resize(Np, NZ, NR);
BR.resize(Np, NZ, NR); 		BR = 0;
BPHI.resize(Np, NZ, NR);	BPHI = 0;
BZ.resize(Np, NZ, NR);		BZ = 0;

for (i=1;i<=Ncoils;i++)
{
	ss << setw(3) << setfill('0') << i; s = ss.str();

	varname = "br_" + LA_STRING(s.c_str());
	chk = nc_inq_varid(ncid, varname, &varid);
	chk = nc_get_var_double(ncid, varid, input.data());
	BR += extcur(i) * input.copy();

	varname = "bp_" + LA_STRING(s.c_str());
	chk = nc_inq_varid(ncid, varname, &varid);
	chk = nc_get_var_double(ncid, varid, input.data());
	BPHI += extcur(i) * input.copy();

	varname = "bz_" + LA_STRING(s.c_str());
	chk = nc_inq_varid(ncid, varname, &varid);
	chk = nc_get_var_double(ncid, varid, input.data());
	BZ += extcur(i) * input.copy();

	// reset stringstream
	ss.str("");
	ss.clear();
}
// make B-field arrays: (NR, NZ, Np)
BR.transposeSelf(thirdDim, secondDim, firstDim);
BPHI.transposeSelf(thirdDim, secondDim, firstDim);
BZ.transposeSelf(thirdDim, secondDim, firstDim);

//prep_interpolation();
}

//--------------------- prep_interpolation --------------------------------------------------------------------------------
// get the C coefficients for the interpolation
void MGRID::prep_interpolation(void)
{
int i,j,k;
Array<double,2> slice, slicedR, slicedZ, sliced2;
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
Range all = Range::all();

dBRdR.resize(NR, NZ, Np);
dBPHIdR.resize(NR, NZ, Np);
dBZdR.resize(NR, NZ, Np);

dBRdZ.resize(NR, NZ, Np);
dBPHIdZ.resize(NR, NZ, Np);
dBZdZ.resize(NR, NZ, Np);

d2BR.resize(NR, NZ, Np);
d2BPHI.resize(NR, NZ, Np);
d2BZ.resize(NR, NZ, Np);

CaBR.resize(NR-1,NZ-1,Np,4,4);
CaBPHI.resize(NR-1,NZ-1,Np,4,4);
CaBZ.resize(NR-1,NZ-1,Np,4,4);

// Get the derivatives
for(i=1;i<=Np;i++)
{
	slice.reference(BR(all,all,i));
	slicedR.reference(dBRdR(all,all,i));
	slicedZ.reference(dBRdZ(all,all,i));
	sliced2.reference(d2BR(all,all,i));
	bcuderiv(slice,dR,dZ,slicedR,slicedZ,sliced2);

	slice.reference(BPHI(all,all,i));
	slicedR.reference(dBPHIdR(all,all,i));
	slicedZ.reference(dBPHIdZ(all,all,i));
	sliced2.reference(d2BPHI(all,all,i));
	bcuderiv(slice,dR,dZ,slicedR,slicedZ,sliced2);

	slice.reference(BZ(all,all,i));
	slicedR.reference(dBZdR(all,all,i));
	slicedZ.reference(dBZdZ(all,all,i));
	sliced2.reference(d2BZ(all,all,i));
	bcuderiv(slice,dR,dZ,slicedR,slicedZ,sliced2);
}

// Get the c's for bcuint, as done by bcucof
for(i=1;i<NR;i++)
{
	for(j=1;j<NZ;j++)
	{
		for(k=1;k<=Np;k++)
		{
			y_sq(1) = BR(i,j,k); y_sq(2) = BR(i+1,j,k); y_sq(3) = BR(i+1,j+1,k); y_sq(4) = BR(i,j+1,k);
			y1_sq(1) = dBRdR(i,j,k); y1_sq(2) = dBRdR(i+1,j,k); y1_sq(3) = dBRdR(i+1,j+1,k); y1_sq(4) = dBRdR(i,j+1,k);
			y2_sq(1) = dBRdZ(i,j,k); y2_sq(2) = dBRdZ(i+1,j,k); y2_sq(3) = dBRdZ(i+1,j+1,k); y2_sq(4) = dBRdZ(i,j+1,k);
			y12_sq(1) = d2BR(i,j,k); y12_sq(2) = d2BR(i+1,j,k); y12_sq(3) = d2BR(i+1,j+1,k); y12_sq(4) = d2BR(i,j+1,k);
			slice.reference(CaBR(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,slice);

			y_sq(1) = BPHI(i,j,k); y_sq(2) = BPHI(i+1,j,k); y_sq(3) = BPHI(i+1,j+1,k); y_sq(4) = BPHI(i,j+1,k);
			y1_sq(1) = dBPHIdR(i,j,k); y1_sq(2) = dBPHIdR(i+1,j,k); y1_sq(3) = dBPHIdR(i+1,j+1,k); y1_sq(4) = dBPHIdR(i,j+1,k);
			y2_sq(1) = dBPHIdZ(i,j,k); y2_sq(2) = dBPHIdZ(i+1,j,k); y2_sq(3) = dBPHIdZ(i+1,j+1,k); y2_sq(4) = dBPHIdZ(i,j+1,k);
			y12_sq(1) = d2BPHI(i,j,k); y12_sq(2) = d2BPHI(i+1,j,k); y12_sq(3) = d2BPHI(i+1,j+1,k); y12_sq(4) = d2BPHI(i,j+1,k);
			slice.reference(CaBPHI(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,slice);

			y_sq(1) = BZ(i,j,k); y_sq(2) = BZ(i+1,j,k); y_sq(3) = BZ(i+1,j+1,k); y_sq(4) = BZ(i,j+1,k);
			y1_sq(1) = dBZdR(i,j,k); y1_sq(2) = dBZdR(i+1,j,k); y1_sq(3) = dBZdR(i+1,j+1,k); y1_sq(4) = dBZdR(i,j+1,k);
			y2_sq(1) = dBZdZ(i,j,k); y2_sq(2) = dBZdZ(i+1,j,k); y2_sq(3) = dBZdZ(i+1,j+1,k); y2_sq(4) = dBZdZ(i,j+1,k);
			y12_sq(1) = d2BZ(i,j,k); y12_sq(2) = d2BZ(i+1,j,k); y12_sq(3) = d2BZ(i+1,j+1,k); y12_sq(4) = d2BZ(i,j+1,k);
			slice.reference(CaBZ(i,j,k,all,all));
			bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,slice);
		}
	}
}
}

//------------------------- interpolate_B ---------------------------------------------------------------------------------
// gets br(R,Z), bphi(R,Z), bz(R,Z) in the phi-plane (k-1)*dp through bicubic spline interpolation
void MGRID::interpolate_B(double r, double z, int k, double& br, double& bphi, double& bz)
{
double dummy;
Array<double,4> slice;
Range all = Range::all();

slice.reference(CaBR(all,all,k,all,all));
bcuint(R,Z,slice,dR,dZ,r,z,br,dummy,dummy);

slice.reference(CaBPHI(all,all,k,all,all));
bcuint(R,Z,slice,dR,dZ,r,z,bphi,dummy,dummy);

slice.reference(CaBZ(all,all,k,all,all));
bcuint(R,Z,slice,dR,dZ,r,z,bz,dummy,dummy);
}

//------------------------ End of Class MGRID -----------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Begin Class POTENTIAL -----------------------------------------------------------------------------------------
class POTENTIAL
{
private:
	// Member Variables
	double epsabs; 			// absolute accuracy tolerance in adaptive Simpson, default = 1e-6
	double epsrel; 			// relative accuracy tolerance in adaptive Simpson, default = 1e-4
	int maxRecDepth; 		// max recursion in adaptive Simpson, defaut = 14
	MGRID mgrid;
	VMEC wout;
	double ctor_fac;		// 1e-7 * wout.ctor * Raxis
	double Raxis;			// average Raxis of VMEC
	double Zaxis;			// average Zaxis of VMEC

	double du,dv;
	Array<double,1> U;
	Array<double,1> V;
	Array<double,4> CaRs, CaZs, CaPot;

	//int count, count2;	// counter for Simpson integration

	// Member-Functions
	Array<double,1> biotSavartIntegs(double v2, double Rs, double v, double Zs);
	Array<double,1> integ_v(double u, int N);
	Array<double,1> integ_u(double v, int N);
	Array<double,1> normal(double dRdu, double dRdv, double dZdu, double dZdv);
	double green(double Rs, double v, double Zs);
	Array<double,1> integs(double u, double v, int N);
	Array<double,1> adaptiveSimpson(int flag, double args[], int N, double a, double b, double epsabs, double epsrel, int maxRecursionDepth);
	Array<double,1> adaptiveSimpsonsAux(int flag, double args[], int N, double a, double b, double epsabs, double epsrel, Array<double,1>&  S,
										Array<double,1>& fa, Array<double,1>& fb, Array<double,1>& fc, int bottom);
	void prep_plasmaCurrentB(void);							// compute plasmaCurrentB on MGRID and add to MGRID B-field
	void prep_sInterpolation(int Nu = 200, int Nv = 200);	// get s = 1 surface and prepare interpolation
	void interpolate_RZ(double u, double v, double& Rs, double& Zs, double& dRdu, double& dRdv, double& dZdu, double& dZdv);
	double interpolate_pot(double u, double v);

public:
	// Member Variables
	double R;		// Point where to evaluate
	double phi;
	double Z;

	// Constructors
	POTENTIAL();						// Default Constructor
	POTENTIAL(const POTENTIAL& pot);	// Copy Constructor
	POTENTIAL(VMEC woutin);				// Standard Constructor, uses defaults
	POTENTIAL(VMEC woutin, double epsabsin, double epsrelin, int maxRecDepthin); // set all Constructor

	// Member-Operators
	POTENTIAL& operator =(const POTENTIAL& pot);	// Operator =

	// Member-Functions
	void init(VMEC woutin);												// load mgrid and prepare all interpolations
	void setup_simps(double epsabsin = 1e-6, double epsrelin = 1e-4, int maxRecDepthin = 14);	// set control parameter for adaptive Simpson
	double ev(double R, double phi, double Z);							// evaluate scalar potential
	Array<double,1> grad(double Rin, double phiin, double Zin);			// evaluate Bplasma = gradPHI; returns Bplasma = (BR,Bphi,BZ)
	Array<double,1> get_vacuumB(double Rs, double v, double Zs);		// returns vacuum B-field Bvac = Bmgrid + Bip = (BR,Bphi,BZ)
	Array<double,1> get_plasmaCurrentB(double Rs, double v, double Zs);	// returns B-field from plasma current Bip = (BR,Bphi,BZ)
}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
POTENTIAL::POTENTIAL()
{
R = 0;
phi = 0;
Z = 0;
setup_simps();
}

// Copy Constructor
POTENTIAL::POTENTIAL(const POTENTIAL& pot)
{
*this = pot;
}

// Standard Constructor
POTENTIAL::POTENTIAL(VMEC woutin)
{
R = 0;
phi = 0;
Z = 0;
setup_simps();	// default values
init(woutin);
}

// set all Constructor
POTENTIAL::POTENTIAL(VMEC woutin, double epsabsin, double epsrelin, int maxRecDepthin)
{
R = 0;
phi = 0;
Z = 0;
setup_simps(epsabsin, epsrelin, maxRecDepthin);
init(woutin);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
POTENTIAL& POTENTIAL::operator =(const POTENTIAL& pot)
{
if (this == &pot) return(*this);	    // if: x=x
epsabs = pot.epsabs;
epsrel = pot.epsrel;
maxRecDepth = pot.maxRecDepth;
wout = pot.wout;
mgrid = pot.mgrid;
ctor_fac = pot.ctor_fac;
Raxis = pot.Raxis;
Zaxis = pot.Zaxis;
R = pot.R;
phi = pot.phi;
Z = pot.Z;

du = pot.du;
dv = pot.dv;
U.reference(pot.U);
V.reference(pot.V);
CaRs.reference(pot.CaRs);
CaZs.reference(pot.CaZs);
CaPot.reference(pot.CaPot);

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- read ---------------------------------------------------------------------------------------------------------------
void POTENTIAL::init(VMEC woutin)
{
int i;
int N = 50;
double Ra, Za, v;
wout = woutin;
mgrid.read(wout.mgrid_file, wout.nextcur, wout.extcur);

// get average magnetic axis for VMEC
Raxis = 0; Zaxis = 0;
for(i=0;i<N;i++)
{
	v = i*pi2/N;
	wout.get_axis(v, Ra, Za);
	Raxis += Ra;
	Zaxis += Za;
}
Raxis /= N;
Zaxis /= N;
ctor_fac = 1e-7 * wout.ctor * Raxis;	// factor for currentB -field

// compute plasmaCurrentB on MGRID and add to MGRID B-field
//cout << "prepare plasmaCurrentB" << endl;
prep_plasmaCurrentB();

// get s = 1 surface and prepare interpolation
//cout << "prepare Interpolations" << endl;
prep_sInterpolation();
}

// --- set_simps ----------------------------------------------------------------------------------------------------------
void POTENTIAL::setup_simps(double epsabsin, double epsrelin, int maxRecDepthin)
{
epsabs = epsabsin;
epsrel = epsrelin;
maxRecDepth = maxRecDepthin;
}

// --- ev -----------------------------------------------------------------------------------------------------------------
double POTENTIAL::ev(double Rin, double phiin, double Zin)
{
Array<double,1> integ(2);
double args[1];
R = Rin; phi = phiin; Z = Zin; 		// load into member variables
//****************
//count = 0; count2 = 0;
//****************
// integ = adaptiveSimpson(1, args, 2, 0, pi2, epsabs, epsrel, maxRecDepth);
integ = adaptiveSimpson(4, args, 2, 0, pi2, epsabs, epsrel, maxRecDepth);
//****************
//cout << "total function calls: " << count << "\t" << "1D integrals: " << count2 << endl;
//****************
return -(integ(0) + integ(1))/pi2;
}

// --- grad ---------------------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::grad(double Rin, double phiin, double Zin)
{
Array<double,1> integ(3);
double args[1];
R = Rin; phi = phiin; Z = Zin; 		// load into member variables
//****************
//count = 0; count2 = 0;
//****************
// integ = adaptiveSimpson(1, args, 3, 0, pi2, epsabs, epsrel, maxRecDepth);
integ = adaptiveSimpson(4, args, 3, 0, pi2, epsabs, epsrel, maxRecDepth);
integ(1) /= R;
integ /= -pi2;
//****************
//cout << "total function calls: " << count << "\t" << "1D integrals: " << count2 << endl;
//****************
return integ;
}

// --- get_vacuumB --------------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::get_vacuumB(double Rs, double v, double Zs)
{
double dphi = pi2 / mgrid.Np;
int k = (int(round(v/dphi)) % mgrid.Np) + 1;	// nearest neighbor approximation
Array<double,1> B(3);
double br, bphi,bz;

mgrid.interpolate_B(Rs, Zs, k, br, bphi, bz);
B(0) = br;
B(1) = bphi;
B(2) = bz;
return B;
}

// --- get_plasmaCurrentB -------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::get_plasmaCurrentB(double Rs, double v, double Zs)
{
Array<double,1> B(3), integ(2);
double args[3] = {Rs,v,Zs};
integ = adaptiveSimpson(2, args, 2, 0, pi2, epsabs, epsrel, maxRecDepth);
B(0) = ctor_fac * (Zs - Zaxis)*integ(0);
B(1) = 0;
B(2) = ctor_fac * (Raxis*integ(1) - Rs*integ(0));
return B;
}

//--------------------- Private Member Functions --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- integ_u/v ----------------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::integ_v(double u, int N)
{
//****************
//count2 += 1;
//****************
double args[1] = {u};
return adaptiveSimpson(0, args, N, 0, pi2, epsabs, epsrel, maxRecDepth);
}

// --- integ_u ------------------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::integ_u(double v, int N)
{
//****************
//count2 += 1;
//****************
double args[1] = {v};
return adaptiveSimpson(3, args, N, 0, pi2, epsabs, epsrel, maxRecDepth);
}

// --- integs -------------------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::integs(double u, double v, int N)
{
double Rs, Zs, G, pot;
double nx,nB,G3,G3p,nxG5p_nBG3, cospv, sinpv;
double dZdu, dZdv, dRdu, dRdv;
Array<double,1> out(N);
Array<double,1> n(3), B(3);
//***************
//count +=1;
//***************

// point on s = 1 surface
//double dummy;
//Rs = wout.rmn.ev(1.0, u, v, dummy, dRdu, dRdv);
//Zs = wout.zmn.ev(1.0, u, v, dummy, dZdu, dZdv);
interpolate_RZ(u, v, Rs, Zs, dRdu, dRdv, dZdu, dZdv);

// normal vector at this point
n = normal(dRdu, dRdv, dZdu, dZdv);

// evaluate Green's function
G = green(Rs, v, Zs);

// scalar potential on s = 1 surface
//pot = wout.pot(u, v);
pot = interpolate_pot(u, v);

// vacuum magnetic field on s = 1 surface
B = get_vacuumB(Rs, v, Zs);

sinpv = sin(phi-v);
cospv = cos(phi-v);
nx = n(0)*(R*cospv - Rs) + n(1)*R*sinpv + n(2)*(Z - Zs);
nB = n(0)*B(0) + n(1)*B(1) + n(2)*B(2);
G3 = G*G*G;
G3p = G3*pot;

if(N == 2)
{
	out(0) = nx * G3p;
	out(1) = nB * G;
}
if(N == 3)
{
	nxG5p_nBG3 = -1.5*nx*G3p*G*G - 0.5*nB*G3;
	out(0) = (n(0)*cospv + n(1)*sinpv)*G3p + 2*(R - Rs*cospv)*nxG5p_nBG3;
	out(1) = R*(-n(0)*sinpv + n(1)*cospv)*G3p + 2*R*Rs*sinpv*nxG5p_nBG3;
	out(2) = n(2)*G3p + 2*(Z - Zs)*nxG5p_nBG3;
}
return out;
}

// --- normal -------------------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::normal(double drdu, double drdv, double dzdu, double dzdv)
{
Array<double,1> n(3);
n(0) = -dzdu;
n(1) = dzdu*drdv - drdu*dzdv;
n(2) = drdu;

n /= sqrt(n(0)*n(0) + n(1)*n(1) + n(2)*n(2));
return n;
}

// --- green --------------------------------------------------------------------------------------------------------------
double POTENTIAL::green(double Rs, double v, double Zs)
{
return 1.0/sqrt(R*R + Rs*Rs - 2*R*Rs*cos(phi-v) + (Z - Zs)*(Z - Zs));
}

// --- biotSavartIntegs ---------------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::biotSavartIntegs(double v2, double Rs, double v, double Zs)
{
double cosv, Zm, d2, d3;
Array<double,1> out(2);
Zm = Zs - Zaxis;
cosv = cos(v - v2);
d2 = Rs*Rs + Raxis*Raxis - 2*Rs*Raxis*cosv + Zm*Zm;
d3 = d2*sqrt(d2);
out(0) = cosv/d3;
out(1) = 1.0/d3;
return out;
}

//--- Adaptive Simpson's Rule Recursor ------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::adaptiveSimpsonsAux(int flag, double args[], int N, double a, double b, double epsabs, double epsrel,
									Array<double,1>& S, Array<double,1>& fa, Array<double,1>& fb, Array<double,1>& fc, int bottom)
{
double c = (a + b)/2.0;
double h = (b - a)/12.0;
double d = (a + c)/2.0;
double e = (c + b)/2.0;

Array<double,1> fd(N), fe(N), out(N), Sleft(N), Sright(N), Snew(N);
switch(flag)
{
case 0:	// integs(u, x)
	fd = integs(args[0], d, N);
	fe = integs(args[0], e, N);
	break;
case 1: // integ_v(x)
	fd = integ_v(d, N);
	fe = integ_v(e, N);
	break;
case 2:	// biotSavartIntegs(x, Rs, v, Zs)
	fd = biotSavartIntegs(d, args[0], args[1], args[2]);
	fe = biotSavartIntegs(e, args[0], args[1], args[2]);
	break;
case 3:	// integs(x, v)
	fd = integs(d, args[0], N);
	fe = integs(e, args[0], N);
	break;
case 4: // integ_u(x)
	fd = integ_u(d, N);
	fe = integ_u(e, N);
	break;
}

Sleft = h*(fa + 4*fd + fc);
Sright = h*(fc + 4*fe + fb);
Snew = Sleft + Sright;

if (bottom <= 0 || (max(fabs(Snew - S)) <= 15*epsabs || max(fabs(Snew - S)/fabs(Snew)) <= 15*epsrel))
{
	out = Snew + (Snew - S)/15.0;
	return out;
}
out = adaptiveSimpsonsAux(flag, args, N, a, c, epsabs/2, epsrel/2, Sleft, fa, fc, fd, bottom-1) + adaptiveSimpsonsAux(flag, args, N, c, b, epsabs/2, epsrel/2, Sright, fc, fb, fe, bottom-1);
return out;
}

//--- Adaptive Simpson's Rule ---------------------------------------------------------------------------------------------
Array<double,1> POTENTIAL::adaptiveSimpson(int flag, double args[], int N, double a, double b, double epsabs, double epsrel, int maxRecursionDepth)
{
double c = (a + b)/2;
double h = (b - a)/6.0;

Array<double,1> fa(N),fb(N),fc(N),out(N),S(N);
switch(flag)
{
case 0:	// integs(u, x)
	fa = integs(args[0], a, N);
	fb = integs(args[0], b, N);
	fc = integs(args[0], c, N);
	break;
case 1: // integ_v(x)
	fa = integ_v(a, N);
	fb = integ_v(b, N);
	fc = integ_v(c, N);
	break;
case 2:	// biotSavartIntegs(x, Rs, v, Zs)
	fa = biotSavartIntegs(a, args[0], args[1], args[2]);
	fb = biotSavartIntegs(b, args[0], args[1], args[2]);
	fc = biotSavartIntegs(c, args[0], args[1], args[2]);
	break;
case 3:	// integs(x, v)
	fa = integs(a, args[0], N);
	fb = integs(b, args[0], N);
	fc = integs(c, args[0], N);
	break;
case 4: // integ_u(x)
	fa = integ_u(a, N);
	fb = integ_u(b, N);
	fc = integ_u(c, N);
	break;
}
S = h*(fa + 4*fc + fb);

return adaptiveSimpsonsAux(flag, args, N, a, b, epsabs, epsrel, S, fa, fb, fc, maxRecursionDepth);
}

// --- prep_plasmaCurrentB ------------------------------------------------------------------------------------------------
void POTENTIAL::prep_plasmaCurrentB(void)
{
Array<double,3> BR(Range(1,mgrid.NR), Range(1,mgrid.NZ), Range(1,mgrid.Np));
Array<double,3> BPHI(Range(1,mgrid.NR), Range(1,mgrid.NZ), Range(1,mgrid.Np));
Array<double,3> BZ(Range(1,mgrid.NR), Range(1,mgrid.NZ), Range(1,mgrid.Np));
//#pragma omp parallel
//{
	Array<double,1> B(3);
//	#pragma omp for collapse(3)
	for(int i=1;i<=mgrid.NR;i++)
		for(int j=1;j<=mgrid.NZ;j++)
			for(int k=1;k<=mgrid.Np;k++)
			{
				B = get_plasmaCurrentB(mgrid.R(i), (k-1)*pi2/mgrid.Np, mgrid.Z(j));
				BR(i,j,k) = B(0);
				BPHI(i,j,k) = B(1);
				BZ(i,j,k) = B(2);
			}
//}
mgrid.BR += BR;
mgrid.BPHI += BPHI;
mgrid.BZ += BZ;
mgrid.prep_interpolation();
}

// --- prep_sInterpolation ------------------------------------------------------------------------------------------------
void POTENTIAL::prep_sInterpolation(int Nu, int Nv)
{
int i,j;
double dummy;
Array<double,2> Rs(Range(1,Nu), Range(1,Nv)), Zs(Range(1,Nu), Range(1,Nv)), Pot(Range(1,Nu), Range(1,Nv));
Array<double,2> dRdu(Range(1,Nu), Range(1,Nv)), dZdu(Range(1,Nu), Range(1,Nv));
Array<double,2> dRdv(Range(1,Nu), Range(1,Nv)), dZdv(Range(1,Nu), Range(1,Nv));
Array<double,2> d2R(Range(1,Nu), Range(1,Nv)), d2Z(Range(1,Nu), Range(1,Nv));
Array<double,2> dPdu(Range(1,Nu), Range(1,Nv)), dPdv(Range(1,Nu), Range(1,Nv)), d2P(Range(1,Nu), Range(1,Nv));
Range all = Range::all();

TinyVector <int,1> index(1);
TinyVector <int,4> index4(1,1,1,1);
U.resize(Nu);	U.reindexSelf(index);
V.resize(Nv);	V.reindexSelf(index);
CaRs.resize(Nu-1,Nv-1,4,4);		CaRs.reindexSelf(index4);
CaZs.resize(Nu-1,Nv-1,4,4);		CaZs.reindexSelf(index4);
CaPot.resize(Nu-1,Nv-1,4,4);	CaPot.reindexSelf(index4);

// s = 1 surface
du = pi2/(Nu-1);
dv = pi2/(Nv-1);
for(i=1;i<=Nu;i++) U(i) = (i-1)*du;
for(j=1;j<=Nv;j++) V(j) = (j-1)*dv;

//#pragma omp parallel for private(i,j) collapse(2)
for(i=1;i<=Nu;i++)
	for(j=1;j<=Nv;j++)
	{
		Rs(i,j) = wout.rmn.ev(1.0, U(i), V(j), dummy, dRdu(i,j), dRdv(i,j), d2R(i,j));
		Zs(i,j) = wout.zmn.ev(1.0, U(i), V(j), dummy, dZdu(i,j), dZdv(i,j), d2Z(i,j));
		Pot(i,j) = wout.pot(U(i), V(j), dPdu(i,j), dPdv(i,j), d2P(i,j));
	}

// get the C's
Array<double,2> slice;
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
for(i=1;i<Nu;i++)
{
	for(j=1;j<Nv;j++)
	{
		y_sq(1) = Rs(i,j); y_sq(2) = Rs(i+1,j); y_sq(3) = Rs(i+1,j+1); y_sq(4) = Rs(i,j+1);
		y1_sq(1) = dRdu(i,j); y1_sq(2) = dRdu(i+1,j); y1_sq(3) = dRdu(i+1,j+1); y1_sq(4) = dRdu(i,j+1);
		y2_sq(1) = dRdv(i,j); y2_sq(2) = dRdv(i+1,j); y2_sq(3) = dRdv(i+1,j+1); y2_sq(4) = dRdv(i,j+1);
		y12_sq(1) = d2R(i,j); y12_sq(2) = d2R(i+1,j); y12_sq(3) = d2R(i+1,j+1); y12_sq(4) = d2R(i,j+1);
		slice.reference(CaRs(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,du,dv,slice);

		y_sq(1) = Zs(i,j); y_sq(2) = Zs(i+1,j); y_sq(3) = Zs(i+1,j+1); y_sq(4) = Zs(i,j+1);
		y1_sq(1) = dZdu(i,j); y1_sq(2) = dZdu(i+1,j); y1_sq(3) = dZdu(i+1,j+1); y1_sq(4) = dZdu(i,j+1);
		y2_sq(1) = dZdv(i,j); y2_sq(2) = dZdv(i+1,j); y2_sq(3) = dZdv(i+1,j+1); y2_sq(4) = dZdv(i,j+1);
		y12_sq(1) = d2Z(i,j); y12_sq(2) = d2Z(i+1,j); y12_sq(3) = d2Z(i+1,j+1); y12_sq(4) = d2Z(i,j+1);
		slice.reference(CaZs(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,du,dv,slice);

		y_sq(1) = Pot(i,j); y_sq(2) = Pot(i+1,j); y_sq(3) = Pot(i+1,j+1); y_sq(4) = Pot(i,j+1);
		y1_sq(1) = dPdu(i,j); y1_sq(2) = dPdu(i+1,j); y1_sq(3) = dPdu(i+1,j+1); y1_sq(4) = dPdu(i,j+1);
		y2_sq(1) = dPdv(i,j); y2_sq(2) = dPdv(i+1,j); y2_sq(3) = dPdv(i+1,j+1); y2_sq(4) = dPdv(i,j+1);
		y12_sq(1) = d2P(i,j); y12_sq(2) = d2P(i+1,j); y12_sq(3) = d2P(i+1,j+1); y12_sq(4) = d2P(i,j+1);
		slice.reference(CaPot(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,du,dv,slice);
	}
}
}

//------------------------- interpolate_RZ --------------------------------------------------------------------------------
// gets Rs(u,v), Zs(u,v) on the s = 1 surface through bicubic spline interpolation
void POTENTIAL::interpolate_RZ(double u, double v, double& Rs, double& Zs, double& dRdu, double& dRdv, double& dZdu, double& dZdv)
{
bcuint(U,V,CaRs,du,dv,u,v,Rs,dRdu,dRdv);
bcuint(U,V,CaZs,du,dv,u,v,Zs,dZdu,dZdv);
}

//------------------------- interpolate_pot -------------------------------------------------------------------------------
// gets potential(u,v) on the s = 1 surface through bicubic spline interpolation
double POTENTIAL::interpolate_pot(double u, double v)
{
double dummy, pot;
bcuint(U,V,CaPot,du,dv,u,v,pot,dummy,dummy);
return pot;
}

//------------------------ End of Class POTENTIAL -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  EXTENDER_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

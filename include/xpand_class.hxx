// Class calculates B-field outside of VMEC boundary
// A.Wingen						9.10.14


// Define
//--------
#ifndef XPAND_CLASS_INCLUDED
#define XPAND_CLASS_INCLUDED

// Include
//--------
#ifdef USE_MPI
#ifdef USE_MPICH
	#include <mpi.h>
#else
	#include <openmpi/ompi/mpi/cxx/mpicxx.h>
#endif
#endif
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <netcdf.h>
#include <andi.hxx>
#include <splines.hxx>
#include <vmec_class.hxx>
using namespace blitz;

// Prototypes
//-----------
void biot_savart(Array<double,2>& xc, Array<double,2>& dv, double I, Array<double,1>& x, Array<double,1>& b);
int read_extcur(Array<double,1>& extcur, LA_STRING input_extension);

// Golbal Parameters
//------------------


//--------- Begin Class MGRID ---------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Calculate vacuum magnetic field from MGRID file
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class MGRID
{
private:
	// Member Variables
	LA_STRING mgrid_file;		// pathname of MGRID file
	Array<double,1> extcur;		// coil currents 1->Ncoils
	double dR, dZ, dp;
	Array<double,3> dBRdR, dBPHIdR, dBZdR;
	Array<double,3> dBRdZ, dBPHIdZ, dBZdZ;
	Array<double,3> d2BR, d2BPHI, d2BZ;
	Array<double,5> CaBR, CaBPHI, CaBZ;
	Array<double,4> br, bp, bz;

	// Member-Functions

public:
	// Member Variables
	int Ncoils;					// number of coils in mgrid
	int NR, NZ, Np;
	Array<double,1> R;
	Array<double,1> Z;
	Array<double,3> BR;
	Array<double,3> BPHI;
	Array<double,3> BZ;

	// Constructors
	MGRID();								// Default Constructor
	MGRID(const MGRID& mgrid);				// Copy Constructor

	// Member-Operators
	MGRID& operator =(const MGRID& mgrid);	// Operator =

	// Member-Functions
	void read(LA_STRING file, int N);
	void read(LA_STRING file, int N, Array<double,1> cur);
	void makeB(Array<double,1> cur);
	void get_nearest(double r, double v, double z, int& ir, int& k, int& iz);
	void getB(int ir, int k, int iz, int icoil, double& Br, double& Bp, double& Bz);
	void prep_interpolation(void);	// initiate the interpolations
	void interpolate_B(double r, double z, int k, double& br, double& bphi, double& bz);	// interpolate all 3 vacuum B-field components on the mgrid grid
	Array<double,1> get_vacuumB(double Rs, double v, double Zs, bool n0only);
}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
MGRID::MGRID()
{
TinyVector <int,1> index(1);
TinyVector <int,3> index3(1,1,1);
TinyVector <int,4> index4(1,1,1,1);
TinyVector <int,5> index5(1,1,1,1,1);

mgrid_file = "None";
Ncoils = 36;
extcur.resize(Ncoils); extcur.reindexSelf(index);

NR = 129;
NZ = 129;
Np = 48;

dR = 0.001;
dZ = 0.001;
dp = pi2/Np;				// 0 -> 2*pi

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

br.resize(Np, NZ, NR, Ncoils);		br.reindexSelf(index4);
bp.resize(Np, NZ, NR, Ncoils);		bp.reindexSelf(index4);
bz.resize(Np, NZ, NR, Ncoils);		bz.reindexSelf(index4);
}

// Copy Constructor
// makes a true copy
MGRID::MGRID(const MGRID& mgrid)
{
mgrid_file = mgrid.mgrid_file;
Ncoils = mgrid.Ncoils;
extcur.reference(mgrid.extcur.copy());

NR = mgrid.NR;
NZ = mgrid.NZ;
Np = mgrid.Np;

dR = mgrid.dR;
dZ = mgrid.dZ;
dp = mgrid.dp;

R.reference(mgrid.R.copy());
Z.reference(mgrid.Z.copy());

BR.reference(mgrid.BR.copy());
BPHI.reference(mgrid.BPHI.copy());
BZ.reference(mgrid.BZ.copy());

dBRdR.reference(mgrid.dBRdR.copy());
dBPHIdR.reference(mgrid.dBPHIdR.copy());
dBZdR.reference(mgrid.dBZdR.copy());

dBRdZ.reference(mgrid.dBRdZ.copy());
dBPHIdZ.reference(mgrid.dBPHIdZ.copy());
dBZdZ.reference(mgrid.dBZdZ.copy());

d2BR.reference(mgrid.d2BR.copy());
d2BPHI.reference(mgrid.d2BPHI.copy());
d2BZ.reference(mgrid.d2BZ.copy());

CaBR.reference(mgrid.CaBR.copy());
CaBPHI.reference(mgrid.CaBPHI.copy());
CaBZ.reference(mgrid.CaBZ.copy());

br.reference(mgrid.br.copy());
bp.reference(mgrid.bp.copy());
bz.reference(mgrid.bz.copy());
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

br.reference(mgrid.br);
bp.reference(mgrid.bp);
bz.reference(mgrid.bz);

return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- read ---------------------------------------------------------------------------------------------------------------
void MGRID::read(LA_STRING file, int N)
{
// Variables
int i,chk, ncid, varid;
int k;
stringstream ss;
string s;
LA_STRING varname;
Range all = Range::all();

// assign private members
mgrid_file = file;
Ncoils = N;

// open file
chk = nc_open(mgrid_file, NC_NOWRITE, &ncid);
if(chk != 0) {cout << "MGRID: Unable to open file " << mgrid_file << endl; EXIT;}

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
dR = (rmax - rmin)/double(NR-1);
dZ = (zmax - zmin)/double(NZ-1);
dp = pi2/double(Np);	// Np + 1 surface is the phi = 2pi surface
for(k=1;k<=NR;k++) R(k) = rmin + (k-1)*dR;
for(k=1;k<=NZ;k++) Z(k) = zmin + (k-1)*dZ;

// read B-field data
Array<double,3> input;
input.resize(Np, NZ, NR);
br.resize(Np, NZ, NR, Ncoils);
bp.resize(Np, NZ, NR, Ncoils);
bz.resize(Np, NZ, NR, Ncoils);

for (i=1;i<=Ncoils;i++)
{
	ss << setw(3) << setfill('0') << i; s = ss.str();

	varname = "br_" + LA_STRING(s.c_str());
	chk = nc_inq_varid(ncid, varname, &varid);
	chk = nc_get_var_double(ncid, varid, input.data());
	br(all,all,all,i) = input.copy();

	varname = "bp_" + LA_STRING(s.c_str());
	chk = nc_inq_varid(ncid, varname, &varid);
	chk = nc_get_var_double(ncid, varid, input.data());
	bp(all,all,all,i) = input.copy();

	varname = "bz_" + LA_STRING(s.c_str());
	chk = nc_inq_varid(ncid, varname, &varid);
	chk = nc_get_var_double(ncid, varid, input.data());
	bz(all,all,all,i) = input.copy();

	// reset stringstream
	ss.str("");
	ss.clear();
}
}

// ---------------------------------------------------------
void MGRID::read(LA_STRING file, int N, Array<double,1> cur)
{
read(file,N);
makeB(cur);
}

// --- makeB --------------------------------------------------------------------------------------------------------------
void MGRID::makeB(Array<double,1> cur)
{
// Variables
int i;
Range all = Range::all();

extcur.resize(Ncoils);
extcur.reference(cur);

BR.resize(Np, NZ, NR); 		BR = 0;
BPHI.resize(Np, NZ, NR);	BPHI = 0;
BZ.resize(Np, NZ, NR);		BZ = 0;

for (i=1;i<=Ncoils;i++)
{
	if(extcur(i) == 0) continue;
	BR += extcur(i) * br(all,all,all,i);
	BPHI += extcur(i) * bp(all,all,all,i);
	BZ += extcur(i) * bz(all,all,all,i);
}

// make B-field arrays: (NR, NZ, Np)
BR.transposeSelf(thirdDim, secondDim, firstDim);
BPHI.transposeSelf(thirdDim, secondDim, firstDim);
BZ.transposeSelf(thirdDim, secondDim, firstDim);

//prep_interpolation();
}

// --- getB ---------------------------------------------------------------------------------------------------------------
void MGRID::get_nearest(double r, double v, double z, int& ir, int& k, int& iz)
{
ir = int((r - R(1))/dR) + 1;
if (fabs(R(ir+1) - r) < fabs(R(ir) - r)) ir += 1;
iz = int((z - Z(1))/dZ) + 1;
if (fabs(Z(iz+1) - z) < fabs(Z(iz) - z)) iz += 1;
k = int(v/dp) + 1;
if (fabs(k*dp - v) < fabs((k-1)*dp -v)) k += 1;
k = ((k-1) % Np) + 1;
}

void MGRID::getB(int ir, int k, int iz, int icoil, double& Br, double& Bp, double& Bz)
{
Br = br(k,iz,ir,icoil);
Bp = bp(k,iz,ir,icoil);
Bz = bz(k,iz,ir,icoil);
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


// --- get_vacuumB --------------------------------------------------------------------------------------------------------
Array<double,1> MGRID::get_vacuumB(double Rs, double v, double Zs, bool n0only)
{
Array<double,1> B(3);
double br, bphi, bz;
double bru, bphiu, bzu;
double t, dphi;
int k,ku;

if(n0only)	// compute B for all k*dphi planes and average
{
	B = 0;
	for(k=1;k<=Np;k++)
	{
		interpolate_B(Rs, Zs, k, br, bphi, bz);
		B(0) += br;
		B(1) += bphi;
		B(2) += bz;
	}
	B /= Np;
}
else	// use linear interpolation between k and k+1 planes
{
	v = modulo2pi(v);	// make sure v in [0, 2pi]
	dphi = pi2 / Np;
	k = int(v/dphi) + 1;	// k*dphi plane
	//k = (int(round(v/dphi)) % mgrid.Np) + 1;	// nearest neighbor approximation
	t = v/dphi - k + 1;

	interpolate_B(Rs, Zs, k, br, bphi, bz);
	B(0) = br;
	B(1) = bphi;
	B(2) = bz;
	if(t > 0)	// not exactly in the v-plane (k-1)*dphi
	{
		if(k == Np) ku = 1;
		else ku = k+1;
		interpolate_B(Rs, Zs, ku, bru, bphiu, bzu);

		// linear interpolation of B-field in v
		B(0) += t*(bru - br);
		B(1) += t*(bphiu - bphi);
		B(2) += t*(bzu - bz);
	}
}
return B;
}

//------------------------ End of Class MGRID -----------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//--------- Begin Class BFIELDVC ------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Calculate magnetic field outside of VMEC s = 1 using a virtual casing principle
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class BFIELDVC
{
private:
	// Member Variables
	double epsabs; 			// absolute accuracy tolerance, default = 1e-6
	double epsrel; 			// relative accuracy tolerance, default = 1e-4
	int maxRecDepth; 		// max recursion in adaptive Simpson, defaut = 14
	bool noVacuumField;		// True: ignore any vacuum fields -> Bvac = 0, False: use vacuum fields from Mgrid
	MGRID mgrid;
	VMEC wout;

	double du,dv;
	Array<double,1> U;
	Array<double,1> V;
	Array<double,4> CaRs, CaZs, CaKR, CaKp, CaKZ, CanB;

	// Member-Functions
	Array<double,1> integ_v(double u);
	Array<double,1> integ_u(double v);
	Array<double,1> normal(double r, double sinv, double cosv, double drdu, double drdv, double dzdu, double dzdv);	// carthesian
	double areal(double r, double drdu, double drdv, double dzdu, double dzdv);
	double green(Array<double,1>& x);
	Array<double,1> adaptiveSimpson(int flag, double args[], double a, double b, double epsabs, double epsrel, int maxRecursionDepth);
	Array<double,1> adaptiveSimpsonsAux(int flag, double args[], double a, double b, double epsabs, double epsrel, Array<double,1>&  S,
										Array<double,1>& fa, Array<double,1>& fb, Array<double,1>& fc, int bottom);
	void prep_sInterpolation(int Nu = 300, int Nv = 300);	// get s = 1 surface and prepare interpolation
	void interpolate_RZ(double u, double v, double& Rs, double& Zs);
	void interpolate_VCcur(double u, double v, double& KR, double& Kp, double& KZ);
	double interpolate_dipole(double u, double v);

public:
	// Member Variables
	double R;		// Point where to evaluate
	double phi;
	double Z;

	// Constructors
	BFIELDVC();							// Default Constructor
	BFIELDVC(const BFIELDVC& bvc);		// Copy Constructor
	BFIELDVC(VMEC woutin, LA_STRING mgrid_file = "None", bool useVacuum = true);				// Standard Constructor, uses defaults
	BFIELDVC(VMEC woutin, double epsabsin, double epsrelin, int maxRecDepthin, LA_STRING mgrid_file = "None", bool useVacuum = true); // set all Constructor

	// Member-Operators
	BFIELDVC& operator =(const BFIELDVC& bvc);	// Operator =

	// Member-Functions
	void init(VMEC woutin, LA_STRING mgrid_file = "None");				// load mgrid and prepare all interpolations
	void setup_accuracy(double epsabsin = 1e-6, double epsrelin = 1e-4, int maxRecDepthin = 14);	// set control parameter for adaptive integration
	void use_Vacuum(bool useVacuum);									// True: use Vacuum fields, False: ignore all Vacuum fields -> Bvac = 0
	Array<double,1> ev(double R, double phi, double Z);					// evaluate magentic field by virtual casing, uses adaptive Simpson integration
	Array<double,1> get_vacuumB(double Rs, double v, double Zs, bool n0only = false);		// returns vacuum B-field Bvac = Bmgrid = (BR,Bphi,BZ)
	Array<double,1> integs(double u, double v);							// integrad for all three integrals: BR, Bphi and BZ
}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
BFIELDVC::BFIELDVC()
{
R = 0;
phi = 0;
Z = 0;
noVacuumField = false;
setup_accuracy();
}

// Copy Constructor
BFIELDVC::BFIELDVC(const BFIELDVC& bvc)
{
*this = bvc;
}

// Standard Constructor
BFIELDVC::BFIELDVC(VMEC woutin, LA_STRING mgrid_file, bool useVacuum)
{
R = 0;
phi = 0;
Z = 0;
use_Vacuum(useVacuum);
setup_accuracy();	// default values
init(woutin, mgrid_file);
}

// set all Constructor
BFIELDVC::BFIELDVC(VMEC woutin, double epsabsin, double epsrelin, int maxRecDepthin, LA_STRING mgrid_file, bool useVacuum)
{
R = 0;
phi = 0;
Z = 0;
use_Vacuum(useVacuum);
setup_accuracy(epsabsin, epsrelin, maxRecDepthin);
init(woutin, mgrid_file);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
BFIELDVC& BFIELDVC::operator =(const BFIELDVC& bvc)
{
if (this == &bvc) return(*this);	    // if: x=x
epsabs = bvc.epsabs;
epsrel = bvc.epsrel;
maxRecDepth = bvc.maxRecDepth;
noVacuumField = bvc.noVacuumField;
wout = bvc.wout;
mgrid = bvc.mgrid;
R = bvc.R;
phi = bvc.phi;
Z = bvc.Z;

du = bvc.du;
dv = bvc.dv;
U.reference(bvc.U);
V.reference(bvc.V);
CaRs.reference(bvc.CaRs);
CaZs.reference(bvc.CaZs);
CaKR.reference(bvc.CaKR);
CaKp.reference(bvc.CaKp);
CaKZ.reference(bvc.CaKZ);
CanB.reference(bvc.CanB);

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- init ---------------------------------------------------------------------------------------------------------------
void BFIELDVC::init(VMEC woutin, LA_STRING mgrid_file)
{
#ifdef USE_MPI
int mpi_rank = MPI::COMM_WORLD.Get_rank();
#else
int mpi_rank = 0;
#endif
Array<double,1> extcur;
int nextcur;
wout = woutin;
if(not noVacuumField)
{
	if(strcmp(mgrid_file,"None") == 0) mgrid_file = wout.mgrid_file;
	if(not wout.lfreeb)
	{
		if(strcmp(mgrid_file,"None") == 0) {if(mpi_rank == 0) cout << "XPAND: provide mgrid file" << endl; EXIT;}
		nextcur = read_extcur(extcur, wout.input_extension);
		wout.nextcur = nextcur;
		wout.extcur.reference(extcur);
	}
	mgrid.read(mgrid_file, wout.nextcur, wout.extcur);
	mgrid.prep_interpolation();
}
// get s = 1 surface and prepare interpolation
prep_sInterpolation();
}


// --- set_accuracy -------------------------------------------------------------------------------------------------------
void BFIELDVC::setup_accuracy(double epsabsin, double epsrelin, int maxRecDepthin)
{
epsabs = epsabsin;
epsrel = epsrelin;
maxRecDepth = maxRecDepthin;	// Simpson only
}

// --- use_Vacuum -------------------------------------------------------------------------------------------------------
void BFIELDVC::use_Vacuum(bool useVacuum)
{
if(useVacuum) noVacuumField = false;
else noVacuumField = true;
}

// --- ev -----------------------------------------------------------------------------------------------------------------
Array<double,1> BFIELDVC::ev(double Rin, double phiin, double Zin)
{
Array<double,1> integ(3);
double args[1];
R = Rin; phi = phiin; Z = Zin; 		// load into member variables

integ = adaptiveSimpson(4, args, 0, pi2, epsabs, epsrel, maxRecDepth);
integ /= -4*pi;	// mu0/4pi and mu0*H = B; because K = n x H originally, so unit of K is actually mu0*current;   the negative sign is because VC returns the negative plasma response field

double Bx, By;
Bx = integ(0); By = integ(1);
integ(0) =  Bx*cos(phiin) + By*sin(phiin);
integ(1) = -Bx*sin(phiin) + By*cos(phiin);
return integ;
}


// --- get_vacuumB --------------------------------------------------------------------------------------------------------
Array<double,1> BFIELDVC::get_vacuumB(double Rs, double v, double Zs, bool n0only)
{
Array<double,1> B(3);
if(noVacuumField) B = 0;
else B = mgrid.get_vacuumB(Rs, v, Zs, n0only);
return B;
}


//--------------------- Private Member Functions --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- integ_v ------------------------------------------------------------------------------------------------------------
Array<double,1> BFIELDVC::integ_v(double u)
{
double args[1] = {u};
return adaptiveSimpson(0, args, 0, pi2, epsabs, epsrel, maxRecDepth);
}


// --- integ_u ------------------------------------------------------------------------------------------------------------
Array<double,1> BFIELDVC::integ_u(double v)
{
double args[1] = {v};
return adaptiveSimpson(3, args, 0, pi2, epsabs, epsrel, maxRecDepth);
}


// --- integs -------------------------------------------------------------------------------------------------------------
Array<double,1> BFIELDVC::integs(double u, double v)
{
double Rs, Zs, G, G3, nB;
Array<double,1> out(3);
Array<double,1> K(3), x(3);

// point on last closed flux surface (LCFS
interpolate_RZ(u, v, Rs, Zs);

// dipole moment
nB = interpolate_dipole(u, v);

// virtual casing current density on s = 1; K = n x B
interpolate_VCcur(u, v, K(0), K(1), K(2));

// \vec r - \vec r';  carthesian
x(0) = R*cos(phi) - Rs*cos(v);
x(1) = R*sin(phi) - Rs*sin(v);
x(2) = Z - Zs;

// virtual casing integrals
out(0) = K(1)*x(2) - K(2)*x(1) + nB*x(0);
out(1) = K(2)*x(0) - K(0)*x(2) + nB*x(1);
out(2) = K(0)*x(1) - K(1)*x(0) + nB*x(2);

// evaluate Green's function
G = green(x);
G3 = G*G*G;
out *= G3;

return out;
}


// --- normal -------------------------------------------------------------------------------------------------------------
// normal vector, \vec n = \partial/(\partial u) \vec r x \partial/(\partial v) \vec r, in carthesian coordinates
Array<double,1> BFIELDVC::normal(double r, double sinv, double cosv, double drdu, double drdv, double dzdu, double dzdv)
{
Array<double,1> n(3);
n(0) = drdu*dzdv*sinv - drdv*dzdu*sinv - dzdu*r*cosv;
n(1) = drdv*dzdu*cosv - dzdu*r*sinv - drdu*dzdv*cosv;
n(2) = drdu * r;

//n /= sqrt(n(0)*n(0) + n(1)*n(1) + n(2)*n(2));
// the length of n is also the area element, so set return from 'areal' to 1 and omit the normalization here
return n;
}


// --- areal -------------------------------------------------------------------------------------------------------------
double BFIELDVC::areal(double r, double drdu, double drdv, double dzdu, double dzdv)
{
return 1;
}


// --- green --------------------------------------------------------------------------------------------------------------
double BFIELDVC::green(Array<double,1>& x)
{
return 1.0/sqrt(x(0)*x(0) + x(1)*x(1) + x(2)*x(2));
}


//--- Adaptive Simpson's Rule Recursor ------------------------------------------------------------------------------------
Array<double,1> BFIELDVC::adaptiveSimpsonsAux(int flag, double args[], double a, double b, double epsabs, double epsrel,
									Array<double,1>& S, Array<double,1>& fa, Array<double,1>& fb, Array<double,1>& fc, int bottom)
{
int N = 3;
double c = (a + b)/2.0;
double h = (b - a)/12.0;
double d = (a + c)/2.0;
double e = (c + b)/2.0;

Array<double,1> fd(N), fe(N), out(N), Sleft(N), Sright(N), Snew(N);
switch(flag)
{
case 0:	// integs(u, x)
	fd = integs(args[0], d);
	fe = integs(args[0], e);
	break;
case 1: // integ_v(x)
	fd = integ_v(d);
	fe = integ_v(e);
	break;
case 3:	// integs(x, v)
	fd = integs(d, args[0]);
	fe = integs(e, args[0]);
	break;
case 4: // integ_u(x)
	fd = integ_u(d);
	fe = integ_u(e);
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
out = adaptiveSimpsonsAux(flag, args, a, c, epsabs/2, epsrel/2, Sleft, fa, fc, fd, bottom-1) + adaptiveSimpsonsAux(flag, args, c, b, epsabs/2, epsrel/2, Sright, fc, fb, fe, bottom-1);
return out;
}


//--- Adaptive Simpson's Rule ---------------------------------------------------------------------------------------------
Array<double,1> BFIELDVC::adaptiveSimpson(int flag, double args[], double a, double b, double epsabs, double epsrel, int maxRecursionDepth)
{
int N = 3;
double c = (a + b)/2;
double h = (b - a)/6.0;

Array<double,1> fa(N),fb(N),fc(N),out(N),S(N);
switch(flag)
{
case 0:	// integs(u, x)
	fa = integs(args[0], a);
	fb = integs(args[0], b);
	fc = integs(args[0], c);
	break;
case 1: // integ_v(x)
	fa = integ_v(a);
	fb = integ_v(b);
	fc = integ_v(c);
	break;
case 3:	// integs(x, v)
	fa = integs(a, args[0]);
	fb = integs(b, args[0]);
	fc = integs(c, args[0]);
	break;
case 4: // integ_u(x)
	fa = integ_u(a);
	fb = integ_u(b);
	fc = integ_u(c);
	break;
}
S = h*(fa + 4*fc + fb);

return adaptiveSimpsonsAux(flag, args, a, b, epsabs, epsrel, S, fa, fb, fc, maxRecursionDepth);
}


// --- prep_sInterpolation ------------------------------------------------------------------------------------------------
void BFIELDVC::prep_sInterpolation(int Nu, int Nv)
{
int i,j,k;
double dummy;
double bsupu,bsupv;
double sinv, cosv, BR, Bphi;
Array<double,1> sinuv, cosuv;
Array<double,1> n(3), B(3), Bvac(3);

Array<double,2> Rs(Range(1,Nu), Range(1,Nv)), Zs(Range(1,Nu), Range(1,Nv));
Array<double,2> dRdu(Range(1,Nu), Range(1,Nv)), dZdu(Range(1,Nu), Range(1,Nv));
Array<double,2> dRdv(Range(1,Nu), Range(1,Nv)), dZdv(Range(1,Nu), Range(1,Nv));
Array<double,2> d2R(Range(1,Nu), Range(1,Nv)), d2Z(Range(1,Nu), Range(1,Nv));

Array<double,2> KR(Range(1,Nu), Range(1,Nv)), Kp(Range(1,Nu), Range(1,Nv)), KZ(Range(1,Nu), Range(1,Nv));
Array<double,2> dKRdu(Range(1,Nu), Range(1,Nv)), dKpdu(Range(1,Nu), Range(1,Nv)), dKZdu(Range(1,Nu), Range(1,Nv));
Array<double,2> dKRdv(Range(1,Nu), Range(1,Nv)), dKpdv(Range(1,Nu), Range(1,Nv)), dKZdv(Range(1,Nu), Range(1,Nv));
Array<double,2> d2KR(Range(1,Nu), Range(1,Nv)), d2Kp(Range(1,Nu), Range(1,Nv)), d2KZ(Range(1,Nu), Range(1,Nv));

Array<double,2> nB(Range(1,Nu), Range(1,Nv));
Array<double,2> dnBdu(Range(1,Nu), Range(1,Nv));
Array<double,2> dnBdv(Range(1,Nu), Range(1,Nv));
Array<double,2> d2nB(Range(1,Nu), Range(1,Nv));

Range all = Range::all();

TinyVector <int,1> index(1);
TinyVector <int,4> index4(1,1,1,1);
U.resize(Nu);	U.reindexSelf(index);
V.resize(Nv);	V.reindexSelf(index);
CaRs.resize(Nu-1,Nv-1,4,4);		CaRs.reindexSelf(index4);
CaZs.resize(Nu-1,Nv-1,4,4);		CaZs.reindexSelf(index4);
CaKR.resize(Nu-1,Nv-1,4,4);		CaKR.reindexSelf(index4);
CaKp.resize(Nu-1,Nv-1,4,4);		CaKp.reindexSelf(index4);
CaKZ.resize(Nu-1,Nv-1,4,4);		CaKZ.reindexSelf(index4);
CanB.resize(Nu-1,Nv-1,4,4);		CanB.reindexSelf(index4);

// last closed flux surface
double s = wout.Shalf(wout.nshalf);
//cout << "VC shell is at s = " << s << endl;
du = pi2/(Nu-1);
dv = pi2/(Nv-1);
for(i=1;i<=Nu;i++) U(i) = (i-1)*du;
for(j=1;j<=Nv;j++) V(j) = (j-1)*dv;

#ifdef USE_MPI
// Prepare parallel execution
int mpi_rank = MPI::COMM_WORLD.Get_rank();
int mpi_size = MPI::COMM_WORLD.Get_size();
int mpi_parts = Nu*Nv / mpi_size;
int kstart = 1 + mpi_rank*mpi_parts;
int kend = (mpi_rank + 1)*mpi_parts;
if (mpi_rank == mpi_size - 1) kend = Nu*Nv;
int istart = int((kstart-1)/Nv) + 1;
int jstart = (kstart-1)%Nv + 1;
//cout << "Proc: " << mpi_rank << "\t" << kstart << "\t" << kend << endl;
#else
int kstart = 1;
int kend = Nu*Nv;
#endif

// Compute virtual current on LCFS
//for(i=1;i<=Nu;i++)
	//for(j=1;j<=Nv;j++)
for(k=kstart;k<=kend;k++)
	{
		i = int((k-1)/Nv) + 1;
		j = (k-1)%Nv + 1;

		wout.get_sincos(U(i), V(j), sinuv, cosuv);
		sinv = sin(V(j));
		cosv = cos(V(j));

		// Surface
		Rs(i,j) = wout.rmn.ev(s, U(i), V(j), dummy, dRdu(i,j), dRdv(i,j), d2R(i,j), sinuv, cosuv);
		Zs(i,j) = wout.zmn.ev(s, U(i), V(j), dummy, dZdu(i,j), dZdv(i,j), d2Z(i,j), sinuv, cosuv);

		// normal vector at this point
		n = normal(Rs(i,j), sinv, cosv, dRdu(i,j), dRdv(i,j), dZdu(i,j), dZdv(i,j));	// carthesian

		// Bvmec field components
		if(wout.use_nyq) wout.get_sincos(U(i), V(j), sinuv, cosuv, true);
		bsupu = wout.bsupumn.ev(s, U(i), V(j), sinuv, cosuv);
		bsupv = wout.bsupvmn.ev(s, U(i), V(j), sinuv, cosuv);

		// Bvmec
		BR = dRdu(i,j)*bsupu + dRdv(i,j)*bsupv;
		Bphi = Rs(i,j) * bsupv;
		B(2) = dZdu(i,j)*bsupu + dZdv(i,j)*bsupv;	// BZ

		// make cartesian
		B(0) = BR*cosv - Bphi*sinv;	// Bx
		B(1) = BR*sinv + Bphi*cosv;	// By

		// Bplasma = Bvmec - Vacuum B-field
		// **** the toroidal interpolation error is much larger than the numerical noise due to including the vacuum field ****
		//Bvac = get_vacuumB(Rs(i,j), V(j), Zs(i,j));
		//BR = Bvac(0); Bphi = Bvac(1);
		//Bvac(0) = BR*cosv - Bphi*sinv;	// Bx
		//Bvac(1) = BR*sinv + Bphi*cosv;	// By
		//B -= Bvac;

		// dipole moment
		nB(i,j) = n(0)*B(0) + n(1)*B(1) + n(2)*B(2);

		// virtual casing current density on s = 1; K = n x B
		KR(i,j) = n(1)*B(2) - n(2)*B(1);	// Kx
		Kp(i,j) = n(2)*B(0) - n(0)*B(2);	// Ky
		KZ(i,j) = n(0)*B(1) - n(1)*B(0);
	}

#ifdef USE_MPI
// send results to every node
double *buffer;
buffer = Rs.dataZero() + Rs.stride(firstDim)*istart + Rs.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, Rs.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = Zs.dataZero() + Zs.stride(firstDim)*istart + Zs.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, Zs.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = KR.dataZero() + KR.stride(firstDim)*istart + KR.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, KR.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = Kp.dataZero() + Kp.stride(firstDim)*istart + Kp.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, Kp.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = KZ.dataZero() + KZ.stride(firstDim)*istart + KZ.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, KZ.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = nB.dataZero() + nB.stride(firstDim)*istart + nB.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, nB.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = dRdu.dataZero() + dRdu.stride(firstDim)*istart + dRdu.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, dRdu.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = dRdv.dataZero() + dRdv.stride(firstDim)*istart + dRdv.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, dRdv.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = d2R.dataZero() + d2R.stride(firstDim)*istart + d2R.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, d2R.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = dZdu.dataZero() + dZdu.stride(firstDim)*istart + dZdu.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, dZdu.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = dZdv.dataZero() + dZdv.stride(firstDim)*istart + dZdv.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, dZdv.dataFirst(), mpi_parts, MPI::DOUBLE);
buffer = d2Z.dataZero() + d2Z.stride(firstDim)*istart + d2Z.stride(secondDim)*jstart;
MPI::COMM_WORLD.Allgather(buffer, mpi_parts, MPI::DOUBLE, d2Z.dataFirst(), mpi_parts, MPI::DOUBLE);

int origin;
if((mpi_size*mpi_parts) < (Nu*Nv))
{
	kstart = mpi_size * mpi_parts + 1;
	istart = int((kstart-1)/Nv) + 1;
	jstart = (kstart-1)%Nv + 1;
	mpi_parts = Nu*Nv - kstart + 1;
	origin = mpi_size-1;

	buffer = Rs.dataZero() + Rs.stride(firstDim)*istart + Rs.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = Zs.dataZero() + Zs.stride(firstDim)*istart + Zs.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = KR.dataZero() + KR.stride(firstDim)*istart + KR.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = Kp.dataZero() + Kp.stride(firstDim)*istart + Kp.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = KZ.dataZero() + KZ.stride(firstDim)*istart + KZ.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = nB.dataZero() + nB.stride(firstDim)*istart + nB.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = dRdu.dataZero() + dRdu.stride(firstDim)*istart + dRdu.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = dRdv.dataZero() + dRdv.stride(firstDim)*istart + dRdv.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = d2R.dataZero() + d2R.stride(firstDim)*istart + d2R.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = dZdu.dataZero() + dZdu.stride(firstDim)*istart + dZdu.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = dZdv.dataZero() + dZdv.stride(firstDim)*istart + dZdv.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
	buffer = d2Z.dataZero() + d2Z.stride(firstDim)*istart + d2Z.stride(secondDim)*jstart;
	MPI::COMM_WORLD.Bcast(buffer, mpi_parts, MPI::DOUBLE, origin);
}
#endif

// get the derivatives
bcuderiv(KR, du, dv, dKRdu, dKRdv, d2KR);
bcuderiv(Kp, du, dv, dKpdu, dKpdv, d2Kp);
bcuderiv(KZ, du, dv, dKZdu, dKZdv, d2KZ);
bcuderiv(nB, du, dv, dnBdu, dnBdv, d2nB);

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

		y_sq(1) = KR(i,j); y_sq(2) = KR(i+1,j); y_sq(3) = KR(i+1,j+1); y_sq(4) = KR(i,j+1);
		y1_sq(1) = dKRdu(i,j); y1_sq(2) = dKRdu(i+1,j); y1_sq(3) = dKRdu(i+1,j+1); y1_sq(4) = dKRdu(i,j+1);
		y2_sq(1) = dKRdv(i,j); y2_sq(2) = dKRdv(i+1,j); y2_sq(3) = dKRdv(i+1,j+1); y2_sq(4) = dKRdv(i,j+1);
		y12_sq(1) = d2KR(i,j); y12_sq(2) = d2KR(i+1,j); y12_sq(3) = d2KR(i+1,j+1); y12_sq(4) = d2KR(i,j+1);
		slice.reference(CaKR(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,du,dv,slice);

		y_sq(1) = Kp(i,j); y_sq(2) = Kp(i+1,j); y_sq(3) = Kp(i+1,j+1); y_sq(4) = Kp(i,j+1);
		y1_sq(1) = dKpdu(i,j); y1_sq(2) = dKpdu(i+1,j); y1_sq(3) = dKpdu(i+1,j+1); y1_sq(4) = dKpdu(i,j+1);
		y2_sq(1) = dKpdv(i,j); y2_sq(2) = dKpdv(i+1,j); y2_sq(3) = dKpdv(i+1,j+1); y2_sq(4) = dKpdv(i,j+1);
		y12_sq(1) = d2Kp(i,j); y12_sq(2) = d2Kp(i+1,j); y12_sq(3) = d2Kp(i+1,j+1); y12_sq(4) = d2Kp(i,j+1);
		slice.reference(CaKp(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,du,dv,slice);

		y_sq(1) = KZ(i,j); y_sq(2) = KZ(i+1,j); y_sq(3) = KZ(i+1,j+1); y_sq(4) = KZ(i,j+1);
		y1_sq(1) = dKZdu(i,j); y1_sq(2) = dKZdu(i+1,j); y1_sq(3) = dKZdu(i+1,j+1); y1_sq(4) = dKZdu(i,j+1);
		y2_sq(1) = dKZdv(i,j); y2_sq(2) = dKZdv(i+1,j); y2_sq(3) = dKZdv(i+1,j+1); y2_sq(4) = dKZdv(i,j+1);
		y12_sq(1) = d2KZ(i,j); y12_sq(2) = d2KZ(i+1,j); y12_sq(3) = d2KZ(i+1,j+1); y12_sq(4) = d2KZ(i,j+1);
		slice.reference(CaKZ(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,du,dv,slice);

		y_sq(1) = nB(i,j); y_sq(2) = nB(i+1,j); y_sq(3) = nB(i+1,j+1); y_sq(4) = nB(i,j+1);
		y1_sq(1) = dnBdu(i,j); y1_sq(2) = dnBdu(i+1,j); y1_sq(3) = dnBdu(i+1,j+1); y1_sq(4) = dnBdu(i,j+1);
		y2_sq(1) = dnBdv(i,j); y2_sq(2) = dnBdv(i+1,j); y2_sq(3) = dnBdv(i+1,j+1); y2_sq(4) = dnBdv(i,j+1);
		y12_sq(1) = d2nB(i,j); y12_sq(2) = d2nB(i+1,j); y12_sq(3) = d2nB(i+1,j+1); y12_sq(4) = d2nB(i,j+1);
		slice.reference(CanB(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,du,dv,slice);
	}
}
}


//------------------------- interpolate_RZ --------------------------------------------------------------------------------
// gets Rs(u,v), Zs(u,v) on the s = 1 surface through bicubic spline interpolation
void BFIELDVC::interpolate_RZ(double u, double v, double& Rs, double& Zs)
{
double dummy;
bcuint(U,V,CaRs,du,dv,u,v,Rs,dummy,dummy);
bcuint(U,V,CaZs,du,dv,u,v,Zs,dummy,dummy);
}


//------------------------- interpolate_VCcur -----------------------------------------------------------------------------
// gets KR(u,v), Kp(u,v) and KZ(u,v) with K = n x Bplasma on the s = 1 surface through bicubic spline interpolation
void BFIELDVC::interpolate_VCcur(double u, double v, double& KR, double& Kp, double& KZ)
{
double dummy;
bcuint(U,V,CaKR,du,dv,u,v,KR,dummy,dummy);
bcuint(U,V,CaKp,du,dv,u,v,Kp,dummy,dummy);
bcuint(U,V,CaKZ,du,dv,u,v,KZ,dummy,dummy);
}


//------------------------- interpolate_dipole ----------------------------------------------------------------------------
// gets dipole_moment(u,v) = n · Bplasma on the s = 1 surface through bicubic spline interpolation
double BFIELDVC::interpolate_dipole(double u, double v)
{
double nB,dummy;
bcuint(U,V,CanB,du,dv,u,v,nB,dummy,dummy);
return nB;
}

//------------------------ End of Class BFIELDVC --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//------------------------- biot_savart -----------------------------------------------------------------------------------
// Biot-Savart integrator for closed loop by segments in cartesian coordinates
// The numerical method of Hanson & Hirshman, Phys. Plasmas vol.9 (2002) 4410 is used
// xc(i,j), i=1,2,3 j=1...ns are (x,y,z) Cartesian position
//      vectors. There must be j=1,ns path-defining loop-point
//      position vectors.
// dv(i,j), i=1,2,3,4 j=1...ns are Cartesian direction vectors (NOT normalized!)
//      (3 elements) and segment lengths (4th element).
//      There must be j=1,ns path-defining loop-point position vectors.
//      The start point of segment j is col j of xc, its end point is
//      col j+1 of xc, except the end point of segment ns is col 1 of xc.
// I = Electric current (amp) in the closed polygon loop.
// x(i), i=1,2,3 is (x,y,z) Cartesian position vector of field point.
// b(i), i=1,2,3 is (Bx,By,Bz) Cartesian magnetic induction vector
//      (tesla) of the closed current polygon at x.
void biot_savart(Array<double,2>& xc, Array<double,2>& dv, double I, Array<double,1>& x, Array<double,1>& b)
{
int j,jend;
int ns = xc.cols();
double Rimag,Rfmag,Rsum,dbf;
Array<double,1> Ri(Range(1,3)), Rf(Range(1,3)), cross(Range(1,3));
Range all = Range::all();

b = 0;
for(j=1;j<=ns;j++)
{
    if(j < ns) jend = j+1;	// segment end point for all segments but last
    else jend = 1;			// segment end point for last segment

    Ri = x - xc(all,j);		// start point of segment j
    Rf = x - xc(all,jend);	// end point of segment j

    Rimag = sqrt(Ri(1)*Ri(1) + Ri(2)*Ri(2) + Ri(3)*Ri(3));
    Rfmag = sqrt(Rf(1)*Rf(1) + Rf(2)*Rf(2) + Rf(3)*Rf(3));
    Rsum  = Rimag + Rfmag;

    // calculate cross product between segment direction vector and Ri
    cross(1) = dv(2,j)*Ri(3) - dv(3,j)*Ri(2);
    cross(2) = dv(3,j)*Ri(1) - dv(1,j)*Ri(3);
    cross(3) = dv(1,j)*Ri(2) - dv(2,j)*Ri(1);

    // calculate Hansen & Hirshman's factor (the "missing" factor L is implicitly in dv, because dv(all,j) is not unit length, but length L)
    dbf = (2*Rsum/Rimag/Rfmag)/(Rsum*Rsum - dv(4,j)*dv(4,j));

    // the "geometric B" (the vector components of B without the common scale factors) of segment j is:
	cross *= dbf;

    // sum to get total "geometric B"
	b += cross;
}
b *= 1e-7*I;	// mu0/4pi * I
}


//------------------------- read_extcur -----------------------------------------------------------------------------------
// read coil currents from separate file, if VMEC is fixed boundary mode
// checks VMEC input file first, if not then checks file "mgrid_coil_currents.in"
// expects extcur input syntax as in a free boundary VMEC input file
int read_extcur(Array<double,1>& extcur, LA_STRING input_extension)
{
#ifdef USE_MPI
int mpi_rank = MPI::COMM_WORLD.Get_rank();
#else
int mpi_rank = 0;
#endif

string line,word,equal,index;
double value;
int pos1,i,nextcur = 0;
TinyVector <int,1> index1(1);
extcur.resize(100); extcur.reindexSelf(index1);

LA_STRING inputfile = LA_STRING("input.") + input_extension;
ifstream in;
in.open(inputfile);
if(in.fail() == 1)
{
	in.clear();
	in.open("mgrid_coil_currents.in");
	if(in.fail() == 1)  {if(mpi_rank == 0) cout << "XPAND: No coil currents found. Run XPAND with -V option." << endl; EXIT;}
	else {if(mpi_rank == 0) cout << "XPAND: EXTCUR found in: mgrid_coil_currents.in" << endl;}
}
else {if(mpi_rank == 0) cout << "XPAND: EXTCUR found in: " << inputfile << endl;}
while(getline(in, line))
{
	if(line.length() < 1) continue; 	// blank lines anywhere don't matter
	if(line.find("EXTCUR") != std::string::npos)
	{
		stringstream ss(line);
		while(ss >> word >> equal >> value)
		{
			pos1 = word.find("(");
			index = word.substr(pos1 + 1, 2);
			i = atoi(index.c_str());	// in c++ 11 (needs -std=c++11 flag for compiler) you can do:  stoi(index);
			extcur(i) = value;
			nextcur += 1;
		}
	}
}
extcur.resizeAndPreserve(nextcur);
return nextcur;
}


#endif //  XPAND_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

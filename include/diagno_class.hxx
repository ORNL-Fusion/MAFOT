// Class provides B-field from VMEC (inside s = 1) and Diagno (outside s = 1)
// A.Wingen						10.10.14


// Define
//--------
#ifndef DIAGNO_CLASS_INCLUDED
#define DIAGNO_CLASS_INCLUDED

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


//--------- Begin Class DIAGNO -------------------------------------------------------------------------------------------
class DIAGNO
{
private:
	// Member Variables
	double dR, dZ, dp;
	Array<double,1> Ra,Za;
	Array<double,3> dBRdR, dBPHIdR, dBZdR;
	Array<double,3> dBRdZ, dBPHIdZ, dBZdZ;
	Array<double,3> d2BR, d2BPHI, d2BZ;
	Array<double,5> CaBR, CaBPHI, CaBZ;

	// Member-Functions
	void prep_interpolation(void);			// initiate the interpolations
	void interpolate_B(double R, double Z, int k, double& br, double& bphi, double& bz);	// interpolate all 3 B-field components on the (R,phi[k],Z) grid

public:
	// Member Variables
	int NR, NZ, Np;
	double Rmin, Rmax, Zmin, Zmax;

	Array<double,3> BR;
	Array<double,3> BPHI;
	Array<double,3> BZ;

	// Constructors
	DIAGNO();								// Default Constructor

	// Member-Operators
	DIAGNO& operator =(const DIAGNO& spec);	// Operator =

	// Member-Functions
	void read(LA_STRING filename);			// read file with B-field on 3-D grid (R,phi,Z)
	void get_B(double R, double phi, double Z, double& br, double& bphi, double& bz);		// evaluate B at any location (R,phi,Z)

}; //end of class

//------------------------ End of Class DIAGNO ----------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
DIAGNO::DIAGNO()
{
TinyVector <int,1> index(1);
TinyVector <int,3> index3(1,1,1);
TinyVector <int,5> index5(1,1,1,1,1);

NR = 129;
NZ = 129;
Np = 48;

Rmin = 1; 		Rmax = 2.4;
Zmin = -1.5; 	Zmax = 1.5;

dR = (Rmax - Rmin) / (NR-1);
dZ = (Zmax - Zmin) / (NZ-1);
dp = pi2/(Np-1);				// 0 -> 2*pi

Ra.resize(NR);				Ra.reindexSelf(index);
Za.resize(NZ);				Za.reindexSelf(index);

BR.resize(NR, NZ, Np);		BR.reindexSelf(index3);
BPHI.resize(NR, NZ, Np);	BPHI.reindexSelf(index3);
BZ.resize(NR, NZ, Np);		BZ.reindexSelf(index3);

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
DIAGNO& DIAGNO::operator =(const DIAGNO& dia)
{
if (this == &dia) return(*this);	    // if: x=x

NR = dia.NR;
NZ = dia.NZ;
Np = dia.Np;

Rmin = dia.Rmin; 	Rmax = dia.Rmax;
Zmin = dia.Zmin; 	Zmax = dia.Zmax;

dR = dia.dR;
dZ = dia.dZ;
dp = dia.dp;

Ra.reference(dia.Ra.copy());
Za.reference(dia.Za.copy());

BR.reference(dia.BR.copy());
BPHI.reference(dia.BPHI.copy());
BZ.reference(dia.BZ.copy());

dBRdR.reference(dia.dBRdR.copy());
dBPHIdR.reference(dia.dBPHIdR.copy());
dBZdR.reference(dia.dBZdR.copy());

dBRdZ.reference(dia.dBRdZ.copy());
dBPHIdZ.reference(dia.dBPHIdZ.copy());
dBZdZ.reference(dia.dBZdZ.copy());

d2BR.reference(dia.d2BR.copy());
d2BPHI.reference(dia.d2BPHI.copy());
d2BZ.reference(dia.d2BZ.copy());

CaBR.reference(dia.CaBR.copy());
CaBPHI.reference(dia.CaBPHI.copy());
CaBZ.reference(dia.CaBZ.copy());

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------------------- read ---------------------------------------------------------------------------------------
// read B-field file and set member variables
void DIAGNO::read(LA_STRING filename)
{
// Variables
int i,j,k;
string stdummy;
double dummy;
Range all = Range::all();
LA_STRING line;

// Input
ifstream in;
in.open(filename);
if(in.fail()==1) {cout << "Unable to open file " << filename << endl; EXIT;}

// skip first line
in >> line;

// Read the dimensions, 3 integers
in >> stdummy;
in >> NR;
in >> Np;
in >> NZ;
Np += 1;	// increase Np by 1 to make periodic

// skip third line
in >> line;

// resize data: (1 -> Np, 1 -> NZ, 1 -> NR)
Array<double,3> R(Range(1,NR), Range(1,NZ), Range(1,Np));
Array<double,3> Z(Range(1,NR), Range(1,NZ), Range(1,Np));

BR.resize(NR, NZ, Np);
BPHI.resize(NR, NZ, Np);
BZ.resize(NR, NZ, Np);

// Read data
for(i=2;i<=Np;i++)	// start at 2, since Np is one larger; first point = last point with first point missing
{
	for(j=1;j<=NZ;j++)
	{
		for(k=1;k<=NR;k++)
		{
			in >> R(k,j,i);
			in >> dummy;	// phi, but not necessary
			in >> Z(k,j,i);
			in >> BR(k,j,i);
			in >> BPHI(k,j,i);
			in >> BZ(k,j,i);
		}
	}
}
in.close();

// make periodic
R(all,all,1) = R(all,all,Np);
Z(all,all,1) = Z(all,all,Np);
BR(all,all,1) = BR(all,all,Np);
BPHI(all,all,1) = BPHI(all,all,Np);
BZ(all,all,1) = BZ(all,all,Np);

// set grid parameters
Rmin = min(R);
Rmax = max(R);
Zmin = min(Z);
Zmax = max(Z);

dR = (Rmax - Rmin) / (NR-1);
dZ = (Zmax - Zmin) / (NZ-1);
dp = pi2/(Np-1);				// 0 -> 2*pi

Ra.resize(NR);
for(i=1;i<=NR;i++) Ra(i) = Rmin + (i-1)*dR;
Za.resize(NZ);
for(i=1;i<=NZ;i++) Za(i) = Zmin + (i-1)*dZ;

// prepare the bicubic splines
prep_interpolation();
}


//-------------------------------- get_B ----------------------------------------------------------------------------------
// evaluate br, bphi, bz at any arbitrary location (R, phi, Z)
void DIAGNO::get_B(double R, double phi, double Z, double& br, double& bphi, double& bz)
{
int k;
double t;
double bru, bphiu, bzu;

// locate neighboring phi planes: lower (or matching) phi = (k-1)*dp; upper phi = k*dp
phi = modulo2pi(phi);	// make sure phi in [0, 2pi]
k = int(phi/dp) + 1;
t = phi/dp - k + 1;

interpolate_B(R, Z, k, br, bphi, bz);
//cout << k << "\t" << t << "\t" << R << "\t" << Z << "\t" << br << "\t" << bphi << "\t" << bz << endl;
if(t > 0)	// not exactly in the phi-plane (k-1)*dp
{
	interpolate_B(R, Z, k+1, bru, bphiu, bzu);

	// linear interpolation of B-field in phi
	br += t*(bru - br);
	bphi += t*(bphiu - bphi);
	bz += t*(bzu - bz);
}
}


//----------------------- Private Member Functions ------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------------------- prep_interpolation --------------------------------------------------------------------------------
// get the C coefficients for the interpolation
void DIAGNO::prep_interpolation(void)
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
void DIAGNO::interpolate_B(double R, double Z, int k, double& br, double& bphi, double& bz)
{
double dummy;
Array<double,4> slice;
Range all = Range::all();

slice.reference(CaBR(all,all,k,all,all));
bcuint(Ra,Za,slice,dR,dZ,R,Z,br,dummy,dummy);

slice.reference(CaBPHI(all,all,k,all,all));
bcuint(Ra,Za,slice,dR,dZ,R,Z,bphi,dummy,dummy);

slice.reference(CaBZ(all,all,k,all,all));
bcuint(Ra,Za,slice,dR,dZ,R,Z,bz,dummy,dummy);
}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



#endif //  DIAGNO_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

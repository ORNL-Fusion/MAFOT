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
using namespace blitz;

// Prototypes
//-----------
double* adaptiveSimpson(double (*f)(double), double a, double b, double epsilon, int maxRecursionDepth);
double* adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon, double S, double fa, double fb, double fc, int bottom);
void spline(Array<double,1>& x, Array<double,1>& y, int n, double yp1, double ypn, Array<double,1>& y2);
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, int n, double x, double& y, double& yx);

// Golbal Parameters
//------------------
extern ofstream ofs2;


//--------- Begin Class POTENTIAL -------------------------------------------------------------------------------------
class POTENTIAL
{
private:
	// Member Variables
	double eps; 				// accuracy tolerance in adaptive Simpson, default = 1e-8
	double maxRecDepth; 		// max recursion in adaptive Simpson, defaut = 10

public:
	// Member Variables
	VMEC wout;
	MGRID mgrid;

	double R;		// Point where to evaluate
	double phi;
	double Z;

	// Constructors
	POTENTIAL();								// Default Constructor

	// Member-Operators
	POTENTIAL& operator =(const POTENTIAL& spec);	// Operator =

	// Member-Functions
	double ev(double R, double phi, double Z);
	double integ_v(double u);
	double* integs(double u, double v);
	Array<double,1> normal(double dRdu, double dRdv, double dZdu, double dZdv);
	double green(double Rs, double v, double Zs);
	Array<double,1> get_vacuumB(double Rs, double v, double Zs);
	Array<double,1> get_plasmaCurrentB(double Rs, double v, double Zs);
	void read_VMEC_input(VMEC wout);


}; //end of class

//------------------------ End of Class VMEC_SPECTRAL----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
POTENTIAL::POTENTIAL()
{

}

//--------- Operator = ----------------------------------------------------------------------------------------------------
POTENTIAL& POTENTIAL::operator =(const POTENTIAL& spec)
{
if (this == &spec) return(*this);	    // if: x=x

ymnc.reference(spec.ymnc.copy());

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- ev -----------------------------------------------------------------------------------
double POTENTIAL::ev(double Rin, double phiin, double Zin)
{
double integ[2];
// load into member variables
R = Rin; phi = phiin; Z = Zin;
integ = adaptiveSimpson(integ_v, 0, pi2, eps, maxRecDepth);
return -(integ[0] + integ[1]) / pi2
}

// --- integ_v -----------------------------------------------------------------------------------
double POTENTIAL::integ_v(double u)
{
auto f = [] (v) {integs(u,v)}
return adaptiveSimpson(f, 0, pi2, eps, maxRecDepth);
}

// --- integs -----------------------------------------------------------------------------------
double* POTENTIAL::integs(double u, double v)
{
double Rs, Zs, G, pot;
double dummy, dZdu, dZdv, dRdu, dRdv;
double out[2];
Array<double,1> n, B;

// point on s = 1 surface
Rs = wout.rmn.ev(1.0, u, v, dummy, dRdu, dRdv);
Zs = wout.zmn.ev(1.0, u, v, dummy, dZdu, dZdv);

// normal vector at this point
n = normal(dRdu, dRdv, dZdu, dZdv);

// evaluate Green's function
G = green(Rs, v, Zs);

// scalar potential on s = 1 surface
pot = wout.pot(u, v);

// vacuum magbetic field on s = 1 surface
B = get_vacuumB(Rs, v, Zs);


I1 = (n(0) * (R - Rs) + n(1) * (phi - v) + n(2) * (Z - Zs)) * G*G*G * pot;
I2 = (n(0)*B(0) + n(1)*B(1) + n(2)*B(2)) * G;

out = {I1, I2};
return out;
}

// --- normal -----------------------------------------------------------------------------------
Array<double,1> POTENTIAL::normal(double drdu, double drdv, double dzdu, double dzdv)
{
Array<double,1> n(3);
n(0) = dzdu;
n(1) = dzdu*drdv - drdu*dzdv;
n(2) = drdu;

n /= sqrt(n(0)*n(0) + n(1)*n(1) + n(2)*n(2));
return n;
}

// --- green -----------------------------------------------------------------------------------
double POTENTIAL::green(double Rs, double v, double Zs)
{
return 1.0/sqrt((R - Rs)*(R - Rs) + (Z - Zs)*(Z - Zs) + (phi - v)*(phi - v));
}


// --- get_vacuumB -----------------------------------------------------------------------------------
Array<double,1> POTENTIAL::get_vacuumB(double Rs, double v, double Zs)
{
double dphi = pi2 / kp;
int k = int(round(v/dphi)) % kp;
Array<double,1> B(3);
double br, bphi,bz;

mgrid.interpolate_B(Rs, Zs, k, br, bphi, bz);
B(0) = br;
B(1) = bphi;
B(2) = bz;

B += get_plasmaCurrentB(Rs, v, Zs);

return B;
}








//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------










//--------- Begin Class MGRID -------------------------------------------------------------------------------------
class MGRID
{
private:
	// Member Variables
	LA_STRING mgrid_file;		// pathname of MGRID file
	Array<double,1> extcur;		// coil currents

public:
	// Member Variables
	Array<double,1> R
	Array<double,1> Z
	int kp;						// number of toroidal planes in mgrid
	Array<double,3> BR;
	Array<double,3> BPHI;
	Array<double,3> BZ;

	// Constructors
	MGRID();								// Default Constructor

	// Member-Operators
	MGRID& operator =(const MGRID& spec);	// Operator =

	// Member-Functions
	void read(void);
	void prep_interpolation(void);	// initiate the interpolations
	void interpolate_B(double R, double Z, int k, double& br, double& bphi, double& bz);	// interpolate all 3 vacuum B-field components on the mgrid grid


}; //end of class

//------------------------ End of Class VMEC_SPECTRAL----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
MGRID::MGRID()
{

}

//--------- Operator = ----------------------------------------------------------------------------------------------------
MGRID& MGRID::operator =(const MGRID& spec)
{
if (this == &spec) return(*this);	    // if: x=x

ymnc.reference(spec.ymnc.copy());

return(*this);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


// --- prep_vacuumB -----------------------------------------------------------------------------------
void MGRID::read(void)
{
// Variables
int chk, ncid, varid;
firstIndex k;

// open file
chk = nc_open(mgrid_file, NC_NOWRITE, &ncid);
if(chk != 0) {cout << "Unable to open file " << mgrid_file << endl; EXIT;}

// read grid data
int ir, jz;
chk = nc_inq_varid(ncid, "ir", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &ir);	// read variable
chk = nc_inq_varid(ncid, "jz", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &jz);	// read variable
chk = nc_inq_varid(ncid, "kp", &varid);	// get variable id
chk = nc_get_var_int(ncid, varid, &kp);	// read variable

double rmin, rmax, zmin, zmax;
chk = nc_inq_varid(ncid, "rmin", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &rmin);		// read
chk = nc_inq_varid(ncid, "rmax", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &rmax);		// read
chk = nc_inq_varid(ncid, "zmin", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &zmin);		// read
chk = nc_inq_varid(ncid, "zmax", &varid);			// get variable id
chk = nc_get_var_double(ncid, varid, &zmax);		// read

R.resize(ir);
Z.resize(jz);
R = rmin + k*(rmax - rmin)/double(ir);
Z = zmin + k*(zmax - zmin)/double(jz);

// read B-field data
Array<double,2> input;
input.resize(kp, jz, ir);
BR.resize(kp, jz, ir); 		BR = 0;
BPHI.resize(kp, jz, ir);	BPHI = 0;
BZ.resize(kp, jz, ir);		BZ = 0;

for (i=1;i<=36;i++)
{
	chk = nc_inq_varid(ncid, "br_"+format(i,'03d'), &varid);			// get variable id
	chk = nc_get_var_double(ncid, varid, input.data());	// read
	BR += input.copy();								// move into place

	// BPHI & BZ too
}
}

//--------------------- prep_interpolation --------------------------------------------------------------------------------
// get the 2. derivative in s-direction, which is input to splint_2D, for all three fields, all u's in all v-planes
void POTENTIAL::prep_interpolation(void)
{
int i,j,k;
Array<double,2> slice, sliceds, slicedu, sliced2;
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
Range all = Range::all();

dBRds.resize(Ns,Nu,Nv);
dBPHIds.resize(Ns,Nu,Nv);
dBZds.resize(Ns,Nu,Nv);

dBRdu.resize(Ns,Nu,Nv);
dBPHIdu.resize(Ns,Nu,Nv);
dBZdu.resize(Ns,Nu,Nv);

d2BR.resize(Ns,Nu,Nv);
d2BPHI.resize(Ns,Nu,Nv);
d2BZ.resize(Ns,Nu,Nv);

CaBR.resize(Ns-1,Nu-1,Nv-1,4,4);
CaBPHI.resize(Ns-1,Nu-1,Nv-1,4,4);
CaBZ.resize(Ns-1,Nu-1,Nv-1,4,4);

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
}


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
		}
	}
}
}





//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------













double* adaptiveSimpsonsAux(double (*f)(double), double a, double b, double epsilon,
                         double S, double fa, double fb, double fc, int bottom)
{
double c = (a + b)/2, h = b - a;
double d = (a + c)/2, e = (c + b)/2;
double fd = f(d), fe = f(e);
double Sleft = (h/12)*(fa + 4*fd + fc);
double Sright = (h/12)*(fc + 4*fe + fb);
double S2 = Sleft + Sright;
if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
	return S2 + (S2 - S)/15;
return adaptiveSimpsonsAux(f, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
	 adaptiveSimpsonsAux(f, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

//
// Adaptive Simpson's Rule
//
double* adaptiveSimpson(double (*f)(double),   // ptr to function
                           double a, double b,  // interval [a,b]
                           double epsilon,  // error tolerance
                           int maxRecursionDepth) {   // recursion cap
double c = (a + b)/2, h = b - a;
double fa = f(a), fb = f(b), fc = f(c);
double S = (h/6)*(fa + 4*fc + fb);
return adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}

#endif //  EXTENDER_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

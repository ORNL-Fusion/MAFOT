// EFIT Class handles all EFIT data read from g-files
// provides interpolation member functions for all data sets in g-file
// used by all ITER drift programs
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						5.08.10


// Define
//--------
#ifndef EFIT_CLASS_INCLUDED
#define EFIT_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;

// Prototypes  
void spline(Array<double,1>& x, Array<double,1>& y, int n, double yp1, double ypn, Array<double,1>& y2);
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, int n, double x, double& y, double& yx);
void splint_2D(Array<double,1>& x1a, Array<double,1>& x2a, Array<double,2>& ya, Array<double,2>& d2ydx1, int n1, int n2,
				double x1, double x2, double& y, double& dydx1, double& dydx2);

void bcuderiv(Array<double,2>& y, double d1, double d2, Array<double,2>& y1, Array<double,2>& y2, Array<double,2>& y12);
void bcucof(Array<double,1>& y, Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12, double d1, double d2, 
			Array<double,2>& c);
int bcuint(Array<double,1>& Ra, Array<double,1>& Za, Array<double,4>& Ca, double dR, double dZ,
			double R, double Z, double& y, double& y1, double& y2);

// Golbal Parameters 
#if defined(windows)		// Windows system
	#if defined(ITER)
	const LA_STRING gFilePath = "C:\\C++\\ITER\\gfiles\\";	// double backslash causes \ to appear in string
	#elif defined(NSTX)
	const LA_STRING gFilePath = "C:\\C++\\NSTX\\gfiles\\";
	#else
	const LA_STRING gFilePath = "C:\\C++\\D3D\\d3gfiles\\";
	#endif
#elif defined(genatomix)	// GA u-drive
	#if defined(ITER)
	const LA_STRING gFilePath = "/u/wingen/d3gfiles/iter/";
	#elif defined(NSTX)
	const LA_STRING gFilePath = "/u/wingen/d3gfiles/nstx/";
	#else
	const LA_STRING gFilePath = "/u/wingen/d3gfiles/";
	#endif
#elif defined(mac)			// Mac Book
	#if defined(ITER)
	const LA_STRING gFilePath = "/Users/wingen/c++/iter/gfiles/";
	#elif defined(NSTX)
	const LA_STRING gFilePath = "/Users/wingen/c++/nstx/gfiles/";
	#else
	const LA_STRING gFilePath = "/Users/wingen/c++/d3d/gfiles/";
	#endif
#else						// default linux system
	#if defined(ITER)
	const LA_STRING gFilePath = "/home/wingen/c++/iter/gfiles/";
	#elif defined(NSTX)
	const LA_STRING gFilePath = "/home/wingen/c++/nstx/gfiles/";
	#else
	const LA_STRING gFilePath = "/home/wingen/c++/d3d/gfiles/";
	#endif
#endif

//--------- Begin Class EFIT ----------------------------------------------------------------------------------------------
class EFIT
{
private:
// Member Variables
	int helicity_adjust;	// default = 1; set to -1 if helicity is wrong -> modifies dpsidR und dpsidZ

public: 
// Member Variables
	LA_STRING Shot;		// Shot number
	LA_STRING Time;		// Time slice in ms
	LA_STRING Path;		// Path to g-file

	// To be read from g-file
	int NR;	// Number of Points in R-direction
	int NZ;	// Number of Points in Z-direction

	double Xdim;	// Distance from inner edge to outer edge covered by EFIT, in [m]
	double Zdim;	// Distance from bottom to top (Z_axis) covered by EFIT, equally spaced around midplane, in [m]
	double R0;		// Major Radius of torus in [m]
	double R1;		// Position of inner edge on radial scale in [m]
	double Zmid;	// Position of midplane on z-scale in [m]
	double RmAxis;	// R-position of magnetic Axis in [m]
	double ZmAxis;	// Z-position of magnetic Axis in [m]
	double psiAxis;	// poloidal Flux at magnetic Axis
	double psiSep;	// poloidal Flux at Separatrix
	double Bt0;		// toroidal magnetic field at Major Radius in [T]
	double Ip;		// Plasma current in [A]

	Array<double,1> Fpol;		// Fpol = Fpol(psi_norm) = Btor * R 
	Array<double,1> Pres;		// Pres = Plasma pressure profile
	Array<double,1> FFprime;	// FFprime = Fpol * dFpol/dpsi
	Array<double,1> Pprime;		// Pprime = Pressure derivative dp/dpsi
	Array<double,1> qpsi;		// q(psi) = Safety Factor

	Array<double,2> psiRZ;	// psi(R,Z) = NOT-normalized poloidal Flux		i=1:Nr rows, j=1:Nz columns  

	int Nlcfs;	// number of points for the Last Closed Flux Surface (LCFS), Array has twice the size!
	int Nwall;	// number of points for the Wall, Array has twice the size!	

	Array<double,1> lcfs;	// Position of the LCFS in R,Z plane; R,Z coordinates are alternating: R1,Z1,R2,Z2,R3,Z3,...
	Array<double,1> wall;	// Position of the Wall in R,Z plane; R,Z coordinates are alternating: R1,Z1,R2,Z2,R3,Z3,...
	Array<double,1> lcfs_th;// poloidal angles of lcfs

	// Calculated in ReadData
	double dR;		// grid distance in R direction
	double dZ;		// grid distance in Z direction
	double dpsi;	// grid distance in Z direction

	Array<double,1> R;		// R grid
	Array<double,1> Z;		// Z grid
	Array<double,1> psi;	// Psi_norm grid

	Array<double,1> d2Fpol;		// Spline of Fpol in psi_norm
	Array<double,1> d2Pres;		// Spline of Pres in psi_norm
	Array<double,1> d2FFprime;	// Spline of FFprime in psi_norm
	Array<double,1> d2Pprime;	// Spline of Pprime in psi_norm
	Array<double,1> d2qpsi;		// Spline of qpsi in psi_norm

	Array<double,2> d2psi;		// Spline of psi in R-direction
	Array<double,2> dpsidR;		// Derivative of psi in R-direction
	Array<double,2> dpsidZ;		// Derivative of psi in Z-direction
	Array<double,2> d2psidRdZ;	// Cross-derivative of psi

	Array<double,4> Ca;		// Coefficients for bcuint

	int helicity;	// +1 = right handed equlibrium		-1 = left handed

// Konstruktoren
	EFIT();		// Default Constructor

// Member-Operatoren
	EFIT& operator =(const EFIT& EQD);		// Operator =

// Member-Funktionen
	void ReadData(const LA_STRING ShotNr, const LA_STRING ShotTime, double Raxis = 0, double Zaxis = 0);	// Reads EFIT data from g-file

	// Interpolates psi and derivatives; flag: 0 = 2D splines, 1 = bicubic interpolation; norm specifies if psi is to be normalized
	int get_psi(const double x1, const double x2, double& y, double& dy1, double& dy2, int flag=1, bool norm=true);		// 0: ok	-1: Point outside of EFIT grid
	
	double get_Fpol(const double x);		// Spline interpolates Fpol
	double get_Pres(const double x);		// Spline interpolates Pres
	double get_FFprime(const double x);		// Spline interpolates FFprime
	double get_Pprime(const double x);		// Spline interpolates Pprime
	double get_q(const double x);		// Spline interpolates qpsi
	int lcfs_RZ_nn(const double th, double& r, double& z);

}; //end of class

//------------------------ End of Class -----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Default Constructor -------------------------------------------------------------------------------------------
EFIT::EFIT()
{
TinyVector <int,1> index(1);
TinyVector <int,2> index2(1,1);
TinyVector <int,4> index4(1,1,1,1);

Shot = "none";
Time = "00000";
Path = gFilePath;

NR = 129;
NZ = 129;

Fpol.resize(NR);	Fpol.reindexSelf(index);
Pres.resize(NR);	Pres.reindexSelf(index);
FFprime.resize(NR);	FFprime.reindexSelf(index);
Pprime.resize(NR);	Pprime.reindexSelf(index);
qpsi.resize(NR);	qpsi.reindexSelf(index);

psiRZ.resize(NR,NZ);	psiRZ.reindexSelf(index2);

Nlcfs = 90;
Nwall = 87;
lcfs.resize(2*Nlcfs);	lcfs.reindexSelf(index);
wall.resize(2*Nwall);	wall.reindexSelf(index);
lcfs_th.resize(Nlcfs);	lcfs_th.reindexSelf(index);

R.resize(NR);	R.reindexSelf(index);
Z.resize(NZ);	Z.reindexSelf(index);
psi.resize(NR);	psi.reindexSelf(index);

d2Fpol.resize(NR);		d2Fpol.reindexSelf(index);
d2Pres.resize(NR);		d2Pres.reindexSelf(index);
d2FFprime.resize(NR);	d2FFprime.reindexSelf(index);
d2Pprime.resize(NR);	d2Pprime.reindexSelf(index);
d2qpsi.resize(NR);		d2qpsi.reindexSelf(index);

d2psi.resize(NR,NZ);		d2psi.reindexSelf(index2);
dpsidR.resize(NR,NZ);		dpsidR.reindexSelf(index2);
dpsidZ.resize(NR,NZ);		dpsidZ.reindexSelf(index2);
d2psidRdZ.resize(NR,NZ);	d2psidRdZ.reindexSelf(index2);

Ca.resize(NR-1,NZ-1,4,4);	Ca.reindexSelf(index4);

helicity = 1;
helicity_adjust = 1;
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
EFIT& EFIT::operator =(const EFIT& EQD)
{
if (this == &EQD) return(*this);	    // if: x=x

Shot = EQD.Shot;
Time = EQD.Time;
Path = EQD.Path;

NR = EQD.NR;
NZ = EQD.NZ;

Xdim = EQD.Xdim;
Zdim = EQD.Zdim;
R0 = EQD.R0;
R1 = EQD.R1;
Zmid = EQD.Zmid;
RmAxis = EQD.RmAxis;
ZmAxis = EQD.ZmAxis;
psiAxis = EQD.psiAxis;
psiSep = EQD.psiSep;
Bt0 = EQD.Bt0;
Ip = EQD.Ip;

Fpol.reference(EQD.Fpol.copy());
Pres.reference(EQD.Pres.copy());
FFprime.reference(EQD.FFprime.copy());
Pprime.reference(EQD.Pprime.copy());
qpsi.reference(EQD.qpsi.copy());

psiRZ.reference(EQD.psiRZ.copy());

Nlcfs = EQD.Nlcfs;
Nwall = EQD.Nwall;
lcfs.reference(EQD.lcfs.copy());
wall.reference(EQD.wall.copy());
lcfs_th.reference(EQD.lcfs_th.copy());

dR = EQD.dR;
dZ = EQD.dZ;
dpsi = EQD.dpsi;

R.reference(EQD.R.copy());
Z.reference(EQD.Z.copy());
psi.reference(EQD.psi.copy());

d2Fpol.reference(EQD.d2Fpol.copy());
d2Pres.reference(EQD.d2Pres.copy());
d2FFprime.reference(EQD.d2FFprime.copy());
d2Pprime.reference(EQD.d2Pprime.copy());
d2qpsi.reference(EQD.d2qpsi.copy());

d2psi.reference(EQD.d2psi.copy());
dpsidR.reference(EQD.dpsidR.copy());
dpsidZ.reference(EQD.dpsidZ.copy());
d2psidRdZ.reference(EQD.d2psidRdZ.copy());

Ca.reference(EQD.Ca.copy());

helicity = EQD.helicity;
helicity_adjust = EQD.helicity_adjust;

return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------------- ReadEFITdata -------------------------------------------------------------------------------------
void EFIT::ReadData(LA_STRING ShotNr, LA_STRING ShotTime, double Raxis, double Zaxis)
{
int i,j;
string stdummy;
double dummy;

int NZ2;
firstIndex k;

double d1,dn;
Range all = Range::all();
Array<double,1> slice_1,slice_2;	//Slice-array

// Set Shot and Time
Shot = ShotNr;
Time = ShotTime;

// Open File
for(i=1;i<=5-int(ShotTime.len());i++) ShotTime.insert('0');		//insert zeros until ShotTime has a length of 5
LA_STRING filename = "g" + ShotNr + "." + ShotTime;
ifstream in;
in.open(Path + filename);
if(in.fail()==1) {cout << "Unable to open file " << Path + filename << endl; exit(0);}

// Read tags and dimension
in >> stdummy;
while (stdummy[0] != '#') in >> stdummy;
in >> stdummy;
in >> dummy;	
in >> NR;	// Number of Points in R-direction
in >> NZ;	// Number of Points in Z-direction
//in >> stdummy;	// useless string: endHead 

// Rezize Arrays (if size = 129 nothing is done!); All Arrays start with index 1 
Fpol.resize(NR);
Pres.resize(NR);
FFprime.resize(NR);
Pprime.resize(NR);
psiRZ.resize(NR,NZ);
qpsi.resize(NR);

// Read Parameters
in >> Xdim;		// Distance from inner edge to outer edge covered by EFIT, in [m]
in >> Zdim;		// Distance from bottom to top (Z_axis) covered by EFIT, equally spaced around midplane, in [m]
in >> R0;		// Major Radius of torus in [m]
in >> R1;		// Position of inner edge on radial scale in [m]
in >> Zmid;		// Position of midplane on z-scale in [m]
in >> RmAxis;	// R-position of magnetic Axis in [m]
in >> ZmAxis;	// Z-position of magnetic Axis in [m]
in >> psiAxis;	// poloidal Flux at magnetic Axis
in >> psiSep;	// poloidal Flux at Separatrix
in >> Bt0;		// toroidal magnetic field at Major Radius in [T]
in >> Ip;		// Plasma current in [A]
for (i=1;i<=9;i++) in >> dummy;		// Unused parameters

// Read Arrays
for(i=1;i<=NR;i++) in >> Fpol(i);	// Fpol = Fpol(psi_norm) = Btor * R; Fpol(1) = Fpol(psi=0), Fpol(NR) = Fpol(psi=1), same for the others
for(i=1;i<=NR;i++) in >> Pres(i);	// Pres = pressure profile
for(i=1;i<=NR;i++) in >> FFprime(i);	// FFprime = ?
for(i=1;i<=NR;i++) in >> Pprime(i);	// Pprime = dPres/dpsi
for(j=1;j<=NZ;j++) 
	for (i=1;i<=NR;i++) in >> psiRZ(i,j);	// psi(R,Z) = NOT-normalized poloidal Flux		i=1:Nr rows, j=1:Nz columns  
for(i=1;i<=NR;i++) in >> qpsi(i);	// q(psi) = Safety Factor

// Read additional stuff
in >> Nlcfs;	// twice the number of points for the Last Closed Flux Surface (LCFS)
in >> Nwall;	// twice the number of points for the Wall
lcfs.resize(2*Nlcfs);
wall.resize(2*Nwall);
for(i=1;i<=2*Nlcfs;i++) in >> lcfs(i);	// Position of the LCFS in R,Z plane; R,Z coordinates are alternating: R1,Z1,R2,Z2,R3,Z3,...
for(i=1;i<=2*Nwall;i++) in >> wall(i);	// Position of the Wall in R,Z plane; R,Z coordinates are alternating: R1,Z1,R2,Z2,R3,Z3,...

in.close();

// shift R,Z grid
if(Raxis > 0)
{
	double dRmAxis = Raxis - RmAxis;
	double dZmAxis = Zaxis - ZmAxis;

	RmAxis += dRmAxis;
	ZmAxis += dZmAxis;
	R1 += dRmAxis;
	Zmid += dZmAxis;
	for(i=1;i<=2*Nlcfs;i+=2)
	{
		lcfs(i) += dRmAxis;
		lcfs(i+1) += dZmAxis;
	}
}

// Set R, Z and psi Arrays
R.resize(NR);
Z.resize(NZ);
psi.resize(NR);

dR = Xdim/double(NR-1);
R = R1 + (k-1)*dR;

dpsi = 1.0/double(NR-1);
psi = (k-1)*dpsi;

dZ = Zdim/double(NZ-1);
NZ2 = int(floor(0.5*NZ));
Z = Zmid + (k-NZ2-1)*dZ;	// Z(NZ2+1) = Zmid

// Prepare Bicubic Spline interpolation of psiRZ
d2psi.resize(NR,NZ);

for(i=1;i<=NZ;i++) 
{
	d1 = (psiRZ(2,i)-psiRZ(1,i))/dR;
	dn = (psiRZ(NR,i)-psiRZ(NR-1,i))/dR;
	slice_1.reference(psiRZ(all,i));
	slice_2.reference(d2psi(all,i));
	spline(R,slice_1,NR,d1,dn,slice_2);
}

// Prepare Bicubic Interpolation of psiRZ -> get gradients and cross-derivative
dpsidR.resize(NR,NZ);
dpsidZ.resize(NR,NZ);
d2psidRdZ.resize(NR,NZ);

bcuderiv(psiRZ,dR,dZ,dpsidR,dpsidZ,d2psidRdZ);

// Get the c's for bcuint, as done by bcucof
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
Array<double,2> slice;
Ca.resize(NR-1,NZ-1,4,4);
for(i=1;i<NR;i++)
{
	for(j=1;j<NZ;j++)
	{
		y_sq(1) = psiRZ(i,j); y_sq(2) = psiRZ(i+1,j); y_sq(3) = psiRZ(i+1,j+1); y_sq(4) = psiRZ(i,j+1);
		y1_sq(1) = dpsidR(i,j); y1_sq(2) = dpsidR(i+1,j); y1_sq(3) = dpsidR(i+1,j+1); y1_sq(4) = dpsidR(i,j+1);
		y2_sq(1) = dpsidZ(i,j); y2_sq(2) = dpsidZ(i+1,j); y2_sq(3) = dpsidZ(i+1,j+1); y2_sq(4) = dpsidZ(i,j+1);
		y12_sq(1) = d2psidRdZ(i,j); y12_sq(2) = d2psidRdZ(i+1,j); y12_sq(3) = d2psidRdZ(i+1,j+1); y12_sq(4) = d2psidRdZ(i,j+1);

		slice.reference(Ca(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,slice); 
	}
}

// Prepare Spline interpolation of Fpol
d2Fpol.resize(NR);
d1 = (Fpol(2)-Fpol(1))/dpsi;
dn = (Fpol(NR)-Fpol(NR-1))/dpsi;

spline(psi,Fpol,NR,d1,dn,d2Fpol);

// Prepare Spline interpolation of Pres
d2Pres.resize(NR);
d1 = (Pres(2)-Pres(1))/dpsi;
dn = (Pres(NR)-Pres(NR-1))/dpsi;

spline(psi,Pres,NR,d1,dn,d2Pres);

// Prepare Spline interpolation of FFprime
d2FFprime.resize(NR);
d1 = (FFprime(2)-FFprime(1))/dpsi;
dn = (FFprime(NR)-FFprime(NR-1))/dpsi;

spline(psi,FFprime,NR,d1,dn,d2FFprime);

// Prepare Spline interpolation of Pprime
d2Pprime.resize(NR);
d1 = (Pprime(2)-Pprime(1))/dpsi;
dn = (Pprime(NR)-Pprime(NR-1))/dpsi;

spline(psi,Pprime,NR,d1,dn,d2Pprime);

// Prepare Spline interpolation of qpsi
d2qpsi.resize(NR);
d1 = (qpsi(2)-qpsi(1))/dpsi;
dn = (qpsi(NR)-qpsi(NR-1))/dpsi;

spline(psi,qpsi,NR,d1,dn,d2qpsi);

// Get helicity of Equilibrium
// B_tor is parallel to Ip if line trajectory derivative dZ/dphi < 0 at <-> "right-handed"
// outer midplane. Otherwise B_tor is anti-parallel to Ip <-> "left-handed".
int chk;
double psi_hel,dpsidr_hel;
double R_hel = max(lcfs(Range(1,toEnd,2))) - 10*dR;	// max of R coordinates from lcfs, just inside the separatrix ...
chk = get_psi(R_hel,0,psi_hel,dpsidr_hel,dummy);	// ... at midplande
helicity = sign(dpsidr_hel/get_Fpol(psi_hel));		// = sign(-dZ/dphi)
helicity_adjust = 1;

// Adjust helicity according to Ip and B_tor
if((Ip*Bt0 > 0 && helicity == -1) || (Ip*Bt0 < 0 && helicity == 1))
{
	helicity *= -1;
	helicity_adjust = -1;
}

// Identify the sign of toroidal and poloidal fields for use elsewhere.
//    btSign  	 of g-file Bt, rel to phi in cylindrical (R,phi,Z) 
//    bpSign 	 of g-file Iplasma, rel to phi in cylindrical (R,phi,Z)
// Standard EFIT writes flux to g-file with FLUX ALWAYS 
// INCREASING away from the magnetic axis (grad psi > 0)
// and the bpSign is used to give flux it correct physical sign.
// NON STANDARD EFITS MAY NOT FOLLOW THIS CONVENTION.
//btSign = sign(Bt0)	// +1.0 or -1.0 (real)
//bpSign = sign(Ip)	// +1.0 or -1.0 (real)

lcfs_th.resize(Nlcfs);
for(i=1;i<=Nlcfs;i++) lcfs_th(i) = polar_phi(lcfs(2*i-1) - RmAxis, lcfs(2*i) - ZmAxis);


return;
}

//-------------- get_psi -------------------------------------------------------------------------------------------------
int EFIT::get_psi(const double x1, const double x2, double& y, double& dy1, double& dy2, int flag, bool norm)
{
int chk = 0;
if(x1>R(NR) || x1<R(1) || x2>Z(NZ) || x2<Z(1))	{return -1;}

// get normalized poloidal Flux psi (should be chi in formulas!)
if(flag==0) splint_2D(R,Z,psiRZ,d2psi,NR,NZ,x1,x2,y,dy1,dy2);
else chk = bcuint(R,Z,Ca,dR,dZ,x1,x2,y,dy1,dy2);
if(chk == -1) {return -1;}

// normalize psi
if(norm==true) y = (y-psiAxis) / (psiSep-psiAxis);

// change helicity, if necessary
if(helicity_adjust == -1)
{
	dy1 *= -1;
	dy2 *= -1;
}

return 0;
}

//-------------- get_Fpol -------------------------------------------------------------------------------------------------
double EFIT::get_Fpol(const double x)
{
double y,dy;
splint(psi,Fpol,d2Fpol,NR,x,y,dy);
return y;
}

//-------------- get_Pres -------------------------------------------------------------------------------------------------
double EFIT::get_Pres(const double x)
{
double y,dy;
splint(psi,Pres,d2Pres,NR,x,y,dy);
return y;
}

//-------------- get_FFprime ----------------------------------------------------------------------------------------------
double EFIT::get_FFprime(const double x)
{
double y,dy;
splint(psi,FFprime,d2FFprime,NR,x,y,dy);
return y;
}

//-------------- get_Pprime -----------------------------------------------------------------------------------------------
double EFIT::get_Pprime(const double x)
{
double y,dy;
splint(psi,Pprime,d2Pprime,NR,x,y,dy);
return y;
}

//-------------- get_q ----------------------------------------------------------------------------------------------------
double EFIT::get_q(const double x)
{
double y,dy;
splint(psi,qpsi,d2qpsi,NR,x,y,dy);
return y;
}

//-------------- lcfs_RZ_nn -----------------------------------------------------------------------------------------------
// return nearest neighbor lcfs point to poloidal angle th
int EFIT::lcfs_RZ_nn(const double th, double& r, double& z)
{
double tmp = 10;
int i,idx;
for(i=1;i<=Nlcfs;i++)
{
	if (fabs(lcfs_th(i) - th) < tmp)
	{
		idx = i;
		tmp = fabs(lcfs_th(i) - th);
	}
}
r = lcfs(2*i-1);
z = lcfs(2*i);
return idx;
}


//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ------------------------ spline (mit Blitz-Arrays) -------------------------------------------------------------------
//Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
//x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative of the interpolating
//function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
//the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
//ypn are equal to 1e30 or larger, the routine is signaled to set the corresponding boundary
//condition for a natural spline, with zero second derivative on that boundary.
// Aufruf mit z.B.: spline(I,q1,NumberOfPoints,dq1,dqN,qr2);
void spline(Array<double,1>& x, Array<double,1>& y, int n, double yp1, double ypn, Array<double,1>& y2)
{
int i,k;
double p,qn,sig,un;
Array<double,1> u(n);

if (yp1 > 0.99e30) y2(1)=u(1)=0.0;	//The lower boundary condition is set either to be "natural"
else	//or else to have a specified first derivative.
{ 
	y2(1) = -0.5;
	u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1);
}

//This is the decomposition loop of the tridiagonal algorithm.
//y2 and u are used for temporary storage of the decomposed factors.
for (i=2;i<=n-1;i++)
{ 
	sig=(x(i)-x(i-1))/(x(i+1)-x(i-1));
	p=sig*y2(i-1)+2.0;
	y2(i)=(sig-1.0)/p;
	u(i)=(y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1));
	u(i)=(6.0*u(i)/(x(i+1)-x(i-1))-sig*u(i-1))/p;
}

if (ypn > 0.99e30) qn=un=0.0;	//The upper boundary condition is set either to be "natural" 
else	//or else to have a specified first derivative.
{ 
	qn=0.5;
	un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)));
}

y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0);

//This is the backsubstitution loop of the tridiagonal algorithm.
for (k=n-1;k>=1;k--)  y2(k)=y2(k)*y2(k+1)+u(k);
}

//--------------------- splint (mit Blitz-Arrays) ---------------------------------------------------------------------------
// Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai�s in order),
// and given the array y2a[1..n], which is the output from spline above, and given a value of
// x, this routine returns a cubic-spline interpolated value y(x) and dy/dx(x)
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, int n, double x, double& y, double& yx)
{
int klo,khi,k;
double h,b,a;

klo=1; 
khi=n;

while (khi-klo>1) 
{
	k=(khi+klo) >> 1;
	if (xa(k)>x) khi=k;
	else klo=k;
} //klo and khi now bracket the input value of x.

h=xa(khi)-xa(klo);
if (h==0.0) cout << "Bad xa input to routine splint" << endl; //The xa�s must be distinct.
a=(xa(khi)-x)/h;
b=(x-xa(klo))/h; //Cubic spline polynomial is now evaluated.

y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6.0;
yx=(ya(khi)-ya(klo))/h - ((3.0*a*a-1.0)*y2a(klo)-(3.0*b*b-1.0)*y2a(khi))*h/6.0;
}

//----------------- splint-2D (Blitz-Arrays) ------------------------------------------------------------------
//Given x1a, x2a, ya, n1, n2 and d2ydx1 as 1D-spline in x1 direction and
//given a desired interpolating point x1,x2; this routine returns an interpolated function value y
//by bicubic spline interpolation.
void splint_2D(Array<double,1>& x1a, Array<double,1>& x2a, Array<double,2>& ya, Array<double,2>& d2ydx1, int n1, int n2,
				double x1, double x2, double& y, double& dydx1, double& dydx2)
{
int j;
double d1,dn,dx2a;
double d2ydx1dx2;
Array<double,1> ytmp(n2+1),yytmp(n2+1),dyytmp(n2+1),dytmp(n2+1);
Range all = Range::all();
Array<double,1> slice_1,slice_2;	//Slice-array

//Perform n1 evaluations of the row splines (r-direction), using the one-dimensional spline evaluator splint.
for (j=1;j<=n2;j++) 
{
	slice_1.reference(ya(all,j));
	slice_2.reference(d2ydx1(all,j));
	splint(x1a,slice_1,slice_2,n1,x1,yytmp(j),dyytmp(j));
}

//Construct the one-dimensional column spline (z-direction) and evaluate it.
dx2a=x2a(2)-x2a(1);	// equidistant grid!!!!
d1=(yytmp(2)-yytmp(1))/dx2a;	// estimate of first derivative at boundary
dn=(yytmp(n2)-yytmp(n2-1))/dx2a;
spline(x2a,yytmp,n2,d1,dn,ytmp);  
splint(x2a,yytmp,ytmp,n2,x2,y,dydx2); 

// Construct the one-dimensional column spline of the derivative dydx1
d1=(dyytmp(2)-dyytmp(1))/dx2a;	// estimate of first derivative at boundary; does NOT have to be exact
dn=(dyytmp(n2)-dyytmp(n2-1))/dx2a;
spline(x2a,dyytmp,n2,d1,dn,dytmp);  
splint(x2a,dyytmp,dytmp,n2,x2,dydx1,d2ydx1dx2); 
}

//------------------ bcuderiv ---------------------------------------------------------------------------------------------
//Calculates the gradients y1 and y2 and the cross-derivative y12 numerically
//second order accuracy for interior and first order for boundary
//needs function y and grid-distances d1 and d2 (equidistant grid required)
void bcuderiv(Array<double,2>& y, double d1, double d2, Array<double,2>& y1, Array<double,2>& y2, Array<double,2>& y12)
{
int i,j;
const int N1 = y.rows();
const int N2 = y.cols();
const double d1d2 = d1*d2;

// Interior
for(i=2;i<N1;i++)
{
	for(j=2;j<N2;j++)
	{
		y1(i,j) = 0.5*(y(i+1,j)-y(i-1,j))/d1;
		y2(i,j) = 0.5*(y(i,j+1)-y(i,j-1))/d2;
		y12(i,j) = 0.25*(y(i+1,j+1)-y(i-1,j+1)-y(i+1,j-1)+y(i-1,j-1))/d1d2;
	}
}

// Boundaries
for(i=2;i<N1;i++) 
{
	// lower x2-boundary
	y1(i,1) = 0.5*(y(i+1,1)-y(i-1,1))/d1;
	y2(i,1) = (y(i,2)-y(i,1))/d2;
	y12(i,1) = 0.5*(y(i+1,2)-y(i-1,2)-y(i+1,1)+y(i-1,1))/d1d2;

	// upper x2-boundary
	y1(i,N2) = 0.5*(y(i+1,N2)-y(i-1,N2))/d1;
	y2(i,N2) = (y(i,N2)-y(i,N2-1))/d2;
	y12(i,N2) = 0.5*(y(i+1,N2)-y(i-1,N2)-y(i+1,N2-1)+y(i-1,N2-1))/d1d2;
}

for(j=2;j<N2;j++)
{
	// lower x1-boundary
	y1(1,j) = (y(2,j)-y(1,j))/d1;
	y2(1,j) = 0.5*(y(1,j+1)-y(1,j-1))/d2;
	y12(1,j) = 0.5*(y(2,j+1)-y(1,j+1)-y(2,j-1)+y(1,j-1))/d1d2;

	// upper x1-boundary
	y1(N1,j) = (y(N1,j)-y(N1-1,j))/d1;
	y2(N1,j) = 0.5*(y(N1,j+1)-y(N1,j-1))/d2;
	y12(N1,j) = 0.5*(y(N1,j+1)-y(N1-1,j+1)-y(N1,j-1)+y(N1-1,j-1))/d1d2;
}

// 4 corner points of the grid
y1(1,1) = (y(2,1)-y(1,1))/d1;
y2(1,1) = (y(1,2)-y(1,1))/d2;
y12(1,1) = (y(2,2)-y(1,2)-y(2,1)+y(1,1))/d1d2;

y1(1,N2) = (y(2,N2)-y(1,N2))/d1;
y2(1,N2) = (y(1,N2)-y(1,N2-1))/d2;
y12(1,N2) = (y(2,N2)-y(1,N2)-y(2,N2-1)+y(1,N2-1))/d1d2;

y1(N1,1) =(y(N1,1)-y(N1-1,1))/d1;
y2(N1,1) = (y(N1,2)-y(N1,1))/d2;
y12(N1,1) = (y(N1,2)-y(N1-1,2)-y(N1,1)+y(N1-1,1))/d1d2;

y1(N1,N2) = (y(N1,N2)-y(N1-1,N2))/d1;
y2(N1,N2) = (y(N1,N2)-y(N1,N2-1))/d2;
y12(N1,N2) = (y(N1,N2)-y(N1-1,N2)-y(N1,N2-1)+y(N1-1,N2-1))/d1d2;
}

//----------------- bcucof ------------------------------------------------------------------------------------------------
//Given arrays y[1..4], y1[1..4], y2[1..4], and y12[1..4], containing the function, gradients,
//and cross derivative at the four grid points of a rectangular grid cell (numbered counterclockwise
//from the lower left), and given d1 and d2, the length of the grid cell in the 1- and
//2-directions, this routine returns the table c[1..4][1..4] that is used by routine bcuint
//for bicubic interpolation.
void bcucof(Array<double,1>& y, Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12, double d1, double d2, 
			Array<double,2>& c)
{
static int wt[16][16]= {
{ 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
{-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
{2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
{0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
{0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
{-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
{9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
{-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
{2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
{-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
{4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1} };
int l,k,j,i;
double xx,d1d2,cl[16],x[16];
d1d2=d1*d2;

// Pack a temporary vector x.
for (i=1;i<=4;i++) 
{ 
	x[i-1]=y(i);
	x[i+3]=y1(i)*d1;
	x[i+7]=y2(i)*d2;
	x[i+11]=y12(i)*d1d2;
}
// Matrix multiply by the stored table.
for (i=0;i<=15;i++) 
{ 
	xx=0.0;
	for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
	cl[i]=xx;
}
l=0;

// Unpack the result into the output table.
for (i=1;i<=4;i++) 
	for (j=1;j<=4;j++) c(i,j)=cl[l++];
}

//------------------ bcuint ------------------------------------------------------------------------------------------------
//Bicubic interpolation. Input quantities are Ra and Za containing the grid coordinates,
//Ca contains the bcucof parameters for each grid square, stored with respect to the lower left corner.
//R and Z are the coordinates of the desired point for
//the interpolation. The interpolated function value is returned as y, and the interpolated
//gradient values as y1 and y2. This routine calls bcucof.
int bcuint(Array<double,1>& Ra, Array<double,1>& Za, Array<double,4>& Ca, double dR, double dZ,
			double R, double Z, double& y, double& y1, double& y2)
{
int i,j,k;
int NR = Ra.rows();
int NZ = Za.rows();
double t,u;
Array<double,2> c;
Range all = Range::all();

// Determine grid square where (R,Z) is in
// square includes lower and left boundary, but not upper or right boundary -> in next box included
j = int((R-Ra(1))/dR) + 1;
k = int((Z-Za(1))/dZ) + 1;
if(j == NR) j -= 1;	// exception: add outermost right boundary to square one to the left
if(k == NZ) k -= 1;	// exception: add outermost top boundary to square one down
if(j>NR || j<1 || k>NZ || k<1)	{cout << "bcuint: Point outside of grid" << endl; return -1;}

// Get the c�s.
c.reference(Ca(j,k,all,all));

// Interpolate
t=(R-Ra(j))/dR;
u=(Z-Za(k))/dZ;

y = y1 = y2 = 0.0;
for (i=4;i>=1;i--) 
{ 
	y = t*y + ((c(i,4)*u + c(i,3))*u + c(i,2))*u + c(i,1);
	y2 = t*y2 + (3.0*c(i,4)*u + 2.0*c(i,3))*u + c(i,2);
	y1 = u*y1 + (3.0*c(4,i)*t + 2.0*c(3,i))*t + c(2,i);
}

y1 /= dR;
y2 /= dZ;
return 0;
}

#endif //  EFIT_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

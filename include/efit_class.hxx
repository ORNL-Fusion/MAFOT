// EFIT Class handles all EFIT data read from g-files
// provides interpolation member functions for all data sets in g-file
// used by all MAFOT programs
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
#include <splines.hxx>
using namespace blitz;

// Prototypes  

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

#endif //  EFIT_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

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
	Array<double,1> Swall;	// length along the wall in m, Swall = 0 is the inner midplane; Swall goes ccw
	double Swall_max;		// total length of wall in m

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

	// external 3D wall
	bool use_3Dwall;		// default is false
	Array<int,1> Nwall3D;	// number of points for the Wall at each full integer toroidal angle phi = 0,...,359; each respective Array has twice the size!
	Array<double,2> wall3D;	// for each toroidal angle: Position of the Wall in R,Z plane; R,Z coordinates are alternating: R1,Z1,R2,Z2,R3,Z3,...

	// radial electric field
	bool hasEfield;				// E-field is available, default: False
	bool addScalPot;			// add scalar potential to GAMMA, default: False
	int NE;						// Number of Points in Er profile
	Array<double,1> Epsi;		// psi of Er profile
	Array<double,1> Eradial;	// (minor) radial electric field profile Er(psi) in [kV/m]
	Array<double,1> d2Eradial;	// Spline of Eradial in Epsi
	Array<double,1> scalPot;	// scalar potential Phi(psi) in [kV]
	Array<double,1> d2scalPot;	// Spline of scalar potential in Epsi

	// Temperature Profile
	bool hasTprofile;			// T profile is available, default: False
	int NT;						// Number of Points in T profile
	Array<double,1> Tpsi;		// psi of T profile
	Array<double,1> Tprofile;	// Temperature profile T(psi) in [keV]
	Array<double,1> d2Tprofile;	// Spline of Tprofile in Tpsi

// Konstruktoren
	EFIT();		// Default Constructor

// Member-Operatoren
	EFIT& operator =(const EFIT& EQD);		// Operator =

// Member-Funktionen
	void ReadData(const LA_STRING ShotNr, const LA_STRING ShotTime, double Raxis = 0, double Zaxis = 0);	// Reads EFIT data from g-file
	void ReadEfield(LA_STRING file);	// Read E-field profile from file
	void ReadTprofile(LA_STRING file);	// Read T profile from file

	// Interpolates psi and derivatives; flag: 0 = 2D splines, 1 = bicubic interpolation; norm specifies if psi is to be normalized
	int get_psi(const double x1, const double x2, double& y, double& dy1, double& dy2, int flag=1, bool norm=true);		// 0: ok	-1: Point outside of EFIT grid

	double get_Fpol(const double x);		// Spline interpolates Fpol
	double get_Pres(const double x);		// Spline interpolates Pres
	double get_FFprime(const double x);		// Spline interpolates FFprime
	double get_Pprime(const double x);		// Spline interpolates Pprime
	double get_q(const double x);			// Spline interpolates qpsi
	int lcfs_RZ_nn(const double th, double& r, double& z);
	void set3Dwall(LA_STRING wall_file);	// read 3D wall and set arrays
	double getEfield(const double x); 		// Spline interpolates Eradial
	double getscalPot(const double x); 		// Spline interpolates scalar Potential
	double getTprofile(const double x); 	// Spline interpolates Tprofile
	double wallDistance(double R, double Z);// returns disntance of R,Z from wall in m

}; //end of class

//------------------------ End of Class -----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Default Constructor -------------------------------------------------------------------------------------------
EFIT::EFIT()
{
TinyVector <int,1> index(1);
TinyVector <int,2> index2(1,1);
TinyVector <int,2> index2_10(1,0);
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
Swall.resize(Nwall);	Swall.reindexSelf(index);

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

use_3Dwall = false;
Nwall3D.resize(360);	// Nwall3D starts at index 0, which corresponds to angle phi = 0
wall3D.resize(500,360);		wall3D.reindexSelf(index2_10);	// second index (= angle) starts with 0 again

hasEfield = false;
addScalPot = false;
NE = 129;
Epsi.resize(NE);	Epsi.reindexSelf(index);
Eradial.resize(NE);	Eradial.reindexSelf(index);
d2Eradial.resize(NE);	d2Eradial.reindexSelf(index);
scalPot.resize(NE);		scalPot.reindexSelf(index);
d2scalPot.resize(NE);	d2scalPot.reindexSelf(index);

hasTprofile = false;
NT = 129;
Tpsi.resize(NE);	Tpsi.reindexSelf(index);
Tprofile.resize(NE);	Tprofile.reindexSelf(index);
d2Tprofile.resize(NE);	d2Tprofile.reindexSelf(index);

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
Swall.reference(EQD.Swall.copy());
Swall_max = EQD.Swall_max;

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

use_3Dwall = EQD.use_3Dwall;
Nwall3D.reference(EQD.Nwall3D.copy());
wall3D.reference(EQD.wall3D.copy());

hasEfield = EQD.hasEfield;
addScalPot = EQD.addScalPot;
NE = EQD.NE;
Epsi.reference(EQD.Epsi.copy());
Eradial.reference(EQD.Eradial.copy());
d2Eradial.reference(EQD.d2Eradial.copy());
scalPot.reference(EQD.scalPot.copy());
d2scalPot.reference(EQD.d2scalPot.copy());

hasTprofile = EQD.hasTprofile;
NT = EQD.NT;
Tpsi.reference(EQD.Tpsi.copy());
Tprofile.reference(EQD.Tprofile.copy());
d2Tprofile.reference(EQD.d2Tprofile.copy());
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
if(Path.right(1) == '/') Path = Path + filename;				// Path should be the full pathname of the gfile. If it is only the path (last char is a /), then add the conventional gfile name
ifstream in;
in.open(Path);
if(in.fail()==1) {cout << "Unable to open file " << Path << endl; exit(0);}

// Read dimensions, last two enties in first line
getline(in, stdummy);
vector<string> words;
words = split(stdummy);
NR = atoi(words[words.size()-2].c_str());	// Number of Points in R-direction
NZ = atoi(words[words.size()-1].c_str());	// Number of Points in Z-direction

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

// check if wall closes back into itself. If not, fix it
if((wall(2*Nwall-1) == wall(1)) && (wall(2*Nwall) == wall(2))) dummy = 0;	// just a dummy statement here. Only the else clause is relevant
else
{
	Nwall += 1;
	wall.resizeAndPreserve(2*Nwall);
	wall(2*Nwall-1) = wall(1);
	wall(2*Nwall) = wall(2);
}

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

//Commented / Modified by TL for HEAT (small tokamaks need a small multiplier or else R_hel is way inside core)
//double R_hel = max(lcfs(Range(1,toEnd,2))) - 10*dR;	// max of R coordinates from lcfs, just inside the separatrix ...
double R_hel = max(lcfs(Range(1,toEnd,2))) - 1*dR;	// max of R coordinates from lcfs, just inside the separatrix ...

chk = get_psi(R_hel,0,psi_hel,dpsidr_hel,dummy);	// ... at midplande
helicity = sign(dpsidr_hel/get_Fpol(psi_hel));		// = sign(-dZ/dphi)
helicity_adjust = 1;

// Adjust helicity according to Ip and B_tor
if((Ip*Bt0 > 0 && helicity == -1) || (Ip*Bt0 < 0 && helicity == 1))
{
	helicity *= -1;
	helicity_adjust = -1;

	//Added cout by TL for HEAT troubleshooting
	cout << "Flipping Helicity!" << endl;
	//cout << "Ip = " << Ip << endl;
	//cout << "Bt0 = " << Bt0 << endl;
	//cout << "helicity = "<< helicity << endl;
	//cout << "dpsidr_hel = " << dpsidr_hel << endl;
	//cout << "Fpol(psi_hel) = " << get_Fpol(psi_hel) << endl;
	//cout << "psi_hel = "<< psi_hel << endl;
	//cout << "R_hel = " << R_hel << endl;
	//cout << "Maximum LCFS Point = " << max(lcfs(Range(1,toEnd,2))) << endl;
	//cout << "dR = " << dR << endl;

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

// poloidal angle of LCFS
lcfs_th.resize(Nlcfs);
for(i=1;i<=Nlcfs;i++) lcfs_th(i) = polar_phi(lcfs(2*i-1) - RmAxis, lcfs(2*i) - ZmAxis);

// length along the wall; Swall = 0 is the inner midplane; Swall goes ccw
Swall.resize(Nwall);
int dir = 1;
double t,S0;
Swall(1) = sqrt((wall(1)-wall(2*Nwall-1))*(wall(1)-wall(2*Nwall-1)) + (wall(2)-wall(2*Nwall))*(wall(2)-wall(2*Nwall)));
if (Swall(1) > 0)
{
	S0 = wall(2*Nwall)/(wall(2*Nwall) - wall(2))*Swall(1);
	if (wall(2) < wall(2*Nwall)) dir = 1;	// ccw
	else dir = -1;							// cw
}
for(i=2;i<=Nwall;i++)
{
	Swall(i) = Swall(i-1) + sqrt((wall(2*i-1)-wall(2*(i-1)-1))*(wall(2*i-1)-wall(2*(i-1)-1)) + (wall(2*i)-wall(2*(i-1)))*(wall(2*i)-wall(2*(i-1))));	//length of curve in m
	if ((wall(2*i)*wall(2*i-2) <= 0) && (wall(2*i-1) < R0))
	{
		t = wall(2*i-2)/(wall(2*i-2) - wall(2*i));
		S0 = Swall(i-1) + t*(Swall(i) - Swall(i-1));
		if (wall(2*i) < wall(2*i-2)) dir = 1;	// ccw
		else dir = -1;							// cw
	}
}
Swall_max = Swall(Nwall);

// set direction and Swall = 0 location
for(i=1;i<=Nwall;i++)
{
	Swall(i) = dir*(Swall(i) - S0);
	if (Swall(i) < 0) Swall(i) += Swall_max;
	if (Swall(i) > Swall_max) Swall(i) -= Swall_max;
	if (fabs(Swall(i)) < 1e-12) Swall(i) = 0;
}

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

//-------------- set3Dwall ------------------------------------------------------------------------------------------------
// Read 3D wall file
void EFIT::set3Dwall(LA_STRING wall_file)
{
// Variables
int i,k;
LA_STRING line;
int count = 0;
int N,Nmax,Nangles;
vector<string> words;

// Input
ifstream in;
in.open(wall_file);
if(in.fail()==1) {cout << "Unable to open file " << wall_file << endl; exit(0);}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}

in.close();	// Important to start reading from the beginning of the file
in.clear(); // Important to clear EOF flag

in.open(wall_file);	// Open file again
for(i=1;i<=count;i++)	// Skip IO data rows
{
	in >> line;
}

// Read data
in >> line;
words = split(string(line));
Nmax = atoi(words[0].c_str());
if(words.size() > 1) Nangles = atoi(words[1].c_str());
else Nangles = 360;
Nwall3D.resize(Nangles);
wall3D.resize(2*Nmax,Nangles);

for(k=0;k<Nangles;k++)
{
	in >> N;
	Nwall3D(k) = N;
	for(i=1;i<=2*N;i++) in >> wall3D(i,k);
}

in.close();

// use the 3D wall now
use_3Dwall = true;
}

//---------------------- ReadEfield ---------------------------------------------------------------------------------------
void EFIT::ReadEfield(LA_STRING name)
{
int i;
LA_STRING line;
int count = 0;
int rows = 0;
double d1,dn,dpsi;
vector<string> words;
int column = 1;

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; return;}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}

// count number of columns in file
words = split(string(line));
column = words.size();
//cout << "Columns found in Er-file: " << column << endl;

// count the number of rows with data
while(in.eof()==0) // Last row is read twice --- can't be changed --- -> rows starts with 0 and is actual number of data rows in file
{
	in >> line;
	rows += 1;
}

in.close();	// Important to start reading from the beginning of the file
in.clear(); // Important to clear EOF flag

// resize
NE = rows;
Epsi.resize(NE);
Eradial.resize(NE);
scalPot.resize(NE);

in.open(name);	// Open file again
for(i=1;i<=count;i++) in >> line;	// Skip IO data rows

// Read data
if(column == 2)	// file has only E-field
{
	for(i=1;i<=rows;i++)
	{
		in >> Epsi(i);
		in >> Eradial(i);
	}
	scalPot = 0;
}
else	// file has E-field and scalar potential
{
	for(i=1;i<=rows;i++)
	{
		in >> Epsi(i);
		in >> Eradial(i);
		in >> scalPot(i);
	}
}
in.close();

// Prepare Spline interpolation of Eradial
d2Eradial.resize(NE);
dpsi = Epsi(2) - Epsi(1);
d1 = (Eradial(2)-Eradial(1))/dpsi;
dn = (Eradial(NE)-Eradial(NE-1))/dpsi;
spline(Epsi,Eradial,NE,d1,dn,d2Eradial);

// Prepare Spline interpolation of scalar potential
d2scalPot.resize(NE);
if(column == 2)
{
	d2scalPot = 0;
}
else
{
	dpsi = Epsi(2) - Epsi(1);
	d1 = (scalPot(2)-scalPot(1))/dpsi;
	dn = (scalPot(NE)-scalPot(NE-1))/dpsi;
	spline(Epsi,scalPot,NE,d1,dn,d2scalPot);
	addScalPot = true;
}

hasEfield = true;
}

//-------------- getEfield ------------------------------------------------------------------------------------------------
// spline interpolates Er(x)
double EFIT::getEfield(const double x)
{
double y,dy;
splint(Epsi,Eradial,d2Eradial,NE,x,y,dy);
return y;
}

//-------------- getscalPot -----------------------------------------------------------------------------------------------
// spline interpolates Er(x)
double EFIT::getscalPot(const double x)
{
double y,dy;
splint(Epsi,scalPot,d2scalPot,NE,x,y,dy);
return y;
}

//---------------------- ReadTprofile -------------------------------------------------------------------------------------
void EFIT::ReadTprofile(LA_STRING name)
{
int i;
LA_STRING line;
int count = 0;
int rows = 0;
double d1,dn,dpsi;
vector<string> words;
int column = 1;

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; return;}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}

// count number of columns in file
words = split(string(line));
column = words.size();
//cout << "Columns found in Er-file: " << column << endl;

// count the number of rows with data
while(in.eof()==0) // Last row is read twice --- can't be changed --- -> rows starts with 0 and is actual number of data rows in file
{
	in >> line;
	rows += 1;
}

in.close();	// Important to start reading from the beginning of the file
in.clear(); // Important to clear EOF flag

// resize
NT = rows;
Tpsi.resize(NT);
Tprofile.resize(NT);

in.open(name);	// Open file again
for(i=1;i<=count;i++) in >> line;	// Skip IO data rows

// Read data
for(i=1;i<=rows;i++)
{
	in >> Tpsi(i);
	in >> Tprofile(i);
}

in.close();

// Prepare Spline interpolation of Eradial
d2Tprofile.resize(NT);
dpsi = Tpsi(2) - Tpsi(1);
d1 = (Tprofile(2)-Tprofile(1))/dpsi;
dn = (Tprofile(NT)-Tprofile(NT-1))/dpsi;
spline(Tpsi,Tprofile,NT,d1,dn,d2Tprofile);

hasTprofile = true;
}

//-------------- getTprofile ---------------------------------------------------------------------------------------------
// spline interpolates T(x)
double EFIT::getTprofile(const double x)
{
double y,dy;
splint(Tpsi,Tprofile,d2Tprofile,NT,x,y,dy);
return y;
}

//-------------- wallDistance --------------------------------------------------------------------------------------------
// returns distance of point R,Z from the wall in m
// assumes wall closes back into itself, which is ensured in ReadData
double EFIT::wallDistance(double R, double Z)
{
int n;
double dR,dZ,s,d1,d2,d;
double dmin = 1e+16;
for(n=1;n<Nwall;n++)
{
	dR = wall(2*n+1) - wall(2*n-1);
	dZ = wall(2*n+2) - wall(2*n);
	if ((dR == 0) && (dZ == 0)) continue;
	s = ((R-wall(2*n-1))*dR + (Z-wall(2*n))*dZ)/(dR*dR + dZ*dZ);
	if ((s < 0) || (s > 1))
	{
		d1 = sqrt((R-wall(2*n-1))*(R-wall(2*n-1)) + (Z-wall(2*n))*(Z-wall(2*n)));
		d2 = sqrt((R-wall(2*n+1))*(R-wall(2*n+1)) + (Z-wall(2*n+2))*(Z-wall(2*n+2)));
		if (d1 < d2) d = d1;
		else d = d2;
	}
	else d = sqrt((R-wall(2*n-1)-s*dR)*(R-wall(2*n-1)-s*dR) + (Z-wall(2*n)-s*dZ)*(Z-wall(2*n)-s*dZ));
	if (d < dmin) dmin = d;
}
return dmin;
}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  EFIT_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

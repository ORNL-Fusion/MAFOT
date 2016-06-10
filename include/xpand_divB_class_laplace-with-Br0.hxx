// Class minimizes divB for B from xpand
// A.Wingen						1.4.16


// Define
//--------
#ifndef XPAND_DIVB_CLASS_INCLUDED
#define XPAND_DIVB_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <andi.hxx>
#include <vmec_class.hxx>
using namespace blitz;

// Prototypes
//-----------


// Golbal Parameters
//------------------


//--------- Begin Class divB ----------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Calculate divergenc B for B-field from xpand file or xpand-like grid
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class divB
{
private:
	// Member Variables
	Array<double,2> xpand_data;
	vector<LA_STRING> header;
	int header_lines;
	int N;

	// Member-Functions

public:
	// Member Variables
	int NR;
	int Np;
	int NZ;

	double dR;
	double dZ;
	double dphi;

	Array<double,1> R;
	Array<double,1> phi;
	Array<double,1> Z;
	Array<double,3> BR;
	Array<double,3> Bphi;
	Array<double,3> BZ;
	Array<double,3> div;	// divergence B; lives on the half grid in R and Z and on the full grid in phi
							// i.e. div(k,j,i) = divB_i+1/2_j+1/2_k in script

	// Constructors
	divB(int nr = 128, int np = 48, int nz = 128);								// Default Constructor
	divB(LA_STRING file);

	// Member-Operators
	divB& operator =(const divB& d);	// Operator =

	// Member-Functions
	void read(LA_STRING file);
	void ev(int order = 2);
	void write_xpand(LA_STRING file);
	void write_div(LA_STRING file);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
divB::divB(int nr, int np, int nz)
{
TinyVector <int,1> index(1);
TinyVector <int,3> index3(1,1,1);

NR = nr;
NZ = nz;
Np = np;

R.resize(NR);				R.reindexSelf(index);
phi.resize(NZ);				phi.reindexSelf(index);
Z.resize(NZ);				Z.reindexSelf(index);

BR.resize(Np, NZ, NR);		BR.reindexSelf(index3);
Bphi.resize(Np, NZ, NR);	Bphi.reindexSelf(index3);
BZ.resize(Np, NZ, NR);		BZ.reindexSelf(index3);

div.resize(Np, NZ-1, NR-1); div.reindexSelf(index3);
}

//----------------------------------------------------
// Read file Constructor
divB::divB(LA_STRING file)
{
TinyVector <int,1> index(1);
TinyVector <int,3> index3(1,1,1);

NR = 2;
NZ = 2;
Np = 1;

R.resize(NR);				R.reindexSelf(index);
phi.resize(NZ);				phi.reindexSelf(index);
Z.resize(NZ);				Z.reindexSelf(index);

BR.resize(Np, NZ, NR);		BR.reindexSelf(index3);
Bphi.resize(Np, NZ, NR);	Bphi.reindexSelf(index3);
BZ.resize(Np, NZ, NR);		BZ.reindexSelf(index3);

div.resize(Np, NZ-1, NR-1); div.reindexSelf(index3);

read(file);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
divB& divB::operator = (const divB& d)
{
if (this == &d) return(*this);	    // if: x=x

NR = d.NR;
NZ = d.NZ;
Np = d.Np;

dR = d.dR;
dZ = d.dZ;
dphi = d.dphi;

header = d.header;
header_lines = d.header_lines;
N = d.N;

xpand_data.reference(d.xpand_data);
R.reference(d.R);
phi.reference(d.phi);
Z.reference(d.Z);
BR.reference(d.BR);
Bphi.reference(d.Bphi);
BZ.reference(d.BZ);
div.reference(d.div);
return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- read ---------------------------------------------------------------------------------------------------------------
void divB::read(LA_STRING file)
{
// Variables
int i,j,k;
int idx;
int NPNZ = 0;
int NRNZ = 0;
int NRNP = 0;

header_lines = readFileHeader(file, header);
N = count_column(file);
readfile(file, N, xpand_data);

int M = xpand_data.rows();
for(i=1;i<=M;i++)
{
	if(xpand_data(i,1) == xpand_data(1,1)) NPNZ += 1;
	if(xpand_data(i,2) == xpand_data(1,2)) NRNZ += 1;
	if(xpand_data(i,3) == xpand_data(1,3)) NRNP += 1;
}
NZ = int(sqrt(NRNZ*NPNZ / NRNP) + 0.5);
NR = NRNZ / NZ;
Np = NPNZ / NZ;
cout << "NR = " << NR << " , Nphi = " << Np << " , NZ = " << NZ << endl;

R.resize(NR);
for(i=1;i<=NR;i++) R(i) = xpand_data(i,1);
dR = (R(NR) - R(1))/double(NR-1);

Z.resize(NZ);
for(j=1;j<=NZ;j++) Z(j) = xpand_data((j-1)*NR + 1,3);
dZ = (Z(NZ) - Z(1))/double(NZ-1);

phi.resize(Np);
if(Np > 1)
{
	for(k=1;k<=Np;k++) phi(k) = xpand_data((k-1)*NR*NZ + 1,2);
	dphi = (max(phi) - min(phi))/double(Np-1);
}
else
{
	phi(1) = xpand_data(1,2);
	dphi = 1;
}

BR.resize(Np, NZ, NR);
Bphi.resize(Np, NZ, NR);
BZ.resize(Np, NZ, NR);
for(k=1;k<=Np;k++)
	for(j=1;j<=NZ;j++)
		for(i=1;i<=NR;i++)
		{
			idx = i + (j-1)*NR + (k-1)*NR*NZ;
			BR(k,j,i) = xpand_data(idx,4);
			Bphi(k,j,i) = xpand_data(idx,5);
			BZ(k,j,i) = xpand_data(idx,6);
		}
}

//--- ev -------------------------------------------------------------------------------------------------------------------
void divB::ev(int order)
{
int i,j,k;
double BR_half,BZ_half, R_half, BR_half_plus1,BZ_half_plus1;
double Bp_half_m1,Bp_half_p1;
div.resize(Np, NZ-1, NR-1);

switch(order)
{
case 2:
	for(k=1;k<=Np;k++)
		for(j=1;j<=NZ-1;j++)
			for(i=1;i<=NR-1;i++)
			{
				BR_half_plus1 = 0.5*(BR(k,j+1,i+1) + BR(k,j,i+1));	// BR(k, j+0.5, i+1)
				BZ_half_plus1 = 0.5*(BZ(k,j+1,i+1) + BZ(k,j+1,i));	// BZ(k, j+1, i+0.5)
				BR_half = 0.5*(BR(k,j+1,i) + BR(k,j,i));	// BR(k, j+0.5, i)
				BZ_half = 0.5*(BZ(k,j,i+1) + BZ(k,j,i));	// BZ(k, j, i+0.5)
				R_half = 0.5*(R(i+1) + R(i));

				div(k,j,i) = (R(i+1)*BR_half_plus1 - R(i)*BR_half)/dR/R_half + (BZ_half_plus1 - BZ_half)/dZ;
			}
	break;
case 4:
	double BR_half_plus1half,BR_half_minushalf,BZ_half_plus1half,BZ_half_minushalf;
	double R_half_plus1, R_half_minus1;
	div = 0;
	for(k=1;k<=Np;k++)
		for(j=2;j<=NZ-2;j++)
			for(i=2;i<=NR-2;i++)
			{
				BR_half_plus1 = 0.5*(BR(k,j+1,i+1) + BR(k,j,i+1));	// BR(k, j+0.5, i+1)
				BZ_half_plus1 = 0.5*(BZ(k,j+1,i+1) + BZ(k,j+1,i));	// BZ(k, j+1, i+0.5)
				BR_half = 0.5*(BR(k,j+1,i) + BR(k,j,i));	// BR(k, j+0.5, i)
				BZ_half = 0.5*(BZ(k,j,i+1) + BZ(k,j,i));	// BZ(k, j, i+0.5)
				R_half = 0.5*(R(i+1) + R(i));
				R_half_plus1 = 0.5*(R(i+2) + R(i+1));
				R_half_minus1 = 0.5*(R(i) + R(i-1));
				BR_half_plus1half = 0.5*(0.5*(BR(k,j+1,i+2) + BR(k,j,i+2)) + BR_half_plus1);
				BR_half_minushalf = 0.5*(BR_half + 0.5*(BR(k,j+1,i-1) + BR(k,j,i-1)));
				BZ_half_plus1half = 0.5*(0.5*(BZ(k,j+2,i+1) + BZ(k,j+2,i)) + BZ_half_plus1);
				BZ_half_minushalf = 0.5*(BZ_half + 0.5*(BZ(k,j-1,i+1) + BZ(k,j-1,i)));

				div(k,j,i) = (-R_half_plus1*BR_half_plus1half + 8*R(i+1)*BR_half_plus1 - 8*R(i)*BR_half + R_half_minus1*BR_half_minushalf)/(6*dR)/R_half
							+ (-BZ_half_plus1half + 8*BZ_half_plus1 - 8*BZ_half + BZ_half_minushalf)/(6*dZ);
			}
	break;
default:
	cout << "Order not implemented. Available options: 2, 4" << endl;
}

// keep phi derivative at 2nd order, since it is almost zero anyway
if(Np > 1)
{
	if(k == 1) Bp_half_m1 = 0.25*(Bphi(Np,j+1,i) + Bphi(Np,j,i) + Bphi(Np,j+1,i+1) + Bphi(Np,j,i+1));
	else Bp_half_m1 = 0.25*(Bphi(k-1,j+1,i) + Bphi(k-1,j,i) + Bphi(k-1,j+1,i+1) + Bphi(k-1,j,i+1));
	if(k == Np) Bp_half_p1 = 0.25*(Bphi(1,j+1,i) + Bphi(1,j,i) + Bphi(1,j+1,i+1) + Bphi(1,j,i+1));
	else Bp_half_p1 = 0.25*(Bphi(k+1,j+1,i) + Bphi(k+1,j,i) + Bphi(k+1,j+1,i+1) + Bphi(k+1,j,i+1));

	div(k,j,i) += 0.5*(Bp_half_p1 - Bp_half_m1)/dphi/R_half;
}
}

//--- write_div ------------------------------------------------------------------------------------------------------------
void divB::write_div(LA_STRING file)
{
int i,j,k;
double R_half, Z_half;
ofstream out(file);
out.precision(16);
out << "# NR = " << NR << endl;
out << "# Np = " << Np << endl;
out << "# NZ = " << NZ << endl;
out << "# R[m]      \t phi[rad]   \t Z[m]      \t BR[T]      \t Bphi[T]    \t BZ[T]     \t divB " << endl;

for(k=1;k<=Np;k++)
	for(j=1;j<=NZ-1;j++)
		for(i=1;i<=NR-1;i++)
		{
			R_half = 0.5*(R(i+1) + R(i));
			Z_half = 0.5*(Z(j+1) + Z(j));
			out << R_half << "\t" << phi(k) << "\t" << Z_half << "\t" << BR(k,j,i) << "\t" << Bphi(k,j,i) << "\t" << BZ(k,j,i) << "\t" << div(k,j,i) << endl;
		}
}

//--- write_xpand ----------------------------------------------------------------------------------------------------------
void divB::write_xpand(LA_STRING file)
{
int i,j,k,idx,n;
double R_half, Z_half;
ofstream out(file);
out.precision(16);
for(n=0;n<header_lines;n++) out << header[n] << endl;

for(k=1;k<=Np;k++)
	for(j=1;j<=NZ;j++)
		for(i=1;i<=NR;i++)
		{
			idx = i + (j-1)*NR + (k-1)*NR*NZ;
			out << R(i) << "\t" << phi(k) << "\t" << Z(j) << "\t" << BR(k,j,i) << "\t" << Bphi(k,j,i) << "\t" << BZ(k,j,i);
			for(n=7;n<=N;n++) out << "\t" << xpand_data(idx,n);
			out << endl;
		}
}
//------------------------ End of Class divB ------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//--------- Begin Class POT -----------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Calculate a magnetic scalar potential that corrects divergence of B
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class POT
{
private:
	// Member Variables
	divB& DIVB;				// Only a Reference, not a copy
	Array<double,1> Raxis;
	Array<double,1> Zaxis;
	Array<double,3> p;		// intermediate variable for conjugate gradient method,  lives on the same grid as phi
	Array<double,3> v;		// laplace(p), intermediate variable for conjugate gradient method,  lives on the same grid as phi
	Array<double,3> Ft;		// support for iteration
	Array<double,3> Ftplus1;// support for iteration

	// Member-Functions

public:
	// Member Variables
	Array<double,3> phi;	// scalar potential, NOT toroidal angle; lives on the same grid as divB.div
	Array<double,3> lap;	// laplace(phi),  lives on the same grid as phi

	// Constructors
	POT(divB& d, VMEC& wout);

	// Member-Operators
	POT& operator =(const POT& P);	// Operator =

	// Member-Functions
	void laplace(Array<double,3>& phi, Array<double,3>& lap);
	void update_divB(int order = 2);
	void iterate(int tmax, double ftol = 1e-12, bool restart = false);
	void write(LA_STRING file);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
POT::POT(divB& d, VMEC& wout): DIVB(d)
{
double angle;
TinyVector <int,1> index1(1);
TinyVector <int,3> index3(1,1,1);
Raxis.resize(d.Np);		Raxis.reindexSelf(index1);
Zaxis.resize(d.Np);		Zaxis.reindexSelf(index1);
phi.resize(d.Np, d.NZ-1, d.NR-1);		phi.reindexSelf(index3);
lap.resize(d.Np, d.NZ-1, d.NR-1);		lap.reindexSelf(index3);
p.resize(d.Np, d.NZ-1, d.NR-1);			p.reindexSelf(index3);
v.resize(d.Np, d.NZ-1, d.NR-1);			v.reindexSelf(index3);
Ft.resize(d.Np, d.NZ-1, d.NR-1);		Ft.reindexSelf(index3);
Ftplus1.resize(d.Np, d.NZ-1, d.NR-1);	Ftplus1.reindexSelf(index3);

phi = 0;
lap = 0;
for(int k=1;k<=d.Np;k++)
{
	angle = (k-1)*pi2/d.Np;
	wout.get_axis(angle, Raxis(k), Zaxis(k));
}
}



//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
POT& POT::operator = (const POT& P)
{
if (this == &P) return(*this);	    // if: x=x
// References to divB remains unchanged! There is only one divB
Raxis.reference(P.Raxis);
Zaxis.reference(P.Zaxis);
phi.reference(P.phi);
lap.reference(P.lap);
p.reference(P.p);
v.reference(P.v);
Ft.reference(P.Ft);
Ftplus1.reference(P.Ftplus1);

return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// --- laplace ------------------------------------------------------------------------------------------------------------
void POT::laplace(Array<double,3>& phi, Array<double,3>& lap)
{
int i,j,k;
double phi1p, phi2p, phi1m, phi2m, R_half, phi1, phi2;
double Z_half,factor,phi_jm_i,phi_jm_ip,phi_jm_im;
double R_half_m,R_half_p,Z_half_m,Z_half_p;
double phi_jp_i,phi_jp_ip,phi_jp_im;
double phi_j_im,phi_j_ip;

// inside
// assume dphi/dvarphi = 0;  varphi = toroidal angle
for(k=1;k<=DIVB.Np;k++)
	for(j=2;j<=DIVB.NZ-2;j++)
		for(i=2;i<=DIVB.NR-2;i++)
		{
			R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));
			phi1p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j,i+1) + phi(k,j-1,i+1));
			phi1m = 0.25*(phi(k,j+1,i-1) + 2*phi(k,j,i-1) + phi(k,j-1,i-1));
			phi2p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j+1,i) + phi(k,j+1,i-1));
			phi2m = 0.25*(phi(k,j-1,i+1) + 2*phi(k,j-1,i) + phi(k,j-1,i-1));
			phi1 = 0.25*(phi(k,j+1,i) + 2*phi(k,j,i) + phi(k,j-1,i));
			phi2 = 0.25*(phi(k,j,i+1) + 2*phi(k,j,i) + phi(k,j,i-1));

			lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
						+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;
		}

// boundary
// enforce Br = 0  <=>  (R-Raxis)*dphi/dR = -(Z-Zaxis)*dphi/dZ
// 0.5*(R-Raxis)*(phi(k,j,i+1) - phi(k,j,i-1))/dR = -0.5*(Z-Zaxis)*(phi(k,j+1,i) - phi(k,j-1,i))/dZ
// d[i] = 0.5*(y[i+1] - y[i-1]) / dx				symmetric, 2nd order
// d[i] = 0.5*(-y[i+2] + 4*y[i+1] - 3*y[i]) / dx	upwards asymmetric, 2nd order
// d[i] = 0.5*(y[i-2] - 4*y[i-1] + 3*y[i]) / dx		downwards asymmetric, 2nd order
for(k=1;k<=DIVB.Np;k++)
{
	// lower Z boundary: set phi(k,j-1,i)
	j = 1;
	for(i=2;i<=DIVB.NR-2;i++)
	{
		R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));

		Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
		factor = DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k));
		phi_jm_i = (R_half - Raxis(k))*factor*(phi(k,j,i+1) - phi(k,j,i-1)) + phi(k,j+1,i);
		R_half_m = 0.5*(DIVB.R(i) + DIVB.R(i-1));
		if(i == 2) phi_jm_im = (R_half_m - Raxis(k))*factor*(-phi(k,j,i+1) + 4*phi(k,j,i) - 3*phi(k,j,i-1)) + phi(k,j+1,i-1);
		else phi_jm_im = (R_half_m - Raxis(k))*factor*(phi(k,j,i) - phi(k,j,i-2)) + phi(k,j+1,i-1);
		R_half_p = 0.5*(DIVB.R(i+2) + DIVB.R(i+1));
		if(i == DIVB.NR-2) phi_jm_ip = (R_half_p - Raxis(k))*factor*(phi(k,j,i-1) - 4*phi(k,j,i) + 3*phi(k,j,i+1)) + phi(k,j+1,i+1);
		else phi_jm_ip = (R_half_p - Raxis(k))*factor*(phi(k,j,i+2) - phi(k,j,i)) + phi(k,j+1,i+1);

		phi1p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j,i+1) + phi_jm_ip);
		phi1m = 0.25*(phi(k,j+1,i-1) + 2*phi(k,j,i-1) + phi_jm_im);
		phi2p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j+1,i) + phi(k,j+1,i-1));
		phi2m = 0.25*(phi_jm_ip + 2*phi_jm_i + phi_jm_im);
		phi1 = 0.25*(phi(k,j+1,i) + 2*phi(k,j,i) + phi_jm_i);
		phi2 = 0.25*(phi(k,j,i+1) + 2*phi(k,j,i) + phi(k,j,i-1));

		lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
					+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;
	}

	// upper Z boundary: set phi(k,j+1,i)
	j = DIVB.NZ-1;
	for(i=2;i<=DIVB.NR-2;i++)
	{
		R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));

		Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
		factor = DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k));
		phi_jp_i = -(R_half - Raxis(k))*factor*(phi(k,j,i+1) - phi(k,j,i-1)) + phi(k,j-1,i);
		R_half_m = 0.5*(DIVB.R(i) + DIVB.R(i-1));
		if(i == 2) phi_jp_im = -(R_half_m - Raxis(k))*factor*(-phi(k,j,i+1) + 4*phi(k,j,i) - 3*phi(k,j,i-1)) + phi(k,j-1,i-1);
		else phi_jp_im = -(R_half_m - Raxis(k))*factor*(phi(k,j,i) - phi(k,j,i-2)) + phi(k,j-1,i-1);
		R_half_p = 0.5*(DIVB.R(i+2) + DIVB.R(i+1));
		if(i == DIVB.NR-2) phi_jp_ip = -(R_half_p - Raxis(k))*factor*(phi(k,j,i-1) - 4*phi(k,j,i) + 3*phi(k,j,i+1)) + phi(k,j-1,i+1);
		else phi_jp_ip = -(R_half_p - Raxis(k))*factor*(phi(k,j,i+2) - phi(k,j,i)) + phi(k,j-1,i+1);

		phi1p = 0.25*(phi_jp_ip + 2*phi(k,j,i+1) + phi(k,j-1,i+1));
		phi1m = 0.25*(phi_jp_im + 2*phi(k,j,i-1) + phi(k,j-1,i-1));
		phi2p = 0.25*(phi_jp_ip + 2*phi_jp_i + phi_jp_im);
		phi2m = 0.25*(phi(k,j-1,i+1) + 2*phi(k,j-1,i) + phi(k,j-1,i-1));
		phi1 = 0.25*(phi_jp_i + 2*phi(k,j,i) + phi(k,j-1,i));
		phi2 = 0.25*(phi(k,j,i+1) + 2*phi(k,j,i) + phi(k,j,i-1));

		lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
					+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;
	}

	// min R boundary: set phi(k,j,i-1)
	i = 1;
	for(j=2;j<=DIVB.NZ-2;j++)
	{
		R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));

		Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
		factor = DIVB.dR/DIVB.dZ/(R_half - Raxis(k));
		phi_j_im = (Z_half - Zaxis(k))*factor*(phi(k,j+1,i) - phi(k,j-1,i)) + phi(k,j,i+1);
		Z_half_m = 0.5*(DIVB.Z(j) + DIVB.Z(j-1));
		if(j == 2) phi_jm_im = (Z_half_m - Zaxis(k))*factor*(-phi(k,j+1,i) + 4*phi(k,j,i) - 3*phi(k,j-1,i)) + phi(k,j-1,i+1);
		else phi_jm_im = (Z_half_m - Zaxis(k))*factor*(phi(k,j,i) - phi(k,j-2,i)) + phi(k,j-1,i+1);
		Z_half_p = 0.5*(DIVB.Z(j+2) + DIVB.Z(j+1));
		if(j == DIVB.NZ-2) phi_jp_im = (Z_half_p - Zaxis(k))*factor*(phi(k,j-1,i) - 4*phi(k,j,i) + 3*phi(k,j+1,i)) + phi(k,j+1,i+1);
		else phi_jp_im = (Z_half_p - Zaxis(k))*factor*(phi(k,j+2,i) - phi(k,j,i)) + phi(k,j+1,i+1);

		phi1p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j,i+1) + phi(k,j-1,i+1));
		phi1m = 0.25*(phi_jp_im + 2*phi_j_im + phi_jm_im);
		phi2p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j+1,i) + phi_jp_im);
		phi2m = 0.25*(phi(k,j-1,i+1) + 2*phi(k,j-1,i) + phi_jm_im);
		phi1 = 0.25*(phi(k,j+1,i) + 2*phi(k,j,i) + phi(k,j-1,i));
		phi2 = 0.25*(phi(k,j,i+1) + 2*phi(k,j,i) + phi_j_im);

		lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
					+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;
	}

	// max R boundary: set phi(k,j,i+1)
	i = DIVB.NR-1;
	for(j=2;j<=DIVB.NZ-2;j++)
	{
		R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));

		Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
		factor = DIVB.dR/DIVB.dZ/(R_half - Raxis(k));
		phi_j_ip = -(Z_half - Zaxis(k))*factor*(phi(k,j+1,i) - phi(k,j-1,i)) + phi(k,j,i-1);
		Z_half_m = 0.5*(DIVB.Z(j) + DIVB.Z(j-1));
		if(j == 2) phi_jm_ip = -(Z_half_m - Zaxis(k))*factor*(-phi(k,j+1,i) + 4*phi(k,j,i) - 3*phi(k,j-1,i)) + phi(k,j-1,i-1);
		else phi_jm_ip = -(Z_half_m - Zaxis(k))*factor*(phi(k,j,i) - phi(k,j-2,i)) + phi(k,j-1,i-1);
		Z_half_p = 0.5*(DIVB.Z(j+2) + DIVB.Z(j+1));
		if(j == DIVB.NZ-2) phi_jp_ip = -(Z_half_p - Zaxis(k))*factor*(phi(k,j-1,i) - 4*phi(k,j,i) + 3*phi(k,j+1,i)) + phi(k,j+1,i-1);
		else phi_jp_ip = -(Z_half_p - Zaxis(k))*factor*(phi(k,j+2,i) - phi(k,j,i)) + phi(k,j+1,i-1);

		phi1p = 0.25*(phi_jp_ip + 2*phi_j_ip + phi_jm_ip);
		phi1m = 0.25*(phi(k,j+1,i-1) + 2*phi(k,j,i-1) + phi(k,j-1,i-1));
		phi2p = 0.25*(phi_jp_ip + 2*phi(k,j+1,i) + phi(k,j+1,i-1));
		phi2m = 0.25*(phi_jm_ip + 2*phi(k,j-1,i) + phi(k,j-1,i-1));
		phi1 = 0.25*(phi(k,j+1,i) + 2*phi(k,j,i) + phi(k,j-1,i));
		phi2 = 0.25*(phi_j_ip + 2*phi(k,j,i) + phi(k,j,i-1));

		lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
					+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;
	}

	// min R, lower Z corner: set phi(k,j-1,i-1)
	i=1; j=1;
	R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));
	Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
	phi_j_im = (Z_half - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(-phi(k,j+2,i) + 4*phi(k,j+1,i) - 3*phi(k,j,i)) + phi(k,j,i+1);
	phi_jm_i = (R_half - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(-phi(k,j,i+2) + 4*phi(k,j,i+1) - 3*phi(k,j,i)) + phi(k,j+1,i);
	R_half_p = 0.5*(DIVB.R(i+2) + DIVB.R(i+1));
	phi_jm_ip = (R_half_p - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(phi(k,j,i+2) - phi(k,j,i)) + phi(k,j+1,i+1);
	Z_half_p = 0.5*(DIVB.Z(j+2) + DIVB.Z(j+1));
	phi_jp_im = (Z_half_p - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j+2,i) - phi(k,j,i)) + phi(k,j+1,i+1);
	R_half_m = R_half - DIVB.dR;
	phi_jm_im = (R_half_m - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(-phi(k,j,i+1) + 4*phi(k,j,i) - 3*phi_j_im) + phi_jp_im;

	phi1p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j,i+1) + phi_jm_ip);
	phi1m = 0.25*(phi_jp_im + 2*phi_j_im + phi_jm_im);
	phi2p = 0.25*(phi(k,j+1,i+1) + 2*phi(k,j+1,i) + phi_jp_im);
	phi2m = 0.25*(phi_jm_ip + 2*phi_jm_i + phi_jm_im);
	phi1 = 0.25*(phi(k,j+1,i) + 2*phi(k,j,i) + phi_jm_i);
	phi2 = 0.25*(phi(k,j,i+1) + 2*phi(k,j,i) + phi_j_im);

	lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
				+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;

	// min R, upper Z corner: set phi(k,j+1,i-1)
	i=1; j=DIVB.NZ-1;
	R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));
	Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
	phi_j_im = (Z_half - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j-2,i) - 4*phi(k,j-1,i) + 3*phi(k,j,i)) + phi(k,j,i+1);
	phi_jp_i = -(R_half - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(-phi(k,j,i+2) + 4*phi(k,j,i+1) - 3*phi(k,j,i)) + phi(k,j-1,i);
	Z_half_m = 0.5*(DIVB.Z(j) + DIVB.Z(j-1));
	phi_jm_im = (Z_half_m - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j,i) - phi(k,j-2,i)) + phi(k,j-1,i+1);
	R_half_p = 0.5*(DIVB.R(i+2) + DIVB.R(i+1));
	phi_jp_ip = -(R_half_p - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(phi(k,j,i+2) - phi(k,j,i)) + phi(k,j-1,i+1);
	Z_half_p = Z_half + DIVB.dZ;
	phi_jp_im = (Z_half_p - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j-1,i) - 4*phi(k,j,i) + 3*phi_jp_i) + phi_jp_ip;

	phi1p = 0.25*(phi_jp_ip + 2*phi(k,j,i+1) + phi(k,j-1,i+1));
	phi1m = 0.25*(phi_jp_im + 2*phi_j_im + phi_jm_im);
	phi2p = 0.25*(phi_jp_ip + 2*phi_jp_i + phi_jp_im);
	phi2m = 0.25*(phi(k,j-1,i+1) + 2*phi(k,j-1,i) + phi_jm_im);
	phi1 = 0.25*(phi_jp_i + 2*phi(k,j,i) + phi(k,j-1,i));
	phi2 = 0.25*(phi(k,j,i+1) + 2*phi(k,j,i) + phi_j_im);

	lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
				+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;

	// max R, lower Z corner: set phi(k,0,NR)
	i=DIVB.NR-1; j=1;
	R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));
	Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
	phi_jm_i = (R_half - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(phi(k,j,i-2) - 4*phi(k,j,i-1) + 3*phi(k,j,i)) + phi(k,j+1,i);
	phi_j_ip = -(Z_half - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(-phi(k,j+2,i) + 4*phi(k,j+1,i) - 3*phi(k,j,i)) + phi(k,j,i-1);
	R_half_m = 0.5*(DIVB.R(i) + DIVB.R(i-1));
	phi_jm_im = (R_half_m - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(phi(k,j,i) - phi(k,j,i-2)) + phi(k,j+1,i-1);
	Z_half_p = 0.5*(DIVB.Z(j+2) + DIVB.Z(j+1));
	phi_jp_ip = -(Z_half_p - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j+2,i) - phi(k,j,i)) + phi(k,j+1,i-1);
	R_half_p = R_half + DIVB.dR;
	phi_jm_ip = (R_half_p - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(phi(k,j,i-1) - 4*phi(k,j,i) + 3*phi_j_ip) + phi_jp_ip;

	phi1p = 0.25*(phi_jp_ip + 2*phi_j_ip + phi_jm_ip);
	phi1m = 0.25*(phi(k,j+1,i-1) + 2*phi(k,j,i-1) + phi_jm_im);
	phi2p = 0.25*(phi_jp_ip + 2*phi(k,j+1,i) + phi(k,j+1,i-1));
	phi2m = 0.25*(phi_jm_ip + 2*phi_jm_i + phi_jm_im);
	phi1 = 0.25*(phi(k,j+1,i) + 2*phi(k,j,i) + phi_jm_i);
	phi2 = 0.25*(phi_j_ip + 2*phi(k,j,i) + phi(k,j,i-1));

	lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
				+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;

	// max R, upper Z corner: set phi(k,NZ,NR)
	i=DIVB.NR-1; j=DIVB.NZ-1;
	R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));
	Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
	phi_j_ip = -(Z_half - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j-2,i) - 4*phi(k,j-1,i) + 3*phi(k,j,i)) + phi(k,j,i-1);
	phi_jp_i = -(R_half - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(phi(k,j,i-2) - 4*phi(k,j,i-1) + 3*phi(k,j,i)) + phi(k,j-1,i);
	Z_half_m = 0.5*(DIVB.Z(j) + DIVB.Z(j-1));
	phi_jm_ip = -(Z_half_m - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j,i) - phi(k,j-2,i)) + phi(k,j-1,i-1);
	R_half_m = 0.5*(DIVB.R(i) + DIVB.R(i-1));
	phi_jp_im = -(R_half_m - Raxis(k))*DIVB.dZ/DIVB.dR/(Z_half - Zaxis(k))*(phi(k,j,i) - phi(k,j,i-2)) + phi(k,j-1,i-1);
	Z_half_p = Z_half + DIVB.dZ;
	phi_jp_ip = -(Z_half_p - Zaxis(k))*DIVB.dR/DIVB.dZ/(R_half - Raxis(k))*(phi(k,j-1,i) - 4*phi(k,j,i) + 3*phi_jp_i) + phi_jp_im;

	phi1p = 0.25*(phi_jp_ip + 2*phi_j_ip + phi_jm_ip);
	phi1m = 0.25*(phi_jp_im + 2*phi(k,j,i-1) + phi(k,j-1,i-1));
	phi2p = 0.25*(phi_jp_ip + 2*phi_jp_i + phi_jp_im);
	phi2m = 0.25*(phi_jm_ip + 2*phi(k,j-1,i) + phi(k,j-1,i-1));
	phi1 = 0.25*(phi_jp_i + 2*phi(k,j,i) + phi(k,j-1,i));
	phi2 = 0.25*(phi_j_ip + 2*phi(k,j,i) + phi(k,j,i-1));

	lap(k,j,i) = (DIVB.R(i+1)*(phi1p - phi1) - DIVB.R(i)*(phi1 - phi1m))/R_half/DIVB.dR/DIVB.dR
				+ (phi2p - 2*phi2 + phi2m)/DIVB.dZ/DIVB.dZ;
}
}

// --- update_divB --------------------------------------------------------------------------------------------------------
void POT::update_divB(int order)
{
int i,j,k;
double phi1, phi2, phi1m, phi2m;

for(k=1;k<=DIVB.Np;k++)
	for(j=2;j<=DIVB.NZ-1;j++)
		for(i=2;i<=DIVB.NR-1;i++)
		{
			phi1 = 0.5*(phi(k,j,i) + phi(k,j-1,i));
			phi1m = 0.5*(phi(k,j,i-1) + phi(k,j-1,i-1));
			phi2 = 0.5*(phi(k,j,i) + phi(k,j,i-1));
			phi2m = 0.5*(phi(k,j-1,i) + phi(k,j-1,i-1));

			DIVB.BR(k,j,i) += (phi1 - phi1m)/DIVB.dR;
			DIVB.BZ(k,j,i) += (phi2 - phi2m)/DIVB.dZ;
		}
DIVB.ev(order);
}

// --- iterate --------------------------------------------------------------------------------------------------------
void POT::iterate(int tmax, double ftol, bool restart)
{
int i,j,k,t;
double rms,a,g,a1,g1;

if(restart)
{
	cout << "Restarting iterations" << endl;
}
else
{
	//phi = 0; 			// initial guess; already done in Constructor
	Ft = -DIVB.div;		// residual for iteration t;  use all points, instead of: Ft = where(abs(DIVB.div) < 1, -DIVB.div, 0);
	p = Ft;				// initial search direction
}

a1 = sum(Ft*Ft);		// same as the dot product of the two arrays
rms = sqrt(a1/Ft.size());
cout << "Iteration: 0" << "\t" << "RMS divB = " << rms << endl;
for(t=1;t<=tmax;t++)
{
	// update solution
	laplace(p, v);
	a = a1/sum(p*v);
	phi += a*p;

	// update residual & check convergence
	Ftplus1 = Ft - a*v;	// same as: laplace(phi, lap); Ftplus1 = -DIVB.div - lap;
	g1 = sum(Ftplus1*Ftplus1);
	rms = sqrt(g1/Ftplus1.size());
	if(rms <= ftol)
	{
		cout << "Iteration: " << t << "\t" << "RMS divB = " << rms << endl;
		break;
	}

	// update search direction
	g = g1/a1;
	p = Ftplus1 + g*p;

	// prepare next loop
	Ft = Ftplus1;
	a1 = g1;
	if(t%300 == 0) cout << "Iteration: " << t << "\t" << "RMS divB = " << rms << endl;
}

laplace(phi, lap);	// Update laplace(phi)
}

//--- write ---------------------------------------------------------------------------------------------------------------
void POT::write(LA_STRING file)
{
int i,j,k;
double R_half, Z_half;
ofstream out(file);
out.precision(16);
out << "# NR = " << DIVB.NR << endl;
out << "# Np = " << DIVB.Np << endl;
out << "# NZ = " << DIVB.NZ << endl;
out << "# R[m]      \t phi[rad]   \t Z[m]       \t scalar Potential " << endl;

for(k=1;k<=DIVB.Np;k++)
	for(j=1;j<=DIVB.NZ-1;j++)
		for(i=1;i<=DIVB.NR-1;i++)
		{
			R_half = 0.5*(DIVB.R(i+1) + DIVB.R(i));
			Z_half = 0.5*(DIVB.Z(j+1) + DIVB.Z(j));
			out << R_half << "\t" << DIVB.phi(k) << "\t" << Z_half << "\t" << phi(k,j,i) << "\t" << lap(k,j,i) << endl;
		}
}

//------------------------ End of Class POT -------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  XPAND_DIVB_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

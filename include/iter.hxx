// Header-File for the ITER Programs 
// Only Machine specific subroutines
// uses Nate Ferraro's M3D-C1 plasma response code output, fixed filename: C1.h5
// Plasma response can be for Equilibrium, or I-coils, or both
// C-coils and F-coils are not yet included in Plasma response
// ++++++ IMPORTANT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Libraries for the M3D-C1 routines only exist in Nate's u-drive account at GA
// use -Dm3dc1 when compiling -> this define activates this part of the code
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						22.06.11

// Define
//--------
#ifndef ITER_INCLUDED
#define ITER_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
int get_target(EFIT& EQD, IO& PAR);
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					   EFIT& EQD, IO& PAR, PARTICLE& FLT);
void point_along_wall(double swall, Array<double,1>& p, EFIT& EQD);	// defined in mafot.hxx

// ------------ Set Parameters for fortran --------------------------------------------------------------------------------
const int mxbands = 3;
const int mxloops = 9;
const int mxsegs = 52;
const int nturns = 1;

// Output of iterigeom_, set in: prep_perturbation()
int kuse[mxbands][mxloops];
int nbands;
int nloops[mxbands];
int nsegs[mxbands][mxloops];
double xs[mxbands][mxloops][mxsegs][3];
double dvs[mxbands][mxloops][mxsegs][4];
double curntw[mxbands][mxloops];

// ------------------- Fortran Common Blocks ------------------------------------------------------------------------------
extern "C" 
{
	extern struct{double pi,twopi,cir,rtd,dtr;} consts_;
	extern struct{double Iadj[mxbands][3];
				  double Icur[mxbands][mxloops];} currents_;
}

// ----------------- Fortran Routines -------------------------------------------------------------------------------------
extern "C" 
{
	void iterigeom_(int *kuse, const int *nbands, int nloops[], int *nsegs, double *xs, double *dvs, double *curntw);
	void polygonb_(const int *loopsdim, const int *segsdim, const int *nloops, int nsegs[], int kuse[],
					double *xs, double *dvs, double curnt[], 
					double *x, double *y, double *z, double *bx, double *by, double *bz);
}

// -------------- global Parameters ---------------------------------------------------------------------------------------
double bndy[7] = {3.95, 8.45, -4.6, 4.75, 6, -1, 2};	// Boundary Box: Rmin, Rmax, Zmin, Zmax		and line parameter: Rstart, Zdown, Zup
double bndy2[4] = { (bndy[2]-bndy[5])/(bndy[4]-bndy[1]),	// slope of line 1
					(bndy[4]*bndy[5]-bndy[1]*bndy[2])/(bndy[4]-bndy[1]),	// ordinate of line 1
					(bndy[3]-bndy[6])/(bndy[4]-bndy[1]),	// slope of line 2
					(bndy[4]*bndy[6]-bndy[1]*bndy[3])/(bndy[4]-bndy[1]) };	// ordinate of line 2

// nex box				old box
//	__					 ___
// |  \	<- line 2		|   |
// |   |				|	|
// |__/	<- line 1		|___|

// Real Wall
// Now read directly from EFIT file
//double wall_data[] = {4.05460000,   -2.50630000,	4.05460000,   -1.50000000,
//						4.05460000,   -0.4836000,	4.05460000,   0.53280000,
//						4.05460000,   1.54920000,	4.05460000,   2.56560000,
//						4.05460000,   3.58200000,	4.32000000,   4.32400000,
//						4.91280000,   4.71150000,	5.76290000,   4.53230000,
//						6.59610000,   3.89340000,	7.47630000,   3.08330000,
//						7.94290000,   2.40240000,	8.27940000,   1.68140000,
//						8.40350000,   0.63290000,	8.31540000,   -0.42150000,
//						7.90780000,   -1.34170000,	7.29200000,   -2.25700000,
//						6.27560000,   -3.04610000,	6.16250000,   -3.23970000,
//						5.97380000,   -3.28690000,	5.80690000,   -3.38700000,
//						5.67620000,   -3.53110000,	5.59300000,   -3.70700000,
//						5.56640000,   -3.89850000,	5.56440000,   -3.99970000,
//						5.55740000,   -4.09980000,	5.55740000,   -4.42530000,
//						5.55740000,   -4.42810000,	5.55736000,   -4.55894000,
//						5.54940000,   -4.55070000,	5.45480000,   -4.45610000,
//						5.36020000,   -4.38150000,	5.26550000,   -4.28680000,
//						5.24580000,   -3.98890000,	5.14250000,   -3.84210000,
//						4.99140000,   -3.74540000,	4.81490000,   -3.71300000,
//						4.63930000,   -3.75000000,	4.48570000,   -3.91300000,
//						4.38470000,  -3.90500000,	4.28380000,   -3.89710000,
//						4.18280000,   -3.88910000,	4.17395000,   -3.88824000,
//						4.22630000,   -3.77750000,	4.29550000,   -3.63990000,
//						4.37140000,   -3.48420000,	4.40040000,   -3.40880000,
//						4.44610000,  -3.28470000,	4.50950000,   -3.11880000,
//						4.50020000,  -2.94610000,	4.43480000,  -2.78610000,
//						4.31980000,  -2.65690000,	4.16650000,   -2.57300000};
//int Nwall = 54;
//Array<double,2> wall(wall_data,shape(Nwall,2),neverDeleteData);

// extern
#ifdef USE_SIESTA
	extern SIESTA SIES;
#endif
#ifdef USE_XFIELD
	extern XFIELD XPND;
#endif
#ifdef m3dc1
	extern M3DC1 M3D;
#endif

extern Array<double,4> field;
extern fakeIsland FISLD;

extern ofstream ofs2;

//---------------- getBfield ----------------------------------------------------------------------------------------------
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR)
{
int chk;
double X,Y,bx,by,bz;
double B_X,B_Y;
double sinp,cosp;

B_R = 0; B_phi = 0; B_Z = 0;

sinp = sin(phi);
cosp = cos(phi);

X = R*cosp;
Y = R*sinp;

chk = getBfield_general(R,Z,phi,B_R,B_Z,B_phi,EQD,PAR);
if(chk==-1) {return -1;}

B_X = 0;	B_Y = 0;
// I-coil perturbation field
if(PAR.useIcoil==1) 
{
	for(int i=0;i<nbands;i++)
	{
		bx = 0;	by = 0;	bz = 0;

		polygonb_(&mxloops, &mxsegs, &nloops[i], &nsegs[i][0], &kuse[i][0],
					&xs[i][0][0][0], &dvs[i][0][0][0], &curntw[i][0],
					&X, &Y, &Z, &bx, &by, &bz);
		B_X += bx;
		B_Y += by;
		B_Z += bz;
	}
}

// Transform B_perturbation = (B_X, B_Y, B_Z) to cylindrical coordinates and add
B_R += B_X*cosp + B_Y*sinp;
B_phi += -B_X*sinp + B_Y*cosp;

return 0;
}

//---------- prep_perturbation --------------------------------------------------------------------------------------------
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank, LA_STRING supPath)
{
int i,j;
int chk;
LA_STRING line;	// entire line is read by ifstream
ifstream in;

// Set common blocks parameters
consts_.pi = pi;
consts_.twopi = pi2;
consts_.cir = 360.0;
consts_.rtd = 360.0/pi2;
consts_.dtr = 1.0/consts_.rtd;

#ifdef m3dc1
	// Prepare loading M3D-C1
	if(PAR.response_field >= 0) chk = M3D.read_m3dc1sup(supPath);
	else chk = 0;
#else
	chk = 0;
#endif

// Read itersub.in file, if coils or M3D-C1 are on
if(PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	in.open(supPath + "itersup.in");
	if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open itersup.in file " << endl; EXIT;}

	for(i=1;i<=5;i++) in >> line;	// Skip 5 lines
	for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) in >> currents_.Icur[i][j];}		// Read coil currents

	in >> line;	// Skip line
	for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) in >> currents_.Iadj[i][j];}		// Read coil adjustments

	in.close();	// close file
	in.clear();	// reset ifstream for next use
}

#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{
		if(chk == -1) M3D.scale_from_coils(currents_.Icur[0], mxloops, mxloops*mxbands, currents_.Iadj[0][2]);	// no m3dc1sup.in file found -> scale from itersup.in, uses the top row
		M3D.load(PAR, mpi_rank);
	}
	else
	{
		if(mpi_rank < 1) cout << "Using g-file!" << endl;
		ofs2 << "Using g-file!" << endl;
	}
#else
	if(mpi_rank < 1) cout << "Using g-file!" << endl;
	ofs2 << "Using g-file!" << endl;
#endif

if(mpi_rank < 1) ofs2 << "Helicity = " << EQD.helicity << endl;
if(mpi_rank < 1) cout << "ITER-coil (0 = off, 1 = on): " << PAR.useIcoil << endl << endl;
ofs2 << "ITER-coil (0 = off, 1 = on): " << PAR.useIcoil << endl;

// Write Currents to log files (Check if corretly read in)
ofs2 << "Currents:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) ofs2 << currents_.Icur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "Adjustments:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) ofs2 << currents_.Iadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << endl;

// Set I-coil geometry
if(PAR.useIcoil==1) iterigeom_(&kuse[0][0],&nbands,&nloops[0],&nsegs[0][0],&xs[0][0][0][0],&dvs[0][0][0][0],&curntw[0][0]);
}

//---------------- get_target ---------------------------------------------------------------------------------------------
// finds target plate corners in EQD.wall data
int get_target(EFIT& EQD, IO& PAR)
{
int i;
double R = 0;
double Z = 0;
switch(PAR.which_target_plate)
{
case 1:			// Corner linear inner target and curve, t = 0 here
	R = 4.1787;
	Z = -3.8815;
	break;
case 2:			// Corner linear outer target and curve, t = 0 here
	R = 5.5563;		
	Z = -4.5318;
	break;
}
int endpt = 2*EQD.Nwall - 1;	// R-Index of last wall point
for(i=endpt;i>0;i-=2)
{
	if(fabs(EQD.wall(i) - R) < 1e-2 && fabs(EQD.wall(i+1) - Z) < 1e-2) return i;
}
cout << "Could not find target plate coordinates -> Abort!" << endl;
ofs2 << "Could not find target plate coordinates -> Abort!" << endl;
EXIT;
return -1;
}

//---------------- start_on_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t in cm
// inner target: t posivie on linear (right) part, negative on curve
// outer target: t positive on curve, negative on linear (left) part
// t = 0 specifies corner point where linear part and curve connect
// points of targets are directly taken from EFIT file
// Position (R0,Z0) of magnetic axis is required
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					   EFIT& EQD, IO& PAR, PARTICLE& FLT)
{
int i_p = 0;
int i_phi = 0;
int N = Np*Nphi;
int index;
double dp,dphi,t,x;
Array<double,1> p1(Range(1,2)),p2(Range(1,2)),p(Range(1,2)),d(Range(1,2));

// Grid stepsizes and t
if(Np == 1) tmax = tmin;
if(Nphi == 1) phimax = phimin;
if(N<=1) {dp = 0; dphi = 0;}
else
{
	dp = (tmax-tmin)/(N-1);
	dphi = (phimax-phimin)/(N-1);
}
if(dp==0) i_phi = i-1;
if(dphi==0) i_p = i-1;
if(dp!=0 && dphi!=0) 
{
	dp = (tmax-tmin)/double(Np-1);
	dphi = (phimax-phimin)/double(Nphi-1);
	i_phi = (i-1)%Nphi;
	i_p = int(double(i-1)/double(Nphi));
}
t = tmin + i_p*dp;	// t in cm

// Postion of Target-Plate
double R1 = 0, Z1 = 0;	// linear part left Corner
double R2 = 0, Z2 = 0;	// linear part right Corner
Array<double,1> R,Z,S;	// Curve
index = get_target(EQD,PAR);	// Get target plate corner positions
switch(PAR.which_target_plate)
{
case 1:	// inner target plate, linear part length 30.7752 cm, curved part length 165.6735 cm
	R1 = EQD.wall(index);	Z1 = EQD.wall(index+1);	// left Corner (4.1787, -3.8815)
	R2 = EQD.wall(index+2);	Z2 = EQD.wall(index+3);	// right Corner (4.4857, -3.903)
	
	N = 13;
	R.resize(Range(1,N));
	Z.resize(Range(1,N));
	S.resize(Range(1,N));
	R = EQD.wall(index), EQD.wall(index-2), EQD.wall(index-4), EQD.wall(index-6), EQD.wall(index-8), EQD.wall(index-10), 
		EQD.wall(index-12), EQD.wall(index-14), EQD.wall(index-16), EQD.wall(index-18), EQD.wall(index-20), 
		EQD.wall(index-22), EQD.wall(index-24);
	Z = EQD.wall(index+1), EQD.wall(index-1), EQD.wall(index-3), EQD.wall(index-5), EQD.wall(index-7), EQD.wall(index-9), 
		EQD.wall(index-11), EQD.wall(index-13), EQD.wall(index-15), EQD.wall(index-17), EQD.wall(index-19), 
		EQD.wall(index-21), EQD.wall(index-23);

	if(tmin < -165.6735 || tmax > 30.7752) ofs2 << "start_on_target: Warning, t out of range" << endl; 
	break;
case 2:	// outer target plate, linear part length 40.0237 cm, curved part length 175.0177 cm
	R2 = EQD.wall(index-2);	Z2 = EQD.wall(index-1);	// left Corner (5.2655, -4.2568) 
	R1 = EQD.wall(index);	Z1 = EQD.wall(index+1);	// right Corner (5.5563, -4.5318)

	N = 14;
	R.resize(Range(1,N));
	Z.resize(Range(1,N));
	S.resize(Range(1,N));
	R = EQD.wall(index), EQD.wall(index+2), EQD.wall(index+4), EQD.wall(index+6), EQD.wall(index+8), EQD.wall(index+10), 
		EQD.wall(index+12), EQD.wall(index+14), EQD.wall(index+16), EQD.wall(index+18), EQD.wall(index+20), 
		EQD.wall(index+22), EQD.wall(index+24), EQD.wall(index+26);
	Z = EQD.wall(index+1), EQD.wall(index+3), EQD.wall(index+5), EQD.wall(index+7), EQD.wall(index+9), EQD.wall(index+11), 
		EQD.wall(index+13), EQD.wall(index+15), EQD.wall(index+17), EQD.wall(index+19), EQD.wall(index+21), 
		EQD.wall(index+23), EQD.wall(index+25), EQD.wall(index+27);

	if(tmin < -40.0237 || tmax > 175.0177) ofs2 << "start_on_target: Warning, t out of range" << endl; 
	t *= -1;	// reverse t
	break;
case 0:
	point_along_wall(t, p, EQD);
	break;
default:
	ofs2 << "start_on_target: No target specified" << endl;
	EXIT;
	break;
}

if(target > 0)
{
	p1(1) = R1;		p1(2) = Z1;
	p2(1) = R2;		p2(2) = Z2;
	d = p2 - p1;	// positive R direction for inner target, negative otherwise;	t would have to be dimensionless in [0,1]
	d *= 0.01/sqrt(d(1)*d(1)+d(2)*d(2));	// d is now scaled for t in cm

	// Coordinates
	if(t>=0)
	{
		p = p1 + t*d;
		FLT.R = p(1);
		FLT.Z = p(2);
	}
	else	// t < 0
	{
		S(1) = 0;
		index = 1;
		for(int i=2;i<=N;i++)
		{
			S(i) = S(i-1) + sqrt((R(i)-R(i-1))*(R(i)-R(i-1)) + (Z(i)-Z(i-1))*(Z(i)-Z(i-1)));	//length of curve in m
			if(S(i) < -0.01*t) index = i;
			else break;
		}
		p1(1) = R(index);		p1(2) = Z(index);
		p2(1) = R(index+1);		p2(2) = Z(index+1);
		d = p2 - p1;
		x = (-0.01*t - S(index))/sqrt(d(1)*d(1)+d(2)*d(2));	// rescale t in m (like S); x is dimensionless in [0,1]
		p = p1 + x*d;
		FLT.R = p(1);
		FLT.Z = p(2);
	}
}
else
{
	FLT.R = p(1);
	FLT.Z = p(2);
}
FLT.phi = (phimin + dphi*i_phi)*rTOd;	// phi in deg
FLT.get_psi(p(1),p(2),FLT.psi);

FLT.Lc = 0;
FLT.psimin = 10;

if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}

if(PAR.which_target_plate == 2) t *= -1;	// undo reverse t
return t;
}
#endif // ITER_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

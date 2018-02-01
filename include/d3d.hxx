// Header-File for the DIII-D Programs 
// Only Machine specific subroutines
// uses Nate Ferraro's M3D-C1 plasma response code output, fixed filename: C1.h5
// Plasma response can be for Equilibrium, or I-coils, or both
// C-coils and F-coils are not yet included in Plasma response
// ++++++ IMPORTANT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Libraries for the M3D-C1 routines only exist in Nate's u-drive account at GA
// use -Dm3dc1 when compiling -> this define activates this part of the code
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						3.5.12

// Define
//--------
#ifndef D3D_INCLUDED
#define D3D_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
void prep_Bcoil_shiftTilt(void);
void Bcoil_shiftTilt_field(double X, double Y, double Z, double& Bx, double& By, double& Bz, EFIT& EQD);
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT);

// ------------ Set Parameters for fortran --------------------------------------------------------------------------------
const int nFc = 18;
const int nFlps = 2*nFc;
const int nIloops = 12;
const int nIsegs = 14;
const int nCloops = 6;
const int nCsegs = 10;
const int ntflimits = 4;
const int nmaxBusloops = 30;
const int nmaxBussegs = 30;
const int nBusloops = 11;

// Global Variables: have to be known during integration for perturbations, set in: prep_perturbation()
int kuseF[nFlps];
int kuseC[nCloops];
int kuseI[nIloops];
int nccsegs[nCloops];
int nicsegs[nIloops];
struct{double xs[nmaxBusloops][nmaxBussegs][3];
	   double dvs[nmaxBusloops][nmaxBussegs][4];
	   double curnt[nmaxBusloops];
	   int nsegs[nmaxBusloops];
	   int kbus[nmaxBusloops];
	   int nloops;} d3bus;	// nloops is a dummy and not used; see nBusloops

// ------------------- Fortran Common Blocks ------------------------------------------------------------------------------
// carefull, the order of the variables in the struct matters!
extern "C" 
{
	extern struct{double pi,twopi,cir,rtd,dtr;} consts_;
	extern struct{double fcshft[nFc],fashft[nFc],fctilt[nFc],fatilt[nFc],fcur[nFc];} d3pfer_;
	extern struct{double amat[nFlps][3][3];
				  double origin[nFlps][3]; 
				  double rcoil[nFlps], curlps[nFlps];} d3pflps_;
	extern struct{double dsbp,alfsbp,dthbp,alfthbp;
				  int iplasbp,ipdir,lbpol;} eqpol_;

	extern struct{double curntIc[nIloops]; 
				  double addanglIU, addanglIL, scaleIU, scaleIL;} d3icoil_;
	extern struct{double xis[nIloops][nIsegs][3]; 
				  double divs[nIloops][nIsegs][4];} d3iloops_;

	extern struct{double curntC[nCloops], curntw[nCloops]; 
				  double addanglC, scaleC;} d3ccoil_;
	extern struct{double xcs[nCloops][nCsegs][3];
				  double dcvs[nCloops][nCsegs][4];} d3cloops_;
	extern struct{double tflimits[ntflimits];
				  double rma, bma, dsbtc, alfsbtc, dthbtc, alfthbtc;
				  int ntfturns, ibcoil, lbtc, iripple;
				  double bcoil;
				  int usebcoilmds;} tfcoil_;
}

// ----------------- Fortran Routines -------------------------------------------------------------------------------------
extern "C" 
{
	void d3pfgeom_(int kuse[]);
	void d3igeom_(int kuse[]);
	void d3cgeom_(int kuse[]);
	void d3pferrb_(int kuse[], double *x, double *y, double *z, double *bxf, double *byf, double *bzf);
	void d3busnew2geom_(int kbus[], int *nloops, int nsegs[], double *xs, double *dvs, double curnt[]);
	void polygonb_(const int *loopsdim, const int *segsdim, const int *nloops, int nsegs[], int kuse[],
					double *xs, double *dvs, double curnt[], 
					double *x, double *y, double *z, double *bx, double *by, double *bz);
}

// -------------- global Parameters ---------------------------------------------------------------------------------------
double bndy[4] = {1.0, 2.4, -1.367, 1.36};	// Boundary

Array<double,2> Bcoil_Tilt_Matrix(Range(1,3),Range(1,3));
Array<double,1> Bcoil_Shift_Vector(Range(1,3));

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

// ***************************************************************
// NOTE: There are machine specific settings in IO_CLASS as well
// ***************************************************************

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
// F-coil perturbation field
if(PAR.useFcoil==1)
{
	bx = 0;	by = 0;	bz = 0;
	d3pferrb_(&kuseF[0], &X, &Y, &Z, &bx, &by, &bz);
	B_X += bx;
	B_Y += by;
	B_Z += bz;
}

// C-coil perturbation field
if(PAR.useCcoil==1)
{
	bx = 0;	by = 0;	bz = 0;
	polygonb_(&nCloops, &nCsegs, &nCloops, &nccsegs[0], &kuseC[0],
						  &d3cloops_.xcs[0][0][0], &d3cloops_.dcvs[0][0][0], &d3ccoil_.curntw[0], 
						  &X, &Y, &Z, &bx, &by, &bz);
	B_X += bx;
	B_Y += by;
	B_Z += bz;
}

// I-coil perturbation field
if(PAR.useIcoil==1)
{
	bx = 0;	by = 0;	bz = 0;
	polygonb_(&nIloops, &nIsegs, &nIloops, &nicsegs[0], &kuseI[0],
						  &d3iloops_.xis[0][0][0], &d3iloops_.divs[0][0][0], &d3icoil_.curntIc[0], 
						  &X, &Y, &Z, &bx, &by, &bz);
	B_X += bx;
	B_Y += by;
	B_Z += bz;
}

// Buswork perturbation field
if(PAR.useBuswork==1)
{
	bx = 0;	by = 0;	bz = 0;
	polygonb_(&nmaxBusloops, &nmaxBussegs, &nBusloops, &d3bus.nsegs[0], &d3bus.kbus[0],
						  &d3bus.xs[0][0][0], &d3bus.dvs[0][0][0], &d3bus.curnt[0],
						  &X, &Y, &Z, &bx, &by, &bz);
	B_X += bx;
	B_Y += by;
	B_Z += bz;
}

// shifted&tilted Bcoil perturbation field
if(PAR.useBcoil==1)
{
	bx = 0;	by = 0;	bz = 0;
	Bcoil_shiftTilt_field(X,Y,Z,bx,by,bz,EQD);
	B_X += bx;
	B_Y += by;
	B_Z += bz;
	B_phi -= EQD.Bt0*EQD.R0/R;	// subtract original B_phi component; new Bphi is in Bx & By
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

eqpol_.dsbp = 0.0;		// no shift !!!
eqpol_.dthbp = 0.0;		// no tilt !!!
eqpol_.alfsbp = 85.0*consts_.dtr;
eqpol_.alfthbp = 110.0*consts_.dtr;
eqpol_.ipdir = 1;		// positive Ip direction

d3ccoil_.scaleC = 1.0;
d3ccoil_.addanglC = 0.0;

d3icoil_.addanglIU = 0.0;
d3icoil_.addanglIL = 0.0;
d3icoil_.scaleIU = 1.0;
d3icoil_.scaleIL = 1.0;

tfcoil_.rma = EQD.R0;
tfcoil_.ntfturns = 144;
tfcoil_.bcoil = 0.5*1e+7*tfcoil_.rma*EQD.Bt0/tfcoil_.ntfturns;	// Magn. field at torus major radius: B = mu0 I*n/(2pi*R)
// the rest of this extern struct is not used in any way by MAFOT
for(i=0;i<ntflimits;i++) tfcoil_.tflimits[i] = 0;
tfcoil_.ibcoil = 1;
tfcoil_.usebcoilmds = 0;
tfcoil_.bma = 0;
tfcoil_.dsbtc = 0;
tfcoil_.alfsbtc = 0;
tfcoil_.dthbtc = 0;
tfcoil_.alfthbtc = 0;
tfcoil_.lbtc = 0;
tfcoil_.iripple = 0;

#ifdef m3dc1
	// Prepare loading M3D-C1
	if(PAR.response_field >= 0) chk = M3D.read_m3dc1sup(supPath);
	else chk = 0;
#else
	chk = 0;
#endif

// Read diiidsub.in file, if coils or M3D-C1 are on
if(PAR.useFcoil == 1 || PAR.useCcoil == 1 || PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	in.open(supPath + "diiidsup.in");
	if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open diiidsup.in file " << endl; EXIT;}

	for(i=1;i<=4;i++) in >> line;	// Skip 4 lines
	for(i=0;i<nFc;i++) in >> d3pfer_.fcur[i];		// Read F-coil currents

	in >> line;	// Skip line
	for(i=0;i<nCloops;i++) in >> d3ccoil_.curntC[i];		// Read C-coil currents

	in >> line;	// Skip line
	for(i=0;i<nIloops;i++) in >> d3icoil_.curntIc[i];		// Read I-coil currents

	in.close();	// close file
	in.clear();	// reset ifstream for next use
}

#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{

		if(chk == -1) M3D.scale_from_coils(d3icoil_.curntIc, nIloops, nIloops);	// no m3dc1sup.in file found -> scale from diiidsup.in file
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

if(mpi_rank < 1) cout << "F-coil: " << PAR.useFcoil << "\t" << "C-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << "\t" << "Bus Error: " << PAR.useBuswork << "\t" << "B-coil Error: " << PAR.useBcoil << endl << endl;
ofs2 << "F-coil: " << PAR.useFcoil << "\t" << "C-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << "\t" << "Bus Error: " << PAR.useBuswork << "\t" << "B-coil Error: " << PAR.useBcoil << endl << endl;

// Write I-coil currents to log files (Check if corretly read in)
ofs2 << "I-coil currents:" << endl;
for(i=0;i<nIloops/2;i++) ofs2 << d3icoil_.curntIc[i] << "\t";
ofs2 << endl;
for(i=nIloops/2;i<nIloops;i++) ofs2 << d3icoil_.curntIc[i] << "\t";
ofs2 << endl;

// Set F-coil geometry
if(PAR.useFcoil==1)
{
	for(i=0;i<nFlps;i++) kuseF[i]=0;	// kuseF is set inside the subroutine !?!
	d3pfgeom_(&kuseF[0]);	
}

// Set C-coil geometry
if(PAR.useCcoil==1)
{
	for(i=0;i<nCloops;i++) {if(d3ccoil_.curntC[i] != 0) kuseC[i]=1; else kuseC[i]=0;}
	for(i=0;i<nCloops;i++) nccsegs[i] = nCsegs;
	d3cgeom_(&kuseC[0]);
}

// Set I-coil geometry
if(PAR.useIcoil==1)
{
	for(i=0;i<nIloops;i++) {if(d3icoil_.curntIc[i] != 0) kuseI[i]=1; else kuseI[i]=0;}
	for(i=0;i<nIloops;i++) nicsegs[i] = nIsegs;
	d3igeom_(&kuseI[0]);
}

// Set buswork geometry
if(PAR.useBuswork==1)
{
	for(i=0;i<nBusloops;i++) d3bus.kbus[i] = 1;	// all loops are on
	d3busnew2geom_(&d3bus.kbus[0], &d3bus.nloops, &d3bus.nsegs[0], &d3bus.xs[0][0][0], &d3bus.dvs[0][0][0], &d3bus.curnt[0]);	// latest geometry only, since 2006
	ofs2 << "B-coil current [A] for Bus error field: " << tfcoil_.bcoil << endl;
}

// set Bcoil shift&tilt
if(PAR.useBcoil==1)
{
	prep_Bcoil_shiftTilt();
}

}

//---------- prep_Bcoil_shiftTilt --------------------------------------------------------------------------------------------
void prep_Bcoil_shiftTilt(void)
{
// inverse shift&tilt
// we want the new coordinates, within the shifted&tilted system, of the same point from the old system
// the forward shift&tilt would give the old coordinates (in the old system) of the shifted and tilted new system
const double shiftR = -5.7e-3;  			// in m;  INVERSE!
const double shift_tor_angle = -71 /rTOd;   // in deg, rhs
const double tilt_angle = -0.06 /rTOd;   	// in deg, from z-axis;  INVERSE!
const double tilt_tor_angle = -250 /rTOd;   // in deg, rhs

// tilting
Bcoil_Tilt_Matrix(1,1) = cos(tilt_tor_angle)*cos(tilt_tor_angle)*cos(tilt_angle) + sin(tilt_tor_angle)*sin(tilt_tor_angle);
Bcoil_Tilt_Matrix(1,2) = sin(tilt_tor_angle)*cos(tilt_tor_angle)*(cos(tilt_angle)-1);
Bcoil_Tilt_Matrix(1,3) = cos(tilt_tor_angle)*sin(tilt_angle);

Bcoil_Tilt_Matrix(2,1) = Bcoil_Tilt_Matrix(1,2);
Bcoil_Tilt_Matrix(2,2) = sin(tilt_tor_angle)*sin(tilt_tor_angle)*cos(tilt_angle) + cos(tilt_tor_angle)*cos(tilt_tor_angle);
Bcoil_Tilt_Matrix(2,3) = sin(tilt_tor_angle)*sin(tilt_angle);

Bcoil_Tilt_Matrix(3,1) = -Bcoil_Tilt_Matrix(1,3);
Bcoil_Tilt_Matrix(3,2) = -Bcoil_Tilt_Matrix(2,3);
Bcoil_Tilt_Matrix(3,3) = cos(tilt_angle);

// shifting
Bcoil_Shift_Vector(1) = shiftR*cos(shift_tor_angle);
Bcoil_Shift_Vector(2) = shiftR*sin(shift_tor_angle);
Bcoil_Shift_Vector(3) = 0;
}

//---------- Bcoil_shiftTilt_field --------------------------------------------------------------------------------------------
void Bcoil_shiftTilt_field(double X, double Y, double Z, double& Bx, double& By, double& Bz, EFIT& EQD)
{
double Rnew,Xnew,Ynew,phinew;
double B_phi_new,Bx_new,By_new;

Xnew = Bcoil_Tilt_Matrix(1,1)*X + Bcoil_Tilt_Matrix(1,2)*Y + Bcoil_Tilt_Matrix(1,3)*Z + Bcoil_Shift_Vector(1);
Ynew = Bcoil_Tilt_Matrix(2,1)*X + Bcoil_Tilt_Matrix(2,2)*Y + Bcoil_Tilt_Matrix(2,3)*Z + Bcoil_Shift_Vector(2);

Rnew = sqrt(Xnew*Xnew + Ynew*Ynew);
phinew = polar_phi(Xnew, Ynew);

// B_phi in the new system
B_phi_new = EQD.Bt0*EQD.R0/Rnew;

// B_phi_new in cartesian
Bx_new = -B_phi_new*sin(phinew);
By_new = B_phi_new*cos(phinew);
//Bz_new = 0;

// Bphi in the old coordinate system
Bx = Bcoil_Tilt_Matrix(1,1)*Bx_new + Bcoil_Tilt_Matrix(1,2)*By_new;
By = Bcoil_Tilt_Matrix(2,1)*Bx_new + Bcoil_Tilt_Matrix(2,2)*By_new;
Bz = Bcoil_Tilt_Matrix(3,1)*Bx_new + Bcoil_Tilt_Matrix(3,2)*By_new;
}

//---------------- start_on_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t=[0 1]
// end points of target are explicitly defined here
// Position (R0,Z0) of magnetic axis is required
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT)
{
int i_p=0;
int i_phi=0;
int N=Np*Nphi;
int target;
double dp,dphi,t;
Array<double,1> p1(Range(1,2)),p2(Range(1,2)),p(Range(1,2)),d(Range(1,2));
Array<double,1> R,Z,S;	// Curve
int idx;
double x,Smax;

// Magnetic Axis
//const double R0=EQD.RmAxis;
//const double Z0=EQD.ZmAxis;

// Grid stepsizes and t
if(Np == 1) tmax = tmin;
if(Nphi == 1) phimax = phimin;
if(N<=1) {dp=0; dphi=0;}
else
{
	dp=(tmax-tmin)/(N-1);
	dphi=(phimax-phimin)/(N-1);
}
if(dp==0) i_phi=i-1;
if(dphi==0) i_p=i-1;
if(dp!=0 && dphi!=0) 
{
	dp=(tmax-tmin)/double(Np-1);
	dphi=(phimax-phimin)/double(Nphi-1);
	i_phi=(i-1)%Nphi;
	i_p=int(double(i-1)/double(Nphi));
}
t=tmin+i_p*dp;
if(PAR.which_target_plate==1 && t<0) target=0;
else target=PAR.which_target_plate;

// Postion of Target-Plate
double R1 = 0,Z1 = 0;	// upper or left Point
double R2 = 0,Z2 = 0;	// lower or right Point
switch(target)
{
case 0:	// 19.59cm (same length as inner target) vertical wall above inner target, t = 0 -> -1, t=0 <=> P1 at inner target
	R1=1.016;	Z1=-1.223;
	R2=1.016;	Z2=-1.0271;
	break;
case 1:	// inner target plate
	R1=1.016;	Z1=-1.223;
	R2=1.153;	Z2=-1.363;
	break;
case 2:	// outer target plate
	R1=1.153;	Z1=-1.363;
	R2=1.372;	Z2=-1.363;	// normally R2=1.42, but that is inside the pump
	break;
case 3:	// 21.9cm (same length as outer target) horizontal shelf above pump to outer target
	R1=1.372;	Z1=-1.25;
	R2=1.591;	Z2=-1.25;
	break;
case 4:	// SAS divertor at upper outer divertor;  here t is dimensionless length along the wall;  t = 1 is same as Smax = 1.01693189 m; t = 0 is at the upper pump exit
	N = 36;
	R.resize(Range(1,N));
	Z.resize(Range(1,N));
	S.resize(Range(1,N));

	R = 1.372  ,  1.37167,  1.37003,  1.36688,  1.36719,  1.37178,
	        1.37224,  1.38662,  1.38708,  1.40382,  1.41127,  1.41857,
	        1.421  ,  1.48663,  1.4973 ,  1.49762,  1.49745,  1.49275,
	        1.4926 ,  1.49261,  1.49279,  1.4934 ,  1.4947 ,  1.49622,
	        1.47981,  1.48082,  1.48149,  1.48646,  1.49095,  1.50305,
	        1.59697,  1.6255 ,  1.63752,  1.647  ,  1.785  ,  2.07;
	Z = 1.31   ,  1.29238,  1.28268,  1.25644,  1.22955,  1.19576,
	        1.19402,  1.16487,  1.16421,  1.15696,  1.1573 ,  1.16132,
	        1.164  ,  1.2405 ,  1.23458,  1.23428,  1.23174,  1.2133 ,
	        1.21061,  1.20486,  1.20214,  1.19642,  1.18511,  1.1607 ,
	        1.12426,  1.12256,  1.12138,  1.11692,  1.11439,  1.11244,
	        1.09489,  1.0853 ,  1.07988,  1.077  ,  1.077  ,  1.04;
	Smax = 1.01693189;

	if(tmin < 0 || tmax > 1) ofs2 << "start_on_target: Warning, t out of range" << endl;
	S(1) = 0;
	idx = 1;
	for(int i=2;i<=N;i++)
	{
		S(i) = S(i-1) + sqrt((R(i)-R(i-1))*(R(i)-R(i-1)) + (Z(i)-Z(i-1))*(Z(i)-Z(i-1)));	//length of curve in m
		if(S(i) < Smax*t) idx = i;
		else break;
	}
	p1(1) = R(idx);		p1(2) = Z(idx);
	p2(1) = R(idx+1);		p2(2) = Z(idx+1);
	d = p2 - p1;
	x = (Smax*t - S(idx))/sqrt(d(1)*d(1)+d(2)*d(2));	// rescale t in m (like S); x is dimensionless in [0,1]
	p = p1 + x*d;
	break;
default:
	ofs2 << "No target specified" << endl;
	EXIT;
	break;
}

if(target <= 3)
{
	p1(1) = R1;	 p1(2) = Z1;
	p2(1) = R2;	 p2(2) = Z2;
	d = p2 - p1;
	if(target == 0) d(2) *= -1;

	// Coordinates
	p = p1 + t*d;
}

FLT.R = p(1);
FLT.Z = p(2);
FLT.phi = (phimin + dphi*i_phi)*rTOd;	// phi in deg
FLT.get_psi(p(1),p(2),FLT.psi);

FLT.Lc = 0;
FLT.psimin = 10;

if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}
return t;
}

#endif // D3D_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

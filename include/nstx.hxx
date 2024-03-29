// Header-File for the NSTX Drift Programs 
// Only Machine specific subroutines
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						22.6.11

// Define
//--------
#ifndef NSTX_INCLUDED
#define NSTX_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD);
void reorder_coil_currents(int mpi_rank=0);

// ------------ Set Parameters for fortran --------------------------------------------------------------------------------
const int mxbands = 1;
const int mxloops = 6;
const int mxsegs = 26;
const int nturns = 2;

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
	extern struct{double ECadj[mxbands][3];
				  double ECcur[mxbands][mxloops];} currents_;
}

// ----------------- Fortran Routines -------------------------------------------------------------------------------------
extern "C" 
{
	void nstxecgeom_(int *kuse, const int *nbands, int nloops[], int *nsegs, double *xs, double *dvs, double *curntw);
	void polygonb_(const int *loopsdim, const int *segsdim, const int *nloops, int nsegs[], int kuse[],
					double *xs, double *dvs, double curnt[], 
					double *x, double *y, double *z, double *bx, double *by, double *bz);
}

// -------------- global Parameters ---------------------------------------------------------------------------------------
double bndy[4] = {0.28, 1.6, -1.63, 1.63};	// Boundary Box: Rmin, Rmax, Zmin, Zmax

// extern
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

// Read nstxsub.in file, if coils or M3D-C1 are on
if(PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	in.open(supPath + "nstxsup.in");
	if(in.fail()==1) {if(mpi_rank < 1) cout << "Unable to open nstxsup.in file " << endl; EXIT;}

	for(i=1;i<=5;i++) {in >> line;} // Skip 5 lines
	for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) in >> currents_.ECcur[i][j];}		// Read coil currents

	in >> line;	// Skip line
	for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) in >> currents_.ECadj[i][j];}		// Read coil adjustments

	in >> line;	// last line, can be '&END'
	if (line.includes("use_new_coil_order"))
		if(line.includes("1"))
			reorder_coil_currents(mpi_rank);

	in.close();	// close file
	in.clear();	// reset ifstream for next use
}

#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{
		if(chk == -1) M3D.scale_from_coils(currents_.ECcur[0], mxloops, mxloops*mxbands, currents_.ECadj[0][2]);	// no m3dc1sup.in file found -> scale from nstxsup.in file, there is only one band
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

if(mpi_rank < 1) cout << "NSTX-coil (0 = off, 1 = on): " << PAR.useIcoil << endl << endl;
ofs2 << "NSTX-coil (0 = off, 1 = on): " << PAR.useIcoil << endl;

// Write Currents to log files (Check if corretly read in)
ofs2 << "Currents:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<mxloops;j++) ofs2 << currents_.ECcur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "Adjustments:" << endl;
for(i=0;i<mxbands;i++) {for(j=0;j<3;j++) ofs2 << currents_.ECadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << endl;

// Set ECoil geometry
if(PAR.useIcoil==1) nstxecgeom_(&kuse[0][0],&nbands,&nloops[0],&nsegs[0][0],&xs[0][0][0][0],&dvs[0][0][0][0],&curntw[0][0]);
}

//---------------- reorder_coil_currents -------------------------------------------------------------------------------------
// reorder coil currents to NSTX naming comvention: coil 1 is on top and numbering is clockwise
// otherwise coil oder is: first coil in centered around -30 deg and numbering is mathematical = counter-clockwise
void reorder_coil_currents(int mpi_rank)
{
if(mpi_rank < 1) cout << "Currents are reordered to NSTX naming convention:" << endl;
ofs2 << "Currents are reordered to NSTX naming convention:" << endl;
int i,j;
int order[6] = {2,1,0,5,4,3};
double temp[mxbands][mxloops];
for(i=0;i<mxbands;i++)
	for(j=0;j<mxloops;j++)
		temp[i][j] = currents_.ECcur[i][j];

for(i=0;i<mxbands;i++)
	for(j=0;j<mxloops;j++)
		currents_.ECcur[i][j] = temp[i][order[j]] ;
}

//---------------- point_along_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t in m
// inner target: t = Z
// outer target: t = R
// points of targets are set to fixed values here
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD)
{
double R1 = 0,Z1 = 0;	// upper or left Point
double R2 = 0,Z2 = 0;	// lower or right Point

if(target == 4 && t > 0.5712) target = 5;
if(target == 2 && t > 0.5712) target = 6;
if(target == 40 && t > 0.5712) target = 5;
if(target == 20 && t > 0.5712) target = 6;
if(target == 10 && t < 1.27) target = 80;
if(target == 30 && t > -1.27) target = 70;

switch(target)
{
case 1:	// upper inner target plate, length 40.66 cm
	R1 = 0.2794;	Z1 = 1.1714;
	R2 = 0.2794;	Z2 = 1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 2:	// upper outer target plate, length 27.33 cm
	R1 = 0.2979;	Z1 = 1.6034;
	R2 = 0.5712;	Z2 = 1.6034;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 3:	// lower inner target plate, length 40.66 cm
	R1 = 0.2794;	Z2 = -1.1714;
	R2 = 0.2794;	Z1 = -1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 4:	// lower outer target plate, length 27.33 cm
	R1 = 0.2979;	Z1 = -1.6034;
	R2 = 0.5712;	Z2 = -1.6034;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 5:	// lower outer inclined wall at R > 0.6 and small horizontal piece before that
	R1 = 0.617;		Z1 = -1.628;
	R2 = 1.0433;	Z2 = -1.4603;
	p(1) = t;
	if(t < R1) p(2) = Z1;
	else p(2) = (Z2-Z1)/(R2-R1)*t + (Z1*R2-R1*Z2)/(R2-R1);
	if(t < 0.5712 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
case 6:	// upper outer declined wall at R > 0.6
	R1 = 0.617;		Z1 = 1.628;
	R2 = 1.0433;	Z2 = 1.4603;
	p(1) = t;
	if(t < R1) p(2) = Z1;
	else p(2) = (Z2-Z1)/(R2-R1)*t + (Z1*R2-R1*Z2)/(R2-R1);
	if(t < 0.5712 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
case 10:	// NSTX-U upper inner target plate
	R1 = 0.4150;	Z1 = 1.27;
	R2 = 0.4150;	Z2 = 1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 20:	// NSTX-U upper outer target plate
	R1 = 0.4350;	Z1 = 1.6234;
	R2 = 0.5712;	Z2 = 1.6234;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 30:	// NSTX-U lower inner target plate
	R1 = 0.4150;	Z2 = -1.27;
	R2 = 0.4150;	Z1 = -1.578;
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = R1;	p(2) = t;
	break;
case 40:	// NSTX-U lower outer target plate
	R1 = 0.4350;	Z1 = -1.6234;
	R2 = 0.5712;	Z2 = -1.6234;
	if(t < R1 || t > R2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	p(1) = t;	p(2) = Z1;
	break;
case 70:	// NSTX-U lower inner inclined wall at R < 0.4
	R1 = 0.415;		Z1 = -1.27;
	R2 = 0.315;		Z2 = -1.05;
	p(2) = t;
	p(1) = (R2-R1)/(Z2-Z1)*t + (R1*Z2-Z1*R2)/(Z2-Z1);
	if(t < Z1 || t > Z2){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
case 80:	// NSTX-U upper inner inclined wall at R < 0.4
	R1 = 0.415;		Z1 = 1.27;
	R2 = 0.315;		Z2 = 1.05;
	p(2) = t;
	p(1) = (R2-R1)/(Z2-Z1)*t + (R1*Z2-Z1*R2)/(Z2-Z1);
	if(t < Z2 || t > Z1){ofs2 << "start_on_target: Warning, Coordinates out of range" << endl; EXIT;};
	break;
default:
	ofs2 << "start_on_target: No target specified" << endl;
	EXIT;
	break;
}
}

#endif // NSTX_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

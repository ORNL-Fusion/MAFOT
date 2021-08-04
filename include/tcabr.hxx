// Header-File for the TCABR Programs
// Only Machine specific subroutines
// uses Nate Ferraro's M3D-C1 plasma response code output, fixed filename: C1.h5
// Plasma response can be for Equilibrium, or I-coils, or both
// ++++++ IMPORTANT +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// uses separate Libraries for the M3D-C1 routines by Nate Ferraro
// use -Dm3dc1 when compiling -> this define activates this part of the code
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						03.08.21

// Define
//--------
#ifndef TCABR_INCLUDED
#define TCABR_INCLUDED

// Include
//--------

// --------------- Prototypes ---------------------------------------------------------------------------------------------
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);	// declared here, defined in mafot.hxx
int getBfield(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR);
void prep_perturbation(EFIT& EQD, IO& PAR, int mpi_rank=0, LA_STRING supPath="./");
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD);

// ------------ Set Parameters for fortran --------------------------------------------------------------------------------
// Internal coils on outer wall, like DIII-D I-coils, henceforth the I-coils
const int mxIbands = 3;
const int mxIloops = 24;
const int mxIsegs = 52;
const int niturns = 10;

// Centerpost coils (here still called EC, but also referred to as CP)
const int mxECbands = 3;
const int mxECloops = 24;
const int mxECsegs = 52;
const int ncturns = 10;

// Output of tcabrigeom_, set in: prep_perturbation()
int kuseI[mxIbands][mxIloops];
int nIbands;
int nIloops[mxIbands];
int nIsegs[mxIbands][mxIloops];
double xsI[mxIbands][mxIloops][mxIsegs][3];
double dvsI[mxIbands][mxIloops][mxIsegs][4];
double curntwI[mxIbands][mxIloops];

// Output of tcabrcpgeom_, set in: prep_perturbation()
int kuseEC[mxECbands][mxECloops];
int nECbands;
int nECloops[mxECbands];
int nECsegs[mxECbands][mxECloops];
double xsEC[mxECbands][mxECloops][mxECsegs][3];
double dvsEC[mxECbands][mxECloops][mxECsegs][4];
double curntwEC[mxECbands][mxECloops];

// ------------------- Fortran Common Blocks ------------------------------------------------------------------------------
extern "C" 
{
	extern struct{double pi,twopi,cir,rtd,dtr;} consts_;
	extern struct{double Iadj[mxIbands][3];
				  double Icur[mxIbands][mxIloops];} icurrents_;
	extern struct{double ECadj[mxECbands][3];
				  double ECcur[mxECbands][mxECloops];} eccurrents_;
}

// ----------------- Fortran Routines -------------------------------------------------------------------------------------
extern "C" 
{
	void tcabrigeom_(int *kuse, const int *nbands, int nloops[], int *nsegs, double *xs, double *dvs, double *curntw);
	void tcabrcpgeom_(int *kuse, const int *nbands, int nloops[], int *nsegs, double *xs, double *dvs, double *curntw);
	void polygonb_(const int *loopsdim, const int *segsdim, const int *nloops, int nsegs[], int kuse[],
					double *xs, double *dvs, double curnt[], 
					double *x, double *y, double *z, double *bx, double *by, double *bz);
}

// -------------- global Parameters ---------------------------------------------------------------------------------------
double bndy[4] = {0.43, 0.815, -0.245, 0.245};	// Boundary; EFIT boundary = {0.25, 1, -0.5, 0.5}

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
// CP-coil perturbation field
if(PAR.useCcoil==1)
{
	for(int i=0;i<nECbands;i++)
	{
		bx = 0;	by = 0;	bz = 0;

		polygonb_(&mxECloops, &mxECsegs, &nECloops[i], &nECsegs[i][0], &kuseEC[i][0],
					&xsEC[i][0][0][0], &dvsEC[i][0][0][0], &curntwEC[i][0],
					&X, &Y, &Z, &bx, &by, &bz);
		B_X += bx;
		B_Y += by;
		B_Z += bz;
	}
}

// I-coil perturbation field
if(PAR.useIcoil==1)
{
	for(int i=0;i<nIbands;i++)
	{
		bx = 0;	by = 0;	bz = 0;

		polygonb_(&mxIloops, &mxIsegs, &nIloops[i], &nIsegs[i][0], &kuseI[i][0],
					&xsI[i][0][0][0], &dvsI[i][0][0][0], &curntwI[i][0],
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
//LA_STRING line;	// entire line is read by ifstream
ifstream file;
string line;
vector<string> words;
int ECband,ECrow,Iband,Irow;

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

// Read tcabr.in file, if coils or M3D-C1 are on
if(PAR.useCcoil == 1 || PAR.useIcoil == 1 || (PAR.response_field > 0 && chk == -1))
{
	file.open(supPath + "tcabr.in");
	if(file.fail()==1) {if(mpi_rank < 1) cout << "Unable to open tcabr.in file " << endl; EXIT;}

	ECband = -1;
	Iband = -1;
	ECrow = -1;
	Irow = -1;

	while(getline(file, line))
	{
		if (line.length() == 0) continue;	// skip empty lines
		if (line[0] == '!') continue; 		// process and skip header lines
		words = split(line);
		if (words.size() == 0) continue;	// line just has blanks

		if (words[0].find("ECcur") != std::string::npos) {ECband += 1; ECrow = 0;}
		if (ECband == 3)
		{
			for(i=0;i<mxECbands;i++) {for(j=0;j<mxECloops;j++) eccurrents_.ECcur[i][j] = 0;}
		}
		if (ECband >3) {cout << "Error in reading CP-coil currents in tcabr.in file " << endl; break;}
		if (ECrow == 2)
		{
			for(j=12;j<18;j++) eccurrents_.ECcur[ECband][j] = atof(words[j-12].c_str());
			ECrow = 3;
		}
		if (ECrow == 1)
		{
			for(j=6;j<12;j++) eccurrents_.ECcur[ECband][j] = atof(words[j-6].c_str());
			ECrow = 2;
		}
		if (ECrow == 0)
		{
			for(j=0;j<6;j++) eccurrents_.ECcur[ECband][j] = atof(words[j+2].c_str());
			ECrow = 1;
		}
		if (words[0].find("ECadj") != std::string::npos)
		{
			for(j=0;j<3;j++) eccurrents_.ECadj[ECband][j] = atof(words[j+2].c_str());
		}

		if (words[0].find("Icur") != std::string::npos) {Iband += 1; Irow = 0;}
		if (Iband == 3)
		{
			for(i=0;i<mxIbands;i++) {for(j=0;j<mxIloops;j++) icurrents_.Icur[i][j] = 0;}
		}
		if (Iband >3) {cout << "Error in reading I-coil currents in tcabr.in file " << endl; break;}
		if (Irow == 2)
		{
			for(j=12;j<18;j++) icurrents_.Icur[Iband][j] = atof(words[j-12].c_str());
			Irow = 3;
		}
		if (Irow == 1)
		{
			for(j=6;j<12;j++) icurrents_.Icur[Iband][j] = atof(words[j-6].c_str());
			Irow = 2;
		}
		if (Irow == 0)
		{
			for(j=0;j<6;j++) icurrents_.Icur[Iband][j] = atof(words[j+2].c_str());
			Irow = 1;
		}
		if (words[0].find("Iadj") != std::string::npos) {for(j=0;j<3;j++) icurrents_.Iadj[Iband][j] = atof(words[j+2].c_str());}
	}
	file.close();	// close file
	file.clear();	// reset ifstream for next use
}

#ifdef m3dc1
	// Read C1.h5 file
	if(PAR.response_field >= 0)
	{
		if(chk == -1) M3D.scale_from_coils(icurrents_.Icur[1], mxIloops, mxIloops*mxIbands, icurrents_.Iadj[1][3]);	// no m3dc1sup.in file found -> scale from tcabr.in file, use only the lower band for scaling
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

if(mpi_rank < 1) cout << "CP-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << endl << endl;
ofs2 << "CP-coil: " << PAR.useCcoil << "\t" << "I-coil: " << PAR.useIcoil << endl << endl;

// Write Currents to log files (Check if corretly read in)
ofs2 << "CP-Currents:" << endl;
for(i=0;i<mxECbands;i++) {for(j=0;j<mxECloops;j++) ofs2 << eccurrents_.ECcur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "CP-Adjustments:" << endl;
for(i=0;i<mxECbands;i++) {for(j=0;j<3;j++) ofs2 << eccurrents_.ECadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << "I-Currents:" << endl;
for(i=0;i<mxIbands;i++) {for(j=0;j<mxIloops;j++) ofs2 << icurrents_.Icur[i][j] << "\t"; ofs2 << endl;}
ofs2 << "I-Adjustments:" << endl;
for(i=0;i<mxIbands;i++) {for(j=0;j<3;j++) ofs2 << icurrents_.Iadj[i][j] << "\t"; ofs2 << endl;}
ofs2 << endl;

// Set EC-coil geometry
if(PAR.useCcoil==1) tcabrcpgeom_(&kuseEC[0][0],&nECbands,&nECloops[0],&nECsegs[0][0],&xsEC[0][0][0][0],&dvsEC[0][0][0][0],&curntwEC[0][0]);

// Set I-coil geometry
if(PAR.useIcoil==1) tcabrigeom_(&kuseI[0][0],&nIbands,&nIloops[0],&nIsegs[0][0],&xsI[0][0][0][0],&dvsI[0][0][0][0],&curntwI[0][0]);
}

//---------------- point_along_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// t parametrizes the target with constant Phi and t=[0 1]
// end points of target are explicitly defined here
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
void point_along_target(int target, double t, Array<double,1>& p, EFIT& EQD)
{
Array<double,1> p1(Range(1,2)),p2(Range(1,2)),d(Range(1,2));
double R1 = 0,Z1 = 0;	// upper or left Point
double R2 = 0,Z2 = 0;	// lower or right Point

switch(target)
{
case 1:	// inner target plate
	R1 = 0.435;	Z1 = -0.24;
	R2 = 0.435;	Z2 = -0.14;
	break;
case 2:	// outer target plate
	R1 = 0.435;	Z1 = -0.24;
	R2 = 0.535;	Z2 = -0.24;
	break;
default:
	ofs2 << "No target specified" << endl;
	EXIT;
	break;
}

p1(1) = R1;	 p1(2) = Z1;
p2(1) = R2;	 p2(2) = Z2;
d = p2 - p1;

// Coordinates
p = p1 + t*d;
}

#endif // TCABR_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

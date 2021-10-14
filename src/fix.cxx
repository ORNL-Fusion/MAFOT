// Program searches for periodic fixed points in Poincare-Plot
// include drifts and time dependent perturbations
// derivatives are calculated numerically using a 5-point stencil
// 2-dimensional Newton method is used
// Fortran subroutines are used for perturbations
// A.Wingen						16.06.11

// Input: 1: Parameterfile	2: period of fixed point	3: praefix (optional)
// Output:	fixed points
//			log-file

// Define
//--------
#if defined(ITER)
	#define program_name "iterfix"
#elif defined(NSTX)
	#define program_name "nstxfix"
#elif defined(MAST)
	#define program_name "mastfix"
#elif defined(CMOD)
	#define program_name "cmodfix"
#elif defined(TCABR)
	#define program_name "tcabrfix"
#else
	#define program_name "dtfix"
#endif

// Include
//--------
#include <mafot.hxx>
#include <unistd.h>

// Prototypes  
//-----------
int gcd(int a, int b);
bool inside(PARTICLE& FLT, IO& PAR);
int newton2D(PARTICLE& FLT, double phistart, int periode);
int mapit_J(PARTICLE& FLT,  double J[], int itt, int Map=1);

// Switches
//----------

// Golbal Parameters
//------------------

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// Variables
int i,j,chk;
double x,y;
double dx,dy,xmin,ymin;
EFIT EQD;

// Period of fixed point
int periode;

// Use system time as seed(=idum) for random numbers
double now=zeit();
long idum=long(now);

// defaults
LA_STRING woutfile = "wout.nc";
LA_STRING xpandfile = "xpand.dat";
LA_STRING siestafile = "siesta.dat";
LA_STRING islandfile = "fakeIslands.in";
LA_STRING TprofileFile = "prof_t.dat";
LA_STRING NprofileFile = "prof_n.dat";
double f = 0;  // ratio of impurity to hydrogen ions
double zbar = 2;  // average over impurity ion charge states
double rc = 1;
double mc = 1;
bool use_collision = false;

// Command line input parsing
int c;
opterr = 0;
while ((c = getopt(argc, argv, "hX:V:S:I:")) != -1)
switch (c)
{
case 'h':
	cout << "usage: dtfix [-h] [-I island] [-S siesta] [-V wout] [-X xpand] file period [tag]" << endl << endl;
	cout << "Calculate the position of periodic points, like x-points and o-points." << endl << endl;
	cout << "positional arguments:" << endl;
	cout << "  file          Contol file (starts with '_')" << endl;
	cout << "  period        periodicity of point" << endl;
	cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
	cout << endl << "optional arguments:" << endl;
	cout << "  -h            show this help message and exit" << endl;
	cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
	cout << "  -S            filename for SIESTA; default, see below" << endl;
	cout << "  -V            filename for VMEC; default, see below" << endl;
	cout << "  -X            filename for XPAND; default, see below" << endl;
	cout << endl << "Examples:" << endl;
	cout << "  dtfix _fix.dat 1 blabla" << endl;
	cout << endl << "Infos:" << endl;
	cout << "  To use B-field from M3DC1, set response_field >= 0, and provide file in cwd:" << endl;
	cout << "    m3dc1sup.in    ->  location and scale factor for M3DC1 output C1.h5" << endl;
	cout << "  To use B-field from XPAND, set response_field = -3, and provide files in cwd:" << endl;
	cout << "    xpand.dat      ->  B-field on 3D grid from XPAND; use option -X to specify other filename" << endl;
	cout << "    wout.nc        ->  VMEC output; use option -V to specify other filename" << endl;
	cout << "  To use B-field from SIESTA, set response_field = -2, and provide file in cwd:" << endl;
	cout << "    siesta.dat     ->  B-field on 3D grid; use option -S to specify other filename" << endl;
	cout << "  To use B-field for mock-up islands, set response_field = -10, and provide file in cwd:" << endl;
	cout << "    fakeIslands.in ->  each line gives: Amplitude, pol. mode m, tor. mode n, phase [rad]" << endl;
	cout << "                       use option -I to specify other filename" << endl;
	cout << endl << "Current MAFOT version is: " << MAFOT_VERSION << endl;
	return 0;
case 'I':
	islandfile = optarg;
	break;
case 'S':
	siestafile = optarg;
	break;
case 'V':
	woutfile = optarg;
	break;
case 'X':
	xpandfile = optarg;
	break;
case '?':
	if (optopt == 'c')
		fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	else if (isprint (optopt))
		fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	else
		fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
	exit(0);
default:
	exit(0);
}
// Input file names
LA_STRING basename;
LA_STRING praefix = "";
if(argc==optind+3) praefix = "_" + LA_STRING(argv[optind+2]);
if(argc>=optind+2)
{
	periode = atoi(argv[optind+1]);
	basename = LA_STRING(argv[optind]);
}
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// log file
ofs2.open("log_" + LA_STRING(program_name) + "_" + LA_STRING(periode) + praefix + ".dat");
ofs2.precision(16);

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// additional parameters for IO
PAR.pv[6].name = "phistart";		PAR.pv[6].wert = PAR.phistart;
PAR.pv[7].name = "MapDirection";	PAR.pv[7].wert = PAR.MapDirection;
PAR.pv[8].name = "Ekin";			PAR.pv[8].wert = PAR.Ekin;
PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;

// set search grid
const int Nx = PAR.Nr;	// always use Nr = sqrt(N), because NR or NZ or others are given by separate positions in the _fix.dat file
const int Ny = PAR.Nth;	// same here
switch(PAR.create_flag)
{
case 3: case 4:
	PAR.pv[0].name = "psi-grid";		PAR.pv[0].wert = Nx;
	PAR.pv[1].name = "theta-grid";		PAR.pv[1].wert = Ny;
	PAR.pv[2].name = "psimin";			PAR.pv[2].wert = PAR.rmin;
	PAR.pv[3].name = "psimax";			PAR.pv[3].wert = PAR.rmax;
	PAR.pv[4].name = "thmin";			PAR.pv[4].wert = PAR.thmin;
	PAR.pv[5].name = "thmax";			PAR.pv[5].wert = PAR.thmax;
	dx = (PAR.rmax-PAR.rmin)/double(Nx-1);
	dy = (PAR.thmax-PAR.thmin)/double(Ny-1);
	xmin = PAR.rmin;
	ymin = PAR.thmin;
	break;
case 0:
	PAR.pv[0].name = "r-grid";			PAR.pv[0].wert = Nx;
	PAR.pv[1].name = "theta-grid";		PAR.pv[1].wert = Ny;
	PAR.pv[2].name = "rmin";			PAR.pv[2].wert = PAR.rmin;
	PAR.pv[3].name = "rmax";			PAR.pv[3].wert = PAR.rmax;
	PAR.pv[4].name = "thmin";			PAR.pv[4].wert = PAR.thmin;
	PAR.pv[5].name = "thmax";			PAR.pv[5].wert = PAR.thmax;
	dx = (PAR.rmax-PAR.rmin)/double(Nx-1);
	dy = (PAR.thmax-PAR.thmin)/double(Ny-1);
	xmin = PAR.rmin;
	ymin = PAR.thmin;
	break;
case 5:
	PAR.pv[0].name = "R-grid";			PAR.pv[0].wert = Nx;
	PAR.pv[1].name = "Z-grid";			PAR.pv[1].wert = Ny;
	PAR.pv[2].name = "Rmin";			PAR.pv[2].wert = PAR.Rmin;
	PAR.pv[3].name = "Rmax";			PAR.pv[3].wert = PAR.Rmax;
	PAR.pv[4].name = "Zmin";			PAR.pv[4].wert = PAR.Zmin;
	PAR.pv[5].name = "Zmax";			PAR.pv[5].wert = PAR.Zmax;
	dx = (PAR.Rmax-PAR.Rmin)/double(Nx-1);
	dy = (PAR.Zmax-PAR.Zmin)/double(Ny-1);
	xmin = PAR.Rmin;
	ymin = PAR.Zmin;
	break;
default:
	cout << "Unknown create points option" << endl;
	EXIT;
}

// Read EFIT-data
double Raxis = 0, Zaxis = 0;
#ifdef USE_XFIELD
if(PAR.response_field == -3)
{
	cout << "Read VMEC file" << endl;
	ofs2 << "Read VMEC file" << endl;
	vmec.read(woutfile);
	vmec.n0only = true;
	vmec.get_axis(PAR.phistart/rTOd,Raxis,Zaxis);	// need axisymmetric axis only
	vmec.n0only = false;
}
#endif

#ifdef m3dc1
if(PAR.response_field == 0 || PAR.response_field == 2)
{
	M3D.read_m3dc1sup();
	M3D.open_source(PAR.response, PAR.response_field, -1);
	Raxis = M3D.RmAxis;
	Zaxis = M3D.ZmAxis;
	M3D.unload();
}
#endif

EQD.ReadData(EQD.Shot,EQD.Time,Raxis,Zaxis);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prepare_common_perturbations(EQD,PAR,0,siestafile,xpandfile,islandfile);
prep_perturbation(EQD,PAR);

// Prepare collisions
COLLISION COL;
if (use_collision) COL.init(TprofileFile, NprofileFile, f, zbar, PAR.Zq, PAR.Mass, rc, mc);
// Prepare particles
PARTICLE FLT(EQD,PAR,COL);

// Use Boundary Box
simpleBndy = 1;

// Output
LA_STRING filenameout = "fix_" + LA_STRING(periode) + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(6);
var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "period";  var[3] = "psi";  var[4] = "theta[rad]";  var[5] = "r[m]";
PAR.writeiodata(out,bndy,var);

// get rational surface
double qrat = 1;
int nrat = 1;
int mrat = periode;
int mngcd;
if(periode > 1)
{
	qrat = EQD.get_q(0.5*(PAR.rmin + PAR.rmax));
	nrat = int(mrat/qrat + 0.5);
	ofs2 << "Searching at rational surface q = " << mrat << "/" << nrat << endl;
	cout << "Searching at rational surface q = " << mrat << "/" << nrat << endl;
	mngcd = gcd(mrat,nrat);
	if(mngcd > 1)
	{
		periode /= mngcd;
		ofs2 << "Periodic points are heteroclinic. Use new period = " << periode << endl;
		cout << "Periodic points are heteroclinic. Use new period = " << periode << endl;
	}
}

Array<double,1> Rall(2*mrat), Zall(2*mrat);
int idx_all = 0;
Array<double,1> R(periode),Z(periode),psi(periode),theta(periode),r(periode);
double dist;
ofs2 << Nx << " rows, done:" << endl;
for(i=0;i<Nx;i++)	//r or R
{
	x = xmin + i*dx;
	for(j=0;j<Ny;j++)	//theta or Z
	{
		y = ymin + j*dy;
		switch(PAR.create_flag)
		{
		case 0:		// r-theta grid
			FLT.convertRZ(y,x);
			break;
		case 3: 	// psi-theta grid
			FLT.psi = x;
			FLT.getRZ(x, y, FLT.R, FLT.Z);
			break;
		case 4:		// psi-theta random
			FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,2);
			break;
		case 5:		// R,Z
			FLT.R = x;
			FLT.Z = y;
			break;
		}

		chk = newton2D(FLT,PAR.phistart,periode);
		if(chk == -1) continue;

		if(inside(FLT,PAR))
		{
			if(periode == 1)
			{
				out << FLT.R << "\t" << FLT.Z << "\t" << -periode << "\t" << FLT.psi << "\t" << FLT.get_theta() << "\t" << FLT.get_r() << endl;
				cout << "Program terminated normally" << endl;
				ofs2 << "Program terminated normally" << endl;
				#ifdef m3dc1
				if(PAR.response_field >= 0) M3D.unload();
				#endif
				return 0; 
			}

			// check if point is already in list
			chk = 0;
			for(int k=0;k<PAR.N;k++)
			{
				if(Rall(k) == 0) break;
				dist = sqrt((Rall(k)-FLT.R)*(Rall(k)-FLT.R) + (Zall(k)-FLT.Z)*(Zall(k)-FLT.Z));
				if(dist < 1e-5) {chk = -1; break;}	// same point already in list
			}
			if(chk < 0) continue;

			R(0) = FLT.R; Z(0) = FLT.Z; psi(0) = FLT.psi; r(0) = FLT.get_r(); theta(0) = FLT.get_theta();

			// map out all intermediate points and check if periodic point actually has the period periode or a smaller one
			for(int k=1;k<periode;k++)
			{
				chk = FLT.mapstep();
				if(chk < 0) {break;}	// particle has left system
				dist = sqrt((R(0)-FLT.R)*(R(0)-FLT.R) + (Z(0)-FLT.Z)*(Z(0)-FLT.Z));
				if(dist < 1e-5) {chk = -1; break;}	// same point with k < periode, so not a periodic point we are looking for
				R(k) = FLT.R; Z(k) = FLT.Z; psi(k) = FLT.psi; r(k) = FLT.get_r(); theta(k) = FLT.get_theta();
			}
			if(chk < 0) continue;

			// Check if periodic point is an O-point or X-point
			FLT.R = R(0) + 1e-5; FLT.Z = Z(0) + 1e-5;	// perturb the periodic point
			chk = FLT.mapit(10*periode);
			dist = sqrt((R(0)-FLT.R)*(R(0)-FLT.R) + (Z(0)-FLT.Z)*(Z(0)-FLT.Z));
			if(dist > 1e-2) chk = -1;		// point has moved away by more than a cm
			if(chk >= 0) chk = 1;

			// chk == -1: X-point; chk == 1: O-point
			for(int k=0;k<periode;k++)
			{
				Rall(idx_all) = R(k); Zall(idx_all) = Z(k);
				idx_all += 1;
				out << R(k) << "\t" << Z(k) << "\t" << chk*periode << "\t" << psi(k) << "\t" << theta(k) << "\t" << r(k) << endl;
			}
		}

		if(idx_all >= 2*mrat)
		{
			ofs2 << "All possible periodic points have been found." << endl;
			cout << "All possible periodic points have been found." << endl;
			ofs2 << "Program terminated normally" << endl;
			cout << "Program terminated normally" << endl;

			#ifdef m3dc1
			if(PAR.response_field >= 0) M3D.unload();
			#endif

			return 0;
		}
	}
	ofs2 << i+1 << "\t" << flush;
}
ofs2 << endl;
ofs2 << "End of search grid. Not all periodic points have been found." << endl;
cout << "End of search grid. Not all periodic points have been found." << endl;
ofs2 << "Program terminated normally" << endl;
cout << "Program terminated normally" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------- gcd -------------------------------------------------------
// find greatest common divisor
int gcd(int a, int b)
{
   if (b == 0) return a;
   return gcd(b, a % b);

   /*
   if (a == 0 || b == 0)
      return 0;
   else if (a == b)
      return a;
   else if (a > b)
      return gcd(a-b, b);
   else return gcd(a, b-a);
   */
}

//----------- inside ----------------------------------------------------
// checks if current FLT point is inside the search grid
bool inside(PARTICLE& FLT, IO& PAR)
{
double r,th;
switch(PAR.create_flag)
{
case 0: case 1:
	r = FLT.get_r();
	th = FLT.get_theta();
	if(r >= PAR.rmin && r <= PAR.rmax && th >= PAR.thmin && th <= PAR.thmax) return true;
	else return false;
case 3: case 4:
	th = FLT.get_theta();
	if(FLT.psi >= PAR.rmin && FLT.psi <= PAR.rmax && th >= PAR.thmin && th <= PAR.thmax) return true;
	else return false;
case 5:
	if(FLT.R >= PAR.Rmin && FLT.R <= PAR.Rmax && FLT.Z >= PAR.Zmin && FLT.Z <= PAR.Zmax) return true;
	else return false;
}
return false;
}

//----------- newton2dim ------------------------------------------------
// all variables with 'r' represent R and all with 't' represent Z!!!
// J[0] = dR(i+1)/dR(i),  J[1] = dR(i+1)/dZ(i),  J[2] = dZ(i+1)/dR(i),  J[3] = dZ(i+1)/dZ(i)
int newton2D(PARTICLE& FLT, double phistart, int periode)	//0: ok		-1: Error
{
double fr,ft,dr,dt,det,length;
int i,chk;

const int imax = 100;
const double delta = 1e-12;

// Vectors
double J[4];

// Search
double R = FLT.R;
double Z = FLT.Z;

for(i=0;i<=imax;i++)
{
	FLT.phi = phistart;
	chk = mapit_J(FLT,J,periode);
	if(chk<0){ofs2 << "No convergence " << chk << endl; return -1;}
	
	fr = FLT.R - R;
	ft = FLT.Z - Z;

	det = (J[0]-1)*(J[3]-1) - J[1]*J[2];
	dr = ((J[3]-1)*fr-J[1]*ft)/det;
	dt = ((J[0]-1)*ft-J[2]*fr)/det;

	length = sqrt(dr*dr + dt*dt);
	if(length<delta) 
	{
		return 0;	// convergence
	} 

	R -= dr;
	Z -= dt;
	FLT.R = R;
	FLT.Z = Z;
}

ofs2 << "No convergence " <<  R << "\t" << Z << "\t" << dr << "\t" << dt << "\t" << length << endl;
return -1;
}

//--------- mapit_J -----------------------------------------------------------------
int mapit_J(PARTICLE& FLT, double J[], int itt, int Map)
{
int chk;
double Ralt,Zalt,phialt;
const double dr = 0.00001;

double R_stencil[4],Z_stencil[4];	// 0:R+dr 1:R-dr 2:Z+dr 3:Z-dr

Ralt = FLT.R;
Zalt = FLT.Z;
phialt = FLT.phi;

//right: R+dr
FLT.R = Ralt + dr;
FLT.Z = Zalt;
FLT.phi = phialt;

chk = FLT.mapit(itt,Map);
if(chk<0) return -1;
R_stencil[0] = FLT.R;
Z_stencil[0] = FLT.Z;

//left: R-dr
FLT.R = Ralt - dr;
FLT.Z = Zalt;
FLT.phi = phialt;

chk = FLT.mapit(itt,Map);
if(chk<0) return -1;
R_stencil[1] = FLT.R;
Z_stencil[1] = FLT.Z;

//up: Z+dr
FLT.R = Ralt;
FLT.Z = Zalt + dr;
FLT.phi = phialt;

chk = FLT.mapit(itt,Map);
if(chk<0) return -1;
R_stencil[2] = FLT.R;
Z_stencil[2] = FLT.Z;

//down: Z-dr
FLT.R = Ralt;
FLT.Z = Zalt - dr;
FLT.phi = phialt;

chk = FLT.mapit(itt,Map);
if(chk<0) return -1;
R_stencil[3] = FLT.R;
Z_stencil[3] = FLT.Z;

// derivatives
// J[0]=dR(i+1)/dR(i) J[1]=dR(i+1)/dZ(i) J[2]=dZ(i+1)/dR(i) J[3]=dZ(i+1)/dZ(i)
J[0] = 0.5*(R_stencil[0]-R_stencil[1])/dr;
J[1] = 0.5*(R_stencil[2]-R_stencil[3])/dr;
J[2] = 0.5*(Z_stencil[0]-Z_stencil[1])/dr;
J[3] = 0.5*(Z_stencil[2]-Z_stencil[3])/dr;

// center
FLT.R = Ralt;
FLT.Z = Zalt;
FLT.phi = phialt;

chk = FLT.mapit(itt,Map);
if(chk<0) return -1;

return 0;
}

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



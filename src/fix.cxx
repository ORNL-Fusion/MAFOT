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
#else
	#define program_name "dtfix"
#endif

// Include
//--------
#include <mafot.hxx>
#if defined(ITER)
	#if defined(m3dc1)
		#include <iter_m3dc1.hxx>
	#else
		#include <iter.hxx>
	#endif
#elif defined(NSTX)
	#if defined(m3dc1)
		#include <nstx_m3dc1.hxx>
	#else
		#include <nstx.hxx>
	#endif
#elif defined(MAST)
	#if defined(m3dc1)
		#include <mast_m3dc1.hxx>
	#else
		#include <mast.hxx>
	#endif
#else
	#if defined(m3dc1)
		#include <d3d_m3dc1.hxx>
	#else
		#include <d3d.hxx>
	#endif
#endif

// Prototypes  
//-----------
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

// Input filenames
LA_STRING basename;
LA_STRING praefix = "";
if(argc==4) praefix = "_" + LA_STRING(argv[3]);
if(argc>=3) 
{
	periode = atoi(argv[2]);
	basename = LA_STRING(argv[1]); 
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
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prep_perturbation(EQD,PAR);

// Prepare particles
PARTICLE FLT(EQD,PAR);

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

ofs2 << Nx << " rows, done:" << endl;
for(i=0;i<Nx;i++)	//r or R
{
	x = xmin + i*dx;
	for(j=0;j<Ny;j++)	//theta or Z
	{
		y = ymin + j*dy;
		switch(PAR.create_flag)
		{
		case 0:
			FLT.convertRZ(y,x);
			break;
		case 5:
			FLT.R = x;
			FLT.Z = y;
			break;
		}

		chk = newton2D(FLT,PAR.phistart,periode);
		if(chk==-1) continue;

		if(periode==1)
		{
			if(inside(FLT,PAR))
			{
				out << FLT.R << "\t" << FLT.Z << "\t" << periode << "\t" << FLT.psi << "\t" << FLT.get_theta() << "\t" << FLT.get_r() << endl;
				cout << "Program terminated normally" << endl;
				ofs2 << "Program terminated normally" << endl;
				#ifdef m3dc1
				if(PAR.response_field >= 0) M3D.unload();
				#endif
				return 0; 
			}
			else continue;
		}
		else out << FLT.R << "\t" << FLT.Z << "\t" << periode << "\t" << FLT.psi << "\t" << FLT.get_theta() << "\t" << FLT.get_r() << endl;
	}
	ofs2 << i+1 << "\t" << flush;
}
ofs2 << endl;
ofs2 << "Program terminated normally" << endl;
cout << "Program terminated normally" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------- inside ----------------------------------------------------
// checks if current FLT point is inside the search grid
bool inside(PARTICLE& FLT, IO& PAR)
{
double r,th;
switch(PAR.create_flag)
{
case 0:
	r = FLT.get_r();
	th = FLT.get_theta();
	if(r >= PAR.rmin && r <= PAR.rmax && th >= PAR.thmin && th <= PAR.thmax) return true;
	else return false;
case 5:
	if(FLT.R >= PAR.Rmin && FLT.R <= PAR.Rmax && FLT.Z >= PAR.Zmin && FLT.Z <= PAR.Zmax) return true;
	else return false;
}
}

//----------- newton2dim ------------------------------------------------
// all variables with 'r' represent R and all with 't' represent Z!!!
// J(1) = dR(i+1)/dR(i),  J(2) = dR(i+1)/dZ(i),  J(3) = dZ(i+1)/dR(i),  J(4) = dZ(i+1)/dZ(i)
int newton2D(PARTICLE& FLT, double phistart, int periode)	//0: ok		-1: Fehler
{
double fr,ft,dr,dt,det,length;
int i,chk;

const int imax = 100;
const double delta = 1e-12;

// Vectors
//Array<double,1> J(Range(1,4));
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
	//if(i%20==0){cout << dr << "\t" << dt << "\t" << length <<  endl;	getchar();}
	if(length<delta) 
	{
		return 0;	// convergency
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

// no Array and it dies here
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



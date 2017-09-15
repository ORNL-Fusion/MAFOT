// Program calculates path of field lines in 10 deg steps. Used for 3-D Visualization
// for particle-drift with Time dependent perturbations
// Fortran subroutines are used for perturbations
// A.Wingen						20.06.11

// Input: 1: Parameterfile	2: Initial value file (optional; if not, enter x instead) 3: praefix (optional)
// Output:	paths of individual field lines in 10 deg steps between the target plates
//			log-file


// Define
//--------
#if defined(ITER)
	#define program_name "iterstructure"
#elif defined(NSTX)
	#define program_name "nstxstructure"
#elif defined(MAST)
	#define program_name "maststructure"
#else
	#define program_name "dtstructure"
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

// namespaces
//-----------
using namespace blitz;

// Prototypes 
//-----------
inline double modulo(double x, double y);

// Switches
//---------

// Golbal Parameters 
//------------------

// Main Program
//--------------
int main(int argc, char *argv[])
{
// Variables
int i,j,k,chk;
int usePointfile;
//double dummy;
EFIT EQD;

// Use system time as seed(=idum) for random numbers
double now = zeit();
//long idum = long(now);

// Input file names
LA_STRING basename,pointname;
LA_STRING praefix = "";
if(argc==4) praefix = "_" + LA_STRING(argv[3]);
if(argc>=3) {basename = LA_STRING(argv[1]); pointname = LA_STRING(argv[2]);}
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

if(pointname == 'x') usePointfile = 0;
else
{
	usePointfile = 1;
	pointname = checkparfilename(pointname);
	pointname = pointname + ".dat";
}

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + ".dat");
ofs2.precision(16);

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
int nstep = 10;					// Number of dpinit steps
double dphi = nstep*dpinit;		// in deg;  dpinit = 1.0 (default)
double alpha = PAR.verschieb;	// Scales parabolic deformation of line between start and end point, alpha = 0: no deformation, alpha = 1: max deformation equals distance between points
bool angleInDeg = true;			// phi angle in output file is in degrees (left-handed machine angle), else in radiants (right-handed angle)

// additional parameters for IO
PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
PAR.pv[1].name = "phistart";		PAR.pv[1].wert = PAR.phistart;
PAR.pv[2].name = "MapDirection";	PAR.pv[2].wert = PAR.MapDirection;
PAR.pv[3].name = "energy ratio lambda";	PAR.pv[3].wert = PAR.lambda;
PAR.pv[4].name = "Ekin";			PAR.pv[4].wert = PAR.Ekin;
PAR.pv[5].name = "Number of Points";PAR.pv[5].wert = PAR.N;
PAR.pv[6].name = "Rmin";			PAR.pv[6].wert = PAR.Rmin;
PAR.pv[7].name = "Zmin";			PAR.pv[7].wert = PAR.Zmin;
PAR.pv[8].name = "Rmax";			PAR.pv[8].wert = PAR.Rmax;
PAR.pv[9].name = "Zmax";			PAR.pv[9].wert = PAR.Zmax;

// Read pointfile 
Array<double,2> initial;
double dR,dZ,dN;
int inputPoints_columns;
if(usePointfile == 1) 
{
	// points in file are:    R[m]		phi[deg] (left-handed machine angle)	Z[m]
	inputPoints_columns = count_column(pointname);
	readfile(pointname,inputPoints_columns,initial);
	PAR.N = initial.rows();
	PAR.pv[5].wert = PAR.N;
}
else	// or construct initial points from straight line between (Rmin,Zmin) and (Rmax,Zmax) at constant phi angle
{
	initial.resize(Range(1,PAR.N),Range(1,2));	// secondIndex = 1: R	secondIndex = 2: Z

	if(PAR.N>1) dN = 1.0/double(PAR.N-1);
	else dN = 0;
	dR = (PAR.Rmax-PAR.Rmin)*dN;
	dZ = (PAR.Zmax-PAR.Zmin)*dN;

	for(i=1;i<=PAR.N;i++) 
	{
		initial(i,1) = PAR.Rmin + (i-1)*(dR - 4*alpha*((i-1)*dN-1)*dZ);
		initial(i,2) = -PAR.phistart;	// phistart is degrees, but right-handed, so initial is left-handed machine angle
		initial(i,3) = PAR.Zmin + (i-1)*(dZ + 4*alpha*((i-1)*dN-1)*dR);
	}
}
//cout << initial << endl;

// Prepare Perturbation
prep_perturbation(EQD,PAR);

// Prepare particles
PARTICLE FLT(EQD,PAR);

// Output
LA_STRING filenameout = "struct" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "X[m]";  var[1] = "Y[m]";  var[2] = "Z[m]";  var[3] = "R[m]";  var[4] = "phi[deg]";
if(not angleInDeg) var[4] = "phi[rad]";
PAR.writeiodata(out,bndy,var);

cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
cout << "Start Tracer for " << PAR.N << " points ... " << endl;
ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
ofs2 << "Start Tracer for " << PAR.N << " points ... " << endl;

// Follow the field lines
int size = PAR.itt*int(360.0/double(dphi));
Array<double,2> data(Range(1,size),Range(1,3));
for(i=1;i<=PAR.N;i++)
{
	// Set initial conditions
	FLT.R = initial(i,1);
	PAR.phistart = -initial(i,2);
	FLT.Z = initial(i,3);
	FLT.phi = PAR.phistart;
	FLT.get_psi(FLT.R,FLT.Z,FLT.psi);
	if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}

	// negative direction
	if(PAR.MapDirection <= 0)
	{
		for(j=1;j<=size;j++)
		{
			chk = FLT.mapstep(-1,nstep);
			if(chk==-1) break;
			if(fabs(FLT.phi + j*dpinit*dphi - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi + j*dpinit*dphi - PAR.phistart) << endl;
			FLT.phi = -j*dpinit*dphi + PAR.phistart;

			//Store Values
			data(j,1) = FLT.R;	data(j,2) = FLT.Z;	data(j,3) = FLT.phi;
		}
		// Write stored values in reverse direction
		if(angleInDeg) {for(k=j-1;k>=1;k--) out << data(k,1)*cos(data(k,3)/rTOd) << "\t" << data(k,1)*sin(data(k,3)/rTOd) << "\t" << data(k,2) << "\t" << data(k,1) << "\t" << modulo(360.0 - data(k,3), 360.0) << endl;}
		else {for(k=j-1;k>=1;k--) out << data(k,1)*cos(data(k,3)/rTOd) << "\t" << data(k,1)*sin(data(k,3)/rTOd) << "\t" << data(k,2) << "\t" << data(k,1) << "\t" << data(k,3)/rTOd << endl;}
	}

	// Restore start values and write them
	FLT.R = initial(i,1);
	PAR.phistart = -initial(i,2);
	FLT.Z = initial(i,3);
	FLT.phi = PAR.phistart;
	FLT.get_psi(FLT.R,FLT.Z,FLT.psi);
	if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}
	if(angleInDeg) {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << modulo(360.0 - FLT.phi, 360.0) << endl;}
	else {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << FLT.phi/rTOd << endl;}

	//positive direction
	if(PAR.MapDirection >= 0)
	{
		for(j=1;j<=size;j++)
		{
			chk = FLT.mapstep(1,nstep);
			if(chk==-1) break;
			if(fabs(FLT.phi - j*dpinit*dphi - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi - j*dpinit*dphi - PAR.phistart) << endl;
			FLT.phi = j*dpinit*dphi + PAR.phistart;

			if(angleInDeg) {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << modulo(360.0 - FLT.phi, 360.0) << endl;}
			else {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << FLT.phi/rTOd << endl;}
		}
	}
}

// Copy results to filament.in file
// interactive
LA_STRING input;
int FileNr;
double Current;

// Close previous output file and clear ofstream
out.close();
out.clear();

// Get interactive input
cout << "Save as filament.in file? (y;n): "; cin >> input;
if(input[1] == 'n') 
{
	double now2 = zeit();
	cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
	ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;
	return 0;
}

cout << "Enter File Number: "; cin >> FileNr;
cout << "Enter Current[A]: "; cin >> Current;

// Output
LA_STRING filamentname = "filament" + LA_STRING(FileNr) + ".in";
out.open(filamentname);
out.precision(16);

// Write Header
out << "# Copy of file: " << filenameout << endl;
out << "#-------------------------------------------------" << endl;
out << "### Current[A]: " << Current << endl;
out << "#-------------------------------------------------" << endl;
out << "### Data:" << endl;
out << "# ";
for(i=0;i<int(var.size());i++) out << var[i] << "     ";
out << endl;
out << "#" << endl;

// Read Data
data.free();
readfile(filenameout,5,data);

// Write Data
if(angleInDeg) {for(j=1;j<=data.rows();j++) out << data(j,1) << "\t" << data(j,2) << "\t" << data(j,3) << "\t" << data(j,4) << "\t" << (360 - data(j,5))/rTOd << endl;}	// read in angle is in lhs degree, but filament.in wants rhs radiants
else {for(j=1;j<=data.rows();j++) out << data(j,1) << "\t" << data(j,2) << "\t" << data(j,3) << "\t" << data(j,4) << "\t" << data(j,5) << endl;}

double now2 = zeit();
cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

return 0;
} //end of main
 

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

inline double modulo(double x, double y)
{
double z = fmod(x,y);
if(z < 0) z += y;
return z;
}

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



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
#elif defined(CMOD)
	#define program_name "cmodstructure"
#elif defined(TCABR)
	#define program_name "tcabrstructure"
#elif defined(HEAT)
	#define program_name "heatstructure"
#elif defined(ANYM)
	#define program_name "anystructure"
#else
	#define program_name "dtstructure"
#endif

// Include
//--------
#include <mafot.hxx>
#include <unistd.h>

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
int i,j,k,chk,chk2,skip_connect;
//double dummy;
EFIT EQD;

// Use system time as seed(=idum) for random numbers
double now = zeit();
//long idum = long(now);

// defaults
LA_STRING pointname;
bool filaments = false;
int nstep = 10;					// Number of dpinit steps
bool angleInDeg = false;			// phi angle in output file is in degrees (left-handed machine angle), else in radiants (right-handed angle)
bool usePointfile = false;
LA_STRING woutfile = "wout.nc";
LA_STRING xpandfile = "None";
LA_STRING siestafile = "siesta.dat";
LA_STRING islandfile = "fakeIslands.in";
bool use_3Dwall = false;
LA_STRING wall_file = "none";
bool checkFistStep = false;
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
while ((c = getopt(argc, argv, "habfd:P:X:V:S:I:W:c")) != -1)
switch (c)
{
case 'h':
	cout << "usage: dtstructure [-h] [-a] [-b] [-c] [-d step] [-f] [-I island] [-P points] [-S siesta] [-V wout] [-W wall] [-X xpand] file [tag]" << endl << endl;
	cout << "Trace field lines in 3D and return the path every step d in toroidal angle." << endl << endl;
	cout << "positional arguments:" << endl;
	cout << "  file          Contol file (starts with '_')" << endl;
	cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
	cout << endl << "optional arguments:" << endl;
	cout << "  -h            show this help message and exit" << endl;
	cout << "  -a            output angle in DIII-D angle (left-handed) in degrees, default = radiants and right-handed," << endl;
	cout << "  -b            use simple boundary instead of g-file wall as limiter" << endl;
	cout << "  -c            check if field line crosses wall during first step; only with 3D wall, default = No" << endl;
	cout << "  -d            step size for output, default = 10 degrees" << endl;
	cout << "  -f            create filament.in file, default = No" << endl;
	cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
	cout << "  -P            use separate input file for initial conditions; argument is the file name; default is None" << endl;
	cout << "                File format of columns: R [m], phi [deg, right-handed coord.], Z [m]" << endl;
	cout << "                Header lines start with '#'; no comment lines between/after data possible" << endl;
	cout << "  -S            filename for SIESTA; default, see below" << endl;
	cout << "  -V            filename for VMEC; default, see below" << endl;
	cout << "  -W            use separate 3D Wall-File; default is 2D wall from EFIT file" << endl;
	cout << "  -X            filename for XPAND; default is None" << endl;
	cout << endl << "Examples:" << endl;
	cout << "  dtstructure _struct.dat blabla" << endl;
	cout << "  dtstructure _struct.dat testrun -P startpoints.dat -d 2 -a" << endl;
	cout << endl << "Infos:" << endl;
	cout << "  To use B-field from M3DC1, set response_field >= 0, and provide file in cwd:" << endl;
	cout << "    m3dc1sup.in    ->  location and scale factor for M3DC1 output C1.h5" << endl;
	cout << "  To use B-field from XPAND, set response_field = -3, and provide files in cwd:" << endl;
	cout << "    xpand.dat      ->  B-field on 3D grid from XPAND; use option -X to specify a filename (default is None -> inside VMEC only)" << endl;
	cout << "    wout.nc        ->  VMEC output; use option -V to specify other filename" << endl;
	cout << "  To use B-field from VMEC, inside only (no xpand file given), set response_field = -3, and provide file in cwd:" << endl;
	cout << "    wout.nc        ->  VMEC output; use option -V to specify other filename" << endl;
	cout << "  To use B-field from SIESTA, set response_field = -2, and provide file in cwd:" << endl;
	cout << "    siesta.dat     ->  B-field on 3D grid; use option -S to specify other filename" << endl;
	cout << "  To use B-field for mock-up islands, set response_field = -10, and provide file in cwd:" << endl;
	cout << "    fakeIslands.in ->  each line gives: Amplitude, pol. mode m, tor. mode n, phase [rad]" << endl;
	cout << "                       use option -I to specify other filename" << endl;
	cout << endl << "Current MAFOT version is: " << MAFOT_VERSION << endl;
	return 0;
case 'a':
	angleInDeg = true;
	break;
case 'b':
	simpleBndy = 1;
	break;
case 'c':
	checkFistStep = true;
	break;
case 'd':
	nstep = atoi(optarg);
	break;
case 'f':
	filaments = true;
	break;
case 'P':
	usePointfile = true;
	pointname = optarg;
	break;
case 'I':
	islandfile = optarg;
	break;
case 'S':
	siestafile = optarg;
	break;
case 'V':
	woutfile = optarg;
	break;
case 'W':
	use_3Dwall = true;
	wall_file = optarg;
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
	EXIT;
default:
	EXIT;
}

// Input file names
LA_STRING basename;
LA_STRING praefix = "";
if(argc==optind+2) praefix = "_" + LA_STRING(argv[optind+1]);
if(argc>=optind+1) basename = LA_STRING(argv[optind]);
else {cout << "No Input files -> Abort!" << endl; EXIT;}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + ".dat");
ofs2.precision(16);

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

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

EQD.ReadData(EQD.gFile,Raxis,Zaxis);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.gFile << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.gFile << endl;

if(simpleBndy == 1)	// set boundary to the EFIT limit.
{
	bndy[0] = min(EQD.R);
	bndy[1] = max(EQD.R);
	bndy[2] = min(EQD.Z);
	bndy[3] = max(EQD.Z);
	cout << "Bounding box: Rmin = " << bndy[0] << "  Rmax = " << bndy[1] << "  Zmin = " << bndy[2] << "  Zmax = " << bndy[3] << endl;
	ofs2 << "Bounding box: Rmin = " << bndy[0] << "  Rmax = " << bndy[1] << "  Zmin = " << bndy[2] << "  Zmax = " << bndy[3] << endl;
}

// Read 3D wall file and add to EQD
if(use_3Dwall)
{
	cout << "Using 3D wall from file: " << wall_file << endl;
	ofs2 << "Using 3D wall from file: " << wall_file << endl;
	EQD.set3Dwall(wall_file);
}
else checkFistStep = false;

// Set starting parameters
double dphi = nstep*dpinit;		// in deg;  dpinit = 1.0 (default)
double alpha = PAR.verschieb;	// Scales parabolic deformation of line between start and end point, alpha = 0: no deformation, alpha = 1: max deformation equals distance between points

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

PAR.output_step_size = nstep;

// Read pointfile
Array<double,2> initial;
double dR,dZ,dN;
int inputPoints_columns;
if(usePointfile)
{
	// points in file are:    R[m]		phi[deg] (right-handed angle)	Z[m]
	inputPoints_columns = count_column(pointname);
	if(inputPoints_columns < 3) {cout << "Less than 3 columns of data in points file. -> Abort!" << endl; EXIT;}
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
		initial(i,2) = PAR.phistart;	// phistart is degrees, right-handed
		initial(i,3) = PAR.Zmin + (i-1)*(dZ + 4*alpha*((i-1)*dN-1)*dR);
	}
}
//cout << initial << endl;

// Prepare Perturbation
prepare_common_perturbations(EQD,PAR,0,siestafile,xpandfile,islandfile);
prep_perturbation(EQD,PAR);

// Prepare collisions
COLLISION COL;
if (use_collision) COL.init(TprofileFile, NprofileFile, f, zbar, PAR.Zq, PAR.Mass, rc, mc);
// Prepare particles
PARTICLE FLT(EQD,PAR,COL);

// Output
LA_STRING filenameout = "struct" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(8);
var[0] = "X[m]";  var[1] = "Y[m]";  var[2] = "Z[m]";  var[3] = "R[m]";  var[4] = "phi[deg]";  var[5] = "psi";  var[6] = "Lc[m]";  var[7] = "dpsidLc[1/m]";
if(not angleInDeg) var[4] = "phi[rad]";
PAR.writeiodata(out,bndy,var);

cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
cout << "Start Tracer for " << PAR.N << " points ... " << endl;
ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
ofs2 << "Start Tracer for " << PAR.N << " points ... " << endl;

// Follow the field lines
#ifdef HEAT
//Modified by TL 20191117
//Using 360 in the line below constrains the minimum field trace to be 360deg
//Now, we just make the minimum itt
int size = PAR.itt;
#else
int size = PAR.itt*int(360.0/double(dphi));
#endif
Array<double,2> data(Range(1,size),Range(1,6));
for(i=1;i<=PAR.N;i++)
{
	// Set initial conditions
	FLT.R = initial(i,1);
	PAR.phistart = initial(i,2);
	FLT.Z = initial(i,3);
	FLT.phi = PAR.phistart;
	FLT.Lc = 0;
	FLT.psimin = 10;
	FLT.steps = 0;
	FLT.dpsidLc = 0;
	FLT.get_psi(FLT.R,FLT.Z,FLT.psi);
	if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}


	// check if particle/field line left the system during the first integration step
	// this does not update/change any FLT parameters
	skip_connect = 0;
	if(checkFistStep)
	{
		if(PAR.MapDirection == 0)
		{
			chk2 = FLT.mapstep(1, 1, true);		// one step in forward MapDirection
			chk = FLT.mapstep(-1, 1, true);		// one step in backward MapDirection
			if(chk + chk2 == -2) skip_connect = 1;
		}
		else
		{
			chk = FLT.mapstep(PAR.MapDirection, 1, true);	// one step in backward MapDirection
			if(chk == -1) skip_connect = 1;
		}
		if ((fabs(initial(i,1) - FLT.R) > 0) || (fabs(initial(i,3) - FLT.Z) > 0) || (fabs(initial(i,2) - FLT.phi) > 0))
			cout << initial(i,1) - FLT.R << "\t" << initial(i,3) - FLT.Z << "\t" << initial(i,2) - FLT.phi << endl;
		cout << "Field line exits during first step (0 = No, 1 = Yes): " << skip_connect << endl;
	}


	// negative direction
	if(PAR.MapDirection <= 0)
	{
		for(j=1;j<=size;j++)
		{
			chk = FLT.mapstep(-1,nstep);
			if(fabs(FLT.phi + j*dpinit*dphi - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi + j*dpinit*dphi - PAR.phistart) << endl;
			FLT.phi = -j*dpinit*dphi + PAR.phistart;

			//Store Values
			data(j,1) = FLT.R;	data(j,2) = FLT.Z;	data(j,3) = FLT.phi;	data(j,4) = FLT.psi;	data(j,5) = FLT.Lc;	data(j,6) = FLT.dpsidLc;
			if(chk==-1) break;
		}
		// Write stored values in reverse direction
		if(j > size) j = size;
		if(angleInDeg) {for(k=j;k>=1;k--) out << data(k,1)*cos(data(k,3)/rTOd) << "\t" << data(k,1)*sin(data(k,3)/rTOd) << "\t" << data(k,2) << "\t" << data(k,1) << "\t" << modulo(360.0 - data(k,3), 360.0) << "\t" << data(k,4) << "\t" << FLT.Lc-data(k,5) << "\t" << data(k,6) << endl;}
		else {for(k=j;k>=1;k--) out << data(k,1)*cos(data(k,3)/rTOd) << "\t" << data(k,1)*sin(data(k,3)/rTOd) << "\t" << data(k,2) << "\t" << data(k,1) << "\t" << data(k,3)/rTOd << "\t" << data(k,4) << "\t" << FLT.Lc-data(k,5) << "\t" << data(k,6) << endl;}
		FLT.dpsidLc = data(1,6);
		FLT.Lcmin = FLT.Lc - FLT.Lcmin;
	}

	// Restore start values and write them, keep Lc and others though
	FLT.R = initial(i,1);
	PAR.phistart = initial(i,2);
	FLT.Z = initial(i,3);
	FLT.phi = PAR.phistart;
	FLT.get_psi(FLT.R,FLT.Z,FLT.psi);
	if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}
	if(angleInDeg) {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << modulo(360.0 - FLT.phi, 360.0) << "\t" << FLT.psi << "\t" << FLT.Lc << "\t" << FLT.dpsidLc << endl;}
	else {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << FLT.phi/rTOd << "\t" << FLT.psi << "\t" << FLT.Lc << "\t" << FLT.dpsidLc  << endl;}

	//positive direction
	if(PAR.MapDirection >= 0)
	{
		for(j=1;j<=size;j++)
		{
			chk = FLT.mapstep(1,nstep);
			if(fabs(FLT.phi - j*dpinit*dphi - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi - j*dpinit*dphi - PAR.phistart) << endl;
			FLT.phi = j*dpinit*dphi + PAR.phistart;

			if(angleInDeg) {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << modulo(360.0 - FLT.phi, 360.0) << "\t" << FLT.psi << "\t" << FLT.Lc << "\t" << FLT.dpsidLc << endl;}
			else {out << FLT.R*cos(FLT.phi/rTOd) << "\t" << FLT.R*sin(FLT.phi/rTOd) << "\t" << FLT.Z << "\t" << FLT.R << "\t" << FLT.phi/rTOd << "\t" << FLT.psi << "\t" << FLT.Lc << "\t" << FLT.dpsidLc << endl;}
			if(chk==-1) break;
		}
	}

	// Write field line averaged statistics
	//commented by TL (prints too much to HEAT)
	//cout << "psimin = " << FLT.psimin << "\t" << "Lc at psimin = " << FLT.Lcmin << "\t" << "<dpsi/dLc> = " << FLT.get_dpsidLc_average() << "\t" << "Lc = " << FLT.Lc << "\t" << "steps = " << FLT.steps << endl;
	ofs2 << "psimin = " << FLT.psimin << "\t" << "Lc at psimin = " << FLT.Lcmin << "\t" << "<dpsi/dLc> = " << FLT.get_dpsidLc_average() << "\t" << "Lc = " << FLT.Lc << "\t" << "steps = " << FLT.steps << endl;
}

// Close previous output file and clear ofstream
out.close();
out.clear();

// If filament.in file is not requested, end here!
if(not filaments)
{
	double now2 = zeit();
	cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
	ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;
	#ifdef m3dc1
	if(PAR.response_field >= 0) M3D.unload();
	#endif
	return 0;
}

// Copy results to filament.in file
// interactive
int FileNr;
double Current;
double phistep = nstep*dpinit;
if(not angleInDeg) phistep /= rTOd;
Array<double,1> jstart(PAR.N),jend(PAR.N);

// Read Data
data.free();
readfile(filenameout,8,data);

// find individual flux tubes
jstart(0) = 1;
jend(PAR.N-1) = data.rows();
i = 1;
j = 2;
while((i < PAR.N) && (j < data.rows()))
{
	if(fabs(data(j,5) - data(j-1,5)) > 2*phistep)
	{
		jstart(i)  = j;
		jend(i-1) = j-1;
		i += 1;
	}
	j += 1;
}

for(i=0;i<PAR.N;i++)
{
	// Get interactive input
	cout << endl << "Creating filament.in file" << endl;
	cout << "Enter File Number: "; cin >> FileNr;
	cout << "Enter Current[A]: "; cin >> Current;

	// Output
	LA_STRING filamentname = "filament" + LA_STRING(FileNr) + ".in";
	out.open(filamentname);
	out.precision(16);

	// Write Header
	out << "# Single flux tube from file: " << filenameout << endl;
	out << "#-------------------------------------------------" << endl;
	out << "### Current[A]: " << Current << endl;
	out << "#-------------------------------------------------" << endl;
	out << "### Data:" << endl;
	out << "# ";
	var[4] = "phi[rad]";
	for(int k=0;k<int(var.size());k++) out << var[k] << "     ";
	out << endl;
	out << "#" << endl;

	// Write Data
	if(angleInDeg) {for(j=jstart(i);j<=jend(i);j++) out << data(j,1) << "\t" << data(j,2) << "\t" << data(j,3) << "\t" << data(j,4) << "\t" << (360 - data(j,5))/rTOd << endl;}	// read in angle is in lhs degree, but filament.in wants rhs radiants
	else {for(j=jstart(i);j<=jend(i);j++) out << data(j,1) << "\t" << data(j,2) << "\t" << data(j,3) << "\t" << data(j,4) << "\t" << data(j,5) << endl;}
	out.close();
}
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

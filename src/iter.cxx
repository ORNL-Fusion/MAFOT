// Program iterates the first fixed point from a file forward using the map machinery
// Takes a fixed point file (output from dtfix) and iterates forward N times
// A.Wingen						17.01.26

// Input: 1: Parameterfile	2: Fixed point file	3: Number of iterations	4: praefix(optional)
// Output:	iterated points
//			log-file

// Define
//--------
#if defined(ITER)
	#define program_name "iteriter"
#elif defined(NSTX)
	#define program_name "nstxiter"
#elif defined(MAST)
	#define program_name "mastiter"
#elif defined(CMOD)
	#define program_name "cmoditer"
#elif defined(TCABR)
	#define program_name "tcabrfiter"
#elif defined(ANYM)
	#define program_name "anyiter"
#else
	#define program_name "dtiter"
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

// Switches
//----------

// Golbal Parameters
//------------------

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// Variables
int i,chk;
int iterations = 1;
int periode = 1;
double Rstart, Zstart;
EFIT EQD;

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
	cout << "usage: dtiter [-h] [-I island] [-S siesta] [-V wout] [-X xpand] file fix_file iterations [tag]" << endl << endl;
	cout << "Iterate a fixed point forward using the map machinery." << endl << endl;
	cout << "positional arguments:" << endl;
	cout << "  file          Control file (starts with '_')" << endl;
	cout << "  fix_file      Fixed point file (output from dtfix)" << endl;
	cout << "  iterations    Number of times to apply the map" << endl;
	cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
	cout << endl << "optional arguments:" << endl;
	cout << "  -h            show this help message and exit" << endl;
	cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
	cout << "  -S            filename for SIESTA; default, see below" << endl;
	cout << "  -V            filename for VMEC; default, see below" << endl;
	cout << "  -X            filename for XPAND; default, see below" << endl;
	cout << endl << "Examples:" << endl;
	cout << "  dtiter _iter.dat fix_1.dat 5 test" << endl;
	cout << endl << "Infos:" << endl;
	cout << "  To use B-field from M3DC1, set response_field >= 0, and provide file in cwd:" << endl;
	cout << "    m3dc1sup.in    ->  location and scale factor for M3DC1 output C1.h5" << endl;
	cout << "  To use B-field from XPAND, set response_field = -3, and provide files in cwd:" << endl;
	cout << "    xpand.dat      ->  B-field on 3D grid from XPAND; use option -X to specify a filename" << endl;
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
LA_STRING basename,fixname;
LA_STRING praefix = "";
if(argc==optind+4) praefix = "_" + LA_STRING(argv[optind+3]);
if(argc>=optind+3)
{
	basename = LA_STRING(argv[optind]);
	fixname = LA_STRING(argv[optind+1]);
	iterations = atoi(argv[optind+2]);
}
if(argc<=optind+2)
{
	cout << "Not enough input arguments -> Abort!" << endl; 
	exit(0);
}

basename = checkparfilename(basename);
fixname = checkparfilename(fixname);
LA_STRING parfilename = "_" + basename + ".dat";
LA_STRING name = fixname + ".dat";

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// read fixed points
Array<double,2> data;
readfile(name,6,data);

if(data.rows() < 1)
{
	cout << "No fixed points in file -> Abort!" << endl;
	exit(0);
}

// Use only the first fixed point
Rstart = data(1,1);
Zstart = data(1,2);
periode = int(fabs(data(1,3)));

// Output
ofstream ofs2,out;
ofs2.precision(16);
out.precision(16);

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + ".dat");
ofs2 << "Read Parameterfile " << parfilename << endl;

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

// Prepare Perturbation
prepare_common_perturbations(EQD,PAR,0,siestafile,xpandfile,islandfile);
prep_perturbation(EQD,PAR);

// Prepare collisions
COLLISION COL;
if (use_collision) COL.init(TprofileFile, NprofileFile, f, zbar, PAR.Zq, PAR.Mass, rc, mc);

// Create particle for iteration
PARTICLE FLT(EQD,PAR,COL);

// Output file
LA_STRING filenameout = "iter_" + LA_STRING(iterations) + "_" + LA_STRING(periode) + praefix + ".dat";
outputtest(filenameout);
out.open(filenameout);

// Write header
out << "# Iteration of fixed point from file " << name << endl;
out << "# Starting point: R = " << Rstart << ", Z = " << Zstart << ", phi = " << PAR.phistart << endl;
out << "# Period = " << periode << ", MapDirection = " << PAR.MapDirection << endl;
out << "# Iterations = " << iterations << endl;
out << "# Columns: iteration, R[m], Z[m], psi, theta[rad], r[m], phi[deg]" << endl;

ofs2 << "Starting point: R = " << Rstart << ", Z = " << Zstart << endl;
ofs2 << "Period = " << periode << ", Iterations = " << iterations << endl;

// Initialize particle at fixed point
FLT.R = Rstart;
FLT.Z = Zstart;
FLT.phi = PAR.phistart;

// Output initial point
out << 0 << "\t" << FLT.R << "\t" << FLT.Z << "\t" << FLT.psi << "\t" 
	<< FLT.get_theta() << "\t" << FLT.get_r() << "\t" << FLT.phi*rTOd << endl;

// Iterate forward
for(i=1; i<=iterations; i++)
{
	// Apply the map
	chk = FLT.mapit(1, PAR.MapDirection);
	
	if(chk < 0)
	{
		ofs2 << "Error during iteration " << i << ": chk = " << chk << endl;
		cout << "Error during iteration " << i << ": chk = " << chk << endl;
		break;
	}
	
	// Output iterated point
	out << i << "\t" << FLT.R << "\t" << FLT.Z << "\t" << FLT.psi << "\t" 
		<< FLT.get_theta() << "\t" << FLT.get_r() << "\t" << FLT.phi*rTOd << endl;
	
	ofs2 << "Iteration " << i << ": R = " << FLT.R << ", Z = " << FLT.Z << ", psi = " << FLT.psi << endl;
}

out.close();
ofs2 << "Program terminates normally" << endl;
cout << "Program terminates normally" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

return 0;
} //end of main

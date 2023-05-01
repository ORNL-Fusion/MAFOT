// Program calculates unstable manifold of a hyp. fixed point for particle-Drift and time dependent perturbation
// For stable manifold use reverse Integration (see MapDirection)
// For left or right hand-sided set sign of 'verschieb' in Parameterfile to - or + respectively
// Fortran subroutines are used for the perturbations
// A.Wingen						16.06.11

// Input: 1: Parameterfile	2: File with fixed points	3: praefix(optional)
//			fixed points have to be in toroidal coordinates, not cylindrical!!!
// Output:	one file for each of the fixed points specified in the input-file, giving the respective manifold
//			log-file

// Define
//--------
#if defined(ITER)
	#define program_name "iterman"
#elif defined(NSTX)
	#define program_name "nstxman"
#elif defined(MAST)
	#define program_name "mastman"
#elif defined(CMOD)
	#define program_name "cmodman"
#elif defined(TCABR)
	#define program_name "tcabrman"
#elif defined(ANYM)
	#define program_name "anyman"
#else
	#define program_name "dtman"
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
void find_start_values(Array<double,1>& xs, Array<double,1>& d, int periode, double& verschieb, 
					   double phistart, int MapDirection, PARTICLE& FLT);
inline double abstand(PARTICLE& FLT1, PARTICLE& FLT2);
inline double abstand(Array<double,1>& x1, Array<double,1>& x2);

// Switches
//----------

// Golbal Parameters
//------------------

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// Variables
int i,k,periode,chk,plotchk,end,variate;
int skipattempt,outside;
double xfix,yfix;
double t,dt,talt,dist;
EFIT EQD;

// Arrays
Array<double,1> xa(Range(1,2));
Array<double,1> xsa(Range(1,2));
Array<double,1> da(Range(1,2));

// defaults
int trytoskip = 1;			// 0: stop after first wall contact		1: try to continue, may cause errors!!!
int skipmax = 3;			// number of skip attempts, if trytoskip == 1 (usually 3)
int preventSmallSteps = 0;	// 0: code stops if step size dt < 1e-14	1: code continues with dt = 1e-10 as long as step size controll would reduce dt below 1e-10
int bndy_type = 0;			// 0: vessel wall   1: box around vessel wall   2: EFIT grid boundary
double minabs = 0.001;		//0.001
double maxabs = 0.005;		//0.005
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
while ((c = getopt(argc, argv, "hsm:b:n:x:X:V:S:I:")) != -1)
switch (c)
{
case 'h':
	cout << "usage: dtman [-h] [-s] [-m max] [-b bndy] [-N MIN] [-X MAX] [-I island] [-S siesta] [-V wout] [-X xpand] file fix_file [tag]" << endl << endl;
	cout << "Calculate stable and unstable manifolds of the given x-point." << endl << endl;
	cout << "positional arguments:" << endl;
	cout << "  file          Contol file (starts with '_')" << endl;
	cout << "  fix_file      Output from dtfix (starts with 'fix')" << endl;
	cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
	cout << endl << "optional arguments:" << endl;
	cout << "  -h            show this help message and exit" << endl;
	cout << "  -s            stop at first boundary contact, default = No" << endl;
	cout << "  -m            max number of skip-boundary attempts, default = 3, ignored with -s option" << endl;
	cout << "  -b            boundary to stop at: " << endl;
	cout << "                    wall = vacuum vessel wall (default), " << endl;
	cout << "                    box  = simple rectangular box around vessel wall, " << endl;
	cout << "                    efit = limit of EFIT grid, forces -s option" << endl;
	cout << "  -n            Step-Size-Control: min step size, default = 0.001" << endl;
	cout << "  -x            Step-Size-Control: max step size, default = 0.005" << endl;
	cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
	cout << "  -S            filename for SIESTA; default, see below" << endl;
	cout << "  -V            filename for VMEC; default, see below" << endl;
	cout << "  -X            filename for XPAND; default, see below" << endl;
	cout << endl << "Examples:" << endl;
	cout << "  dtman _fix.dat fix_1_lower.dat blabla" << endl;
	cout << "  dtman -b efit _fix.dat fix_1_lower.dat testrun" << endl;
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
case 's':
	trytoskip = 0;
	break;
case 'm':
	skipmax = atoi(optarg);
	break;
case 'n':
	minabs = atof(optarg);
	break;
case 'x':
	maxabs = atof(optarg);
	break;
case 'b':
	if(LA_STRING(optarg) == "box") bndy_type = 1;
	if(LA_STRING(optarg) == "efit") {bndy_type = 2; trytoskip = 0;}
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
LA_STRING basename,startname;
LA_STRING praefix = "";
if(argc==optind+3) praefix = "_" + LA_STRING(argv[optind+2]);
if(argc>=optind+2)
{
	basename = LA_STRING(argv[optind]);
	startname = LA_STRING(argv[optind+1]);
}
if(argc<=optind+1)
{
	cout << "No Input files -> Abort!" << endl; 
	exit(0);
}
basename = checkparfilename(basename);
startname = checkparfilename(startname);
LA_STRING parfilename = "_" + basename + ".dat";
LA_STRING name = startname + ".dat";

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// Set starting parameters
const int kstart = 1;
const int kend = 30;

// read fixed points
Array<double,2> data;
readfile(name,6,data);

// Output
LA_STRING filenameout,type,dir;
if(PAR.MapDirection==1) type = "_unst";
else type = "_st";
if(PAR.verschieb>0) dir = "r";
else dir = "l";
ofstream out;
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "psi";  var[3] = "theta[rad]";  var[4] = "r[m]";

// log file
ofs2.open("log_" + LA_STRING(program_name) + type + dir + praefix + ".dat");
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

EQD.ReadData(EQD.Shot,EQD.Time,Raxis,Zaxis);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.Path << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.Path << endl;

// set boundary

if(bndy_type == 0) simpleBndy = 0;   // use real wall as boundary
else
{
	simpleBndy = 1;		// use simple boundary box
	if(bndy_type == 2) 	// extend boundary to EFIT boundary
	{
		bndy[0] = min(EQD.R);
		bndy[1] = max(EQD.R);
		bndy[2] = min(EQD.Z);
		bndy[3] = max(EQD.Z);
	}
}
cout << "Boundary (0 = Wall, 1 = Box, 2 = EFIT): " << bndy_type << endl;
cout << "Box limits: " << bndy[0] << "\t" << bndy[1] << "\t" << bndy[2] << "\t" << bndy[3] << endl;

// Prepare Perturbation
prepare_common_perturbations(EQD,PAR,0,siestafile,xpandfile,islandfile);
prep_perturbation(EQD,PAR);

// Prepare collisions
COLLISION COL;
if (use_collision) COL.init(TprofileFile, NprofileFile, f, zbar, PAR.Zq, PAR.Mass, rc, mc);
// Prepare particles
PARTICLE FLT(EQD,PAR,COL);
PARTICLE FLTold(EQD,PAR,COL);

// loop for all fixed points
for(i=1;i<=data.rows();i++)
{
	xfix = data(i,1);	 yfix = data(i,2);	periode = int(fabs(data(i,3)));
	xsa(1) = xfix;	 xsa(2) = yfix;	

	find_start_values(xsa,da,periode,PAR.verschieb,PAR.phistart,PAR.MapDirection,FLT); 
	
	ofs2 << "Period: " << periode << "\t" << "Shift: " << PAR.verschieb << endl;
	ofs2 << "fixed point: R = " << xfix << "\t" << "Z = " << yfix << endl;
	ofs2 << "k goes from k = " << kstart << " to k= " << kend << endl;
	ofs2 << "Distances: Min = " << minabs << "\t" << "Max= " << maxabs << endl;
	ofs2 << endl;

	// additional parameters for IO
	PAR.pv[0].name = "Period";				PAR.pv[0].wert = periode;
	PAR.pv[1].name = "fixed point R";		PAR.pv[1].wert = xfix;
	PAR.pv[2].name = "fixed point Z";		PAR.pv[2].wert = yfix;
	PAR.pv[3].name = "Shift in R";			PAR.pv[3].wert = PAR.verschieb;
	PAR.pv[4].name = "Minimal distance";	PAR.pv[4].wert = minabs;
	PAR.pv[5].name = "Maximal distance";	PAR.pv[5].wert = maxabs;
	PAR.pv[6].name = "phistart";			PAR.pv[6].wert = PAR.phistart;
	PAR.pv[7].name = "MapDirection";		PAR.pv[7].wert = PAR.MapDirection;
	PAR.pv[8].name = "Ekin";				PAR.pv[8].wert = PAR.Ekin;
	PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;

	// Output
	filenameout = "man" + type + dir + LA_STRING(periode) + "_" + LA_STRING(i) + praefix + ".dat";
	outputtest(filenameout);
	out.open(filenameout);
	PAR.writeiodata(out,bndy,var);

	// calculate Manifold
	end = 0;  skipattempt = 0;  outside = 0;
	plotchk = 0;
	variate = 0;
	xa = xsa + da;
	FLTold.R = xa(1);
	FLTold.Z = xa(2);
	dt = 0.1;

	for(k=kstart;k<=kend;k++)
	{
		t = 0;
		talt = t;
		//dt = 0.1;

		while(t<1)
		{
			t += dt;
			if(t>1) 
			{
				dt = 1 - talt;
				t = 1;
			}
			xa = xsa + t*da;
			//xitta = xa;

			FLT.R = xa(1);
			FLT.Z = xa(2);
			FLT.phi = PAR.phistart;
			chk = FLT.mapit(k*2*periode,PAR.MapDirection);
			if(chk<0) 	// outside wall, try to skip outside part, no step size management
			{
				skipattempt = 1; 
				if(outside >= skipmax) {end = 1; ofs2 << "final wall hit" << endl; break;} // stop after skipmax skip attempts
				if(trytoskip==1) continue;
				else {end = 1; ofs2 << "wall hit" << endl; break;}
			}
			if(skipattempt==1) {skipattempt = 0; outside += 1; variate = 1;}

			// step size management
			dist = abstand(FLT,FLTold);
			if(variate==0 && dist>maxabs)
			{
				t = talt;
				dt *= 0.5;
				if(dt < 1e-8 && preventSmallSteps == 1)	// forced to execute one step with dt = 1e-6
					{dt = 1e-6; variate = 1; ofs2 << "Step size too small -> try to skip" << endl;}
				if(dt < 1e-14){end = 1; ofs2 << "Step size cannot be further decreased" << endl; break;}
				continue;
			}

			// Output 
			if(plotchk==0) {ofs2 << "Start plotting..." << endl; plotchk = 1;}
			out << FLT.R << "\t" << FLT.Z << "\t" << FLT.psi << "\t" << FLT.get_theta() << "\t" << FLT.get_r() << endl;
		
			// step size management part 2
			if(dist<minabs)
			{
				dt *= 2.5;
			}

			talt = t;
			FLTold = FLT;
			variate = 0;
		} // end while

		if(end==1) break;
		ofs2 << "k= " << k << "\t" << flush; 
	} // end for k
	ofs2 << endl;
	out.close();
} // end for i
ofs2 << "Program terminates normally" << endl;
cout << "Program terminates normally" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

return 0;
} //end of main

//------------------------ End of Main -------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------

//---------- find_start_values ----------------
void find_start_values(Array<double,1>& xsa, Array<double,1>& da, int periode, double& verschieb, 
					   double phistart, int MapDirection, PARTICLE& FLT)
{
int chk;
double laenge = 0,scalar;

Array<double,1> d2a(Range(1,2));
Array<double,1> fixa(Range(1,2));
Array<double,1> xsitta(Range(1,2));

fixa = xsa;
xsitta = xsa;

da(2) = fabs(verschieb)/sqrt(2); da(1)= verschieb/sqrt(2);	// shift in R and Z, Z always positive
//da(2) = 0; da(1) = verschieb;	// shift in R

// adjust shift that way, that xsa and xsitta are between 0.0005 and 0.001
int end = 0;
while(end<50)
{
	xsa = fixa + da;
	FLT.R = xsa(1);
	FLT.Z = xsa(2);
	FLT.phi = phistart;
	chk = FLT.mapit(2*periode,MapDirection);
	if(chk<0) {end += 1; da *= 0.5; continue;}
	xsitta(1) = FLT.R;
	xsitta(2) = FLT.Z;


	//ofs2 << xsa(1) << "\t" << xsa(2) << "\t" << xsitta(1) << "\t" << xsitta(2) << "\t" << da(1) << "\t" << da(2) << endl;

	laenge = abstand(xsa,xsitta);
	end += 1;
	if(laenge > 0.001) {da *= 0.5; continue;}
	if(laenge < 0.0005){da *= 1.1; continue;}
	break;
}
if(end==50) ofs2 << "Shift adjustment not successfull! " << laenge << endl;
verschieb = sqrt(da(1)*da(1) + da(2)*da(2))*sign(da(1));

// adjust direction of d
end = 0;  scalar = 0;
while(1-scalar > 1e-5)
{
	// set d
	da = xsitta - fixa;
	//if(fabs(da(1)) > 4) {da(1) -= sign(da(1))*pi2;}	// possible 2pi jump between xsitt and fix.
	da *= fabs(verschieb) / sqrt(da(1)*da(1)+da(2)*da(2));		// set length of d to shift

	// check new direction
	xsa = fixa + da;
	FLT.R = xsa(1);
	FLT.Z = xsa(2);
	FLT.phi = phistart;
	chk = FLT.mapit(2*periode,MapDirection);
	if(chk<0) {ofs2 << "Error during direction adjustment: " << chk << endl; exit(0);}
	xsitta(1) = FLT.R;
	xsitta(2) = FLT.Z;

	// terminate condition.
	d2a = xsitta - xsa;
	//if(fabs(d2a(1)) > 4) {d2a(1) -= sign(d2a(1))*pi2;}
	scalar = (da(1)*d2a(1)+da(2)*d2a(2)) / sqrt(da(1)*da(1)+da(2)*da(2)) / sqrt(d2a(1)*d2a(1)+d2a(2)*d2a(2));

	end += 1;
	if(end>50) {ofs2 << "Direction adjustment not successfull" << endl; break;}

	//ofs2 << xsa(1) << "\t" << xsa(2) << "\t" << xsitta(1) << "\t" << xsitta(2) << "\t" << da(1) << "\t" << da(2) << endl;
}
ofs2 << "Effort of direction adjustment: " << end << endl;

// set d to be vector from xs to xsitt
da = xsitta - xsa;
//if(fabs(da(1)) > 4) {da(1) -= sign(da(1))*pi2;}	// possible 2pi jump between xsitt and fix.
ofs2 << "Length of d: " << sqrt(da(1)*da(1)+da(2)*da(2)) << endl;

}

//---------- abstand ----------------------
// Calculates distance between (R,Z) positions of two PARTICLE objects
double abstand(PARTICLE& FLT1, PARTICLE& FLT2)
{
const double dR = FLT1.R - FLT2.R;
const double dZ = FLT1.Z - FLT2.Z;
return sqrt(dR*dR + dZ*dZ);
}

//---------- abstand ----------------------
// Calculates distance between two (R,Z) arrays
double abstand(Array<double,1>& x1, Array<double,1>& x2)
{
const double dR = x1(1) - x2(1);
const double dZ = x1(2) - x2(2);
return sqrt(dR*dR + dZ*dZ);
}

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

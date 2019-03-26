// Program calculates last closed flux surface for particle-drift with time dependent perturbations
// Fortran subroutines for Perturbation are used
// A.Wingen						19.3.19

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	Poincare particle drift data file
//			log-file

// Define
//--------
#if defined(ITER)
	#define program_name "iterlcfs"
#elif defined(NSTX)
	#define program_name "nstxlcfs"
#elif defined(MAST)
	#define program_name "mastlcfs"
#elif defined(CMOD)
	#define program_name "cmodlcfs"
#else
	#define program_name "dtlcfs"
#endif

// Include
//--------
#include <mafot.hxx>
#include <unistd.h>

// namespaces
//-----------
using namespace blitz;

// Switches
//----------

// Golbal Parameters
//------------------

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// MPI initialize
int mpi_rank = 0;

// Variables
int i,j,n,chk;
EFIT EQD;
Range all = Range::all();
double s,u;

// Use system time as seed(=idum) for random numbers
double now=zeit();

// defaults
LA_STRING woutfile = "wout.nc";
LA_STRING xpandfile = "None";
LA_STRING siestafile = "siesta.dat";
LA_STRING islandfile = "fakeIslands.in";
LA_STRING ErProfileFile = "None";
bool use_ErProfile = false;
LA_STRING TprofileFile = "None";
bool use_Tprofile = false;

// Command line input parsing
int c;
opterr = 0;
while ((c = getopt(argc, argv, "hX:V:S:I:E:T:")) != -1)
switch (c)
{
case 'h':
	cout << "usage: dtlcfs [-h] [-E ErProfile] [-I island] [-S siesta] [-T Tprofile] [-V wout] [-X xpand] file [tag]" << endl << endl;
	cout << "Calculate the last closed flux surface in a Poincare Plot." << endl << endl;
	cout << "positional arguments:" << endl;
	cout << "  file          Contol file (starts with '_')" << endl;
	cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
	cout << endl << "optional arguments:" << endl;
	cout << "  -h            show this help message and exit" << endl;
	cout << "  -E            use electric field with particle drifts. Filename of Er(psi) profile." << endl;
	cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
	cout << "  -S            filename for SIESTA; default, see below" << endl;
	cout << "  -T            use temperature profile with particle drifts. Filename of T(psi) profile." << endl;
	cout << "  -V            filename for VMEC; default, see below" << endl;
	cout << "  -X            filename for XPAND; default is None" << endl;
	cout << endl << "Examples:" << endl;
	cout << "  dtlcfs _plot.dat blabla" << endl;
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
case 'E':
	ErProfileFile = optarg;
	use_ErProfile = true;
	break;
case 'T':
	TprofileFile = optarg;
	use_Tprofile = true;
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
if(argc==optind+2) praefix = "_" + LA_STRING(argv[optind+1]);
if(argc>=optind+1) basename = LA_STRING(argv[optind]);
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + ".dat");
ofs2.precision(16);

// Output
LA_STRING filenameout = "lcfs" + praefix + ".dat";
outputtest(filenameout);

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,9);

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

// Read E-field data
if(use_ErProfile)
{
	EQD.ReadEfield(ErProfileFile);
	PAR.useErProfile = 1;
	cout << "Read Er profile: " << ErProfileFile << endl;
	ofs2 << "Read Er profile: " << ErProfileFile << endl;
}

// Read T-profile data
if(use_Tprofile)
{
	EQD.ReadTprofile(TprofileFile);
	PAR.useTprofile = 1;
	cout << "Ignoring Ekin, read T profile: " << TprofileFile << endl;
	ofs2 << "Ignoring Ekin, read T profile: " << TprofileFile << endl;
}

// Prepare particles
PARTICLE FLT(EQD,PAR);

// additional parameters for IO
PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
PAR.pv[1].name = "Refinements";		PAR.pv[1].wert = PAR.N;
PAR.pv[5].name = "phistart";		PAR.pv[5].wert = PAR.phistart;
PAR.pv[6].name = "MapDirection";	PAR.pv[6].wert = PAR.MapDirection;
PAR.pv[7].name = "Ekin";			PAR.pv[7].wert = PAR.Ekin;
PAR.pv[8].name = "energy ratio lambda";	PAR.pv[8].wert = PAR.lambda;

PAR.pv[2].name = "psimin";			PAR.pv[2].wert = PAR.rmin;
PAR.pv[3].name = "psimax";			PAR.pv[3].wert = PAR.rmax;
PAR.pv[4].name = "theta";			PAR.pv[4].wert = PAR.thmin;

// Output
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(6);
var[0] = "theta[rad]"; var[1] = "r[m]"; var[2] = "phi[deg]"; var[3] = "psi"; var[4] = "R[m]"; var[5] = "Z[m]";
if(PAR.response_field == -2) {var[0] = "u"; var[3] = "s";}
PAR.writeiodata(out,bndy,var);

cout << "Helicity = " << EQD.helicity << endl;
ofs2 << "Helicity = " << EQD.helicity << endl;

// Result array:	 Column Number,  Values
Array<double,2> results(Range(1,6),Range(1,PAR.itt));
Array<double,2> temp(Range(1,6),Range(1,PAR.itt));


// Prepare Perturbation
prepare_common_perturbations(EQD,PAR,0,siestafile,xpandfile,islandfile);
prep_perturbation(EQD,PAR);

// each process gets a different seed
long idum=long(now);

cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;

// Get last closed flux surface
int Nmax = PAR.N;			// number of refinements
double psi0 = PAR.rmin;
double psi1 = PAR.rmax;
double th0 = pi - 0.1;
double th1 = pi + 0.1;
double psix,psimin,psimax;
for(n=1;n<=Nmax;n++)
{
	// Set initial values in FLT
	psix = 0.5*(psi1 + psi0);
	FLT.set(1,1,psix,psix,PAR.thmin,PAR.thmin,1,2);

	// use psi as control parameter and set to -1
	temp(4,all) = -1;
	psimin = 1;
	psimax = 0;

	// Integrate
	for(i=1;i<=PAR.itt;i++)
	{
		chk = FLT.mapstep(PAR.MapDirection);
		if(chk<0) {break;}	// particle has left system

		//if(fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) << endl;
		FLT.phi=PAR.MapDirection*i*dpinit*ilt + PAR.phistart;

		// Store results
#ifdef USE_SIESTA
		if(PAR.response_field == -2)
		{
			SIES.get_su(FLT.R, FLT.phi/rTOd, FLT.Z, s, u);
			temp(1,i) = u;
			temp(4,i) = s;
		}
		else
		{
			temp(1,i) = FLT.get_theta();
			temp(4,i) = FLT.psi;
		}
#else
		temp(1,i) = FLT.get_theta();
		temp(4,i) = FLT.psi;
#endif
		temp(2,i) = FLT.get_r();
		temp(3,i) = FLT.phi;
		temp(5,i) = FLT.R;
		temp(6,i) = FLT.Z;

		// check for spread in psi
		if ((temp(1,i) >= th0) && (temp(1,i) <= th1))
		{
			if (FLT.psi < psimin) psimin = FLT.psi;
			if (FLT.psi > psimax) psimax = FLT.psi;
		}
		if ((psimin < 1) && (psimax > 0))
		{
			if((psimax - psimin) > 0.03) {chk = -2; break;}
		}
	} // end for i

	if (chk < 0) psi1 = psix;
	else
	{
		psi0 = psix;
		results = temp.copy();
	}

	cout << "Step: " << n << "/" << Nmax << endl;
	ofs2 << "----------------------------------------- Step: " << n << "/" << Nmax << endl;
	ofs2 << psix << "\t" << i << "\t" << chk << endl;
} // end for n

// Output
for(j=1;j<=PAR.itt;j++)	if(results(4,j) > 0) out << results(1,j) << "\t" << results(2,j) << "\t" << results(3,j) << "\t" << results(4,j) << "\t" << results(5,j) << "\t" << results(6,j) << endl;


double now2=zeit();
cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

return 0;
} //end main

//----------------------- End of Main -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

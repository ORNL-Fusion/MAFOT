// Program calculates connection length and penetration depth for D3D inside the plasma volume
// for D3D-Drift with Time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						20.06.11

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	2d connection length data for colored contour plot
//			log-file


// Parallel computation possible by slicing the Z-direction -> easy reconnection of files

// Define
//--------
//#define BZ_DEBUG
#define program_name "dtlaminar"

// Include
//--------
#include <mafot.hxx>
#include <d3d.hxx>

// Prototypes  

// Switches
const int spare_interior = 0;	// 0: all points are calculated		1: inside psi=0.95 results are set to fixed values (code runs faster)

// Golbal Parameters 

// Main Program
//--------------
int main(int argc, char *argv[])
{
// Variables
int i,j;
int chk;
double Rout,Zout;
double ntor,length,psimin;
EFIT EQD;

// Use system time as seed(=idum) for random numbers
double now = zeit();
long idum = long(now);

// Input file names
LA_STRING basename;
LA_STRING praefix = "";
if(argc==3) praefix = "_" + LA_STRING(argv[2]);
if(argc>=2) basename = LA_STRING(argv[1]);
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

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,11);

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
int N = PAR.NR*PAR.NZ;

// additional parameters for IO
PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
PAR.pv[1].name = "R-grid";			PAR.pv[1].wert = PAR.NR;
PAR.pv[2].name = "Z-grid";			PAR.pv[2].wert = PAR.NZ;
PAR.pv[3].name = "Rmin";			PAR.pv[3].wert = PAR.Rmin;
PAR.pv[4].name = "Rmax";			PAR.pv[4].wert = PAR.Rmax;
PAR.pv[5].name = "Zmin";			PAR.pv[5].wert = PAR.Zmin;
PAR.pv[6].name = "Zmax";			PAR.pv[6].wert = PAR.Zmax;
PAR.pv[7].name = "phistart";		PAR.pv[7].wert = PAR.phistart;
PAR.pv[8].name = "MapDirection";	PAR.pv[8].wert = PAR.MapDirection;
PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;
PAR.pv[10].name = "Ekin";			PAR.pv[10].wert = PAR.Ekin;

// Prepare Perturbation
prep_perturbation(EQD,PAR);

// Prepare particles
PARTICLE FLT(EQD,PAR);

// Output
LA_STRING filenameout = "lam" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";
PAR.writeiodata(out,bndy,var);

cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
cout << "Start Tracer for " << N << " points ... " << endl;
ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << PAR.MapDirection << endl;
ofs2 << "Start Tracer for " << N << " points ... " << endl;
for(i=1;i<=N;i++)
{
	// Set and store initial condition
	FLT.set(i,N,PAR.Rmin,PAR.Rmax,PAR.Zmin,PAR.Zmax,PAR.NZ);
	Rout = FLT.R;
	Zout = FLT.Z;

	// Spare the calculation of the interior
	if(spare_interior == 1 && FLT.psi <= 0.95 && FLT.Z > -1.25) 
	{
		ntor = 2*PAR.itt;
		length = 4000.0;
		psimin = 0.95;
	}
	else chk = FLT.connect(ntor,length,psimin,PAR.itt,PAR.MapDirection);

	out << Rout << "\t" << Zout << "\t" << ntor << "\t" << length/1000.0 << "\t" << psimin << endl;

	if(i%100==0) ofs2 << "Trax: " << i << endl;
	//if(i%100==0) {cout << "Trax: " << i << endl;}
}

double now2 = zeit();
cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

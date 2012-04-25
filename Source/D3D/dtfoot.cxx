// Program calculates connection length and penetration depth for D3D
// include drifts and time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						20.06.11

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	2d footprint data for colored contour plot
//			log-file


// Parallel computation possible by slicing the t-direction -> easy reconnection of files


// Define
//--------
//#define BZ_DEBUG
#define program_name "dtfoot"

// Include
//--------
#include <mafot.hxx>
#include <d3d.hxx>

// Prototypes  

// Switches

// Golbal Parameters 

// Main Program
//--------------
int main(int argc, char *argv[])
{
// Variables
int i;
int chk;
double phiout;
double t,ntor,length,psimin;
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

// Read Parameterfile
cout << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// Set target type for output-filename
LA_STRING type;
switch(PAR.which_target_plate)
{
case 0:
	type = "_cp";
	break;
case 1:
	type = "_in";
	break;
case 2:
	type = "_out";
	break;
case 3:
	type = "_shelf";
	break;
default:
	type = "";
	break;
}

// log file
ofs2.open("log_" + LA_STRING(program_name) + type + praefix + ".dat");
ofs2.precision(16);
ofs2 << "Read Parameterfile " << parfilename << endl;

// Set starting parameters
int N = PAR.Nt*PAR.Nphi;

// Use Boundary Box and extend lower boundary to prevent horizontal plate to be outside of boundary
simpleBndy = 1;
bndy[2] = -1.4;	 // originally  bndy[2] = -1.367 and plate at -1.3664-Z0

// additional parameters for IO
PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
PAR.pv[1].name = "t-grid";			PAR.pv[1].wert = PAR.Nt;
PAR.pv[2].name = "phi-grid";		PAR.pv[2].wert = PAR.Nphi;
PAR.pv[3].name = "tmin";			PAR.pv[3].wert = PAR.tmin;
PAR.pv[4].name = "tmax";			PAR.pv[4].wert = PAR.tmax;
PAR.pv[5].name = "phimin";			PAR.pv[5].wert = PAR.phimin;
PAR.pv[6].name = "phimax";			PAR.pv[6].wert = PAR.phimax;
PAR.pv[7].name = "MapDirection";	PAR.pv[7].wert = PAR.MapDirection;
PAR.pv[8].name = "Ekin";			PAR.pv[8].wert = PAR.Ekin;
PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prep_perturbation(EQD,PAR);

// Prepare particles
PARTICLE FLT(EQD,PAR);

// Output
LA_STRING filenameout = "foot" + type + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "phi[rad]";  var[1] = "length t";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";
PAR.writeiodata(out,bndy,var);

cout << "Target (0=CP, 1=inner, 2=outer, 3=shelf): " << PAR.which_target_plate << endl;
cout << "Start Tracer for " << N << " points ... " << endl;
ofs2 << "Target (0=CP, 1=inner, 2=outer, 3=shelf): " << PAR.which_target_plate << endl;
ofs2 << "Start Tracer for " << N << " points ... " << endl;
for(i=1;i<=N;i++)
{
	// Set and store initial condition
	t = start_on_target(i,PAR.Nt,PAR.Nphi,PAR.tmin,PAR.tmax,PAR.phimin,PAR.phimax,EQD,PAR,FLT);
	phiout = FLT.phi/rTOd;	//phi in rad;

	chk = FLT.connect(ntor,length,psimin,PAR.itt,PAR.MapDirection);

	out << phiout << "\t" << t << "\t" << ntor << "\t" << length/1000.0 << "\t" << psimin << endl;
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

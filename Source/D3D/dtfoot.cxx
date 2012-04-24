// Program calculates connection length and penetration depth for D3D
// include drifts and time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						12.01.09

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
#include <andi.hxx>
#include <efit_class.hxx>
#include <d3d-drift.hxx>

// Prototypes  

// Switches
const int useparfile = 1;	// 0: additional parameters are set in the code		1: All parameters are read from file

// Golbal Parameters 
EFIT EQD;
double GAMMA;
double eps0;
double Ix;

// Main Program
//--------------
int main(int argc, char *argv[])
{
// Variables
int i;
int chk;
int MapDirection;	// set within the Code, depending on which_target_plate; 1: positive phi-direction	-1: negative phi-direction
double phi,phistart,psimin;
double t,ntor,length;
double omegac;
Array<double,1> xa(Range(1,2));

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
vector<double> startvec;
cout << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);

// Set target type for output-filename
LA_STRING type;
switch(which_target_plate)
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
int itt = 15;
double tmin = 0.76;
double tmax = 0.85;
double phimin = 0;
double phimax = pi2;
int Np = 50;
int Nphi = 20;

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

if(useparfile==1)
{
	cout << "All parameters are read from file" << endl;
	ofs2 << "All parameters are read from file" << endl;
	itt = int(startvec[1]);
	tmin = startvec[2];
	tmax = startvec[3];
	phimin = startvec[4];
	phimax = startvec[5];
	Np = int(startvec[6]);
	Nphi = int(startvec[0]);
	Ekin = startvec[18];
	lambda = startvec[19];
}
int N = Np*Nphi;
if(which_target_plate==2 || which_target_plate==3) MapDirection = 1;
else  MapDirection = -1;

// Extend lower boundary to prevent horizontal plate to be outside of boundary
bndy[2] = -1.4;	 // originally  bndy[2] = -1.367 and plate at -1.3664-Z0

// additional parameters for IO
const int psize=10;
parstruct * parvec = new parstruct[psize];
parvec[0].name = "Max. Iterations";	parvec[0].wert = itt;
parvec[1].name = "t-grid";			parvec[1].wert = Np;
parvec[2].name = "phi-grid";		parvec[2].wert = Nphi;
parvec[3].name = "tmin";			parvec[3].wert = tmin;
parvec[4].name = "tmax";			parvec[4].wert = tmax;
parvec[5].name = "phimin";			parvec[5].wert = phimin;
parvec[6].name = "phimax";			parvec[6].wert = phimax;
parvec[7].name = "MapDirection";	parvec[7].wert = MapDirection;
parvec[8].name = "Ekin";			parvec[8].wert = Ekin;
parvec[9].name = "energy ratio lambda";	parvec[9].wert = lambda;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prep_perturbation();

// Prepare particles
if(sigma==0)
{
	cout << "Field lines are calculated" << endl;
	ofs2 << "Field lines are calculated" << endl;
	GAMMA = 1;	omegac = 1;	 eps0 = 1;	Ix = 1;
}
else
{
	if(Zq>=1) // Ions
	{
		GAMMA = 1 + Ekin/(E0p*Massnumber);	// relativistic gamma factor 1/sqrt(1-v^2/c^2)
		omegac = e*EQD.Bt0/(mp*Massnumber);	// normalized gyro frequency (SI-System)
		cout << "Ions are calculated" << endl;
		ofs2 << "Ions are calculated" << endl;
	}
	else // Electrons
	{
		Zq = -1;	// default!
		GAMMA = 1 + Ekin/(E0e*Massnumber);	// relativistic gamma factor 1/sqrt(1-v^2/c^2)
		omegac = e*EQD.Bt0/(me*Massnumber);	// normalized gyro frequency (SI-System)
		cout << "Electrons are calculated" << endl;
		ofs2 << "Electrons are calculated" << endl;
	}
	eps0 = c*c/omegac/omegac/EQD.R0/EQD.R0;	// normalized rest energy
	Ix = -0.5/double(Zq)*eps0*((lambda*(GAMMA-1)+1)*(lambda*(GAMMA-1)+1)-1);
	cout << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
	ofs2 << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
}

// Output
LA_STRING filenameout = "foot" + type + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "phi[rad]";  var[1] = "length t";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";
writeiodata(out,var,parvec,psize,parfilename);

cout << "Target (0=vertical, 1=45°, 2=horizontal, 3=shelf): " << which_target_plate << endl;
cout << "Start Tracer for " << N << " points ... " << endl;
ofs2 << "Target (0=vertical, 1=45°, 2=horizontal, 3=shelf): " << which_target_plate << endl;
ofs2 << "Start Tracer for " << N << " points ... " << endl;
for(i=1;i<=N;i++)
{
	t = start_on_target(i,Np,Nphi,tmin,tmax,phimin,phimax,xa(2),xa(1),phistart);
	phi = phistart*rTOd;	//phi in deg;

	chk = connect(xa,phi,itt,MapDirection,ntor,length,psimin);

	out << phistart << "\t" << t << "\t" << ntor << "\t" << length/1000.0 << "\t" << psimin << endl;
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

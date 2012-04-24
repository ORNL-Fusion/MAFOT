// Program calculates connection length and penetration depth for D3D inside the plasma volume
// for D3D-Drift with Time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						13.01.09

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
#include <andi.hxx>
#include <efit_class.hxx>
#include <d3d-drift.hxx>

// Prototypes  
int follow(double& R, double& Z, double& phi, double& psi, int MapDirection);

// Switches
const int useparfile = 1;	// 0: additional parameters are set in the code		1: All parameters are read from file
const int spare_interior = 0;	// 0: all points are calculated		1: inside psi=0.95 results are set to fixed values (code runs faster)

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
int i,j;
int chk;
double R,Z,phi,psi,psimin;
double ntor,length,dummy;
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

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + ".dat");
ofs2.precision(16);

// Read parameter file
vector<double> startvec;
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
int itt = 40;
double Rmin = 1.016 - EQD.RmAxis;
double Rmax = 1.45 - EQD.RmAxis;
double Zmin = -1.367 - EQD.ZmAxis;
double Zmax = -0.5 - EQD.ZmAxis;
int NR = 200;
int NZ = 300;
double phistart = 0;
int MapDirection = 0;	// 1: positive phi-direction	-1: negative phi-direction	0: both directions

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

if(useparfile==1)
{
	cout << "All parameters are read from file" << endl;
	ofs2 << "All parameters are read from file" << endl;
	itt = int(startvec[1]);
	Rmin = startvec[2] - EQD.RmAxis;
	Rmax = startvec[3] - EQD.RmAxis;
	Zmin = startvec[4] - EQD.ZmAxis;
	Zmax = startvec[5] - EQD.ZmAxis;
	NR = int(startvec[6]);
	NZ = int(startvec[0]);
	phistart = startvec[7];
	MapDirection = int(startvec[8]);
	Ekin = startvec[18];
	lambda = startvec[19];
}
int N = NR*NZ;

// additional parameters for IO
const int psize = 11;
parstruct * parvec = new parstruct[psize];
parvec[0].name = "Max. Iterations";	parvec[0].wert = itt;
parvec[1].name = "R-grid";			parvec[1].wert = NR;
parvec[2].name = "Z-grid";			parvec[2].wert = NZ;
parvec[3].name = "Rmin";			parvec[3].wert = Rmin + EQD.RmAxis;
parvec[4].name = "Rmax";			parvec[4].wert = Rmax + EQD.RmAxis;
parvec[5].name = "Zmin";			parvec[5].wert = Zmin + EQD.ZmAxis;
parvec[6].name = "Zmax";			parvec[6].wert = Zmax + EQD.ZmAxis;
parvec[7].name = "phistart";		parvec[7].wert = phistart;
parvec[8].name = "MapDirection";	parvec[8].wert = MapDirection;
parvec[9].name = "energy ratio lambda";	parvec[9].wert = lambda;
parvec[10].name = "Ekin";			parvec[10].wert = Ekin;

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
LA_STRING filenameout = "lam" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "N_toroidal";  var[3] = "connection length [km]";  var[4] = "psimin (penetration depth)";
writeiodata(out,var,parvec,psize,parfilename);

cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;
cout << "Start Tracer for " << N << " points ... " << endl;
ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;
ofs2 << "Start Tracer for " << N << " points ... " << endl;
for(i=1;i<=N;i++)
{
	set(i,N,Zmin,Zmax,Rmin,Rmax,Z,R,NZ);	// swap x and y, because matlab requires R to vary first
	xa(1) = atan(Z/R);	// theta
	if(R<0) xa(1) += pi;
	if(R>=0 && Z<0) xa(1) += pi2;
	xa(2) = sqrt(R*R+Z*Z);	// r
	phi = phistart;

	// Spare the calculation of the interior
	EQD.get_psi(R+EQD.RmAxis,Z+EQD.ZmAxis,psi,dummy,dummy); // get psi
	if(spare_interior == 1 && psi <= 0.95 && Z+EQD.ZmAxis > -1.25) 
	{
		ntor = 2*itt;
		length = 4000.0;
		psimin = 0.95;
	}
	else chk = connect(xa,phi,itt,MapDirection,ntor,length,psimin);

	out << R + EQD.RmAxis << "\t" << Z + EQD.ZmAxis << "\t" << ntor << "\t" << length/1000.0 << "\t" << psimin << endl;

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

//----------------- follow ----------------------------
int follow(double& R, double& Z, double& phi, double& psi, int MapDirection)
{
int chk;
double dummy;
double dphi = MapDirection*dpinit/rTOd;
double phi_rad = phi/rTOd;	// phi in rad

Array<double,1> y(2);
y(0) = R;	// R
y(1) = Z;	// Z

chk = rkint(nvar,10,dphi,y,phi_rad);
if(chk<0) return -1;	// particle has left system

// return coordinates
R = y(0); // R 
Z = y(1); // Z 

// phi back in deg
phi = phi_rad*rTOd;

// Get normalized flux
EQD.get_psi(y(0),y(1),psi,dummy,dummy);

return 0;
}


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

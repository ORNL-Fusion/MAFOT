// Program calculates path of field lines in 10 deg steps. Used for 3-D Visualization
// for D3D-Drift with Time dependent perturbations
// Fortran Subroutines are used for perturbations
// A.Wingen						4.02.09

// Input: 1: Parameterfile	2: Initial value file (optional; if not, enter x instead) 3: praefix (optional)
// Output:	paths of individual field lines in 10 deg steps between the target plates
//			log-file


// Define
//--------
//#define BZ_DEBUG
#define program_name "dtstructure"

// Include
//--------
#include <andi.hxx>
#include <efit_class.hxx>
#include <d3d-drift.hxx>

// Prototypes 
int follow(double& R, double& Z, double& phi, double& psi, int MapDirection);

// Switches
const int useparfile=1;	// 0: additional parameters are set in the code		1: All parameters are read from file

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
int i,j,k,chk;
double Z,R,phi,psi;
double omegac;
int usePointfile;

// Use system time as seed(=idum) for random numbers
double now = zeit();
long idum = long(now);

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
vector<double> startvec;
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Set starting parameters
double dphi = 10;	// in deg
int itt = 40;
double Rmin = 1.016;	// Startpoint of line
double Zmin = -1.367;
double Rmax = 1.45;	// Endpoint of line
double Zmax = -0.5;
int N = 10;		// Number of points on line, including start and end points
double phistart = 0;
int MapDirection = 0;	// 1: positive phi-direction	-1: negative phi-direction	0: both directions
double alpha = 0;	//Scales parabolic deformation of line between start and end point, alpha = 0: no deformation, alpha = 1: max deformation equals distance between points 

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

if(useparfile==1)
{
	cout << "All parameters are read from file" << endl;
	ofs2 << "All parameters are read from file" << endl;
	itt = int(startvec[1]);
	Rmin = startvec[2];
	Zmin = startvec[4];
	Rmax = startvec[3];
	Zmax = startvec[5];
	N = int(startvec[6]);
	phistart = startvec[7];
	MapDirection = int(startvec[8]);
	Ekin = startvec[18];
	lambda = startvec[19];
	alpha = startvec[0];
}

// additional parameters for IO
int psize = 10;
parstruct * parvec = new parstruct[psize];
parvec[0].name = "Max. Iterations";	parvec[0].wert = itt;
parvec[1].name = "phistart";		parvec[1].wert = phistart;
parvec[2].name = "MapDirection";	parvec[2].wert = MapDirection;
parvec[3].name = "energy ratio lambda";	parvec[3].wert = lambda;
parvec[4].name = "Ekin";			parvec[4].wert = Ekin;
parvec[5].name = "Number of Points";parvec[5].wert = N;
parvec[6].name = "Rmin";			parvec[6].wert = Rmin;
parvec[7].name = "Zmin";			parvec[7].wert = Zmin;
parvec[8].name = "Rmax";			parvec[8].wert = Rmax;
parvec[9].name = "Zmax";			parvec[9].wert = Zmax;

// Read pointfile 
Array<double,2> initial;
double dR,dZ,dN;
if(usePointfile == 1) 
{
	psize = 6;	// only the parameters 0 - 5 are written to file header
	readfile(pointname,2,initial);
	N = initial.rows();
	parvec[5].wert = N;
}
else	// or construct initial points from staight line between (Rmin,Zmin) and (Rmax,Zmax)
{
	initial.resize(Range(1,N),Range(1,2));	// secondIndex = 1: R	secondIndex = 2: Z

	if(N>1) dN = 1.0/double(N-1);
	else dN = 0;
	dR = (Rmax-Rmin)*dN;
	dZ = (Zmax-Zmin)*dN;

	for(i=1;i<=N;i++) 
	{
		initial(i,1) = Rmin + (i-1)*(dR - 4*alpha*((i-1)*dN-1)*dZ);	
		initial(i,2) = Zmin + (i-1)*(dZ + 4*alpha*((i-1)*dN-1)*dR);
	}
}
//cout << initial << endl;

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
LA_STRING filenameout = "struct" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "X[m]";  var[1] = "Y[m]";  var[2] = "Z[m]";  var[3] = "R[m]";  var[4] = "phi[rad]";
writeiodata(out,var,parvec,psize,parfilename);

cout << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;
cout << "Start Tracer for " << N << " points ... " << endl;
ofs2 << "MapDirection(0=both, 1=pos.phi, -1=neg.phi): " << MapDirection << endl;
ofs2 << "Start Tracer for " << N << " points ... " << endl;

// Follow the field lines
int size = itt*int(360.0/double(dphi));
Array<double,2> data(Range(1,size),Range(1,3));
for(i=1;i<=N;i++)
{
	R = initial(i,1);
	Z = initial(i,2);
	phi = phistart;

	// negative direction
	for(j=1;j<=size;j++)
	{
		chk = follow(R,Z,phi,psi,-1);	
		if(fabs(phi+j*dpinit*dphi-phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(phi+j*dpinit*dphi-phistart) << endl;
		phi = -j*dpinit*dphi+phistart;

		//Store Values
		data(j,1) = R;	data(j,2) = Z;	data(j,3) = phi;
		if(chk==-1) break;
	}
	// Write stored values in reverse direction
	for(k=j-1;k>=1;k--) out << data(k,1)*cos(data(k,3)/rTOd) << "\t" << data(k,1)*sin(data(k,3)/rTOd) << "\t" << data(k,2) << "\t" << data(k,1) << "\t" << data(k,3)/rTOd << endl;

	// Restore start values and write them
	R = initial(i,1);	Z = initial(i,2);	phi = phistart;
	out << R*cos(phi/rTOd) << "\t" << R*sin(phi/rTOd) << "\t" << Z << "\t" << R << "\t" << phi/rTOd << endl;

	//positive direction
	for(j=1;j<=size;j++)
	{
		chk = follow(R,Z,phi,psi,1);	
		if(fabs(phi-j*dpinit*dphi-phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(phi-j*dpinit*dphi-phistart) << endl;
		phi = j*dpinit*dphi+phistart;

		out << R*cos(phi/rTOd) << "\t" << R*sin(phi/rTOd) << "\t" << Z << "\t" << R << "\t" << phi/rTOd << endl;
		if(chk==-1) break;
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
for(j=1;j<=data.rows();j++) out << data(j,1) << "\t" << data(j,2) << "\t" << data(j,3) << "\t" << data(j,4) << "\t" << data(j,5) << endl;

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



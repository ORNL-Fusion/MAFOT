// Program searches for periodic fixed points in DIII-D Poincare-Plot
// include drifts and time dependent perturbations
// derivatives are calculated numerically using a 5-point stencil
// 2-dimensional Newton method is used
// Fortran Subroutines are used for perturbations
// A.Wingen						12.01.09

// Input: 1: Parameterfile	2: period of fixed point	3: praefix (optional)
// Output:	fixed points
//			log-file

// Define
//--------
//#define BZ_DEBUG
#define program_name "dtfix"

// Include
//--------
#include <andi.hxx>
#include <efit_class.hxx>
#include <d3d-drift.hxx>

// Prototypes  
int check_boundary(Array<double,1> x);
int newton2D(Array<double,1>& x, double phistart, double& psi, int periode);
int mapit_J(Array<double,1>& xa, double& phi, double& psi, Array<double,1>& J, int itt, int Map=1);

// Switches
const int useparfile = 1;		// 0: additional parameters are set in the code		1: All parameters are read from file

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
int i,j,chk;
double r,theta,psi;
double omegac;

Array<double,1> fix(Range(1,2));

// Period of fixed point
int periode;

// Input filenames
LA_STRING basename;
LA_STRING praefix = "";
if(argc==4) praefix = "_" + LA_STRING(argv[3]);
if(argc>=3) 
{
	periode = atoi(argv[2]);
	basename = LA_STRING(argv[1]); 
}
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// log file
ofs2.open("log_" + LA_STRING(program_name) + "_" + LA_STRING(periode) + praefix + ".dat");
ofs2.precision(16);

// Search grid
vector<double> startvec;
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);

double xmin = 0;		//theta
double xmax = pi2;
int Nx = 20;

double ymin = 0.4;	//r
double ymax = 0.7;
int Ny = 20;

double phistart = 0;
int MapDirection = 1;

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

if(useparfile==1)
{
	cout << "All parameters are read from file" << endl;
	ofs2 << "All parameters are read from file" << endl;
	ymin = startvec[2];
	ymax = startvec[3];
	xmin = startvec[4];
	xmax = startvec[5];
	Nx = Ny = int(sqrt(startvec[6]));
	phistart = startvec[7];
	Ekin = startvec[18];
	lambda = startvec[19];
}
const double dx = (xmax-xmin)/double(Nx-1);
const double dy = (ymax-ymin)/double(Ny-1);

// additional parameters for IO
const int psize = 10;
parstruct * parvec = new parstruct[psize];
parvec[0].name = "r-grid";			parvec[0].wert = Ny;
parvec[1].name = "theta-grid";		parvec[1].wert = Nx;
parvec[2].name = "rmin";			parvec[2].wert = ymin;
parvec[3].name = "rmax";			parvec[3].wert = ymax;
parvec[4].name = "thmin";			parvec[4].wert = xmin;
parvec[5].name = "thmax";			parvec[5].wert = xmax;
parvec[6].name = "phistart";		parvec[6].wert = phistart;
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
LA_STRING filenameout = "fix_" + LA_STRING(periode) + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(4);
var[0] = "theta[rad]";  var[1] = "r[m]";  var[2] = "period";  var[3] = "psi";
writeiodata(out,var,parvec,psize,parfilename);

// Extend Boundary close to efit boundary 0.84 2.54 -1.6 1.6
bndy[0] = 0.88;  bndy[1] = 2.5;  bndy[2] = -1.56;  bndy[3] = 1.56;	

ofs2 << Ny << " rows, done:" << endl;
for(i=0;i<Ny;i++)	//r
{
	r = ymin + i*dy;
	for(j=0;j<Nx;j++)	//theta
	{
		theta = xmin + j*dx;
		fix(1) = theta;
		fix(2) = r;

		chk = newton2D(fix,phistart,psi,periode);
		if(chk==-1) continue;

		if(periode==1)
		{
			if(fix(2)>1)
			{
				out << fix(1) << "\t" << fix(2) << "\t" << periode << "\t" << psi << endl;
				ofs2 << "Program terminated normally" << endl;
				return 0; 
			}
			else continue;
		}
		else out << fix(1) << "\t" << fix(2) << "\t" << periode << "\t" << psi << endl;
	}
	ofs2 << i+1 << "\t" << flush;
}
ofs2 << endl;
ofs2 << "Program terminated normally" << endl;
return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- check_boundary ---------------------------------------------------------
int check_boundary(Array<double,1> x)
{
double R,Z;
R = EQD.RmAxis + x(2)*cos(x(1));
Z = EQD.ZmAxis + x(2)*sin(x(1));

if(R<bndy[0] || R>bndy[1] || Z<bndy[2] || Z>bndy[3]) 
{
	ofs2 << "Boundary crossed" << endl;
	return -1;
}
//else cout << R << "\t" << Z << endl;
return 0;
}

//----------- newton2dim ------------------------------------------------
// all variables with 'r' represent theta and all with 't' represent r!!!
// x(1) = theta, x(2) = r
// J(1) = dth(i+1)/dth(i),  J(2) = dth(i+1)/dr(i),  J(3) = dr(i+1)/dth(i),  J(4) = dr(i+1)/dr(i)
int newton2D(Array<double,1>& x, double phistart, double& psi, int periode)	//0: ok		-1: Fehler
{
double fr,ft,dr,dt,det,length;
int i,chk;
double phi;

const int imax = 100;
const double delta = 1e-12;

// Vectors
Array<double,1> J(Range(1,4));		
Array<double,1> xn(Range(1,2));

// Search
xn = x;

for(i=0;i<=imax;i++)
{
	phi = phistart;
	chk = check_boundary(xn);
	if(chk==-1) return -1;
	chk = mapit_J(xn,phi,psi,J,periode);
	if(chk<0){ofs2 << "No convergency " << chk << endl; return -1;}

	
	fr = xn(1) - x(1);
	ft = xn(2) - x(2);

	det = (J(1)-1)*(J(4)-1) - J(2)*J(3);
	dr = ((J(4)-1)*fr-J(2)*ft)/det;
	dt = ((J(1)-1)*ft-J(3)*fr)/det;

	length = sqrt(dr*dr + dt*dt);
	//if(i%20==0){cout << dr << "\t" << dt << "\t" << length <<  endl;	getchar();}
	if(length<delta) 
	{
		return 0;	// convergency
	} 

	x(1) -= dr;
	x(2) -= dt;
	//x(1) = modulo2pi(x(1));
	xn = x;
}

ofs2 << "No convergency " <<  xn(1) << "\t" << xn(2) << "\t" << dr << "\t" << dt << "\t" << length << endl;
return -1;
}

//--------- mapit_J -----------------------------------------------------------------
// tracer needs angle in degrees! --> transform theta!!!!!
int mapit_J(Array<double,1>& xa, double& phi, double& psi, Array<double,1>& J, int itt, int Map)
{
int chk;
double r,theta;
double ralt,thetaalt,phialt;
const double dr = 0.00001;
const double dth = 0.0001;	// in Rad

Array<double,1> r_stencil(Range(1,4)),th_stencil(Range(1,4));	// 1:r+dr 2:r-dr 3:th+dth 4:th-dth

thetaalt = xa(1);	//in Rad
ralt = xa(2);
phialt = phi;

//up: r+dr
r = ralt + dr;
theta = thetaalt;
phi = phialt;

chk = mapit(theta,r,phi,psi,itt,Map);
if(chk<0) return -1;
th_stencil(1) = theta;
r_stencil(1) = r;

//low: r-dr
r = ralt - dr;
theta = thetaalt;
phi = phialt;

chk = mapit(theta,r,phi,psi,itt,Map);
if(chk<0) return -1;
th_stencil(2) = theta;
r_stencil(2) = r;

//right: th+dth
r = ralt;
theta = thetaalt + dth;
phi = phialt;

chk = mapit(theta,r,phi,psi,itt,Map);
if(chk<0) return -1;
th_stencil(3) = theta;
r_stencil(3) = r;

//left: th-dth
r = ralt;
theta = thetaalt - dth;
phi = phialt;

chk = mapit(theta,r,phi,psi,itt,Map);
if(chk<0) return -1;
th_stencil(4) = theta;
r_stencil(4) = r;

// derivatives
// J(1)=dth(i+1)/dth(i) J(2)=dth(i+1)/dr(i) J(3)=dr(i+1)/dth(i) J(4)=dr(i+1)/dr(i)
J(1) = 0.5*(th_stencil(3)-th_stencil(4))/dth;
J(2) = 0.5*(th_stencil(1)-th_stencil(2))/dr;
J(3) = 0.5*(r_stencil(3)-r_stencil(4))/dth;
J(4) = 0.5*(r_stencil(1)-r_stencil(2))/dr;

// center
r = ralt;
theta = thetaalt;
phi = phialt;

chk = mapit(theta,r,phi,psi,itt,Map);
if(chk<0) return -1;

xa(1) = theta;
xa(2) = r;

return 0;
}

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



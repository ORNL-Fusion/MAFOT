// Program calculates Poincaré Plot for ITER-Drift with Time dependent perturbations
// Fortran Subroutines for Perturbation are used
// A.Wingen						5.08.10

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	Poincaré particel drift data file
//			log-file



// Define
//--------
#define program_name "iterplot"
//#define BZ_DEBUG		// Debug Blitz-Arrays

// Include
//--------
#include <andi.hxx>
#include <efit_class_iter.hxx>
#include <iter-drift.hxx>

// Prototypes

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
int i,n,chk;
double r,theta,phi,psi;
double omegac;

// Arrays
Array<double,1> y(nvar); // Array to set initial conditions

// Use system time as seed(=idum) for random numbers
double now=zeit();
long idum=long(now);

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

// Read parameter file
vector<double> startvec;
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);


// Set starting parameters
int itt = 100;
double rmin = 1.0;
double rmax = 1.9;
double thmin = 0;
double thmax = 0;
int N = 31;
double phistart = 0;	// in deg
int MapDirection = 1;

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

if(useparfile==1)
{
	cout << "All parameters are read from file" << endl << endl;
	ofs2 << "All parameters are read from file" << endl << endl;
	itt = int(startvec[1]);
	rmin = startvec[2];
	rmax = startvec[3];
	thmin = startvec[4];
	thmax = startvec[5];
	N = int(startvec[6]);
	phistart = startvec[7];
	MapDirection = int(startvec[8]);
	Ekin = startvec[16];
	lambda = startvec[17];
}

// additional parameters for IO
const int psize = 10;
parstruct * parvec = new parstruct[psize];
parvec[0].name = "Max. Iterations";	parvec[0].wert = itt;
parvec[1].name = "Points";			parvec[1].wert = N;
parvec[2].name = "rmin";			parvec[2].wert = rmin;
parvec[3].name = "rmax";			parvec[3].wert = rmax;
parvec[4].name = "thmin";			parvec[4].wert = thmin;
parvec[5].name = "thmax";			parvec[5].wert = thmax;
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
LA_STRING filenameout = "plot" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(6);
var[0] = "theta[rad]"; var[1] = "r[m]"; var[2] = "phi[deg]"; var[3] = "psi"; var[4] = "R[m]"; var[5] = "Z[m]";
writeiodata(out,var,parvec,psize,parfilename);

cout << "Create (0=set, 1=random, 2=target): " << create_flag << endl; 
if(create_flag==2) cout << "Target (0=vertical, 1=45°, 2=horizontal): " << which_target_plate << endl;
cout << endl << "Start Tracer for " << N << " points ... " << endl;
ofs2 << "Create (0=grid, 1=random, 2=target): " << create_flag << endl; 
if(create_flag==2) ofs2 << "Target (0=vertical, 1=45°, 2=horizontal): " << which_target_plate << endl;
ofs2 << endl << "Start Tracer for " << N << " points ... " << endl;

// Get Poincaré section
for(n=1;n<=N;n++)
{
	// Set initial values y(phistart)
	switch(create_flag) 
	{
		//case 2:
		//	start_on_target(n,N,1,rmin,rmax,phistart,phistart,R_target,Z_target,S_target,r,theta,phi);
		//	break;
		case 1:
			create(idum,rmin,rmax,thmin,thmax,r,theta);
			break;
		default:
			set(n,N,rmin,rmax,thmin,thmax,r,theta);
			break;
	}
	phi = phistart;
	y(0) = r*cos(theta) + EQD.RmAxis;	// R
	y(1) = r*sin(theta) + EQD.ZmAxis;	// Z

	// Integrate
	for(i=1;i<=itt;i++)
	{
		chk = mapstep(y,theta,r,phi,psi,MapDirection);
		if(chk<0) {ofs2 << "mapit: wall hit" << endl; break;}	// particle has left system

		if(fabs(phi-MapDirection*i*dpinit*ilt-phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(phi-MapDirection*i*dpinit*ilt-phistart) << endl;
		phi=MapDirection*i*dpinit*ilt+phistart;

		// Output
		out << theta << "\t" << r << "\t" << phi << "\t" << psi << "\t" << y(0) << "\t" << y(1) << endl;
	} // end for i
	//if(n%50==0) cout << "Trax: " << n << "\t" << "Steps: " << i-1 << endl;
	ofs2 << "Trax: " << n << "\t" << "Steps: " << i-1 << endl;
} // end for n

double now2=zeit();
cout << endl;
cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;
return 0;
} //end main

//----------------------- End of Main -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

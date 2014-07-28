// Subroutines of MAFOT programs for D3D-Drift with Time dependent perturbations
// Fortran Subroutines for Perturbation are used
//
//
//void dtplot (double R[], double Z[], double Rout[], double Zout[], double phiout[], double psiout[],
//			 double phistart, double phiend, int N, int shot, int time,
//			 int quiet = 0, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1, int response = 0,
//			 int sigma = 0, int charge = -1, double Ekin = 0, double lambda = 0.1, int useFilament = 0, int write2file = 1,
//			 char* supPath = "./", char* gPath = gFilePath, char* outPath = "./", char* char_praefix = "");
// A.Wingen						8.10.13
//
// Input:
//			Arrays: R in [m], Z in [m] 			initial conditions
//			phistart in [deg] 					initial toroidal angle
//			double phiend 						toroidal angle for end of integration, in [deg],
//			int N								number of elements in arrays
//			int shot, time						which discharge and time, values passed as integer!
//		Flags: (integer)
//  		quiet								1: no infos written to stdout, 0: several stdout infos (default)
//			useFcoil, useCcoil, useIcoil		0: off;	1: on (default)
//			useM3DC1 							-1: g-file vacuum (default); 0: equilibrium only; 1: perturbation only; 2: total field
//			response							0: no (default), 1: yes
//			int sigma							0: field-lines (default), 1: passing particles, -1: counter-passing
//			int charge							-1: electrons (default), >=1: ions with this charge number
//			double Ekin							kinetic particle Energy in keV
//			double lambda						perp. to parallel energy ratio (default = 0.1)
//			int useFilament						0: off (default), >=1: Number of current filaments to include
//			int write2file						save output to file 1: yes (default), 0: no
//		Paths: (character pointer)
//			supPath								path to diiidsup.in file (default is the current working directory)
//			gPath								path to g-file (default is gFilePath, which is defined in efit_class.hxx)
//			outPath								path to write output file (default is the current working directory)
//		File Name tag: (character pointer)
//			char_praefix						arbitrary tag which is added to the output file name, default is an empty string
// Output:
//			Arrays: Rout, Zout, phiout and psiout
//				output occurs every 360 deg toroidally
//				IMPORTANT: initial array size must be sufficient to hold N*itt values
//
//--------------------------------------------------------------------------------------------------------------------------
// follow_field_lines (double R[], double Z[], double phi[], int N, double phiend, int shot, int time,
//					   int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
//					   char* supPath = "./", char* gPath = gFilePath)
// A.Wingen						21.6.12
//
// Input:
//			Arrays: R in [m], Z in [m], phi in [deg] 	initial conditions
//			int N										number of elements in arrays
//			double phiend 								toroidal angle for end of integration, in [deg],
//			int shot, time								which discharge and time, values passed as integer!
//		Flags: (integer)
//  		quiet										1: no infos written to stdout (default);		0: several stdout infos
//			useFcoil, useCcoil, useIcoil				0: off;	1: on (default)
//			useM3DC1									-1: g-file vacuum (default); 0: equilibrium only; 1: perturbation only; 2: total field
//		Paths: (character pointer)
//			supPath										path to diiidsup.in file (default is the current working directory)
//			gPath										path to g-file (default is gFilePath, which is defined in efit_class.hxx)
// Output:
//			Arrays: R, Z, phi contain the final position of each initial condition at phiend (check: phi = phiend should be true for all elements)
//			initial arrays are overwritten
//
//--------------------------------------------------------------------------------------------------------------------------
// get_B_field (double R[], double Z[], double phi[], int N, int shot, int time,
//					   int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
//					   char* supPath = "./", char* gPath = gFilePath)
// A.Wingen						27.7.13
//
// Input:
//			Arrays: R in [m], Z in [m], phi in [deg] 	initial conditions
//			int N										number of elements in arrays
//			int shot, time								which discharge and time, values passed as integer!
//		Flags: (integer)
//  		quiet										1: no infos written to stdout (default);		0: several stdout infos
//			useFcoil, useCcoil, useIcoil				0: off;	1: on (default)
//			useM3DC1									-1: g-file vacuum (default); 0: equilibrium only; 1: perturbation only; 2: total field
//		Paths: (character pointer)
//			supPath										path to diiidsup.in file (default is the current working directory)
//			gPath										path to g-file (default is gFilePath, which is defined in efit_class.hxx)
// Output:
//			Arrays: R, Z, phi contain the magnetic field components B_R, B_Z and B_phi respectively
//			initial arrays are overwritten


// Define
//--------
#define program_name "mafot"
//#define BZ_DEBUG		// Debug Blitz-Arrays

// Include
//--------
#include <mafot.hxx>
#ifdef m3dc1
#include <d3d_m3dc1.hxx>
#else
#include <d3d.hxx>
#endif

// Prototypes
extern "C"	// this makes the function name in the dylib to be as it is here, otherwise the name is mangled (contains lots of strange characters, standard c++ style)
{
void dtplot (double R[], double Z[], double Rout[], double Zout[], double phiout[], double psiout[],
			 double phistart, double phiend, int N, int shot, int time,
			 int quiet = 0, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1, int response = 0,
			 int sigma = 0, int charge = -1, double Ekin = 0, double lambda = 0.1, int useFilament = 0, int write2file = 1,
			 char* supPath = "./", char* gPath = gFilePath, char* outPath = "./", char* char_praefix = "");

void follow_field_lines (double R[], double Z[], double phi[], int N, double phiend, int shot, int time,
						 int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
						 char* supPath = "./", char* gPath = gFilePath);

void get_B_field (double R[], double Z[], double phi[], int N, int shot, int time,
				  int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
				  char* supPath = "./", char* gPath = gFilePath);
}

// Switches

// Golbal Parameters



//-------------------------------------------------------------------------------------------------------------------------
// dtplot Program  ->  dtplot.cxx as a callable subroutine
//-------------------------------------------------------------------------------------------------------------------------
void dtplot (double R[], double Z[], double Rout[], double Zout[], double phiout[], double psiout[],
			 double phistart, double phiend, int N, int shot, int time,
			 int quiet, int useFcoil, int useCcoil, int useIcoil, int useM3DC1, int response,
			 int sigma, int charge, double Ekin, double lambda, int useFilament, int write2file,
			 char* supPath, char* gPath, char* outPath, char* char_praefix)
{
// Variables
int i,n,chk,index;
EFIT EQD;
IO PAR(EQD);
LA_STRING praefix;
if(strlen(char_praefix) > 0) praefix = "_" + LA_STRING(char_praefix);

// Use system time as seed(=idum) for random numbers
double now = zeit();

// log file
if(write2file == 1) ofs2.open(LA_STRING(outPath) + "log_dtplot" + praefix + ".dat");

// Parameters
EQD.Shot = LA_STRING(shot);
EQD.Time = LA_STRING(time);
EQD.Path = LA_STRING(gPath);

// R grid
PAR.N = N;
PAR.NR = 1;
PAR.Rmin = 0;
PAR.Rmax = 0;

// Z grid
PAR.NZ = 1;
PAR.Zmin = 0;
PAR.Zmax = 0;

// Mapping
PAR.phistart = phistart;
PAR.MapDirection = sign(phiend - PAR.phistart);
PAR.itt = abs(int(round((phiend - PAR.phistart)/360)));

// Particle Parameters
PAR.Ekin = Ekin;
PAR.lambda = lambda;
PAR.verschieb = 0;

// Set switches
PAR.which_target_plate = 0;
PAR.create_flag = 0;
PAR.useFcoil = useFcoil;
PAR.useCcoil = useCcoil;
PAR.useIcoil = useIcoil;
PAR.sigma = sigma;
PAR.Zq = charge;
PAR.useFilament = useFilament;
PAR.useTprofile = 0;
PAR.response = response;
PAR.response_field = useM3DC1;	// -1: use vacuum field from g-file (= M3D-C1 off, here: default); 0: equilibrium only; 1: perturbation only; 2: total field
if(quiet == 0) cout << "use M3D-C1: " << PAR.response_field << endl;

// additional parameters for IO
PAR.set_pv(6);
PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
PAR.pv[1].name = "Points";			PAR.pv[1].wert = PAR.N;
PAR.pv[2].name = "phistart";		PAR.pv[2].wert = PAR.phistart;
PAR.pv[3].name = "MapDirection";	PAR.pv[3].wert = PAR.MapDirection;
PAR.pv[4].name = "Ekin";			PAR.pv[4].wert = PAR.Ekin;
PAR.pv[5].name = "energy ratio lambda";	PAR.pv[5].wert = PAR.lambda;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
if(quiet == 0) cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
if(write2file == 1) ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prep_perturbation(EQD,PAR,quiet,supPath);

// Prepare particles
PARTICLE FLT(EQD,PAR,quiet);

// Output
ofstream out;
if(write2file == 1)
{
	LA_STRING filenameout = LA_STRING(outPath) + "plot" + praefix + ".dat";
	outputtest(filenameout);
	out.open(filenameout);
	out.precision(16);
	vector<LA_STRING> var(6);
	var[0] = "theta[rad]"; var[1] = "r[m]"; var[2] = "phi[deg]"; var[3] = "psi"; var[4] = "R[m]"; var[5] = "Z[m]";
	PAR.writeiodata(out,bndy,var);

	ofs2 << "MapDirection: " << PAR.MapDirection << "\t" << "toroidal turns: " << PAR.itt << endl;
	ofs2 << "Helicity: " << EQD.helicity << endl;
	ofs2 << endl << "Start Tracer for " << PAR.N << " points ... " << endl;
}

if(quiet == 0)
{
	cout << "MapDirection: " << PAR.MapDirection << "\t" << "toroidal turns: " << PAR.itt << endl;
	cout << "Helicity = " << EQD.helicity << endl;
	cout << endl << "Start Tracer for " << PAR.N << " points ... " << endl;
}

// Follow all orbits/field lines
index = 0;
for(n=0;n<PAR.N;n++)
{
	// Set initial values in FLT
	FLT.R = R[n];
	FLT.Z = Z[n];
	FLT.phi = PAR.phistart;

	// Integrate
	for(i=1;i<=PAR.itt;i++)
	{
		chk = FLT.mapstep(PAR.MapDirection);
		if(chk < 0 && write2file == 1) {ofs2 << "mapit: wall hit" << endl; break;}	// particle has left system

		if(fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) > 1e-10 && quiet == 0) cout << "MAFOT Warning: wrong toroidal angle" << endl;
		FLT.phi = PAR.MapDirection*i*dpinit*ilt + PAR.phistart;

		// Output
		Rout[index] = FLT.R;
		Zout[index] = FLT.Z;
		phiout[index] = FLT.phi;
		psiout[index] = FLT.psi;
		index += 1;

		if(write2file == 1) out << FLT.get_theta() << "\t" << FLT.get_r() << "\t" << FLT.phi << "\t" << FLT.psi << "\t" << FLT.R << "\t" << FLT.Z << endl;

	} // end for i
	if(write2file == 1) ofs2 << "Trax: " << n << "\t" << "Steps: " << i-1 << endl;
} // end for n

double now2 = zeit();
if(quiet == 0) cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
if(write2file == 1) ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;
ofs2.close();

#ifdef m3dc1
if(PAR.response_field >= 0) m3dc1_unload_file_();
#endif

return;
} //end dtplot
//-------------------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------------------
// follow_field_lines Program  ->  wrapper around mapit from particle_class.hxx
//-------------------------------------------------------------------------------------------------------------------------
void follow_field_lines (double R[], double Z[], double phi[], int N, double phiend, int shot, int time,
						 int quiet, int useFcoil, int useCcoil, int useIcoil, int useM3DC1,
						 char* supPath, char* gPath)
{
// Variables
int n,chk;
EFIT EQD;
IO PAR(EQD);

// Parameters
EQD.Shot = LA_STRING(shot);
EQD.Time = LA_STRING(time);
EQD.Path = LA_STRING(gPath);

if(quiet == 0)
{
	cout << N << "\t" << R[0] << "\t" << R[1] << endl;
	cout << phi[0] << "\t" << phiend << endl;
	cout << shot << "\t" << time << endl;
}

// R grid
PAR.N = N;
PAR.NR = 1;
PAR.Rmin = 0;
PAR.Rmax = 0;

// Z grid
PAR.NZ = 1;
PAR.Zmin = 0;
PAR.Zmax = 0;

// Particle Parameters
PAR.Ekin = 0;
PAR.lambda = 0;
PAR.verschieb = 0;

// Set switches
PAR.which_target_plate = 0;
PAR.create_flag = 0;
PAR.useFcoil = useFcoil;
PAR.useCcoil = useCcoil;
PAR.useIcoil = useIcoil;
PAR.sigma = 0;
PAR.Zq = 1;
PAR.useFilament = 0;
PAR.useTprofile = 0;
PAR.response = 1;
PAR.response_field = useM3DC1;	// -1: use vacuum field from g-file (= M3D-C1 off, here: default); 0: equilibrium only; 1: perturbation only; 2: total field
if(quiet == 0) cout << "use M3D-C1: " << PAR.response_field << endl;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);

// Prepare Perturbation
prep_perturbation(EQD,PAR,quiet,supPath);	// quiet = 1 causes function to be quiet

// Prepare particles
PARTICLE FLT(EQD,PAR,quiet);		// quiet = 1 causes function to be quiet

// Follow all field lines
for(n=0;n<PAR.N;n++)
{
	// Set initial values in FLT
	FLT.R = R[n];
	FLT.Z = Z[n];
	FLT.phi = phi[n];

	// Map Parameters
	PAR.phistart = phi[n];
	PAR.MapDirection = sign(phiend - PAR.phistart);
	PAR.itt = abs(int(round(phiend - PAR.phistart)));
	if(quiet == 0) cout << PAR.MapDirection << "\t" << PAR.itt << endl;


	// Integrate
	chk = FLT.mapstep(PAR.MapDirection, PAR.itt);
	if(chk == -1 && quiet == 0) cout << "MAFOT: field line crossed wall" << endl;

	// Output
	R[n] = FLT.R;
	Z[n] = FLT.Z;
	phi[n] = FLT.phi;
	if(int(round(FLT.phi - phiend)) != 0 && chk == 0 && quiet == 0) cout << "MAFOT Warning: wrong toroidal angle" << endl;
} // end for n

#ifdef m3dc1
if(PAR.response_field >= 0) m3dc1_unload_file_();
#endif

return;
} //end follow_field_lines
//-------------------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------------------
// get_B_field Program  ->  wrapper around getBfield from d3d.hxx
//-------------------------------------------------------------------------------------------------------------------------
void get_B_field (double R[], double Z[], double phi[], int N, int shot, int time,
				  int quiet, int useFcoil, int useCcoil, int useIcoil, int useM3DC1,
				  char* supPath, char* gPath)
{
// Variables
int n;
EFIT EQD;
IO PAR(EQD);

// Parameters
EQD.Shot = LA_STRING(shot);
EQD.Time = LA_STRING(time);
EQD.Path = LA_STRING(gPath);

if(quiet == 0)
{
	cout << N << "\t" << R[0] << "\t" << R[1] << endl;
	cout << shot << "\t" << time << endl;
}

// R grid
PAR.N = N;
PAR.NR = 1;
PAR.Rmin = 0;
PAR.Rmax = 0;

// Z grid
PAR.NZ = 1;
PAR.Zmin = 0;
PAR.Zmax = 0;

// Particle Parameters
PAR.Ekin = 0;
PAR.lambda = 0;
PAR.verschieb = 0;

// Set switches
PAR.which_target_plate = 0;
PAR.create_flag = 0;
PAR.useFcoil = useFcoil;
PAR.useCcoil = useCcoil;
PAR.useIcoil = useIcoil;
PAR.sigma = 0;
PAR.Zq = 1;
PAR.useFilament = 0;
PAR.useTprofile = 0;
PAR.response = 1;
PAR.response_field = useM3DC1;	// -1: use vacuum field from g-file (= M3D-C1 off, here: default); 0: equilibrium only; 1: perturbation only; 2: total field
if(quiet == 0) cout << "use M3D-C1: " << PAR.response_field << endl;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);

// Prepare Perturbation
prep_perturbation(EQD,PAR,quiet,supPath);	// quiet = 1 causes function to be quiet

double B_R, B_Z, B_phi;

// Follow all field lines
for(n=0;n<PAR.N;n++)
{
	// Get magnetic field
	getBfield(R[n], Z[n], phi[n], B_R, B_Z, B_phi, EQD, PAR);

	// Output
	R[n] = B_R;
	Z[n] = B_Z;
	phi[n] = B_phi;
} // end for n

#ifdef m3dc1
if(PAR.response_field >= 0) m3dc1_unload_file_();
#endif

return;
} //end get_B_field
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// Subroutine to follow field lines for D3D-Drift with Time dependent perturbations
// Fortran Subroutines for Perturbation are used
// A.Wingen						21.6.12

// follow_field_lines (double R[], double Z[], double phi[], int N, double phiend, int shot, int time,
//					   int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
//					   char* supPath = "./", char* gPath = gFilePath)
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



// Define
//--------
#define program_name "follow_field_lines"
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
void follow_field_lines (double R[], double Z[], double phi[], int N, double phiend, int shot, int time,
						 int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
						 char* supPath = "./", char* gPath = gFilePath);
}

// Switches

// Golbal Parameters

// Main Program
//--------------
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
} //end main

//----------------------- End of Main -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

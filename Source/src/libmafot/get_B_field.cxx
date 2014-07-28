// Subroutine to return the total magnetic field for D3D-Drift with Time dependent perturbations
// Fortran Subroutines for Perturbation are used
// A.Wingen						26.7.13

// get_B_field (double R[], double Z[], double phi[], int N, int shot, int time,
//					   int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
//					   char* supPath = "./", char* gPath = gFilePath)
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
//			Arrays: R, Z, phi contain the total field B_R, B_Z and B_phi
//			initial arrays are overwritten



// Define
//--------
#define program_name "get_B_field"
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
void get_B_field (double R[], double Z[], double phi[], int N, int shot, int time,
						 int quiet = 1, int useFcoil = 1, int useCcoil = 1, int useIcoil = 1, int useM3DC1 = -1,
						 char* supPath = "./", char* gPath = gFilePath);
}

// Switches

// Golbal Parameters

// Main Program
//--------------
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
} //end main

//----------------------- End of Main -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

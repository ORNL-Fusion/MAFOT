// IO Class contains all control parameters read from file
// needs input variables from EFIT class
// member functions to read and write parameters are defined in Machine Header file
// used by all D3D, ITER and NSTX drift programs
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						7.06.11


// Define
//--------
#ifndef IO_CLASS_INCLUDED
#define IO_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <efit_class.hxx>
using namespace blitz;

// Prototypes  

// Typedef
typedef struct {string name; double wert;} parstruct;

// Golbal Parameters 

//--------- Begin Class PARTICLE ----------------------------------------------------------------------------------------------
class IO
{
private:
// Member Variables
	EFIT& EQDr;			// Only a Reference, not a copy
	LA_STRING filename;	// Name of Parameter file
	int psize;			// length of parstruct array

public: 
// Member Variables
	int itt;									// toroidal iterations
	double phistart;							// toroidal start-angle
	int MapDirection;							// +1 forward, -1 backwards, 0 both
	int Nt; double tmin; double tmax;			// t grid (footprints only)
	int Nphi; double phimin; double phimax;		// phi grid (footprints only)
	int NR; double Rmin; double Rmax;			// R grid (laminar only)
	int NZ; double Zmin; double Zmax;			// Z grid (laminar only)
	int N; int Nr; double rmin; double rmax;	// r grid
	int Nth; double thmin; double thmax;		// theta grid
	int which_target_plate;						// 0: center post(above inner target), 1: inner target(45deg), 2: outer target(horizontal), 3: shelf (right of outer target)
	int create_flag;							// 0: fixed grid	1: random numbers	2: Start on target
	int useIcoil; int useCcoil; int useFcoil;	// Perturbation coils:  0: off, 1: on
	int useFilament;							// Current filaments: 0: off, >=1: Number of filaments to use
	int Zq;										// Charge number: 1: ions are calculated	-1: electrons are calculated
	int sigma;									// 1: co-passing particles		-1: count-passing particles		0: field lines only
	int useTprofile;							// Temperature profile: 0: constant Energy, 1: according to Position
	double Ekin;								// Kinetic Energy in keV
	double lambda;								// fraction of Energy in radial direction
	double verschieb;							// parameter for mainfolds
	int response;								// 1: use M3D-C1 plasma response solution; 0: no plasma response included -> vacuum field as used in M3D-C1
	int response_field;							// -1: use vacuum field from g-file (= M3D-C1 off); 0: equilibrium only; 1: perturbation only; 2: total field
	int output_step_size;

	parstruct* pv;								// parstruct array

// Constructors
	IO(EFIT& EQD);												// Default Constructor
	IO(EFIT& EQD, LA_STRING name, int size=0, int mpi_rank=0);	// Constructor with readiodata

// Member-Operators
	IO& operator =(const IO& PAR);											// Operator =
	friend std::ostream& operator <<(std::ostream& out, const IO& PAR);		// Operator <<

// Member-Functions
	void readiodata(char* name, int mpi_rank=0);
	void writeiodata(ofstream& out, double bndy[], vector<LA_STRING>& var);
	void set_pv(int size);

}; //end of class

//------------------------ End of Class -----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Constructor list necessary to set EQDr since references cannot be empty
IO::IO(EFIT& EQD): EQDr(EQD)										// Default
{
	filename = "";
	psize = 0;
	pv = 0;
	output_step_size = 360;
}					
IO::IO(EFIT& EQD, LA_STRING name, int size, int mpi_rank): EQDr(EQD)	// with readiodata
{
	psize = size;
	readiodata(name, mpi_rank);
	if(psize > 0) pv = new parstruct[psize];
	else pv = 0;
	output_step_size = 360;
}	

//--------- Operator = ----------------------------------------------------------------------------------------------------
IO& IO::operator =(const IO& PAR)
{
if (this == &PAR) return(*this);	    // if: x=x
// Private Member Variables
// Reference to EQD remains unchanged! There is only one EQD
filename = PAR.filename;
psize = PAR.psize;

// Map Parameters
itt = PAR.itt;
phistart = PAR.phistart;
MapDirection = PAR.MapDirection;

// t grid (footprints only)
Nt = PAR.Nt; 
tmin = PAR.tmin; 
tmax = PAR.tmax;			

// phi grid (footprints only)
Nphi = PAR.Nphi; 
phimin = PAR.phimin; 
phimax = PAR.phimax;		

// R grid (laminar only)
NR = PAR.NR; 
Rmin = PAR.Rmin; 
Rmax = PAR.Rmax;		

// Z grid (laminar only)
NZ = PAR.NZ; 
Zmin = PAR.Zmin; 
Zmax = PAR.Zmax;		

// r grid
N = PAR.N; 
Nr = PAR.Nr; 
rmin = PAR.rmin; 
rmax = PAR.rmax;		

// theta grid
Nth = PAR.Nth; 
thmin = PAR.thmin; 
thmax = PAR.thmax;		

// Particle Parameters
Ekin = PAR.Ekin;
lambda = PAR.lambda;
verschieb = PAR.verschieb;

// Set switches
which_target_plate = PAR.which_target_plate;
create_flag = PAR.create_flag;
useFcoil = PAR.useFcoil;
useCcoil = PAR.useCcoil;
useIcoil = PAR.useIcoil;
sigma = PAR.sigma;
Zq = PAR.Zq;
useFilament = PAR.useFilament;
useTprofile = PAR.useTprofile;

// M3D-C1 parameter
response = PAR.response;
response_field = PAR.response_field;

// Parstruct array
if(pv) delete[] pv;
pv = new parstruct[psize];
for(int i=0;i<psize;i++) 
{
	pv[i].name = PAR.pv[i].name;
	pv[i].wert = PAR.pv[i].wert;
}

return(*this);
}

//--------- Operator << ---------------------------------------------------------------------------------------------------
ostream& operator <<(ostream& out, const IO& PAR)
{
out << "--- Data from file: " << PAR.filename << " ----------------------" << endl;
out << "--- Map Parameters ---" << endl;
out << "itt = " << PAR.itt << endl;
out << "phistart = " << PAR.phistart << endl;
out << "MapDirection = " << PAR.MapDirection << endl;
out << endl;

out << "--- t Grid (footprints only) ---" << endl;
out << "Nt = " << PAR.Nt << endl; 
out << "tmin = " << PAR.tmin << endl; 
out << "tmax = " << PAR.tmax << endl;			
out << endl;

out << "--- phi Grid (footprints only) ---" << endl;
out << "Nphi = " << PAR.Nphi << endl; 
out << "phimin = " << PAR.phimin << endl; 
out << "phimax = " << PAR.phimax << endl;		
out << endl;

out << "--- R Grid (laminar only) ---" << endl;
out << "NR = " << PAR.NR << endl; 
out << "Rmin = " << PAR.Rmin << endl; 
out << "Rmax = " << PAR.Rmax << endl;		
out << endl;

out << "--- Z Grid (laminar only) ---" << endl;
out << "NZ = " << PAR.NZ << endl; 
out << "Zmin = " << PAR.Zmin << endl; 
out << "Zmax = " << PAR.Zmax << endl;		
out << endl;

out << "--- r Grid ---" << endl;
out << "N = " << PAR.N << endl; 
out << "Nr = " << PAR.Nr << endl; 
out << "rmin = " << PAR.rmin << endl; 
out << "rmax = " << PAR.rmax << endl;		
out << endl;

out << "--- theta Grid ---" << endl;
out << "Nth = " << PAR.Nth << endl; 
out << "thmin = " << PAR.thmin << endl; 
out << "thmax = " << PAR.thmax << endl;		
out << endl;

out << "--- Particle Parameters ---" << endl;
out << "Ekin = " << PAR.Ekin << endl;
out << "lambda = " << PAR.lambda << endl;
out << "verschieb = " << PAR.verschieb << endl;
out << endl;

out << "--- Switches ---" << endl;
out << "which_target_plate = " << PAR.which_target_plate << endl;
out << "create_flag = " << PAR.create_flag << endl;
out << "useFcoil = " << PAR.useFcoil << endl;
out << "useCcoil = " << PAR.useCcoil << endl;
out << "useIcoil = " << PAR.useIcoil << endl;
out << "sigma = " << PAR.sigma << endl;
out << "Zq = " << PAR.Zq << endl;
out << "useFilament = " << PAR.useFilament << endl;
out << "useTprofile = " << PAR.useTprofile << endl;
out << endl;

#ifdef m3dc1
out << "--- M3D-C1 ---" << endl;
out << "response = " << PAR.response << endl;
out << "response_field = " << PAR.response_field << endl;
out << endl;
#endif

if(PAR.psize>0) 
{
	out << "--- Parstruct array ---" << endl;
	for(int i=0;i<PAR.psize;i++) out << PAR.pv[i].name << " = " << PAR.pv[i].wert << endl;
}
return out;
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// readiodata and writeiodata are defined in Machine specific header file

// ------------------- set_pv ---------------------------------------------------------------------------------------------
void IO::set_pv(int size)
{
psize = size;
pv = new parstruct[psize];
}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  IO_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

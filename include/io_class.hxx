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
void readparfile_new(char* name, vector<double>& vec, int& shot, int& time, LA_STRING& path);

// Typedef
typedef struct {string name; double wert;} parstruct;

// Golbal Parameters 
extern int simpleBndy;
extern const double dpinit;

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
	int useBuswork;								// Bus work error field:  0: off, 1: on
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
useBuswork = PAR.useBuswork;
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
out << "useBuswork = " << PAR.useBuswork << endl;
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

// ----------------- readiodata -------------------------------------------------------------------------------------------
void IO::readiodata(char* name, int mpi_rank)
{
// Get Parameters
vector<double> vec;
int shot,time;
readparfile_new(name,vec,shot,time,EQDr.Path);
EQDr.Shot = LA_STRING(shot);
EQDr.Time = LA_STRING(time);
//if(vec.size()<22) {if(mpi_rank < 1) cout << "Fail to read all parameters, File incomplete" << endl; EXIT;}

// private Variables
filename = name;

// Map Parameters
itt = int(vec[1]);
phistart = vec[7];
MapDirection = int(vec[8]);

// t grid (footprints only)
Nt = int(vec[6]);
tmin = vec[2];
tmax = vec[3];

// phi grid (footprints only)
Nphi = int(vec[0]);
phimin = vec[4];
phimax = vec[5];

// R grid (laminar only)
NR = int(vec[6]);
Rmin = vec[2];
Rmax = vec[3];

// Z grid (laminar only)
NZ = int(vec[0]);
Zmin = vec[4];
Zmax = vec[5];

// r grid
N = int(vec[6]);
Nr = int(sqrt(N));
rmin = vec[2];
rmax = vec[3];

// theta grid
Nth = int(sqrt(N));
thmin = vec[4];
thmax = vec[5];

verschieb = vec[0];

#if defined(ITER)
	// Particle Parameters
	Ekin = vec[16];
	lambda = vec[17];

	// Set switches
	which_target_plate = int(vec[9]);
	create_flag = int(vec[10]);
	useFcoil = int(vec[11]);
	useCcoil = int(vec[11]);
	useIcoil = int(vec[11]);
	useFilament = int(vec[12]);
	useTprofile = int(vec[13]);
	sigma = int(vec[14]);
	Zq = int(vec[15]);

	// External fields parameter
	response = int(vec[18]);
	response_field = int(vec[19]);
	useTprofile = 0;
#elif defined(NSTX)
	// Particle Parameters
	Ekin = vec[16];
	lambda = vec[17];

	// Set switches
	which_target_plate = int(vec[9]);
	create_flag = int(vec[10]);
	useFcoil = int(vec[11]);
	useCcoil = int(vec[11]);
	useIcoil = int(vec[11]);
	useFilament = int(vec[12]);
	useTprofile = int(vec[13]);
	sigma = int(vec[14]);
	Zq = int(vec[15]);

	// External fields parameter
	response = int(vec[18]);
	response_field = int(vec[19]);
	useTprofile = 0;
#elif defined(MAST)
	// Particle Parameters
	Ekin = vec[18];
	lambda = vec[19];

	// Set switches
	which_target_plate = int(vec[11]);
	create_flag = int(vec[12]);
	useCcoil = int(vec[13]);
	useIcoil = int(vec[14]);
	sigma = int(vec[16]);
	Zq = int(vec[17]);
	useFilament = int(vec[15]);

	// External fields parameter
	response = int(vec[9]);
	response_field = int(vec[10]);

	// Set unused Parameters to defaults
	useFcoil = 0;
	useTprofile = 0;
#elif defined(CMOD)
	// Particle Parameters
	Ekin = vec[18];
	lambda = vec[19];

	// Set switches
	which_target_plate = int(vec[11]);
	create_flag = int(vec[12]);
	useFcoil = 0;
	useCcoil = 0;
	useIcoil = 0;
	sigma = int(vec[16]);
	Zq = int(vec[17]);
	useFilament = int(vec[20]);

	// External fields parameter
	response = int(vec[9]);
	response_field = int(vec[10]);

	useTprofile = 0;
#else
	// Particle Parameters
	Ekin = vec[18];
	lambda = vec[19];

	// Set switches
	which_target_plate = int(vec[11]);
	create_flag = int(vec[12]);
	useFcoil = int(vec[13]);
	useCcoil = int(vec[14]);
	useIcoil = int(vec[15]);
	sigma = int(vec[16]);
	Zq = int(vec[17]);
	useFilament = int(vec[20]);

	// External fields parameter
	response = int(vec[9]);
	response_field = int(vec[10]);

	if(vec[21] > 1) useTprofile = 0;
	else useTprofile = int(vec[21]);
	if(vec.size() < 23) useBuswork = 0;
	else
	{
		if(vec[22] > 1) useBuswork = 0;
		else useBuswork = int(vec[22]);
	}
#endif
}


// ------------------- writeiodata ----------------------------------------------------------------------------------------
void IO::writeiodata(ofstream& out, double bndy[], vector<LA_STRING>& var)
{
int i;
out << "# " << program_name << endl;
out << "#-------------------------------------------------" << endl;
out << "### Parameterfile: " << filename << endl;
out << "# Shot: " << EQDr.Shot << endl;
out << "# Time: " << EQDr.Time << endl;
#ifdef m3dc1
	out << "#-------------------------------------------------" << endl;
	out << "### M3D-C1:" << endl;
	out << "# Plasma response (0=no, >1=yes): " << response << endl;
	out << "# Field (-1=M3D-C1 off, 0=Eq, 1=I-coil, 2=both): " << response_field << endl;
#endif
out << "#-------------------------------------------------" << endl;
out << "### Switches:" << endl;
#if defined(ITER)
	out << "# I-coil active (0=no, 1=yes): " << useIcoil << endl;
#elif defined(NSTX)
	out << "# EC-coil active (0=no, 1=yes): " << useIcoil << endl;
#elif defined(MAST)
	out << "# ECC-coil active (0=no, 1=yes): " << useCcoil << endl;
	out << "# I-coil active (0=no, 1=yes): " << useIcoil << endl;
#elif defined(CMOD)
#else
	out << "# F-coil active (0=no, 1=yes): " << useFcoil << endl;
	out << "# C-coil active (0=no, 1=yes): " << useCcoil << endl;
	out << "# I-coil active (0=no, 1=yes): " << useIcoil << endl;
	out << "# Bus work error field active (0=no, 1=yes): " << useBuswork << endl;
#endif
out << "# No. of current filaments (0=none): " << useFilament << endl;
out << "# Use Temperature Profile (0=off, 1=on): " << useTprofile << endl;
out << "# Target (0=cp, 1=inner, 2=outer, 3=shelf): " << which_target_plate << endl;
out << "# Create Points (0=r-grid, 1=r-random, 2=target, 3=psi-grid, 4=psi-random, 5=RZ-grid): " << create_flag << endl;
out << "# Direction of particles (1=co-pass, -1=count-pass, 0=field lines): " << sigma << endl;
out << "# Charge number of particles (=-1:electrons, >=1:ions): " << Zq << endl;
out << "# Boundary (0=Wall, 1=Box): " << simpleBndy << endl;
out << "#-------------------------------------------------" << endl;
out << "### Global Parameters:" << endl;
out << "# Steps till Output (ilt): " << output_step_size << endl;
out << "# Step size (dpinit): " << dpinit << endl;
out << "# Boundary Rmin: " << bndy[0] << endl;
out << "# Boundary Rmax: " << bndy[1] << endl;
out << "# Boundary Zmin: " << bndy[2] << endl;
out << "# Boundary Zmax: " << bndy[3] << endl;
out << "# Magnetic Axis: R0: " << EQDr.RmAxis << endl;
out << "# Magnetic Axis: Z0: " << EQDr.ZmAxis << endl;
out << "#-------------------------------------------------" << endl;
out << "### additional Parameters:" << endl;
for(i=0;i<psize;++i)
{
	out << "# " << pv[i].name << ": " << pv[i].wert << endl;
}
out << "#-------------------------------------------------" << endl;
out << "### Data:" << endl;
out << "# ";
for(i=0;i<int(var.size());i++) out << var[i] << "     ";
out << endl;
out << "#" << endl;
}

// ------------------- set_pv ---------------------------------------------------------------------------------------------
void IO::set_pv(int size)
{
psize = size;
pv = new parstruct[psize];
}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ------------- readparfile ---------------------------------------------------------------------------------------------
void readparfile_new(char* name, vector<double>& vec, int& shot, int& time, LA_STRING& path)
{
// Variables
int i;
int shotfound,timefound;
size_t found;
double val;
string parname;
string line,word;
vector<string> words;

// Read file
ifstream file;
file.open(name);
if(file.fail()==1) {cout << "Unable to open file " << name << endl; exit(0);}

// Read data
while(getline(file, line))
{
	if (line[0] == '#') // process and skip header lines
	{
		// split line into words
		found = 0;
		while ((found = line.find_first_of("\t ")) != std::string::npos)
		{
		    word = line.substr(0, found);
		    line.erase(0, found + 1);
		    if (word.length() == 0) continue;
		    words.push_back(word);
		}
		words.push_back(line); // last word in line

		// search for key words
		shotfound = -1;
		timefound = -1;
		for(i=0;i<words.size();i++)
		{
			if (int(words[i].find("Shot")) > -1) shotfound = i;
			if (int(words[i].find("Time")) > -1) timefound = i;
			if ((shotfound > -1) && (timefound > -1))
			{
				shot = stoi(words[shotfound + 1]);
				word = words[timefound + 1];
				if (int(word.find("ms")) > -1) word.erase(word.find("ms")); // erase all from beginning of 'ms' to end of string
				time = stoi(word);
				break;
			}
			if (int(words[i].find("Path")) > -1)
			{
				word = words[i+1];
				found = word.find_last_of("/");
				if (found < word.length() -1) word += "/";
				path = LA_STRING(word.c_str());
				break;
			}
		}

		words.clear();
		continue;
	}

    if (int(line.find_first_of('#')) > -1)	// does the line include a comment after the data
    {
    	found = line.find_first_of('#');
    	line = line.substr(0,found);
    }

    found = line.find_last_of("=");
    parname = line.substr(0,found);
    val = stod(line.substr(found+1));
    vec.push_back(val);
}
file.close();
}

#endif //  IO_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

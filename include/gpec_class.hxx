// Class loads and sets GPEC input
// A.Wingen						12.11.24


// Define
//--------
#ifndef GPEC_CLASS_INCLUDED
#define GPEC_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <andi.hxx>
#include <io_class.hxx>
#include <splines.hxx>

using namespace blitz;

// Prototypes
//-----------

// Golbal Parameters
//------------------
extern ofstream ofs2;

//--------- Begin Class GPEC ---------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Open and load GPEC data files, return B field
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class GPEC
{
private:
	// Parameter
	int N;		// number of toroidal modes

	// Member Variables
	vector<LA_STRING> filenames;
	vector<double> scale;
	vector<double> phase;

public:
	// Member Variables
	vector<GPEC_mode> Modes;

	// Constructors
	GPEC();								// Default Constructor

	// Member-Functions
	int getB(double R, double phi, double Z, double& Br, double& Bp, double& Bz);
	int read_gpecsup(LA_STRING supPath="./");
	void load(IO& PAR, int mpi_rank, LA_STRING supPath="./");

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
GPEC::GPEC()
{
N = 0;
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- getB -------------------------------------------------------------------------------------------------
// returns the magnetic field at location (R,phi,Z)
int GPEC::getB(double R, double phi, double Z, double& Br, double& Bp, double& Bz)
{
int n,chk;
double Bnr, Bnp, Bnz;

Br = 0; Bp = 0; Bz = 0;
for(n=0;n<N;n++)
{
	 chk = Modes[n].getBn(R, phi + phase[n], Z, Bnr, Bnp, Bnz);
	 if(chk < 0) {return -1;}
	 Br += scale[n] * Bnr;
	 Bp += scale[n] * Bnp;
	 Bz += scale[n] * Bnz;
}

return 0;
}

// ----------------- read_gpecsup ----------------------------------------------------------------------------------------
// Prepare loading M3D-C1
int GPEC::read_gpecsup(LA_STRING supPath)
{
string name, line;
double sc,ph;

ifstream file;
file.open(supPath + "gpecsup.in");
if(file.fail()==1) {cout << "Unable to open gpecsup.in file" << endl; exit(0);}

// m3dc1 control file found
N = 0;
while(getline(file, line))
{
	if(line.length() < 1) continue; 	// blank lines anywhere don't matter
	stringstream ss(line);
	if( ss >> name >> sc >> ph)	// assigns any one of file, scale and ph if possible
	{
		if(fabs(ph) > pi) ph /= rTOd;	// convert phase to radiants
		phase.push_back(ph);
	}
	else	// ph assignment was not possible, others are set though
	{
		phase.push_back(0);
	}
	filenames.push_back(name.c_str());
	scale.push_back(sc);
	N += 1;
}

file.close();
return 0;
}


// ----------------- load -------------------------------------------------------------------------------------------------
void GPEC::load(IO& PAR, int mpi_rank, LA_STRING supPath)
{
int j,chk;

if(PAR.response_field == 5)		// only if GPEC is used
{
	if(mpi_rank < 1) cout << "Loading GPEC output" << endl;
	ofs2 << "Loading GPEC output" << endl;

	read_gpecsup(supPath);

	//if(mpi_rank < 1) cout << "Plasma response (0 = off, 1 = on): " << PAR.response << "\t" << "GPEC perturbation (0 = off, 1 = on): " << PAR.response_field << endl;
	//ofs2 << "Plasma response (0 = off, 1 = on): " << PAR.response << "\t" << "GPEC perturbation (0 = off, 1 = on): " << PAR.response_field << endl;





	for(j=0;j<N;j++)
	{
		if(mpi_rank < 1) cout << "GPEC file: " << filenames[j] <<  " -> perturbation scaling factor: " << scale[j] << "  and phase: " << phase[j]*rTOd << endl;
		ofs2 << "GPEC file: " << filenames[j] <<  " -> perturbation scaling factor: " << scale[j] << "  and phase: " << phase[j]*rTOd << endl;
		Modes[j] = ??? // Problem of copy and destruction of GPEC_mode class object
	}

	if(mpi_rank < 1) cout << "Coils turned off: perturbation (only) already included in GPEC output" << endl;
	ofs2 << "Coils turned off: perturbation (only) already included in GPEC output" << endl;

	PAR.useFcoil = 0;
	PAR.useCcoil = 0;
	PAR.useIcoil = 0;
	PAR.useBuswork = 0;
	PAR.useBcoil = 0;
}
}





//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//------------------------ End of Class GPEC ------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Begin Class GPEC_mode -----------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Loads a single n-mode and returns Bn(R,Z)
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class GPEC_mode
{
private:
	// Parameter
	int NR;		// Number of points in R
	int NZ;		// Number of points in Z

public:
	// Member Variables
	int n;		// toroidal mode number

	Array<double,1> R;		// R grid
	Array<double,1> Z;		// Z grid

	GPEC_var L;		// ???;  i=1:Nr rows, j=1:Nz columns
	GPEC_var ReBr;	// Real(B_r) (R,Z);  i=1:Nr rows, j=1:Nz columns
	GPEC_var ImBr;	// Imag(B_r) (R,Z);  i=1:Nr rows, j=1:Nz columns
	GPEC_var ReBp;	// Real(B_phi) (R,Z);  i=1:Nr rows, j=1:Nz columns
	GPEC_var ImBp;	// Imag(B_phi) (R,Z);  i=1:Nr rows, j=1:Nz columns
	GPEC_var ReBz;	// Real(B_z) (R,Z);  i=1:Nr rows, j=1:Nz columns
	GPEC_var ImBz;	// Imag(B_z) (R,Z);  i=1:Nr rows, j=1:Nz columns

	// Constructors
	GPEC_mode();								// Default Constructor

	// Member-Functions
	int getBn(double r, double phi, double z, double& Bnr, double& Bnp, double& Bnz);
	void ReadData(LA_STRING filename);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
GPEC_Nmode::GPEC_mode()
{
TinyVector <int,1> index(1);

NR = 129;
NZ = 129;

R.resize(NR);	R.reindexSelf(index);
Z.resize(NZ);	Z.reindexSelf(index);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- getB -------------------------------------------------------------------------------------------------
// returns the magnetic field at location (R,phi,Z) for toroidal mode n
int GPEC_mode::getBn(double r, double phi, double z, double& Bnr, double& Bnp, double& Bnz)
{
int i;
int chk = 0;
double rebr,imbr,rebp,imbp,rebz,imbz,dummy;

chk += ReBr.ev(r,z,rebr,dummy,dummy);
chk += ImBr.ev(r,z,imbr,dummy,dummy);
chk += ReBp.ev(r,z,rebp,dummy,dummy);
chk += ImBp.ev(r,z,imbp,dummy,dummy);
chk += ReBz.ev(r,z,rebz,dummy,dummy);
chk += ImBz.ev(r,z,imbz,dummy,dummy);
if(chk < 0) {return -1;}

cosnp = cos(n*phi);
sinnp = sin(n*phi);

Bnr = 2*(rebr * cosnp - imbr * sinnp);
Bnp = 2*(rebp * cosnp - imbp * sinnp);
Bnz = 2*(rebz * cosnp - imbz * sinnp);

return chk;
}


//---------------------- ReadData -------------------------------------------------------------------------------------
void GPEC_mode::ReadData(LA_STRING filename)
{
int i;
int j = 0;
int HeaderDone = 0;
string line,word;
vector<string> words;
Array<int,1> nums(3);

// Open File
ifstream file;
file.open(filename);
if(file.fail()==1) {cout << "Unable to open file " << filename << endl; exit(0);}

// Read dimensions, last two enties in first line
while(HeaderDone == 0)
{
	getline(file, line);
	if(line.length() < 1) continue; 	// blank lines anywhere don't matter
	words = split(line);
	if (words.size() == 0) continue;	// line only has blanks and no words
	for(i=0;i<words.size();i++)
	{
		if (int(words[i].find("real(b_r)")) > -1) HeaderDone = 1;
		if (int(words[i].find("=")) > -1) {nums(j) = atoi(words[i + 1].c_str()); j++;}
	}
}
n = nums(0);	// toroidal mode number
NR = nums(1);	// Number of Points in R-direction
NZ = nums(2);	// Number of Points in Z-direction

// Rezize Arrays; All Arrays start with index 1
R.resize(NR);	//  if size = 129 nothing is done!
Z.resize(NZ);	//  if size = 129 nothing is done!
L.resize(NR,NZ);	// This calls resize in GPEC_var
ReBr.resize(NR,NZ);
ImBr.resize(NR,NZ);
ReBp.resize(NR,NZ);
ImBp.resize(NR,NZ);
ReBz.resize(NR,NZ);
ImBz.resize(NR,NZ);

// Read Data
int rows = NR*NZ;
int row;	// row = (i-1)*NZ + j
for(row=1;row<=rows;row++)
{
	i = row/NZ + 1;		// uses integer division
	j = row - (i-1)*NZ;

	file >> L[i,j];
	file >> R(i);
	file >> Z(j);
	file >> ReBr[i,j];
	file >> ImBr[i,j];
	file >> ReBz[i,j];
	file >> ImBz[i,j];
	file >> ReBp[i,j];
	file >> ImBp[i,j];
}
file.close();

// Set up interpolation
ReBr.set(R,Z);
ImBr.set(R,Z);
ReBp.set(R,Z);
ImBp.set(R,Z);
ReBz.set(R,Z);
ImBz.set(R,Z);
}
//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//--------- End Class GPEC_mode -------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//--------- Begin Class GPEC_var ------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// A single field component and its interpolation
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class GPEC_var
{
private:
	// Parameter
	int NR;		// Number of points in R
	int NZ;		// Number of points in Z

	// Member Variables
	double dR;		// grid distance in R direction
	double dZ;		// grid distance in Z direction

	Array<double,1> R;		// R grid
	Array<double,1> Z;		// Z grid

	Array<double,2> dBdR;		// Derivative in R-direction
	Array<double,2> dBdZ;		// Derivative in Z-direction
	Array<double,2> d2BdRdZ;	// Cross-derivative

	Array<double,4> Ca;			// Coefficients for bcuint

public:
	// Member Variables
	Array<double,2> B;		// i=1:Nr rows, j=1:Nz columns

	// Constructors & Operator
	GPEC_var();										// Default Constructor
	const double &operator[] (int i, int j) const;	// Get element
	double &operator[](int i, int j);				// Assign element
	GPEC_var& operator =(const GPEC_var& var);		// Operator =

	// Member-Functions
	int ev(const double x1, const double x2, double& y, double& dy1, double& dy2);
	void set(Array<double,1>& Rin, Array<double,1>& Zin);


}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
GPEC_var::GPEC_var()
{
NR = 0;
NZ = 0;
dR = 0;
dZ = 0;
}

//--------- [] Operator --------------------------------------------------------------------------------------------------
const double &GPEC_var::operator[](int i, int j) const		//  Get element
{
	return B(i,j);
}

double &GPEC_var::operator[](int i, int j)		//  Assign element
{
	return B(i,j);
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
GPEC_var& GPEC_var::operator =(const GPEC_var& var)
{
if (this == &var) return(*this);	    // if: x=x
NR = var.NR;
NZ = var.NZ;
dR = var.dR;
dZ = var.dZ;
B.reference(var.B.copy());
R.reference(var.R.copy());
Z.reference(var.Z.copy());
dBdR.reference(var.dBdR.copy());
dBdZ.reference(var.dBdZ.copy());
d2BdRdZ.reference(var.d2BdRdZ.copy());
Ca.reference(var.Ca.copy());
return(*this);
}


//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- resize --------------------------------------------------------------------------------------------------
void GPEC_var::resize(int nr, int nz)
{
TinyVector <int,1> index(1);
TinyVector <int,2> index2(1,1);
TinyVector <int,4> index4(1,1,1,1);

NR = nr;
NZ = nz;

B.resize(NR,NZ);	B.reindexSelf(index2);
dBdR.resize(NR,NZ);	dBdR.reindexSelf(index2);
dBdZ.resize(NR,NZ);	dBdZ.reindexSelf(index2);
d2BdRdZ.resize(NR,NZ);	d2BdRdZ.reindexSelf(index2);
Ca.resize(NR-1,NZ-1,4,4);	Ca.reindexSelf(index4);
}

// ----------------- ev -------------------------------------------------------------------------------------------------
// returns the magnetic field at location (R,phi,Z)
int GPEC_var::ev(const double x1, const double x2, double& y, double& dy1, double& dy2)
{
int chk = 0;
if(x1>R(NR) || x1<R(1) || x2>Z(NZ) || x2<Z(1))	{return -1;}
chk = bcuint(R,Z,Ca,dR,dZ,x1,x2,y,dy1,dy2);
if(chk == -1) {return -1;}
return 0;
}

//---------------------- set -------------------------------------------------------------------------------------
void GPEC_var::set(Array<double,1>& Rin, Array<double,1>& Zin)
{
int i,j;
Range all = Range::all();

R = Rin.reference();
Z = Zin.reference();
dR = (R(NR) - R(1))/double(NR-1);
dZ = (Z(NZ) - Z(1))/double(NR-1);

// Prepare Bicubic Interpolation of B -> get gradients and cross-derivative
dBdR.resize(NR,NZ);
dBdZ.resize(NR,NZ);
d2BdRdZ.resize(NR,NZ);

bcuderiv(B,dR,dZ,dBdR,dBdZ,d2BdRdZ);

// Get the c's for bcuint, as done by bcucof
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
Array<double,2> slice;
Ca.resize(NR-1,NZ-1,4,4);
for(i=1;i<NR;i++)
{
	for(j=1;j<NZ;j++)
	{
		y_sq(1) = B(i,j); y_sq(2) = B(i+1,j); y_sq(3) = B(i+1,j+1); y_sq(4) = B(i,j+1);
		y1_sq(1) = dBdR(i,j); y1_sq(2) = dBdR(i+1,j); y1_sq(3) = dBdR(i+1,j+1); y1_sq(4) = dBdR(i,j+1);
		y2_sq(1) = dBdZ(i,j); y2_sq(2) = dBdZ(i+1,j); y2_sq(3) = dBdZ(i+1,j+1); y2_sq(4) = dBdZ(i,j+1);
		y12_sq(1) = d2BdRdZ(i,j); y12_sq(2) = d2BdRdZ(i+1,j); y12_sq(3) = d2BdRdZ(i+1,j+1); y12_sq(4) = d2BdRdZ(i,j+1);

		slice.reference(Ca(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,slice);
	}
}
}

//----------------------- End of Member Functions -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
//--------- End Class GPEC_var --------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  GPEC_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

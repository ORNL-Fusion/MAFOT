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

	// Member Variables

	// Member-Functions

public:
	// Member Variables

	// Constructors
	GPEC();								// Default Constructor

	// Member-Functions
	int getB(double R, double phi, double Z, double& Br, double& Bp, double& Bz)
	int read_gpecsup(LA_STRING supPath="./");
	void load(IO& PAR, int mpi_rank);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
GPEC::GPEC()
{
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- getB -------------------------------------------------------------------------------------------------
// returns the magnetic field at location (R,phi,Z)
int GPEC::getB(double R, double phi, double Z, double& Br, double& Bp, double& Bz)
{
int i,chk;


return chk;
}

// ----------------- read_gpecsup ----------------------------------------------------------------------------------------
// Prepare loading M3D-C1
int GPEC::read_gpecsup(LA_STRING supPath)
{
//freopen("/dev/null", "w", stderr);	// suppress stderr output
int i;
string file, line;
ifstream in;
double ph;

in.open(supPath + "m3dc1sup.in");
if(in.fail()==1) // no m3dc1 control file found -> use default: scale by coils and filename = "C1.h5"
{
	nfiles = 1;
	filenames[0] = "C1.h5";
	scale[0] = 1;
	phase[0] = 0;
	in.close();
	return -1;
}
else	// m3dc1 control file found
{
	i = 0;
	while(getline(in, line))
	{
		if(line.length() < 1) continue; 	// blank lines anywhere don't matter
		stringstream ss(line);
		if( ss >> file >> scale[i] >> ph)	// assigns any one of file, scale and ph if possible
		{
			if(fabs(ph) > pi) ph /= rTOd;	// convert phase to radiants
			phase[i] = ph;
		}
		else	// ph assignment was not possible, others are set though
		{
			phase[i] = 0;
		}
		filenames[i] = file.c_str();
		i += 1;
	}
	nfiles = i;
}
in.close();
return 0;
}


// ----------------- load -------------------------------------------------------------------------------------------------
void GPEC::load(IO& PAR, int mpi_rank)
{
int j,chk;

if(mpi_rank < 1) cout << "Loading M3D-C1 output file C1.h5" << endl;
ofs2 << "Loading M3D-C1 output file C1.h5" << endl;

chk = open_source(PAR.response, PAR.response_field);	// load C1.h5
if(chk != 0) {if(mpi_rank < 1) cout << "Error loding C1.h5 file" << endl; EXIT;}
if(nonlinear) PAR.response_field = 2;	// nonlinear runs should only use the full field

if(mpi_rank < 1) cout << "Plasma response (0 = off, 1 = on): " << PAR.response << "\t" << "Field (0 = Eq, 1 = I-coil, 2 = total): " << PAR.response_field << endl;
ofs2 << "Plasma response (0 = off, 1 = on): " << PAR.response << "\t" << "Field (0 = Eq, 1 = I-coil, 2 = total): " << PAR.response_field << endl;

if(PAR.response_field > 0)		// Perturbation already included in M3D-C1 output
{
	for(j=0;j<nfiles;j++)
	{
		if(mpi_rank < 1) cout << "M3D-C1 file: " << filenames[j] <<  " -> perturbation scaling factor: " << scale[j] << "  and phase: " << phase[j]*rTOd << endl;
		ofs2 << "M3D-C1 file: " << filenames[j] <<  " -> perturbation scaling factor: " << scale[j] << "  and phase: " << phase[j]*rTOd << endl;
	}

	if(mpi_rank < 1) cout << "Coils turned off: perturbation (only) already included in M3D-C1 output" << endl;
	ofs2 << "Coils turned off: perturbation (only) already included in M3D-C1 output" << endl;

	PAR.useFcoil = 0;
	PAR.useCcoil = 0;
	PAR.useIcoil = 0;
	PAR.useBuswork = 0;
	PAR.useBcoil = 0;
}
}



//--------------------- Private Member Functions --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



// ----------------- make_psi ---------------------------------------------------------------------------------------------
int GPEC::make_psi(void)
{
int ierr,i;


return ierr;
}

// ----------------- axis0 ------------------------------------------------------------------------------------------------
int GPEC::axis0(void)
{
int ir0, iz0, ierr;


return ierr;
}

// ----------------- axis -------------------------------------------------------------------------------------------------
int GPEC::axis(void)
{
int chk,i;

return chk;
}



//------------------------ End of Class GPEC -----------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Begin Class GPEC_Nmode ---------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Loads a single n-mode and returns Bn(R,Z)
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class GPEC_Nmode
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
	Array<double,2> ReBr;	// Real(B_r) (R,Z);  i=1:Nr rows, j=1:Nz columns
	Array<double,2> ImBr;	// Imag(B_r) (R,Z);  i=1:Nr rows, j=1:Nz columns
	Array<double,2> ReBp;	// Real(B_phi) (R,Z);  i=1:Nr rows, j=1:Nz columns
	Array<double,2> ImBp;	// Imag(B_phi) (R,Z);  i=1:Nr rows, j=1:Nz columns
	Array<double,2> ReBz;	// Real(B_z) (R,Z);  i=1:Nr rows, j=1:Nz columns
	Array<double,2> ImBz;	// Imag(B_z) (R,Z);  i=1:Nr rows, j=1:Nz columns

	Array<double,2> dReBrdR;		// Derivative in R-direction
	Array<double,2> dReBrdZ;		// Derivative in Z-direction
	Array<double,2> d2ReBrdRdZ;		// Cross-derivative
	Array<double,4> Ca_ReBr;		// Coefficients for bcuint

	Array<double,4> Ca_ImBr;
	Array<double,4> Ca_ReBp;
	Array<double,4> Ca_ImBp;
	Array<double,4> Ca_ReBz;
	Array<double,4> Ca_ImBz;

	// Member-Functions

public:
	// Member Variables
	int n;		// toroidal mode number

	// Constructors
	GPEC_Nmode();								// Default Constructor

	// Member-Functions
	int getBn(double R, double Z, double& Bnr, double& Bnp, double& Bnz)
	void ReadData(LA_STRING file)

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
GPEC_Nmode::GPEC_Nmode()
{
TinyVector <int,1> index(1);
TinyVector <int,2> index2(1,1);
TinyVector <int,2> index2_10(1,0);
TinyVector <int,4> index4(1,1,1,1);

NR = 129;
NZ = 129;

ReBr.resize(NR,NZ);	ReBr.reindexSelf(index2);
ImBr.resize(NR,NZ);	ImBr.reindexSelf(index2);
ReBp.resize(NR,NZ);	ReBp.reindexSelf(index2);
ImBp.resize(NR,NZ);	ImBp.reindexSelf(index2);
ReBz.resize(NR,NZ);	ReBz.reindexSelf(index2);
ImBz.resize(NR,NZ);	ImBz.reindexSelf(index2);

Ca_ReBr.resize(NR-1,NZ-1,4,4);	Ca_ReBr.reindexSelf(index4);
Ca_ImBr.resize(NR-1,NZ-1,4,4);	Ca_ImBr.reindexSelf(index4);
Ca_ReBp.resize(NR-1,NZ-1,4,4);	Ca_ReBp.reindexSelf(index4);
Ca_ImBp.resize(NR-1,NZ-1,4,4);	Ca_ImBp.reindexSelf(index4);
Ca_ReBz.resize(NR-1,NZ-1,4,4);	Ca_ReBz.reindexSelf(index4);
Ca_ImBz.resize(NR-1,NZ-1,4,4);	Ca_ImBz.reindexSelf(index4);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- getB -------------------------------------------------------------------------------------------------
// returns the magnetic field at location (R,phi,Z)
int GPEC_Nmode::getBn(double R, double Z, double& Bnr, double& Bnp, double& Bnz)
{
int i,chk;


return chk;
}


//---------------------- ReadData -------------------------------------------------------------------------------------
void GPEC_Nmode::ReadData(LA_STRING file)
{
int i,j;

// Set R, Z and psi Arrays
R.resize(NR);
Z.resize(NZ);




// Prepare Bicubic Interpolation of psiRZ -> get gradients and cross-derivative
dpsidR.resize(NR,NZ);
dpsidZ.resize(NR,NZ);
d2psidRdZ.resize(NR,NZ);

bcuderiv(psiRZ,dR,dZ,dpsidR,dpsidZ,d2psidRdZ);

// Get the c's for bcuint, as done by bcucof
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
Array<double,2> slice;
Ca.resize(NR-1,NZ-1,4,4);
for(i=1;i<NR;i++)
{
	for(j=1;j<NZ;j++)
	{
		y_sq(1) = psiRZ(i,j); y_sq(2) = psiRZ(i+1,j); y_sq(3) = psiRZ(i+1,j+1); y_sq(4) = psiRZ(i,j+1);
		y1_sq(1) = dpsidR(i,j); y1_sq(2) = dpsidR(i+1,j); y1_sq(3) = dpsidR(i+1,j+1); y1_sq(4) = dpsidR(i,j+1);
		y2_sq(1) = dpsidZ(i,j); y2_sq(2) = dpsidZ(i+1,j); y2_sq(3) = dpsidZ(i+1,j+1); y2_sq(4) = dpsidZ(i,j+1);
		y12_sq(1) = d2psidRdZ(i,j); y12_sq(2) = d2psidRdZ(i+1,j); y12_sq(3) = d2psidRdZ(i+1,j+1); y12_sq(4) = d2psidRdZ(i,j+1);

		slice.reference(Ca(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,slice);
	}
}


}


//--------- Begin Class GPEC_var ---------------------------------------------------------------------------------------------
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

	Array<double,2> B;		// i=1:Nr rows, j=1:Nz columns

	Array<double,2> dBdR;		// Derivative in R-direction
	Array<double,2> dBdZ;		// Derivative in Z-direction
	Array<double,2> d2BdRdZ;	// Cross-derivative
	Array<double,4> Ca;			// Coefficients for bcuint


	// Member-Functions

public:
	// Member Variables
	int n;		// toroidal mode number

	// Constructors
	GPEC_var();								// Default Constructor

	// Member-Functions


}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
GPEC_var::GPEC_var()
{
TinyVector <int,1> index(1);
TinyVector <int,2> index2(1,1);
TinyVector <int,2> index2_10(1,0);
TinyVector <int,4> index4(1,1,1,1);

NR = 129;
NZ = 129;

ReBr.resize(NR,NZ);	ReBr.reindexSelf(index2);
ImBr.resize(NR,NZ);	ImBr.reindexSelf(index2);
ReBp.resize(NR,NZ);	ReBp.reindexSelf(index2);
ImBp.resize(NR,NZ);	ImBp.reindexSelf(index2);
ReBz.resize(NR,NZ);	ReBz.reindexSelf(index2);
ImBz.resize(NR,NZ);	ImBz.reindexSelf(index2);

Ca_ReBr.resize(NR-1,NZ-1,4,4);	Ca_ReBr.reindexSelf(index4);
Ca_ImBr.resize(NR-1,NZ-1,4,4);	Ca_ImBr.reindexSelf(index4);
Ca_ReBp.resize(NR-1,NZ-1,4,4);	Ca_ReBp.reindexSelf(index4);
Ca_ImBp.resize(NR-1,NZ-1,4,4);	Ca_ImBp.reindexSelf(index4);
Ca_ReBz.resize(NR-1,NZ-1,4,4);	Ca_ReBz.reindexSelf(index4);
Ca_ImBz.resize(NR-1,NZ-1,4,4);	Ca_ImBz.reindexSelf(index4);
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- getB -------------------------------------------------------------------------------------------------
// returns the magnetic field at location (R,phi,Z)
int GPEC_var::eval(double R, double Z, double& val)
{
int i,chk;


return chk;
}


//---------------------- ReadData -------------------------------------------------------------------------------------
void GPEC_var::set(double& R, double& Z, double& Val)
{
int i,j;

NR = R.size();
NZ = Z.size();
B = Val;



// Prepare Bicubic Interpolation of psiRZ -> get gradients and cross-derivative
dpsidR.resize(NR,NZ);
dpsidZ.resize(NR,NZ);
d2psidRdZ.resize(NR,NZ);

bcuderiv(psiRZ,dR,dZ,dpsidR,dpsidZ,d2psidRdZ);

// Get the c's for bcuint, as done by bcucof
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
Array<double,2> slice;
Ca.resize(NR-1,NZ-1,4,4);
for(i=1;i<NR;i++)
{
	for(j=1;j<NZ;j++)
	{
		y_sq(1) = psiRZ(i,j); y_sq(2) = psiRZ(i+1,j); y_sq(3) = psiRZ(i+1,j+1); y_sq(4) = psiRZ(i,j+1);
		y1_sq(1) = dpsidR(i,j); y1_sq(2) = dpsidR(i+1,j); y1_sq(3) = dpsidR(i+1,j+1); y1_sq(4) = dpsidR(i,j+1);
		y2_sq(1) = dpsidZ(i,j); y2_sq(2) = dpsidZ(i+1,j); y2_sq(3) = dpsidZ(i+1,j+1); y2_sq(4) = dpsidZ(i,j+1);
		y12_sq(1) = d2psidRdZ(i,j); y12_sq(2) = d2psidRdZ(i+1,j); y12_sq(3) = d2psidRdZ(i+1,j+1); y12_sq(4) = d2psidRdZ(i,j+1);

		slice.reference(Ca(i,j,all,all));
		bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,slice);
	}
}


}

#endif //  GPEC_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

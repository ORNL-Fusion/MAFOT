// Class loads and sets M3DC1 input
// A.Wingen						2.7.15


// Define
//--------
#ifndef M3DC1_CLASS_INCLUDED
#define M3DC1_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <andi.hxx>
#include <fusion_io_defs.h>
#include <fusion_io_c.h>

using namespace blitz;

// Prototypes
//-----------

// Golbal Parameters
//------------------
extern ofstream ofs2;

//--------- Begin Class M3DC1 ---------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Open and load C1.h5 files
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class M3DC1
{
private:
	// Parameter
	static const int nfiles_max = 10;	// max number of files = 10

	// Member Variables
	int isrc[nfiles_max];				// handle for sources
	LA_STRING filenames[nfiles_max];	// M3DC1 files
	double scale[nfiles_max];			// scaling factor for linear runs

	// Member-Functions
	int open_source(int response, int response_field);
	int make_psi(void);

public:
	// Member Variables
	int nfiles;				// actual number of files
	bool nonlinear;			// flag if run is linear or not
	int imag[nfiles_max];	// handle for magnetic fields
	int ia;					// handle for the vector potential
	double psi_axis;		// poloidal flux on axis
	double psi_lcfs;		// poloidal flux at lscf

	// Constructors
	M3DC1();								// Default Constructor

	// Member-Functions
	int read_m3dc1sup(LA_STRING supPath);
	void scale_from_coils(double currents[], int loops, int turns, double adj = 1);
	void load(IO& PAR, int mpi_rank);
	void unload(void);


}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
M3DC1::M3DC1()
{
nfiles = 1;				// default: at least one file
nonlinear = true;		// to be on the safe side
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- read_m3dc1sup ----------------------------------------------------------------------------------------
// Prepare loading M3D-C1
int M3DC1::read_m3dc1sup(LA_STRING supPath)
{
int i;
string file;
ifstream in;

in.open(supPath + "m3dc1sup.in");
if(in.fail()==1) // no m3dc1 control file found -> use default: scale by coils and filename = "C1.h5"
{
	nfiles = 1;
	filenames[0] = "C1.h5";
	in.close();
	return -1;
}
else	// m3dc1 control file found
{
	i = 0;
	while(in.eof()==0) // Last row is read twice --- can't be changed --- -> i-1 is actual number of rows in file
	{
		in >> file;
		in >> scale[i];
		filenames[i] = file.c_str();
		i += 1;
	}
	nfiles = i - 1;
}
in.close();
return 0;
}

// ----------------- scale_from_coils -------------------------------------------------------------------------------------
// no m3dc1sup.in file found -> scale from coil currents
void M3DC1::scale_from_coils(double currents[], int loops, int turns, double adj)
{
int i;
// Scale perturbation in M3D-C1 according to coil currents (current = scale * 1kA)
scale[0] = 0;
for(i=0;i<loops;i++) scale[0] += fabs(currents[i]*adj);
scale[0] /= 1000.0*turns;
scale[0] *= -sign(currents[0]);		// proper phasing: first current < 0 -> 60° phase, else 0° phase
}

// ----------------- load -------------------------------------------------------------------------------------------------
void M3DC1::load(IO& PAR, int mpi_rank)
{
int j,chk;

if(mpi_rank < 1) cout << "Loading M3D-C1 output file C1.h5" << endl;
ofs2 << "Loading M3D-C1 output file C1.h5" << endl;

chk = open_source(PAR.response, PAR.response_field);	// load C1.h5
if(chk != 0) {if(mpi_rank < 1) cout << "Error loding C1.h5 file" << endl; EXIT;}

if(mpi_rank < 1) cout << "Plasma response (0 = off, 1 = on): " << PAR.response << "\t" << "Field (0 = Eq, 1 = I-coil, 2 = total): " << PAR.response_field << endl;
ofs2 << "Plasma response (0 = off, 1 = on): " << PAR.response << "\t" << "Field (0 = Eq, 1 = I-coil, 2 = total): " << PAR.response_field << endl;

if(PAR.response_field > 0)		// Perturbation already included in M3D-C1 output
{
	for(j=0;j<nfiles;j++)
	{
		if(mpi_rank < 1) cout << "M3D-C1 file: " << filenames[j] <<  " -> perturbation scaling factor: " << scale[j] << endl;
		ofs2 << "M3D-C1 file: " << filenames[j] <<  " -> perturbation scaling factor: " << scale[j] << endl;
	}

	if(mpi_rank < 1) cout << "Coils turned off: perturbation (only) already included in M3D-C1 output" << endl;
	ofs2 << "Coils turned off: perturbation (only) already included in M3D-C1 output" << endl;

	PAR.useFcoil = 0;
	PAR.useCcoil = 0;
	PAR.useIcoil = 0;
}
}

// ----------------- unload -----------------------------------------------------------------------------------------------
void M3DC1::unload(void)
{
int i, ierr;
ierr = fio_close_field(ia);
for(i=0;i<nfiles;i++)
{
	ierr = fio_close_field(imag[i]);
	ierr = fio_close_source(isrc[i]);
}
}

//--------------------- Private Member Functions --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// ----------------- open_source ------------------------------------------------------------------------------------------
int M3DC1::open_source(int response, int response_field)
{
int i, ierr, ierr2, ierr3;

// open file(s)
ierr =  fio_open_source(FIO_M3DC1_SOURCE, filenames[0], &(isrc[0]));

// Set options appropriate to this source
ierr2 = fio_get_options(isrc[0]);
ierr2 += fio_set_int_option(FIO_TIMESLICE, response);	// response = 1: plasma response solution; response = 0: vacuum field;

switch(response_field)
{
case 0: 	// M3D-C1: equilibrium field only
	ierr2 += fio_set_int_option(FIO_PART, FIO_EQUILIBRIUM_ONLY);
	break;

case 1:		// M3D-C1: I-coil perturbation field only
    ierr2 += fio_set_real_option(FIO_LINEAR_SCALE, scale[0]);	// Scale I-coil perturbation in M3D-C1 according to diiidsup.in (current = scale * 1kA)
	ierr2 += fio_set_int_option(FIO_PART, FIO_PERTURBED_ONLY);
	break;

case 2:		// M3D-C1: total field
    ierr2 += fio_set_real_option(FIO_LINEAR_SCALE, scale[0]);	// Scale I-coil perturbation in M3D-C1 according to diiidsup.in (current = scale * 1kA)
	ierr2 += fio_set_int_option(FIO_PART, FIO_TOTAL);
	break;
}

// set field handle
ierr2 += fio_get_field(isrc[0], FIO_MAGNETIC_FIELD, &(imag[0]));

// further sources only contribute the perturbation, since the equilibrium is already in source 0
for(i=1;i<nfiles;i++)
{
	// open file(s)
	ierr +=  fio_open_source(FIO_M3DC1_SOURCE, filenames[i], &(isrc[i]));

    // Set options appropriate to this source
    ierr2 += fio_get_options(isrc[i]);
    ierr2 += fio_set_int_option(FIO_TIMESLICE, response);	// response = 1: plasma response solution; response = 0: vacuum field;

    // M3D-C1: I-coil perturbation field only
    ierr2 += fio_set_real_option(FIO_LINEAR_SCALE, scale[i]);	// Scale I-coil perturbation in M3D-C1 according to diiidsup.in (current = scale * 1kA)
    ierr2 += fio_set_int_option(FIO_PART, FIO_PERTURBED_ONLY);

    // set field handle
    ierr2 += fio_get_field(isrc[i], FIO_MAGNETIC_FIELD, &(imag[i]));
}

// prepare psi eval
ierr3 = make_psi();

return ierr+ierr2+ierr3;
}

// ----------------- make_psi ---------------------------------------------------------------------------------------------
int M3DC1::make_psi(void)
{
int ierr, ierr2;
int ipsi_axis, ipsi_lcfs;

// Set options appropriate to this source
ierr2 = fio_get_options(isrc[0]);
//ierr2 += fio_set_int_option(FIO_PART, FIO_EQUILIBRIUM_ONLY);	// works only in a linear run, where equilibrium and perturbation are separable; in a non-linear run one has to evaluate the vector-potential at several toroidal locations and average it to get an (almost) axisymmetric psi
nonlinear = true;

// set field handle
ierr2 += fio_get_field(isrc[0], FIO_VECTOR_POTENTIAL, &ia);

// get psi_axis and psi_sep
ierr = fio_get_series(isrc[0], FIO_MAGAXIS_PSI, &ipsi_axis);
ierr += fio_get_series(isrc[0], FIO_LCFS_PSI, &ipsi_lcfs);
ierr += fio_eval_series(ipsi_axis, 0., &psi_axis);
ierr += fio_eval_series(ipsi_lcfs, 0., &psi_lcfs);
ierr += fio_close_series(ipsi_axis);
ierr += fio_close_series(ipsi_lcfs);

return ierr+ierr2;
}

//------------------------ End of Class M3DC1 -----------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  M3DC1_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

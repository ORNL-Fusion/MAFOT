// Class loads and sets M3DC1 input
// A.Wingen						2.7.15

/* Details on how to scale the perturbation
The field at any point is the sum of the equilibrium part (from the “equilibrium” time slice) and the perturbed part (from the “time_###” time slice).  
In linear calculations the perturbed part is usually calculated assuming 1 kA in the I-coils, so that part needs to be scaled to the experimental 
value to compare with the experimental data.  You just need to multiply the perturbed part by the appropriate factor. 
The appropriate factor will just be the experimental amplitude of the current in the I-coils, times a geometric factor.  
The factor is different from 1 and different for each toroidal mode number, because of the toroidal discreteness of the coils.

n = 3: the geometric factor is 4/pi 

To match inner footprints in practice:
n = 3	scale factor: 1.116	phase: -39 deg
n = 1	scale factor: 4		phase: -72 deg	(C-coils + F-coil error)

Furthermore, if the M3D-C1 run was done using ProbeG data, then the perturbation needs to be multiplied by a factor 2.
The reason is that for ProbeG, only positive toroidal modes are calculated, to account for the negative n modes, multiply by 2. 
Unfortunatelky, there is no way to tell from an M3D-C1 output file, if ProbeG was used.
*/

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
#include <io_class.hxx>

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
	static const int Nphi = 12;			// number of toroidal angles to evaluate

	// Member Variables
	int locB;							// last finite element where B was found
	int locA;							// last finite element where A was found
	int isrc[nfiles_max];				// handle for sources
	LA_STRING filenames[nfiles_max];	// M3DC1 files
	double scale[nfiles_max];			// scaling factor for linear runs
	double RX[Nphi];					// R array of X-point at various phi
	double ZX[Nphi];					// Z array of X-point at various phi
	double RmAxis_a[Nphi];				// R array of magnetic axis at various phi
	double ZmAxis_a[Nphi];				// Z array of magnetic axis at various phi

	// Member-Functions
	void chk_linear(int response);
	int make_A(int response);
	int make_psi(void);
	int axis0(void);
	int axis(void);
	int Xpoint(double RX0 = 1.3, double ZX0 = -1.2);
	int newton2D(double& R, double phi, double& Z);
	int get_dA(double R, double phi, double Z, double& a, double& dadr, double& dadz, Array<double,1>& J);

public:
	// Member Variables
	int Np;						// public counterpart of Nphi
	int nfiles;					// actual number of files
	bool nonlinear;				// flag if run is linear or not
	int imag[nfiles_max];		// handle for magnetic fields
	int ia;						// handle for the vector potential
	double RmAxis;				// average R of magnetic axis
	double ZmAxis;				// average Z of magnetic axis
	double psi_axis_a[Nphi];	// array of poloidal flux on axis for various toroidal angles
	double psi_lcfs_a[Nphi];	// array of poloidal flux on lcfs for various toroidal angles
	double phase[nfiles_max];	// perturbation phase shift for linear runs

	// Constructors
	M3DC1();								// Default Constructor

	// Member-Functions
	int getB(double R, double phi, double Z, double& Br, double& Bp, double& Bz);
	int getA(double R, double phi, double Z, double& Ar, double& Ap, double& Az);
	int read_m3dc1sup(LA_STRING supPath="./");
	void scale_from_coils(double currents[], int loops, int turns, double adj = 1);
	void load(IO& PAR, int mpi_rank);
	void unload(void);
	int open_source(int response, int response_field, int flag = 0);
	void show_m3dc1sup_data(void);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
M3DC1::M3DC1()
{
nfiles = 1;				// default: at least one file
nonlinear = true;		// to be on the safe side
Np = Nphi;
RmAxis = 0;
ZmAxis = 0;
locA = 0;
locB = 0;
}

//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- getB -------------------------------------------------------------------------------------------------
// returns the magnetic field at location (R,phi,Z)
int M3DC1::getB(double R, double phi, double Z, double& Br, double& Bp, double& Bz)
{
int i,chk;
double coord[3], B[3];
coord[0] = R; coord[1] = phi; coord[2] = Z;

chk = 0;
Br = 0; Bp = 0; Bz = 0;
for(i=0;i<nfiles;i++)
{
	coord[1] = phi + phase[i];
	#if defined(m3dc1_newfio)
		chk += fio_eval_field(imag[i], coord, B, &locB);
	#else
		chk += fio_eval_field(imag[i], coord, B);
	#endif
	Br += B[0]; Bp += B[1]; Bz += B[2];
}
return chk;
}

// ----------------- getA -------------------------------------------------------------------------------------------------
// returns the vector potential of the first source at location (R,phi,Z)
int M3DC1::getA(double R, double phi, double Z, double& Ar, double& Ap, double& Az)
{
int chk;
double coord[3], A[3];
coord[0] = R; coord[1] = phi; coord[2] = Z;
#if defined(m3dc1_newfio)
	chk = fio_eval_field(ia, coord, A, &locA);
#else
	chk = fio_eval_field(ia, coord, A);
#endif
Ar = A[0]; Ap = A[1]; Az = A[2];
return chk;
}

// ----------------- read_m3dc1sup ----------------------------------------------------------------------------------------
// Prepare loading M3D-C1
int M3DC1::read_m3dc1sup(LA_STRING supPath)
{
freopen("/dev/null", "w", stderr);	// suppress stderr output
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

// ----------------- open_source ------------------------------------------------------------------------------------------
int M3DC1::open_source(int response, int response_field, int flag)
{
int i, ierr, ierr2, ierr3;

// open file(s)
ierr =  fio_open_source(FIO_M3DC1_SOURCE, filenames[0], &(isrc[0]));
if(ierr != 0) {cout << filenames[0] << " not found" << endl;}

// Set options appropriate to this source
chk_linear(response);
ierr2 = fio_get_options(isrc[0]);
ierr2 += fio_set_int_option(FIO_TIMESLICE, response);	// response = 1: plasma response solution; response = 0: vacuum field;

if(nonlinear) response_field = 2;	// nonlinear runs should only use the full field

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
if(ierr2 != 0) {cout << "M3DC1: setting options failed" << endl;}

// prepare vector potential
ierr3 = make_A(response);
if(ierr3 != 0) {cout << "M3DC1: Vector Potential setup failed" << endl;}

// get magnetic axis
if(RmAxis == 0 || not nonlinear) ierr3 += axis();
if(ierr3 != 0) {cout << "M3DC1: Magn. Axis setup failed" << endl;}
//ofs2 << endl << "Axis locations:" << endl;
//for(i=0;i<Nphi;i++) ofs2 << i*pi2/Nphi << "\t" << RmAxis_a[i] << "\t" << ZmAxis_a[i] << endl;

// get X-point location
//if(flag == 0 && nonlinear) ierr3 += Xpoint();
//if(ierr3 != 0) {cout << "M3DC1: Locating Xpoint failed" << endl;}
//ofs2 << endl << "X-point locations:" << endl;
//for(i=0;i<Nphi;i++) ofs2 << i*pi2/Nphi << "\t" << RX[i] << "\t" << ZX[i] << endl;

// prepare psi eval
if(flag == 0) ierr3 += make_psi();
if(ierr3 != 0) {cout << "M3DC1: poloidal flux setup failed" << endl;}

return ierr+ierr2+ierr3;
}

// ----------------- show_m3dc1sup_data -----------------------------------------------------------------------------------
void M3DC1::show_m3dc1sup_data(void)
{
int i;

cout << "M3DC1 Files found: " << nfiles << endl;
for(i=0;i<nfiles;i++)
{
	cout << "File: " << filenames[i] << ",  Scale: " << scale[i] << ",  Phase: " << phase[i] << endl;
}
}

//--------------------- Private Member Functions --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- chk_linear -------------------------------------------------------------------------------------------
void M3DC1::chk_linear(int response)
{
	if(response > 1) nonlinear = true;
	else nonlinear = false;
}

// ----------------- make_A -----------------------------------------------------------------------------------------------
int M3DC1::make_A(int response)
{
int ierr;

// Set options appropriate to this source
ierr = fio_get_options(isrc[0]);
ierr += fio_set_int_option(FIO_TIMESLICE, response);
if(not nonlinear) ierr += fio_set_int_option(FIO_PART, FIO_EQUILIBRIUM_ONLY);	// works only in a linear run, where equilibrium and perturbation are separable; in a non-linear run one has to evaluate the vector-potential at several toroidal locations and average it to get an (almost) axisymmetric psi

// set field handle
ierr += fio_get_field(isrc[0], FIO_VECTOR_POTENTIAL, &ia);

return ierr;
}

// ----------------- make_psi ---------------------------------------------------------------------------------------------
int M3DC1::make_psi(void)
{
int ierr,i;
int ipsi_axis, ipsi_lcfs;
//double a, dummy, phi;

// get psi_axis and psi_sep
// in a nonlinear run, this returns roughly a toroidally averaged value
ierr = fio_get_series(isrc[0], FIO_MAGAXIS_PSI, &ipsi_axis);
ierr += fio_get_series(isrc[0], FIO_LCFS_PSI, &ipsi_lcfs);
ierr += fio_eval_series(ipsi_axis, 0., &psi_axis_a[0]);
ierr += fio_eval_series(ipsi_lcfs, 0., &psi_lcfs_a[0]);
ierr += fio_close_series(ipsi_axis);
ierr += fio_close_series(ipsi_lcfs);

// just make sure that the other array elements are something usefull, in case they are still used elsewhere
if(nonlinear)
{
	//ierr = 0;
	for(i=1;i<Nphi;i++)
	{
		//phi = i*pi2/Nphi;
		//ierr += getA(RmAxis_a[i],phi,ZmAxis_a[i],dummy,a,dummy);
		//psi_axis_a[i] = a*RmAxis_a[i];
		psi_axis_a[i] = psi_axis_a[0];

		//ierr += getA(RX[i],phi,ZX[i],dummy,a,dummy);
		//psi_lcfs_a[i] = a*RX[i];
		psi_lcfs_a[i] = psi_lcfs_a[0];
	}
}

return ierr;
}

// ----------------- axis0 ------------------------------------------------------------------------------------------------
int M3DC1::axis0(void)
{
int ir0, iz0, ierr;

ierr = fio_get_series(isrc[0], FIO_MAGAXIS_R, &ir0);
ierr += fio_get_series(isrc[0], FIO_MAGAXIS_Z, &iz0);

ierr += fio_eval_series(ir0, 0., &RmAxis_a[0]);
ierr += fio_eval_series(iz0, 0., &ZmAxis_a[0]);

ierr += fio_close_series(ir0);
ierr += fio_close_series(iz0);
return ierr;
}

// ----------------- axis -------------------------------------------------------------------------------------------------
int M3DC1::axis(void)
{
int chk,i;
double phi;
chk = axis0();

if(nonlinear)
{
	RmAxis = RmAxis_a[0]; ZmAxis = ZmAxis_a[0];
	for(i=1;i<Nphi;i++)
	{
		RmAxis_a[i] = RmAxis_a[0];
		ZmAxis_a[i] = ZmAxis_a[0];
		phi = i*pi2/Nphi;
		chk += newton2D(RmAxis_a[i],phi,ZmAxis_a[i]);
		RmAxis += RmAxis_a[i];
		ZmAxis += ZmAxis_a[i];
	}
	RmAxis /= Nphi;
	ZmAxis /= Nphi;
}
else
{
	RmAxis = RmAxis_a[0];
	ZmAxis = ZmAxis_a[0];
}
return chk;
}

// ----------------- Xpoint -----------------------------------------------------------------------------------------------
int M3DC1::Xpoint(double RX0, double ZX0)
{
int i, chk;
double phi;
RX[0] = RX0;
ZX[0] = ZX0;
chk = newton2D(RX[0],0,ZX[0]);

for(i=1;i<Nphi;i++)
{
	RX[i] = RX[i-1];
	ZX[i] = ZX[i-1];
	phi = i*pi2/Nphi;
	chk += newton2D(RX[i],phi,ZX[i]);
}
return chk;
}

// ----------- newton2D ---------------------------------------------------------------------------------------------------
// J(1) = d2A/dR^2,  J(2) = d2A/dRdZ = J(3),  J(4) = d2A/dZ^2
int M3DC1::newton2D(double& R, double phi, double& Z)	//0: ok		-1: Fehler
{
double dr,dz,det,length;
double a, dadr, dadz, fr, fz;
int i,chk;

const int imax = 100;
const double delta = 1e-10;

// Vectors
Array<double,1> J(Range(1,4));
Array<double,1> dda(Range(1,4));

// Search
for(i=0;i<=imax;i++)
{
	chk = get_dA(R,phi,Z,a,dadr,dadz,dda);
	if(chk != 0){ofs2 << "No convergence " << chk << endl; return -1;}

	fr = a + R*dadr;
	fz = R*dadz;

	J(1) = 2*dadr + R*dda(1);
	J(4) = R*dda(4);
	J(2) = dadz + R*dda(2);
	J(3) = J(2);

	det = J(1)*J(4) - J(2)*J(3);
	dr = (J(4)*fr - J(2)*fz)/det;
	dz = (J(1)*fz - J(3)*fr)/det;

	length = sqrt(dr*dr + dz*dz);
	//if(i%20==0){cout << dr << "\t" << dt << "\t" << length <<  endl;	getchar();}
	if(length < delta)
	{
		return 0;	// convergence
	}

	R -= dr;
	Z -= dz;
}

ofs2 << "No convergence " <<  R << "\t" << Z << "\t" << dr << "\t" << dz << "\t" << length << endl;
return -1;
}

// J(1) = d2A/dR^2,  J(2) = d2A/dRdZ = J(3),  J(4) = d2A/dZ^2
int M3DC1::get_dA(double R, double phi, double Z, double& a, double& dadr, double& dadz, Array<double,1>& J)
{
const double h = 0.0001;
int chk;
double apdr,amdr,apdz,amdz;
double apdrpdz,amdrpdz,apdrmdz,amdrmdz;
double dummy;

chk  = getA(R,phi,Z,dummy,a,dummy);
chk += getA(R+h,phi,Z,dummy,apdr,dummy);
chk += getA(R-h,phi,Z,dummy,amdr,dummy);
chk += getA(R,phi,Z+h,dummy,apdz,dummy);
chk += getA(R,phi,Z-h,dummy,amdz,dummy);
chk += getA(R+h,phi,Z+h,dummy,apdrpdz,dummy);
chk += getA(R-h,phi,Z+h,dummy,amdrpdz,dummy);
chk += getA(R+h,phi,Z-h,dummy,apdrmdz,dummy);
chk += getA(R-h,phi,Z-h,dummy,amdrmdz,dummy);
if(chk != 0) return chk;

dadr = 0.5*(apdr - amdr)/h;
dadz = 0.5*(apdz - amdz)/h;

J(1) = (apdr - 2*a + amdr)/h/h;	// d2A/dR^2
J(4) = (apdz - 2*a + amdz)/h/h;	// d2A/dZ^2
J(2) = 0.25*(apdrpdz - amdrpdz - apdrmdz + amdrmdz)/h/h;	// d2A/dRdZ
J(3) = J(2);
return 0;
}

//------------------------ End of Class M3DC1 -----------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#endif //  M3DC1_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <fstream>
#include <andi.hxx>
#include <la_string.hxx>
#include <splines.hxx>
#include <blitz/array.h>

using namespace blitz;

//---------- begin class COLLISION
class COLLISION 
{	
public:
<<<<<<< HEAD
	double last_coll;  // connection length at last collision
	int num_colls;
	bool occurs(double Lc, const double x, double & meanfreepath, double & probability);
	void collide(double& R, double& Z, double modB, double x, double Lc);
	void reset();
	
// initializer	
	void init(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar, int Zq, int Mass, double rc=1, double mc=1, int mpi_rank=0); // set up collision module
=======
	bool occurs(double Lc, const double x);
	void collide(double& R, double& Z, double modB, double x);
	void init(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar, int Zq, int Mass, int mpi_rank=0); // set up collision module

>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
// default constructor
	COLLISION();

private:
// member variables
	bool use_me;
<<<<<<< HEAD
	double zeff, stdev, te_const, mfp_const, rho_const, sq2, rnum;
=======
	double zeff, stdev, te_const, mfp_const, rho_const, sq2;
	double last_coll;  // connection length at last collision
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
	double rho_coeff, mfp_coeff;
	long seed;
	
	// Temperature
	int NT;  // number of rows in T profile
	
	Array<double, 2> Tdata;  // 2 dimensional array of temperature profile data
	Array<double, 1> d2Tprofile;  // T profile spline
	
	// Density
	int ND;  // number of rows in N profile
	Array<double, 2> Ndata;  // 2 dimensional array of temperature profile data
	Array<double, 1> d2Nprofile;  // N profile spline
	
	
	
// member functions
	double meanFreePath(const double x);  // calculates and returns mean free path
	double getRho(double modB, double x);  // calculates and returns the larmor radius
	void readProfile(LA_STRING file, int& N, Array<double, 2>& Data, Array<double, 1>& d2Profile);  // reads profiles from external file, prepares splines
	double getProfile(int N, Array<double, 2>& Data, Array<double, 1>& d2Profile, const double x);  // evaluates splines of profiles
	
}; //end of class

// default constructor...
COLLISION::COLLISION() 
{
	// initialize everything to zero
	use_me = false;
	zeff = 0;
	stdev = 0;
	te_const = 0;
	mfp_const = 0;
	rho_const = 0;
	sq2 = 0;
	last_coll = 0;
	seed = 0;
	NT = 0;
	ND = 0;	
<<<<<<< HEAD
	rnum = 0;
	num_colls = 0;
=======
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
	
	rho_coeff = 0;
	mfp_coeff = 0;
}

//----------- initializer
<<<<<<< HEAD
void COLLISION::init(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar, int Zq, int Mass, double rc, double mc, int mpi_rank)
=======
void COLLISION::init(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar, int Zq, int Mass, int mpi_rank)
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
{
	use_me = true;
	last_coll = 0;
	
<<<<<<< HEAD
	rho_coeff = rc;
	mfp_coeff = 1.0/mc;
=======
	rho_coeff = 1;
	mfp_coeff = 1;
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
	
	// determine whether particle is electron or ion
	if (Zq < 0) 
	{
		rho_const = 1.07e-4;
		zeff = 1;
	} else {
		rho_const = 4.57e-3 * sqrt(Mass);
		zeff = (f + pow(zbar, 2)) / (f + zbar);  // this is actually the effective charge squared	
	}
		
	// define constants
	double m_e = 9.10938e-31;
	double e_0 = 8.854187817e-12;
	double e = 1.60217662e-19;
	
<<<<<<< HEAD
	//seed = (long) time(NULL) + mpi_rank; // add mpi_rank to give each process a unique seed
	seed = 1625981776;
	rnum = ran0(seed);
	if (mpi_rank == 0) 
	{
		std::cout << "#Using collision class: Tprofile= " << filename_te << " Nprofile= " << filename_ne << " Zeff=" << sqrt(zeff) << endl;
		std::cout << "#Seed: " << seed << endl;
		std::cout << "#rho_coeff=" << rho_coeff << " mfp_coeff=" << mfp_coeff << endl;
=======
	//seed = (long) time(NULL);
	seed = 1625616217;
	
	if (mpi_rank == 0) 
	{
		std::cout << "Using collision class: Tprofile= " << filename_te << " Nprofile= " << filename_ne << " Zeff=" << sqrt(zeff) << endl;
		std::cout << "Seed: " << seed << endl;
		std::cout << "rho_coeff=" << rho_coeff << " mfp_coeff=" << mfp_coeff << endl;
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
	}
	
	te_const = (1.09e16) / zeff;  // Wesson Tokamaks p729
	mfp_const = sqrt(1000 * e / m_e);
	sq2 = sqrt(2);
	
	// read the profiles
	readProfile(filename_te, NT, Tdata, d2Tprofile);
	readProfile(filename_ne, ND, Ndata, d2Nprofile);
<<<<<<< HEAD
	
	std::cout << "flag last_coll\tmfp\tprob\tR\tZ\tPhi\tPsi\tLc" << endl;
}

void COLLISION::reset() 
{
	last_coll = 0;
=======
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
}

//----------- mean free path
double COLLISION::meanFreePath(const double x)
{
	double Tprof = getProfile(NT, Tdata, d2Tprofile, x);
	double Nprof = getProfile(ND, Ndata, d2Nprofile, x);
		
	// constants te_const and mfp_const calculated in constructor for efficiency
	double lnA = 15.2 - (0.5)*log(Nprof) + log(Tprof);  // coulomb logarithm
	double tau_e = te_const * pow(Tprof, 1.5) / (Nprof * 1e20 * lnA);  //electron collision time
	double mfp = mfp_const * sqrt(Tprof) * tau_e;  // mean free path
	return mfp;
}

//----------- getRho
// params: modB - magnitude of magnetic field, x - flux of particle
double COLLISION::getRho(double modB, double x) 
{
	double temp = getProfile(NT, Tdata, d2Tprofile, x);
<<<<<<< HEAD
	double rho = rho_const * sqrt(temp) / modB;
	//std::cout << "Psi: " << x << endl;
	//std::cout << "Temp (kev): " << temp << endl;
	//std::cout << "Rho: " << rho << endl;
	return rho;
}

//------------ check for collision
bool COLLISION::occurs(double Lc, const double x, double & meanfreepath, double & probability) 
=======
	return rho_const * sqrt(temp) / modB;
}

//------------ check for collision
bool COLLISION::occurs(double Lc, const double x) 
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
{
	if (not use_me) return false;
	
	double l = Lc - last_coll;  // distance since last collision
<<<<<<< HEAD
	if (l == 0) 
	{
		std::cout << "# early exit" << endl;
		return false;  // don't collide on first integration step
	}
	
	double mfp = mfp_coeff * meanFreePath(x);
	meanfreepath = mfp;
	stdev = 0.2 * mfp;
	double prob = 0.5 + 0.5 * erf((l - mfp)/(stdev * sq2)); // sqrt(2) calculated in constructor for efficiency
	probability = prob;
=======
	if (l == 0) return false;  // don't collide on first integration step
	
	double mfp = mfp_coeff * meanFreePath(x);
	stdev = 0.2 * mfp;
    double rnum = ran0(seed);
	double prob = 0.5 + 0.5 * erf((l - mfp)/(stdev * sq2)); // sqrt(2) calculated in constructor for efficiency
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
	bool occured = (rnum < prob);
	
	if (occured)
	{
<<<<<<< HEAD
		std::cout << "#/----------/" << endl;
		std::cout << "#Lc: " << Lc << endl;
		std::cout << "#last_coll: " << last_coll << endl;
		std::cout << "#l: " << l << endl;
		std::cout << "#stdev: " << stdev << endl;
		std::cout << "#mean free path: " << mfp << endl;
		std::cout << "#prob: " << prob << endl;
	}
=======
		std::cout << "/----------/" << endl;
		std::cout << "Lc: " << Lc << endl;
		std::cout << "last_coll: " << last_coll << endl;
		std::cout << "l: " << l << endl;
		std::cout << "stdev: " << stdev << endl;
		std::cout << "mean free path: " << mfp << endl;
		std::cout << "prob: " << prob << endl;
	
		last_coll = Lc;
	}
	
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
	return occured;
}

//----------- collide
<<<<<<< HEAD
void COLLISION::collide(double& R, double& Z, double modB, double x, double Lc)
{
	last_coll = Lc;
	rnum = ran0(seed);	
	double theta = ran0(seed) * 2 * M_PI;
	double rho = rho_coeff * getRho(modB, x);
	std::cout << "#Angle: " << theta << " Rho: " << rho << endl;
	R += rho * cos(theta);
	Z += rho * sin(theta);
	num_colls++;
=======
void COLLISION::collide(double& R, double& Z, double modB, double x)
{
	double theta = ran0(seed) * 2 * M_PI;
	double rho = rho_coeff * getRho(modB, x);
	std::cout << "Angle: " << theta << " Rho: " << rho << endl;
	R += rho * cos(theta);
	Z += rho * sin(theta);
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
}

//----------- readProfile
// params:
<<<<<<< HEAD
// file - the file to read from, N - # columns in file, Data - 2D array for file data, 
=======
// file - the file to read from, N - # columns in file, File - 2D array for file data, 
>>>>>>> 7f9abe6e0aa4455a8df5f01b76db834000aeb6db
// Psi - 1D array for psi values, Profile - 1D array for profile values, d2Profile - array for spline
void COLLISION::readProfile(LA_STRING file, int& N, Array<double, 2>& Data, Array<double, 1>& d2Profile)
{
	double d1, dn, dpsi;
	Array <double, 1> Psi, Profile;
	
	readfile(file, 2, Data);
	N = Data.rows();
	Psi.reference(Data(Range::all(), 1));
	Profile.reference(Data(Range::all(), 2));
	
	// prepare splines
	d2Profile.resize(N);
	dpsi = Psi(2) - Psi(1);
	d1 = (Profile(2) - Profile(1)) / dpsi;
	dn = (Profile(N) - Profile(N - 1)) / dpsi;
	spline(Psi, Profile, N, d1, dn, d2Profile);	
}

//----------- getProfile
// params:
// N - # columns in file, Psi - 1D array for psi values, Profile - 1D array for profile values
// d2Profile - 1D array of splines calculated in readProfile(), x - psi for which to calculate profile
double COLLISION::getProfile(int N, Array <double, 2>& Data, Array<double, 1>& d2Profile, const double x) 
{
	Array<double, 1> Psi, Profile;
	Psi.reference(Data(Range::all(), 1));
	Profile.reference(Data(Range::all(), 2));
		
	if (x > Psi(N)) return Profile(N);
	double y, dy;
	splint(Psi, Profile, d2Profile, N, x, y, dy);
	return y;
}
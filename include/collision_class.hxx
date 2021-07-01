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
	bool occurs(double Lc, const double x);
	void collide(double& R, double& Z, double modB, double x);

// default constructor
	COLLISION();
// constructor
	COLLISION(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar); // add inputs for profiles etc.

private:
// member variables
	bool use_me;
	double zeff, stdev, te_const, mfp_const, sq2;
	double last_coll;  // connection length at last collision
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
	sq2 = 0;
	last_coll = 0;
	seed = 0;
	NT = 0;
	ND = 0;	
}

//----------- actual constructor
COLLISION::COLLISION(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar)
{
	use_me = true;
	last_coll = 0;
	zeff = (f + pow(zbar, 2)) / (f + zbar);  // this is actually the effective charge squared
	
	// define constants
	double m_e = 9.10938e-31;
	double e_0 = 8.854187817e-12;
	double e = 1.60217662e-19;
	
	seed = (long) time(NULL);
	
	te_const = (1.09e16) / zeff;  // Wesson Tokamaks p729
	mfp_const = sqrt(1000 * e / m_e);
	sq2 = sqrt(2);
	
	stdev = 10;
	
	
	// read the profiles
	readProfile(filename_te, NT, Tdata, d2Tprofile);  // rename file array..?
	readProfile(filename_ne, ND, Ndata, d2Nprofile);
	//readTprofile("profile_te");
	//readNprofile("profile_ne");
}

//----------- mean free path
// for psi=0.99, n=0.084(1e20/m^-3), T=0.165keV, lnA=14.64, tau_e=5.95e-6s, mfp=32.1m
double COLLISION::meanFreePath(const double x)
{
	double Tprof = getProfile(NT, Tdata, d2Tprofile, x);
	//std::cout << "Tprof: " << Tprof << endl;
	double Nprof = getProfile(ND, Ndata, d2Nprofile, x);
	//std::cout << "Nprof: " << Nprof << endl;
		
	// constants te_const and oosme calculated in constructor for efficiency
	double lnA = 15.2 - (0.5)*log(Nprof) + log(Tprof);  // coulomb logarithm
	//std::cout << "lnA: " << lnA << endl;
	//std::cout << "te_const: " << te_const << endl;
	double tau_e = te_const * pow(Tprof, 1.5) / (Nprof * 1e20 * lnA);  //electron collision time
	//std::cout << "tau_e: " << tau_e << endl;
	double mfp = mfp_const * sqrt(Tprof) * tau_e;  // mean free path
	return mfp;
}

//----------- getRho
// params: modB - magnitude of magnetic field, x - flux of particle
double COLLISION::getRho(double modB, double x) 
{
	double temp = getProfile(NT, Tdata, d2Tprofile, x);
	return (1.07e-4) * sqrt(temp) / modB;
}

//------------ check for collision
bool COLLISION::occurs(double Lc, const double x) 
{
	if (not use_me) return false;
	
	double l = Lc - last_coll;  // distance since last collision
	if (l == 0) return false;  // don't collide on first integration step
	
	double mfp = meanFreePath(x);
	//std::cout << "mean free path: " << mfp << endl;
	
    double rnum = ran0(seed);
	//std::cout << "rnum: " << rnum << endl;
	double prob = 0.5 + 0.5 * erf((l - mfp)/((double)stdev * sq2)); // sqrt(2) calculated in constructor for efficiency
	//std::cout << "prob: " << prob << endl;
	bool occured = (rnum < prob);
	
	if (occured) last_coll = Lc;
	
	return occured;
}

//----------- collide
void COLLISION::collide(double& R, double& Z, double modB, double x)
{
	double theta = ran0(seed) * 2 * M_PI;
	double rho = getRho(modB, x);
	R += rho * cos(theta);
	Z += rho * sin(theta);
}

//----------- readProfile
// params:
// file - the file to read from, N - # columns in file, File - 2D array for file data, 
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

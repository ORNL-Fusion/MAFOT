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
	void collide(double& R, double& Z);

// default constructor
	COLLISION();
// constructor
	COLLISION(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar); // add inputs for profiles etc.

private:
// member variables
	double zeff, stdev, m_e, e_0, e, kb, te_const, ln1e20, mfp_const, sq2;
	double last_coll;  // connection length at last collision
	long seed;
	double rho;  // larmor radius
	
	// Temperature
	
	int NT;  // number of rows in T profile
	
	Array<double, 2> Tfile;  // 2 dimensional array of temperature profile data
	Array<double, 1> d2Tprofile;  // T profile spline
	
	// Density
	int ND;  // number of rows in N profile
	Array<double, 2> Nfile;  // 2 dimensional array of temperature profile data
	Array<double, 1> d2Nprofile;  // N profile spline
	
	
	
// member functions
	double meanFreePath(const double x);  // calculates mean free path
	
	// General
	void readProfile(LA_STRING file, int& N, Array<double, 2>& File, Array<double, 1>& d2Profile);
	double getProfile(int N, Array<double, 2> File, Array<double, 1>& d2Profile, const double x);
	
}; //end of class

// default constructor...
COLLISION::COLLISION() 
{
	last_coll = 0;
}

//----------- actual constructor
COLLISION::COLLISION(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar)
{
	last_coll = 0;
	zeff = (f + pow(zbar, 2)) / (f + zbar);  // this is actually the effective charge squared
	rho = 1;
	
	// define constants
	m_e = 9.10938e-31;
	e_0 = 8.854187817e-12;
	e = 1.60217662e-19;
	kb = 1.38064852e-23;
	
	seed = (long) time(NULL);
	
	te_const = (1.09e16) / zeff;  // Wesson Tokamaks p729
	//ln1e20 = log(1e20);
	mfp_const = sqrt(1000 * e / m_e);
	sq2 = sqrt(2);
	
	stdev = 10;
	
	
	// read the profiles
	readProfile(filename_te, NT, Tfile, d2Tprofile);  // rename file array..?
	readProfile(filename_ne, ND, Nfile, d2Nprofile);
	//readTprofile("profile_te");
	//readNprofile("profile_ne");
}

//----------- mean free path
// for psi=0.99, n=0.084(1e20/m^-3), T=0.165keV, lnA=14.64, tau_e=5.95e-6s, mfp=32.1m
double COLLISION::meanFreePath(const double x)
{
	double Tprof = getProfile(NT, Tfile, d2Tprofile, x);
	std::cout << "Tprof: " << Tprof << endl;
	double Nprof = getProfile(ND, Nfile, d2Nprofile, x);
	std::cout << "Nprof: " << Nprof << endl;
		
	// constants te_const and oosme calculated in constructor for efficiency
	double lnA = 15.2 - (0.5)*log(Nprof) + log(Tprof);  // coulomb logarithm
	//std::cout << "lnA: " << lnA << endl;
	//std::cout << "te_const: " << te_const << endl;
	double tau_e = te_const * pow(Tprof, 1.5) / (Nprof * 1e20 * lnA);  //electron collision time
	//std::cout << "tau_e: " << tau_e << endl;
	double mfp = mfp_const * sqrt(Tprof) * tau_e;  // mean free path
	return mfp;
}

//------------ check for collision (eventually...)
bool COLLISION::occurs(double Lc, const double x) 
{
	double l = Lc - last_coll;  // distance since last collision
	double mfp = meanFreePath(x);
	//double mfp = 1;
	std::cout << "mean free path: " << mfp << endl;
	
    double rnum = ran0(seed);
	std::cout << "rnum: " << rnum << endl;
	double prob = 0.5 + 0.5 * erf((l - mfp)/((double)stdev * sq2)); // sqrt(2) calculated in constructor for efficiency
	std::cout << "prob: " << prob << endl;
	bool occured = (rnum < prob);
	
	if (occured) {
		std::cout << "collision\n";
		last_coll = Lc;
	} else {
		std::cout << "no collision\n";
	}
	
	// some more testing
	//std::cout << "psi:\tT:\tN:\n" << endl;
	/*for (double d = 0.0; d < 1.2; d += 0.005) 
	{
		std::cout << d << "\t" << getProfile(NT, Tfile, d2Tprofile, d) << "\t" << getProfile(ND, Nfile, d2Nprofile, d) << endl;
	}*/

	return occured;
}

//----------- collide
void COLLISION::collide(double& R, double& Z)  // add modB as input
{
	double theta = ran0(seed) * 2 * M_PI;
	R += rho * cos(theta);
	Z += rho * sin(theta);
}

// getrho

//----------- readProfile
// params:
// file - the file to read from, N - # columns in file, File - 2D array for file data, 
// Psi - 1D array for psi values, Profile - 1D array for profile values, d2Profile - array for spline
void COLLISION::readProfile(LA_STRING file, int& N, Array<double, 2>& File, Array<double, 1>& d2Profile)
{
	double d1, dn, dpsi;
	Array <double, 1> Psi, Profile;
	
	readfile(file, 2, File);
	N = File.rows();
	Psi.reference(File(Range::all(), 1));
	Profile.reference(File(Range::all(), 2));
	
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
double COLLISION::getProfile(int N, Array <double, 2> File, Array<double, 1>& d2Profile, const double x) 
{
	Array<double, 1> Psi, Profile;
	Psi.reference(File(Range::all(), 1));
	Profile.reference(File(Range::all(), 2));
		
	if (x > Psi(N)) return Profile(N);
	double y, dy;
	splint(Psi, Profile, d2Profile, N, x, y, dy);
	return y;
}

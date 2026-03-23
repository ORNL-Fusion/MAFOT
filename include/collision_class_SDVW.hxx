#include <cmath>
#include <random>
#include <utility>
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
	enum PDFModel { PDF_GAUSSIAN = 0, PDF_POISSON = 1 };

	double last_coll;  // connection length at last collision
	bool occurs(double Lc, const double x, double & meanfreepath, double & probability);  // determines probability of collision occurring, checks to see if it does
	void collide(double& R, double& Z, double B_R, double B_phi, double B_Z, double& Ekin, double x, double Lc, double mu, int& sigma);  // displaces particle to model collision with velocity space scattering (can modify sigma)
	void reset(double x);  // gets module ready for next package
	void initializeMFP(double x);  // initializes drawn_mfp with current particle position (called after particle is created)
	
// initializer	
	void init(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar, int Zq, int Mass, double x, double rc=1, double mc=1, int mpi_rank=0); // set up collision module
	void setPDFModel(PDFModel m) { pdf_model = m; }
// default constructor
	COLLISION();

private:
// member variables
	bool use_me;
	double zeff, stdev, te_const, mfp_const, rho_const, sq2, rnum;
	double rho_coeff, mfp_coeff;
	double particle_mass;  // stored for use in collide() method
	double drawn_mfp;  // randomly drawn mean free path for next collision
	long seed;
	PDFModel pdf_model; // which probability density to use when deciding collisions
	mt19937 generator;
	
	// Temperature
	int NT;  // number of rows in T profile
	Array<double, 2> Tdata;  // 2 dimensional array of temperature profile data
	Array<double, 1> d2Tprofile;  // T profile spline
	
	// Density
	int ND;  // number of rows in N profile
	Array<double, 2> Ndata;  // 2 dimensional array of temperature profile data
	Array<double, 1> d2Nprofile;  // N profile spline
	
	// Type alias for basis vectors
	using Vec3 = std::array<double,3>;

// member functions
	double meanFreePath(const double x);  // calculates and returns mean free path
	double ChandrasekharFunction(double x);  // calculates and returns the Chandrasekhar function for a given x
	double getRho(double modB, double x);  // calculates and returns the larmor radius
	void readProfile(LA_STRING file, int& N, Array<double, 2>& Data, Array<double, 1>& d2Profile);  // reads profiles from external file, prepares splines
	double getProfile(int N, Array<double, 2>& Data, Array<double, 1>& d2Profile, const double x);  // evaluates splines of profiles
	pair<double, double> draw_maxwellian(double T, double m);  // draws a random velocity from a maxwellian distribution with temperature T and mass m
	tuple<Vec3, Vec3, Vec3> generate_basis(double Br, double Bphi, double Bz);  // generates b-hat and random perpendicular unit vector
}; //end of class

// default constructor
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
	rnum = 0;
	particle_mass = 0;
	drawn_mfp = 0;
	
	rho_coeff = 0;
	mfp_coeff = 0;
}

//----------- initializer
void COLLISION::init(LA_STRING filename_te, LA_STRING filename_ne, double f, double zbar, int Zq, int Mass, double x, double rc, double mc, int mpi_rank)
{
	use_me = true;
	last_coll = 0;
	
	rho_coeff = rc;
	mfp_coeff = mc;
	
	// determine whether particle is electron or ion
	if (Zq < 0) 
	{
		rho_const = 1.07e-4;
		zeff = 1;
		particle_mass = 9.10938e-31;  // electron mass in kg
	} else {
		rho_const = 4.57e-3 * sqrt(Mass);
		zeff = (f + pow(zbar, 2)) / (f + zbar);  // this is actually the effective charge squared
		particle_mass = Mass * 1.67262192e-27;  // convert amu to kg
	}
		
	// define constants
	double m_e = 9.10938e-31;
	double e_0 = 8.854187817e-12;
	double e = 1.60217662e-19;
	double base_freq = zeff * pow(e, 2) / (4 * M_PI * pow(e_0 * particle_mass, 2)); // base frequency, Helander and Sigmar
	
	seed = (long) time(NULL) + mpi_rank; // add mpi_rank to give each process a unique seed
	//seed = 1625981776;
	generator.seed(seed);  // seed the random number generator for maxwellian draws
	uniform_real_distribution<double> U(0.0,1.0);
	rnum = U(generator);

	// select PDF model from environment variable COLLISION_PDF ("gaussian" or "poisson"). Default: gaussian
	const char* pdf_env = getenv("COLLISION_PDF");
	if(pdf_env != nullptr) {
		LA_STRING pdfs = LA_STRING(pdf_env);
		if(strcasecmp((const char*)pdfs, "poisson") == 0) pdf_model = PDF_POISSON;
		else pdf_model = PDF_GAUSSIAN;
	} else {
		pdf_model = PDF_GAUSSIAN;
	}
	
	if (mpi_rank == 0) 
	{
		std::cout << "#Using collision class: Tprofile= " << filename_te << " Nprofile= " << filename_ne << " Zeff=" << sqrt(zeff) << endl;
		std::cout << "#Seed: " << seed << endl;
		std::cout << "#rho_coeff=" << rho_coeff << " mfp_coeff=" << mfp_coeff << endl;
	}
	
	te_const = (1.09e16) / zeff;  // Wesson Tokamaks p729
	mfp_const = sqrt(1000 * e / m_e);  // Wesson Tokamaks ch14
	sq2 = sqrt(2);
	
	// read the profiles
	readProfile(filename_ne, ND, Ndata, d2Nprofile);
	readProfile(filename_te, NT, Tdata, d2Tprofile);
	
	// drawn_mfp will be initialized in reset() when particle position is available
	
	//std::cout << "flag last_coll\tmfp\tprob\tR\tZ\tPhi\tPsi\tLc" << endl;
}

//----------- reset
void COLLISION::reset(double x) 
{
	last_coll = 0;
	// Draw a new random mean free path for the next package
	double mfp = mfp_coeff * meanFreePath(x);  // use current position for mfp
	if (pdf_model == PDF_GAUSSIAN) {
		normal_distribution<double> normal_dist(mfp, 0.2 * mfp);
		do {
			drawn_mfp = normal_dist(generator);
		} while (drawn_mfp <= 0.0);
	} else {
		exponential_distribution<double> exp_dist(1.0 / mfp);
		drawn_mfp = exp_dist(generator);
	}
}

//----------- initializeMFP
// Initialize drawn_mfp when particle is created (called from particle constructor)
void COLLISION::initializeMFP(double x)
{
	double mfp = mfp_coeff * meanFreePath(x);  // use current particle position for mfp
	if (pdf_model == PDF_GAUSSIAN) {
		normal_distribution<double> normal_dist(mfp, 0.2 * mfp);
		do {
			drawn_mfp = normal_dist(generator);
		} while (drawn_mfp <= 0.0);
	} else {
		exponential_distribution<double> exp_dist(1.0 / mfp);
		drawn_mfp = exp_dist(generator);
	}
}

//----------- Decide ion or electron bath particle
bool COLLISION::isElectron(double x){
	auto [ion_mfp, elec_mfp] = meanFreePath(x);
	elec_prop = pow(ion_mfp,-1) / (pow(ion_mfp,-1) + pow(elec_mfp,-1));

	bernoulli_distribution choose_a(elec_prop);
    return choose_a(generator);
}

//----------- mean free path
pair<double,double> COLLISION::meanFreePath(const double x)
{
	double Teprof = getProfile(NT, Tdata, d2Tprofile, x); // Bath temperature in keV
	double Neprof = getProfile(ND, Ndata, d2Nprofile, x); // Bath density in 10^20 m^-3
	double Tiprof = getProfile(NT, Tdata, d2Tprofile, x); // Bath temperature in keV
	double Niprof = getProfile(ND, Ndata, d2Nprofile, x); // Bath density in 10^20 m^-3
		
	// constants te_const and mfp_const calculated in constructor for efficiency
	double lnA = 14.9 - (0.5)*log(Nprof) + log(Tprof);  // coulomb logarithm Wesson 4th ed.
	double tau_e = te_const * pow(Tprof, 1.5) / (Nprof * 1e20 * lnA);  //electron collision time
	double mfp = mfp_const * sqrt(Tprof) * tau_e;  // mean free path
	return {mfp, mfp};
}

//----------- Chandrasekhar function
double COLLISION::ChandrasekharFunction(double x)
{
	if (x < 1e-5){
		// Use asymptotic expansion for small x to avoid numerical issues
		return (2.0 * x/sqrt(M_PI)/3.0);
	} else {
		return (erf(x) - x*2/sqrt(M_PI)*exp(-x*x)) / (2.0*x*x);
	}
}

//----------- getRho
// params: modB - magnitude of magnetic field, x - flux of particle
double COLLISION::getRho(double modB, double x) 
{
	double temp = getProfile(NT, Tdata, d2Tprofile, x);
	double rho = rho_const * sqrt(temp) / modB;
	//std::cout << "Psi: " << x << endl;
	//std::cout << "Temp (kev): " << temp << endl;
	//std::cout << "Rho: " << rho << endl;
	return rho;
}

//------------ check for collision
// params: Lc - connection length, x - flux of particle, meanfreepath - mfp reference, 
// probability - prob reference (now unused but kept for interface compatibility)
bool COLLISION::occurs(double Lc, const double x, double & meanfreepath, double & probability) 
{
	if (not use_me) return false;
	
	double l = Lc - last_coll;  // distance since last collision
	if (l == 0) 
	{
		//std::cout << "# early exit" << endl;
		return false;  // don't collide on first integration step
	}
	
	// calculate mean free path at this location
	double mfp = mfp_coeff * meanFreePath(x);
	meanfreepath = mfp;
	stdev = 0.2 * mfp;
	
	// Check if distance traveled exceeds or equals the drawn mean free path
	bool occured = (l >= drawn_mfp);
	
	// For interface compatibility, set probability based on comparison
	probability = occured ? 1.0 : 0.0;
	
	/*if (occured)
	{
		std::cout << "#/----------/" << endl;
		std::cout << "#Lc: " << Lc << endl;
		std::cout << "#last_coll: " << last_coll << endl;
		std::cout << "#l: " << l << endl;
		std::cout << "#drawn_mfp: " << drawn_mfp << endl;
		std::cout << "#mfp: " << mfp << endl;
	}*/
	return occured;
}

//----------- collide with velocity space scattering
// params: R - particle R value, Z - particle Z value, B_R, B_phi, B_Z - magnetic field components,
// Ekin - particle kinetic energy (modified in place), x - particle flux, Lc - connection length,
// mu - pitch angle, sigma - particle direction (modified if velocity reverses during collision)
void COLLISION::collide(double& R, double& Z, double B_R, double B_phi, double B_Z, double& Ekin, double x, double Lc, double mu, int& sigma)
{
	last_coll = Lc;
	// Generate magnetic field basis vectors
	auto [b_hat, p_hat, q_hat] = generate_basis(B_R, B_phi, B_Z);

	// Test Particle Velocity
	double v = sqrt(2 * Ekin * 1000 * 1.60217662e-19 / particle_mass);  // convert keV to J for velocity calculation
	
	// Store old parallel velocity component sign to detect direction reversal
	double v_old_par = v * mu;
	
	// Update spatial position
	rnum = ran0(seed);	
	double theta = ran0(seed) * 2 * M_PI;  // get random angle
	double modB = sqrt(B_R*B_R + B_phi*B_phi + B_Z*B_Z);
	double rho = rho_coeff * getRho(modB, x);
	R += rho * cos(theta);
	Z += rho * sin(theta);
	
	// Velocity space scattering using temperature profile and magnetic field orientation
	double T = getProfile(NT, Tdata, d2Tprofile, x);  // temperature at this location in keV
	
	// Draw thermal velocity from Maxwellian
	auto [vpar_new, vbath_new] = draw_maxwellian(T, particle_mass);
	
	// Velocity magnitude in new frame using pitch angle mu from particle initialization
	double vnorm_perp = vbath_new * sqrt(1.0 - mu*mu);
	double vnorm_par = vbath_new * mu;
	
	// Check if parallel component reversed direction
	if ((v_old_par > 0 && vnorm_par < 0) || (v_old_par < 0 && vnorm_par > 0)) {
		sigma = -sigma;  // flip direction if velocity reversed
	}
	double dR = p_hat[0] * vnorm_perp + b_hat[0] * vnorm_par;
	double dphi = p_hat[1] * vnorm_perp + b_hat[1] * vnorm_par;
	double dZ = p_hat[2] * vnorm_perp + b_hat[2] * vnorm_par;
	
	// Update kinetic energy from thermal speed
	double e = 1.60217662e-19;
	Ekin = 0.5 * particle_mass * (dR*dR + dphi*dphi + dZ*dZ) / e / 1000.0;  // keV
	
	// Draw new random mean free path for next collision
	double mfp_local = mfp_coeff * meanFreePath(x);
	if (pdf_model == PDF_GAUSSIAN) {
		normal_distribution<double> normal_dist(mfp_local, 0.2 * mfp_local);
	do {
		drawn_mfp = normal_dist(generator);
	} while (drawn_mfp <= 0.0);
	} else {
		exponential_distribution<double> exp_dist(1.0 / mfp_local);
		drawn_mfp = exp_dist(generator);
	}
}

//----------- readProfile
// params:
// file - the file to read from, N - # columns in file, Data - 2D array for file data, 
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
	d2Profile.resize(Range(1, N));
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
	
	// evaluate splines
	if (x > Psi(N)) return Profile(N);
	double y, dy;
	splint(Psi, Profile, d2Profile, N, x, y, dy);
	return y;
}

//------------ Draw from Maxwellian to get v and vpar
// params: T - temperature in keV
// m: mass of particle
pair<double,double> COLLISION::draw_maxwellian(double T, double m)
{
	double e = 1.60217662e-19;
	double T_J = T * 1000 * e;  // convert keV to J
  	normal_distribution<double> N(0.0, sqrt(T_J / m));
	double ux = N(generator), uy = N(generator), uz = N(generator);
	double vbath = sqrt(ux*ux + uy*uy + uz*uz);
	double vpar  = ux;

	return {vpar, vbath};
}

//----------- Generate b-hat and random perpendicular unit vector
// params: Br, Bphi, Bz - components of magnetic field at location of particle
tuple<COLLISION::Vec3, COLLISION::Vec3, COLLISION::Vec3> COLLISION::generate_basis(double Br, double Bphi, double Bz)
{
    const double invB = 1.0 / sqrt(Br*Br + Bphi*Bphi + Bz*Bz);
    const double br = Br * invB, bphi = Bphi * invB, bz = Bz * invB;

	// Generate random vector r with normal distribution
    normal_distribution<double> N(0.0, 1.0);
    const double rx = N(generator), ry = N(generator), rz = N(generator);

    // Project out b: p = r - (r·b) b
    const double rdotb = rx*br + ry*bphi + rz*bz;
    double px = rx - rdotb*br;
    double py = ry - rdotb*bphi;
    double pz = rz - rdotb*bz;
    const double p2 = px*px + py*py + pz*pz;
    if (p2 < 1e-30) {
        const double rx2 = N(generator), ry2 = N(generator), rz2 = N(generator);
        const double rdotb2 = rx2*br + ry2*bphi + rz2*bz;
        px = rx2 - rdotb2*br;
        py = ry2 - rdotb2*bphi;
        pz = rz2 - rdotb2*bz;
    }
    const double invP = 1.0 / sqrt(px*px + py*py + pz*pz);
    px *= invP; py *= invP; pz *= invP;

	// Binormal: q = b x p
	double qx =  bphi*pz - bz*py;
	double qy =  bz*px   - br*pz;
	double qz =  br*py   - bphi*px;

    return { Vec3{br, bphi, bz}, Vec3{px, py, pz}, Vec3{qx, qy, qz} };
}

COLLISION::Vec3 COLLISION::delta_u(double u) // Taken from Helander and Sigmar pg. 26
{
	// Central Limit Theorem for Many Small Angle Collisions
	constexpr double sigma = M_PI / 2.0;
    normal_distribution<double> gauss(0.0, sigma);
	uniform_real_distribution<double> U(0.0,1.0);

    double x = gauss(generator);
    double y = gauss(generator);

	double alpha = sqrt(x*x + y*y);
	double azimuth = 2 * M_PI * U(generator);

	double du_x = alpha * (cos(azimuth) - 1);
	double du_y = u * sin(alpha)*cos(azimuth);
	double du_z = u * sin(alpha)*sin(azimuth);
	return Vec3{du_x, du_y, du_z};
}
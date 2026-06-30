#include <cmath>
#include <algorithm>
#include <random>
#include <utility>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <andi.hxx>
#include <la_string.hxx>
#include <splines.hxx>
#include <blitz/array.h>

using namespace blitz;

// Physical constants
namespace {
	constexpr double c    = 299792458.0;        // speed of light, m/s
	constexpr double e    = 1.60217662e-19;     // elementary charge, C
	constexpr double e_0  = 8.854187817e-12;    // vacuum permittivity, F/m
	constexpr double me   = 9.10938e-31;        // electron mass, kg
	constexpr double mp   = 1.67262192e-27;     // proton mass, kg
	constexpr double amu  = 1.66053906660e-27;  // atomic mass unit, kg
}

//---------- begin class COLLISION
class COLLISION
{
public:
	enum PDFModel { NONE = -1, PDF_GAUSSIAN = 0, PDF_POISSON = 1, NANBU = 2, PDF_NANBU = NANBU };

	PDFModel pdf_model; // which probability density to use when deciding collisions

	double last_coll;   // connection length at last collision
	double last_t_coll; // transit time at last collision [s]
	bool occurs(double Lc, const double x, double & meanfreepath, double & probability, double Ekin=0);  // determines probability of collision occurring, checks to see if it does
	void collide(double& R, double& Z, double B_R, double B_phi, double B_Z, double& Ekin, double x, double Lc, double t, double& mu, int& sigma,
	             double J_R=0.0, double J_phi=0.0, double J_Z=0.0);  // displaces particle to model collision with velocity space scattering (can modify sigma)
	void reset(double x, double Ekin);  // gets module ready for next package (Ekin in keV)
	void initializeMFP(double x, double Ekin);  // initializes drawn_mfp with current particle position (called after particle is created)
	void updateAverageTimestep(double t, bool initialize=false); // updates the shared Nanbu timestep average from orbit time
	//bool isElectron(double x, double Ekin); // decide whether to treat bath as electron or ion using local mfp

// initializer
	// init() sets up the test-particle properties and loads the center bath species into slot 0.
	// The slot-0 density profile is the electron/quasi-neutral density profile.
	// Additional bath species can be added via addBathSpecies() after init().
	void init(LA_STRING filename_te, LA_STRING filename_ne, LA_STRING filename_vt, double f, double zbar, int Zq, int Mass,
	          double rc=1, double mc=1, int mpi_rank=0,
	          double bath_mass_amu=5.4858e-4, double bath_charge_e=-1.0); // filename_vt is optional toroidal bulk velocity profile data
	void setPDFModel(PDFModel m) { pdf_model = m; }
	// Bath species management.
	// The ORDER in which species are added determines the Strang splitting:
	//   bath_species[0] is always the CENTER (full dt) and uses the electron/quasi-neutral density profile.
	//   Additional species are outer half-steps, from outermost (last added) inward.
	// mass_amu : bath particle mass in atomic mass units (use 5.486e-4 for electrons)
	// charge_e : bath particle charge in units of e (use -1 for electrons, +1 for protons, +2 for He, etc.)
	void addBathSpecies(double mass_amu, double charge_e, LA_STRING te_file, LA_STRING ne_file, int mpi_rank=0);
	void clearBathSpecies();   // remove all bath species (including the center species from init())
// default constructor
	COLLISION();

private:
// member variables
	bool use_me;
	int Zq;  // test particle charge number (stored from init() for use in PerezScatteringAngle)
	double zeff, stdev, te_const, mfp_const, rho_const, sq2, rnum, r0;
	double rho_coeff, mfp_coeff;
	double particle_mass;  // stored for use in collide() method
	double E0; // rest energy of this particle in keV (m c^2 / e / 1000)
	double drawn_mfp;  // randomly drawn mean free path for next collision
	double active_dt; // timestep currently used by Nanbu collisions
	double average_dt; // running average timestep used by Nanbu collisions
	long average_dt_count;
	long seed;
	// pdf_model moved to public section
	mt19937 generator;

	// Optional global velocity profile data.
	int VD;
	double Vconst;
	Array<double,2> Vdata;
	Array<double,1> d2Vprofile;

	// Unified bath species storage.
	// bath_species[0] = center species loaded by init(); it is always the Strang-splitting center.
	// bath_species[1..N-1] = additional species (added via addBathSpecies()), outermost-in order.
	struct BathSpecies
	{
		double mass   = 0;   // kg
		double charge = 0;   // C (negative for electrons)
		int NT = 0, ND = 0;
		Array<double,2> Tdata, Ndata;
		Array<double,1> d2Tprofile, d2Nprofile;
	};
	std::vector<BathSpecies> bath_species;
	static const std::array<double, 200> Ztable; // pre-computed table for Maxwell Juttner distribution

	// Type alias for basis vectors
	using Vec3 = mafot::Vec3;

// member functions
	double meanFreePath(const double x, double& Ekin);  // calculates and returns mean free path
	double getRho(double modB, double mu, double Ekin);  // calculates and returns the larmor radius
	void readProfile(LA_STRING file, int& N, Array<double, 2>& Data, Array<double, 1>& d2Profile);  // reads profiles from external file, prepares splines
	double getProfile(int N, Array<double, 2>& Data, Array<double, 1>& d2Profile, const double x);  // evaluates splines of profiles
	double getBulkVelocity(double x);
	// Returns T or N profile for bath species k (falls back to species 0 if k out of range)
	double getSpeciesT(int k, double x) { const int i = (k < (int)bath_species.size()) ? k : 0; return getProfile(bath_species[i].NT, bath_species[i].Tdata, bath_species[i].d2Tprofile, x); }
	double getSpeciesN(int k, double x) { const int i = (k < (int)bath_species.size()) ? k : 0; return getProfile(bath_species[i].ND, bath_species[i].Ndata, bath_species[i].d2Nprofile, x); }
	pair<double, double> draw_maxwellian(double T, double m);  // draws a random velocity from a maxwellian distribution with temperature T and mass m
	tuple<Vec3, Vec3, Vec3> generate_basis(double Br, double Bphi, double Bz);  // generates b-hat and random perpendicular unit vector
	double PerezScatteringAngle(const double x, double p1, double p2, double pstar1, double mb, double qb, double GAM1star, double GAM2star, double GAMcm, double dt, int sp_idx=0);
	double cos_chi(double s12);
	tuple<Vec3, Vec3, double> toCoMFrame(Vec3 p1, Vec3 p2, double mb);
	Vec3 CoM_pf(double chi, Vec3 p);
	tuple<double, double> coulombLog(double T, double N, double p_hat, double pTe, double Zq);
	double draw_maxwell_juttner(double T, double m);

	}; //end of class

// default constructor
COLLISION::COLLISION()
{
	use_me = false;
	Zq = 0;
	zeff = 0;
	stdev = 0;
	te_const = 0;
	mfp_const = 0;
	rho_const = 0;
	sq2 = 0;
	last_coll = 0;
	last_t_coll = 0;
	active_dt = 0;
	average_dt = 0;
	average_dt_count = 0;
	seed = 0;
	rnum = 0;
	particle_mass = 0;
	drawn_mfp = 0;
	pdf_model = NONE;
	VD = 0;
	Vconst = 0;
	bath_species.clear();
}

//----------- Bath species management
// remove all bath species; caller must add at least one before using NANBU collisions
void COLLISION::clearBathSpecies()
{
	bath_species.clear();
}

void COLLISION::addBathSpecies(double mass_amu, double charge_e, LA_STRING te_file, LA_STRING ne_file, int mpi_rank)
{
	BathSpecies bs;
	bs.mass   = mass_amu * amu;
	bs.charge = charge_e * e;
	readProfile(te_file, bs.NT, bs.Tdata, bs.d2Tprofile);
	readProfile(ne_file, bs.ND, bs.Ndata, bs.d2Nprofile);
	bath_species.push_back(std::move(bs));
	if (mpi_rank == 0)
		cout << "#Bath species " << bath_species.size()-1
		          << " (m=" << mass_amu << " amu, q=" << charge_e << " e): "
		          << "T=" << (const char*)te_file << " N=" << (const char*)ne_file << endl;
}

//----------- initializer
void COLLISION::init(LA_STRING filename_te, LA_STRING filename_ne, LA_STRING filename_vt, double f, double zbar, int Zq, int Mass,
                     double rc, double mc, int mpi_rank,
                     double bath_mass_amu, double bath_charge_e)
{
	use_me = true;
	last_coll = 0;
	last_t_coll = 0;
	active_dt = 0;
	average_dt = 0;
	average_dt_count = 0;
	this->Zq = Zq;

	rho_coeff = rc;
	mfp_coeff = mc;
	if(mfp_coeff <= 0.0)
	{
		if(mpi_rank == 0) cout << "Collision mfp_coeff must be positive." << endl;
		EXIT;
	}

	// determine whether particle is electron or ion
	if (Zq < 0)
	{
		zeff = 1;
		particle_mass = 9.10938e-31;  // electron mass in kg
	} else {
		zeff = (f + pow(zbar, 2)) / (f + zbar);  // this is actually the effective charge squared
		particle_mass = Mass * 1.67262192e-27;  // convert amu to kg
	}

	// define constants
	double base_freq = zeff * pow(e, 2) / (4 * M_PI * pow(e_0 * particle_mass, 2)); // base frequency, Helander and Sigmar

	// store particle rest energy (keV) for later use
	E0 = particle_mass*c*c/e/1000.0; 	// Rest energy in [keV]

	// relativistic gyroradius prefactor: rho = rho_const * p_hat_perp / B
	rho_const = particle_mass * c / (abs(Zq) * e);

	seed = (long) time(NULL) + mpi_rank; // add mpi_rank to give each process a unique seed

	generator.seed(seed);  // seed the random number generator for maxwellian draws
	uniform_real_distribution<double> U(0.0,1.0);
	rnum = U(generator);

	// select PDF model from environment variable COLLISION_PDF ("gaussian" or "poisson"). Default: gaussian
	const char* pdf_env = getenv("COLLISION_PDF");
	if(pdf_env != nullptr) {
		LA_STRING pdfs = LA_STRING(pdf_env);
		if(strcasecmp((const char*)pdfs, "gaussian") == 0) pdf_model = PDF_GAUSSIAN;
		else if(strcasecmp((const char*)pdfs, "nanbu") == 0) pdf_model = NANBU;
		else {
			pdf_model = PDF_POISSON;
			if(mpi_rank == 0) std::cout << "Warning: Unrecognized COLLISION_PDF environment variable value '" << pdf_env << "'. Defaulting to Poisson." << std::endl;
		}
	} else {
		pdf_model = PDF_POISSON; // default to Poisson if environment variable is not set
	}

	VD = 0;
	Vconst = 0.0;
	if (filename_vt != "")
	{
		char* end = nullptr;
		const char* text = filename_vt;
		Vconst = strtod(text, &end);
		while(end != nullptr && *end != '\0' && isspace((unsigned char)*end)) end++;
		if(!(text != end && end != nullptr && *end == '\0'))
		{
			Vconst = 0.0;
			readProfile(filename_vt, VD, Vdata, d2Vprofile);
		}
	}

	if (mpi_rank == 0)
	{
		std::cout << "#Using collision class: Tprofile= " << filename_te << " Nprofile= " << filename_ne << " Vprofile= " << filename_vt << endl;
		std::cout << " Zeff=" << sqrt(zeff) << endl;
		std::cout << "#Seed: " << seed << endl;
		std::cout << "#rho_coeff=" << rho_coeff << " mfp_coeff=" << mfp_coeff << endl;
	}

	r0 = pow(e, 2) / (4 * M_PI * e_0 * particle_mass * pow(c, 2)); // classical radius
	te_const = 4 * M_PI * c * pow(r0, 2);  // Hesslow 2018
	sq2 = sqrt(2);

	// Load the Strang-splitting center (slot 0); its density profile is the electron/quasi-neutral density.
	bath_species.clear();
	{
		BathSpecies bs;
		bs.mass   = bath_mass_amu * amu;
		bs.charge = bath_charge_e * e;
		readProfile(filename_te, bs.NT, bs.Tdata, bs.d2Tprofile);
		readProfile(filename_ne, bs.ND, bs.Ndata, bs.d2Nprofile);
		bath_species.push_back(std::move(bs));
		if (mpi_rank == 0)
			std::cout << "#Bath species 0 (Strang center): m=" << bath_mass_amu << " amu, q=" << bath_charge_e
			          << " e, T=" << filename_te << " N=" << filename_ne << std::endl;
	}

	// drawn_mfp will be initialized in reset() when particle position is available

	//std::cout << "flag last_coll\tmfp\tprob\tR\tZ\tPhi\tPsi\tLc" << endl;
}

//----------- reset
void COLLISION::reset(double x, double Ekin)
{
	last_coll = 0;
	last_t_coll = 0;

	if (pdf_model != PDF_NANBU) {
		// Draw a new random mean free path for the next package
		double mfp = mfp_coeff * meanFreePath(x, Ekin);  // use current position for mfp
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
}

//----------- initializeMFP
// Initialize drawn_mfp when particle is created (called from particle constructor)
void COLLISION::initializeMFP(double x, double Ekin)
{
	last_coll = 0;
	last_t_coll = 0;
	double mfp = mfp_coeff * meanFreePath(x, Ekin);  // use current particle position for mfp
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

//----------- updateAverageTimestep
void COLLISION::updateAverageTimestep(double t, bool initialize)
{
	if(initialize)
	{
		last_t_coll = t;
		if(average_dt > 0.0) active_dt = average_dt;
		return;
	}

	double dt = t - last_t_coll;
	last_t_coll = t;
	if(dt <= 0.0) return;

	average_dt_count++;
	average_dt += (dt - average_dt)/double(average_dt_count);
}

//----------- Coulomb logarithm
tuple<double, double> COLLISION::coulombLog(double T, double N, double p_hat, double pTe, double Zq) {
	double lnA = 14.9 - 0.5 * log(N) + log(T);  // coulomb logarithm Wesson 4th ed.

	if (Zq < 0){
		// corrected coulomb logarithm for electron-electron and electron-ion collisions Hesslow 2018
		return { lnA + 0.2 * log(1 + pow((2 * p_hat / pTe), 5)), lnA + 0.2 * log(1 + pow((2 * p_hat / pTe), 5)) };
	} else {
		// Even energetic ions have small gamma
		return { lnA, lnA };
	}
}

//----------- mean free path
// Sums collision frequencies across all bath species.
double COLLISION::meanFreePath(const double x, double& Ekin)
{
	const double eps   = Ekin / E0;
	const double p_hat = sqrt(eps * (eps + 2.0));
	const double GAMMA = 1.0 + eps;
	const double v     = c * p_hat / GAMMA;

	double nu_total = 0.0;
	for (int k = 0; k < (int)bath_species.size(); k++)
	{
		double Tprof = getSpeciesT(k, x);  // keV
		double Nprof = getSpeciesN(k, x);  // 10^20 m^-3
		const double pTk = sqrt(Tprof / E0);
		double lnA   = 14.9 - 0.5 * log(Nprof) + log(Tprof); // coulomb logarithm Wesson 4th ed.
		double lnA_k = lnA + 0.2 * log(1 + pow((2 * p_hat / pTk), 5)); // corrected coulomb logarithm Hesslow 2018
		nu_total += te_const * GAMMA / pow(p_hat, 3) * (Nprof * 1e20 * lnA_k); // electron-ion collision frequency Hesslow 2018
	}
	if (nu_total == 0.0) return 0.0;
	return v / nu_total;
}

//----------- getRho
// params: modB - magnitude of magnetic field, x - flux of particle, mu - pitch angle, Ekin - particle kinetic energy
double COLLISION::getRho(double modB, double mu, double Ekin)
{
	// use stored rest energy E0 (keV)
	const double eps = Ekin / E0;
	const double p_hat = sqrt(eps * (eps + 2.0));
	const double p_perp_hat = p_hat * sqrt(max(0.0, 1.0 - mu*mu));
	return rho_const * p_perp_hat / modB;
}

//------------ boost frame
mafot::Vec3 boost(const mafot::Vec3& p0, const mafot::Vec3& b_hat, double beta_d)
{
    beta_d = std::max(-0.999, std::min(0.999, beta_d));
    const double gamma_d = 1.0 / sqrt(1.0 - beta_d*beta_d);

    const double gamma0 = sqrt(1.0 + blitz::dot(p0, p0));
    const double ppar0 = blitz::dot(p0, b_hat);
    const mafot::Vec3 pperp0 = p0 - ppar0 * b_hat;

    const double ppar_lab = gamma_d * (ppar0 + beta_d * gamma0);
    return pperp0 + ppar_lab * b_hat;
}

//------------ check for collision
// params: Lc - connection length, x - flux of particle, meanfreepath - mfp reference,
// probability - prob reference (now unused but kept for interface compatibility)
bool COLLISION::occurs(double Lc, const double x, double& meanfreepath, double& probability, double Ekin)
{
	if (not use_me) return false;

	double l = Lc - last_coll;  // distance since last collision
	if (l == 0)
	{
		//std::cout << "# early exit" << endl;
		return false;  // don't collide on first integration step
	}

	// calculate mean free path at this location
	double mfp = mfp_coeff * meanFreePath(x, Ekin);
	meanfreepath = mfp;
	stdev = 0.2 * mfp;

	// Check if distance traveled exceeds or equals the drawn mean free path
	bool occured = (l >= drawn_mfp);

	// For interface compatibility, set probability based on comparison
	probability = occured ? 1.0 : 0.0;

	return occured;
}

//----------- collide with velocity space scattering
// params: R - particle R value, Z - particle Z value, B_R, B_phi, B_Z - magnetic field components,
// Ekin - particle kinetic energy (modified in place), x - particle flux, Lc - connection length,
// mu - pitch angle cosine, sigma - particle direction (modified if velocity reverses during collision)
void COLLISION::collide(double& R, double& Z, double B_R, double B_phi, double B_Z, double& Ekin, double x, double Lc, double t, double& mu, int& sigma,
                        double J_R, double J_phi, double J_Z)
{
	double dt = t - last_t_coll;  // elapsed time since last collision [s]
	last_coll = Lc;
	if (pdf_model == PDF_NANBU) {
		if (active_dt > 0.0) dt = active_dt;
		else last_t_coll = t;
	}
	else last_t_coll = t;
	if (dt <= 0.0) return;  // should not happen, but guard PerezScatteringAngle against zero dt
	uniform_real_distribution<double> U(-1.0,1.0); // Bath pitch angle isotropic in 3D, so uniform in cos(theta)
	// Velocity space scattering using temperature profile and magnetic field orientation
	if (pdf_model == PDF_NANBU)
	{
		// Generate magnetic field basis vectors
		auto [b_hat, p_hat_basis, q_hat] = generate_basis(B_R, B_phi, B_Z);
		const Vec3 J_vec(J_R, J_phi, J_Z);
		const Vec3 ion_velocity_vec(0.0, getBulkVelocity(x), 0.0);
		const double ne = getSpeciesN(0, x) * 1e20; // slot 0 carries the electron/quasi-neutral density profile [m^-3]

		// Normalized test particle momentum magnitude
		double eps = Ekin / E0;
		double p_hat_mag = sqrt(eps * (eps + 2.0));
		double p_old_par = p_hat_mag * mu;
		double p_perp_par = p_hat_mag * sqrt(max(0.0, 1.0 - mu*mu));
		double phi_rand = U(generator) * M_PI;  // random gyrophase
		Vec3 p_vec = p_old_par * b_hat + p_perp_par * (cos(phi_rand) * p_hat_basis + sin(phi_rand) * q_hat);

		// Draw bath momenta for all species from their respective T profiles
		const int N_bath = (int)bath_species.size();
		std::vector<Vec3> pbath(N_bath);
		for (int k = 0; k < N_bath; k++)
		{
			double T_k    = getSpeciesT(k, x);
			double pb_mag = draw_maxwell_juttner(T_k, bath_species[k].mass);
			double mu_k   = U(generator);
			double az_k   = U(generator) * M_PI;
			Vec3 p_bath_rest = pb_mag * mu_k * b_hat
			                  + pb_mag * sqrt(max(0.0, 1.0 - mu_k*mu_k))
			                    * (cos(az_k) * p_hat_basis + sin(az_k) * q_hat);

			Vec3 drift_velocity(0.0, 0.0, 0.0);
			if (bath_species[k].charge > 0.0)
			{
				drift_velocity = ion_velocity_vec;
			}
			else if (bath_species[k].charge < 0.0 && ne > 0.0)
			{
				// EFIT J is conventional current. The optional bulk profile fixes
				// the ion velocity gauge; electron drift follows from J = e ne (vi - ve).
				drift_velocity = ion_velocity_vec - J_vec/(e * ne);
			}
			const double drift_speed = mafot::norm(drift_velocity);
			pbath[k] = (drift_speed == 0.0) ? p_bath_rest : boost(p_bath_rest, drift_velocity/drift_speed, drift_speed/c);
		}

		// ---- Telescoping Strang splitting ----
		// bath_species[0] is ALWAYS the center (full dt).
		// Outer species (indices N-1 down to 1) each get dt/2 as symmetric half-steps.
		// For a single species this reduces to just: species[0](dt).
		// For two species [A, B]: B(dt/2)  A(dt)  B(dt/2)
		// For three species [A, B, C]: C(dt/2) B(dt/2) A(dt) B(dt/2) C(dt/2)

		// Helper: one Perez scatter step
		auto scatter_step = [&](Vec3 p_in, int k, double step_dt) -> Vec3
		{
			const Vec3& p_bath = pbath[k];
			const double mb    = bath_species[k].mass;
			const double qb    = bath_species[k].charge;
			auto [p1cm, Bcm, GAMcm] = toCoMFrame(p_in, p_bath, mb);
			double p1cm_mag = sqrt(blitz::dot(p1cm,p1cm));
			double GAMMA1   = sqrt(1.0 + blitz::dot(p_in,p_in));
			double GAMMA2b  = sqrt(1.0 + blitz::dot(p_bath,p_bath));
			double GAM1star = GAMcm * (GAMMA1 - blitz::dot(Bcm, p_in));
			double GAM2star = GAMcm * (GAMMA2b - blitz::dot(Bcm, p_bath));
			double chi = PerezScatteringAngle(x, sqrt(blitz::dot(p_in,p_in)), sqrt(blitz::dot(p_bath,p_bath)),
				p1cm_mag, mb, qb, GAM1star, GAM2star, GAMcm, step_dt, k);
			Vec3 CMpf = CoM_pf(chi, p1cm);
			return CMpf + Bcm * ((GAMcm*GAMcm) / (GAMcm + 1) * blitz::dot(CMpf, Bcm) + GAM1star * GAMcm);
		};

		// Forward outer half-steps: outermost (N-1) down to 1
		Vec3 pf = p_vec;
		for (int k = N_bath - 1; k >= 1; k--)
			pf = scatter_step(pf, k, 0.5*dt);
		// Center: species 0, full dt
		pf = scatter_step(pf, 0, dt);
		// Reverse outer half-steps: 1 up to N-1
		for (int k = 1; k < N_bath; k++)
			pf = scatter_step(pf, k, 0.5*dt);

		// Update kinetic energy from thermal speed
		Ekin = E0 * (sqrt(1 + blitz::dot(pf,pf)) - 1);  // keV

		// Update pitch angle cosine
		mu = blitz::dot(pf, b_hat) / sqrt(blitz::dot(pf,pf));

		// Update gyrocenter
		double modB = sqrt(B_R*B_R + B_phi*B_phi + B_Z*B_Z);
		Vec3 p_par_vec = blitz::dot(pf, b_hat) * b_hat;
		Vec3 p_perp_vec = pf - p_par_vec;
		Vec3 p_perp_old_vec = p_vec - blitz::dot(p_vec, b_hat) * b_hat;
		Vec3 drho = sqrt(rho_coeff) * c*particle_mass/(Zq*e) * mafot::cross(p_perp_vec - p_perp_old_vec, b_hat) / modB;
		R += drho[0];
		Z += drho[2];

		// Check if parallel component reversed direction
		double p_new_par = blitz::dot(pf, b_hat);
		if ((p_new_par > 0 && p_old_par < 0) || (p_new_par < 0 && p_old_par > 0)) {
			sigma = -sigma;  // flip direction if velocity reversed
		}
	} else {
		uniform_real_distribution<double> U(-1.0,1.0);

		// Basic pitch angle scattering model (no energy change)
		mu = U(generator);  // randomized pitch angle at each step
		int sign = (mu > 0) - (mu < 0);  // sign of mu for direction of scattering
		sigma = sign;  // update particle direction based on new pitch angle
		double theta = U(generator) * M_PI;  // get random angle
		double rho = sqrt(rho_coeff) * getRho(sqrt(B_R*B_R + B_phi*B_phi + B_Z*B_Z), mu, Ekin);
		R += rho * cos(theta);
		Z += rho * sin(theta);

		// Draw new random mean free path for next collision
		double mfp_local = mfp_coeff * meanFreePath(x, Ekin);
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
}

//----------- readProfile
// params:
// file - the file to read from, N - # columns in file, Data - 2D array for file data,
// Psi - 1D array for psi values, Profile - 1D array for profile values, d2Profile - array for spline
void COLLISION::readProfile(LA_STRING file, int& N, Array<double, 2>& Data, Array<double, 1>& d2Profile)
{
	double d1, dn, dpsi;
	Array <double, 1> Psi, Profile;
	vector<pair<double, double>> rows;
	ifstream in((const char*)file);
	if(in.fail()==1) {cout << "Unable to open file " << file << endl; exit(0);}

	string line;
	bool extraColumnsWarned = false;
	while(getline(in, line))
	{
		size_t first = line.find_first_not_of(" \t\r\n");
		if(first == string::npos || line[first] == '#') continue;

		istringstream iss(line);
		vector<double> values;
		double value;
		while(iss >> value)
		{
			if(isnan(value))
			{
				cout << "profiles have nan" << endl;
				exit(0);
			}
			values.push_back(value);
		}
		if(values.size() > 2 && !extraColumnsWarned)
		{
			cout << values.size() << " columns present, reading only the first two" << endl;
			extraColumnsWarned = true;
		}
		if(values.size() >= 2) rows.push_back({values[0], values[1]});
	}
	in.close();

	N = rows.size();
	if(N < 2) {cout << "Profile file needs at least two data rows: " << file << endl; exit(0);}

	Data.resize(Range(1,N), Range(1,2));
	for(int i=1; i<=N; i++)
	{
		Data(i,1) = rows[i-1].first;
		Data(i,2) = rows[i-1].second;
	}
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
	if (x < Psi(1)) return Profile(1);
	if (x > Psi(N)) return Profile(N);
	double y, dy;
	splint(Psi, Profile, d2Profile, N, x, y, dy);
	return y;
}

double COLLISION::getBulkVelocity(double x)
{
	if(VD > 0) return getProfile(VD, Vdata, d2Vprofile, x);
	return Vconst;
}

// ---------- Helper to convert to center of momentum frame. mb=bath particle mass
tuple<COLLISION::Vec3, COLLISION::Vec3, double> COLLISION::toCoMFrame(Vec3 p1, Vec3 p2, double mb) {
	// Lorentz factor for test particle and bath particle
	double GAMMA1 = sqrt(1.0 + blitz::dot(p1,p1));
	double GAMMAb = sqrt(1.0 + blitz::dot(p2,p2));
	double GAMrel = GAMMA1 * GAMMAb - blitz::dot(p1, p2);

	// CM velocity and Lorentz factor
	double Mcm2 = particle_mass * particle_mass + mb * mb + 2 * mb * particle_mass * GAMrel;
	double Mcm = sqrt(Mcm2);
	double mGAM = particle_mass * GAMMA1 + mb * GAMMAb;
	double GAMcm = mGAM / Mcm;
	Vec3 Bcm = (particle_mass * p1 + mb * p2) / (mGAM);  // CM Beta

	// Lorentz transform normalized test particle momentum to CM frame
	Vec3 p1cm = p1 + ((GAMcm*GAMcm)/(GAMcm + 1) * blitz::dot(Bcm, p1) - GAMcm * GAMMA1) * Bcm;
	return {p1cm, Bcm, GAMcm};
}

//----------- Determine deflection angle according to Perez/Nanbu model
// params: p1 magnitude of test particle momentum, p2 magnitude of bath particle momentum both in lab frame
// mb bath particle mass
// qb bath particle charge
// GAMstar Lorentz factor in CoM frame for test particle
double COLLISION::cos_chi(double s12) {
	uniform_real_distribution<double> U(0.0,1.0);

	// Inversion sampling Perez 2012, eqn 11
	if (s12 < 0.1) {
		return 1.0 + s12 * log(U(generator));
	} else if (s12 < 3.0) {
		double s = s12;
		double s2 = s * s;
		double s3 = s2 * s;
		double s4 = s3 * s;
		double s5 = s4 * s;

		double Ainv =
			0.0056958
		+ 0.9560202 * s
		- 0.508139  * s2
		+ 0.47913906 * s3
		- 0.12788975 * s4
		+ 0.02389567 * s5;

		double A = 1.0 / Ainv;

		return (1.0 / A) * log(exp(-A) + 2.0 * U(generator) * sinh(A));

	} else if (s12 < 6.0) {
		double A = 3.0 * exp(-s12);
		return (1.0 / A) * log(exp(-A) + 2.0 * U(generator) * sinh(A));

	} else {
		return 2.0 * U(generator) - 1.0;
	}
}

double COLLISION::PerezScatteringAngle(const double x, double p1, double p2,
	double pstar1, double mb, double qb, double GAM1star, double GAM2star,
	double GAMcm, double dt, int sp_idx)
{
	double Tprof = getSpeciesT(sp_idx, x);  // keV for bath species sp_idx
	double Nprof = getSpeciesN(sp_idx, x);  // 10^20 m^-3 for bath species sp_idx

	double GAM1 = sqrt(1.0 + p1*p1);
	double GAM2 = sqrt(1.0 + p2*p2);
	double mGAM = particle_mass * GAM1 + mb * GAM2;
	double emFactor = 4 * M_PI * e_0*e_0 * c*c*c*c * particle_mass * GAM1 * mb * GAM2;
	double pstar1_phys = pstar1 * particle_mass * c;  // convert to physical momentum for use in scattering parameter

	// Calculate Coulomb logarithm, two distinct values for electron test, only base for ion test
	auto [lnA_i, lnA_e] = coulombLog(Tprof, Nprof, p1, sqrt(Tprof / E0), Zq);

	// Calculate scattering parameter (Perez 2012, eqn 9)
	double stb = Nprof*1e20 * dt * Zq*e*qb*Zq*e*qb / (emFactor) * GAMcm * pstar1_phys / mGAM
				* (particle_mass * GAM1star * mb * GAM2star / (pstar1_phys*pstar1_phys) * c*c +1)
				* (particle_mass * GAM1star * mb * GAM2star / (pstar1_phys*pstar1_phys) * c*c +1);

	if (qb < 0) {
		stb = stb * lnA_e;
	} else {
		stb = stb * lnA_i;
	}
	stb /= mfp_coeff;
	double cchi = cos_chi(stb);
	cchi = max(-1.0, min(1.0, cchi));
	return acos(cchi);
}

//----------- Generate b-hat and random perpendicular unit vector not used for basic pitch angle scattering model
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
	const Vec3 b_hat = Vec3(br, bphi, bz);
	const Vec3 p_hat = Vec3(px, py, pz);
	const Vec3 q_hat = mafot::cross(b_hat, p_hat);

    return { b_hat, p_hat, q_hat };
}

COLLISION::Vec3 COLLISION::CoM_pf(double chi, Vec3 p) // Taken from Perez et. al 2012
{
	// Isotropic azimuth
	uniform_real_distribution<double> U(0.0,1.0);
	double azimuth = 2 * M_PI * U(generator);

	// Momentum components
	double p1x = p[0], p1y = p[1], p1z = p[2];
	double pmag = sqrt(p1x*p1x + p1y*p1y + p1z*p1z);

	// Trig arguments
	double sc_cp = sin(chi)*cos(azimuth);
	double sc_sp = sin(chi)*sin(azimuth);
	double cc = cos(chi);
	double p_perp = sqrt(p1x*p1x + p1y*p1y);

	// Rare singular case: incoming CM momentum is nearly along z.
    if (p_perp < 1e-14 * pmag) {
        return Vec3(
            pmag * sc_cp,
            pmag * sc_sp,
            p1z * cc
        );
    }

	// Momentum update in CoM frame. Perez et. al 2012
	double dp_x = p1x*p1z/p_perp*sc_cp - p1y*pmag/p_perp*sc_sp + p1x*cc;
	double dp_y = p1y*p1z/p_perp*sc_cp + p1x*pmag/p_perp*sc_sp + p1y*cc;
	double dp_z = -p_perp*sc_cp + p1z*cc;
	return Vec3(dp_x, dp_y, dp_z);
}

//------------ Draw from Maxwellian to get v and vpar
// params: T - temperature in keV
// m: mass of particle
double COLLISION::draw_maxwell_juttner(double T, double m)
{
	if(!isfinite(T) || !isfinite(m) || T <= 0.0 || m <= 0.0)
	{
		cout << "Maxwell-Juttner sampler received invalid input" << endl;
		exit(0);
	}
	double T_J = T * 1000 * e;  // convert keV to J

	uniform_real_distribution<double> U(0.0,1.0);
	gamma_distribution<double> maxwell(1.5, 2*T_J / m);

	double E0m = m*c*c/e/1000.0; 	// Rest energy in [keV]
	double tt = T / E0m; // dimensionless temperature
	double dY = 0.19095477386934673; // Y step size for interpolation in log-log space
	double gamma1 = 0.0; // initialization for while loop, more than 0 for safety
	double Z = 0; // initialize Z for interpolation
	const int max_attempts = 1000000;

	if (tt < 0.1) {
		// Non-relativistic Maxwellian: v^2 ~ Gamma(3/2, 2T/m)
		double v2 = maxwell(generator);

		for(int attempt=0; v2 >= c*c && attempt<max_attempts; attempt++) {
			v2 = maxwell(generator);
		}
		if(v2 >= c*c) {cout << "Maxwell-Juttner sampler failed" << endl; exit(0);}

		const double beta2 = v2 / (c*c);

		return sqrt(beta2 / (1.0 - beta2));
	}

	for(int attempt=0; attempt<max_attempts; attempt++) {
		double X = log1p(-U(generator)) - 1/tt + log1p(1/tt+0.5/tt/tt);
		if(!isfinite(X) || X >= 0.0) continue;

		if (log(-X) < -26.0){
			Z = log(-6*X)/3.0; // X -> 0^- limit
		}else if (log(-X) > 12.0){
			Z = log(-X + 2*log(-X)); // X -> -infinity limit
		} else {
			double Y = log(-X);
			double t = (Y+26.0)/dY;
			int i = int(t);

			if (i < 0) i = 0;
			else if (i > 198) i = 198; // keep i in bounds for interpolation

			double frac = t - i;
			Z = (1.0 -frac) * COLLISION::Ztable[i] + frac * COLLISION::Ztable[i+1];
		}

		double u0 = exp(Z);

		// one Newton step
		double Hu = -u0 + log(1.0 + u0 + 0.5*u0*u0);
		double Hp = -(u0*u0) / (2.0 + 2.0*u0 + u0*u0);

		gamma1 = (u0 - (Hu - X) / Hp) * tt;
		if(!isfinite(gamma1) || gamma1 <= 1.0) continue;
		if (U(generator) < sqrt((1-1/gamma1)*(1+1/gamma1))) {
			return sqrt((gamma1-1.0)*(gamma1+1.0));
		}
	}
	cout << "Maxwell-Juttner sampler failed" << endl;
	exit(0);
}

const std::array<double, 200> COLLISION::Ztable = {
	// Table for Z = log(u) used in interpolation step for Maxwell-Juttner
	-8.0693338890372797e+00,
	-8.0056853949234590e+00,
	-7.9420254247720914e+00,
	-7.8783615739995927e+00,
	-7.8147060375542567e+00,
	-7.7510479473089031e+00,
	-7.6873895420476979e+00,
	-7.6237304371830383e+00,
	-7.5600692751948495e+00,
	-7.4964094194694457e+00,
	-7.4327486440437243e+00,
	-7.3690873211492045e+00,
	-7.3054263368963372e+00,
	-7.2417634114379217e+00,
	-7.1781006233835098e+00,
	-7.1144358916468571e+00,
	-7.0507715476406556e+00,
	-6.9871055058730454e+00,
	-6.9234388008275713e+00,
	-6.8597711059868542e+00,
	-6.7961023163631546e+00,
	-6.7324322582037821e+00,
	-6.6687609989533874e+00,
	-6.6050885739203373e+00,
	-6.5414146986940391e+00,
	-6.4777394236236470e+00,
	-6.4140625286551201e+00,
	-6.3503840693343898e+00,
	-6.2867037655391611e+00,
	-6.2230215601534136e+00,
	-6.1593373972378211e+00,
	-6.0956510541061482e+00,
	-6.0319624329597232e+00,
	-5.9682713633593512e+00,
	-5.9045777253022740e+00,
	-5.8408812985980063e+00,
	-5.7771819330526206e+00,
	-5.7134794053629934e+00,
	-5.6497735413256418e+00,
	-5.5860641041881420e+00,
	-5.5223508507295822e+00,
	-5.4586335424229295e+00,
	-5.3949119020207510e+00,
	-5.3311856512534552e+00,
	-5.2674544829435979e+00,
	-5.2037180719284457e+00,
	-5.1399760734165136e+00,
	-5.0762281172647352e+00,
	-5.0124738086607330e+00,
	-4.9487127278550025e+00,
	-4.8849444289777813e+00,
	-4.8211684331178413e+00,
	-4.7573842305949130e+00,
	-4.6935912789630772e+00,
	-4.6297889970212545e+00,
	-4.5659767678723426e+00,
	-4.5021539303714810e+00,
	-4.4383197815916162e+00,
	-4.3744735699266073e+00,
	-4.3106144936877042e+00,
	-4.2467416977016699e+00,
	-4.1828542688003925e+00,
	-4.1189512324530044e+00,
	-4.0550315484781763e+00,
	-3.9910941059665861e+00,
	-3.9271377191848189e+00,
	-3.8631611217330031e+00,
	-3.7991629611690110e+00,
	-3.7351417928599040e+00,
	-3.6710960735477194e+00,
	-3.6070241544097641e+00,
	-3.5429242735354300e+00,
	-3.4787945476895490e+00,
	-3.4146329639579776e+00,
	-3.3504373700737857e+00,
	-3.2862054645689218e+00,
	-3.2219347858253617e+00,
	-3.1576227003975537e+00,
	-3.0932663903679685e+00,
	-3.0288628397534696e+00,
	-2.9644088197681993e+00,
	-2.8999008728892930e+00,
	-2.8353352957009896e+00,
	-2.7707081202417414e+00,
	-2.7060150938807013e+00,
	-2.6412516574666531e+00,
	-2.5764129216668046e+00,
	-2.5114936413228834e+00,
	-2.4464881875782774e+00,
	-2.3813905176466186e+00,
	-2.3161941419605356e+00,
	-2.2508920884572001e+00,
	-2.1854768637372364e+00,
	-2.1199404108157678e+00,
	-2.0542740631119796e+00,
	-1.9884684943616298e+00,
	-1.9225136640457894e+00,
	-1.8563987579239087e+00,
	-1.7901121232254817e+00,
	-1.7236411979947248e+00,
	-1.6569724340672436e+00,
	-1.5900912131015110e+00,
	-1.5229817550517371e+00,
	-1.4556270184276221e+00,
	-1.3880085916446561e+00,
	-1.3201065747410721e+00,
	-1.2518994506891397e+00,
	-1.1833639455284479e+00,
	-1.1144748765202352e+00,
	-1.0452049875354323e+00,
	-9.7552477091432455e-01,
	-9.0540227509219939e-01,
	-8.3480289737670654e-01,
	-7.6368916140742171e-01,
	-6.9202047902791786e-01,
	-6.1975289658132138e-01,
	-5.4683882601454414e-01,
	-4.7322676166761990e-01,
	-3.9886098425598615e-01,
	-3.2368125436387929e-01,
	-2.4762249875697168e-01,
	-1.7061449407544310e-01,
	-9.2581553965416868e-02,
	-1.3442227503326652e-02,
	6.6890981105981309e-02,
	1.4851185913568921e-01,
	2.3152068306040940e-01,
	3.1602434917892519e-01,
	4.0213641569321179e-01,
	4.8997703156000499e-01,
	5.7967272436402240e-01,
	6.7135601691476876e-01,
	7.6516484064057289e-01,
	8.6124171366801949e-01,
	9.5973265333315216e-01,
	1.0607857974217343e+00,
	1.1645497163074363e+00,
	1.2711714098564344e+00,
	1.3807939987472817e+00,
	1.4935541395706089e+00,
	1.6095792160252107e+00,
	1.7289843833619327e+00,
	1.8518695678989439e+00,
	1.9783165453007352e+00,
	2.1083862373686530e+00,
	2.2421163743597887e+00,
	2.3795196658576310e+00,
	2.5205826065413492e+00,
	2.6652650138798770e+00,
	2.8135003545669863e+00,
	2.9651968688204824e+00,
	3.1202394511738194e+00,
	3.2784921983915845e+00,
	3.4398014947441449e+00,
	3.6039994761879637e+00,
	3.7709077005020797e+00,
	3.9403408506917050e+00,
	4.1121103126778848e+00,
	4.2860274926466646e+00,
	4.4619067707276105e+00,
	4.6395680219580679e+00,
	4.8188386691529104e+00,
	4.9995552624826782e+00,
	5.1815646053753106e+00,
	5.3647244648874972e+00,
	5.5489039168785705e+00,
	5.7339833827307292e+00,
	5.9198544159528632e+00,
	6.1064192949112082e+00,
	6.2935904732778276e+00,
	6.4812899335906726e+00,
	6.6694484824107754e+00,
	6.8580050185769394e+00,
	7.0469057994388473e+00,
	7.2361037239723425e+00,
	7.4255576464893691e+00,
	7.6152317302958723e+00,
	7.8050948470997401e+00,
	7.9951200251639056e+00,
	8.1852839470471075e+00,
	8.3755664961812268e+00,
	8.5659503504016250e+00,
	8.7564206197878800e+00,
	8.9469645257053330e+00,
	9.1375711176965702e+00,
	9.3282310247973754e+00,
	9.5189362378992985e+00,
	9.7096799199117427e+00,
	9.9004562406633490e+00,
	1.0091260233701437e+01,
	1.0282087672384556e+01,
	1.0472934962902935e+01,
	1.0663799052098078e+01,
	1.0854677348178248e+01,
	1.1045567652639297e+01,
	1.1236468101896349e+01,
	1.1427377117311627e+01,
	1.1618293362466384e+01,
	1.1809215706670679e+01,
	1.2000143193835402e+01
};

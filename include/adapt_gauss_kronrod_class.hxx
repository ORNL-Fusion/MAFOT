// Class solves integral though adaptive stepsize Gauss-Kronrod quadrature
// Based off Quadpack++
// A.Wingen						26.01.15

/* Gauss-Kronrod quadrature rule(s) plus error-estimate data for
 adaptive quadrature routines.

 The Gauss-Kronrod abscissae consist of 2m+1 points x_1 < ... < x_{2m+1}
 in the interval (-1, 1) used for a low-order and a high-order quadrature rule:
 Q_m^G f = \sum_{k=1}^m a_k f(x_{2k}), Q_m^{GK} f = \sum_{k=1}^{2m+1} b_k f(x_k).

 The weights and abscissae are available through member functions, however
 they are stored according to compact QUADPACK convention. Due to symmetry,
 the positive abscissae x_{2m+1}, x_{2m},..., x_m are returned
 as values xgk(0), xgk(1),..., xgk(m+1) respectively
 of the member function xgk(). Note the reverse order. The corresponding
 weights b_{2m+1}, b_{2m},..., b_m are returned by the
 respective values of wgk(). The weights a_m, a_{m-1},...,
 corresponding to the even-indexed x_{2m}, x_{2m-2}, ...., are
 returned by the values of wg() in their reverse order.

 Computational details:

 The even-indexed abscissae x_2, ..., x_{2m} are the zeros of the m-th
 Legendre polynomial P_m. The odd indexed points are zeros of a polynomial
 that is represented as a Chebyshev sum,
 E_{m+1} = T_{m+1} + c_{m-1}T_{m-1} + c_{m-3} T_{m-3} + ...,
 whose coefficients are defined by explicit formulae in
 - Giovanni Monegato, "Some remarks on the construction of extended Gaussian
 quadrature rules," Math. Comp., Vol. 32 (1978) pp. 247-252.
 http://www.jstor.org/stable/2006272.

 The zeros of both of these polynomials are computed by Newton's method. Upper
 bounds for their round-off errors, as functions of machine epsilon, are
 incorporated in the stopping criteria for the root finders.

 The weights a_1, ..., a_m are Gauss-Legendre weights.
 The b_k are given by the formulae
 b_{2k} = a_k + \frac{2 p_m}{(2m+1)t_{m+1} P_m'(x_{2k}) E_{m+1}(x_{2k})}, k = 1,...,m,
 and
 b_{2k+1} = \frac{2 p_m}{(2m+1) t_{m+1} P_m(x_{2k+1}) E_{m+1}'(x_{2k+1})}, k = 0,..., m,
 where p_m and t_{m+1} are the leading coefficients of the
 polynomials P_m and T_{m+1} respectively. These are from
 - Giovanni Monegato, "A note on extended Gaussian quadrature rules,"
 Math. Comp., Vol. 30 (1976) pp. 812-817.
 http://www.jstor.org/stable/2005400
 */


// Define
//--------
#ifndef ADPAT_GAUSS_KRONROD_CLASS_INCLUDED
#define ADPAT_GAUSS_KRONROD_CLASS_INCLUDED

// Include
//--------
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <andi.hxx>
using namespace blitz;

// Prototypes
//-----------
extern ofstream logfile;

// Golbal Parameters
//------------------
const double EPSILON = 2.22045e-16;	// Machine epsilon
const double XMIN = 2.22507e-308;	// Underflow threshold



//--------- Begin Class FtnBase -------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
/* Basis for a function to get passed into the algorithm
User-defined functions for numerical routines may require additional parameters
depending on the application. However function-pointers, when passed as
arguments to C-routines, depend on their "signature" type. To accomodate a
universal signature, only the base type, FtnBase, is passed as an argument to routines.

A user defined class which inherits FtnBase and overloads FtnBase::operator() can then be used */
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class FtnBase
{
private:
	Array<double,1> out;
public:
	FtnBase(int N = 1) { out.resize(N); };
	virtual Array<double,1> operator() (double x) { out = 0; return out; };
};
//------------------------ End of Class FtnBase ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//--------- Begin Class GaussKronrod --------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Solves the integral int_a^b \vec f(x) dx
// for all N components of \vec f simultaneously
// using Gauss-Kronrod quadrature.
// The weights and abscissae are precalculated for grad (2m+1) Legendre polynoms
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class GaussKronrod
{
private:
	// Member Variables
	int N_;	 // Number of integrals to solve simultaneously
	int m_;  // Gauss-Legendre degree
	int n_;  // size of Gauss-Kronrod arrays
	Array<double,1> xgk_;  // Gauss-Kronrod abscissae
	Array<double,1> wg_;   // Gauss-Legendre weights
	Array<double,1> wgk_;  // Gauss-Kronrod weights

	Array<double,1> coefs;  // Chebyshev coefficients
	Array<double,1> zeros;  // zeros of Legendre polynomial

	Array<double,2> fv1, fv2;  // scratch space for error estimator

	// Member-Functions
	double rescale_error(double err, const double result_abs,const double result_asc);

	void legendre_zeros();
	void chebyshev_coefs();
	void gauss_kronrod_abscissae();
	void gauss_kronrod_weights();

	double legendre_err(int deg, double x, double& err);
	double legendre_deriv(int deg, double x);
	double chebyshev_series(double x, double& err);
	double chebyshev_series_deriv(double x);


public:
	// Member Variables

	// Constructors
	GaussKronrod(int m = 10, int N = 1);

	// Member-Operators

	// Member-Functions
	void qk(FtnBase& f, double a, double b, Array<double,1>& result, Array<double,1>& abserr, Array<double,1>& resabs, Array<double,1>& resasc);
	int size() { return n_; };	// Size of arrays of Gauss-Kronrod abscissae and weights.
	double xgk(int k) { return (0 <= k && k < n_) ? xgk_(k) : 0.0; };	// Array of Gauss-Kronrod abscissae in (0, 1); QUADPACK convention.
	double wgk(int k) { return (0 <= k && k < n_) ? wgk_(k) : 0.0; };	// Array of corresponding Gauss-Kronrod weights; QUADPACK convention.
	double wg(int k) { return (0 <= k && k < n_/2) ? wg_(k) : 0.0; };	// Gauss-Legendre weights for odd indexed abscissae; QUADPACK convention.
}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
GaussKronrod::GaussKronrod(int m, int N)
{
	N_ = N;
	m_ = m;
	n_ = m_ + 1;
	xgk_.resize(n_);
	wg_.resize(n_ / 2);
	wgk_.resize(n_);
	coefs.resize(n_ + 1);
	zeros.resize(m_ + 2);
	fv1.resize(n_, N_);
	fv2.resize(n_, N_);

	legendre_zeros();
	chebyshev_coefs();
	gauss_kronrod_abscissae();
	gauss_kronrod_weights();
}


//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// evaluate the integral of f(x) on the interval [a,b]
void GaussKronrod::qk(FtnBase& f, double a, double b, Array<double,1>& result, Array<double,1>& abserr, Array<double,1>& resabs, Array<double,1>& resasc)
{
	int jtw, jtwm1;
	double abscissa;
	const double center = (a + b) / 2;
	const double half_length = (b - a) / 2;
	const double abs_half_length = fabs(half_length);
	Array<double,1> f_center(N_), result_gauss(N_), result_kronrod(N_), result_abs(N_), result_asc(N_), fval1(N_), fval2(N_), fsum(N_), err(N_), err_scaled(N_);
	Range all = Range::all();

	f_center = f(center);
	result_gauss = 0.0;
	result_kronrod = f_center * wgk_(n_ - 1);
	result_abs = fabs(result_kronrod);
	result_asc = 0.0; err = 0.0;

	int j;

	if(n_ % 2 == 0) result_gauss = f_center * wg_(n_/2 - 1);

	for(j = 0; j < (n_ - 1) / 2; j++)
	{
		jtw = j * 2 + 1;        // j=1,2,3 jtw=2,4,6
		abscissa = half_length * xgk_(jtw);
		fval1 = f(center - abscissa);
		fval2 = f(center + abscissa);
		fsum = fval1 + fval2;
		fv1(jtw,all) = fval1;
		fv2(jtw,all) = fval2;
		result_gauss += wg_(j) * fsum;
		result_kronrod += wgk_(jtw) * fsum;
		result_abs += wgk_(jtw) * (fabs(fval1) + fabs(fval2));
	}

	for(j = 0; j < n_ / 2; j++)
	{
		jtwm1 = j * 2;
		abscissa = half_length * xgk_(jtwm1);
 		fval1 = f(center - abscissa);
		fval2 = f(center + abscissa);
		fv1(jtwm1,all) = fval1;
		fv2(jtwm1,all) = fval2;
		result_kronrod += wgk_(jtwm1) * (fval1 + fval2);
		result_abs += wgk_(jtwm1) * (fabs(fval1) + fabs(fval2));
	}

	result_asc = wgk_(n_ - 1) * fabs(f_center - 0.5*result_kronrod);

	for(j = 0; j < n_ - 1; j++)
	{
		result_asc += wgk_(j) * (fabs(fv1(j,all) - 0.5*result_kronrod) + fabs(fv2(j,all) - 0.5*result_kronrod));
	}

	// scale by the width of the integration region
	err = (result_kronrod - result_gauss) * half_length;

	result_kronrod *= half_length;
	result_abs *= abs_half_length;
	result_asc *= abs_half_length;
	for(j=0;j<N_;j++) err_scaled(j) = rescale_error(err(j), result_abs(j), result_asc(j));

	result.resize(N_);
	resabs.resize(N_);
	resasc.resize(N_);
	abserr.resize(N_);
	result = result_kronrod;
	resabs = result_abs;
	resasc = result_asc;
	abserr = err_scaled;
}

//--------------------- Private Member Functions --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------------------
// Computes the zeros of the Legendre polynomial P_m. Upon exit, these
// are stored consecutively as elements of the array zeros[] indexed by 1,..., m.
void GaussKronrod::legendre_zeros()
{
	Array<double,1> temp(m_+1);
	zeros(0) = -1.0;
	zeros(1) = 1.0;
	double delta, epsilon;

	for(int k = 1; k <= m_; ++k)
	{
		// loop to locate zeros of P_k interlacing z_0,...,z_k
		for(int j = 0; j < k; ++j)
		{
			// Newton's method for P_k :
			// initialize solver at midpoint of (z_j, z_{j+1})
			delta = 1;
			double x_j = (zeros(j) + zeros(j+1)) / 2;
			double P_k = legendre_err(k, x_j, epsilon);
			while(fabs(P_k) > epsilon && fabs(delta) > EPSILON)
			{
				delta = P_k / legendre_deriv(k, x_j);
				x_j -= delta;
				P_k = legendre_err(k, x_j, epsilon);
			}
			temp(j) = x_j;
		}

		// copy roots tmp_0,...,tmp_{k-1} to z_1,...z_k:
		zeros(k+1) = zeros(k);
		for(int j = 0; j < k; ++j)
			zeros(j+1) = temp(j);

	}
}


//-------------------------------------------------------------------------------------------------------------------------
// Computes coefficients of polynomial E_{m+1} in the array coefs[].
void GaussKronrod::chebyshev_coefs()
{
	int ell = (m_ + 1)/2;
	Array<double,1> alpha(ell+1);
	Array<double,1> f(ell+1);

	// Care must be exercised in initalizing the constants in the definitions.
	// Compilers interpret expressions like "(2*k + 1.0)/(k + 2.0)" as floating
	// point precision, before casting to double.
	f(1) = double(m_+1)/double(2*m_ + 3);
	alpha(0) = 1.0; // coefficient of T_{m+1}
	alpha(1) = -f(1);

	for(int k = 2; k <= ell; ++k)
	{
		f(k) = f(k-1) * (2*k - 1) * (m_ + k) / (k * (2*m_ + 2*k + 1));
		alpha(k) = -f(k);
		for(int i = 1; i < k; ++i)
			alpha(k) -= f(i) * alpha(k-i);
	}

	for(int k = 0; k <= ell; ++k)
	{
		coefs(m_ + 1 - 2*k) = alpha(k);
		if(m_ >= 2*k) coefs(m_ - 2*k) = 0.0;
	}
}


//-------------------------------------------------------------------------------------------------------------------------
// Computes Gauss-Legendre weights wg_[] and Gauss-Kronrod weights wgk_[].
void GaussKronrod::gauss_kronrod_weights()
{
	double err;
	// Gauss-Legendre weights:
	for(int k = 0; k < n_ / 2; ++k)
	{
		double x = xgk_(2*k + 1);
		wg_(k) = (double(-2) / ((m_ + 1) * legendre_deriv(m_, x) * legendre_err(m_+1, x, err)));
	}

	// The ratio of leading coefficients of P_n and T_{n+1} is computed from the recursive formulae for the respective polynomials.
	double F_m = 2.0 / double(2*m_ + 1);
	for(int k = 1; k <= m_; ++k)
		F_m *= (double(2*k) / double(2*k - 1));

	// Gauss-Kronrod weights:
	for(int k = 0; k < n_; ++k)
	{
		double x = xgk_(k);
		if (k % 2 == 0) wgk_(k) = F_m / (legendre_err(m_, x, err) * chebyshev_series_deriv(x));
		else wgk_(k) = (wg_(k/2) + F_m / (legendre_deriv(m_, x) * chebyshev_series(x, err)));
	}
}


//-------------------------------------------------------------------------------------------------------------------------
// Computes the zeros of the polynomial E_{m+1}, using the fact that these
// interlace the zeros of the Legendre polynomial P_m, which are stored in
// the array zeros[]. Appropriate elements of zeros[] are then copied into xgk_[].
void GaussKronrod::gauss_kronrod_abscissae()
{
	double delta, epsilon;

	for(int k = 0; k < n_ / 2; ++k)
	{
		delta = 1;
		// Newton's method for E_{n+1} :
		double x_k = (zeros(m_-k) + zeros(m_+1-k))/2.0;
		double E = chebyshev_series(x_k, epsilon);
		while(fabs(E) > epsilon && fabs(delta) > EPSILON)
		{
			delta = E / chebyshev_series_deriv(x_k);
			x_k -= delta;
			E = chebyshev_series(x_k, epsilon);
		}
		xgk_(2*k) = x_k;
		// copy adjacent Legendre-zero into the array:
		if(2*k+1 < n_) xgk_(2*k+1) = zeros(m_-k);
	}
}


//-------------------------------------------------------------------------------------------------------------------------
// Recursive definition of the Legendre polynomials (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x),
// and estimate of the rounding error, E_{k+1} = \frac{(2k+1)|x|E_k + kE_{k-1}}{2(k+1)},
// are from the routine gsl_sf_legendre_Pl_e distributed with GSL.
double GaussKronrod::legendre_err(int n, double x, double& err)
{
	if(n == 0)
	{
		err = 0.0;
		return 1.0;
	}
	else if(n == 1)
	{
		err = 0.0;
		return x;
	}

	double P0 = 1.0, P1 = x, P2;
	double E0 = EPSILON;
	double E1 = fabs(x) * EPSILON;
	for(int k = 1; k < n; ++k)
	{
		P2 = ((2*k + 1) * x * P1 - k * P0) / (k + 1);
		err = ((2*k + 1) * fabs(x) * E1 + k * E0) / (2*(k + 1));
		P0 = P1; P1 = P2;
		E0 = E1; E1 = err;
	}
	return P2;
}


//-------------------------------------------------------------------------------------------------------------------------
// Three-term recursion identity for the Legendre derivatives:
// P_{k+1}'(x) = (2k+1) P_k(x) + P_{k-1}'(x).
double GaussKronrod::legendre_deriv(int n, double x)
{
	if(n == 0) return 0.0;
	else if(n == 1) return 1.0;

	double P0 = 1.0, P1 = x, P2;
	double dP0 = 0.0, dP1 = 1.0, dP2;
	for(int k = 1; k < n; ++k)
	{
		P2 = ((2*k + 1) * x * P1 - k * P0) / double(k + 1);
		dP2 = (2*k + 1) * P1 + dP0;
		P0 = P1; P1 = P2;
		dP0 = dP1; dP1 = dP2;
	}
	return dP2;
}


//-------------------------------------------------------------------------------------------------------------------------
// Evaluation of the polynomial E_{m+1} is using the Clenshaw method and
// (truncation) error estimate is taken from the routine
// gsl_cheb_eval_err distributed with GSL.
double GaussKronrod::chebyshev_series(double x, double& err)
{
	double d1(0), d2(0);
	double absc = fabs(coefs(0)); // final term for truncation error
	double y2 = 2 * x; // linear term for Clenshaw recursion

	for(int k = n_; k >= 1; --k)
	{
		double temp = d1;
		d1 = y2 * d1 - d2 + coefs(k);
		d2 = temp;
		absc += fabs(coefs(k));
	}

	err = absc * EPSILON;
	return x * d1 - d2 + coefs(0)/2;
}


//-------------------------------------------------------------------------------------------------------------------------
// Derivatives of Chebyshev polynomials satisfy the identity T_n' = nU_{n-1},
// where the U_k are Chebyshev polynomials of the second kind. The derivative
// of the polynomial E_{m+1} is implemented using the Clenshaw algorithm for
// the latter polynomials.
double GaussKronrod::chebyshev_series_deriv(double x)
{
	double d1(0), d2(0);
	double y2 = 2 * x; // linear term for Clenshaw recursion

	for(int k = n_; k >= 2; --k)
	{
		double temp = d1;
		d1 = y2 * d1 - d2 + k * coefs(k);
		d2 = temp;
	}

	return y2 * d1 - d2 + coefs(1);
}


//-------------------------------------------------------------------------------------------------------------------------
// QUADPACK's nonlinear formula for the absolute error.
double GaussKronrod::rescale_error(double err, const double result_abs, const double result_asc)
{
	err = fabs(err);

	if(result_asc != 0.0 && err != 0.0)
	{
		double scale = pow((200 * err / result_asc), 1.5);

		if(scale < 1.0) err = result_asc * scale ;
		else err = result_asc ;
	}

	if(result_abs > XMIN / (50 * EPSILON))
	{
		double min_err = 50 * EPSILON * result_abs ;
		if (min_err > err) err = min_err ;
	}
	return err ;
}

//------------------------ End of Class GaussKronrod ----------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



//--------- Begin Class AdaptiveGK ----------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// Solves the integral int_a^b \vec f(x) dx
// by adaptively refining the interval [a,b]
// using the Gauss-Kronrod quadrature on each subinterval
// until the overall error estimate is below the tolerances epsabs or epsrel.
// Since \vec f is an array, the refinement uses the largest of all errors
// within a subinterval as criterion.
// The results and errors of a subinterval are stored in lists. The lists are sorted,
// based on the largest max error. The subinterval on top gets refined. Then repeat...
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class AdaptiveGK : public GaussKronrod
{
private:
	// Member Variables
	int N_;	 // Number of integrals to solve simultaneously
	int limit;		// maximum number of refinements
	int size;
	int nrmax;
	int i_work;
	int maximum_level;
	int done_iter;

	Array<double,1> alist;
	Array<double,1> blist;
	Array<double,1> elist;
	Array<int,1> order;
	Array<int,1> level;

	Array<double,2> rlist;
	Array<double,2> EList;

	// Member-Functions
	void  initialise(double a, double b);
	void  set_initial_result(Array<double,1> result, Array<double,1> error);
	void  qpsrt();
	void  retrieve(double& a, double& b, Array<double,1>& r, Array<double,1>& e);
	int   subinterval_too_small (double a1, double a2, double b2);
	Array<double,1>  sum_results();
	void  update(double a1, double b1, Array<double,1>& area1, Array<double,1>& error1, double a2, double b2, Array<double,1>& area2, Array<double,1>& error2);

public:
	// Member Variables

	// Constructors
	AdaptiveGK(int limit = 1, int m = 10, int N = 1);	// Initialize for limit refinement intervals and (2m+1)-point quadrature.

	// Member-Operators

	// Member-Functions
	int qag(FtnBase& f, double a, double b, double epsabs, double epsrel, Array<double,1>& result, double& abserr);
	int qag(FtnBase& f, double a, double b, double epsabs, double epsrel, double& result, double& abserr);	// Overload for special case N_ = 1
	int get_iterations();

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Standard Constructor
AdaptiveGK::AdaptiveGK(int n, int m, int N) : GaussKronrod(m, N)
{
	N_ = N;
	if (n == 0) n = 1;	// ensure that workspace has positive size
	alist.resize(n);
	blist.resize(n);
	elist.resize(n);
	order.resize(n);
	level.resize(n);

	rlist.resize(n,N);
	EList.resize(n,N);

	size = 0;
	limit = n;
	maximum_level = 0;
	done_iter = 0;
}


//--------------------- Public Member Functions ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//------------------ qag --------------------------------------------------------------------------------------------------
// computes the overall integral by consecutively bisecting [a,b] and iterating through all subintervals
int AdaptiveGK::qag(FtnBase& f, double a, double b, double epsabs, double epsrel, Array<double,1>& result, double& abserr)
{
	Array<double,1> area(N_), result0(N_), abserr0(N_), resabs0(N_), resasc0(N_), errsum(N_), tolerance(N_), round_off(N_);
	int iteration = 0, error_type = 0;
	Array<int,1> roundoff_type1(N_), roundoff_type2(N_);
	roundoff_type1 = 0; roundoff_type2 = 0;

	// Initialize results
	initialise (a, b);

	result.resize(N_);
	result = 0.0;
	abserr = 0.0;

	if(epsabs <= 0.0 && epsrel < 50 * EPSILON)
	{
		logfile << "tolerance cannot be acheived with given epsabs and epsrel" << endl;
		return -1;
	}

	// perform the first integration
	this->qk(f, a, b, result0, abserr0, resabs0, resasc0);
	set_initial_result(result0, abserr0);

	// Test on accuracy
	for(int i = 0;i<N_;i++) tolerance(i) = max(epsabs, epsrel * fabs(result0(i)));

	// need IEEE rounding here to match original quadpack behavior
	round_off = 50 * EPSILON * resabs0;

	if(any(abserr0 <= round_off && abserr0 > tolerance))
	{
		result = result0;
		abserr = max(abserr0);

		logfile << "cannot reach tolerance because of roundoff error on first attempt" << endl;
		return -1;
	}
	else if(all((abserr0 <= tolerance && abserr0 != resasc0) || abserr0 == 0.0))
	{
		result = result0;
		abserr = max(abserr0);

		return 0;
	}
	else if(limit == 1)
	{
		result = result0;
		abserr = max(abserr0);

		logfile << "a maximum of one iteration was insufficient" << endl;
		return -1;
	}

	area = result0;
	errsum = abserr0;
	iteration = 1;

	double a1, b1, a2, b2;
	double a_i, b_i;
	Array<double,1> area1(N_), area2(N_), area12(N_), r_i(N_), delta(N_), e_i(N_);
	Array<double,1> error1(N_), error2(N_), error12(N_);
	Array<double,1> resasc1(N_), resasc2(N_);
	Array<double,1> resabs1(N_), resabs2(N_);

	do
	{
		// Bisect the subinterval with the largest error estimate
		retrieve(a_i, b_i, r_i, e_i);

		a1 = a_i;
		b1 = (a_i + b_i) / 2.0;
		a2 = b1;
		b2 = b_i;

		this->qk(f, a1, b1, area1, error1, resabs1, resasc1);
		this->qk(f, a2, b2, area2, error2, resabs2, resasc2);

		area12 = area1 + area2;
		error12 = error1 + error2;

		errsum += error12 - e_i;
		area += area12 - r_i;

		for(int i=0;i<N_;i++)
		{
			if(resasc1(i) != error1(i) && resasc2(i) != error2(i))
			{
				delta(i) = r_i(i) - area12(i);

				if(fabs(delta(i)) <= 1.0e-5 * fabs(area12(i)) && error12(i) >= 0.99 * e_i(i))
				{
					roundoff_type1(i)++;
				}
				if(iteration >= 10 && error12(i) > e_i(i)) roundoff_type2(i)++;
			}
			tolerance(i) = max(epsabs, epsrel * fabs(area(i)));
		}

		if(any(errsum > tolerance))
		{
			if(any(roundoff_type1 >= 6 || roundoff_type2 >= 20)) error_type = 2;   // round off error

			//set error flag in the case of bad integrand behaviour at a point of the integration range
			if(subinterval_too_small(a1, a2, b2)) error_type = 3;
		}

		update(a1, b1, area1, error1, a2, b2, area2, error2);
		iteration++;
	}
	while(iteration < limit && !error_type && any(errsum > tolerance));

	result = sum_results();
	abserr = max(errsum);
	done_iter += iteration;

	if(all(errsum <= tolerance)) return 0;
	else if(error_type == 2)
	{
		logfile << "roundoff error prevents tolerance from being achieved" << endl;
		return -1;
	}
	else if(error_type == 3)
	{
		logfile << "bad integrand behavior found in the integration interval" << endl;
		return -1;
	}
	else if(iteration == limit)
	{
		logfile << "maximum number of subdivisions reached" << endl;
		return -1;
	}
	else
	{
		logfile << "could not integrate function" << endl;
		return -1;
	}
}

// --- overload for special case N_ = 1 -------
int AdaptiveGK::qag(FtnBase& f, double a, double b, double epsabs, double epsrel, double& result, double& abserr)
{
	if(N_ > 1)
	{
		logfile << "AdaptiveGK: Result must be an array, not a double. Check your class constructor and member function input -> Abort" << endl;
		return -1;
	}
	int idx;
	Array<double,1> result0(N_);
	idx = qag(f, a, b, epsabs, epsrel, result0, abserr);
	result = result0(0);
	return idx;
}


//------------------ get_iterations ---------------------------------------------------------------------------------------
// return the done_iter member
int AdaptiveGK::get_iterations()
{
	return done_iter;
}

//--------------------- Private Member Functions --------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//------------------ initialise -------------------------------------------------------------------------------------------
// initializes the lists
void AdaptiveGK::initialise(double a, double b)
{
	Range all = Range::all();
	size = 0;
	nrmax = 0;
	i_work = 0;
	alist = 0.0;
	blist = 0.0;
	rlist = 0.0;
	EList = 0.0;
	elist = 0.0;
	order = 0;
	level = 0;

	alist(0) = a;
	blist(0) = b;

	maximum_level = 0;
}


//------------------ qpsrt -----------------------------------------------------------------------------------------------
// sort the lists
void AdaptiveGK::qpsrt()
{
	const int last = size - 1;

	double errmax;
	double errmin;
	int i, k, top;

	int i_nrmax = nrmax;
	int i_maxerr = order(i_nrmax);

	// Check whether the list contains more than two error estimates

	if(last < 2)
	{
		order(0) = 0;
		order(1) = 1;
		i_work = i_maxerr;
		return;
	}

	errmax = elist(i_maxerr);

	/* This part of the routine is only executed if, due to a difficult
	 integrand, subdivision increased the error estimate. In the normal
	 case the insert procedure should start after the nrmax-th largest
	 error estimate. */
	while(i_nrmax > 0 && errmax > elist(order(i_nrmax - 1)))
	{
		order(i_nrmax) = order(i_nrmax - 1);
		i_nrmax--;
	}

	/* Compute the number of elements in the list to be maintained in
	 descending order. This number depends on the number of
	 subdivisions still allowed. */
	if(last < (limit/2 + 2)) top = last;
	else top = limit - last + 1;

	// Insert errmax by traversing the list top-down, starting comparison from the element elist(order(i_nrmax+1)).
	i = i_nrmax + 1;

	// The order of the tests in the following line is important to prevent a segmentation fault
	while(i < top && errmax < elist(order(i)))
	{
		order(i-1) = order(i);
		i++;
	}
	order(i-1) = i_maxerr;

	// Insert errmin by traversing the list bottom-up
	errmin = elist(last);
	k = top - 1;

	while(k > i - 2 && errmin >= elist(order(k)))
	{
		order(k+1) = order(k);
		k--;
	}
	order(k+1) = last;

	// Set i_max and e_max
	i_maxerr = order(i_nrmax);
	i_work = i_maxerr;
	nrmax = i_nrmax;
}


//--------------------- set_initial_result -------------------------------------------------------------------------------
// load first result into lists
void AdaptiveGK::set_initial_result(Array<double,1> result, Array<double,1> error)
{
	Range all = Range::all();
	size = 1;
	rlist(0,all) = result;
	EList(0,all) = error;
	elist(0) = max(error);
}


//---------------------- retrieve -----------------------------------------------------------------------------------------
// get results and errors from lists
void AdaptiveGK::retrieve(double& a, double& b, Array<double,1>& r, Array<double,1>& e)
{
	Range all = Range::all();
	a = alist(i_work);
	b = blist(i_work);
	r.reference(rlist(i_work,all));
	e.reference(EList(i_work,all));
	//e = elist(i_work);
}


//---------------------- sum_results --------------------------------------------------------------------------------------
// sum the results from all subintervals
Array<double,1> AdaptiveGK::sum_results()
{
	int k;
	Array<double,1> result_sum(N_);
	result_sum = 0.0;
	Range all = Range::all();

	for(k = 0; k < size; k++)
	{
      result_sum += rlist(k,all);
	}
	return result_sum;
}


//--------------------- subinterval_too_small -----------------------------------------------------------------------------
// check if refinement reached system accuracy
int AdaptiveGK::subinterval_too_small(double a1, double a2, double b2)
{
	double tmp = (1 + 100 * EPSILON) * (fabs(a2) + 1000 * XMIN);
	int status = (fabs(a1) <= tmp && fabs(b2) <= tmp);
	return status;
}


//------------------------- update ----------------------------------------------------------------------------------------
// load refined results into lists
void AdaptiveGK::update(double a1, double b1, Array<double,1>& area1, Array<double,1>& error1, double a2, double b2, Array<double,1>& area2, Array<double,1>& error2)
{
	const int i_max = i_work;
	const int i_new = size;
	const int new_level = level(i_max) + 1;
	Range all = Range::all();
	double e1 = max(error1);
	double e2 = max(error2);

	// append the newly-created intervals to the list
	if (e2 > e1)
	{
		alist(i_max) = a2;        // blist(maxerr) is already == b2
		rlist(i_max,all) = area2;
		elist(i_max) = e2;
		level(i_max) = new_level;
		EList(i_max,all) = error2;

		alist(i_new) = a1;
		blist(i_new) = b1;
		rlist(i_new,all) = area1;
		elist(i_new) = e1;
		level(i_new) = new_level;
		EList(i_new,all) = error1;
	}
	else
	{
		blist(i_max) = b1;        // alist(maxerr) is already == a1
		rlist(i_max,all) = area1;
		elist(i_max) = e1;
		level(i_max) = new_level;
		EList(i_max,all) = error1;

		alist(i_new) = a2;
		blist(i_new) = b2;
		rlist(i_new,all) = area2;
		elist(i_new) = e2;
		level(i_new) = new_level;
		EList(i_new,all) = error2;
	}
	size++;
	if (new_level > maximum_level) maximum_level = new_level;
	qpsrt ();
}

//------------------------ End of Class AdaptiveGK ------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------






#endif //  ADPAT_GAUSS_KRONROD_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

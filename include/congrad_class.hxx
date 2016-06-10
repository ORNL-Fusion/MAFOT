// Class minimizes Ax - b using conjugate gradient method
// A.Wingen						1.4.16


// Define
//--------
#ifndef CONGRAD_CLASS_INCLUDED
#define CONGRAD_CLASS_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <andi.hxx>
using namespace blitz;

// Prototypes
//-----------
inline double FMAX(double a, double b) {return (a > b ? a : b);}
inline double SIGN(double a, double b) {return (b >= 0.0 ? fabs(a) : -fabs(a));}

// Golbal Parameters
//------------------


//A user defined class which inherits congradFtn and overloads congradFtn::operator() can then be used */
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class congradFtn
{
private:
	Array<double,1> out;
public:
	congradFtn(int N = 1) { out.resize(N); };
	virtual Array<double,1> operator() (Array<double,1>& x) { out = 0; return out; };
};

class congradFtn_f
{
private:
	double out;
public:
	congradFtn_f() { out = 0; };
	virtual double operator() (Array<double,1>& x) { out = 0; return out; };
};
//------------------------ End of Class congradFtn ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//--------- Begin Class congrad -------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// finds array x, that minimizes A*x - b for given operator A and array b
// using a conjugate gradient iteration method
// IMPORTANT: A needs to be symmetric and positive definite
// A trick could be, solving A^T * A * x = A^T * b
// or A * A^T * y = b with then x = A^T * y
// First trick works ok (not great, but ok). Second one lacks a proper initial guess for y
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class congrad
{
private:
	// Member Variables
	int nb;
	Array<double,1> p;		// search direction, intermediate variable for conjugate gradient method, same shape as b
	Array<double,1> v;		// A*p, intermediate variable for conjugate gradient method, same shape as b
	Array<double,1> Ft;		// residual b - A*x at iteration t, same shape as b
	Array<double,1> Ftplus1;// residual b - A*x at iteration t+1, same shape as b


	// Member-Functions

public:
	// Member Variables
	Array<double,1> x;		// solution
	Array<double,1> b;		// rhs

	// Constructors
	congrad();								// Default Constructor
	congrad(Array<double,1>& b0, Array<double,1>& x0);	// set rhs b = b0 and initial guess x = x0

	// Member-Operators
	congrad& operator =(const congrad& cg);	// Operator =

	// Member-Functions
	void iterate(congradFtn& A, int tmax, double ftol, bool restart = false);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
congrad::congrad()
{
TinyVector <int,1> index(1);
nb = 1;
x.resize(nb);	x.reindexSelf(index);
b.resize(nb);	b.reindexSelf(index);

p.resize(nb);	p.reindexSelf(index);
v.resize(nb);	v.reindexSelf(index);
Ft.resize(nb);	Ft.reindexSelf(index);
Ftplus1.resize(nb);	Ftplus1.reindexSelf(index);
}

congrad::congrad(Array<double,1>& b0, Array<double,1>& x0)
{
TinyVector <int,1> index(1);
nb = b0.rows();

x.resize(nb);	x.reindexSelf(index);
b.resize(nb);	b.reindexSelf(index);

p.resize(nb);	p.reindexSelf(index);
v.resize(nb);	v.reindexSelf(index);
Ft.resize(nb);	Ft.reindexSelf(index);
Ftplus1.resize(nb);	Ftplus1.reindexSelf(index);

x = x0.copy();
b = b0.copy();
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
congrad& congrad::operator = (const congrad& cg)
{
if (this == &cg) return(*this);	    // if: x=x
nb = cg.nb;
x.reference(cg.x);
b.reference(cg.b);
p.reference(cg.p);
v.reference(cg.v);
Ft.reference(cg.Ft);
Ftplus1.reference(cg.Ftplus1);
return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// --- iterate ------------------------------------------------------------------------------------------------------------
void congrad::iterate(congradFtn& A, int tmax, double ftol, bool restart)
{
int i,j,k,t;
double rms,a,g,a1,g1;

if(restart)
{
	cout << "Restarting iterations" << endl;
}
else
{
	//x = 0; 			// initial guess; already done in Constructor
	Ft = b - A(x);		// initial residual
	p = Ft;				// initial search direction
}

a1 = sum(Ft*Ft);		// same as the dot product of the two arrays
rms = sqrt(a1/Ft.size());
cout << "Iteration: 0" << "\t" << "Residual = " << rms << endl;
for(t=1;t<=tmax;t++)
{
	// update solution
	v = A(p);
	a = a1/sum(p*v);
	x += a*p;

	// update residual & check convergence
	Ftplus1 = Ft - a*v;	// same as: b - A(x);
	g1 = sum(Ftplus1*Ftplus1);
	rms = sqrt(g1/Ftplus1.size());
	if(rms <= ftol)
	{
		cout << "Iteration: " << t << "\t" << "Residual = " << rms << endl;
		break;
	}

	// update search direction
	g = g1/a1;
	p = Ftplus1 + g*p;

	// prepare next loop
	Ft = Ftplus1;
	a1 = g1;
	if(t%300 == 0) cout << "Iteration: " << t << "\t" << "Residual = " << rms << endl;
	//cout << "Iteration: " << t << "\t" << "Residual = " << rms << endl;
}
}

//------------------------ End of Class congrad ---------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



//--------- Begin Class congradNR -----------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// finds array x, that minimizes sqrt(|A*x - b|) for given operator A and array b
// using a conjugate gradient iteration method
// uses algorithm from numerical recipes
//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
class congradNR
{
private:
	// Member Variables
	int n;
	Array<double,1> xi;		// search direction

	// Member-Functions
	double dlinmin(congradFtn_f& func, congradFtn& dfunc);
	double f1dim(double t, congradFtn_f& func);
	double df1dim(double t, congradFtn& dfunc);
	void mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, congradFtn_f& func);
	double dbrent(double ax, double bx, double cx, congradFtn_f& func, congradFtn& dfunc, double tol, double& xmin);

public:
	// Member Variables
	Array<double,1> x;		// solution

	// Constructors
	congradNR(Array<double,1>& x0);	// initial guess x = x0

	// Member-Operators
	congradNR& operator =(const congradNR& cg);	// Operator =

	// Member-Functions
	void iterate(congradFtn_f& func, congradFtn& dfunc, int itmax, double ftol);

}; //end of class

//------------------------ Contructors & Operator -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//--------- Constructors --------------------------------------------------------------------------------------------------
// Default Constructor
congradNR::congradNR(Array<double,1>& x0)
{
TinyVector <int,1> index(1);
n = x0.rows();
x.resize(n);	x.reindexSelf(index);
xi.resize(n);	xi.reindexSelf(index);
x = x0.copy();
}

//--------- Operator = ----------------------------------------------------------------------------------------------------
// arrays are just referenced; use A.reference(class.A.copy()) for true copy
congradNR& congradNR::operator = (const congradNR& cg)
{
if (this == &cg) return(*this);	    // if: x=x
n = cg.n;
x.reference(cg.x);
xi.reference(cg.xi);
return(*this);
}

//--------------------- Member Functions ----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a
// function func, using its gradient as calculated by a routine dfunc. The convergence tolerance
// on the function value is input as ftol. Returned quantities are p (the location of the minimum),
// iter (the number of iterations that were performed), and fret (the minimum value of the
// function). The routine dlinmin is called to perform line minimizations.
void congradNR::iterate(congradFtn_f& func, congradFtn& dfunc, int itmax, double ftol)
{
int j,its;
double gg,gam,fp,dgg,fret;
Array<double,1> g(Range(1,n)), h(Range(1,n));
const double eps = 1e-15; // eps is a small number to rectify the special case of converging to exactly zero function value

// Initializations.
fp = func(x);
xi = dfunc(x);

g = -xi;
h = g;
xi = h;

cout << "Iteration: 0" << "\t" << "Residual = " << fp << endl;
for (its=1;its<=itmax;its++) // Loop over iterations.
{
	fret = dlinmin(func,dfunc); // now: x is the minimum location along xi, xi is the displacement form previous x to new x, fret is the minimum value
	if (2.0*fabs(fret - fp) <= ftol*(fabs(fret) + fabs(fp) + eps)) // normal return; minimum achieved, when relative change in residual is less than tolerance
	{
		cout << "Iteration: " << its << "\t" << "Residual = " << fret << endl;
		return;
	}
	xi = dfunc(x);
	gg = sum(g*g);
	//dgg = sum(xi*xi);		// This statement for Fletcher-Reeves.
	dgg = sum((xi+g)*xi);	// This statement for Polak-Ribiere.

	if (gg == 0.0) //Unlikely. If gradient is exactly zero then we are already done.
	{
		cout << "Iteration: " << its << "\t" << "Residual = " << fret << endl;
		return;
	}

	gam = dgg/gg;
	g = -xi;
	h = g + gam*h;
	xi = h;
	fp = fret;
	if(its%300 == 0) cout << "Iteration: " << its << "\t" << "Residual = " << fret << endl;
}
}


// Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
// resets p to where the function func(p) takes on a minimum along the direction xi from p,
// and replaces xi by the actual vector displacement that p was moved. Also returns as fret
// the value of func at the returned location p. This is actually all accomplished by calling the
// routines mnbrak and dbrent.
double congradNR::dlinmin(congradFtn_f& func, congradFtn& dfunc)
{
int j;
double fret;
double TOL = 2.0e-4; 	// Tolerance passed to dbrent.
double xx,xmin,fx,fb,fa,bx,ax;

ax = 0.0; // Initial guess for brackets.
xx = 1.0;
mnbrak(ax,xx,bx,fa,fx,fb,func);
fret = dbrent(ax,xx,bx,func,dfunc,TOL,xmin);
xi *= xmin;
x += xi;
return fret;
}


double congradNR::f1dim(double t, congradFtn_f& func)
{
double f;
Array<double,1> xt(Range(1,n));
xt = x + t*xi;
f = func(xt);
return f;
}


double congradNR::df1dim(double t, congradFtn& dfunc)
{
double df1;
Array<double,1> xt(Range(1,n)), df(Range(1,n));
xt = x + t*xi;
df = dfunc(xt);
df1 = sum(df*xi);
return df1;
}


// Given a function func, and given distinct initial points ax and bx, this routine searches in
// the downhill direction (defined by the function as evaluated at the initial points) and returns
// new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
// values at the three points, fa, fb, and fc.
void congradNR::mnbrak(double& ax, double& bx, double& cx, double& fa, double& fb, double& fc, congradFtn_f& func)
{
double ulim,u,r,q,fu,dum;
const double GOLD = 1.618034;	// Here GOLD is the default ratio by which successive intervals are magnified;
const int GLIMIT = 100;			// GLIMIT is the maximum magnification allowed for a parabolic-fit step.
const double TINY = 1e-20;		// TINY is used to prevent any possible division by zero.

fa = f1dim(ax, func);
fb = f1dim(bx, func);

if (fb > fa) // Switch roles of a and b so that we can go downhill in the direction from a to b.
{
	dum = ax; ax = bx; bx = dum;
	dum = fb; fb = fa; fa = dum;
}
cx = bx + GOLD*(bx - ax); // First guess for c.
fc = f1dim(cx, func);
while (fb > fc) // Keep returning here until we bracket.
{
	r = (bx - ax)*(fb - fc); // Compute u by parabolic extrapolation from a, b, c.
	q = (bx - cx)*(fb - fa);
	u = bx - ((bx - cx)*q - (bx - ax)*r)/(2.0*SIGN(FMAX(fabs(q - r),TINY), q - r));
	ulim = bx + GLIMIT*(cx - bx);
	// We won’t go farther than this. Test various possibilities:
	if ((bx - u)*(u - cx) > 0.0) // Parabolic u is between b and c: try it.
	{
		fu = f1dim(u, func);
		if (fu < fc) // Got a minimum between b and c.
		{
			ax = bx;
			bx = u;
			fa = fb;
			fb = fu;
			return;
		}
		else if (fu > fb) // Got a minimum between between a and u.
		{
			cx = u;
			fc = fu;
			return;
		}
		u = cx +GOLD*(cx - bx); //Parabolic fit was no use. Use default magnification.
		fu = f1dim(u, func);
	}
	else if ((cx - u)*(u - ulim) > 0.0) // Parabolic fit is between c and its allowed limit.
	{
		fu = f1dim(u, func);
		if (fu < fc)
		{
			bx = cx; cx = u; u = cx + GOLD*(cx - bx);
			fb = fc; fc = fu; fu = f1dim(u, func);
		}
	}
	else if ((u - ulim)*(ulim - cx) >= 0.0) // Limit parabolic u to maximum allowed value.
	{
		u = ulim;
		fu = f1dim(u, func);
	}
	else // Reject parabolic u, use default magnification.
	{
		u = cx + GOLD*(cx - bx);
		fu = f1dim(u, func);
	}
	// Eliminate oldest point and continue.
	ax = bx; bx = cx; cx = u;
	fa = fb; fb = fc; fc = fu;
}
}


// Given a function f and its derivative function df, and given a bracketing triplet of abscissas ax,
// bx, cx [such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)],
// this routine isolates the minimum to a fractional precision of about tol using a modification of
// Brent’s method that uses derivatives. The abscissa of the minimum is returned as xmin, and
// the minimum function value is returned as dbrent, the returned function value.
double congradNR::dbrent(double ax, double bx, double cx, congradFtn_f& func, congradFtn& dfunc, double tol, double& xmin)
{
int iter,ok1,ok2; // Will be used as flags for whether proposed steps are acceptable or not.
double a,b,d,d1,d2,du,dv,dw,dx,e = 0.0;
double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

a = (ax < cx ? ax : cx);
b = (ax > cx ? ax : cx);
x = w = v = bx;
fw = fv = fx = f1dim(x, func);
dw = dv = dx = df1dim(x, dfunc); // All our housekeeping chores are doubled by the necessity of moving derivative values around as well as function values.

for (iter=1;iter<=100;iter++)
{
	xm = 0.5*(a + b);
	tol1 = tol*fabs(x) + 1e-10;
	tol2 = 2.0*tol1;
	if (fabs(x - xm) <= (tol2 - 0.5*(b - a)))
	{
		xmin = x;
		return fx;
	}
	if (fabs(e) > tol1)
	{
		d1 = 2.0*(b - a); //Initialize these d’s to an out-of-bracket value.
		d2 = d1;
		if (dw != dx) d1 = (w-x)*dx/(dx-dw); // Secant method with one point.
		if (dv != dx) d2 = (v-x)*dx/(dx-dv); // And the other.
		// Which of these two estimates of d shall we take? We will insist that they be within
		// the bracket, and on the side pointed to by the derivative at x:
		u1 = x + d1;
		u2 = x + d2;
		ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
		ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
		olde = e; // Movement on the step before last.
		e = d;
		if (ok1 || ok2) // Take only an acceptable d, and if both are acceptable, then take the smallest one.
		{
			if (ok1 && ok2) d = (fabs(d1) < fabs(d2) ? d1 : d2);
			else if (ok1) d = d1;
			else d = d2;
			if (fabs(d) <= fabs(0.5*olde))
			{
				u = x + d;
				if (u-a < tol2 || b-u < tol2) d = SIGN(tol1, xm-x);
			}
			else // Bisect, not golden section.
			{
				d = 0.5*(e = (dx >= 0.0 ? a-x : b-x)); //Decide which segment by the sign of the derivative.
			}
		}
		else
		{
			d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
		}
	}
	else
	{
		d = 0.5*(e = (dx >= 0.0 ? a-x : b-x));
	}
	if (fabs(d) >= tol1)
	{
		u = x + d;
		fu = f1dim(u, func);
	}
	else
	{
		u = x + SIGN(tol1,d);
		fu = f1dim(u, func);
		if (fu > fx) // If the minimum step in the downhill direction takes us uphill, then we are done.
		{
			xmin = x;
			return fx;
		}
	}
	du = df1dim(u, dfunc); // Now all the housekeeping, sigh.
	if (fu <= fx)
	{
		if (u >= x) a = x;
		else b = x;
		v = w; fv = fw; dv = dw;
		w = x; fw = fx; dw = dx;
		x = u; fx = fu; dx = du;
	}
	else
	{
		if (u < x) a = u;
		else b = u;
		if (fu <= fw || w == x)
		{
			v = w; fv = fw; dv = dw;
			w = u; fw = fu; dw = du;
		}
		else if (fu < fv || v == x || v == w)
		{
			v = u; fv = fu; dv = du;
		}
	}
}
cout << "Too many iterations in routine dbrent" << endl;
return 0.0; // Never get here.
}

//------------------------ End of Class congradNR -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



#endif //  CONGRAD_CLASS_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

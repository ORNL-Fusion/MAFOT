// all spline and bicubic interpolation routines
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						5.08.10

// Define
//--------
#ifndef SPLINES_INCLUDED
#define SPLINES_INCLUDED

// Include
//--------
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
using namespace blitz;

// Prototypes  
void spline(Array<double,1>& x, Array<double,1>& y, int n, double yp1, double ypn, Array<double,1>& y2);
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, int n, double x, double& y, double& yx);
void splint_2D(Array<double,1>& x1a, Array<double,1>& x2a, Array<double,2>& ya, Array<double,2>& d2ydx1, int n1, int n2,
				double x1, double x2, double& y, double& dydx1, double& dydx2);

void bcuderiv(Array<double,2>& y, double d1, double d2, Array<double,2>& y1, Array<double,2>& y2, Array<double,2>& y12);
void bcucof(Array<double,1>& y, Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12, double d1, double d2, 
			Array<double,2>& c);
int bcuint(Array<double,1>& Ra, Array<double,1>& Za, Array<double,4>& Ca, double dR, double dZ,
			double R, double Z, double& y, double& y1, double& y2);


// ------------------------ spline (mit Blitz-Arrays) -------------------------------------------------------------------
//Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with
//x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative of the interpolating
//function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains
//the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
//ypn are equal to 1e30 or larger, the routine is signaled to set the corresponding boundary
//condition for a natural spline, with zero second derivative on that boundary.
// Aufruf mit z.B.: spline(I,q1,NumberOfPoints,dq1,dqN,qr2);
void spline(Array<double,1>& x, Array<double,1>& y, int n, double yp1, double ypn, Array<double,1>& y2)
{
int i,k;
double p,qn,sig,un;
Array<double,1> u(n);

if (yp1 > 0.99e30) y2(1)=u(1)=0.0;	//The lower boundary condition is set either to be "natural"
else	//or else to have a specified first derivative.
{ 
	y2(1) = -0.5;
	u(1)=(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1);
}

//This is the decomposition loop of the tridiagonal algorithm.
//y2 and u are used for temporary storage of the decomposed factors.
for (i=2;i<=n-1;i++)
{ 
	sig=(x(i)-x(i-1))/(x(i+1)-x(i-1));
	p=sig*y2(i-1)+2.0;
	y2(i)=(sig-1.0)/p;
	u(i)=(y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1));
	u(i)=(6.0*u(i)/(x(i+1)-x(i-1))-sig*u(i-1))/p;
}

if (ypn > 0.99e30) qn=un=0.0;	//The upper boundary condition is set either to be "natural" 
else	//or else to have a specified first derivative.
{ 
	qn=0.5;
	un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)));
}

y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0);

//This is the backsubstitution loop of the tridiagonal algorithm.
for (k=n-1;k>=1;k--)  y2(k)=y2(k)*y2(k+1)+u(k);
}

//--------------------- splint (mit Blitz-Arrays) ---------------------------------------------------------------------------
// Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai's in order),
// and given the array y2a[1..n], which is the output from spline above, and given a value of
// x, this routine returns a cubic-spline interpolated value y(x) and dy/dx(x)
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, int n, double x, double& y, double& yx)
{
int klo,khi,k;
double h,b,a;

klo=1; 
khi=n;

while (khi-klo>1) 
{
	k=(khi+klo) >> 1;
	if (xa(k)>x) khi=k;
	else klo=k;
} //klo and khi now bracket the input value of x.

h=xa(khi)-xa(klo);
if (h==0.0) cout << "Bad xa input to routine splint" << endl; //The xa's must be distinct.
a=(xa(khi)-x)/h;
b=(x-xa(klo))/h; //Cubic spline polynomial is now evaluated.

y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6.0;
yx=(ya(khi)-ya(klo))/h - ((3.0*a*a-1.0)*y2a(klo)-(3.0*b*b-1.0)*y2a(khi))*h/6.0;
}

//----------------- splint-2D (Blitz-Arrays) ------------------------------------------------------------------
//Given x1a, x2a, ya, n1, n2 and d2ydx1 as 1D-spline in x1 direction and
//given a desired interpolating point x1,x2; this routine returns an interpolated function value y
//by bicubic spline interpolation.
void splint_2D(Array<double,1>& x1a, Array<double,1>& x2a, Array<double,2>& ya, Array<double,2>& d2ydx1, int n1, int n2,
				double x1, double x2, double& y, double& dydx1, double& dydx2)
{
int j;
double d1,dn,dx2a;
double d2ydx1dx2;
Array<double,1> ytmp(n2+1),yytmp(n2+1),dyytmp(n2+1),dytmp(n2+1);
Range all = Range::all();
Array<double,1> slice_1,slice_2;	//Slice-array

//Perform n1 evaluations of the row splines (r-direction), using the one-dimensional spline evaluator splint.
for (j=1;j<=n2;j++) 
{
	slice_1.reference(ya(all,j));
	slice_2.reference(d2ydx1(all,j));
	splint(x1a,slice_1,slice_2,n1,x1,yytmp(j),dyytmp(j));
}

//Construct the one-dimensional column spline (z-direction) and evaluate it.
dx2a=x2a(2)-x2a(1);	// equidistant grid!!!!
d1=(yytmp(2)-yytmp(1))/dx2a;	// estimate of first derivative at boundary
dn=(yytmp(n2)-yytmp(n2-1))/dx2a;
spline(x2a,yytmp,n2,d1,dn,ytmp);  
splint(x2a,yytmp,ytmp,n2,x2,y,dydx2); 

// Construct the one-dimensional column spline of the derivative dydx1
d1=(dyytmp(2)-dyytmp(1))/dx2a;	// estimate of first derivative at boundary; does NOT have to be exact
dn=(dyytmp(n2)-dyytmp(n2-1))/dx2a;
spline(x2a,dyytmp,n2,d1,dn,dytmp);  
splint(x2a,dyytmp,dytmp,n2,x2,dydx1,d2ydx1dx2); 
}

//------------------ bcuderiv ---------------------------------------------------------------------------------------------
//Calculates the gradients y1 and y2 and the cross-derivative y12 numerically
//second order accuracy for interior and first order for boundary
//needs function y and grid-distances d1 and d2 (equidistant grid required)
void bcuderiv(Array<double,2>& y, double d1, double d2, Array<double,2>& y1, Array<double,2>& y2, Array<double,2>& y12)
{
int i,j;
const int N1 = y.rows();
const int N2 = y.cols();
const double d1d2 = d1*d2;

// Interior
for(i=2;i<N1;i++)
{
	for(j=2;j<N2;j++)
	{
		y1(i,j) = 0.5*(y(i+1,j)-y(i-1,j))/d1;
		y2(i,j) = 0.5*(y(i,j+1)-y(i,j-1))/d2;
		y12(i,j) = 0.25*(y(i+1,j+1)-y(i-1,j+1)-y(i+1,j-1)+y(i-1,j-1))/d1d2;
	}
}

// Boundaries
for(i=2;i<N1;i++) 
{
	// lower x2-boundary
	y1(i,1) = 0.5*(y(i+1,1)-y(i-1,1))/d1;
	y2(i,1) = (y(i,2)-y(i,1))/d2;
	y12(i,1) = 0.5*(y(i+1,2)-y(i-1,2)-y(i+1,1)+y(i-1,1))/d1d2;

	// upper x2-boundary
	y1(i,N2) = 0.5*(y(i+1,N2)-y(i-1,N2))/d1;
	y2(i,N2) = (y(i,N2)-y(i,N2-1))/d2;
	y12(i,N2) = 0.5*(y(i+1,N2)-y(i-1,N2)-y(i+1,N2-1)+y(i-1,N2-1))/d1d2;
}

for(j=2;j<N2;j++)
{
	// lower x1-boundary
	y1(1,j) = (y(2,j)-y(1,j))/d1;
	y2(1,j) = 0.5*(y(1,j+1)-y(1,j-1))/d2;
	y12(1,j) = 0.5*(y(2,j+1)-y(1,j+1)-y(2,j-1)+y(1,j-1))/d1d2;

	// upper x1-boundary
	y1(N1,j) = (y(N1,j)-y(N1-1,j))/d1;
	y2(N1,j) = 0.5*(y(N1,j+1)-y(N1,j-1))/d2;
	y12(N1,j) = 0.5*(y(N1,j+1)-y(N1-1,j+1)-y(N1,j-1)+y(N1-1,j-1))/d1d2;
}

// 4 corner points of the grid
y1(1,1) = (y(2,1)-y(1,1))/d1;
y2(1,1) = (y(1,2)-y(1,1))/d2;
y12(1,1) = (y(2,2)-y(1,2)-y(2,1)+y(1,1))/d1d2;

y1(1,N2) = (y(2,N2)-y(1,N2))/d1;
y2(1,N2) = (y(1,N2)-y(1,N2-1))/d2;
y12(1,N2) = (y(2,N2)-y(1,N2)-y(2,N2-1)+y(1,N2-1))/d1d2;

y1(N1,1) =(y(N1,1)-y(N1-1,1))/d1;
y2(N1,1) = (y(N1,2)-y(N1,1))/d2;
y12(N1,1) = (y(N1,2)-y(N1-1,2)-y(N1,1)+y(N1-1,1))/d1d2;

y1(N1,N2) = (y(N1,N2)-y(N1-1,N2))/d1;
y2(N1,N2) = (y(N1,N2)-y(N1,N2-1))/d2;
y12(N1,N2) = (y(N1,N2)-y(N1-1,N2)-y(N1,N2-1)+y(N1-1,N2-1))/d1d2;
}

//----------------- bcucof ------------------------------------------------------------------------------------------------
//Given arrays y[1..4], y1[1..4], y2[1..4], and y12[1..4], containing the function, gradients,
//and cross derivative at the four grid points of a rectangular grid cell (numbered counterclockwise
//from the lower left), and given d1 and d2, the length of the grid cell in the 1- and
//2-directions, this routine returns the table c[1..4][1..4] that is used by routine bcuint
//for bicubic interpolation.
void bcucof(Array<double,1>& y, Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12, double d1, double d2, 
			Array<double,2>& c)
{
static int wt[16][16]= {
{ 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
{-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
{2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
{0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
{0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
{0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
{-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
{9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
{-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
{2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
{-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
{4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1} };
int l,k,j,i;
double xx,d1d2,cl[16],x[16];
d1d2=d1*d2;

// Pack a temporary vector x.
for (i=1;i<=4;i++) 
{ 
	x[i-1]=y(i);
	x[i+3]=y1(i)*d1;
	x[i+7]=y2(i)*d2;
	x[i+11]=y12(i)*d1d2;
}
// Matrix multiply by the stored table.
for (i=0;i<=15;i++) 
{ 
	xx=0.0;
	for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
	cl[i]=xx;
}
l=0;

// Unpack the result into the output table.
for (i=1;i<=4;i++) 
	for (j=1;j<=4;j++) c(i,j)=cl[l++];
}

//------------------ bcuint ------------------------------------------------------------------------------------------------
//Bicubic interpolation. Input quantities are Ra and Za containing the grid coordinates,
//Ca contains the bcucof parameters for each grid square, stored with respect to the lower left corner.
//R and Z are the coordinates of the desired point for
//the interpolation. The interpolated function value is returned as y, and the interpolated
//gradient values as y1 and y2. This routine calls bcucof.
int bcuint(Array<double,1>& Ra, Array<double,1>& Za, Array<double,4>& Ca, double dR, double dZ,
			double R, double Z, double& y, double& y1, double& y2)
{
int i,j,k;
int NR = Ra.rows();
int NZ = Za.rows();
double t,u;
Array<double,2> c;
Range all = Range::all();

// Determine grid square where (R,Z) is in
// square includes lower and left boundary, but not upper or right boundary -> in next box included
j = int((R-Ra(1))/dR) + 1;
k = int((Z-Za(1))/dZ) + 1;
if(j == NR) j -= 1;	// exception: add outermost right boundary to square one to the left
if(k == NZ) k -= 1;	// exception: add outermost top boundary to square one down
if(j>NR || j<1 || k>NZ || k<1)	{cout << "bcuint: Point outside of grid" << endl; return -1;}

// Get the c's.
c.reference(Ca(j,k,all,all));

// Interpolate
t=(R-Ra(j))/dR;
u=(Z-Za(k))/dZ;

y = y1 = y2 = 0.0;
for (i=4;i>=1;i--) 
{ 
	y = t*y + ((c(i,4)*u + c(i,3))*u + c(i,2))*u + c(i,1);
	y2 = t*y2 + (3.0*c(i,4)*u + 2.0*c(i,3))*u + c(i,2);
	y1 = u*y1 + (3.0*c(4,i)*t + 2.0*c(3,i))*t + c(2,i);
}

y1 /= dR;
y2 /= dZ;
return 0;
}

#endif //  SPLINES_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

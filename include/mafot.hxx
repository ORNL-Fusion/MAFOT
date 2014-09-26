// Header-File for all non-machine-specific MAFOT subroutines
// uses arrays and multiple-arrays from blitz-Library
// interpolated filament fields routines defined here
// used by all D3D, ITER and NSTX drift programs
// A.Wingen						7.6.11

// Define
//--------
#ifndef MAFOT_INCLUDED
#define MAFOT_INCLUDED

// Include
//--------
#include <andi.hxx>				// includes many usefull tools and defines 
#include <efit_class.hxx>			// includes all EFIT data and interpolation routines
#include <io_class.hxx>				// includes Control file data
#include <particle_class.hxx>		// includes all particle/fieldline parameters and Runge-Kutta Integrator

// --------------- Prototypes ---------------------------------------------------------------------------------------------
bool outofRealBndy(double x, double y, EFIT& EQD);

void get_filament_field(double R, double phi, double Z, Array<double,4>& field, double& bx, double& by, double& bz, EFIT& EQD);
void bcuderiv_square(Array<double,2>& y, int j, int k, double d1, double d2, 
					 Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12);
void bcuint_square(Array<double,1>& Ra, Array<double,1>& Za, double dR, double dZ, Array<double,2>& field,
			double R, double Z, double& y, double& y1, double& y2);

//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//------------ outofRealBndy ----------------------------------------------------------------------------------------------
// Check if (x,y) is out of the torus. It uses the jordan curve theorem with 
// additional detection if (x,y) is part of an edge. Edge is defined as inside
// the torus.
bool outofRealBndy(double x, double y, EFIT& EQD)
{
int wn = 0;

double x1,y1;
double x2 = EQD.wall(1);	//R1
double y2 = EQD.wall(2);	//Z1

double a; 

bool startUeber = (y2 >= y) ? 1 : 0;
for(int i=3; i<(2*EQD.Nwall); i=i+2)
{
	// Continue if two wall point are identical
	if(x2 == EQD.wall(i) && y2 == EQD.wall(i+1)) continue;

	x1 = x2;
	y1 = y2;
	x2 = EQD.wall(i);
	y2 = EQD.wall(i+1);

	if((y1==y2) && (y==y1))
	{
      if(((x1<=x)&&(x<=x2)) || ((x2<=x)&&(x<=x1))) return 0;
    } 
    else if((x1==x2) && (x==x1))
    {
      if(((y1<=y)&&(y<=y2)) || ((y2<=y)&&(y<=y1))) return 0;
    } else {
      a = (x-x1)*(y2-y1)-(y-y1)*(x2-x1);
      if((a <= 1e-15) && (a >= -1e-15)) return 0;	//necessary for use in dtfoot
    }
    bool endUeber = (y2 >= y) ? 1 : 0;
    if(startUeber != endUeber) 
    {
      if((y2 - y)*(x2 - x1) <= (y2 - y1)*(x2 - x))
      {
        if(endUeber) { wn++; }
      } else {
        if(!endUeber) { wn--; }
      }
    }
    startUeber = endUeber;
}
return wn == 0;
} 

//------------------ get_filament_field -----------------------------------------------------------------------------------
// determines sum of magnetic fields of all current filaments at point (R,phi,Z) by interpolation
// field contains total field on EFIT grid
void get_filament_field(double R, double phi, double Z, Array<double,4>& field, double& bx, double& by, double& bz, EFIT& EQD)
{
int xyz;
double dummy;
Range all = Range::all();
Array<double,2> slice;
Array<double,1> y(Range(1,3));

// Get phi slice: phi in deg and rounded to integer
int k = int(phi*rTOd+0.5*sign(phi));	// has to be set between 0 and 359
while(k<0) k += 360;	// move k to a positive value
k = k%360;	// k now between 0 and 359

// Get field in carthesian coordinates
for(xyz=1;xyz<=3;xyz++)
{
	slice.reference(field(xyz,k,all,all));
	bcuint_square(EQD.R,EQD.Z,EQD.dR,EQD.dZ,slice,R,Z,y(xyz),dummy,dummy);
}
bx = y(1);
by = y(2);
bz = y(3);
}

//------------------ bcuderiv_square --------------------------------------------------------------------------------------
//Calculates the gradients y1 and y2 and the cross-derivative y12 numerically
//second order accuracy
//needs function y on the entire grid and the indices (j,k) of the lower left corner of the grid square
//and grid-distances d1 and d2 (equidistant grid required)
//efit_class version modified: use only for grid points at a time
//y(0...NR+1,0...NZ+1), y1(1...4), y2(1...4), y12(1...4)
void bcuderiv_square(Array<double,2>& y, int j, int k, double d1, double d2, 
					 Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12)
{
const double d1d2 = d1*d2;

y1(1) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(1) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(1) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;

j += 1;
y1(2) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(2) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(2) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;

k += 1;
y1(3) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(3) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(3) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;

j -= 1;
y1(4) = 0.5*(y(j+1,k)-y(j-1,k))/d1;
y2(4) = 0.5*(y(j,k+1)-y(j,k-1))/d2;
y12(4) = 0.25*(y(j+1,k+1)-y(j-1,k+1)-y(j+1,k-1)+y(j-1,k-1))/d1d2;
}

//------------------ bcuint_square ----------------------------------------------------------------------------------------
//Bicubic interpolation. Input quantities are Ra and Za containing the grid coordinates,
//field contains the function values on the entire grid
//R and Z are the coordinates of the desired point for
//the interpolation. The interpolated function value is returned as y, and the interpolated
//gradient values as y1 and y2. 
//Modified from efit_class version: This routine calls bcucof and bcuderiv_square.
void bcuint_square(Array<double,1>& Ra, Array<double,1>& Za, double dR, double dZ, Array<double,2>& field,
			double R, double Z, double& y, double& y1, double& y2)
{
int i,j,k;
double t,u;
Array<double,2> c(Range(1,4),Range(1,4));
Range all = Range::all();

// Determine grid square where (R,Z) is in
j = int((R-Ra(1))/dR) + 1;
k = int((Z-Za(1))/dZ) + 1;
if(j>Ra.rows() || j<1 || k>Za.rows() || k<1)	{ofs2 << "Point outside of grid" << endl; EXIT;}

// Get derivatives
Array<double,1> y_sq(Range(1,4)),y1_sq(Range(1,4)),y2_sq(Range(1,4)),y12_sq(Range(1,4));
bcuderiv_square(field,j,k,dR,dZ,y1_sq,y2_sq,y12_sq);

// Get the c’s.
y_sq(1) = field(j,k); y_sq(2) = field(j+1,k); y_sq(3) = field(j+1,k+1); y_sq(4) = field(j,k+1);
bcucof(y_sq,y1_sq,y2_sq,y12_sq,dR,dZ,c);

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
}

#endif // MAFOT_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

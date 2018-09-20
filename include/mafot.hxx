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
#include <restart.hxx>
#ifdef USE_SIESTA
	#include <siesta_class_interpolation.hxx>			// includes the SIESTA interface
	//#include <siesta_class_evaluation.hxx>
#endif
#ifdef USE_XFIELD
	#include <xfield_class.hxx>			// includes the XFIELD interface
	#include <vmec_class.hxx>			// includes the VMEC interface
#endif
#include <particle_class.hxx>		// includes all particle/fieldline parameters and Runge-Kutta Integrator
#include <fakeIsland_class.hxx>
#ifdef m3dc1
	#include <m3dc1_class.hxx>
#endif
#if defined(ITER)
	#include <iter.hxx>
#elif defined(NSTX)
	#include <nstx.hxx>
#elif defined(MAST)
	#include <mast.hxx>
#elif defined(CMOD)
	#include <cmod.hxx>
#else
	#include <d3d.hxx>
#endif

// --------------- Prototypes ---------------------------------------------------------------------------------------------
// int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR); // declared in machine header, defined here
void prepare_common_perturbations(EFIT& EQD, IO& PAR, int mpi_rank, LA_STRING siestafile = "siesta.dat", LA_STRING xpandfile = "xpand.dat",
								  LA_STRING islandfile = "fakeIslands.in", LA_STRING filamentfile = "filament_all.in");
bool outofBndy(double phi, double x, double y, EFIT& EQD);
bool outofBndyInBetween(double phi0, double x0, double y0, double phi1, double x1, double y1, EFIT& EQD);
bool outofRealBndy(double phi, double x, double y, EFIT& EQD);
void point_along_wall(double swall, Array<double,1>& p, EFIT& EQD);
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT);
void get_filament_field(double R, double phi, double Z, Array<double,4>& field, double& bx, double& by, double& bz, EFIT& EQD);
void bcuderiv_square(Array<double,2>& y, int j, int k, double d1, double d2, 
					 Array<double,1>& y1, Array<double,1>& y2, Array<double,1>& y12);
void bcuint_square(Array<double,1>& Ra, Array<double,1>& Za, double dR, double dZ, Array<double,2>& field,
			double R, double Z, double& y, double& y1, double& y2);

// -------------- global Parameters ---------------------------------------------------------------------------------------
// Boundary Box
int simpleBndy = 0;		// 0: real wall boundary 	1: simple boundary box

// log file for errors
ofstream ofs2;

// structures to load in
#ifdef USE_SIESTA
	SIESTA SIES;
#endif
#ifdef USE_XFIELD
	VMEC vmec;
	XFIELD XPND;
#endif
#ifdef m3dc1
	M3DC1 M3D;
#endif

Array<double,4> field;	// for current filaments
fakeIsland FISLD;

//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------- getBfield_general --------------------------------------------------------------------------------------
// part of getBfield that is common to all machines
int getBfield_general(double R, double Z, double phi, double& B_R, double& B_Z, double& B_phi, EFIT& EQD, IO& PAR)
{
int i, chk, chk2;
double psi,dpsidr,dpsidz;
double F;
double X,Y,bx,by,bz;
double B_X,B_Y;
double sinp,cosp;
double coord[3], b_field[3];

coord[0] = R; coord[1] = phi; coord[2] = Z;
B_R = 0; B_phi = 0; B_Z = 0;

sinp = sin(phi);
cosp = cos(phi);

X = R*cosp;
Y = R*sinp;

// Equilibrium field
switch(PAR.response_field)
{
#ifdef USE_XFIELD
case -3:
	if (XPND.useit) XPND.get_B(R, phi, Z, B_R, B_phi, B_Z);
	else vmec.get_B(R, phi, Z, B_R, B_phi, B_Z);
	break;
#endif
#ifdef USE_SIESTA
case -2:
	SIES.get_B(R, phi, Z, B_R, B_phi, B_Z);
	break;
#endif
case -1: case 1: case -10:	// Vacuum equilibrium field from g file
	// get normalized poloidal Flux psi (should be chi in formulas!)
	chk = EQD.get_psi(R,Z,psi,dpsidr,dpsidz);
	if(chk==-1) {ofs2 << "Point is outside of EFIT grid" << endl; B_R=0; B_Z=0; B_phi=1; return -1;}	// integration of this point terminates

	// Equilibrium field
	F = EQD.get_Fpol(psi);
	B_R = dpsidz/R;
	B_phi = F/R;	//B_phi = EQD.Bt0*EQD.R0/R;
	B_Z = -dpsidr/R;
	break;

#ifdef m3dc1
case 0: case 2: 	// M3D-C1: equilibrium field or total field
	for(i=0;i<M3D.nfiles;i++)
	{
		chk = fio_eval_field(M3D.imag[i], coord, b_field);
		if(chk != 0) // field eval failed, probably outside of M3DC1 domain -> fall back to g-file equilibrium
		{
			chk2 = EQD.get_psi(R,Z,psi,dpsidr,dpsidz);
			if(chk2 == -1) {ofs2 << "Point is outside of EFIT grid" << endl; B_R=0; B_Z=0; B_phi=1; return -1;}	// integration of this point terminates

			// Equilibrium field
			F = EQD.get_Fpol(psi);
			B_R = dpsidz/R;
			B_phi = F/R;	//B_phi = EQD.Bt0*EQD.R0/R;
			B_Z = -dpsidr/R;
			break;	// break the for loop
		}
		else
		{
			B_R += b_field[0];
			B_phi += b_field[1];
			B_Z += b_field[2];
		}
	}
	break;
#endif
}

#ifdef m3dc1
// M3D-C1: I-coil perturbation field only, coils are turned off in prep_perturbation
if(PAR.response_field == 1)
{
	for(i=0;i<M3D.nfiles;i++)
	{
		coord[1] = phi + M3D.phase[i];
		chk = fio_eval_field(M3D.imag[i], coord, b_field);
		if(chk != 0) {b_field[0] = 0; b_field[1] = 0; b_field[2] = 0; break;}
		B_R += b_field[0];
		B_phi += b_field[1];
		B_Z += b_field[2];
		//coord[1] = phi;
	}
}
#endif

B_X = 0;	B_Y = 0;

// Field of any current filament
bx = 0;	by = 0;	bz = 0;
if(PAR.useFilament>0) get_filament_field(R,phi,Z,field,bx,by,bz,EQD);

B_X += bx;
B_Y += by;
B_Z += bz;

// Transform B_perturbation = (B_X, B_Y, B_Z) to cylindrical coordinates and add
B_R += B_X*cosp + B_Y*sinp;
B_phi += -B_X*sinp + B_Y*cosp;

if(PAR.response_field == -10)
{
	bx = 0;	bz = 0;
	FISLD.get_B(R,phi,Z,bx,bz,EQD);
	B_R += bx;
	B_Z += bz;
}

return 0;
}

// ------------ prepare_common_perturbations ------------------------------------------------------------------------------
// prepare all perturbations that are not machine specific
void prepare_common_perturbations(EFIT& EQD, IO& PAR, int mpi_rank, LA_STRING siestafile, LA_STRING xpandfile, LA_STRING islandfile, LA_STRING filamentfile)
{
int i,j;
int chk;
LA_STRING line;	// entire line is read by ifstream
ifstream in;

// Prepare SIESTA
#ifdef USE_SIESTA
	if(PAR.response_field == -2)
	{
		if(mpi_rank < 1) cout << "Read SIESTA file" << endl;
		ofs2 << "Read SIESTA file" << endl;
		SIES.read(siestafile);
	}
#endif

// Prepare XFIELD
#ifdef USE_XFIELD
	if(PAR.response_field == -3)
	{
		if (xpandfile == "None")
		{
			XPND.useit = false;
			if(mpi_rank < 1) cout << "No XPAND data. Only inside VMEC last-closed-flux-surface available." << endl;
		}
		else
		{
			if(mpi_rank < 1) cout << "Read XPAND file" << endl;
			ofs2 << "Read XPAND file" << endl;
			XPND.read(xpandfile);
			if(mpi_rank < 1) cout << "NR = " << XPND.NR << "\t Nphi = " << XPND.Np-1 << "\t NZ = " << XPND.NZ << endl;
		}
	}
#endif

// Prepare filaments
if(PAR.useFilament>0)
{
	if(mpi_rank < 1) cout << "Interpolated filament field is used" << endl;
	ofs2 << "Interpolated filament field is used" << endl;
	in.open(filamentfile);
	if(in.fail()==1)
	{
		if(mpi_rank == 1) cout << "Unable to open " << filamentfile << " file. Please run fi_prepare." << endl;
		EXIT;
	}
	else	// Read field on grid from file
	{
		// Set field size
		field.resize(Range(1,3),Range(0,359),Range(0,EQD.NR+1),Range(0,EQD.NZ+1));

		// Skip 3 lines
		in >> line;
		if(mpi_rank < 1) cout << line.mid(3) << endl;
		ofs2 << line.mid(3) << endl;
		in >> line;
		if(mpi_rank < 1) cout << line.mid(3) << endl;
		ofs2 << line.mid(3) << endl;
		in >> line;

		// Read data
		for(int k=0;k<360;k++)
		{
			for(i=0;i<=EQD.NR+1;i++)
			{
				for(int j=0;j<=EQD.NZ+1;j++)
				{
					in >> field(1,k,i,j);
					in >> field(2,k,i,j);
					in >> field(3,k,i,j);
				}
			}
		}
		in.close();
	}
	in.clear();
	if(mpi_rank < 1) cout << endl;
	ofs2 << endl;
}

// Prepare fake Islands
if(PAR.response_field == -10)
{
	if(mpi_rank < 1) cout << "Read Fake Islands file" << endl;
	ofs2 << "Read Fake Islands file" << endl;
	FISLD.read(islandfile);
	FISLD.get_surfaces(EQD);
	if(mpi_rank < 1)
	{
		cout << "Amplitude: " << FISLD.A << endl;
		cout << "m: " << FISLD.m << endl;
		cout << "n: " << FISLD.n << endl;
		cout << "Phase: " << FISLD.delta << endl;
		cout << "Location: " << FISLD.psi0 << endl;
	}
}
}

//----------- outofBndy ---------------------------------------------------------------------------------------------------
// Check if (x,y) is out of the torus. Returns 0 if (x,y)
// is in boundary an 1 if (x,y) is out of boundary.
// simpleBndy = 0; use real wall as boundaries
// simpleBndy = 1: use simple boundary box
bool outofBndy(double phi, double x, double y, EFIT& EQD)
{
if(std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) return false;
switch(simpleBndy)
{
case 0:
	return outofRealBndy(phi,x,y,EQD);
	break;
case 1:
	if(x<bndy[0] || x>bndy[1] || y<bndy[2] || y>bndy[3]) return true;	// bndy defined in machine specific heade file
#ifdef ITER
	if(x>bndy[4] && (y<bndy2[0]*x+bndy2[1] || y>bndy2[2]*x+bndy2[3])) return true;
#endif
	break;
default:
    cout << "simpleBndy switch has a wrong value!" << endl;
}
return false;
}

//------------ outofBndyInBetween -----------------------------------------------------------------------------------------
// Check if line between (x0,y0) at phi0 and (x1,y1) at phi1 is out of the torus.
// toroidal angles phi0 and phi1 are continuous (no modulo), in degrees and right-handed
bool outofBndyInBetween(double phi0, double x0, double y0, double phi1, double x1, double y1, EFIT& EQD)
{
int i;
double x,y,phi,t,step;
int Nangle = EQD.Nwall3D.rows();
double dphi = 360.0/double(Nangle);

double delta = phi1 - phi0;
int N = int(fabs(delta)/dphi);	// number of dphi steps in between
int dir = sign(delta);

for(i=1;i<=N;i++)
{
	step = i*dir*dphi;
	t = step/delta;
	phi = phi0 + step;
	if((t<0) || (t>1)) cout << "Wrong R,Z point in outofBndyInBetween" << endl;
	if((phi1 > phi0) && ((phi < phi0) || (phi > phi1))) cout << "Wrong phi angle in outofBndyInBetween: " << phi0 << "\t" << phi << "\t" << phi1 << endl;
	if((phi1 < phi0) && ((phi > phi0) || (phi < phi1))) cout << "Wrong phi angle in outofBndyInBetween: " << phi1 << "\t" << phi << "\t" << phi0 << endl;
	x = x0 + t*(x1-x0);
	y = y0 + t*(y1-y0);
	if(outofBndy(phi,x,y,EQD)) return true;
}
return false;
}

//------------ outofRealBndy ----------------------------------------------------------------------------------------------
// Check if (x,y) is out of the torus. It uses the jordan curve theorem with 
// additional detection if (x,y) is part of an edge. Edge is defined as inside
// the torus.
// toroidal angle phi is in degrees and right-handed
bool outofRealBndy(double phi, double x, double y, EFIT& EQD)
{
int wn = 0;
int Nwall,p,Nangle;
double pd,dphi,phimod;
Array<double,1> wall;
Range all = Range::all();

if(EQD.use_3Dwall) 	// use 3D wall
{
	Nangle = EQD.Nwall3D.rows();
	dphi = 360.0/double(Nangle);
	pd = phi/dphi;
	p = int(pd + sign(pd)*0.5);		// round pd to nearest integer -> nearest neighbor approximation of 3D wall
	p = p % Nangle;					// p now between -Nangle+1 and +Nangle-1
	if(p < 0) p += Nangle;			// p is now between 0 and +Nangle-1

	phimod = fmod(phi,double(360));
	if(phimod < 0) phimod += 360;
	if((p == 0) && (phimod > 360-dphi)) phimod -= 360;
	if(fabs(p*dphi - phimod) > dphi) cout << "Warning, wrong 3D wall plane: " << p << "\t" << phi << endl;

	Nwall = EQD.Nwall3D(p);
	wall.reference(EQD.wall3D(all,p));
}
else				// use 2D wall from EFIT for every phi
{
	Nwall = EQD.Nwall;
	wall.reference(EQD.wall);
}

double x1,y1;
double x2 = wall(1);	//R1
double y2 = wall(2);	//Z1

double a; 

bool startUeber = (y2 >= y) ? 1 : 0;
for(int i=3; i<(2*Nwall); i=i+2)
{
	// Continue if two wall point are identical
	if(x2 == wall(i) && y2 == wall(i+1)) continue;

	x1 = x2;
	y1 = y2;
	x2 = wall(i);
	y2 = wall(i+1);

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

//------------------ point_along_wall -----------------------------------------------------------------------------------
// locates a point (R,Z) along the wall, based on the length along the wall swall.
void point_along_wall(double swall, Array<double,1>& p, EFIT& EQD)
{
int i, idx, idx_jump;
double x,s1,s0;
Array<double,1> p1(Range(1,2)),p2(Range(1,2)),d(Range(1,2));

swall = fmod(swall,EQD.Swall_max);
if(swall < 0) swall += EQD.Swall_max;

// locate discontinuity
s0 = EQD.Swall(EQD.Nwall);
for(i=1;i<=EQD.Nwall;i++)
{
	s1 = EQD.Swall(i);
	if (fabs(s1 - s0) > 0.5*EQD.Swall_max)	idx_jump = i;	// discontinuity between idx_jump and idx_jump-1
	s0 = s1;
}

// locate intervall that brackets swall
if ((swall > max(EQD.Swall)) || (swall < min(EQD.Swall)))
{
	idx = idx_jump;
	if (fabs(swall - EQD.Swall(idx)) > 0.5*EQD.Swall_max)
	{
		if (swall < EQD.Swall(idx)) swall += EQD.Swall_max;
		else swall -= EQD.Swall_max;
	}
}
else
{
	idx = 1;
	s0 = EQD.Swall(EQD.Nwall);
	for(i=1;i<=EQD.Nwall;i++)
	{
		s1 = EQD.Swall(i);

		if (fabs(s1 - s0) > 0.5*EQD.Swall_max) // skip the jump around Swall = 0 point
		{
			s0 = s1;
			continue;
		}

		if((s1 - swall) * (s0 - swall) <= 0)
		{
			idx = i;
			break;
		}
		s0 = s1;
	}
}

// set bracket points
p1(1) = EQD.wall(2*idx-1);		p1(2) = EQD.wall(2*idx);
if (idx == 1)
{
	p2(1) = EQD.wall(2*EQD.Nwall-1);		p2(2) = EQD.wall(2*EQD.Nwall);
}
else
{
	p2(1) = EQD.wall(2*idx-3);		p2(2) = EQD.wall(2*idx-2);
}

// linear interplation between bracket points
d = p2 - p1;
x = fabs(swall - EQD.Swall(idx))/sqrt(d(1)*d(1)+d(2)*d(2));	// x is dimensionless in [0,1]
p = p1 + x*d;
//cout << swall << "\t" << idx  << "\t" << s0 << "\t" << s1 << "\t" << x << endl;
}



//---------------- start_on_target ----------------------------------------------------------------------------------------
// creates initial conditions on the target plate
// generally t parameterizes the wall as length alomg the wall in m
// starting at HFS midplane, advancing counter-clock-wise, also known as s_wall
// or
// t parametrizes a special target with constant Phi and t=[0 1]
// with points on target are explicitly defined here
// in the contrary to 'set', phi (representing the x coordinate) is varied first here, t second
double start_on_target(int i, int Np, int Nphi, double tmin, double tmax, double phimin, double phimax,
					 EFIT& EQD, IO& PAR, PARTICLE& FLT)
{
int i_p = 0;
int i_phi = 0;
int N = Np*Nphi;
double dp,dphi,t;
Array<double,1> p(Range(1,2));

// Grid stepsizes and t
if(Np == 1) tmax = tmin;
if(Nphi == 1) phimax = phimin;
if(N <= 1) {dp = 0; dphi = 0;}
else
{
	dp = (tmax - tmin)/(N-1);
	dphi = (phimax - phimin)/(N-1);
}
if(dp == 0) i_phi = i - 1;
if(dphi == 0) i_p = i - 1;
if(dp != 0 && dphi != 0)
{
	dp = (tmax - tmin)/double(Np-1);
	dphi = (phimax - phimin)/double(Nphi-1);
	i_phi = (i - 1)%Nphi;
	i_p = int(double(i-1)/double(Nphi));
}
t = tmin + i_p*dp;

// get point (R,Z) that matches parameter t
if(PAR.which_target_plate == 0) point_along_wall(t, p, EQD);
else point_along_target(PAR.which_target_plate, t, p, EQD);

// set initial condition
FLT.R = p(1);
FLT.Z = p(2);
FLT.phi = (phimin + dphi*i_phi)*rTOd;	// phi in deg
FLT.get_psi(p(1),p(2),FLT.psi);

FLT.Lc = 0;
FLT.psimin = 10;

if(FLT.sigma != 0 && PAR.useTprofile == 1) {FLT.set_Energy(); FLT.Lmfp_total = get_Lmfp(FLT.Ekin);}
return t;
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

// Get the c�s.
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

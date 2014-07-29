// Program calculates Poincaré Plot for particle-drift with time dependent perturbations
// Fortran subroutines for Perturbation are used
// A.Wingen						7.6.11

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	Poincaré particle drift data file
//			log-file

// Define
//--------
#if defined(ITER)
	#define program_name "iterplot"
#elif defined(NSTX)
	#define program_name "nstxplot"
#elif defined(MAST)
	#define program_name "mastplot"
#else
	#define program_name "dtplot"
#endif

// Include
//--------
#include <mafot.hxx>
#if defined(ITER)
	#if defined(m3dc1)
		#include <iter_m3dc1.hxx>
	#else
		#include <iter.hxx>
	#endif
#elif defined(NSTX)
	#if defined(m3dc1)
		#include <nstx_m3dc1.hxx>
	#else
		#include <nstx.hxx>
	#endif
#elif defined(MAST)
	#if defined(m3dc1)
		#include <mast_m3dc1.hxx>
	#else
		#include <mast.hxx>
	#endif
#else
	#if defined(m3dc1)
		#include <d3d_m3dc1.hxx>
	#else
		#include <d3d.hxx>
	#endif
#endif

// namespaces
//-----------
using namespace blitz;

// Switches
//----------

// Golbal Parameters
//------------------

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// Variables
int i,n,chk;
EFIT EQD;

// Use system time as seed(=idum) for random numbers
double now=zeit();
long idum=long(now);

// Input file names
LA_STRING basename;
LA_STRING praefix = "";
if(argc==3) praefix = "_" + LA_STRING(argv[2]);
if(argc>=2) basename = LA_STRING(argv[1]);
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// log file
ofs2.open("log_" + LA_STRING(program_name) + praefix + ".dat");

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// additional parameters for IO
PAR.pv[0].name = "Max. Iterations";	PAR.pv[0].wert = PAR.itt;
PAR.pv[1].name = "Points";			PAR.pv[1].wert = PAR.N;
PAR.pv[2].name = "rmin";			PAR.pv[2].wert = PAR.rmin;
PAR.pv[3].name = "rmax";			PAR.pv[3].wert = PAR.rmax;
PAR.pv[4].name = "thmin";			PAR.pv[4].wert = PAR.thmin;
PAR.pv[5].name = "thmax";			PAR.pv[5].wert = PAR.thmax;
PAR.pv[6].name = "phistart";		PAR.pv[6].wert = PAR.phistart;
PAR.pv[7].name = "MapDirection";	PAR.pv[7].wert = PAR.MapDirection;
PAR.pv[8].name = "Ekin";			PAR.pv[8].wert = PAR.Ekin;
PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prep_perturbation(EQD,PAR);

// Prepare particles
PARTICLE FLT(EQD,PAR);

// Output
LA_STRING filenameout = "plot" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(6);
var[0] = "theta[rad]"; var[1] = "r[m]"; var[2] = "phi[deg]"; var[3] = "psi"; var[4] = "R[m]"; var[5] = "Z[m]";
PAR.writeiodata(out,bndy,var);

cout << "Create (0=r-grid, 1=r-random, 2=target, 3=psi-grid, 4=psi-random, 5=R-grid): " << PAR.create_flag << endl;
if(PAR.create_flag==2) cout << "Target: " << PAR.which_target_plate << endl;
cout << endl << "Start Tracer for " << PAR.N << " points ... " << endl;
ofs2 << "Create (0=r-grid, 1=r-random, 2=target, 3=psi-grid, 4=psi-random, 5=R-grid): " << PAR.create_flag << endl;
if(PAR.create_flag==2) ofs2 << "Target: " << PAR.which_target_plate << endl;
ofs2 << endl << "Start Tracer for " << PAR.N << " points ... " << endl;

cout << "Helicity = " << EQD.helicity << endl;

// Get Poincaré section
for(n=1;n<=PAR.N;n++)
{
	// Set initial values in FLT
	switch(PAR.create_flag) 
	{
		case 5:		// grid from R, Z
			FLT.set(n,PAR.N,PAR.Rmin,PAR.Rmax,PAR.Zmin,PAR.Zmax,PAR.NZ,0);
			break;
		case 4:		// random from psi, theta
			FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,2);
			break;
		case 3:		// grid from psi, theta
			FLT.set(n,PAR.N,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,PAR.NZ,2);
			break;
		case 2:		// grid on target from t, phi
			start_on_target(n,PAR.Nt,1,PAR.tmin,PAR.tmax,PAR.phistart,PAR.phistart,EQD,PAR,FLT);
			break;
		case 1:		// random from r, theta
			FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,1);
			break;
		default:	// grid from r, theta
			FLT.set(n,PAR.N,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,1,1);
			break;
	}

	// Integrate
	for(i=1;i<=PAR.itt;i++)
	{
		chk = FLT.mapstep(PAR.MapDirection);
		if(chk<0) {ofs2 << "mapit: wall hit" << endl; break;}	// particle has left system
		
		if(fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) > 1e-10) ofs2 << "wrong toroidal angle: " << fabs(FLT.phi - PAR.MapDirection*i*dpinit*ilt - PAR.phistart) << endl;
		FLT.phi=PAR.MapDirection*i*dpinit*ilt + PAR.phistart;

		// Output
		out << FLT.get_theta() << "\t" << FLT.get_r() << "\t" << FLT.phi << "\t" << FLT.psi << "\t" << FLT.R << "\t" << FLT.Z << endl;
	} // end for i
	ofs2 << "Trax: " << n << "\t" << "Steps: " << i-1 << endl;
} // end for n

double now2=zeit();
cout << endl;
cout << "Program terminates normally, Time: " << now2-now  << " s" << endl;
ofs2 << "Program terminates normally, Time: " << now2-now  << " s" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) m3dc1_unload_file_();
#endif

return 0;
} //end main

//----------------------- End of Main -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// Program searches for periodic fixed points in Poincare-Plot
// include drifts and time dependent perturbations
// derivatives are calculated numerically using a 5-point stencil
// 2-dimensional Newton method is used
// Fortran subroutines are used for perturbations
// A.Wingen						16.06.11

// Input: 1: Parameterfile	2: period of fixed point	3: praefix (optional)
// Output:	fixed points
//			log-file

// Define
//--------
#if defined(ITER)
	#define program_name "iterfix"
#elif defined(NSTX)
	#define program_name "nstxfix"
#elif defined(MAST)
	#define program_name "mastfix"
#elif defined(CMOD)
	#define program_name "cmodfix"
#elif defined(TCABR)
	#define program_name "tcabrfix"
#elif defined(ANYM)
	#define program_name "anyfix"
#else
	#define program_name "dtfix"
#endif

// Include
//--------
#include <mafot.hxx>
#include <unistd.h>

// Prototypes  
//-----------
int gcd(int a, int b);
bool inside(PARTICLE& FLT, IO& PAR);
int newton2D(PARTICLE& FLT, double phistart, int periode);
int mapit_J(PARTICLE& FLT, double J[], int itt, int Map=1, int* mask_out=nullptr);
// --- windowed multishoot additions ---
struct BOX2D { double umin, umax, vmin, vmax; };
bool read_windows_2d(const LA_STRING& path, vector<BOX2D>& W, int expected_n,
                     PARTICLE& FLT, IO& PAR);
static inline bool set_state_from_RZ(PARTICLE& FLT, IO& PAR, double phistart, double R, double Z);
static inline bool find_survivor_in_box(PARTICLE& FLT, IO& PAR, double phistart, int map_dir,
                                       int periode,
                                       const BOX2D& B, int Nu, int Nv,
                                       double& useed, double& vseed);
static inline bool uv_to_RZ(PARTICLE& FLT, IO& PAR,
                            double u, double v,
                            double& R, double& Z);
int multishoot_RZ_AL_LM(
  PARTICLE& FLT, IO& PAR, double phistart,
  int periode, const std::vector<BOX2D>& W,
  int map_dir, int box_dir,
  std::vector<double>& R, std::vector<double>& Z,
  double* final_cost
);
static inline int sigma_index(int i, int n, int box_direction);
bool fill_cycle_arrays(
  PARTICLE& FLT, double phistart, int periode,
  int map_dir,
  double R0, double Z0,
  Array<double,1>& R, Array<double,1>& Z,
  Array<double,1>& psi, Array<double,1>& theta, Array<double,1>& r
);
int classify_OX(PARTICLE& FLT, double phistart, int periode, int map_dir, double R0, double Z0);


// Switches
//----------

// Golbal Parameters
//------------------

// Function Definitions
//---------------------
int main(int argc, char *argv[])
{
// Variables
int i,j,chk;
double x,y;
double dx,dy,xmin,ymin;
EFIT EQD;

// Period of fixed point
int periode;

// Use system time as seed(=idum) for random numbers
double now=zeit();
long idum=long(now);

// defaults
LA_STRING woutfile = "wout.nc";
LA_STRING xpandfile = "xpand.dat";
LA_STRING siestafile = "siesta.dat";
LA_STRING islandfile = "fakeIslands.in";
LA_STRING TprofileFile = "prof_t.dat";
LA_STRING NprofileFile = "prof_n.dat";
double f = 0;  // ratio of impurity to hydrogen ions
double zbar = 2;  // average over impurity ion charge states
double rc = 1;
double mc = 1;
bool use_collision = false;
// --- windowed multishoot options ---
bool use_multishoot = false;
LA_STRING boxfile = "";

// Command line input parsing
int c;
opterr = 0;
while ((c = getopt(argc, argv, "hX:V:S:I:MB:")) != -1)
switch (c)
{
case 'h':
	cout << "usage: dtfix [-h] [-I island] [-S siesta] [-V wout] [-X xpand] file period [tag]" << endl << endl;
	cout << "Calculate the position of periodic points, like x-points and o-points." << endl << endl;
	cout << "positional arguments:" << endl;
	cout << "  file          Contol file (starts with '_')" << endl;
	cout << "  period        periodicity of point" << endl;
	cout << "  tag           optional; arbitrary tag, appended to output-file name" << endl;
	cout << endl << "optional arguments:" << endl;
	cout << "  -h            show this help message and exit" << endl;
	cout << "  -I            filename for mock-up island perturbations; default, see below" << endl;
	cout << "  -S            filename for SIESTA; default, see below" << endl;
	cout << "  -V            filename for VMEC; default, see below" << endl;
	cout << "  -X            filename for XPAND; default, see below" << endl;
	cout << "  -M            Enable windowed multishooting (no scan)" << endl;
	cout << "  -B            filename for boxes (period lines; coords follow create_flag: "
        "0=r/theta, 3/4=psi/theta, 5=R/Z)" << endl;
	cout << endl << "Examples:" << endl;
	cout << "  dtfix _fix.dat 1 blabla" << endl;
	cout << endl << "Infos:" << endl;
	cout << "  To use B-field from M3DC1, set response_field >= 0, and provide file in cwd:" << endl;
	cout << "    m3dc1sup.in    ->  location and scale factor for M3DC1 output C1.h5" << endl;
	cout << "  To use B-field from XPAND, set response_field = -3, and provide files in cwd:" << endl;
	cout << "    xpand.dat      ->  B-field on 3D grid from XPAND; use option -X to specify other filename" << endl;
	cout << "    wout.nc        ->  VMEC output; use option -V to specify other filename" << endl;
	cout << "  To use B-field from SIESTA, set response_field = -2, and provide file in cwd:" << endl;
	cout << "    siesta.dat     ->  B-field on 3D grid; use option -S to specify other filename" << endl;
	cout << "  To use B-field for mock-up islands, set response_field = -10, and provide file in cwd:" << endl;
	cout << "    fakeIslands.in ->  each line gives: Amplitude, pol. mode m, tor. mode n, phase [rad]" << endl;
	cout << "                       use option -I to specify other filename" << endl;
	cout << endl << "Current MAFOT version is: " << MAFOT_VERSION << endl;
	return 0;
case 'I':
	islandfile = optarg;
	break;
case 'S':
	siestafile = optarg;
	break;
case 'V':
	woutfile = optarg;
	break;
case 'X':
	xpandfile = optarg;
	break;
case 'M':
  use_multishoot = true;
  break;
case 'B':
  boxfile = optarg;
  break;
case '?':
	if (optopt == 'c')
		fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	else if (isprint (optopt))
		fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	else
		fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
	exit(0);
default:
	exit(0);
}
// Input file names
LA_STRING basename;
LA_STRING praefix = "";
if(argc==optind+3) praefix = "_" + LA_STRING(argv[optind+2]);
if(argc>=optind+2)
{
	periode = atoi(argv[optind+1]);
	basename = LA_STRING(argv[optind]);
}
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// log file
ofs2.open("log_" + LA_STRING(program_name) + "_" + LA_STRING(periode) + praefix + ".dat");
ofs2.precision(16);

// Read parameter file
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// additional parameters for IO
PAR.pv[6].name = "phistart";		PAR.pv[6].wert = PAR.phistart;
PAR.pv[7].name = "MapDirection";	PAR.pv[7].wert = PAR.MapDirection;
PAR.pv[8].name = "Ekin";			PAR.pv[8].wert = PAR.Ekin;
PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;

// set search grid
const int Nx = PAR.Nr;	// always use Nr = sqrt(N), because NR or NZ or others are given by separate positions in the _fix.dat file
const int Ny = PAR.Nth;	// same here
switch(PAR.create_flag)
{
case 3: case 4:
	PAR.pv[0].name = "psi-grid";		PAR.pv[0].wert = Nx;
	PAR.pv[1].name = "theta-grid";		PAR.pv[1].wert = Ny;
	PAR.pv[2].name = "psimin";			PAR.pv[2].wert = PAR.rmin;
	PAR.pv[3].name = "psimax";			PAR.pv[3].wert = PAR.rmax;
	PAR.pv[4].name = "thmin";			PAR.pv[4].wert = PAR.thmin;
	PAR.pv[5].name = "thmax";			PAR.pv[5].wert = PAR.thmax;
	dx = (PAR.rmax-PAR.rmin)/double(Nx-1);
	dy = (PAR.thmax-PAR.thmin)/double(Ny-1);
	xmin = PAR.rmin;
	ymin = PAR.thmin;
	break;
case 0:
	PAR.pv[0].name = "r-grid";			PAR.pv[0].wert = Nx;
	PAR.pv[1].name = "theta-grid";		PAR.pv[1].wert = Ny;
	PAR.pv[2].name = "rmin";			PAR.pv[2].wert = PAR.rmin;
	PAR.pv[3].name = "rmax";			PAR.pv[3].wert = PAR.rmax;
	PAR.pv[4].name = "thmin";			PAR.pv[4].wert = PAR.thmin;
	PAR.pv[5].name = "thmax";			PAR.pv[5].wert = PAR.thmax;
	dx = (PAR.rmax-PAR.rmin)/double(Nx-1);
	dy = (PAR.thmax-PAR.thmin)/double(Ny-1);
	xmin = PAR.rmin;
	ymin = PAR.thmin;
	break;
case 5:
	PAR.pv[0].name = "R-grid";			PAR.pv[0].wert = Nx;
	PAR.pv[1].name = "Z-grid";			PAR.pv[1].wert = Ny;
	PAR.pv[2].name = "Rmin";			PAR.pv[2].wert = PAR.Rmin;
	PAR.pv[3].name = "Rmax";			PAR.pv[3].wert = PAR.Rmax;
	PAR.pv[4].name = "Zmin";			PAR.pv[4].wert = PAR.Zmin;
	PAR.pv[5].name = "Zmax";			PAR.pv[5].wert = PAR.Zmax;
	dx = (PAR.Rmax-PAR.Rmin)/double(Nx-1);
	dy = (PAR.Zmax-PAR.Zmin)/double(Ny-1);
	xmin = PAR.Rmin;
	ymin = PAR.Zmin;
	break;
default:
	cout << "Unknown create points option" << endl;
	EXIT;
}

// Read EFIT-data
double Raxis = 0, Zaxis = 0;
#ifdef USE_XFIELD
if(PAR.response_field == -3)
{
	cout << "Read VMEC file" << endl;
	ofs2 << "Read VMEC file" << endl;
	vmec.read(woutfile);
	vmec.n0only = true;
	vmec.get_axis(PAR.phistart/rTOd,Raxis,Zaxis);	// need axisymmetric axis only
	vmec.n0only = false;
}
#endif

#ifdef m3dc1
if(PAR.response_field == 0 || PAR.response_field == 2)
{
	M3D.read_m3dc1sup();
	M3D.open_source(PAR.response, PAR.response_field, -1);
	Raxis = M3D.RmAxis;
	Zaxis = M3D.ZmAxis;
	M3D.unload();
}
#endif

EQD.ReadData(EQD.gFile,Raxis,Zaxis);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.gFile << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << "\t" << "gFile: " << EQD.gFile << endl;

// Prepare Perturbation
prepare_common_perturbations(EQD,PAR,0,siestafile,xpandfile,islandfile);
prep_perturbation(EQD,PAR);

// Prepare collisions
COLLISION COL;
if (use_collision) COL.init(TprofileFile, NprofileFile, f, zbar, PAR.Zq, PAR.Mass, rc, mc);
// Prepare particles
PARTICLE FLT(EQD,PAR,COL);

// Use Boundary Box
simpleBndy = 1;

// Output
LA_STRING filenameout = "fix_" + LA_STRING(periode) + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(6);
var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "period";  var[3] = "psi";  var[4] = "theta[rad]";  var[5] = "r[m]";
PAR.writeiodata(out,bndy,var);

// -------------------- MULTISHOOT MODE (AL-LM in RZ, uv-box constraints) --------------------
if(use_multishoot)
{
  if(boxfile == ""){ cout << "Multishoot enabled (-M) but no boxes file (-B).\n"; exit(0); }

  vector<BOX2D> W;
  if(!read_windows_2d(boxfile, W, periode, FLT, PAR)){
    cout << "Could not read " << periode << " boxes from " << boxfile << "\n";
    exit(0);
  }

  FLT.phi = PAR.phistart;
  const int map_dir = 1;
  const int Nu = std::max(2, Nx), Nv = std::max(2, Ny);

  auto solve_dir = [&](int box_dir,
                       vector<double>& Rg, vector<double>& Zg,
                       double& cost)->bool
  {
    vector<double> U(periode), V(periode);

    // 1) seed uv in each assigned box
    for(int i=0;i<periode;i++){
      const int wi = sigma_index(i, periode, box_dir);
      if(!find_survivor_in_box(FLT, PAR, PAR.phistart, map_dir, periode, W[wi], Nu, Nv, U[i], V[i])){
        ofs2 << "MS seed fail: wi=" << wi << " grid=" << Nu << "x" << Nv << "\n";
        return false;
      }
    }

    // 2) uv -> initial RZ guess (with EFIT validity gate)
    Rg.assign(periode,0.0); Zg.assign(periode,0.0);
    for(int i=0;i<periode;i++){
      if(!uv_to_RZ(FLT, PAR, U[i], V[i], Rg[i], Zg[i])) return false;
      double psi0;
      if(FLT.get_psi(Rg[i], Zg[i], psi0, PAR.phistart) < 0) return false;
    }

    // 3) AL-LM solve in RZ with uv-box constraints inside solver
    cost = 0.0;
    return (multishoot_RZ_AL_LM(FLT, PAR, PAR.phistart, periode, W,
                                map_dir, box_dir, Rg, Zg, &cost) >= 0);
  };

  vector<double> Rcw,Zcw,Rccw,Zccw;
  double ccw=1e300, cw=1e300;
  const bool ok_cw  = solve_dir(+1, Rcw,  Zcw,  cw);
  const bool ok_ccw = solve_dir(-1, Rccw, Zccw, ccw);

  if(!ok_cw && !ok_ccw){ cout << "Multishoot failed (cw/ccw).\n"; exit(0); }

  const int box_dir = (ok_cw && (!ok_ccw || cw<=ccw)) ? +1 : -1;
  const vector<double>& Rg = (box_dir>0) ? Rcw  : Rccw;
  const vector<double>& Zg = (box_dir>0) ? Zcw  : Zccw;

  Array<double,1> Rrep(periode), Zrep(periode), psirep(periode), thetarep(periode), rrrep(periode);
  fill_cycle_arrays(FLT, PAR.phistart, periode, map_dir, Rg[0], Zg[0], Rrep, Zrep, psirep, thetarep, rrrep);

  double max_node_mismatch = 0.0;
  for(int k=0;k<periode;k++){
    double dR = Rg[k] - Rrep(k);
    double dZ = Zg[k] - Zrep(k);
    max_node_mismatch = std::max(max_node_mismatch, std::sqrt(dR*dR + dZ*dZ));
  }
  ofs2 << "Max mismatch between multishoot nodes and replayed orbit = "
       << max_node_mismatch << "\n";
  ofs2 << "cw ok=" << ok_cw  << " cost=" << cw  << "\n";
  ofs2 << "ccw ok="<< ok_ccw << " cost=" << ccw << "\n";
  ofs2 << "chosen box_dir=" << box_dir << "\n";



  // node 0 seed
  const double R0 = Rg[0], Z0 = Zg[0];
  double psi0;
  if(FLT.get_psi(R0, Z0, psi0, PAR.phistart) < 0){
    cout << "Chosen seed outside EFIT.\n"; exit(0);
  }
  FLT.psi = psi0;

  Array<double,1> R(periode), Z(periode), psi(periode), theta(periode), rr(periode);
  if(!fill_cycle_arrays(FLT, PAR.phistart, periode, map_dir, R0, Z0, R, Z, psi, theta, rr)){
    cout << "Replay failed after multishoot.\n"; exit(0);
  }

  const int ox = classify_OX(FLT, PAR.phistart, periode, map_dir, R0, Z0);
  for(int k=0;k<periode;k++)
    out << R(k) << "\t" << Z(k) << "\t" << ox*periode
        << "\t" << psi(k) << "\t" << theta(k) << "\t" << rr(k) << "\n";

  ofs2 << "Program terminated normally (multishoot mode)\n";
  cout << "Program terminated normally (multishoot mode)\n";

  #ifdef m3dc1
  if(PAR.response_field >= 0) M3D.unload();
  #endif
  return 0;
}
// -------------------- END MULTISHOOT MODE --------------------

// get rational surface
double qrat = 1;
int nrat = 1;
int mrat = periode;
int mngcd;
if(periode > 1)
{
	qrat = EQD.get_q(0.5*(PAR.rmin + PAR.rmax));
	nrat = int(mrat/qrat + 0.5);
	ofs2 << "Searching at rational surface q = " << mrat << "/" << nrat << endl;
	cout << "Searching at rational surface q = " << mrat << "/" << nrat << endl;
	mngcd = gcd(mrat,nrat);
	if(mngcd > 1)
	{
		periode /= mngcd;
		ofs2 << "Periodic points are heteroclinic. Use new period = " << periode << endl;
		cout << "Periodic points are heteroclinic. Use new period = " << periode << endl;
	}
}

Array<double,1> Rall(2*mrat), Zall(2*mrat);
int idx_all = 0;
Array<double,1> R(periode),Z(periode),psi(periode),theta(periode),r(periode);
double dist;
ofs2 << Nx << " rows, done:" << endl;
for(i=0;i<Nx;i++)	//r or R
{
	x = xmin + i*dx;
	for(j=0;j<Ny;j++)	//theta or Z
	{
		y = ymin + j*dy;
		switch(PAR.create_flag)
		{
		case 0:		// r-theta grid
			FLT.convertRZ(y,x);
			break;
		case 3: 	// psi-theta grid
			FLT.psi = x;
			FLT.getRZ(x, y, FLT.R, FLT.Z);
			break;
		case 4:		// psi-theta random
			FLT.create(idum,PAR.rmin,PAR.rmax,PAR.thmin,PAR.thmax,2);
			break;
		case 5:		// R,Z
			FLT.R = x;
			FLT.Z = y;
			break;
		}

		chk = newton2D(FLT,PAR.phistart,periode);
		if(chk == -1) continue;

		if(inside(FLT,PAR))
		{
			if(periode == 1)
			{
				out << FLT.R << "\t" << FLT.Z << "\t" << -periode << "\t" << FLT.psi << "\t" << FLT.get_theta() << "\t" << FLT.get_r() << endl;
				cout << "Program terminated normally" << endl;
				ofs2 << "Program terminated normally" << endl;
				#ifdef m3dc1
				if(PAR.response_field >= 0) M3D.unload();
				#endif
				return 0; 
			}

			// check if point is already in list
			chk = 0;
			for(int k=0;k<PAR.N;k++)
			{
				if(Rall(k) == 0) break;
				dist = sqrt((Rall(k)-FLT.R)*(Rall(k)-FLT.R) + (Zall(k)-FLT.Z)*(Zall(k)-FLT.Z));
				if(dist < 1e-5) {chk = -1; break;}	// same point already in list
			}
			if(chk < 0) continue;

			R(0) = FLT.R; Z(0) = FLT.Z; psi(0) = FLT.psi; r(0) = FLT.get_r(); theta(0) = FLT.get_theta();

			// map out all intermediate points and check if periodic point actually has the period periode or a smaller one
			for(int k=1;k<periode;k++)
			{
				chk = FLT.mapstep();
				if(chk < 0) {break;}	// particle has left system
				dist = sqrt((R(0)-FLT.R)*(R(0)-FLT.R) + (Z(0)-FLT.Z)*(Z(0)-FLT.Z));
				if(dist < 1e-5) {chk = -1; break;}	// same point with k < periode, so not a periodic point we are looking for
				R(k) = FLT.R; Z(k) = FLT.Z; psi(k) = FLT.psi; r(k) = FLT.get_r(); theta(k) = FLT.get_theta();
			}
			if(chk < 0) continue;

			// Check if periodic point is an O-point or X-point
			FLT.R = R(0) + 1e-5; FLT.Z = Z(0) + 1e-5;	// perturb the periodic point
			chk = FLT.mapit(10*periode);
			dist = sqrt((R(0)-FLT.R)*(R(0)-FLT.R) + (Z(0)-FLT.Z)*(Z(0)-FLT.Z));
			if(dist > 1e-2) chk = -1;		// point has moved away by more than a cm
			if(chk >= 0) chk = 1;

			// chk == -1: X-point; chk == 1: O-point
			for(int k=0;k<periode;k++)
			{
				Rall(idx_all) = R(k); Zall(idx_all) = Z(k);
				idx_all += 1;
				out << R(k) << "\t" << Z(k) << "\t" << chk*periode << "\t" << psi(k) << "\t" << theta(k) << "\t" << r(k) << endl;
			}
		}

		if(idx_all >= 2*mrat)
		{
			ofs2 << "All possible periodic points have been found." << endl;
			cout << "All possible periodic points have been found." << endl;
			ofs2 << "Program terminated normally" << endl;
			cout << "Program terminated normally" << endl;

			#ifdef m3dc1
			if(PAR.response_field >= 0) M3D.unload();
			#endif

			return 0;
		}
	}
	ofs2 << i+1 << "\t" << flush;
}
ofs2 << endl;
ofs2 << "End of search grid. Not all periodic points have been found." << endl;
cout << "End of search grid. Not all periodic points have been found." << endl;
ofs2 << "Program terminated normally" << endl;
cout << "Program terminated normally" << endl;

#ifdef m3dc1
if(PAR.response_field >= 0) M3D.unload();
#endif

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------- gcd -------------------------------------------------------
// find greatest common divisor
int gcd(int a, int b)
{
   if (b == 0) return a;
   return gcd(b, a % b);

   /*
   if (a == 0 || b == 0)
      return 0;
   else if (a == b)
      return a;
   else if (a > b)
      return gcd(a-b, b);
   else return gcd(a, b-a);
   */
}

//----------- inside ----------------------------------------------------
// checks if current FLT point is inside the search grid
bool inside(PARTICLE& FLT, IO& PAR)
{
double r,th;
switch(PAR.create_flag)
{
case 0: case 1:
	r = FLT.get_r();
	th = FLT.get_theta();
	if(r >= PAR.rmin && r <= PAR.rmax && th >= PAR.thmin && th <= PAR.thmax) return true;
	else return false;
case 3: case 4:
	th = FLT.get_theta();
	if(FLT.psi >= PAR.rmin && FLT.psi <= PAR.rmax && th >= PAR.thmin && th <= PAR.thmax) return true;
	else return false;
case 5:
	if(FLT.R >= PAR.Rmin && FLT.R <= PAR.Rmax && FLT.Z >= PAR.Zmin && FLT.Z <= PAR.Zmax) return true;
	else return false;
}
return false;
}

//----------- newton2dim ------------------------------------------------
// all variables with 'r' represent R and all with 't' represent Z!!!
// J[0] = dR(i+1)/dR(i),  J[1] = dR(i+1)/dZ(i),  J[2] = dZ(i+1)/dR(i),  J[3] = dZ(i+1)/dZ(i)
int newton2D(PARTICLE& FLT, double phistart, int periode)	//0: ok		-1: Error
{
double fr,ft,dr,dt,det,length;
int i,chk;

const int imax = 100;
const double delta = 1e-12;

// Vectors
double J[4];

// Search
double R = FLT.R;
double Z = FLT.Z;

for(i=0;i<=imax;i++)
{
	FLT.phi = phistart;
	chk = mapit_J(FLT,J,periode);
	if(chk<0){ofs2 << "No convergence " << chk << endl; return -1;}
	
	fr = FLT.R - R;
	ft = FLT.Z - Z;

	det = (J[0]-1)*(J[3]-1) - J[1]*J[2];
	dr = ((J[3]-1)*fr-J[1]*ft)/det;
	dt = ((J[0]-1)*ft-J[2]*fr)/det;

	length = sqrt(dr*dr + dt*dt);
	if(length<delta) 
	{
		return 0;	// convergence
	} 

	R -= dr;
	Z -= dt;
	FLT.R = R;
	FLT.Z = Z;
}

ofs2 << "No convergence " <<  R << "\t" << Z << "\t" << dr << "\t" << dt << "\t" << length << endl;
return -1;
}

//--------- mapit_J -----------------------------------------------------------------
// return 0 if all succeeded; otherwise return negative bitmask (-(mask))
int mapit_J(PARTICLE& FLT, double J[], int itt, int Map, int* mask_out)
{
  int chk;
  int mask = 0;

  double Ralt = FLT.R, Zalt = FLT.Z, phialt = FLT.phi;
  const double dr = 1e-5; 

  double R_stencil[4], Z_stencil[4];

  auto try_map = [&](double R0, double Z0, double& Rout, double& Zout, int bit)->bool{
    FLT.R = R0; FLT.Z = Z0; FLT.phi = phialt;
    int c = FLT.mapit(itt, Map);
    if(c < 0){ mask |= (1<<bit); return false; }
    Rout = FLT.R; Zout = FLT.Z;
    return true;
  };

  // R+dr
  bool ok0 = try_map(Ralt+dr, Zalt, R_stencil[0], Z_stencil[0], 0);
  // R-dr
  bool ok1 = try_map(Ralt-dr, Zalt, R_stencil[1], Z_stencil[1], 1);
  // Z+dr
  bool ok2 = try_map(Ralt, Zalt+dr, R_stencil[2], Z_stencil[2], 2);
  // Z-dr
  bool ok3 = try_map(Ralt, Zalt-dr, R_stencil[3], Z_stencil[3], 3);

  // center
  FLT.R = Ralt; FLT.Z = Zalt; FLT.phi = phialt;
  chk = FLT.mapit(itt, Map);
  if(chk < 0) mask |= (1<<4);

  if(mask_out) *mask_out = mask;

  if(mask != 0){
    // restore
    FLT.R = Ralt; FLT.Z = Zalt; FLT.phi = phialt;
    return -mask;
  }

  // derivatives (only if all four stencils succeeded)
  J[0] = 0.5*(R_stencil[0]-R_stencil[1])/dr;
  J[1] = 0.5*(R_stencil[2]-R_stencil[3])/dr;
  J[2] = 0.5*(Z_stencil[0]-Z_stencil[1])/dr;
  J[3] = 0.5*(Z_stencil[2]-Z_stencil[3])/dr;

  return 0;
}

// -------------------- WINDOWED MULTI-SHOOT (AL-LM in RZ, uv-box constraints) --------------------
// box_direction: +1 => W0,W1,...,Wn-1 ; -1 => W0,Wn-1,...,W1

// set FLT from R,Z and sync psi (EFIT validity gate)
// (ONLY ONE definition!)
// set FLT from R,Z and sync psi (EFIT validity gate)
static inline double unwrap_near(double a, double a0){
  while(a - a0 >  M_PI) a -= 2.0*M_PI;
  while(a - a0 < -M_PI) a += 2.0*M_PI;
  return a;
}

static inline void unwrap_uv_to_box_center(const IO& PAR, const BOX2D& B,
                                           double& u, double& v,
                                           double& umin, double& umax,
                                           double& vmin, double& vmax)
{
  umin = B.umin; umax = B.umax;
  vmin = B.vmin; vmax = B.vmax;

  const bool v_is_theta = (PAR.create_flag == 0 || PAR.create_flag == 3 || PAR.create_flag == 4);
  if(!v_is_theta) return;

  // Use the *current* v as the reference chart.
  const double vc = v;

  // Keep v as-is (unwrap_near(v, v) is a no-op), but it makes intent explicit:
  v    = unwrap_near(v,    vc);
  vmin = unwrap_near(vmin, vc);
  vmax = unwrap_near(vmax, vc);

  if(vmin > vmax) std::swap(vmin, vmax);
}

static inline bool set_state_from_RZ(PARTICLE& FLT, IO&, double phistart, double R, double Z)
{
  FLT.R   = R;
  FLT.Z   = Z;
  FLT.phi = phistart;

  double psi;
  if (FLT.get_psi(R, Z, psi, phistart) < 0) return false;
  FLT.psi = psi;
  return true;
}

static inline int sigma_index(int i, int n, int box_direction)
{
  i %= n; if (i < 0) i += n;
  if (box_direction >= 0) return i;
  if (i == 0) return 0;
  return n - i; // 1..n-1 -> n-1..1
}

static inline bool box_center_halfwidth(const BOX2D& B,
                                        double& uc, double& vc,
                                        double& hu, double& hv)
{
  uc = 0.5*(B.umin + B.umax);
  vc = 0.5*(B.vmin + B.vmax);
  hu = 0.5*(B.umax - B.umin);
  hv = 0.5*(B.vmax - B.vmin);
  return (hu > 0.0 && hv > 0.0);
}

// u,v -> R,Z (matches dtfix scan conventions)
static inline bool uv_to_RZ(PARTICLE& FLT, IO& PAR,
                            double u, double v,
                            double& R, double& Z)
{
  if (PAR.create_flag == 5) { R = u; Z = v; return true; }
  if (PAR.create_flag == 0) { FLT.convertRZ(v, u); R = FLT.R; Z = FLT.Z; return true; }          // u=r, v=theta
  if (PAR.create_flag == 3 || PAR.create_flag == 4) { FLT.getRZ(u, v, R, Z); return true; }     // u=psi, v=theta
  return false;
}

// read u,v (the coordinates boxes are defined in) from the *current* FLT state
static inline bool state_to_uv(PARTICLE& FLT, IO& PAR, double& u, double& v)
{
  if (PAR.create_flag == 5) { u = FLT.R;        v = FLT.Z;         return true; }
  if (PAR.create_flag == 0) { u = FLT.get_r(); v = FLT.get_theta(); return true; }
  if (PAR.create_flag == 3 || PAR.create_flag == 4) { u = FLT.psi; v = FLT.get_theta(); return true; }
  return false;
}

// inequality constraints g<=0 for node i in its assigned box:
// g0=umin-u, g1=u-umax, g2=vmin-v, g3=v-vmax
static inline void box_g(const BOX2D& B, double u, double v, double g[4])
{
  g[0] = B.umin - u;
  g[1] = u - B.umax;
  g[2] = B.vmin - v;
  g[3] = v - B.vmax;
}

static inline bool uv_grad_fd(PARTICLE& FLT, IO& PAR, double phistart,
                              double R, double Z, double epsR, double epsZ,
                              double& du_dR, double& du_dZ,
                              double& dv_dR, double& dv_dZ)
{
  double u0, v0, up, vp, um, vm;

  if (!set_state_from_RZ(FLT, PAR, phistart, R, Z)) return false;
  if (!state_to_uv(FLT, PAR, u0, v0)) return false;

  // Gate: only unwrap if v is an angle (theta)
  const bool v_is_theta = (PAR.create_flag == 0 || PAR.create_flag == 3 || PAR.create_flag == 4);

  // d/dR
  if (!set_state_from_RZ(FLT, PAR, phistart, R + epsR, Z)) return false;
  if (!state_to_uv(FLT, PAR, up, vp)) return false;
  if (!set_state_from_RZ(FLT, PAR, phistart, R - epsR, Z)) return false;
  if (!state_to_uv(FLT, PAR, um, vm)) return false;

  if (v_is_theta) {
    vp = unwrap_near(vp, v0);
    vm = unwrap_near(vm, v0);
  }

  du_dR = (up - um) / (2.0 * epsR);
  dv_dR = (vp - vm) / (2.0 * epsR);

  // d/dZ
  if (!set_state_from_RZ(FLT, PAR, phistart, R, Z + epsZ)) return false;
  if (!state_to_uv(FLT, PAR, up, vp)) return false;
  if (!set_state_from_RZ(FLT, PAR, phistart, R, Z - epsZ)) return false;
  if (!state_to_uv(FLT, PAR, um, vm)) return false;

  if (v_is_theta) {
    vp = unwrap_near(vp, v0);
    vm = unwrap_near(vm, v0);
  }

  du_dZ = (up - um) / (2.0 * epsZ);
  dv_dZ = (vp - vm) / (2.0 * epsZ);

  // restore base
  set_state_from_RZ(FLT, PAR, phistart, R, Z);
  return true;
}

// Dense solver (Gaussian elim w/ partial pivot)
static inline bool solve_dense(std::vector<double>& A, std::vector<double>& b, std::vector<double>& x, int n){
  x.assign(n, 0.0);
  for(int k=0;k<n;k++){
    int piv = k;
    double amax = std::fabs(A[k*n + k]);
    for(int i=k+1;i<n;i++){
      double v = std::fabs(A[i*n + k]);
      if(v > amax){ amax=v; piv=i; }
    }
    if(amax < 1e-30) return false;
    if(piv != k){
      for(int j=k;j<n;j++) std::swap(A[k*n + j], A[piv*n + j]);
      std::swap(b[k], b[piv]);
    }
    const double Akk = A[k*n + k];
    for(int i=k+1;i<n;i++){
      const double f = A[i*n + k] / Akk;
      A[i*n + k] = 0.0;
      for(int j=k+1;j<n;j++) A[i*n + j] -= f * A[k*n + j];
      b[i] -= f * b[k];
    }
  }
  for(int i=n-1;i>=0;i--){
    double s = b[i];
    for(int j=i+1;j<n;j++) s -= A[i*n + j] * x[j];
    x[i] = s / A[i*n + i];
  }
  return true;
}

// -------------------- AL-LM multishoot in RZ --------------------
// Variables: x = (R0,Z0,...,R_{n-1},Z_{n-1})
// Residuals: r_i = x_{i+1} - F(x_i) in RZ (2 per i), F is 1-step map (map_dir)
// Constraints: uv(x_i) must stay inside assigned box W[sigma_index(i,n,box_dir)] via augmented Lagrangian.
int multishoot_RZ_AL_LM(
  PARTICLE& FLT, IO& PAR, double phistart,
  int periode, const std::vector<BOX2D>& W,
  int map_dir, int box_dir,
  std::vector<double>& R, std::vector<double>& Z,
  double* final_cost
){
  const int n   = periode;
  const int dim = 2*n;         // x dimension
  const int m   = 2*n;         // residual dimension

  if ((int)R.size()!=n) R.assign(n,0.0);
  if ((int)Z.size()!=n) Z.assign(n,0.0);
  if (final_cost) *final_cost = 0.0;

  // LM controls (inner)
  const int    lm_imax   = 50;
  double       mu        = 1e-3;
  const double mu_up     = 10.0, mu_dn = 0.2;
  const double tol_step  = 1e-12;
  const double tol_rnorm = 1e-10;
  const double g_band = 1e-8; 


  // AL controls (outer)
  const int    al_imax   = 8;
  double       rho       = 1e3;         // penalty weight
  const double rho_up    = 10.0;
  const double tol_gmax  = 0.5*g_band;

  // guardrail: allow some slack before hard-failing (scaled by halfwidth)
  const double guard_mult = 2.5;

  std::vector<double> lambda(4*n, 0.0); // multipliers for g<=0

  // work buffers
  std::vector<double> rvec(m,0.0), rtry(m,0.0);
  std::vector<double> JTJ(dim*dim,0.0), JTr(dim,0.0);
  std::vector<double> A(dim*dim,0.0), b(dim,0.0), dx(dim,0.0);
  std::vector<double> Rtry(n,0.0), Ztry(n,0.0);

  auto idxR = [&](int i){ return 2*i; };
  auto idxZ = [&](int i){ return 2*i+1; };

  // Evaluate residual r(x) and AL terms. Also returns max constraint violation.
  auto eval_all = [&](const std::vector<double>& Rin,
                      const std::vector<double>& Zin,
                      std::vector<double>& rout,
                      double& gmax,
                      double& phi) -> bool
  {
    if ((int)rout.size()!=m) rout.assign(m,0.0);
    std::fill(rout.begin(), rout.end(), 0.0);
    gmax = 0.0;
    phi  = 0.0;

    const double Rsave=FLT.R, Zsave=FLT.Z, psisave=FLT.psi, phisave=FLT.phi;

    double cost_r = 0.0;

    for (int i=0;i<n;i++){
      const int ip = (i+1)%n;

      // set state at x_i
      if (!set_state_from_RZ(FLT, PAR, phistart, Rin[i], Zin[i])){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }

      // one-step map + DF via mapit_J
      double DF[4];
      const double Ri = Rin[i], Zi = Zin[i];
      if (mapit_J(FLT, DF, 1, map_dir) < 0){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }
      const double FR = FLT.R;
      const double FZ = FLT.Z;

      // ensure x_{i+1} is EFIT-valid too
      if (!set_state_from_RZ(FLT, PAR, phistart, Rin[ip], Zin[ip])){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }

      const double riR = Rin[ip] - FR;
      const double riZ = Zin[ip] - FZ;
      rout[2*i+0] = riR;
      rout[2*i+1] = riZ;
      cost_r += riR*riR + riZ*riZ;

      // constraints at node i in its assigned box (in uv coords)
      const int wi = sigma_index(i, n, box_dir);
      const BOX2D& B = W[wi];

      // compute uv at node i from *its* RZ
      if (!set_state_from_RZ(FLT, PAR, phistart, Ri, Zi)){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }
      double u,v;
	  if (!state_to_uv(FLT, PAR, u, v)) { 
		FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
	  }

	  // unwrap (theta only) and get a consistent interval
	  double umin, umax, vmin, vmax;
	  unwrap_uv_to_box_center(PAR, B, u, v, umin, umax, vmin, vmax);

	  // guardrail (use unwrapped center/halfwidth)
	  {
	  const double uc = 0.5*(umin+umax), vc = 0.5*(vmin+vmax);
	  const double hu = 0.5*(umax-umin), hv = 0.5*(vmax-vmin);
	  if (hu > 0.0 && hv > 0.0){
		  const double xu = std::fabs((u-uc)/hu);
		  const double xv = std::fabs((v-vc)/hv);
		  if (xu > guard_mult || xv > guard_mult){
		  FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
		  return false;
		  }
	    }
	  }

	  // constraints g<=0 using the unwrapped interval
	  double g[4];
	  g[0] = umin - u;
	  g[1] = u - umax;
	  g[2] = vmin - v;
	  g[3] = v - vmax;

      for (int k=0;k<4;k++){
        // banded hinge: gp is zero for g<=-g_band, linear ramp to g for g>=0
		double gp = 0.0;
		if (g[k] >= 0.0) gp = g[k];
		else if (g[k] > -g_band) gp = (g[k] + g_band);   // in (0, g_band)

		gmax = std::max(gmax, gp);
		phi += lambda[4*i+k]*gp + 0.5*rho*gp*gp;
      }
    }

    phi += 0.5*cost_r;
    FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
    return true;
  };

  // Build (JTJ, JTr) for residual r using mapit_J-derived DF blocks
  auto build_residual_normal_eq = [&](const std::vector<double>& Rin,
                                      const std::vector<double>& Zin,
                                      const std::vector<double>& r0) -> bool
  {
    std::fill(JTJ.begin(), JTJ.end(), 0.0);
    std::fill(JTr.begin(), JTr.end(), 0.0);

    const double Rsave=FLT.R, Zsave=FLT.Z, psisave=FLT.psi, phisave=FLT.phi;

    for (int i=0;i<n;i++){
      const int ip = (i+1)%n;

      if (!set_state_from_RZ(FLT, PAR, phistart, Rin[i], Zin[i])){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }

      double DF[4];
      if (mapit_J(FLT, DF, 1, map_dir) < 0){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }

      const double riR = r0[2*i+0];
      const double riZ = r0[2*i+1];

      // JTr
      JTr[idxR(ip)] += riR;
      JTr[idxZ(ip)] += riZ;
      JTr[idxR(i)]  += -(DF[0]*riR + DF[2]*riZ);
      JTr[idxZ(i)]  += -(DF[1]*riR + DF[3]*riZ);

      // JTJ
      const double a00 = DF[0], a01 = DF[1], a10 = DF[2], a11 = DF[3];

      JTJ[idxR(ip)*dim + idxR(ip)] += 1.0;
      JTJ[idxZ(ip)*dim + idxZ(ip)] += 1.0;

      JTJ[idxR(i)*dim + idxR(i)] += a00*a00 + a10*a10;
      JTJ[idxR(i)*dim + idxZ(i)] += a00*a01 + a10*a11;
      JTJ[idxZ(i)*dim + idxR(i)] += a01*a00 + a11*a10;
      JTJ[idxZ(i)*dim + idxZ(i)] += a01*a01 + a11*a11;

      JTJ[idxR(i)*dim + idxR(ip)] += -a00;
      JTJ[idxR(i)*dim + idxZ(ip)] += -a10;
      JTJ[idxZ(i)*dim + idxR(ip)] += -a01;
      JTJ[idxZ(i)*dim + idxZ(ip)] += -a11;

      JTJ[idxR(ip)*dim + idxR(i)] += -a00;
      JTJ[idxR(ip)*dim + idxZ(i)] += -a01;
      JTJ[idxZ(ip)*dim + idxR(i)] += -a10;
      JTJ[idxZ(ip)*dim + idxZ(i)] += -a11;
    }

    FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
    return true;
  };

  // Add AL constraint terms (active constraints only)
  auto add_constraint_terms = [&](const std::vector<double>& Rin,
                                  const std::vector<double>& Zin) -> bool
  {
    const double epsR = 1e-6, epsZ = 1e-6;
    const double Rsave=FLT.R, Zsave=FLT.Z, psisave=FLT.psi, phisave=FLT.phi;

    for (int i=0;i<n;i++){
      const int wi = sigma_index(i, n, box_dir);
      const BOX2D& B = W[wi];

      if (!set_state_from_RZ(FLT, PAR, phistart, Rin[i], Zin[i])){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }
      double u,v;
      if (!state_to_uv(FLT, PAR, u, v)){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }
      // unwrap (theta only) and compute g in that chart
	  double umin, umax, vmin, vmax;
	  unwrap_uv_to_box_center(PAR, B, u, v, umin, umax, vmin, vmax);

	  double g[4];
	  g[0] = umin - u;
	  g[1] = u - umax;
	  g[2] = vmin - v;
	  g[3] = v - vmax;

      double du_dR,du_dZ,dv_dR,dv_dZ;
      if (!uv_grad_fd(FLT, PAR, phistart, Rin[i], Zin[i], epsR, epsZ, du_dR, du_dZ, dv_dR, dv_dZ)){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        return false;
      }

      for (int k=0;k<4;k++){
	    // same banded hinge gp as eval_all
	    double gp = 0.0;
	    double dgp_dg = 0.0;
	    if (g[k] >= 0.0) { gp = g[k]; dgp_dg = 1.0; }
	    else if (g[k] > -g_band) { gp = (g[k] + g_band); dgp_dg = 1.0; }
	    else continue; // fully inactive

	    double gR=0.0, gZ=0.0;
	    if (k==0){ gR = -du_dR; gZ = -du_dZ; }
	    if (k==1){ gR =  du_dR; gZ =  du_dZ; }
	    if (k==2){ gR = -dv_dR; gZ = -dv_dZ; }
	    if (k==3){ gR =  dv_dR; gZ =  dv_dZ; }

	    // chain rule through gp(g)
	    gR *= dgp_dg;
	    gZ *= dgp_dg;

	    const double t = lambda[4*i+k] + rho*gp;  // IMPORTANT: gp, not g

	    JTr[idxR(i)] += gR * t;
	    JTr[idxZ(i)] += gZ * t;

	    JTJ[idxR(i)*dim + idxR(i)] += rho * (gR*gR);
	    JTJ[idxR(i)*dim + idxZ(i)] += rho * (gR*gZ);
	    JTJ[idxZ(i)*dim + idxR(i)] += rho * (gZ*gR);
	    JTJ[idxZ(i)*dim + idxZ(i)] += rho * (gZ*gZ);
	  }
    }

    FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
    return true;
  };

  // ---- Augmented Lagrangian outer loop ----
  double phi=0.0, gmax=0.0;
  if (!eval_all(R, Z, rvec, gmax, phi)) return -1;

  for (int al_it=0; al_it<al_imax; al_it++)
  {
    // ---- LM inner loop ----
    for (int lm_it=0; lm_it<lm_imax; lm_it++)
    {
      double rnorm2 = 0.0;
      for (double v : rvec) rnorm2 += v*v;
      const double rnorm = std::sqrt(rnorm2);

	  if (lm_it == 0 || (lm_it % 5) == 0){
  	    ofs2 << "AL " << al_it
        << "  LM " << lm_it
        << "  rnorm=" << rnorm
        << "  gmax="  << gmax
        << "  phi="   << phi
        << "  mu="    << mu
        << "  rho="   << rho
        << endl;
      }

      if (rnorm < tol_rnorm && gmax < tol_gmax){
        if (final_cost) *final_cost = 2.0*phi;
        return 0;
      }

      if (!build_residual_normal_eq(R, Z, rvec)) return -1;
      if (!add_constraint_terms(R, Z)) return -1;

      A = JTJ;
      for (int k=0;k<dim;k++) A[k*dim + k] += mu;

      for (int k=0;k<dim;k++) b[k] = -JTr[k];

      std::vector<double> Acpy = A, bcpy = b;
      if (!solve_dense(Acpy, bcpy, dx, dim)){
        mu = std::min(1e18, mu*mu_up);
        continue;
      }

      // step norm (full dx)
	  double step2_full = 0.0;
	  for (int i=0;i<n;i++){
	    const double dR = dx[idxR(i)];
	    const double dZ = dx[idxZ(i)];
	    step2_full += dR*dR + dZ*dZ;
 	  }
	  if (std::sqrt(step2_full) < tol_step){
	    if (final_cost) *final_cost = 2.0*phi;
	    break;
	  }

	  // --- backtracking line search ---
	  double alpha = 1.0;
	  bool accepted = false;

	  for (int ls=0; ls<10; ls++){
	    Rtry = R; Ztry = Z;
	    for (int i=0;i<n;i++){
		  Rtry[i] += alpha * dx[idxR(i)];
		  Ztry[i] += alpha * dx[idxZ(i)];
	    }

	    double phi_try=0.0, gmax_try=0.0;
	    if (!eval_all(Rtry, Ztry, rtry, gmax_try, phi_try)){
		  // treat as "bad step": shrink and retry
		  alpha *= 0.5;
		  continue;
	    }

	    if (phi_try < phi){
	  	  // accept
		  R.swap(Rtry);
		  Z.swap(Ztry);
		  rvec.swap(rtry);
		  phi  = phi_try;
		  gmax = gmax_try;

		  mu = std::max(1e-18, mu*mu_dn);
		  accepted = true;
		  break;
	    }

	    alpha *= 0.5;
	  }

	  if (!accepted){
	    // no alpha worked: increase damping
	    mu = std::min(1e18, mu*mu_up);
	    continue;
	  }
    }

    // ---- multiplier update, maybe increase rho ----
	{
	const double Rsave=FLT.R, Zsave=FLT.Z, psisave=FLT.psi, phisave=FLT.phi;

	double gmax_now = 0.0;

	for (int i=0;i<n;i++){
		const int wi = sigma_index(i, n, box_dir);
		const BOX2D& B = W[wi];

		if (!set_state_from_RZ(FLT, PAR, phistart, R[i], Z[i])){
		FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave; 
		return -1;
		}

		double u,v;
		if (!state_to_uv(FLT, PAR, u, v)){
		FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave; 
		return -1;
		}

		// unwrap box interval consistently around current (u,v)
		double umin, umax, vmin, vmax;
		unwrap_uv_to_box_center(PAR, B, u, v, umin, umax, vmin, vmax);

		// constraints with UNWRAPPED interval
		double g[4];
		g[0] = umin - u;
		g[1] = u - umax;
		g[2] = vmin - v;
		g[3] = v - vmax;

		for(int k=0;k<4;k++){
		// classic AL update for inequality g<=0 with lambda>=0
		lambda[4*i+k] = std::max(0.0, lambda[4*i+k] + rho*g[k]);

		double gp = 0.0;
		if (g[k] >= 0.0) gp = g[k];
		else if (g[k] > -g_band) gp = (g[k] + g_band);

		gmax_now = std::max(gmax_now, gp);
		}
	}

	FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;

	gmax = gmax_now;
	if (gmax < tol_gmax) break;

	rho = std::min(1e18, rho * rho_up);
	mu  = std::max(1e-6, mu);
	}

	// recompute merit and residuals cleanly after updates
	if (!eval_all(R, Z, rvec, gmax, phi)) return -1;
  }

  if (final_cost) *final_cost = 2.0*phi;
  return (gmax < tol_gmax) ? 0 : -1;
}

// Find a UV seed in a box that is EFIT-valid and does NOT close early (smaller period).
static inline bool find_survivor_in_box(
  PARTICLE& FLT, IO& PAR, double phistart, int map_dir,
  int periode,
  const BOX2D& B, int Nu, int Nv,
  double& useed, double& vseed
){
  Nu = std::max(2, Nu); Nv = std::max(2, Nv);
  const double du = (B.umax - B.umin) / (Nu - 1);
  const double dv = (B.vmax - B.vmin) / (Nv - 1);
  const double eps2 = 1e-10; // (1e-5)^2

  const double Rsave=FLT.R, Zsave=FLT.Z, psisave=FLT.psi, phisave=FLT.phi;

  for(int iu=0; iu<Nu; iu++){
    const double u = B.umin + iu*du;
    for(int iv=0; iv<Nv; iv++){
      const double v = B.vmin + iv*dv;

      // uv -> RZ (dtfix conventions)
      double R0,Z0;
      if(!uv_to_RZ(FLT, PAR, u, v, R0, Z0)){ ofs2 << "Index (" << iu << "," << iv << ") Failed uv_to_RZ.\n"; continue;}

      // EFIT valid + inside scan region
      if(!set_state_from_RZ(FLT, PAR, phistart, R0, Z0)){ ofs2 << "Index (" << iu << "," << iv << ") Failed set_state_from_RZ.\n"; continue;}
      // if(!inside(FLT, PAR)){ ofs2 << "Index (" << iu << "," << iv << ") Failed inside(FLT, PAR).\n"; continue;}

      // reject smaller period: if you return to seed before 'periode', skip it
      bool ok = true;
      const double Rref = R0, Zref = Z0;

      for(int k=1; k<periode; k++){
        int chk = FLT.mapit(1, map_dir);   // one map step in correct direction
		if(k == 1){
		  if(chk < 0){ ok=false;
			ofs2 << "Index (" << iu << "," << iv << ") Rejecting seed u=" << u << " v=" << v << ": leaves system on first step.\n"; 
			break; } // left system right away
		}else{
		  if(chk < 0){ break; }            // left system later
		}
        const double dR = FLT.R - Rref, dZ = FLT.Z - Zref;
        if(dR*dR + dZ*dZ < eps2){ ok=false; break; } // closed early
      }
      if(!ok){
        FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
        continue;
      }

      useed = u;
      vseed = v;
      FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
      return true;
    }
  }

  FLT.R=Rsave; FLT.Z=Zsave; FLT.psi=psisave; FLT.phi=phisave;
  return false;
}

// -------------------- Boxes file reader --------------------
bool read_windows_2d(const LA_STRING& path, std::vector<BOX2D>& W, int expected_n,
                     PARTICLE& /*FLT*/, IO& /*PAR*/)
{
  W.clear();
  std::ifstream in(path);
  if(!in) return false;

  std::string line;
  while(std::getline(in, line)){
    std::istringstream ts(line);
    std::string first;
    if(!(ts >> first)) continue;
    if(first.size() && first[0] == '#') continue;

    double umin, umax, vmin, vmax;
    { std::istringstream iss(line); if(!(iss >> umin >> umax >> vmin >> vmax)) return false; }
    W.push_back(BOX2D{umin,umax,vmin,vmax});
    if(expected_n > 0 && (int)W.size() == expected_n) break;
  }
  return (expected_n <= 0) ? true : ((int)W.size() == expected_n);
}

bool fill_cycle_arrays(
  PARTICLE& FLT, double phistart, int periode, int map_dir,
  double R0, double Z0,
  Array<double,1>& R, Array<double,1>& Z,
  Array<double,1>& psi, Array<double,1>& theta, Array<double,1>& r
){
  FLT.R = R0; FLT.Z = Z0; FLT.phi = phistart;

  double psival;
  if(FLT.get_psi(FLT.R, FLT.Z, psival, FLT.phi) < 0) return false;
  FLT.psi = psival;

  R(0)=FLT.R; Z(0)=FLT.Z; psi(0)=FLT.psi;
  r(0)=FLT.get_r(); theta(0)=FLT.get_theta();

  for(int k=1;k<periode;k++){
    int chk = FLT.mapit(1, map_dir);
    if(chk < 0) return false;

    if(FLT.get_psi(FLT.R, FLT.Z, psival, FLT.phi) < 0) return false;
    FLT.psi = psival;

    R(k)=FLT.R; Z(k)=FLT.Z; psi(k)=psival;
    r(k)=FLT.get_r(); theta(k)=FLT.get_theta();
  }
  return true;
}

int classify_OX(
  PARTICLE& FLT, double phistart, int periode, int map_dir,
  double R0, double Z0
){
  FLT.R = R0 + 1e-5;
  FLT.Z = Z0 + 1e-5;
  FLT.phi = phistart;

  int chk = FLT.mapit(10*periode, map_dir);
  double dist = std::sqrt((R0-FLT.R)*(R0-FLT.R) + (Z0-FLT.Z)*(Z0-FLT.Z));
  if(dist > 1e-2) chk = -1;
  return (chk >= 0) ? 1 : -1;
}


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



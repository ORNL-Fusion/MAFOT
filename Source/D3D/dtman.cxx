// Program calculates unstable manifold of a hyp. fixed point for D3D-Drift with time dependent perturbation
// For stable manifold use reverse Integration (see MapDirection)
// For left or right hand-sided set sign of 'verschieb' in Parameterfile to - or + respectively
// Fortran Subroutines are used for the perturbations
// A.Wingen						12.01.09

// Input: 1: Parameterfile	2: File with fixed points	3: praefix(optional)
//			fixed points have to be in toroidal coordinates, not cylindrical!!!
// Output:	one file for each of the fixed points specified in the input-file, giving the respective manifold
//			log-file

// Define
//--------
#define program_name "dtman"

// Include
//--------
#include <andi.hxx>
#include <efit_class.hxx>
#include <d3d-drift.hxx>

// Prototypes
void find_start_values(Array<double,1>& xs, Array<double,1>& d, int periode, double& verschieb, 
					   double phistart, int MapDirection);
double distmod2pi(Array<double,1> x, Array<double,1> y);

// Switches
const int useparfile = 1;	// 0: additional parameters are set in the code		1: All parameters are read from file
const int trytoskip = 1;	// 0: stop after first wall contact		1: try to continue, may cause errors!!!
const int preventSmallSteps = 0;	// 0: code stops if step size dt < 1e-14	1: code continues with dt = 1e-10 as long as step size controll would reduce dt below 1e-10

// Golbal Parameters
EFIT EQD;
double GAMMA;
double eps0;
double Ix;

// Main Program
//--------------
int main(int argc, char *argv[])
{
// Variables
int i,k,periode,chk,plotchk,end,variate;
double xfix,yfix;
double t,dt,talt,dist;
double x,y,phi,psi;
int skipattempt,outside;
double omegac;

// Arrays
Array<double,1> xa(Range(1,2));
Array<double,1> xsa(Range(1,2));
Array<double,1> xitta(Range(1,2));
Array<double,1> xittolda(Range(1,2));
Array<double,1> da(Range(1,2));

vector<double> xvec,yvec,periodevec,dummy;

// Input file names
LA_STRING basename,startname;
LA_STRING praefix = "";
if(argc==4) praefix = "_" + LA_STRING(argv[3]);
if(argc>=3) 
{
	basename = LA_STRING(argv[1]);
	startname = LA_STRING(argv[2]);
}
if(argc<=2)
{
	cout << "No Input files -> Abort!" << endl; 
	exit(0);
}
basename = checkparfilename(basename);
startname = checkparfilename(startname);
LA_STRING parfilename = "_" + basename + ".dat";
LA_STRING name = startname + ".dat";

// Read parameter file
vector<double> startvec;
cout << "read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);

// Set starting parameters
const double minabs = 0.005;	//0.001
const double maxabs = 0.01;		//0.005
double verschieb = 1e-4;
double phistart = 0;
int MapDirection = 1;

const int kstart = 1;
const int kend = 30;
const int kstartpkt = kstart;

double Ekin = 10;		// kinetic Energy in [keV]
double lambda = 0.1;	// ratio of kinetic energy in R direction to total kinetic energy, simply estimated; ????? Inluence on results ????? 

if(useparfile==1) 
{
	cout << "All parameters are read from file" << endl;
	verschieb = startvec[0];
	phistart = startvec[7];
	MapDirection = int(startvec[8]);	
	Ekin = startvec[18];
	lambda = startvec[19];
}

// read fixed points
readfile(name,xvec,yvec,periodevec,dummy);

// Output
LA_STRING filenameout,type,dir;
if(MapDirection==1) type = "_unst";
else type = "_st";
if(verschieb>0) dir = "r";
else dir = "l";
ofstream out;
out.precision(16);
vector<LA_STRING> var(3);
var[0] = "theta[rad]";  var[1] = "r[m]";  var[2] = "psi";

// log file
ofs2.open("log_" + LA_STRING(program_name) + type + dir + praefix + ".dat");

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prep_perturbation();

// Prepare particles
if(sigma==0)
{
	cout << "Field lines are calculated" << endl;
	ofs2 << "Field lines are calculated" << endl;
	GAMMA = 1;	omegac = 1;	 eps0 = 1;	Ix = 1;
}
else
{
	if(Zq>=1) // Ions
	{
		GAMMA = 1 + Ekin/(E0p*Massnumber);	// relativistic gamma factor 1/sqrt(1-v^2/c^2)
		omegac = e*EQD.Bt0/(mp*Massnumber);	// normalized gyro frequency (SI-System)
		cout << "Ions are calculated" << endl;
		ofs2 << "Ions are calculated" << endl;
	}
	else // Electrons
	{
		Zq = -1;	// default!
		GAMMA = 1 + Ekin/(E0e*Massnumber);	// relativistic gamma factor 1/sqrt(1-v^2/c^2)
		omegac = e*EQD.Bt0/(me*Massnumber);	// normalized gyro frequency (SI-System)
		cout << "Electrons are calculated" << endl;
		ofs2 << "Electrons are calculated" << endl;
	}
	eps0 = c*c/omegac/omegac/EQD.R0/EQD.R0;	// normalized rest energy
	Ix = -0.5/double(Zq)*eps0*((lambda*(GAMMA-1)+1)*(lambda*(GAMMA-1)+1)-1);
	cout << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
	ofs2 << "kin. Energy: Ekin= " << Ekin << "keV" << "\t" << "rel. gamma-factor: gamma= " << GAMMA << endl;
}

// loop for all fixed points
const int psize = 10;
for(i=0;i<int(xvec.size());i++)
{
	xfix = xvec[i];	 yfix = yvec[i];	periode = int(periodevec[i]);
	xsa(1) = xfix;	 xsa(2) = yfix;	

	find_start_values(xsa,da,periode,verschieb,phistart,MapDirection); 
	
	ofs2 << "Period: " << periode << "\t" << "Shift: " << verschieb << endl;
	ofs2 << "fixed point: theta= " << xfix << "\t" << "r= " << yfix << endl;
	ofs2 << "k goes from k= " << kstart << " to k= " << kend << endl;
	ofs2 << "Distances: Min= " << minabs << "\t" << "Max= " << maxabs << endl;
	ofs2 << endl;

	// additional parameters for IO
	parstruct * parvec = new parstruct[psize];
	parvec[0].name = "Period";				parvec[0].wert = periode;
	parvec[1].name = "fixed point theta";	parvec[1].wert = xfix;
	parvec[2].name = "fixed point r";		parvec[2].wert = yfix;
	parvec[3].name = "Shift in theta";		parvec[3].wert = verschieb;
	parvec[4].name = "Minimal distance";	parvec[4].wert = minabs;
	parvec[5].name = "Maximal distance";	parvec[5].wert = maxabs;
	parvec[6].name = "phistart";			parvec[6].wert = phistart;
	parvec[7].name = "MapDirection";		parvec[7].wert = MapDirection;
	parvec[8].name = "Ekin";				parvec[8].wert = Ekin;
	parvec[9].name = "energy ratio lambda";	parvec[9].wert = lambda;

	// Output
	filenameout = "man" + type + dir + LA_STRING(periode) + praefix + "_" + LA_STRING(i) + ".dat";
	outputtest(filenameout);
	out.open(filenameout);
	writeiodata(out,var,parvec,psize,parfilename);

	// calculate Manifold
	end = 0;  skipattempt = 0;  outside = 0;
	plotchk = 0;
	variate = 0;
	xittolda = xsa + da;
	bndy[0] = 0.85;  bndy[1] = 2.5;  bndy[2] = -1.58;  bndy[3] = 1.58;	// Extend Boundary to efit boundary
	dt = 0.1;

	for(k=kstart;k<=kend;k++)
	{
		t = 0;
		talt = t;
		//dt = 0.1;

		while(t<1)
		{
			t += dt;
			if(t>1) 
			{
				dt = 1 - talt;
				t = 1;
			}
			xa = xsa + t*da;
			xitta = xa;

			phi = phistart;
			chk = mapit(xitta(1),xitta(2),phi,psi,k*2*periode,MapDirection);
			if(chk<0) 	// outside wall, try to skip outside part, no step size management
			{
				skipattempt = 1; 
				if(outside>3) {end = 1; ofs2 << "final wall hit" << endl; break;} // stop after 3 skip attempts
				if(trytoskip==1) continue;
				else {end = 1; ofs2 << "wall hit" << endl; break;}
			}
			if(skipattempt==1) {skipattempt = 0; outside += 1; variate = 1;}

			// step size management
			dist = distmod2pi(xitta,xittolda);
			if(variate==0 && dist>maxabs)
			{
				t = talt;
				dt *= 0.5;
				if(dt < 1e-8 && preventSmallSteps == 1)	// forced to execute one step with dt = 1e-6
					{dt = 1e-6; variate = 1; ofs2 << "Step size too small -> try to skip" << endl;}
				if(dt < 1e-14){end = 1; ofs2 << "Step size cannot be further decreased" << endl; break;}
				continue;
			}

			// Output 
			if(plotchk==0) {ofs2 << "Start plotting..." << endl; plotchk = 1;}
			x = xitta(1);
			y = xitta(2);
			out << x << "\t" << y << "\t" << psi << endl;
		
			// step size management part 2
			if(dist<minabs)
			{
				dt *= 2.5;
			}

			talt = t;
			xittolda = xitta;
			variate = 0;
		} // end while

		if(end==1) break;
		ofs2 << "k= " << k << "\t" << flush; 
	} // end for k
	ofs2 << endl;
	out.close();
	bndy[0] = 1.016;  bndy[1] = 2.4;  bndy[2] = -1.367;  bndy[3] = 1.36;	// Reset Boundary
} // end for i
ofs2 << "Program terminates normally" << endl;
return 0;
} //end of main

//------------------------ End of Main -------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------

//---------- find_start_values ----------------
void find_start_values(Array<double,1>& xsa, Array<double,1>& da, int periode, double& verschieb, 
					   double phistart, int MapDirection)
{
int chk;
double laenge,scalar,phi,psi;

Array<double,1> d2a(Range(1,2));
Array<double,1> fixa(Range(1,2));
Array<double,1> xsitta(Range(1,2));

fixa = xsa;
xsitta = xsa;

da(2) = 0; da(1)= verschieb;	// shift in theta

// adjust shift that way, that xsa and xsitta are between 0.0005 and 0.001
int end = 0;
while(end<50)
{
	xsa = fixa + da;
	xsitta = xsa;
	phi = phistart;
	chk = mapit(xsitta(1),xsitta(2),phi,psi,2*periode,MapDirection);
	if(chk<0) {end += 1; da *= 0.5; continue;}

	laenge=distmod2pi(xsa,xsitta);
	end += 1;
	if(laenge > 0.001) {da *= 0.5; continue;}
	if(laenge < 0.0005){da *= 1.1; continue;}
	break;
}
if(end==50) ofs2 << "Shift adjustment not successfull! " << laenge << endl;
verschieb = da(1);

//ofs2 << xsa(1) << "\t" << xsa(2) << "\t" << xsitta(1) << "\t" << xsitta(2) << "\t" << da(1) << "\t" << da(2) << endl;

// adjust direction of d
end = 0;  scalar = 0;
while(1-scalar > 1e-5)
{
	// set d
	da = xsitta - fixa;
	if(fabs(da(1)) > 4) {da(1) -= sign(da(1))*pi2;}	// possible 2pi jump between xsitt and fix.
	da *= fabs(verschieb) / sqrt(da(1)*da(1)+da(2)*da(2));		// set length of d to shift

	// check new direction
	xsa = fixa + da;
	xsitta = xsa;
	phi = phistart;
	chk = mapit(xsitta(1),xsitta(2),phi,psi,2*periode,MapDirection);
	if(chk<0) {ofs2 << "Error during direction adjustment: " << chk << endl; exit(0);}

	// terminate condition.
	d2a = xsitta - xsa;
	if(fabs(d2a(1)) > 4) {d2a(1) -= sign(d2a(1))*pi2;}
	scalar = (da(1)*d2a(1)+da(2)*d2a(2)) / sqrt(da(1)*da(1)+da(2)*da(2)) / sqrt(d2a(1)*d2a(1)+d2a(2)*d2a(2));

	end += 1;
	if(end>50) {ofs2 << "Direction adjustment not successfull" << endl; break;}

	//ofs2 << xsa(1) << "\t" << xsa(2) << "\t" << xsitta(1) << "\t" << xsitta(2) << "\t" << da(1) << "\t" << da(2) << endl;
}
ofs2 << "Effort of direction adjustment: " << end << endl;

// set d to be vector from xs to xsitt
da = xsitta - xsa;
if(fabs(da(1)) > 4) {da(1) -= sign(da(1))*pi2;}	// possible 2pi jump between xsitt and fix.
ofs2 << "Length of d: " << sqrt(da(1)*da(1)+da(2)*da(2)) << endl;
}


//---------- distmod2pi ----------------------
double distmod2pi(Array<double,1> x, Array<double,1> y)
{
// theta coordinate only determined modulo pi2
// minimize point distance with respect to modulo 2pi
int n;
double dist,distaus;
double mod;

distaus = 1e+10;

for(n=-1;n<=1;n++)
{
	mod = n*pi2;
	//dist = fabs(x-y+mod);
	dist = sqrt((x(1)-y(1)+mod)*(x(1)-y(1)+mod) + (x(2)-y(2))*(x(2)-y(2)));
	if(dist<distaus) distaus = dist;
}

return distaus;
}

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

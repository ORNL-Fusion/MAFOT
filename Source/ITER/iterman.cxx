// Program calculates unstable manifold of a hyp. fixed point for ITER with time dependent perturbation
// For stable manifold use reverse Integration (see MapDirection)
// For left or right hand-sided set sign of 'verschieb' in Parameterfile to - or + respectively
// Fortran Subroutines are used for the perturbations
// A.Wingen						16.06.11

// Input: 1: Parameterfile	2: File with fixed points	3: praefix(optional)
//			fixed points have to be in toroidal coordinates, not cylindrical!!!
// Output:	one file for each of the fixed points specified in the input-file, giving the respective manifold
//			log-file

// Define
//--------
#define program_name "iterman"
#define ITER

// Include
//--------
#include <mafot.hxx>
#include <iter.hxx>

// Prototypes
void find_start_values(Array<double,1>& xs, Array<double,1>& d, int periode, double& verschieb, 
					   double phistart, int MapDirection, PARTICLE& FLT);
inline double abstand(PARTICLE& FLT1, PARTICLE& FLT2);
inline double abstand(Array<double,1>& x1, Array<double,1>& x2);

// Switches
const int trytoskip = 1;	// 0: stop after first wall contact		1: try to continue, may cause errors!!!
const int preventSmallSteps = 0;	// 0: code stops if step size dt < 1e-14	1: code continues with dt = 1e-10 as long as step size controll would reduce dt below 1e-10

// Golbal Parameters

// Main Program
//--------------
int main(int argc, char *argv[])
{
// Variables
int i,k,periode,chk,plotchk,end,variate;
int skipattempt,outside;
double xfix,yfix;
double t,dt,talt,dist;
EFIT EQD;

// Arrays
Array<double,1> xa(Range(1,2));
Array<double,1> xsa(Range(1,2));
Array<double,1> da(Range(1,2));

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
cout << "Read Parameterfile " << parfilename << endl;
ofs2 << "Read Parameterfile " << parfilename << endl;
IO PAR(EQD,parfilename,10);

// Set starting parameters
const double minabs = 0.001;	//0.001
const double maxabs = 0.005;	//0.005
const int kstart = 1;
const int kend = 30;
const int kstartpkt = kstart;

// read fixed points
Array<double,2> data;
readfile(name,6,data);

// Output
LA_STRING filenameout,type,dir;
if(PAR.MapDirection==1) type = "_unst";
else type = "_st";
if(PAR.verschieb>0) dir = "r";
else dir = "l";
ofstream out;
out.precision(16);
vector<LA_STRING> var(5);
var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "psi";  var[3] = "theta[rad]";  var[4] = "r[m]";

// log file
ofs2.open("log_" + LA_STRING(program_name) + type + dir + praefix + ".dat");

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
ofs2 << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;

// Prepare Perturbation
prep_perturbation(EQD,PAR);

// Prepare particles
PARTICLE FLT(EQD,PAR);
PARTICLE FLTold(EQD,PAR);

// loop for all fixed points
for(i=1;i<=data.rows();i++)
{
	xfix = data(i,1);	 yfix = data(i,2);	periode = int(data(i,3));
	xsa(1) = xfix;	 xsa(2) = yfix;	

	find_start_values(xsa,da,periode,PAR.verschieb,PAR.phistart,PAR.MapDirection,FLT); 
	
	ofs2 << "Period: " << periode << "\t" << "Shift: " << PAR.verschieb << endl;
	ofs2 << "fixed point: R = " << xfix << "\t" << "Z = " << yfix << endl;
	ofs2 << "k goes from k = " << kstart << " to k= " << kend << endl;
	ofs2 << "Distances: Min = " << minabs << "\t" << "Max= " << maxabs << endl;
	ofs2 << endl;

	// additional parameters for IO
	PAR.pv[0].name = "Period";				PAR.pv[0].wert = periode;
	PAR.pv[1].name = "fixed point theta";	PAR.pv[1].wert = xfix;
	PAR.pv[2].name = "fixed point r";		PAR.pv[2].wert = yfix;
	PAR.pv[3].name = "Shift in theta";		PAR.pv[3].wert = PAR.verschieb;
	PAR.pv[4].name = "Minimal distance";	PAR.pv[4].wert = minabs;
	PAR.pv[5].name = "Maximal distance";	PAR.pv[5].wert = maxabs;
	PAR.pv[6].name = "phistart";			PAR.pv[6].wert = PAR.phistart;
	PAR.pv[7].name = "MapDirection";		PAR.pv[7].wert = PAR.MapDirection;
	PAR.pv[8].name = "Ekin";				PAR.pv[8].wert = PAR.Ekin;
	PAR.pv[9].name = "energy ratio lambda";	PAR.pv[9].wert = PAR.lambda;

	// Output
	filenameout = "man" + type + dir + LA_STRING(periode) + praefix + "_" + LA_STRING(i) + ".dat";
	outputtest(filenameout);
	out.open(filenameout);
	PAR.writeiodata(out,bndy,var);

	// calculate Manifold
	end = 0;  skipattempt = 0;  outside = 0;
	plotchk = 0;
	variate = 0;
	xa = xsa + da;
	FLTold.R = xa(1);
	FLTold.Z = xa(2);
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
			//xitta = xa;

			FLT.R = xa(1);
			FLT.Z = xa(2);
			FLT.phi = PAR.phistart;
			chk = FLT.mapit(k*2*periode,PAR.MapDirection);
			if(chk<0) 	// outside wall, try to skip outside part, no step size management
			{
				skipattempt = 1; 
				if(outside>3) {end = 1; ofs2 << "final wall hit" << endl; break;} // stop after 3 skip attempts
				if(trytoskip==1) continue;
				else {end = 1; ofs2 << "wall hit" << endl; break;}
			}
			if(skipattempt==1) {skipattempt = 0; outside += 1; variate = 1;}

			// step size management
			dist = abstand(FLT,FLTold);
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
			out << FLT.R << "\t" << FLT.Z << "\t" << FLT.psi << "\t" << FLT.get_theta() << "\t" << FLT.get_r() << endl;
		
			// step size management part 2
			if(dist<minabs)
			{
				dt *= 2.5;
			}

			talt = t;
			FLTold = FLT;
			variate = 0;
		} // end while

		if(end==1) break;
		ofs2 << "k= " << k << "\t" << flush; 
	} // end for k
	ofs2 << endl;
	out.close();
} // end for i
ofs2 << "Program terminates normally" << endl;
cout << "Program terminates normally" << endl;
return 0;
} //end of main

//------------------------ End of Main -------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------

//---------- find_start_values ----------------
void find_start_values(Array<double,1>& xsa, Array<double,1>& da, int periode, double& verschieb, 
					   double phistart, int MapDirection, PARTICLE& FLT)
{
int chk;
double laenge,scalar,phi,psi;

Array<double,1> d2a(Range(1,2));
Array<double,1> fixa(Range(1,2));
Array<double,1> xsitta(Range(1,2));

fixa = xsa;
xsitta = xsa;

da(2) = 0; da(1)= verschieb;	// shift in R

// adjust shift that way, that xsa and xsitta are between 0.0005 and 0.001
int end = 0;
while(end<50)
{
	xsa = fixa + da;
	FLT.R = xsa(1);
	FLT.Z = xsa(2);
	FLT.phi = phistart;
	chk = FLT.mapit(2*periode,MapDirection);
	if(chk<0) {end += 1; da *= 0.5; continue;}
	xsitta(1) = FLT.R;
	xsitta(2) = FLT.Z;

	laenge = abstand(xsa,xsitta);
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
	//if(fabs(da(1)) > 4) {da(1) -= sign(da(1))*pi2;}	// possible 2pi jump between xsitt and fix.
	da *= fabs(verschieb) / sqrt(da(1)*da(1)+da(2)*da(2));		// set length of d to shift

	// check new direction
	xsa = fixa + da;
	FLT.R = xsa(1);
	FLT.Z = xsa(2);
	FLT.phi = phistart;
	chk = FLT.mapit(2*periode,MapDirection);
	if(chk<0) {ofs2 << "Error during direction adjustment: " << chk << endl; exit(0);}
	xsitta(1) = FLT.R;
	xsitta(2) = FLT.Z;

	// terminate condition.
	d2a = xsitta - xsa;
	//if(fabs(d2a(1)) > 4) {d2a(1) -= sign(d2a(1))*pi2;}
	scalar = (da(1)*d2a(1)+da(2)*d2a(2)) / sqrt(da(1)*da(1)+da(2)*da(2)) / sqrt(d2a(1)*d2a(1)+d2a(2)*d2a(2));

	end += 1;
	if(end>50) {ofs2 << "Direction adjustment not successfull" << endl; break;}

	//ofs2 << xsa(1) << "\t" << xsa(2) << "\t" << xsitta(1) << "\t" << xsitta(2) << "\t" << da(1) << "\t" << da(2) << endl;
}
ofs2 << "Effort of direction adjustment: " << end << endl;

// set d to be vector from xs to xsitt
da = xsitta - xsa;
//if(fabs(da(1)) > 4) {da(1) -= sign(da(1))*pi2;}	// possible 2pi jump between xsitt and fix.
ofs2 << "Length of d: " << sqrt(da(1)*da(1)+da(2)*da(2)) << endl;
}

//---------- abstand ----------------------
// Calculates distance between (R,Z) positions of two PARTICLE objects
double abstand(PARTICLE& FLT1, PARTICLE& FLT2)
{
const double dR = FLT1.R - FLT2.R;
const double dZ = FLT1.Z - FLT2.Z;
return sqrt(dR*dR + dZ*dZ);
}

//---------- abstand ----------------------
// Calculates distance between two (R,Z) arrays
double abstand(Array<double,1>& x1, Array<double,1>& x2)
{
const double dR = x1(1) - x2(1);
const double dZ = x1(2) - x2(2);
return sqrt(dR*dR + dZ*dZ);
}

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

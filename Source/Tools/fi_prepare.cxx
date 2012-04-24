// Program precalculates total magnetic field on EFIT grid and 360 toroidal slices of current filaments
// for D3D-Drift with current filament perturbations
// A.Wingen						23.2.09

// Input: 1: Parameterfile	2: praefix(optional)
// Output:	filament field file



// Define
//--------
//#define BZ_DEBUG
#define program_name "fi_prepare"

// Include
//--------
#include <andi.hxx>
#include <efit_class.hxx>
//#include <d3d-drift-fi.hxx>

// Prototypes 
void readiodata(char* name, vector<double>& vec);
void prepare_filament_field(Array<double,4>& field, LA_STRING filename);
double filament_fields_ex(double x, double y, double z, Array<double,2>& data, int N, double I, double dIdt,
					   double& B_x, double& B_y, double& B_z, double& dAdt, double cosp, double sinp);

// Switches
int useFilament = 0;	// 0: no	>= 1: Number of Filaments to be included 

// Golbal Parameters
const double pi = LA_PI;
const double pi2 = 2*pi;
const double rTOd = 180.0/pi;	// rad to deg

EFIT EQD;
Array<double,3> filament_data(Range(1,1),Range(0,1),Range(1,5));	// default size
Array<double,4> field;	// default constructed


// Main Program
//--------------
int main(int argc, char *argv[])
{
// Input file names
LA_STRING basename,filename;
LA_STRING praefix = "";
if(argc==3) praefix = "_" + LA_STRING(argv[2]);
if(argc>=2) basename = LA_STRING(argv[1]);
else	// No Input: Abort
{
	cout << "No Input -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
LA_STRING parfilename = "_" + basename + ".dat";

// Read parameter file
vector<double> startvec;
cout << "Read Parameterfile " << parfilename << endl;
readiodata(parfilename, startvec);
cout << "done" << endl;

// Read EFIT-data
EQD.ReadData(EQD.Shot,EQD.Time);
cout << "Shot: " << EQD.Shot << "\t" << "Time: " << EQD.Time << "ms" << endl;
cout << endl;

// Calculate total field and write to file
filename = "filament" + praefix + ".dat";
prepare_filament_field(field,filename);

return 0;
} //end of main
 

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

// ----------------- readiodata ------------------------------------------------------------
void readiodata(char* name, vector<double>& vec)
{
LA_STRING input;	// !!! LA_STRING reads entire line, string reads only one word !!!

// Get ShotNr and ShotTime from Parameterfile
ifstream in;
in.open(name);
in >> input;	// Skip first line
in >> input;
EQD.Shot = input.mid(9,6);	// 6 characters starting at index 9 of input string
EQD.Time = input.mid(22,4); // 4 characters starting at index 22 of input string

// Get Parameters
readparfile(name,vec);

//// Set switches
//which_target_plate = int(vec[11]);
//create_flag = int(vec[12]);
//
//if(vec[13]>1) 
//{
//	cout << "Coil flags are undefined in file " << name << endl;
//	ofs2 << "Coil flags are undefined in file " << name << endl;
//}
//else
//{
//	useFcoil = int(vec[13]);
//	useCcoil = int(vec[14]);
//	useIcoil = int(vec[15]);
//}
//
//if(vec.size()>=22)
//{
//	sigma = int(vec[16]);
//	Zq = int(vec[17]);
//}
//else {cout << "Fail to read particle parameters" << endl; exit(0);}

if(vec.size()>=23)
{
	useFilament = int(vec[20]);
}
else {cout << "No current filaments included" << endl;}
}

//------------------ prepare_filament_field ----------------------------------------------------
// Reads filament files and calculates total field on EFIT grid
// phi is discretized in 1 deg. steps, according to integrator step size  -->  no interpolation in phi!!!!!
// Bicubic interpolation on R and Z.
void prepare_filament_field(Array<double,4>& field, LA_STRING filename)
{
// Variables
int i,j,k;
ofstream out;

// firstIndex: 1 = B_X, 2 = B_Y, 3 = B_Z;   secondIndex: phi(0,...,359);   thirdIndex: R(0,1,...,NR,NR+1);   forthIndex: Z(0,1,...,NZ,NZ+1); 
// first and last index in R and Z are over the boundary to calculate derivatives at boundary
field.resize(Range(1,3),Range(0,359),Range(0,EQD.NR+1),Range(0,EQD.NZ+1));

// Read current-filament paths
//-----------------------------
ifstream in;
LA_STRING filamentfile,line;
string word;

Array<double,2> data;

if(useFilament>0) 
{
	//write field to file
	out.open(filename);
	out.precision(16);

	cout << "Number of current filaments used: " << useFilament << endl; 
//	ofs2 << "Number of current filaments used: " << useFilament << endl; 
	out << "# Number of current filaments used: " << useFilament << endl; 
	cout << "Filament currents: ";
//	ofs2 << "Filament currents: ";
	out << "# Filament currents: ";
}
for(i=1;i<=useFilament;i++)
{
	filamentfile = "filament" + LA_STRING(i) + ".in";
	readfile(filamentfile,5,data);

	// Store data in filament_data -> Adjust size if necessary 
	if(data.rows()>=filament_data.cols()) filament_data.resizeAndPreserve(useFilament,data.rows()+1,5);
	for(int j=1;j<=data.rows();j++) for(int k=1;k<=5;k++) filament_data(i,j,k) = data(j,k);
	filament_data(i,0,1) = data.rows();	// number of rows in data is stored at col-index 0 and depth-index 1
	data.free();	// data no longer used; prepare data for next use

	// read current from filament file -> store it in filament_data(i,0,2)
	in.open(filamentfile);	// location in the file needed !!!!!!!!!
	in >> line;
	in >> line;
	in >> word;	in >> word;
	in >> filament_data(i,0,2); 
	in.close();
	cout << i << ": " << filament_data(i,0,2) << "A" << "\t";
//	ofs2 << i << ": " << filament_data(i,0,2) << "A" << "\t";
	out << i << ": " << filament_data(i,0,2) << "A" << "\t";
}
if(useFilament>0) 
{
	cout << endl;
//	ofs2 << endl;
	out << endl << "#" << endl;
}

// Calculate total field of filaments on grid
//--------------------------------------------
double R,phi,X,Y,Z;
double bx,by,bz;
double B_X,B_Y,B_Z;
double sinp,cosp;
double dummy,I;
Range all = Range::all();
Array<double,2> slice;
Array<double,2> ddR(Range(1,EQD.NR),Range(1,EQD.NZ)),ddZ(Range(1,EQD.NR),Range(1,EQD.NZ)),d2dRdZ(Range(1,EQD.NR),Range(1,EQD.NZ));

// set grid
cout << "Calculate total field..." << endl;
//ofs2 << "Calculate total field..." << endl;
for(k=0;k<360;k++)
{
	phi = double(k)/rTOd;
	sinp = sin(phi);
	cosp = cos(phi);

	for(i=0;i<=EQD.NR+1;i++)
	{
		if(i==0) R = EQD.R(1) - EQD.dR;
		if(i==EQD.NR+1) R = EQD.R(EQD.NR) + EQD.dR;
		if(i>0 && i<=EQD.NR) R = EQD.R(i);
		X = R*cosp;
		Y = R*sinp;

		for(j=0;j<=EQD.NZ+1;j++) 
		{
			if(j==0) Z = EQD.Z(1) - EQD.dZ;
			if(j==EQD.NZ+1) Z = EQD.Z(EQD.NZ) + EQD.dZ;
			if(j>0 && j<=EQD.NZ) Z = EQD.Z(j);
			B_X = 0;	B_Y = 0;	B_Z = 0;	

			// Sum up fields of any current filament
			for(int l=1;l<=useFilament;l++)	// no filament fields, if useFilament == 0
			{
				bx = 0;	by = 0;	bz = 0;
				//I = get_current(l,0,dummy);	// time is fixed at t = 0; dummy is used for dIdt
				I = filament_data(l,0,2);	// current is stored in here, see prep_perturbation
				// filament_data(l,0,1) gives number of rows of data; index 0 is not used in filament_fields_ex
				slice.reference(filament_data(l,all,all));
				dummy = filament_fields_ex(X,Y,Z,slice,int(filament_data(l,0,1)),I,0,	// dIdt is set to 0 and fixed <- constant current
										   bx,by,bz,dummy,cosp,sinp);			// dummy is used instead of A_phi and dAdt
				B_X += bx;
				B_Y += by;
				B_Z += bz;
			}

			// Store total field
			field(1,k,i,j) = B_X;
			field(2,k,i,j) = B_Y;
			field(3,k,i,j) = B_Z;

			out << B_X << "\t" << B_Y << "\t" << B_Z << endl;

		} // end for j
	} // end for i
	//if(k%10==0) cout << k << "\t" << flush;
} // end for k
//cout << endl;
cout << "done" << endl;
//ofs2 << "done" << endl;
}

//------------- filament_fields_ex ----------------------------------------
// returns toroidal component A_phi of vector potential at position (x, y, z) at time t.
// data = data(1..N, 1..5) contains the pathway of the filament with columns X, Y, Z, R, phi
// Current I(t) in [A] and its time derivative dIdt(t) has to be provided. t is not needed!
// data usually has the form (R, Z, phi, psi) after reading from file -> has to be converted in main program
// before calling this function
// B_x, B_y and B_z are the components of the magnetic field.
// dAdt is the time derivative of A_phi
// Uses exact solution for vector potential, NO Taylor-series
// sinp = sin(phi) and cosp = cos(phi); phi is the toroidal angle of position, phi is not needed
double filament_fields_ex(double x, double y, double z, Array<double,2>& data, int N, double I, double dIdt,
					   double& B_x, double& B_y, double& B_z, double& dAdt, double cosp, double sinp)
{
int i,vor;
double dist,dist2;
double jx,jy,jz;
double dx1,dy1,dz1;
double dx2,dy2,dz2;
double d1,d2;
double dAx,dAy,dAz,A_part;
double A_x,A_y,A_phi;
double f;

B_x = 0;	B_y = 0;	B_z = 0;	
A_x = 0;	A_y = 0;

// calculate fields
for(i=1;i<N;i++)
{
	// Direction of current, without I
	jx = data(i+1,1) - data(i,1);
	jy = data(i+1,2) - data(i,2);
	jz = data(i+1,3) - data(i,3);
	dist = sqrt(jx*jx + jy*jy + jz*jz);
	jx /= dist;
	jy /= dist;
	jz /= dist;

	// Distance between position and filament
	dx1 = x - data(i,1);
	dy1 = y - data(i,2);
	dz1 = z - data(i,3);
	dx2 = x - data(i+1,1);
	dy2 = y - data(i+1,2);
	dz2 = z - data(i+1,3);

	dist = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);	
	dist2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

	// Set help variables
	d1 = dist - (jx*dx1 + jy*dy1 + jz*dz1);
	d2 = dist2 - (jx*dx2 + jy*dy2 + jz*dz2);
	vor = 1;

	// Help variables can get 0 <- not allowed --> Use another antiderivative (differs only on a constant -> same A)
	if(d1 <= 0 || d2 <= 0)
	{
		d1 = dist + (jx*dx1 + jy*dy1 + jz*dz1);
		d2 = dist2 + (jx*dx2 + jy*dy2 + jz*dz2);
		vor = -1;
	}

	// Vector potential, still without I !!!!
	A_part = vor*log(d2 / d1);
	A_x += A_part*jx;
	A_y += A_part*jy;
	//A_z += A_part*jz;

	// magnetic field
	dAx = vor*((dx2/dist2 - jx)/d2 - (dx1/dist - jx)/d1);
	dAy = vor*((dy2/dist2 - jy)/d2 - (dy1/dist - jy)/d1);
	dAz = vor*((dz2/dist2 - jz)/d2 - (dz1/dist - jz)/d1);

	B_x += jz*dAy - jy*dAz;
	B_y += jx*dAz - jz*dAx;
	B_z += jy*dAx - jx*dAy;
}
// transform A_x and A_y to A_phi
A_phi = A_y*cosp - A_x*sinp;

// Multiply with constants: mu0/4pi*I/L, (mu0/4pi = 1e-7)
f = 1e-7;
dAdt = f*dIdt*A_phi;

f *= I;
A_phi *= f;

B_x *= f;
B_y *= f;
B_z *= f;

return A_phi;
}


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



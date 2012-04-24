// Program gets intersections of current filament (struct-file) with Poincaré section.

// Define
//--------
//#define BZ_DEBUG
#define program_name "get_current"

// Include
//--------
#include <andi.hxx>
//#include <efit_class.hxx>
//#include <d3d-drift.hxx>

// Prototypes 
typedef struct {string name; double wert;} parstruct;
void writeiodata(ofstream& out, vector<LA_STRING>& var, parstruct * pv, int psize, char* name);

// Switches

// Golbal Parameters
const double pi = LA_PI;
const double pi2 = 2*pi;
const double rTOd = 180.0/pi;	// rad to deg

//EFIT EQD;

// Main Program
//--------------
int main(int argc, char *argv[])
{
int i;

// Input file names
LA_STRING basename,pointname;
LA_STRING praefix = "";
//if(argc==4) praefix = "_" + LA_STRING(argv[3]);
if(argc>=3) {basename = LA_STRING(argv[1]); pointname = LA_STRING(argv[2]);}
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
basename = checkparfilename(basename);
pointname = checkparfilename(pointname);
praefix = pointname.mid(7);
pointname = pointname + ".dat";
LA_STRING parfilename = "_" + basename + ".dat";

// Read parameter file
vector<double> startvec;
cout << "Read Parameterfile " << parfilename << endl;
readparfile(parfilename, startvec);

double phistart = startvec[7];

// Read struct file
Array<double,2> data;
readfile(pointname,5,data);

// additional parameters for IO
int psize = 2;
parstruct * parvec = new parstruct[psize];
parvec[0].name = "Original struct Filename: " + STRING(pointname);	parvec[0].wert = 0;
parvec[1].name = "phistart";										parvec[1].wert = phistart;

// Output
LA_STRING filenameout = "current" + praefix + ".dat";
outputtest(filenameout);
ofstream out(filenameout);
out.precision(16);
vector<LA_STRING> var(3);
var[0] = "R[m]";  var[1] = "Z[m]";  var[2] = "phi[rad]";
writeiodata(out,var,parvec,psize,parfilename);

// Find intersection points
for(i=1;i<=data.rows();i++)
{
	if(fabs(modulo2pi(data(i,5) - phistart + pi) - pi) < 1e-6) out << data(i,4) << "\t" << data(i,3) << "\t" << data(i,5) << endl;
}

return 0;
} //end of main
 

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

void writeiodata(ofstream& out, vector<LA_STRING>& var, parstruct * pv, int psize, char* name)
{
int i;
out << "# " << program_name << endl;
out << "#-------------------------------------------------" << endl;
out << "### Parameterfile: " << name << endl;
out << "#-------------------------------------------------" << endl;
out << "### Switches:" << endl;
out << "#-------------------------------------------------" << endl;
out << "### Global Parameters:" << endl;
out << "#-------------------------------------------------" << endl;
out << "### additional Parameters:" << endl;
for(i=0;i<psize;++i)
{
	out << "# " << pv[i].name << ": " << pv[i].wert << endl;
}
out << "#-------------------------------------------------" << endl;
out << "### Data:" << endl;
out << "# ";
for(i=0;i<int(var.size());i++) out << var[i] << "     ";
out << endl;
out << "#" << endl;
}
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



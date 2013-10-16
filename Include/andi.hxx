// Helpfull functions used in all programs
// a(ll thats) n(eedful) d(irectly) i(mplemented)
// a(lles) n(ützliche) d(irekt) i(mplementiert)
// Last modified 7.6.11


// Define
//--------
#ifdef USE_MPI
#define EXIT MPI::COMM_WORLD.Abort(0)
#else
#define EXIT exit(0)
#endif

#ifndef program_name		// checks if program_name is defined
#define program_name "default"	// if not, set program_name to default
#endif						// end

// Include
//--------
#include <la_string.hxx>

#ifdef linux
#include <sys/time.h>
#else
#include <time.h>
#include <sys/timeb.h>
#endif

#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>

#include <blitz/array.h>
using namespace blitz;

// Include the rest only if not done yet
#ifndef ANDI_INCLUDED
#define ANDI_INCLUDED

// Global const parameters
//------------------------
const double LA_PI = 3.1415926535897932384626433832795029;
const double LA_PI2 = 2*LA_PI;
const double pi = LA_PI;
const double pi2 = 2*pi;
const double rTOd = 180.0/pi;	// rad to deg

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

// ------------ STRING --------------------------------------------------------------------------------------------------
// Casted alles auf string
template<typename T>
string STRING(const T& t) {
  stringstream s;
  s << t;
  return s.str();
}

// -------------checkparfilename --------------------------------------------------------------------------------------
LA_STRING checkparfilename(LA_STRING name)
{
int i;
i=name.length();
if(name.right(4)==".dat") name=name(1,i-4);
i=name.length();
if(name.left(1)=="_" || name.left(1)=="#") name=name(2,i);
return name;
}

// -------------checkparfilename (string) ------------------------------------------------------------------------------
string checkparfilename(string name)
{
int i;
i=int(name.length());
if(name.substr(i-4,4)==".dat") name=name.substr(0,i-4);
i=int(name.length());
if(name.substr(0,1)=="_" || name.substr(0,1)=="#") name=name.substr(1,i-1);
return name;
}

// ------------- readparfile ---------------------------------------------------------------------------------------------
void readparfile(char* name, vector<double>& vec)
{
// Variablen
int i;
int count=0,chk=0;
double wert;
string parname;
LA_STRING line;

// Einlesen
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; exit(0);}

// Zählen, wieviele Zeilen mit '#' beginnen, diese werden übersprungen
while(1)
{
	in >> line;
	if(line[1]!='#' && chk==0) {cout << "Parameterdatei " << name << " ist nicht im 'Wingen'-Format, bitte 'parfilepipe.exe' ausfuehren." << endl; exit(0);}
	if(line[1]=='#') {chk=1; count+=1; continue;}
	else break;
}
in.close();	// wichtig, damit Einlesen vom Anfang der Datei wieder neu beginnt

in.open(name);	// Erneutes Öffnen der Datei
for(i=1;i<=count;i++)	// Überspringen der IO-Daten Zeilen
{
	in >> line;
}

// Einlesen der Daten 
while(in.eof()==0) // Die letzte Zeile wird immer doppelt eingelesen --- leider nicht zu ändern ---
{
	in >> parname;
	in >> wert;
	vec.push_back(wert);
}
// Löschen der letzten Einträge, um doppeltes Einlesen zu kompensieren
vec.erase(vec.end()-1);

in.close();
}

// ------------- readfile (blitz) ---------------------------------------------------------------------------------------------
// N specifies numger of columns to read
// All indices in vec start with 1 
// First index in vec specifies the row, second one the column, e.g vec(Range::all(),3) are all data in the third column
// This ordering does NOT corresponds to the ordering in memory but is more intuitive
void readfile(char* name, int N, Array<double,2>& vec)
{
// Variables
int i;
int count=0;
LA_STRING line;

// resize vec (Spalte,Zeile)
vec.resize(Range(1,100),Range(1,N));

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; exit(0);}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}
in.close();	// Important to start reading from the beginning of the file

in.open(name);	// Open file again
for(i=1;i<=count;i++)	// Skip IO data rows
{
	in >> line;
}

// Read data
count=0;
while(in.eof()==0) // Last row is read twice --- can't be changed --- -> count-1 is actual number of rows in file
{
	count+=1;
	for(i=1;i<=N;i++) in >> vec(count,i);
	if(count%100==0) vec.resizeAndPreserve(count+100,N);
}

// Resize to actual size (remove dublicate last row as well)
count-=1;
vec.resizeAndPreserve(count,N);

in.close();
}

// ------------- readfile 2 Spalten ---------------------------------------------------------------------------------------------
void readfile(char* name, vector<double>& thetavec, vector<double>& psivec)
{
// Variablen
int i;
int count=0;
double psi,theta;
LA_STRING line;

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; exit(0);}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}
in.close();	// Important to start reading from the beginning of the file

in.open(name);	// Open file again
for(i=1;i<=count;i++)	// Skip IO data rows
{
	in >> line;
}

// Einlesen der Daten (hier für zwei Spalten)
while(in.eof()==0) // Die letzte Zeile wird immer doppelt eingelesen --- leider nicht zu ändern ---
{
	in >> theta;
	in >> psi;
	thetavec.push_back(theta);
	psivec.push_back(psi);
}
// Löschen der letzten Einträge, um doppeltes Einlesen zu kompensieren
thetavec.erase(thetavec.end()-1);
psivec.erase(psivec.end()-1);

in.close();
}

// ------------- readfile 3 Spalten ---------------------------------------------------------------------------------------------
void readfile(char* name, vector<double>& vec1, vector<double>& vec2, vector<double>& vec3)
{
// Variablen
int i;
int count=0;
double A,B,C;
LA_STRING line;

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; exit(0);}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}
in.close();	// Important to start reading from the beginning of the file

in.open(name);	// Open file again
for(i=1;i<=count;i++)	// Skip IO data rows
{
	in >> line;
}

// Einlesen der Daten (hier für drei Spalten)
while(in.eof()==0) // Die letzte Zeile wird immer doppelt eingelesen --- leider nicht zu ändern ---
{
	in >> A;
	in >> B;
	in >> C;
	vec1.push_back(A);
	vec2.push_back(B);
	vec3.push_back(C);
}
// Löschen der letzten Einträge, um doppeltes Einlesen zu kompensieren
vec1.erase(vec1.end()-1);
vec2.erase(vec2.end()-1);
vec3.erase(vec3.end()-1);

in.close();
}

// ------------- readfile 4 Spalten ---------------------------------------------------------------------------------------------
void readfile(char* name, vector<double>& vec1, vector<double>& vec2, vector<double>& vec3, vector<double>& vec4)
{
// Variablen
int i;
int count=0;
double A,B,C,D;
LA_STRING line;

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; exit(0);}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}
in.close();	// Important to start reading from the beginning of the file

in.open(name);	// Open file again
for(i=1;i<=count;i++)	// Skip IO data rows
{
	in >> line;
}

// Einlesen der Daten (hier für vier Spalten)
while(in.eof()==0) // Die letzte Zeile wird immer doppelt eingelesen --- leider nicht zu ändern ---
{
	in >> A;
	in >> B;
	in >> C;
	in >> D;
	vec1.push_back(A);
	vec2.push_back(B);
	vec3.push_back(C);
	vec4.push_back(D);
}
// Löschen der letzten Einträge, um doppeltes Einlesen zu kompensieren
vec1.erase(vec1.end()-1);
vec2.erase(vec2.end()-1);
vec3.erase(vec3.end()-1);
vec4.erase(vec4.end()-1);

in.close();
}

// ------------- readfile 5 Spalten ---------------------------------------------------------------------------------------------
void readfile(char* name, vector<double>& vec1, vector<double>& vec2, vector<double>& vec3, vector<double>& vec4, vector<double>& vec5)
{
// Variablen
int i;
int count=0;
double A,B,C,D,E;
LA_STRING line;

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; exit(0);}

// Count the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#') {count+=1; continue;}
	else break;
}
in.close();	// Important to start reading from the beginning of the file

in.open(name);	// Open file again
for(i=1;i<=count;i++)	// Skip IO data rows
{
	in >> line;
}

// Einlesen der Daten (hier für fünf Spalten)
while(in.eof()==0) // Die letzte Zeile wird immer doppelt eingelesen --- leider nicht zu ändern ---
{
	in >> A;
	in >> B;
	in >> C;
	in >> D;
	in >> E;
	vec1.push_back(A);
	vec2.push_back(B);
	vec3.push_back(C);
	vec4.push_back(D);
	vec5.push_back(E);
}
// Löschen der letzten Einträge, um doppeltes Einlesen zu kompensieren
vec1.erase(vec1.end()-1);
vec2.erase(vec2.end()-1);
vec3.erase(vec3.end()-1);
vec4.erase(vec4.end()-1);
vec5.erase(vec5.end()-1);

in.close();
}

// ------------- countdown 'maxsec' sekunden ---------------------------------------------------------------------------------
void countdown(int maxsec, bool reply = false);	//Prototype with default value
void countdown(int maxsec, bool reply)
{
int sec = 1;
time_t begin,end;

time(&begin);

while(1)
{
        if(begin + sec == time(&end)) {if(reply == true && sec <= maxsec) cout << sec << "\t" << flush; sec += 1;}
        if(sec == maxsec + 2) {cout << "Overwriting File!" << endl; break;}
}
cout << endl;
}

// ------------- outputtest ------------------------------------------------------------------------------------------------
// checks if file already exists. 
// Linux: If yes, main program waits 15 seconds before file is overwritten
// Windows: If yes, asks for permission to overwrite
#ifdef linux
void outputtest(char* name)
{
int chk,sec;
sec = 15;
ifstream in(name);
chk = in.fail();
if(chk==0)
{
        cout << "File " << name << " already exists. Overwrite in " << sec << " seconds" << endl;
        countdown(sec);
}
else return;
}

#else
void outputtest(char* name)
{
int chk;
string frage;
ifstream in(name);
chk = in.fail();
if(chk == 0)
{
	cout << "File " << name << " already exists. Overwrite? (y/n): "; cin >> frage;
	if(frage[0] == 'y') return;
	else exit(0);
}
else return;
}
#endif

// ------------- outputtest mit Rückgabewert -----------------------------------------------------------------------------------
// see above
//int outputtest2(char* name)	// 0: File kann geschrieben werden	1: File existiert und soll nicht überschrieben werden! 
//{
//int chk;
//string frage;
//ifstream in(name);
//chk=in.fail();
//if(chk==0)
//{
//	cout << "File " << name << " already exists. Overwrite? (y/n): "; cin >> frage;
//	if(frage[0]=='j') return 0;
//	else return 1;
//}
//else return 0;
//}

//----------------------- count_rows ---------------------------------------------------------------------------------------
// counts number of rows in file
int count_rows(char* name)
{
int count = 0;
int row = 1;
LA_STRING line;

// read file
ifstream in;
in.open(name);
if(in.fail() == 1) {cout << "Cannot open file " << name << endl; exit(0);}

// Count rows, starting with '#'
while(1)
{
	in >> line;
	if(line[1] == '#') {count += 1; continue;}
	else break;
}

while(in.eof()==0) // Last row is read twice --- can't be changed --- -> row-1 is actual number of rows in file
{
	in >> line;
	row += 1;
}
row -= 1;

in.close();
return row;
}

//----------------------- count_column -----------------------------------------------------------------------------------------
// counts number of columns in file, separated by taps
int count_column(char* name)
{
int i;
int count = 0;
int column = 1;
LA_STRING line;

// read file
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Cannot open file " << name << endl; exit(0);}

// Count rows, starting with '#'
while(1)
{
	in >> line;
	if(line[1] == '#') {count += 1; continue;}
	else break;
}
for(i=1;i<=int(line.length());i++)
{
	if(line[i] == '\t') column += 1;
	//if(a[i]==cout.widen('\n')) break; // doesn't work!???!
}

in.close();	// nessecary to start reading at top of file again
return column;
}

//-------------- modulo2pi ------------------------------------------------------------------------------------------------------
inline double modulo2pi(double x)
{
x = fmod(x,double(LA_PI2));
if(x<0) x += LA_PI2;
return x;
}

//----------- sign -------------------------------------------------------------------------------------------------------------
inline int sign(double x)
{
if(x == 0) return 0;
if(x > 0) return 1;
else return -1;
}

//------------ zeit ------------------------------------------------------------------------------------------------------------
// returns system time in seconds, including milliseconds
#ifdef linux
double zeit(void)
{
double aus;
struct timeval tstruct;

while(gettimeofday(&tstruct,0) == -1);

aus=double( double(tstruct.tv_sec) + int(double(tstruct.tv_usec)*1.e-3+0.5)*1.e-3 );
return aus;
}

#else
double zeit(void)
{
struct _timeb timebuffer;
double aus;

_ftime64_s( &timebuffer );

aus=timebuffer.time+1e-3*timebuffer.millitm;
return aus;
}
#endif

//--------------------- polar_r -------------------------------------
// transforms carthesian (x,y) to polar coordinates (r,phi)
// r coordinate is returned
inline double polar_r(double x, double y)
{
return sqrt(x*x + y*y);
}

//-------------------- polar_phi ------------------------------------
// transforms carthesian (x,y) to polar coordinates (r,phi)
// phi coordinate is returned
inline double polar_phi(double x, double y)
{
double phi = atan(y/x);
if(x<0) phi += LA_PI;
if(x>0 && y<0) phi += LA_PI2;
if(x==0) 
{
	if(y>=0) phi = 0.5*LA_PI;
	else phi = 1.5*LA_PI;
}
return phi;
}

//-------------- ran0 ----------------------------------------------
// Minimal random number generator of Park and Miller. Returns a uniform random deviate
// between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value MASK)
// to initialize the sequence; idum must not be altered between calls for successive deviates in
// a sequence.
double ran0(long& idum)
{
const int IA=16807;
const int IM=2147483647;
const double AM=(1.0/IM);
const int IQ=127773;
const int IR=2836;
const int MASK=123459876;

long k;
double ans;

idum ^= MASK;	//XORing with MASK allows use of zero and other simple bit patterns for idum.
k=idum/IQ;
idum=IA*(idum-k*IQ)-IR*k; //Compute idum=(IA*idum) % IM without overflows by Schrage’s method.
if(idum < 0) idum += IM; 
ans=AM*idum;		//Convert idum to a floating result.
idum ^= MASK;		//Unmask before return.
return ans;
}


#endif // ANDI_INCLUDED



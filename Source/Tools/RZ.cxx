// Program transforms from (theta,r) to (R,Z) coordinates
// First two columns are replaced accordingly, all other columns (up to total number of 5) are copied unchanged
// A. Wingen		6.3.08

// Input:	1. File to transform	2. total number of columns (optional, Warning: automated column count can be wrong) 
// Output:	to (R,Z) coordinates transformed file; it gets the additional praefix _RZ

// Define
//--------
//#define BZ_DEBUG
#define program_name "RZ"

// Include
//--------
#include <andi.hxx>
//#include <d3d.hxx>

// Prototypes  

// Switches

// Golbal Parameters 
const double pi=LA_PI;
const double pi2=2*pi;


int main(int argc, char *argv[])
{
// Variables
int i,j;
int N = 0;
double R,Z;
double R0,Z0;	//Magnetic Axis coordinates
Array<double,2> data;

LA_STRING basename,startname;
if(argc==3) N = atoi(argv[2]);
if(argc>=2) 
{
	startname = LA_STRING(argv[1]);
}
else
{
	cout << "No input files. -> Abort" << endl; 
	exit(0);
}
startname = checkparfilename(startname);
LA_STRING name = startname + ".dat";

// Determine Number of Columns
if(N==0) N = count_column(name);
cout << N << " Columns found" << endl;

// Read data
readfile(name,N,data);
int M = data.rows();

// Output
LA_STRING filenameout = startname + "_RZ" + ".dat";
ofstream out(filenameout);
out.precision(16);

// Read and write IO data; Find and read magnetic Axis coordinates from IO header
int count = 1;
LA_STRING line;
ifstream in(name);
in >> line;	//first line
out << "# Transformed to (R,Z) coordinates. Original file " << name << endl;
while(1)
{
	in >> line;
	//count += 1;
	if(line.length()>18)
	{
		if(line[3] == 'M' && line[4] == 'a' && line[5] == 'g' && line[6] == 'n' && line[7] == 'e' && line[8] == 't' && line[9] == 'i' && line[10] == 'c' && count < 3)
		{
			if(line[18]=='R') R0 = atof(line.mid(22));
			else Z0 = atof(line.mid(22));
			count +=1;
		}
	}
	if(line[1]=='#') out << line << endl;
	else break;
}
if(count!=3) {cout << "R0 und Z0 nicht gefunden -> Abbruch" << endl; exit(0);}
cout << "R0 = " << R0 << "\t" << "Z0 = " << Z0 << endl;

out << "# R[m]" << "\t" << "Z[m]" << "\t" << "see above" << endl;
out << "#" << endl;

// convert Data
for(i=1;i<=M;i++)
{
	R = data(i,2)*cos(data(i,1)) + R0;
	Z = data(i,2)*sin(data(i,1)) + Z0;

	out << R << "\t" << Z << "\t";
	for(j=3;j<N;j++) out << data(i,j) << "\t";
	out << data(i,N) << endl;
}

return 0;
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



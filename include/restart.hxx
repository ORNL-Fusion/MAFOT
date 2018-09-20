// Header-File for restarting MAFOT subroutines
// uses arrays and multiple-arrays from blitz-Library
// A.Wingen						19.9.18

// Define
//--------
#ifndef RESTART_INCLUDED
#define RESTART_INCLUDED

// Include
//--------
#include <la_string.hxx>
#include <fstream>
#include <sstream>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <andi.hxx>
using namespace blitz;

// --------------- Prototypes ---------------------------------------------------------------------------------------------
void readMafotFile(char* name, vector<LA_STRING>& head, Array<double,2>& data);
int restartMafot(ofstream& out, int N, int NoOfPackages, vector<LA_STRING>& head, Array<double,2>& data, bool writeHeader = false);
bool checkOutputFile(char* name);

// -------------- global Parameters ---------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//---------------- checkOutputFile --------------------------------------------------------------------------------------
// checks if file exists and returns true or false
bool checkOutputFile(char* name)
{
int chk;
ifstream in(name);
chk = in.fail();
if(chk==0)
{
		cout << "File " << name << " already exists." << endl;
		return true;
}
else return false;
}

//---------------- restartMafot --------------------------------------------------------------------------------------
// write the data into the file. If writeHeader is true (default is false), also the header is written.
int restartMafot(ofstream& out, int N, int NoOfPackages, vector<LA_STRING>& head, Array<double,2>& data, bool writeHeader)
{
int i,k;
int packageSize = int(N/NoOfPackages);
int donePackages = int((data.rows()-1)/packageSize);	// last row can be garbage

if(data.rows() < N)
{
	// Output
	if(writeHeader)
	{
		for(i=0;i<head.size();i++) out << head[i] << endl;	// write header
	}
	for(k=1;k<=donePackages*packageSize;k++)
	{
		out << data(k,1);
		for(i=2;i<=data.cols();i++) out << "\t" << data(k,i);
		out << endl;
	}
	return donePackages;
}
else return 0;
}

//---------------- readMafotFile --------------------------------------------------------------------------------------
// head contains the header lines as LA_STRING, head[0,1,...,N-1]
// data contains the actual data in the file
void readMafotFile(char* name, vector<LA_STRING>& head, Array<double,2>& data)
{
// Variables
int i,k;
LA_STRING line;
int count = 0;
int rows = 0;
int cols = 1;
vector<string> words;

// Input
ifstream in;
in.open(name);
if(in.fail()==1) {cout << "Unable to open file " << name << endl; EXIT;}

// Read the number of rows starting with #
while(1)
{
	in >> line;
	if(line[1]=='#')
	{
		count+=1;
		head.push_back(line);
	}
	else break;
}

words = split(string(line));
cols = words.size();

// count the number of rows with data
while(in.eof()==0) // Last row is read twice --- can't be changed --- -> rows starts with 0 and is actual number of data rows in file
{
	in >> line;
	rows += 1;
}

in.close();	// Important to start reading from the beginning of the file
in.clear(); // Important to clear EOF flag

// resize data (Spalte,Zeile)
data.resize(Range(1,rows),Range(1,cols));

in.open(name);	// Open file again
for(i=1;i<=count;i++)	// Skip IO data rows
{
	in >> line;
}

string value;
// Read data
for(k=1;k<=rows;k++)
	for(i=1;i<=cols;i++)
	{
		in >> value;						// this can read also nan and inf; direct read into data(k,i) fails for nan or inf!
		data(k,i) = atof(value.c_str());	// returns double, 0 for non numeric values, nan for nan and inf for inf
	}

in.close();
return;
}

#endif // RESTART_INCLUDED
//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

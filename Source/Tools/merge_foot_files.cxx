// Program merges footprint files that have been manually divided for parallel computation
// Files have to have a continuing number at the end of the filename (before the .dat)
// Input filename has to be the first file to be read. All following files are read automatically
// Number of the first file has to be a single digit 
// A. Wingen		15.1.09

// Input:	1. first file to be merged	2. total number of columns (optional, Warning: automated column count can be wrong) 
// Output:	File which contains all data from all files written among one another

// Define
//--------
//#define BZ_DEBUG
#define program_name "merge_foot_files"

// Include
//--------
#include <andi.hxx>
//#include <efit_class.hxx>
//#include <d3d-drift.hxx>

// Prototypes  

// Switches

// Golbal Parameters 

// Main Program
//--------------
int main(int argc, char *argv[])
{
// Variables
int i,j,k;
int count,jstart;
int N = 0;
LA_STRING line;
LA_STRING filename;
ifstream in;

Array<double,2> data;

// Input file names
LA_STRING basename;
if(argc==3) N = atoi(argv[2]);
if(argc>=2) basename = LA_STRING(argv[1]);
else	// No Input: Abort
{
	cout << "No Input files -> Abort!" << endl;
	exit(0);
}
if(basename.right(4) != ".dat") {cout << "Please enter full file name" << endl; exit(0);}

// count columns if necessary
if(N==0) N = count_column(basename);
cout << N << " Columns found" << endl;

// Output
int FileNr = 1;
i = basename.length();
basename = basename(1,i-4);
if(atoi(basename.right(2)) > 0) {FileNr = atoi(basename.right(2)); i -= 6;}	// non-numerical string causes atoi to return 0
else {FileNr = atoi(basename.right(1)); i -= 5;}
basename = basename(1,i);
LA_STRING name = basename + ".dat";
ofstream out(name);
out.precision(16);

i = FileNr;
while(1)
{
	filename = basename + LA_STRING(i) + ".dat";

	// Check if file exists
	in.open(filename);
	if(in.fail()==1) {break;}
	cout << filename << endl;

	// Read and write file header
	if(i==FileNr)	
	{
		while(1) 
		{
			in >> line;
			if(line[1]=='#') out << line << endl;
			else break;
		}
		jstart = 1;
	}
	in.close();

	// Read data
	readfile(filename,N,data);		

	// Find Number of rows for first t value
	if(i>FileNr)
	{
		count = 1;
		for(j=2;j<=data.rows();j++)
		{
			if(data(j,2)==data(1,2)) count += 1;
			else break;
		}
		jstart = count + 1;
	}

	// Write data of other files to output
	for(j=jstart;j<=data.rows();j++)
	{
		for(k=1;k<N;k++) out << data(j,k) << "\t";
		out << data(j,N) << endl;
	}

	// next file
	i += 1;
}


return 0;
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



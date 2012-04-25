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
#define program_name "merge_foot_files_x_sliced"

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
int count,jstart,blocks;
int N = 0;
int M = 0;
LA_STRING line;
LA_STRING filename;
ifstream in;
Range all = Range::all();

Array<double,2> data;
Array<int,1> groesse(Range(1,1));

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
M = count_rows(basename);
cout << M << " Rows found" << endl;
if(N==0) N = count_column(basename);
cout << N << " Columns found" << endl;

Array<double,3> all_data(Range(1,1),Range(1,M),Range(1,N));

// Output
i = basename.length();
basename = basename(1,i-5);
LA_STRING name = basename + ".dat";
ofstream out(name);
out.precision(16);

i = 1;
while(1)
{
	filename = basename + LA_STRING(i) + ".dat";

	// Check if file exists
	in.open(filename);
	if(in.fail()==1) {break;}
	cout << filename << endl;

	// Read and write file header
	if(i==1)	
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

	// Store data
	if(data.rows() > M) M = data.rows();
	all_data.resizeAndPreserve(i,M,N);
	groesse.resizeAndPreserve(i);
	for(j=1;j<=data.rows();j++) all_data(i,j,all) = data(j,all);

	// Get number of rows with same 2. column value for each file -> resolution in y direction
	count = 1;
	for(j=2;j<=data.rows();j++) 
	{
		if(data(j,2)==data(1,2)) count += 1;	
		else break;
	}
	groesse(i) = count;
	blocks = data.rows()/count;	// The same for all files

	// next file
	i += 1;
}

// Write data of all files to output
for(int l=1;l<=blocks;l++)
{
	for(i=1;i<=all_data.rows();i++)
	{
		count = groesse(i);
		jstart = (l-1)*count + 1;
		for(j=jstart;j<jstart+count;j++)
		{
			if(i>1 && j==jstart) continue;
			for(k=1;k<N;k++) out << all_data(i,j,k) << "\t";
			out << all_data(i,j,N) << endl;
		}
	}
}
return 0;
} //end of main
 

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



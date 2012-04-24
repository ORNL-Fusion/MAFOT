// Program calculates targets, wall and boundary box, to plot.
// Necessary input: Position of magnetic axis
// A. Wingen		6.3.08

// Input:	1. Parameterfile (optional)
// Output:	target0.dat		-> vertical shelf above 45° target
//			target1.dat		-> 45° target
//			target2.dat		-> horizontal target
//			target3.dat		-> horizontal shelf above pump to horizontal target
//			wall.dat		-> entire wall 
//			boundary.dat	-> boundary box
// All outputfiles include theta, r, R, Z coordinates

// Define
//--------
//#define BZ_DEBUG
#define program_name "target"

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
double r,theta,R,Z;
double R0,Z0;
int N=101;
LA_STRING name,basename,parfilename;
vector<double> vec;

// Input: Magnetic Axis
if(argc>=2) 
{
	basename = LA_STRING(argv[1]);
	basename=checkparfilename(basename);
	parfilename = "_" + basename + ".dat";
	readparfile(parfilename,vec);
	R0=vec[9];
	Z0=vec[10];
}
else	// no input: ask
{
	cout << "Magnetic Axis Position?" << endl;
	cout << "R0 = "; cin >> R0;
	cout << "Z0 = "; cin >> Z0;
}
cout << "Magnetic Axis Position: R0 = " << R0 << "\t" << "Z0 = " << Z0 << endl;

//-----------------------------------------------------------------------
// Targets
double R1,Z1;	// upper or left Point
double R2,Z2;	// lower or right Point

ofstream out2;

for(j=0;j<4;j++)	// Target (0=vertical wall, 1=45°, 2=horizontal, 3=shelf)
{
	switch(j)
	{
	case 0:
		R1=1.0161-R0;	Z1=-1.22884-Z0;
		R2=1.0161-R0;	Z2=-1.034873-Z0;
		name="target0.dat";
		break;
	case 1:
		R1=1.0161-R0;	Z1=-1.22884-Z0;
		R2=1.15285-R0;	Z2=-1.3664-Z0;
		name="target1.dat";
		break;
	case 2:
		R1=1.15285-R0;	Z1=-1.3664-Z0;
		R2=1.372-R0;	Z2=-1.3664-Z0;	// R2=1.419845-R0;
		name="target2.dat";
		break;
	case 3:
		R1=1.372-R0;	Z1=-1.25-Z0;
		R2=1.59115-R0;	Z2=-1.25-Z0;	
		name="target3.dat";
		break;
	default:
		cout << "No target specified" << endl;
		exit(0);
		break;
	}
	out2.open(name);

	for(i=0;i<N;i++)
	{
		Z=Z1+i*(Z2-Z1)/double(N-1);
		R=R1+i*(R2-R1)/double(N-1);

		r=sqrt(R*R+Z*Z);
		theta=atan(Z/R);
		if(R<0) theta+=pi;
		if(R>=0 && Z<0) theta+=pi2;

		out2 << theta << "\t" << r << "\t" << R+R0 << "\t" << Z+Z0 << endl;
	}
	out2.close();
}

//-----------------------------------------------------------------------
// Wall
name="new_vessel.dat";
vector<double> Rvec,Zvec;
readfile(name,Rvec,Zvec);
N=int(Rvec.size());

ofstream out("wall.dat");

for(i=0;i<N;i++)
{
	R=Rvec[i]/1000.0-R0;
	Z=Zvec[i]/1000.0-Z0;

	r=sqrt(R*R+Z*Z);
	theta=atan(Z/R);
	if(R<0) theta+=pi;
	if(R>=0 && Z<0) theta+=pi2;

	out << theta << "\t" << r << "\t" << R+R0 << "\t" << Z+Z0 << endl;
}

//-----------------------------------------------------------------------
// Boundary Box
const double Rmin=1.016-R0;
const double Rmax=2.4-R0;
const double Zmin=-1.367-Z0;
const double Zmax=1.36-Z0;

ofstream out3("boundary.dat");

for(i=0;i<=50;i++)
{
	Z=i*(Zmax)/50.0;
	R=Rmax;

	r=sqrt(R*R+Z*Z);
	theta=atan(Z/R);
	if(R<0) theta+=pi;
	if(R>=0 && Z<0) theta+=pi2;

	out3 << theta << "\t" << r << "\t" << R+R0 << "\t" << Z+Z0 << endl;
}

for(i=0;i<=100;i++)
{
	Z=Zmax;
	R=Rmax-i*(Rmax-Rmin)/100.0;

	r=sqrt(R*R+Z*Z);
	theta=atan(Z/R);
	if(R<0) theta+=pi;
	if(R>=0 && Z<0) theta+=pi2;

	out3 << theta << "\t" << r << "\t" << R+R0 << "\t" << Z+Z0 << endl;
}

for(i=0;i<=100;i++)
{
	Z=Zmax-i*(Zmax-Zmin)/100.0;
	R=Rmin;

	r=sqrt(R*R+Z*Z);
	theta=atan(Z/R);
	if(R<0) theta+=pi;
	if(R>=0 && Z<0) theta+=pi2;

	out3 << theta << "\t" << r << "\t" << R+R0 << "\t" << Z+Z0 << endl;
}

for(i=0;i<=100;i++)
{
	Z=Zmin;
	R=Rmin+i*(Rmax-Rmin)/100.0;

	r=sqrt(R*R+Z*Z);
	theta=atan(Z/R);
	if(R<0) theta+=pi;
	if(R>=0 && Z<0) theta+=pi2;

	out3 << theta << "\t" << r << "\t" << R+R0 << "\t" << Z+Z0 << endl;
}

for(i=0;i<50;i++)
{
	Z=Zmin+i*(-Zmin)/50.0;
	R=Rmax;

	r=sqrt(R*R+Z*Z);
	theta=atan(Z/R);
	if(R<0) theta+=pi;
	if(R>=0 && Z<0) theta+=pi2;

	out3 << theta << "\t" << r << "\t" << R+R0 << "\t" << Z+Z0 << endl;
}

return 0; 
} //end of main

//------------------------ End of Main ------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------


//----------------------- End of File -------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------



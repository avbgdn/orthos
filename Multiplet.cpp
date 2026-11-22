// File purpose: isotropic hyperfine multiplet from groups of equivalent nuclei: the line positions and relative
// intensities (stick spectrum) that are later used to broaden the spectrum by isotropic hyperfine splitting.
// The class Multiplet itself is declared in BiradicalODF.h (member Mult of MagneticParameters).
// Used by: Discrepancy.cpp (Discrepancy() calls SetFromHFS()); Multiplet::Multiplet(), Multiplet::~Multiplet()
// and CalculateMultiplet() are called implicitly or from this file only.
#include "BiradicalODF.h"


// CalculateMultiplet: number of equivalent combinations (line degeneracy) for each total spin projection of N
// equivalent nuclei, each with Sx2+1 projections (Sx2 is twice the nuclear spin); this is repeated convolution
// of the single-nucleus pattern (Sx2+1 ones) with itself N times.
// mult: output array of Sx2*N+1 integers allocated by the caller; mult[j] is the degeneracy of the j-th line.
// Used by: Multiplet::SetFromHFS() in this file.
void CalculateMultiplet(int Sx2, int N, int *mult)
{
	int *mult2;
	// scratch array for the new multiplet while the old one is still needed
	mult2=(int*)malloc((Sx2*N+1)*sizeof(int));

	// start from a single nucleus-free line of unit weight
	mult[0]=1;
	for(int i=1; i<Sx2*N+1; i++) mult[i]=0;
	// add the equivalent nuclei one by one, each step convolves the multiplet with a flat pattern of Sx2+1 lines
	for(int i=0; i<N; i++)       // for each of the equivalent nuclei:
	{
		for(int j=0; j<Sx2*N+1; j++)		// for each line in the multiplet
		{
			int sum=0;
			for(int k=j-Sx2; k<=j; k++)
				if((k>=0)&&(k<Sx2*N+1)) sum+=mult[k];
			mult2[j]=sum;
		}
		for(int j=0; j<Sx2*N+1; j++)		// refreshing the multiplet array
		{
			mult[j]=mult2[j];
		}
	}

	return;

}

// Commented-out old copy of the class declaration; the live declaration (with a three-column Transitions array) is in BiradicalODF.h
/*
class Multiplet
{
public:
	Array<double, Dynamic, 2> Transitions;
	Array<int, Dynamic, 1> PrincipalSpinProjection;
	int NLines;

	Multiplet();
	~Multiplet();

	void SetFromHFS(double Field0, Hyperfine * HFS);

};
*/

// Multiplet constructor: does nothing, the arrays are sized in SetFromHFS().
// Used by: called implicitly whenever a Multiplet object is created, i.e. for each MagneticParameters object,
// which has a Multiplet member (Mult, declared in BiradicalODF.h).
Multiplet::Multiplet()
{
}

// Multiplet destructor: does nothing.
// Used by: called implicitly whenever a Multiplet object is destroyed (as a member of MagneticParameters).
Multiplet::~Multiplet()
{
}

// SetFromHFS: builds the isotropic hyperfine multiplet of one radical from all of its groups of equivalent nuclei.
// Transitions has NLines rows: column 0 is the line position (field offset, centred on Field0), column 1 the relative
// intensity (the intensities of all lines add up to 1), column 2 the nuclear spin projection of the line, which is
// filled only when the first hyperfine group is a single nucleus (zero otherwise).
// Field0: centre of the multiplet; HFS: array of hyperfine groups (degeneracy, spin, coupling constant);
// NNuclei: number of groups in HFS.
// Used by: Discrepancy() in Discrepancy.cpp (called with Field0=0 for every particle).
void Multiplet::SetFromHFS(double Field0, Hyperfine * HFS, int NNuclei)
{
	// total number of lines is the product of the line counts of all groups
	NLines=1;
	for(int i=0; i<NNuclei; i++)
		NLines*=(HFS[i].Degeneracy*HFS[i].Sx2+1);
//	cout << HFS[0].Degeneracy << " " << HFS[0].Sx2 << " " << NLines << endl;
	// start from a single line at Field0 with unit intensity, all other rows zero
	Transitions.resize(NLines, 3);

	Transitions(0, 0)=Field0;
	Transitions(0, 1)=1.0;
	Transitions(0, 2)=0;
	for(int i=1; i<NLines; i++)
	{
		Transitions(i, 0)=0;
		Transitions(i, 1)=0;
		Transitions(i, 2)=0;
	}

	// add the groups one by one: each existing line is replaced by a whole multiplet of the new group, expanded in place from the end of the array
	int *mult;
	int Completed=1;
	for(int i=NNuclei-1; i>=0; i--)
	{
		mult=(int*)malloc((HFS[i].Sx2*HFS[i].Degeneracy+1)*sizeof(int));
		CalculateMultiplet(HFS[i].Sx2, HFS[i].Degeneracy, mult);
// for each line in the multiplet
		for(int j=HFS[i].Sx2*HFS[i].Degeneracy; j>=0; j--)
		{
			for(int k=Completed-1; k>=0; k--)
			{
				Transitions(j*Completed+k, 1)=Transitions(k, 1)*(double)mult[j]/pow((double)(HFS[i].Sx2+1), HFS[i].Degeneracy);
				Transitions(j*Completed+k, 0)=Transitions(k, 0)-HFS[i].Degeneracy*HFS[i].NuclearSpin*HFS[i].HFConstant+j*HFS[i].HFConstant; 
				if((!i)&&(HFS[i].Degeneracy==1))
					Transitions(j*Completed+k,2)=-HFS[i].NuclearSpin+j;
			}
		}
//		cout << Transitions;
//		cout << endl;
		Completed*=HFS[i].Sx2*HFS[i].Degeneracy+1;

		free(mult);
	}



}

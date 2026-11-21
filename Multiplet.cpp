#include "BiradicalODF.h"


void CalculateMultiplet(int Sx2, int N, int *mult)
{
	int *mult2;
	mult2=(int*)malloc((Sx2*N+1)*sizeof(int));

	mult[0]=1;
	for(int i=1; i<Sx2*N+1; i++) mult[i]=0;
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

Multiplet::Multiplet()
{
}

Multiplet::~Multiplet()
{
}

void Multiplet::SetFromHFS(double Field0, Hyperfine * HFS, int NNuclei)
{
	NLines=1;
	for(int i=0; i<NNuclei; i++)
		NLines*=(HFS[i].Degeneracy*HFS[i].Sx2+1);
//	cout << HFS[0].Degeneracy << " " << HFS[0].Sx2 << " " << NLines << endl;
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

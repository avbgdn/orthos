/*
    This file is a part of Orthos software for Biradical EPR Spectra simulation
    Copyright (C) 2015-2016  Alexey Bogdanov

    The details and updates of this program can be found on its website:
    http://orthos-epr.sourceforge.net/

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "LiouvilleMatrixElementB.h"

// integration net for potential order parameters
// #define INT_NET 4000 
// #define INT_NET_FI 8000

// CLASS RELATING TO MEAN FIELD POTENTIAL MANIPULATION:
// - MF Potential decomposition coefficients storage and variation
// - order parameter calculation and storage
// - XLK coeffients of rotation diffusion superoperator calculation and storage
// - printout

double Potential::Get(int L, int K)
	{
		if (L%2) return 0;
		if (K%2) return 0;
		if (L>MAX_POT) return 0;
		if(K<0) K=-K;
		if (K>L) return 0;
		return A[L*(L+2)/8+K/2];
	}
int Potential::GetVF(int L, int K)
	{
		if (L%2) return 0;
		if (K%2) return 0;
		if (L>MAX_POT) return 0;
		if(K<0) K=-K;
		if (K>L) return 0;
		return VF[L*(L+2)/8+K/2];
	}

int Potential::GetNArg(int L, int K)
	{
		if (L%2) return 0;
		if (K%2) return 0;
		if (L>MAX_POT) return 0;
		if(K<0) K=-K;
		if (K>L) return 0;
		return NArg[L*(L+2)/8+K/2];
	}

void Potential::Add(int L, int K, double Value, int VariateFlag)
	{
		if (L%2) return;
		if (K%2) return;
		if (L>MAX_POT) return;
		if (L>Lmax) Lmax=L;
		A [L*(L+2)/8+K/2]=Value;
		VF[L*(L+2)/8+K/2]=VariateFlag;
		if ((K!=0)&&(Value!=0)) isAxial=false;
	}
void Potential::SetNArg(int L, int K, int NA)
{
		NArg[L*(L+2)/8+K/2]=NA;
}

double Potential::GetSOP(int L, int K, int M)
{
		if (M)   return 0;
		if (L%2) return 0;
		if (K%2) return 0;
		if (L>2*MAX_TRUNC) return 0;
		if(K<0) K=-K;
		if (K>L) return 0;
		return SOP[L*(L+2)/8+K/2];
}

std::complex<double> Potential::GetOP_JUMP_P(int L, int K, int NAxis)
{
		if (NAxis>=MagneticParameterParent->NJumpAxes) return 0;
		if (L%2) return 0;
		if (L>2*MAX_TRUNC) return 0;
		if (K>L) return 0;
		if(K<0) 
		{
			K=-K;
			if(K%2) return -conj(OP_JUMP_P[NAxis][L*L/4+K]);
			else	return  conj(OP_JUMP_P[NAxis][L*L/4+K]);
		}
		else return OP_JUMP_P[NAxis][L*L/4+K];
}

std::complex<double> Potential::GetOP_JUMP_M(int L, int K, int NAxis)
{
		if (NAxis>=MagneticParameterParent->NJumpAxes) return 0;
		if (L%2) return 0;
		if (L>2*MAX_TRUNC) return 0;
		if (K>L) return 0;
		if(K<0) 
		{
			K=-K;
			if(K%2) return -conj(OP_JUMP_M[NAxis][L*L/4+K]);
			else	return  conj(OP_JUMP_M[NAxis][L*L/4+K]);
		}
		else return OP_JUMP_M[NAxis][L*L/4+K];
}

double Potential::GetOP(int L, int K)
{
	if (L%2)
	{ 
//		cout << "GetOP" <<  L << ", " << K << " = " << 0 << endl;
		return 0;
	}
	if (K%2) 
	{ 
//		cout << "GetOP" <<  L << ", " << K << " = " << 0 << endl;
		return 0;
	}
	if (L>2*MAX_TRUNC) 
	{ 
//		cout << "GetOP" <<  L << ", " << K << " = " << 0 << endl;
		return 0;
	}
	if(K<0) K=-K;
	if (K>L) 
	{ 
//		cout << "GetOP" <<  L << ", " << K << " = " << 0 << endl;
		return 0;
	}
//		cout << "GetOP" <<  L << ", " << K << " = " << OP[L*(L+2)/8+K/2] << endl;
		return OP[L*(L+2)/8+K/2];
}

void Potential::CalculateXLK(int MaxL)
{
//	cout << "Calculating XLK\t";

	if (MaxL>2*MAX_POT)
	{
		printf("Error in Potential::CalculateXLK: MaxL>2*MAX_POT!\n");
		return;
	}
	for(int L=0; L<=MaxL; L+=2)
		for(int K=-L; K<=L; K++)
			for(int M=-L; M<=L; M++)
			{
				if(M)
				{
					XLK_Matrix[L/2][MAX_POT+K][MAX_POT+M]=0;
					continue;
				}
				XLK_Matrix[L/2][MAX_POT+K][MAX_POT+M]=XLK_ND(L, K, *this, (MagneticParameterParent->DiffusionR));

			}
//	cout << "Finished\n";

			return;

}

std::complex<double> Potential::GetXLK(int L, int K, int M)
{
	if(L%2)			return 0;
	if(abs(K)>L)	return 0;
	if(M)			return 0;
	return XLK_Matrix[L/2][MAX_POT+K][MAX_POT+M];
}
void Potential::CalculateOP(int MaxL)
{
//	cout << "Calculating potential order parameters\t";
	for (int i=0; i<(2*MAX_TRUNC+2)*(2*MAX_TRUNC+4)/8; i++) 
	{
		OP[i]=0;
		SOP[i]=0;
	}
	if(Lmax==0)
	{
		OP[0]=1.0;
		SOP[0]=1.0;
		return;
	}

	for(double cosbeta=0; cosbeta<=1; cosbeta+=1/((double)INT_NET))
	{
		double beta=acos(cosbeta);

		int int_net_fi=(int)(INT_NET_FI*sin(beta))+1;		// experimental!
//		int int_net_fi=INT_NET_FI;
	  for(double gamma=0; gamma<=PI/2; gamma+=PI/(2*int_net_fi))
	  {
		double ro=0;
		for (int l=0; l<=Lmax; l+=2)
		{
// calculating ro=potential!
			ro+=this->Get(l, 0)*real(ComplexWigner(l, 0, 0, 0, beta, gamma) );
			if(!isAxial)
			  for(int k=2; k<=l; k+=2)
				ro+=this->Get(l, k)*2.0*real(ComplexWigner(l, 0, k, 0, beta, gamma) );
		}
// now ro becomes boltzman orientation distribution function:
			double ro2=exp(ro);
			ro=exp(0.5*ro);
//			cout << ro << " " << ro2 <<endl;
			for( int L=0; L<=2*MaxL; L+=2 )
			  for( int K=0; K<=L; K+=2 )
				if((!isAxial)||(K==0))
				{
					OP[L*(L+2)/8+K/2] += real(ComplexWigner(L, 0, K, 0, beta, gamma))*ro2/(int_net_fi);	// integration step no *sin(beta) because cosine integration!
					SOP[L*(L+2)/8+K/2]+= real(ComplexWigner(L, 0, K, 0, beta, gamma))*ro /(int_net_fi);	// integration step
				}



	  }
	}
	  
// normalization:
//	  cout << MaxL << endl;
	  Normalization_OP=OP[0];
	  for (int i=1; i<(2*MaxL+2)*(2*MaxL+4)/8; i++) 
	  {
		SOP[i]/=SOP[0];
		OP[i] /=OP[0];
	  }
	  SOP[0]=1.0;
	  OP[0]= 1.0;

//	cout << "Finished!\n";
}


Potential::Potential()
{
		isAxial=true;
		Lmax=0;
		for (int i=0; i<(MAX_POT+2)*(MAX_POT+4)/8; i++) 
		{
			A[i]=0;
			VF[i]=0;
		}
		for (int i=0; i<(MAX_TRUNC+2)*(MAX_TRUNC+4)/8; i++) 
		{
			OP[i]=0;
			SOP[i]=0;
		}
}
void Potential::PrintOut(int MaxL)
{
	FILE *potout;
    char buff[MAX_STRLEN];
    sprintf(buff, "potential%d.out", MyCount);

	potout=fopen(buff, "w");
		for(int L=2; L<=Lmax; L+=2)
			for(int K=0; K<=L; K+=2)
				fprintf(potout, "Potential(%d, %d)=%lf; %d\n", L, K, this->Get(L, K), this->GetVF(L, K));
		for(int L=0; L<=2*MaxL; L+=2)
			for(int K=0; K<=L; K+=2)
				fprintf(potout, "OP(%2d, %2d)=%lf +- %e ,\t SOP(%d, %d)=%lf\n", L, K, this->GetOP(L, K), fabs(OP2[L*(L+2)/8+K/2]-this->GetOP(L, K)), L, K, this->GetSOP(L, K, 0));
		for(int L=2; L<=2*Lmax; L+=2)
			for(int K=-L; K<=L; K+=1)
				fprintf(potout, "XLK(%d, %d)=%lf  %lf\n", L, K, real(this->GetXLK(L, K, 0)), imag(this->GetXLK(L, K, 0)));

for(int i=0; i<MagneticParameterParent->NJumpAxes; i++)
{
	fprintf(potout, "\n\n\nJUMP DIFFUSION for axis %d\n\n", i);
		for(int L=0; L<=2*MaxL; L+=2)
			for(int K=0; K<=L; K+=1)
			{
				fprintf(potout, "OP_JUMP_P(%2d, %2d)=%+lf + i*%+lf +- %e ,\t", L, K, real(OP_JUMP_P[i][L*L/4+K]), imag(OP_JUMP_P[i][L*L/4+K]), abs(OP_JUMP_P[i][L*L/4+K]-OP_JUMP_P2[i][L*L/4+K]));
				fprintf(potout, "OP_JUMP_M(%2d, %2d)=%+lf + i*%+lf +- %e \n", L, K, real(OP_JUMP_M[i][L*L/4+K]), imag(OP_JUMP_M[i][L*L/4+K]), abs(OP_JUMP_M[i][L*L/4+K]-OP_JUMP_M2[i][L*L/4+K]));
			}


}

	fclose(potout);
}

Potential::~Potential()
{}






void Potential::CalculateOP_GK(int MaxL)
{
//	cout << "Calculating potential order parameters\t";
	for (int i=0; i<(MAX_TRUNC+2)*(MAX_TRUNC+4)/8; i++) 
	{
		OP[i]=0;
		SOP[i]=0;

		OP2[i]=0;
		SOP2[i]=0;
	}
	if(Lmax==0)
	{
		OP[0]=1.0;
		SOP[0]=1.0;

		OP2[0]=1.0;
		SOP2[0]=1.0;

		return;
	}

	int GaussNodes=50; int KronrodNodes=101; 
	double Node[101], KronWeight[101], GausWeight[101];

		Node[ 50]=  -9.9981190136436471898752764191979778e-01; KronWeight[ 50]=5.0676166803489136800323005679418875e-04; GausWeight[ 50]=0                                       ;
		Node[  0]=  -9.9886640442007105018545944497421851e-01; KronWeight[  0]=1.4201102381663570584892416091048193e-03; GausWeight[  0]=2.9086225531551409584007243428554808e-03;
		Node[ 51]=  -9.9694438701887617830552521284551574e-01; KronWeight[ 51]=2.4231037458205732718680199597014328e-03; GausWeight[ 51]=0                                       ;
		Node[  1]=  -9.9403196943209071258510820042069473e-01; KronWeight[  1]=3.3945908928972374308566272798283577e-03; GausWeight[  1]=6.7597991957454015027788781779850318e-03;
		Node[ 52]=  -9.9016501066968000942392663658159575e-01; KronWeight[ 52]=4.3377036052637357835434236062529196e-03; GausWeight[ 52]=0                                       ;
		Node[  2]=  -9.8535408404800588230900962563248940e-01; KronWeight[  2]=5.2869070843107008149663109052265184e-03; GausWeight[  2]=1.0590548383650969263569681499241022e-02;
		Node[ 53]=  -9.7958729676076942938993428547746899e-01; KronWeight[ 53]=6.2467616528934688052784611929436288e-03; GausWeight[ 53]=0                                       ;
		Node[  3]=  -9.7286438510669207371334410460625205e-01; KronWeight[  3]=7.1958559137583317710477090115739725e-03; GausWeight[  3]=1.4380822761485574419378908927324350e-02;
		Node[ 54]=  -9.6520165106614520049209141719129594e-01; KronWeight[ 54]=8.1275336625593873786803729424595699e-03; GausWeight[ 54]=0                                       ;
		Node[  4]=  -9.5661095524280794299774564415662209e-01; KronWeight[  4]=9.0539068531320007348298463831867097e-03; GausWeight[  4]=1.8115560713489390351259943422354620e-02;
		Node[ 55]=  -9.4709403424493944376335887833533441e-01; KronWeight[ 55]=9.9790805127719727152000689655549135e-03; GausWeight[ 55]=0                                       ;
		Node[  5]=  -9.3665661894487793378087494727249660e-01; KronWeight[  5]=1.0893031149145421391853736320483791e-02; GausWeight[  5]=2.1780243170124792981592069062690341e-02;
		Node[ 56]=  -9.2531354647480188621729501218597896e-01; KronWeight[ 56]=1.1790732100542854676267459488075106e-02; GausWeight[ 56]=0                                       ;
		Node[  6]=  -9.1307855665579189308973564277165709e-01; KronWeight[  6]=1.2678061837012244403571737226388368e-02; GausWeight[  6]=2.5360673570012390440194878385442723e-02;
		Node[ 57]=  -8.9995988756429461789335121594253760e-01; KronWeight[ 57]=1.3557616011559140930458157803339780e-02; GausWeight[ 57]=0                                       ;
		Node[  7]=  -8.8596797952361304863754098246675363e-01; KronWeight[  7]=1.4423324872289546908736058736990278e-02; GausWeight[  7]=2.8842993580535198029906373113232433e-02;
		Node[ 58]=  -8.7111919820410031041116631276794735e-01; KronWeight[ 58]=1.5271463812173062587031639835809499e-02; GausWeight[ 58]=0                                       ;
		Node[  8]=  -8.5542976942994608461136264393475747e-01; KronWeight[  8]=1.6105363478293306328795864483333898e-02; GausWeight[  8]=3.2213728223578016648165827323003953e-02;
		Node[ 59]=  -8.3891258696722439582669453115982172e-01; KronWeight[ 59]=1.6926657875738541082113555463923612e-02; GausWeight[ 59]=0                                       ;
		Node[  9]=  -8.2158207085933594835625411087393954e-01; KronWeight[  9]=1.7731168466700687369451433287130214e-02; GausWeight[  9]=3.5459835615146154160734611000975797e-02;
		Node[ 60]=  -8.0345686805045941104583132344339868e-01; KronWeight[ 60]=1.8516042461777473742240817149347183e-02; GausWeight[ 60]=0                                       ;
		Node[ 10]=  -7.8455583290039926390530519634099120e-01; KronWeight[ 10]=1.9283323521808148253018898011621324e-02; GausWeight[ 10]=3.8568756612587675244770150236385935e-02;
		Node[ 61]=  -7.6489567937235152255213391009405332e-01; KronWeight[ 61]=2.0034045830719434860935189792575789e-02; GausWeight[ 61]=0                                       ;
		Node[ 11]=  -7.4449430222606853826053625268219424e-01; KronWeight[ 11]=2.0765128471719162242080239721667014e-02; GausWeight[ 11]=4.1528463090147697422411978964067018e-02;
		Node[ 62]=  -7.2337276602259251579221941156517307e-01; KronWeight[ 62]=2.1474338201809933451719273439764485e-02; GausWeight[ 62]=0                                       ;
		Node[ 12]=  -7.0155246870682225108954625788365573e-01; KronWeight[ 12]=2.2162984926601566106843498326675211e-02; GausWeight[ 12]=4.4327504338803275492022286830394197e-02;
		Node[ 63]=  -6.7905338940674727147750574623528946e-01; KronWeight[ 63]=2.2831717568874006991633399959920337e-02; GausWeight[ 63]=0                                       ;
		Node[ 13]=  -6.5589646568543936078162486400367982e-01; KronWeight[ 13]=2.3478184078910193004586090503346223e-02; GausWeight[ 13]=4.6955051303948432965633013634987683e-02;
		Node[ 64]=  -6.3210506844506487645226444165317495e-01; KronWeight[ 64]=2.4100625877218305976665760144698111e-02; GausWeight[ 64]=0                                       ;
		Node[ 14]=  -6.0770292718495023918038179639183289e-01; KronWeight[ 14]=2.4699903537325739494760209714837493e-02; GausWeight[ 14]=4.9400938449466314921243580751432729e-02;
		Node[ 65]=  -5.8271281781783653230059341928604869e-01; KronWeight[ 65]=2.5276412574832709146478978607374462e-02; GausWeight[ 65]=0                                       ;
		Node[ 15]=  -5.5715830451465005431552290962580161e-01; KronWeight[ 15]=2.5828336809121795208844863774061424e-02; GausWeight[ 15]=5.1655703069581138489905295840095280e-02;
		Node[ 66]=  -5.3106482452270813956871444365352729e-01; KronWeight[ 66]=2.6354306818690138454725762609577214e-02; GausWeight[ 66]=0                                       ;
		Node[ 16]=  -5.0445814490746420165145913184914119e-01; KronWeight[ 16]=2.6854896362906597368937147115205512e-02; GausWeight[ 16]=5.3710621888996246523458797255664553e-02;
		Node[ 67]=  -4.7736339223304303620238841113348793e-01; KronWeight[ 67]=2.7330331890109643833054043505441732e-02; GausWeight[ 67]=0                                       ;
		Node[ 17]=  -4.4980633497403878914713146777837582e-01; KronWeight[ 17]=2.7779223892439833361406855758959752e-02; GausWeight[ 17]=5.5557744806212517623567425612269498e-02;
		Node[ 68]=  -4.2181415731061309702869639767299001e-01; KronWeight[ 68]=2.8200540579478679146298601877461934e-02; GausWeight[ 68]=0                                       ;
		Node[ 18]=  -3.9341431189756512739422925382381727e-01; KronWeight[ 18]=2.8594668383406057739047317910104991e-02; GausWeight[ 18]=5.7189925647728383723029315065993163e-02;
		Node[ 69]=  -3.6463382886161427721485091875517339e-01; KronWeight[ 69]=2.8961719790501033567571417319566858e-02; GausWeight[ 69]=0                                       ;
		Node[ 19]=  -3.3550024541943735683698825729107170e-01; KronWeight[ 19]=2.9300667035119782336840778651545099e-02; GausWeight[ 19]=5.8600849813222445835122436630848466e-02;
		Node[ 70]=  -3.0604211922918460908352322458238267e-01; KronWeight[ 70]=2.9610784940459955030252135399656593e-02; GausWeight[ 70]=0                                       ;
		Node[ 20]=  -2.7628819377953199032764527852113019e-01; KronWeight[ 20]=2.9892335817039792458804648308739805e-02; GausWeight[ 20]=5.9785058704265457509576405312585231e-02;
		Node[ 71]=  -2.4626694739814402597841233060681906e-01; KronWeight[ 71]=3.0145353774570835194075130209130173e-02; GausWeight[ 71]=0                                       ;
		Node[ 21]=  -2.1600723687604175684728453261710133e-01; KronWeight[ 21]=3.0369133200164873181378495665853924e-02; GausWeight[ 21]=6.0737970841770216031750015384811002e-02;
		Node[ 72]=  -1.8553858172277277270028588015893298e-01; KronWeight[ 72]=3.0563237706310037205460535834084086e-02; GausWeight[ 72]=0                                       ;
		Node[ 22]=  -1.5489058999814590207162862094110950e-01; KronWeight[ 22]=3.0727845669416690205934370165811952e-02; GausWeight[ 22]=6.1455899590316663756406786083915375e-02;
		Node[ 73]=  -1.2409272435916037178300113567216009e-01; KronWeight[ 73]=3.0862934121709973760264323329606162e-02; GausWeight[ 73]=0                                       ;
		Node[ 23]=  -9.3174701560086140854450377639600348e-02; KronWeight[ 23]=3.0968095623668437958542109562984134e-02; GausWeight[ 23]=6.1936067420683243384087509780830689e-02;
		Node[ 74]=  -6.2166564819416169080162369062493662e-02; KronWeight[ 74]=3.1043173996945701068760645698647835e-02; GausWeight[ 74]=0                                       ;
		Node[ 24]=  -3.1098338327188876112328989665949194e-02; KronWeight[ 24]=3.1088287782405219104218184906343066e-02; GausWeight[ 24]=6.2176616655347262321033107360613431e-02;
		Node[ 75]=   3.2073176521106347753686437605817783e-50; KronWeight[ 75]=3.1103366641749575467154644934577086e-02; GausWeight[ 75]=0                                       ;
		Node[ 25]=   3.1098338327188876112328989665949194e-02; KronWeight[ 25]=3.1088287782405219104218184906343066e-02; GausWeight[ 25]=6.2176616655347262321033107360613431e-02;
		Node[ 76]=   6.2166564819416169080162369062493662e-02; KronWeight[ 76]=3.1043173996945701068760645698647835e-02; GausWeight[ 76]=0                                       ;
		Node[ 26]=   9.3174701560086140854450377639600348e-02; KronWeight[ 26]=3.0968095623668437958542109562984134e-02; GausWeight[ 26]=6.1936067420683243384087509780830689e-02;
		Node[ 77]=   1.2409272435916037178300113567216009e-01; KronWeight[ 77]=3.0862934121709973760264323329606162e-02; GausWeight[ 77]=0                                       ;
		Node[ 27]=   1.5489058999814590207162862094110950e-01; KronWeight[ 27]=3.0727845669416690205934370165811952e-02; GausWeight[ 27]=6.1455899590316663756406786083915375e-02;
		Node[ 78]=   1.8553858172277277270028588015893298e-01; KronWeight[ 78]=3.0563237706310037205460535834084086e-02; GausWeight[ 78]=0                                       ;
		Node[ 28]=   2.1600723687604175684728453261710133e-01; KronWeight[ 28]=3.0369133200164873181378495665853924e-02; GausWeight[ 28]=6.0737970841770216031750015384811002e-02;
		Node[ 79]=   2.4626694739814402597841233060681906e-01; KronWeight[ 79]=3.0145353774570835194075130209130173e-02; GausWeight[ 79]=0                                       ;
		Node[ 29]=   2.7628819377953199032764527852113019e-01; KronWeight[ 29]=2.9892335817039792458804648308739805e-02; GausWeight[ 29]=5.9785058704265457509576405312585231e-02;
		Node[ 80]=   3.0604211922918460908352322458238267e-01; KronWeight[ 80]=2.9610784940459955030252135399656593e-02; GausWeight[ 80]=0                                       ;
		Node[ 30]=   3.3550024541943735683698825729107170e-01; KronWeight[ 30]=2.9300667035119782336840778651545099e-02; GausWeight[ 30]=5.8600849813222445835122436630848466e-02;
		Node[ 81]=   3.6463382886161427721485091875517339e-01; KronWeight[ 81]=2.8961719790501033567571417319566858e-02; GausWeight[ 81]=0                                       ;
		Node[ 31]=   3.9341431189756512739422925382381727e-01; KronWeight[ 31]=2.8594668383406057739047317910104991e-02; GausWeight[ 31]=5.7189925647728383723029315065993163e-02;
		Node[ 82]=   4.2181415731061309702869639767299001e-01; KronWeight[ 82]=2.8200540579478679146298601877461934e-02; GausWeight[ 82]=0                                       ;
		Node[ 32]=   4.4980633497403878914713146777837582e-01; KronWeight[ 32]=2.7779223892439833361406855758959752e-02; GausWeight[ 32]=5.5557744806212517623567425612269498e-02;
		Node[ 83]=   4.7736339223304303620238841113348793e-01; KronWeight[ 83]=2.7330331890109643833054043505441732e-02; GausWeight[ 83]=0                                       ;
		Node[ 33]=   5.0445814490746420165145913184914119e-01; KronWeight[ 33]=2.6854896362906597368937147115205512e-02; GausWeight[ 33]=5.3710621888996246523458797255664553e-02;
		Node[ 84]=   5.3106482452270813956871444365352729e-01; KronWeight[ 84]=2.6354306818690138454725762609577214e-02; GausWeight[ 84]=0                                       ;
		Node[ 34]=   5.5715830451465005431552290962580161e-01; KronWeight[ 34]=2.5828336809121795208844863774061424e-02; GausWeight[ 34]=5.1655703069581138489905295840095280e-02;
		Node[ 85]=   5.8271281781783653230059341928604869e-01; KronWeight[ 85]=2.5276412574832709146478978607374462e-02; GausWeight[ 85]=0                                       ;
		Node[ 35]=   6.0770292718495023918038179639183289e-01; KronWeight[ 35]=2.4699903537325739494760209714837493e-02; GausWeight[ 35]=4.9400938449466314921243580751432729e-02;
		Node[ 86]=   6.3210506844506487645226444165317495e-01; KronWeight[ 86]=2.4100625877218305976665760144698111e-02; GausWeight[ 86]=0                                       ;
		Node[ 36]=   6.5589646568543936078162486400367982e-01; KronWeight[ 36]=2.3478184078910193004586090503346223e-02; GausWeight[ 36]=4.6955051303948432965633013634987683e-02;
		Node[ 87]=   6.7905338940674727147750574623528946e-01; KronWeight[ 87]=2.2831717568874006991633399959920337e-02; GausWeight[ 87]=0                                       ;
		Node[ 37]=   7.0155246870682225108954625788365573e-01; KronWeight[ 37]=2.2162984926601566106843498326675211e-02; GausWeight[ 37]=4.4327504338803275492022286830394197e-02;
		Node[ 88]=   7.2337276602259251579221941156517307e-01; KronWeight[ 88]=2.1474338201809933451719273439764485e-02; GausWeight[ 88]=0                                       ;
		Node[ 38]=   7.4449430222606853826053625268219424e-01; KronWeight[ 38]=2.0765128471719162242080239721667014e-02; GausWeight[ 38]=4.1528463090147697422411978964067018e-02;
		Node[ 89]=   7.6489567937235152255213391009405332e-01; KronWeight[ 89]=2.0034045830719434860935189792575789e-02; GausWeight[ 89]=0                                       ;
		Node[ 39]=   7.8455583290039926390530519634099120e-01; KronWeight[ 39]=1.9283323521808148253018898011621324e-02; GausWeight[ 39]=3.8568756612587675244770150236385935e-02;
		Node[ 90]=   8.0345686805045941104583132344339868e-01; KronWeight[ 90]=1.8516042461777473742240817149347183e-02; GausWeight[ 90]=0                                       ;
		Node[ 40]=   8.2158207085933594835625411087393954e-01; KronWeight[ 40]=1.7731168466700687369451433287130214e-02; GausWeight[ 40]=3.5459835615146154160734611000975797e-02;
		Node[ 91]=   8.3891258696722439582669453115982172e-01; KronWeight[ 91]=1.6926657875738541082113555463923612e-02; GausWeight[ 91]=0                                       ;
		Node[ 41]=   8.5542976942994608461136264393475747e-01; KronWeight[ 41]=1.6105363478293306328795864483333898e-02; GausWeight[ 41]=3.2213728223578016648165827323003953e-02;
		Node[ 92]=   8.7111919820410031041116631276794735e-01; KronWeight[ 92]=1.5271463812173062587031639835809499e-02; GausWeight[ 92]=0                                       ;
		Node[ 42]=   8.8596797952361304863754098246675363e-01; KronWeight[ 42]=1.4423324872289546908736058736990278e-02; GausWeight[ 42]=2.8842993580535198029906373113232433e-02;
		Node[ 93]=   8.9995988756429461789335121594253760e-01; KronWeight[ 93]=1.3557616011559140930458157803339780e-02; GausWeight[ 93]=0                                       ;
		Node[ 43]=   9.1307855665579189308973564277165709e-01; KronWeight[ 43]=1.2678061837012244403571737226388368e-02; GausWeight[ 43]=2.5360673570012390440194878385442723e-02;
		Node[ 94]=   9.2531354647480188621729501218597896e-01; KronWeight[ 94]=1.1790732100542854676267459488075106e-02; GausWeight[ 94]=0                                       ;
		Node[ 44]=   9.3665661894487793378087494727249660e-01; KronWeight[ 44]=1.0893031149145421391853736320483791e-02; GausWeight[ 44]=2.1780243170124792981592069062690341e-02;
		Node[ 95]=   9.4709403424493944376335887833533441e-01; KronWeight[ 95]=9.9790805127719727152000689655549135e-03; GausWeight[ 95]=0                                       ;
		Node[ 45]=   9.5661095524280794299774564415662209e-01; KronWeight[ 45]=9.0539068531320007348298463831867097e-03; GausWeight[ 45]=1.8115560713489390351259943422354620e-02;
		Node[ 96]=   9.6520165106614520049209141719129594e-01; KronWeight[ 96]=8.1275336625593873786803729424595699e-03; GausWeight[ 96]=0                                       ;
		Node[ 46]=   9.7286438510669207371334410460625205e-01; KronWeight[ 46]=7.1958559137583317710477090115739725e-03; GausWeight[ 46]=1.4380822761485574419378908927324350e-02;
		Node[ 97]=   9.7958729676076942938993428547746899e-01; KronWeight[ 97]=6.2467616528934688052784611929436288e-03; GausWeight[ 97]=0                                       ;
		Node[ 47]=   9.8535408404800588230900962563248940e-01; KronWeight[ 47]=5.2869070843107008149663109052265184e-03; GausWeight[ 47]=1.0590548383650969263569681499241022e-02;
		Node[ 98]=   9.9016501066968000942392663658159575e-01; KronWeight[ 98]=4.3377036052637357835434236062529196e-03; GausWeight[ 98]=0                                       ;
		Node[ 48]=   9.9403196943209071258510820042069473e-01; KronWeight[ 48]=3.3945908928972374308566272798283577e-03; GausWeight[ 48]=6.7597991957454015027788781779850318e-03;
		Node[ 99]=   9.9694438701887617830552521284551574e-01; KronWeight[ 99]=2.4231037458205732718680199597014328e-03; GausWeight[ 99]=0                                       ;
		Node[ 49]=   9.9886640442007105018545944497421851e-01; KronWeight[ 49]=1.4201102381663570584892416091048193e-03; GausWeight[ 49]=2.9086225531551409584007243428554808e-03;
		Node[100]=   9.9981190136436471898752764191979778e-01; KronWeight[100]=5.0676166803489136800323005679418875e-04; GausWeight[100]=0                                       ;



	for(int ibeta=0; ibeta<KronrodNodes; ibeta++)
	{
		double cosbeta=Node[ibeta];
		double beta=acos(cosbeta);

	  for(int igamma=0; igamma<KronrodNodes; igamma++)
	  {
	    double gamma=(Node[igamma]+1.0)*PI;
		double ro=0;
		for (int l=0; l<=Lmax; l+=2)
		{
// calculating ro=potential!
			ro+=this->Get(l, 0)*real(ComplexWigner(l, 0, 0, 0, beta, gamma) );
			if(!isAxial)
			  for(int k=2; k<=l; k+=2)
				ro+=this->Get(l, k)*2.0*real(ComplexWigner(l, 0, k, 0, beta, gamma) );
		}
// now ro becomes boltzman orientation distribution function:
			double ro2=exp(ro);
			ro=exp(0.5*ro);
//			cout << ro << " " << ro2 <<endl;
			for( int L=0; L<=2*MaxL; L+=2 )
			  for( int K=0; K<=L; K+=2 )
				if((!isAxial)||(K==0))
				{
					OP[L*(L+2)/8+K/2] += real(ComplexWigner(L, 0, K, 0, beta, gamma))*ro2*KronWeight[ibeta]*KronWeight[igamma];	// integration step no *sin(beta) because cosine integration!
					SOP[L*(L+2)/8+K/2]+= real(ComplexWigner(L, 0, K, 0, beta, gamma))*ro *KronWeight[ibeta]*KronWeight[igamma];	// integration step

					OP2[L*(L+2)/8+K/2] += real(ComplexWigner(L, 0, K, 0, beta, gamma))*ro2*GausWeight[ibeta]*GausWeight[igamma];	// integration step no *sin(beta) because cosine integration!
					SOP2[L*(L+2)/8+K/2]+= real(ComplexWigner(L, 0, K, 0, beta, gamma))*ro *GausWeight[ibeta]*GausWeight[igamma];	// integration step

				}



	  }
	}
	  
// normalization:
//	  cout << MaxL << endl;
	  Normalization_OP=OP[0];
	  for (int i=1; i<(2*MaxL+2)*(2*MaxL+4)/8; i++) 
	  {
		SOP[i]/=SOP[0];
		OP[i] /=OP[0];
		SOP2[i]/=SOP2[0];
		OP2[i] /=OP2[0];

	  }
	  SOP[0]=1.0;
	  OP[0]= 1.0;

	  SOP2[0]=1.0;
	  OP2[0]= 1.0;

//	cout << "Finished!\n";
}



void Potential::CalculateOP_Jumps_GK(int MaxL)
{
for(int i=0; i<MagneticParameterParent->NJumpAxes; i++)
{
//	cout << "Calculating potential order parameters\t";
	for (int j=0; j<(MAX_TRUNC)*(MAX_TRUNC); j++) 
	{
		OP_JUMP_P[i][j]=0;
		OP_JUMP_M[i][j]=0;

		OP_JUMP_P2[i][j]=0;
		OP_JUMP_M2[i][j]=0;
	}
	if(Lmax==0)
	{
		OP_JUMP_P[i][0]=1.0;
		OP_JUMP_M[i][0]=1.0;

		OP_JUMP_P2[i][0]=1.0;
		OP_JUMP_M2[i][0]=1.0;

		return;
	}

	int GaussNodes=50; int KronrodNodes=101; 
	double Node[101], KronWeight[101], GausWeight[101];

		Node[ 50]=  -9.9981190136436471898752764191979778e-01; KronWeight[ 50]=5.0676166803489136800323005679418875e-04; GausWeight[ 50]=0                                       ;
		Node[  0]=  -9.9886640442007105018545944497421851e-01; KronWeight[  0]=1.4201102381663570584892416091048193e-03; GausWeight[  0]=2.9086225531551409584007243428554808e-03;
		Node[ 51]=  -9.9694438701887617830552521284551574e-01; KronWeight[ 51]=2.4231037458205732718680199597014328e-03; GausWeight[ 51]=0                                       ;
		Node[  1]=  -9.9403196943209071258510820042069473e-01; KronWeight[  1]=3.3945908928972374308566272798283577e-03; GausWeight[  1]=6.7597991957454015027788781779850318e-03;
		Node[ 52]=  -9.9016501066968000942392663658159575e-01; KronWeight[ 52]=4.3377036052637357835434236062529196e-03; GausWeight[ 52]=0                                       ;
		Node[  2]=  -9.8535408404800588230900962563248940e-01; KronWeight[  2]=5.2869070843107008149663109052265184e-03; GausWeight[  2]=1.0590548383650969263569681499241022e-02;
		Node[ 53]=  -9.7958729676076942938993428547746899e-01; KronWeight[ 53]=6.2467616528934688052784611929436288e-03; GausWeight[ 53]=0                                       ;
		Node[  3]=  -9.7286438510669207371334410460625205e-01; KronWeight[  3]=7.1958559137583317710477090115739725e-03; GausWeight[  3]=1.4380822761485574419378908927324350e-02;
		Node[ 54]=  -9.6520165106614520049209141719129594e-01; KronWeight[ 54]=8.1275336625593873786803729424595699e-03; GausWeight[ 54]=0                                       ;
		Node[  4]=  -9.5661095524280794299774564415662209e-01; KronWeight[  4]=9.0539068531320007348298463831867097e-03; GausWeight[  4]=1.8115560713489390351259943422354620e-02;
		Node[ 55]=  -9.4709403424493944376335887833533441e-01; KronWeight[ 55]=9.9790805127719727152000689655549135e-03; GausWeight[ 55]=0                                       ;
		Node[  5]=  -9.3665661894487793378087494727249660e-01; KronWeight[  5]=1.0893031149145421391853736320483791e-02; GausWeight[  5]=2.1780243170124792981592069062690341e-02;
		Node[ 56]=  -9.2531354647480188621729501218597896e-01; KronWeight[ 56]=1.1790732100542854676267459488075106e-02; GausWeight[ 56]=0                                       ;
		Node[  6]=  -9.1307855665579189308973564277165709e-01; KronWeight[  6]=1.2678061837012244403571737226388368e-02; GausWeight[  6]=2.5360673570012390440194878385442723e-02;
		Node[ 57]=  -8.9995988756429461789335121594253760e-01; KronWeight[ 57]=1.3557616011559140930458157803339780e-02; GausWeight[ 57]=0                                       ;
		Node[  7]=  -8.8596797952361304863754098246675363e-01; KronWeight[  7]=1.4423324872289546908736058736990278e-02; GausWeight[  7]=2.8842993580535198029906373113232433e-02;
		Node[ 58]=  -8.7111919820410031041116631276794735e-01; KronWeight[ 58]=1.5271463812173062587031639835809499e-02; GausWeight[ 58]=0                                       ;
		Node[  8]=  -8.5542976942994608461136264393475747e-01; KronWeight[  8]=1.6105363478293306328795864483333898e-02; GausWeight[  8]=3.2213728223578016648165827323003953e-02;
		Node[ 59]=  -8.3891258696722439582669453115982172e-01; KronWeight[ 59]=1.6926657875738541082113555463923612e-02; GausWeight[ 59]=0                                       ;
		Node[  9]=  -8.2158207085933594835625411087393954e-01; KronWeight[  9]=1.7731168466700687369451433287130214e-02; GausWeight[  9]=3.5459835615146154160734611000975797e-02;
		Node[ 60]=  -8.0345686805045941104583132344339868e-01; KronWeight[ 60]=1.8516042461777473742240817149347183e-02; GausWeight[ 60]=0                                       ;
		Node[ 10]=  -7.8455583290039926390530519634099120e-01; KronWeight[ 10]=1.9283323521808148253018898011621324e-02; GausWeight[ 10]=3.8568756612587675244770150236385935e-02;
		Node[ 61]=  -7.6489567937235152255213391009405332e-01; KronWeight[ 61]=2.0034045830719434860935189792575789e-02; GausWeight[ 61]=0                                       ;
		Node[ 11]=  -7.4449430222606853826053625268219424e-01; KronWeight[ 11]=2.0765128471719162242080239721667014e-02; GausWeight[ 11]=4.1528463090147697422411978964067018e-02;
		Node[ 62]=  -7.2337276602259251579221941156517307e-01; KronWeight[ 62]=2.1474338201809933451719273439764485e-02; GausWeight[ 62]=0                                       ;
		Node[ 12]=  -7.0155246870682225108954625788365573e-01; KronWeight[ 12]=2.2162984926601566106843498326675211e-02; GausWeight[ 12]=4.4327504338803275492022286830394197e-02;
		Node[ 63]=  -6.7905338940674727147750574623528946e-01; KronWeight[ 63]=2.2831717568874006991633399959920337e-02; GausWeight[ 63]=0                                       ;
		Node[ 13]=  -6.5589646568543936078162486400367982e-01; KronWeight[ 13]=2.3478184078910193004586090503346223e-02; GausWeight[ 13]=4.6955051303948432965633013634987683e-02;
		Node[ 64]=  -6.3210506844506487645226444165317495e-01; KronWeight[ 64]=2.4100625877218305976665760144698111e-02; GausWeight[ 64]=0                                       ;
		Node[ 14]=  -6.0770292718495023918038179639183289e-01; KronWeight[ 14]=2.4699903537325739494760209714837493e-02; GausWeight[ 14]=4.9400938449466314921243580751432729e-02;
		Node[ 65]=  -5.8271281781783653230059341928604869e-01; KronWeight[ 65]=2.5276412574832709146478978607374462e-02; GausWeight[ 65]=0                                       ;
		Node[ 15]=  -5.5715830451465005431552290962580161e-01; KronWeight[ 15]=2.5828336809121795208844863774061424e-02; GausWeight[ 15]=5.1655703069581138489905295840095280e-02;
		Node[ 66]=  -5.3106482452270813956871444365352729e-01; KronWeight[ 66]=2.6354306818690138454725762609577214e-02; GausWeight[ 66]=0                                       ;
		Node[ 16]=  -5.0445814490746420165145913184914119e-01; KronWeight[ 16]=2.6854896362906597368937147115205512e-02; GausWeight[ 16]=5.3710621888996246523458797255664553e-02;
		Node[ 67]=  -4.7736339223304303620238841113348793e-01; KronWeight[ 67]=2.7330331890109643833054043505441732e-02; GausWeight[ 67]=0                                       ;
		Node[ 17]=  -4.4980633497403878914713146777837582e-01; KronWeight[ 17]=2.7779223892439833361406855758959752e-02; GausWeight[ 17]=5.5557744806212517623567425612269498e-02;
		Node[ 68]=  -4.2181415731061309702869639767299001e-01; KronWeight[ 68]=2.8200540579478679146298601877461934e-02; GausWeight[ 68]=0                                       ;
		Node[ 18]=  -3.9341431189756512739422925382381727e-01; KronWeight[ 18]=2.8594668383406057739047317910104991e-02; GausWeight[ 18]=5.7189925647728383723029315065993163e-02;
		Node[ 69]=  -3.6463382886161427721485091875517339e-01; KronWeight[ 69]=2.8961719790501033567571417319566858e-02; GausWeight[ 69]=0                                       ;
		Node[ 19]=  -3.3550024541943735683698825729107170e-01; KronWeight[ 19]=2.9300667035119782336840778651545099e-02; GausWeight[ 19]=5.8600849813222445835122436630848466e-02;
		Node[ 70]=  -3.0604211922918460908352322458238267e-01; KronWeight[ 70]=2.9610784940459955030252135399656593e-02; GausWeight[ 70]=0                                       ;
		Node[ 20]=  -2.7628819377953199032764527852113019e-01; KronWeight[ 20]=2.9892335817039792458804648308739805e-02; GausWeight[ 20]=5.9785058704265457509576405312585231e-02;
		Node[ 71]=  -2.4626694739814402597841233060681906e-01; KronWeight[ 71]=3.0145353774570835194075130209130173e-02; GausWeight[ 71]=0                                       ;
		Node[ 21]=  -2.1600723687604175684728453261710133e-01; KronWeight[ 21]=3.0369133200164873181378495665853924e-02; GausWeight[ 21]=6.0737970841770216031750015384811002e-02;
		Node[ 72]=  -1.8553858172277277270028588015893298e-01; KronWeight[ 72]=3.0563237706310037205460535834084086e-02; GausWeight[ 72]=0                                       ;
		Node[ 22]=  -1.5489058999814590207162862094110950e-01; KronWeight[ 22]=3.0727845669416690205934370165811952e-02; GausWeight[ 22]=6.1455899590316663756406786083915375e-02;
		Node[ 73]=  -1.2409272435916037178300113567216009e-01; KronWeight[ 73]=3.0862934121709973760264323329606162e-02; GausWeight[ 73]=0                                       ;
		Node[ 23]=  -9.3174701560086140854450377639600348e-02; KronWeight[ 23]=3.0968095623668437958542109562984134e-02; GausWeight[ 23]=6.1936067420683243384087509780830689e-02;
		Node[ 74]=  -6.2166564819416169080162369062493662e-02; KronWeight[ 74]=3.1043173996945701068760645698647835e-02; GausWeight[ 74]=0                                       ;
		Node[ 24]=  -3.1098338327188876112328989665949194e-02; KronWeight[ 24]=3.1088287782405219104218184906343066e-02; GausWeight[ 24]=6.2176616655347262321033107360613431e-02;
		Node[ 75]=   3.2073176521106347753686437605817783e-50; KronWeight[ 75]=3.1103366641749575467154644934577086e-02; GausWeight[ 75]=0                                       ;
		Node[ 25]=   3.1098338327188876112328989665949194e-02; KronWeight[ 25]=3.1088287782405219104218184906343066e-02; GausWeight[ 25]=6.2176616655347262321033107360613431e-02;
		Node[ 76]=   6.2166564819416169080162369062493662e-02; KronWeight[ 76]=3.1043173996945701068760645698647835e-02; GausWeight[ 76]=0                                       ;
		Node[ 26]=   9.3174701560086140854450377639600348e-02; KronWeight[ 26]=3.0968095623668437958542109562984134e-02; GausWeight[ 26]=6.1936067420683243384087509780830689e-02;
		Node[ 77]=   1.2409272435916037178300113567216009e-01; KronWeight[ 77]=3.0862934121709973760264323329606162e-02; GausWeight[ 77]=0                                       ;
		Node[ 27]=   1.5489058999814590207162862094110950e-01; KronWeight[ 27]=3.0727845669416690205934370165811952e-02; GausWeight[ 27]=6.1455899590316663756406786083915375e-02;
		Node[ 78]=   1.8553858172277277270028588015893298e-01; KronWeight[ 78]=3.0563237706310037205460535834084086e-02; GausWeight[ 78]=0                                       ;
		Node[ 28]=   2.1600723687604175684728453261710133e-01; KronWeight[ 28]=3.0369133200164873181378495665853924e-02; GausWeight[ 28]=6.0737970841770216031750015384811002e-02;
		Node[ 79]=   2.4626694739814402597841233060681906e-01; KronWeight[ 79]=3.0145353774570835194075130209130173e-02; GausWeight[ 79]=0                                       ;
		Node[ 29]=   2.7628819377953199032764527852113019e-01; KronWeight[ 29]=2.9892335817039792458804648308739805e-02; GausWeight[ 29]=5.9785058704265457509576405312585231e-02;
		Node[ 80]=   3.0604211922918460908352322458238267e-01; KronWeight[ 80]=2.9610784940459955030252135399656593e-02; GausWeight[ 80]=0                                       ;
		Node[ 30]=   3.3550024541943735683698825729107170e-01; KronWeight[ 30]=2.9300667035119782336840778651545099e-02; GausWeight[ 30]=5.8600849813222445835122436630848466e-02;
		Node[ 81]=   3.6463382886161427721485091875517339e-01; KronWeight[ 81]=2.8961719790501033567571417319566858e-02; GausWeight[ 81]=0                                       ;
		Node[ 31]=   3.9341431189756512739422925382381727e-01; KronWeight[ 31]=2.8594668383406057739047317910104991e-02; GausWeight[ 31]=5.7189925647728383723029315065993163e-02;
		Node[ 82]=   4.2181415731061309702869639767299001e-01; KronWeight[ 82]=2.8200540579478679146298601877461934e-02; GausWeight[ 82]=0                                       ;
		Node[ 32]=   4.4980633497403878914713146777837582e-01; KronWeight[ 32]=2.7779223892439833361406855758959752e-02; GausWeight[ 32]=5.5557744806212517623567425612269498e-02;
		Node[ 83]=   4.7736339223304303620238841113348793e-01; KronWeight[ 83]=2.7330331890109643833054043505441732e-02; GausWeight[ 83]=0                                       ;
		Node[ 33]=   5.0445814490746420165145913184914119e-01; KronWeight[ 33]=2.6854896362906597368937147115205512e-02; GausWeight[ 33]=5.3710621888996246523458797255664553e-02;
		Node[ 84]=   5.3106482452270813956871444365352729e-01; KronWeight[ 84]=2.6354306818690138454725762609577214e-02; GausWeight[ 84]=0                                       ;
		Node[ 34]=   5.5715830451465005431552290962580161e-01; KronWeight[ 34]=2.5828336809121795208844863774061424e-02; GausWeight[ 34]=5.1655703069581138489905295840095280e-02;
		Node[ 85]=   5.8271281781783653230059341928604869e-01; KronWeight[ 85]=2.5276412574832709146478978607374462e-02; GausWeight[ 85]=0                                       ;
		Node[ 35]=   6.0770292718495023918038179639183289e-01; KronWeight[ 35]=2.4699903537325739494760209714837493e-02; GausWeight[ 35]=4.9400938449466314921243580751432729e-02;
		Node[ 86]=   6.3210506844506487645226444165317495e-01; KronWeight[ 86]=2.4100625877218305976665760144698111e-02; GausWeight[ 86]=0                                       ;
		Node[ 36]=   6.5589646568543936078162486400367982e-01; KronWeight[ 36]=2.3478184078910193004586090503346223e-02; GausWeight[ 36]=4.6955051303948432965633013634987683e-02;
		Node[ 87]=   6.7905338940674727147750574623528946e-01; KronWeight[ 87]=2.2831717568874006991633399959920337e-02; GausWeight[ 87]=0                                       ;
		Node[ 37]=   7.0155246870682225108954625788365573e-01; KronWeight[ 37]=2.2162984926601566106843498326675211e-02; GausWeight[ 37]=4.4327504338803275492022286830394197e-02;
		Node[ 88]=   7.2337276602259251579221941156517307e-01; KronWeight[ 88]=2.1474338201809933451719273439764485e-02; GausWeight[ 88]=0                                       ;
		Node[ 38]=   7.4449430222606853826053625268219424e-01; KronWeight[ 38]=2.0765128471719162242080239721667014e-02; GausWeight[ 38]=4.1528463090147697422411978964067018e-02;
		Node[ 89]=   7.6489567937235152255213391009405332e-01; KronWeight[ 89]=2.0034045830719434860935189792575789e-02; GausWeight[ 89]=0                                       ;
		Node[ 39]=   7.8455583290039926390530519634099120e-01; KronWeight[ 39]=1.9283323521808148253018898011621324e-02; GausWeight[ 39]=3.8568756612587675244770150236385935e-02;
		Node[ 90]=   8.0345686805045941104583132344339868e-01; KronWeight[ 90]=1.8516042461777473742240817149347183e-02; GausWeight[ 90]=0                                       ;
		Node[ 40]=   8.2158207085933594835625411087393954e-01; KronWeight[ 40]=1.7731168466700687369451433287130214e-02; GausWeight[ 40]=3.5459835615146154160734611000975797e-02;
		Node[ 91]=   8.3891258696722439582669453115982172e-01; KronWeight[ 91]=1.6926657875738541082113555463923612e-02; GausWeight[ 91]=0                                       ;
		Node[ 41]=   8.5542976942994608461136264393475747e-01; KronWeight[ 41]=1.6105363478293306328795864483333898e-02; GausWeight[ 41]=3.2213728223578016648165827323003953e-02;
		Node[ 92]=   8.7111919820410031041116631276794735e-01; KronWeight[ 92]=1.5271463812173062587031639835809499e-02; GausWeight[ 92]=0                                       ;
		Node[ 42]=   8.8596797952361304863754098246675363e-01; KronWeight[ 42]=1.4423324872289546908736058736990278e-02; GausWeight[ 42]=2.8842993580535198029906373113232433e-02;
		Node[ 93]=   8.9995988756429461789335121594253760e-01; KronWeight[ 93]=1.3557616011559140930458157803339780e-02; GausWeight[ 93]=0                                       ;
		Node[ 43]=   9.1307855665579189308973564277165709e-01; KronWeight[ 43]=1.2678061837012244403571737226388368e-02; GausWeight[ 43]=2.5360673570012390440194878385442723e-02;
		Node[ 94]=   9.2531354647480188621729501218597896e-01; KronWeight[ 94]=1.1790732100542854676267459488075106e-02; GausWeight[ 94]=0                                       ;
		Node[ 44]=   9.3665661894487793378087494727249660e-01; KronWeight[ 44]=1.0893031149145421391853736320483791e-02; GausWeight[ 44]=2.1780243170124792981592069062690341e-02;
		Node[ 95]=   9.4709403424493944376335887833533441e-01; KronWeight[ 95]=9.9790805127719727152000689655549135e-03; GausWeight[ 95]=0                                       ;
		Node[ 45]=   9.5661095524280794299774564415662209e-01; KronWeight[ 45]=9.0539068531320007348298463831867097e-03; GausWeight[ 45]=1.8115560713489390351259943422354620e-02;
		Node[ 96]=   9.6520165106614520049209141719129594e-01; KronWeight[ 96]=8.1275336625593873786803729424595699e-03; GausWeight[ 96]=0                                       ;
		Node[ 46]=   9.7286438510669207371334410460625205e-01; KronWeight[ 46]=7.1958559137583317710477090115739725e-03; GausWeight[ 46]=1.4380822761485574419378908927324350e-02;
		Node[ 97]=   9.7958729676076942938993428547746899e-01; KronWeight[ 97]=6.2467616528934688052784611929436288e-03; GausWeight[ 97]=0                                       ;
		Node[ 47]=   9.8535408404800588230900962563248940e-01; KronWeight[ 47]=5.2869070843107008149663109052265184e-03; GausWeight[ 47]=1.0590548383650969263569681499241022e-02;
		Node[ 98]=   9.9016501066968000942392663658159575e-01; KronWeight[ 98]=4.3377036052637357835434236062529196e-03; GausWeight[ 98]=0                                       ;
		Node[ 48]=   9.9403196943209071258510820042069473e-01; KronWeight[ 48]=3.3945908928972374308566272798283577e-03; GausWeight[ 48]=6.7597991957454015027788781779850318e-03;
		Node[ 99]=   9.9694438701887617830552521284551574e-01; KronWeight[ 99]=2.4231037458205732718680199597014328e-03; GausWeight[ 99]=0                                       ;
		Node[ 49]=   9.9886640442007105018545944497421851e-01; KronWeight[ 49]=1.4201102381663570584892416091048193e-03; GausWeight[ 49]=2.9086225531551409584007243428554808e-03;
		Node[100]=   9.9981190136436471898752764191979778e-01; KronWeight[100]=5.0676166803489136800323005679418875e-04; GausWeight[100]=0                                       ;

	double fi=		MagneticParameterParent->JumpAxisFi[i];
	double teta=	MagneticParameterParent->JumpAxisTeta[i];
	double alpha=	MagneticParameterParent->JumpAngle[i];


	for(int ibeta=0; ibeta<KronrodNodes; ibeta++)
	{
		double cosbeta=Node[ibeta];
		double beta=acos(cosbeta);

	  for(int igamma=0; igamma<KronrodNodes; igamma++)
	  {
	    double gamma=(Node[igamma]+1.0)*PI;
		std::complex<double> ro_m, ro_p;
		std::complex<double> sum_p=0, sum_m=0;
		for (int l=0; l<=Lmax; l+=2)
		{
// calculating ro=potential!
			for(int k=-l; k<=l; k++)
			{
				ro_p=this->Get(l, k);
				ro_m=this->Get(l, k);

				for(int q=-l; q<=l; q++)
				{
					if(!isAxial)
						for(int p=-l; p<=l; p++)
						{
/*
							std::complex<double> add_p, add_m;

							add_p=this->Get(l, p)*ComplexWigner(l, k, q, fi, teta, alpha)*ComplexWigner(l, q, p, 0, -teta, -fi);
							add_m=this->Get(l, p)*ComplexWigner(l, k, q, fi, teta,-alpha)*ComplexWigner(l, q, p, 0, -teta, -fi);
							ro_p+=add_p;
							ro_m+=add_m;
*/
//							cout << add_p << " " << ro_p << " " << add_m << " " << ro_m << " " << endl;

							ro_p+=this->Get(l, p)*ComplexWigner(l, k, q, fi, teta, alpha)*ComplexWigner(l, q, p, 0, -teta, -fi);
							ro_m+=this->Get(l, p)*ComplexWigner(l, k, q, fi, teta,-alpha)*ComplexWigner(l, q, p, 0, -teta, -fi);
						}
					else
						{
							int p=0;
							ro_p+=this->Get(l, p)*ComplexWigner(l, k, q, fi, teta, alpha)*ComplexWigner(l, q, p, 0, -teta, -fi);
							ro_m+=this->Get(l, p)*ComplexWigner(l, k, q, fi, teta,-alpha)*ComplexWigner(l, q, p, 0, -teta, -fi);
						}
				}
				sum_p+=ro_p*(ComplexWigner(l, 0, k, 0, beta, gamma) );
				sum_m+=ro_m*(ComplexWigner(l, 0, k, 0, beta, gamma) );
			}
		}
// now ro becomes boltzman orientation distribution function:
		ro_p=exp(0.5*sum_p);
		ro_m=exp(0.5*sum_m);
//			cout << ro_p << " " << ro_m <<endl;
			for( int L=0; L<=2*MaxL; L+=2 )
			  for( int K=0; K<=L; K++ )
				{
					std::complex<double> cw=conj(ComplexWigner(L, 0, K, 0, beta, gamma));
					OP_JUMP_P[i][L*L/4+K] += cw*ro_p*KronWeight[ibeta]*KronWeight[igamma];	// integration step no *sin(beta) because cosine integration!
					OP_JUMP_M[i][L*L/4+K] += cw*ro_m*KronWeight[ibeta]*KronWeight[igamma];	// integration step

					OP_JUMP_P2[i][L*L/4+K]+= cw*ro_p*GausWeight[ibeta]*GausWeight[igamma];	// integration step no *sin(beta) because cosine integration!
					OP_JUMP_M2[i][L*L/4+K]+= cw*ro_m*GausWeight[ibeta]*GausWeight[igamma];	// integration step

				}

	  }
	}
	  
// normalization:
//	  cout << MaxL << endl;
	  for (int j=0; j<(2*MaxL+2)*(2*MaxL+2)/4; j++) 
	  {
		OP_JUMP_P [i][j]/=Normalization_OP;
		OP_JUMP_M [i][j]/=Normalization_OP;
		OP_JUMP_P2[i][j]/=Normalization_OP;
		OP_JUMP_M2[i][j]/=Normalization_OP;

	  }

// debug
/*
	int L0=2;
    for(int K0=-L0; K0<=L0; K0++)
	{
	  std::complex<double> sump=0;
	  for (int P=-L0; P<=L0; P++)
	  {
		  for(int Q=-L0; Q<=L0; Q++)
		  {
			  sump+=conj(ComplexWigner(L0, P, Q, fi, teta, -alpha)*ComplexWigner(L0, Q, K0, 0, -teta, -fi))*this->GetOP_JUMP_M(L0, P, i);
		  }
	  }
	  cout << L0 << " " << K0 << " " << sump << " " << this->GetOP_JUMP_P(L0, K0, i) << endl;
	}
*/
//	cout << "Finished!\n";
}
}
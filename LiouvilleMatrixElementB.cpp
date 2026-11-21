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

// Calculates Matrix elements of Liouville superoperator for biradical molecule
//

#include "SphericalTensor.h"
#include "LiouvilleMatrixElementB.h"

//clock_t time_spin;
//clock_t time_diffusion;
//clock_t time_potential;

extern double Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3);

// based on the spin part of the basis vector, calculates the index of a row/column, corresponding to this quantum numbers in the spin matrix.
int SpinIndexB(BasisVectorB &basis)
{
// TODO!!!
	double ns1, nss1, ms1, mss1;
	double ns2, nss2, ms2, mss2;

	ns1 =(basis.qS1+basis.pS1)/2.0;
    nss1=(basis.qS1-basis.pS1)/2.0;

    ms1 =(basis.qI1+basis.pI1)/2.0;
    mss1=(basis.qI1-basis.pI1)/2.0;

	ns2 =(basis.qS2+basis.pS2)/2.0;
    nss2=(basis.qS2-basis.pS2)/2.0;

    ms2 =(basis.qI2+basis.pI2)/2.0;
    mss2=(basis.qI2-basis.pI2)/2.0;

// I2' S2' I1' S1' I2'' S2'' I1'' S1''
//   2   3   2   3    2    3    2    3
//   1   3   6  18   36  108  216  648
	int q=(int)(    (1-ms2) +  3*(0.5- ns2)+  6*(1- ms1)+ 18*(0.5- ns1)
		        +36*(1-mss2)+108*(0.5-nss2)+216*(1-mss1)+648*(0.5-nss1) );
	return q;
	
//	return (int)(18*(0.5-nss1)+6*(1-mss1)+3*(0.5-ns1)+(1-ms1));       // this one is better for understanding!
//	return (int)(17.5-18*nss-6*mss-3*ns-ms);					  // this one is faster for computation
	
}

// M-symmetrized matrix element
std::complex<double> MatrixElementBSM(BasisVectorB basisr, BasisVectorB basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2)
{
	std::complex<double> a, b, c, d, temp1, temp2, temp3, temp4;
	bool invertM1=false, invertM2=false;

//	clock_t time_start=clock();
	if (abs(basisc.L-basisr.L)>2) 
	{
//		time_spin+=(clock()-time_start);
		return 0;
	}

	temp1=MatrixElementBSK1(basisr, basisc, Interactions,NInteractions, psi, J2);
	bool rselfconj=((basisr.M==0)&&(basisr.pI1==0));
	bool cselfconj=((basisc.M==0)&&(basisc.pI1==0));

	if (rselfconj) temp2=temp1;
	else
	{
		basisr.M=-basisr.M;
		basisr.pI1=-basisr.pI1;   
		temp2=MatrixElementBSK1(basisr, basisc, Interactions,NInteractions, psi, J2);
		basisr.M=-basisr.M;
		basisr.pI1=-basisr.pI1;
	}
	if (cselfconj)
		temp3=temp1;
	else
	{
		basisc.M =-basisc.M;
		basisc.pI1=-basisc.pI1;   
		temp3=MatrixElementBSK1(basisr, basisc, Interactions,NInteractions, psi, J2);
		basisc.M =-basisc.M;
		basisc.pI1=-basisc.pI1;
	}
	if (rselfconj) temp4=temp3;
	else if(cselfconj) temp4=temp2;
	else
	{
		basisc.M =-basisc.M;
		basisc.pI1=-basisc.pI1;   
		basisr.M=-basisr.M;
		basisr.pI1=-basisr.pI1;
		temp4=MatrixElementBSK1(basisr, basisc, Interactions,NInteractions, psi, J2);
		basisc.M =-basisc.M;
		basisc.pI1=-basisc.pI1;
		basisr.M=-basisr.M;
		basisr.pI1=-basisr.pI1;		
	}
	if (rselfconj)						a=0.5;
	else								a=0.7071067811865475;
	if (cselfconj)						c=0.5;
	else								c=0.7071067811865475;
	if ((basisr.L+basisr.M)%2)			b=-a;
	else							    b=a;
	if ((basisc.L+basisc.M)%2)			d=-c;
	else							    d=c;

//	cout << basisr.M << " " << basisc.M << " " << a << " " << b << " " << c << " " << d << "\n";
//	cout << temp1 << " " << temp2 << " " << temp3 << " " << temp4 << "\n";
//	time_spin+=(clock()-time_start);


	return a*c*temp1+b*c*temp2+a*d*temp3+b*d*temp4;


}

// calculation of small Wigner d-functions.
double wigner_precalculated_array[5][5];

void precalculate_wigner(double hi, double psi, double ksi)
{
	for(int m=-2; m<=2; m++)
		for(int ms=-2; ms<=2; ms++)
			wigner_precalculated_array[m+2][ms+2]=wigner_d(2, m, ms, psi);
}

double wigner_d_stored(int l, int m, int ms)
{
	if (l==0) return 1.0;
	else      return wigner_precalculated_array[m+2][ms+2];

}

// K-symmetrized matrix elements
std::complex<double> MatrixElementBSK1(BasisVectorB basisr, BasisVectorB basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2)
// Calculates matrix elements of spin part of Liouvillian, in K-symetrized space, with analytic formulae
{
	double jK1, jK2;
	std::complex<double> eK1, eK2;
	
//	double psi=0;
	int j, k;
	j=SpinIndexB(basisr);
	k=SpinIndexB(basisc);

	std::complex<double> temp=0;
	std::complex<double> I(0, 1);

// definition of jK, eK=sqrt(jK);
	if(basisr.K>0)
		jK1=1;
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else
	{
		basisr.K=-basisr.K;
		jK1=-1;
	}

	if(basisc.K>0)
		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;		
		else			jK2=1;
	}
    else
	{
		basisc.K=-basisc.K;
		jK2=-1;
	}

	if(jK1==1)		eK1=1;
	else			eK1=I;

	if(jK2==1)		eK2=1;
	else			eK2=I;

if(basisr.N==basisc.N)
{
	for (int l=0; l<=2; l+=2)
		for (int m=-l; m<=l; m++)
		{
			int ms=basisr.M-basisc.M;
			int mss=basisr.K-basisc.K;
			std::complex<double> temp3=Wigner3jMy(basisr.L, l, basisc.L, -basisr.M, ms, basisc.M);
			if(real(temp3)==0) continue;

			std::complex<double> temp1=0;
//			std::complex<double> temp2=wigner_d(l, m, ms, psi)*sqrt((double)(2*basisc.L+1)*(2*basisr.L+1));
			std::complex<double> temp2=wigner_d_stored(l, m, ms)*sqrt((double)(2*basisc.L+1)*(2*basisr.L+1));
			if(abs(temp2)==0) continue;

			std::complex<double> temp4=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K-basisc.K, basisc.K);
			std::complex<double> temp5=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K+basisc.K, -basisc.K)*jK2;
			if( (basisc.L+basisc.K)%2 ) temp5=-temp5;

			for (int i=0; i<NInteractions; i++)
			{
				if(jK1==jK2)
				{
// according to my formulae:
//					temp1=temp1+( temp4*2.0*real(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*real(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)))*Interactions[i].AuL.GetElement(l, m, j, k);
// according to my MIIF1982, (A42):
//					std::complex<double> temp6=( temp4*2.0*real(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*real(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)));

// if i==NInteractions-1, then the interaction is Lorentz linewidth, which has imaginary matrix elements
					if (i<NInteractions-1)	temp1=temp1+  ( temp4*2.0*real(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*real(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)))*conj(Interactions[i].AuL.GetElement(l, m, j, k));
					else					temp1=temp1+I*( temp4*2.0*real(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*real(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)))*conj(Interactions[i].AuL.GetElement(l, m, j, k));
//					cout << l << " " << m << " " << j << " " << k << " " << temp6 << endl;
				}
				else
				{
// according to my formulae:
//					temp1=temp1+( temp4*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)))*Interactions[i].AuL.GetElement(l, m, j, k);
// according to my MIIF1982, (A42):
//					std::complex<double> temp6=jK1*(temp4*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)));
					if (i<NInteractions-1)	temp1=temp1+  jK1*(temp4*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)))*conj(Interactions[i].AuL.GetElement(l, m, j, k));
					else					temp1=temp1+I*jK1*(temp4*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K-basisc.K))+temp5*2.0*imag(Interactions[i].Fdu.Get(l, basisr.K+basisc.K)))*conj(Interactions[i].AuL.GetElement(l, m, j, k));
//					cout << l << " " << m << " " << j << " " << k << " " << temp6 << endl;
				}
				
			}

//			if((basisr.L==0)&&(basisc.L==2)&&(basisr.M==basisc.M)&&(basisr.pI1==basisc.pI1)&&(basisr.qI1==basisc.qI1))
//			{
//				cout << l << " " << m << " " << ms << " " << wigner_d(l, m, ms, psi) << endl;
//				cout << temp1 << " " << temp2 << " " << temp3<< " " << temp4 << " " << temp5 << endl << endl;
//			}

			temp=temp+temp1*temp2*temp3;
		}

	if((basisr.K+basisr.M)%2) temp=-temp;

	if (basisr.K==0)	temp=temp*0.5;
	else				temp=temp*0.7071067811865475;

	if (basisc.K==0)	temp=temp*0.5;
	else				temp=temp*0.7071067811865475;

//	cout << basisr.K << " " << basisc.K << " " << jK1 << " " << jK2 << " " << eK1 << " " << eK2 << endl;
//	cout << endl;
}
else
{
	if ((basisr.L==basisc.L)&&
		(basisr.M==basisc.M)&&
		(basisr.K==basisc.K)&&
		(jK1==jK2))
	{
		int l=0, m=0;
		int ms=basisr.M-basisc.M;
		int mss=basisr.K-basisc.K;
		std::complex<double> temp3=Wigner3jMy(basisr.L, l, basisc.L, -basisr.M, ms, basisc.M);
		std::complex<double> temp1=0;
		std::complex<double> temp2=wigner_d_stored(l, m, ms)*sqrt((double)(2*basisc.L+1)*(2*basisr.L+1));
		std::complex<double> temp4=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K-basisc.K, basisc.K);
		std::complex<double> temp5=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K+basisc.K, -basisc.K)*jK2;
		if( (basisc.L+basisc.K)%2 ) temp5=-temp5;

		int i=4;   // (only J)
		if((basisr.N==0)&&(basisc.N==2) || (basisr.N==2)&&(basisc.N==0)) temp1=( temp4*2.0000000000000000000000000000000*J2 + temp5*2.0000000000000000000000000000000*J2 )*conj(Interactions[i].AuL.GetElement(l, m, j, k));
		else if (abs(basisr.N-basisc.N)<=2)								 temp1=( temp4*2.8284271247461900976033774484194*J2 + temp5*2.8284271247461900976033774484194*J2 )*conj(Interactions[i].AuL.GetElement(l, m, j, k));
		else															 temp1=0;

		temp=temp1*temp2*temp3;

		if((basisr.K+basisr.M)%2) temp=-temp;

		if (basisr.K==0)	temp=temp*0.5;
		else				temp=temp*0.7071067811865475;

		if (basisc.K==0)	temp=temp*0.5;
		else				temp=temp*0.7071067811865475;


	}


}
	return temp;


}

// auxiliary functions
double Nk(double K1, double K2)
{
	double n=1.0;
	if (K1==0) n *= 0.7071067811865475;
	if (K2==0) n *= 0.7071067811865475;
	return n;
}

double Np(double L, double K)
{
	return sqrt((L-K-1)*(L-K)*(L+K+1)*(L+K+2));
}
double Nm(double L, double K)
{
	return sqrt((L+K-1)*(L+K)*(L-K+1)*(L-K+2));
}
double Mp(double L, double K)
{
	return sqrt((L+K+1)*(L-K));
}
double Mm(double L, double K)
{
	return sqrt((L-K+1)*(L+K));
}
/*
std::complex<double> MatrixElementIsotropicBrownianB_ND_nosym(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R)
// subscipt _ND means that the diffusion tensor need not be diagonal;
{
	std::complex<double> element(0,0);
	std::complex<double> I(0,1);
// definition of jK, eK=sqrt(jK);
	if((basisr.L== basisc.L)&&
	   (basisr.M== basisc.M)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double L1=(double)basisr.L;
		double K1=((double)basisr.K);
		double K2=((double)basisc.K);
			if(K1==K2)
				element=(R(0,0)+R(1,1))*0.5*(L1*(L1+1)-K1*K1)+R(2,2)*K1*K1;
			else if(K1==K2+2)
				element=Np(L1, K2)*(R(0,0)-R(1,1)-2.0*I*R(0,1))*0.25;
			else if(K1==K2-2)
				element=Nm(L1, K2)*(R(0,0)-R(1,1)+2.0*I*R(0,1))*0.25;
			else if(K1==K2+1)
				element=(2.0*K2+1.0)*Mp(L1, K2)*(R(0,2)-I*R(1,2))*0.5;
			else if(K1==K2-1)
				element=(2.0*K2-1.0)*Mm(L1, K2)*(R(0,2)+I*R(1,2))*0.5;

	}
	return element;
}
*/

// THE RELAXATION PART OF THE LIOUVILLE SUPEROPERATOR:
//

// Isotropic Brownian diffusion, diffusion matrix non-diagonal (i. e., diffusion and orientation are non-coaxial).
std::complex<double> MatrixElementIsotropicBrownianB_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, double* Diffusion_Type, double* R_Tau)
// subscipt _ND means that the diffusion tensor need not be diagonal;
{
	std::complex<double> element(0,0);
// definition of jK, eK=sqrt(jK);
	double jK1, jK2;
	if(basisr.K>0)		jK1=1;
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;
	if((basisr.L== basisc.L)&&
	   (basisr.M== basisc.M)&&
	   (basisr.N  ==basisc.N)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double L1=(double)basisr.L;
		double K1=fabs((double)basisr.K);
		double K2=fabs((double)basisc.K);
/*
		if(jK1==jK2)
		{
			std::complex<double> element1, element2;
			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;
			element1=real(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));

			basisr.K=-basisr.K;
			element2=real(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));
			basisr.K=-basisr.K;

			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;

			if((basisr.L+basisr.K)%2) element2=-jK1*element2;
			else element2=jK1*element2;

			element=Nk(K1, K2)*(element1+element2);
		}
		else
		{
			std::complex<double> element1, element2;
			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;
			element1=imag(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));

			basisr.K=-basisr.K;
			element2=imag(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));
			basisr.K=-basisr.K;

			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;

			if((basisr.L+basisr.K)%2) element2=-jK1*element2;
			else element2=jK1*element2;

			element=-jK1*Nk(K1, K2)*(element1+element2);
		}
*/
// Non-Brownian Diffusion:
     if((Diffusion_Type[0]!=0)||(Diffusion_Type[1]!=0)||(Diffusion_Type[2]!=0))
	 {
		if(jK1==jK2)
		{
			if(K1==K2)
			{
				element=(R(0,0)+R(1,1))*0.5*(L1*(L1+1)*pow(1+R_Tau[0]*L1*(L1+1), -Diffusion_Type[0])-K1*K1*pow(1+R_Tau[1]*K1*K1, -Diffusion_Type[1]))+R(2,2)*K1*K1*pow(1+R_Tau[2]*K1*K1, -Diffusion_Type[2]);
			}
		}

	 }
// Brownian diffusion
	 else
	 {
		if(jK1==jK2)
		{
			if(K1==K2)
			{
				element=(R(0,0)+R(1,1))*0.5*(L1*(L1+1)-K1*K1)+R(2,2)*K1*K1;
				if(K1==1)    // особый случай K1=-1, K2=1!
				{
					if(basisr.L%2) element-=-jK1*Nm(L1, K2)*(R(0,0)-R(1,1))*0.25;
					else		   element+=-jK1*Nm(L1, K2)*(R(0,0)-R(1,1))*0.25;
				}
			}
			else if(K1==K2+2)
				element=Np(L1, K2)*(R(0,0)-R(1,1))*0.25/Nk(K1, K2);
			else if(K1==K2-2)
				element=Nm(L1, K2)*(R(0,0)-R(1,1))*0.25/Nk(K1, K2);
			else if(K1==K2+1)
				element=(2.0*K2+1.0)*Mp(L1, K2)*R(0,2)*0.5/Nk(K1, K2);
			else if(K1==K2-1)
				element=(2.0*K2-1.0)*Mm(L1, K2)*R(0,2)*0.5/Nk(K1, K2);
		}
		else	// if (jK1 != jK2)
		{
			if ((K1==1)&&(K2==1))
			{
				 element= Nm(L1, K2)*R(0,1)*0.5;
				 if(basisr.L%2) element=-element;
			}
			if     (K1==K2+2)
				element= jK1*Np(L1, K2)*R(0,1)*0.5/Nk(K1, K2);
			else if(K1==K2-2)
				element=-jK1*Nm(L1, K2)*R(0,1)*0.5/Nk(K1, K2);
			else if(K1==K2+1)
				element= jK1*(2.0*K2+1.0)*Mp(L1, K2)*R(1,2)*0.5/Nk(K1, K2);
			else if(K1==K2-1)
				element=-jK1*(2.0*K2-1.0)*Mm(L1, K2)*R(1,2)*0.5/Nk(K1, K2);

		}
	 }
	

	}

	return element;

}

// isotropic Brownian diffusion
std::complex<double> MatrixElementIsotropicBrownianB(BasisVectorB basisr, BasisVectorB basisc, double Rx, double Ry, double Rz)
{
	std::complex<double> element(0,0);
// definition of jK, eK=sqrt(jK);
	double jK1, jK2;
	if(basisr.K>0)		jK1=1;
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;


	if((basisr.L== basisc.L)&&
	   (basisr.N  ==basisc.N)&&
	   (basisr.M== basisc.M)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (jK1==jK2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double L1=(double)basisr.L;
		double K1=fabs((double)basisr.K);
		double K2=fabs((double)basisc.K);
		if(K1==K2)
			element=(Rx+Ry)*0.5*(L1*(L1+1)-K1*K1)+Rz*K1*K1;
		else if(K1-2==K2)
			element=Np(L1, K1-2)*(Rx-Ry)*0.25/Nk(K1, K1-2);
		else if(K1+2==K2)
			element=Nm(L1, K1+2)*(Rx-Ry)*0.25/Nk(K1, K1+2);
	
	}

	return element;
}


// XLK decomposition elements calculated in the case of non-diagonal rotation diffusion tensor.
std::complex<double> XLK_ND(int L, int K, Potential &pot, Matrix<std::complex<double>, 3, 3> &R)
// changed R matrix from complex<double> to double!
// for non-diagonal diffusion tensor!
{
	std::complex<double> result;
	std::complex<double> sum(0, 0);
	std::complex<double> I(0,1);
	if(L%2) return 0;
	if (abs(K)>L) return 0;
//	if(K%2) return 0;
//	if(K<0) K=-K;
	result=-0.5*pot.Get(L, K+2)*Nm((double)L, (double)(K+2))*0.25*(R(0,0)-R(1,1)+2.0*I*R(0,1));
	result-=0.5*pot.Get(L, K-2)*Np((double)L, (double)(K-2))*0.25*(R(0,0)-R(1,1)-2.0*I*R(0,1));
	result-=0.5*pot.Get(L, K+1)*Mm((double)L, (double)(K+1))*(2.0*(double)K+1)*0.5*(R(0,2)+I*R(1,2));
	result-=0.5*pot.Get(L, K-1)*Mp((double)L, (double)(K-1))*(2.0*(double)K-1)*0.5*(R(0,2)-I*R(1,2));
	result-=0.5*(pot.Get(L, K)*(L*(L+1)-K*K)*(R(0,0)+R(1,1))*0.5+pot.Get(L, K)*K*K*R(2,2));

//	printf("L %d K %d el=%lf\n", L, K, result);


	for(int L1=0; L1<=pot.Lmax; L1+=2)
		for(int L2=0; L2<=pot.Lmax; L2+=2)
			for(int K1=-L1; K1<=L1; K1+=2)
				for(int K2=-L2; K2<=L2; K2+=2)
				{
					std::complex<double> el;
					el=pot.Get(L1, K1)*pot.Get(L2, K2);
					el*=Wigner3jMy(L1, L, L2, 0, 0, 0);
					if(abs(el)==0) continue;
					double temp1, temp2, temp3, temp4, temp5, temp6;
					temp1=Mp((double)L1, (double)K1)*Mp((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2+1);
					temp2=Mm((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1-1, -K, K2-1);
					temp3=Mp((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2-1);
					temp4=												  K1*K2*Wigner3jMy(L1, L, L2, K1,   -K, K2);
					temp5=Mp((double)L1, (double)K1)*						 K2*Wigner3jMy(L1, L, L2, K1+1, -K, K2);
					temp6=Mm((double)L1, (double)K1)*						 K2*Wigner3jMy(L1, L, L2, K1-1, -K, K2);
					el*= 0.25*(R(0,0)-R(1,1)-2.0*I*R(0,1))*temp1
						+0.25*(R(0,0)-R(1,1)+2.0*I*R(0,1))*temp2
						+0.5 *(R(0,0)+R(1,1))			*temp3
						+R(2,2)*temp4
						+(R(0,2)-I*R(1,2))*temp5
						+(R(0,2)+I*R(1,2))*temp6;
//					cout << Wigner3jMy(2, 2, 2, 1, -2, -1) << endl;
//					printf("L %d K %d L1 %d L2 %d K1 %d K2 %d el=%lf\n", L, K, L1, L2, K1, K2, el);
					sum+=el;
				}
	if(L%2) return result+sum*0.25*(2.0*L+1);
	else	return result-sum*0.25*(2.0*L+1);


}

// XLK decomposition components in the case of diagonal rotation diffusion tensor.
double XLK(int L, int K, Potential &pot, double Rx, double Ry, double Rz)
{
	double result, sum=0;
	if(L%2) return 0;
	if(K%2) return 0;
//	if(K<0) K=-K;
	result=-0.5*(pot.Get(L, K+2)*Nm((double)L, (double)(K+2))+pot.Get(L, K-2)*Np((double) L, double(K-2)))*0.25*(Rx-Ry);
	result-=0.5*(pot.Get(L, K)*(L*(L+1)-K*K)*(Rx+Ry)*0.5+pot.Get(L, K)*K*K*Rz);

//	printf("L %d K %d el=%lf\n", L, K, result);


	for(int L1=0; L1<=pot.Lmax; L1+=2)
		for(int L2=0; L2<=pot.Lmax; L2+=2)
			for(int K1=-L1; K1<=L1; K1+=2)
				for(int K2=-L2; K2<=L2; K2+=2)
				{
					double el=pot.Get(L1, K1)*pot.Get(L2, K2);
					el*=Wigner3jMy(L1, L, L2, 0, 0, 0);
					if(el==0) continue;
					double temp1, temp2, temp3, temp4;
					temp1=Mp((double)L1, (double)K1)*Mp((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2+1);
					temp2=Mm((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1-1, -K, K2-1);
					temp3=Mp((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2-1);
					temp4=												  K1*K2*Wigner3jMy(L1, L, L2, K1, -K, K2);

					el*=0.25*(Rx-Ry)*(temp1+temp2)+0.5*(Rx+Ry)*temp3+Rz*temp4;
//					cout << Wigner3jMy(2, 2, 2, 1, -2, -1) << endl;
//					printf("L %d K %d L1 %d L2 %d K1 %d K2 %d el=%lf\n", L, K, L1, L2, K1, K2, el);
					sum+=el;
				}
	return result-sum*0.25*(2.0*L+1);
}	

double Nl(double L1, double L2)
{
	return sqrt((2*L1+1)*(2*L2+1));
}

// Matrix element of Brownian diffusion in anisotropic medium, in the case of non-coaxial orientation and diffusion:
std::complex<double> MatrixElementPotentialB_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, Potential &pot)
// for non-diagonal diffusion tensors!
{
	std::complex<double> element(0,0);
// definition of jK, eK=sqrt(jK);
	double jK1, jK2;
	if(basisr.K>0)		jK1=1;
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;


	if((basisr.M  ==basisc.M)&&
	   (basisr.N  ==basisc.N)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double sum=0;
		int L1, L2, K1, K2, M1, M2;
		L1=basisr.L;
		L2=basisc.L;
		K1=basisr.K;
		K2=basisc.K;
		M1=basisr.M;
		M2=basisc.M;
		double temp1, temp2, temp3;

		for(int L=0; L<=2*pot.Lmax; L+=2)
		{
			temp1=Wigner3jMy(L1, L, L2, M1, 0, -M2);
			if (temp1==0) continue;
			if(jK1==jK2)
			{
//				temp2=real(XLK_ND(L, K1-K2, pot, R))*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
//				temp3=real(XLK_ND(L, K1+K2, pot, R))*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
				temp2=real(pot.GetXLK(L, K1-K2, 0))*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
				temp3=real(pot.GetXLK(L, K1+K2, 0))*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
				if((L2+K2)%2) temp3=-temp3;
			}
			else
			{
//				temp2=imag(XLK_ND(L, K1-K2, pot, R))*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
//				temp3=imag(XLK_ND(L, K1+K2, pot, R))*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
				temp2=imag(pot.GetXLK(L, K1-K2, 0))*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
				temp3=imag(pot.GetXLK(L, K1+K2, 0))*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
				if((L2+K2)%2) temp3=-temp3;
			}
			sum+=temp1*(temp2+temp3);
		}
		sum*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum=-sum;
		if(jK1!=jK2) sum=jK2*sum;
		element=sum;
	}

	return element;


}

// The case of coaxial diffusion and orientation (according to Freed at al.)
std::complex<double> MatrixElementPotentialB(BasisVectorB basisr, BasisVectorB basisc, double Rx, double Ry, double Rz, Potential &pot)
{
	std::complex<double> element(0,0);
// definition of jK, eK=sqrt(jK);
	double jK1, jK2;
	if(basisr.K>0)		jK1=1;
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;


	if((basisr.M== basisc.M)&&
	   (basisr.N  ==basisc.N)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (jK1==jK2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double sum=0;
		int L1, L2, K1, K2, M1, M2;
		L1=basisr.L;
		L2=basisc.L;
		K1=basisr.K;
		K2=basisc.K;
		M1=basisr.M;
		M2=basisc.M;

		for(int L=0; L<=2*pot.Lmax; L+=2)
		{
			double temp1, temp2, temp3;
			temp1=Wigner3jMy(L1, L, L2, M1, 0, -M1);
			if (temp1==0) continue;
			temp2=XLK(L, K1-K2, pot, Rx, Ry, Rz)*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
			temp3=XLK(L, K1+K2, pot, Rx, Ry, Rz)*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
			if((L2+K2)%2) temp3=-temp3;

			sum+=temp1*(temp2+temp3);
		}
		sum*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum=-sum;
		element=sum;
	}

	return element;
}



// M-Symmetrized, for High-field approximation
std::complex<double> MatrixElementPotentialBSM_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, Potential &pot)
// for non-diagonal diffusion tensors!
{
	std::complex<double> element(0,0);
// definition of jK, eK=sqrt(jK);
	double jK1, jK2;
	double jM1, jM2;
//	clock_t time_start=clock();

	if(basisr.K>0)		jK1=1;
//definition of jK1, jK2
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;
//definition of jM1, jM2
	if(basisr.M>0)		jM1=1;
	else if(basisr.M==0)
	{
		if (basisr.L%2) jM1=-1;
		else			jM1=1;
	}
	else				jM1=-1;

	if(basisc.M>0)		jM2=1;
	else if (basisc.M==0)
	{
		if (basisc.L%2) jM2=-1;
		else			jM2=1;
	}
    else				jM2=-1;

	if((basisr.M  ==basisc.M)&&
	   (basisr.N  ==basisc.N)&&
	   (jM1==jM2)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double sum=0;
		int L1, L2, K1, K2, M1, M2;
		L1=basisr.L;
		L2=basisc.L;
		K1=basisr.K;
		K2=basisc.K;
		M1=basisr.M;
		M2=basisc.M;
		double temp1, temp2, temp3;

		for(int L=0; L<=2*pot.Lmax; L+=2)
		{
			temp1=Wigner3jMy(L1, L, L2, M1, 0, -M2);
			if (temp1==0) continue;
			if(jK1==jK2)
			{
				temp2=real(pot.GetXLK(L, K1-K2, 0))*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
				temp3=real(pot.GetXLK(L, K1+K2, 0))*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
				if((L2+K2)%2) temp3=-temp3;
			}
			else
			{
				temp2=imag(pot.GetXLK(L, K1-K2, 0))*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
				temp3=imag(pot.GetXLK(L, K1+K2, 0))*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
				if((L2+K2)%2) temp3=-temp3;
			}
			sum+=temp1*(temp2+temp3);
		}
		sum*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum=-sum;
		if(jK1!=jK2) sum=jK2*sum;
		element=sum;
	}

//	time_potential+=(clock()-time_start);
	return element;


}




std::complex<double> MatrixElementIsotropicBrownianBSM_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, double *Diffusion_Type, double *R_Tau)
// subscipt _ND means that the diffusion tensor need not be diagonal;
// _BSM means that this function is for M-symmetrized basis
{
	std::complex<double> element(0,0);

//	clock_t time_start=clock();


// definition of jK, eK=sqrt(jK);
	double jK1, jK2, jM1, jM2;
	if(basisr.K>0)		jK1=1;
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;

//definition of jM1, jM2
	if(basisr.M>0)		jM1=1;
	else if(basisr.M==0)
	{
		if (basisr.L%2) jM1=-1;
		else			jM1=1;
	}
	else				jM1=-1;

	if(basisc.M>0)		jM2=1;
	else if (basisc.M==0)
	{
		if (basisc.L%2) jM2=-1;
		else			jM2=1;
	}
    else				jM2=-1;


	if((basisr.L== basisc.L)&&
	   (basisr.M== basisc.M)&&
	   (basisr.N  ==basisc.N)&&
	   (jM1==jM2)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double L1=(double)basisr.L;
		double K1=fabs((double)basisr.K);
		double K2=fabs((double)basisc.K);
/*
		if(jK1==jK2)
		{
			std::complex<double> element1, element2;
			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;
			element1=real(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));

			basisr.K=-basisr.K;
			element2=real(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));
			basisr.K=-basisr.K;

			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;

			if((basisr.L+basisr.K)%2) element2=-jK1*element2;
			else element2=jK1*element2;

			element=Nk(K1, K2)*(element1+element2);
		}
		else
		{
			std::complex<double> element1, element2;
			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;
			element1=imag(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));

			basisr.K=-basisr.K;
			element2=imag(MatrixElementIsotropicBrownianB_ND_nosym(basisr, basisc, R));
			basisr.K=-basisr.K;

			if(jK1<0) basisr.K=-basisr.K;
			if(jK2<0) basisc.K=-basisc.K;

			if((basisr.L+basisr.K)%2) element2=-jK1*element2;
			else element2=jK1*element2;

			element=-jK1*Nk(K1, K2)*(element1+element2);
		}
*/
// Non-Brownian Diffusion:
     if((Diffusion_Type[0]!=0)||(Diffusion_Type[1]!=0)||(Diffusion_Type[2]!=0))
	 {
		if(jK1==jK2)
		{
			if(K1==K2)
			{
				element=(R(0,0)+R(1,1))*0.5*(L1*(L1+1)*pow(1+R_Tau[0]*L1*(L1+1), -Diffusion_Type[0])-K1*K1*pow(1+R_Tau[1]*K1*K1, -Diffusion_Type[1]))+R(2,2)*K1*K1*pow(1+R_Tau[2]*K1*K1, -Diffusion_Type[2]);
			}
		}

	 }
// Brownian diffusion
	 else
	 {
		if(jK1==jK2)
		{
			if(K1==K2)
			{
				element=(R(0,0)+R(1,1))*0.5*(L1*(L1+1)-K1*K1)+R(2,2)*K1*K1;
				if(K1==1)    // особый случай K1=-1, K2=1!
				{
					if(basisr.L%2) element-=-jK1*Nm(L1, K2)*(R(0,0)-R(1,1))*0.25;
					else		   element+=-jK1*Nm(L1, K2)*(R(0,0)-R(1,1))*0.25;
				}
			}
			else if(K1==K2+2)
				element=Np(L1, K2)*(R(0,0)-R(1,1))*0.25/Nk(K1, K2);
			else if(K1==K2-2)
				element=Nm(L1, K2)*(R(0,0)-R(1,1))*0.25/Nk(K1, K2);
			else if(K1==K2+1)
				element=(2.0*K2+1.0)*Mp(L1, K2)*R(0,2)*0.5/Nk(K1, K2);
			else if(K1==K2-1)
				element=(2.0*K2-1.0)*Mm(L1, K2)*R(0,2)*0.5/Nk(K1, K2);
		}
		else	// if (jK1 != jK2)
		{
			if ((K1==1)&&(K2==1))
			{
				 element= Nm(L1, K2)*R(0,1)*0.5;
				 if(basisr.L%2) element=-element;
			}
			if     (K1==K2+2)
				element= jK1*Np(L1, K2)*R(0,1)*0.5/Nk(K1, K2);
			else if(K1==K2-2)
				element=-jK1*Nm(L1, K2)*R(0,1)*0.5/Nk(K1, K2);
			else if(K1==K2+1)
				element= jK1*(2.0*K2+1.0)*Mp(L1, K2)*R(1,2)*0.5/Nk(K1, K2);
			else if(K1==K2-1)
				element=-jK1*(2.0*K2-1.0)*Mm(L1, K2)*R(1,2)*0.5/Nk(K1, K2);

		}

	
	 }
	}
//	time_diffusion+=(clock()-time_start);

	return element;

}


/*

std::complex<double> MatrixElementDiscreteJumpsB(BasisVectorB basisr, BasisVectorB basisc, int NSites, double JumpFrequency)
// subscipt _ND means that the diffusion tensor need not be diagonal;
// _BSM means that this function is for M-symmetrized basis
{
	std::complex<double> element(0,0);

	if((basisr.L== basisc.L)&&
	   (basisr.M== basisc.M)&&
	   (basisr.K== basisc.K)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		if (basisr.K % NSites) element=JumpFrequency;
	}

	return element;

}
*/



std::complex<double> MatrixElementJumpsBSM_Pot(BasisVectorB basisr, BasisVectorB basisc, MagneticParameters* MagneticParameter, Potential &pot)
// for non-diagonal diffusion tensors!
{
	std::complex<double> element(0,0);
	std::complex<double> I(0, 1.0);
// definition of jK, eK=sqrt(jK);
	double jK1, jK2;
	double jM1, jM2;
	double sum=0;
	double sum2=0;

//	clock_t time_start=clock();

	if(basisr.K>0)		jK1=1;
//definition of jK1, jK2
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;
//definition of jM1, jM2
	if(basisr.M>0)		jM1=1;
	else if(basisr.M==0)
	{
		if (basisr.L%2) jM1=-1;
		else			jM1=1;
	}
	else				jM1=-1;

	if(basisc.M>0)		jM2=1;
	else if (basisc.M==0)
	{
		if (basisc.L%2) jM2=-1;
		else			jM2=1;
	}
    else				jM2=-1;

	if((basisr.M  ==basisc.M)&&
	   (basisr.N  ==basisc.N)&&
	   (jM1==jM2)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		int L1, L2, K1, K2, M1, M2;
		L1=basisr.L;
		L2=basisc.L;
		K1=basisr.K;
		K2=basisc.K;
		M1=basisr.M;
		M2=basisc.M;
		double temp1, temp2, temp3, temp4, temp5;

		for(int L=abs(L1-L2); L<=L1+L2; L++)
		{
			temp1=Wigner3jMy(L1, L, L2, -M1, 0, M1);
			if (temp1==0) continue;
			temp1 *= (2*L+1);
			if(jK1==jK2)
			{

			for(int i=0; i<MagneticParameter->NJumpAxes; i++)
			{
				sum=0; 
				if ((MagneticParameter->JumpAxisFi[i]==0)&&(MagneticParameter->JumpAxisTeta[i]==0))
				{
					temp4=Wigner3jMy(L1, L, L2, -K1, K1-K2, K2);
					temp5=Wigner3jMy(L1, L, L2, -K1, K1+K2,-K2);

					temp2=2*cos((K2-K1)*MagneticParameter->JumpAngle[i])*pot.GetOP(L, K1-K2)
						  -real(exp(-I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1-K2, i)+exp(I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1-K2, i));

					temp3=2*cos((K2+K1)*MagneticParameter->JumpAngle[i])*pot.GetOP(L, K1+K2)
						  -real(exp(I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1+K2, i)+exp(-I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1+K2, i));

					if((L2+K2)%2) temp3=-temp3;
					temp3 *= jK2;

					sum+=temp2*temp4 + temp3*temp5;

				}
				else
				{

				for (int P=-L; P<=L; P++)
				{
					temp4=Wigner3jMy(L1, L, L2, -K1, K1-K2, K2)*pot.GetOP(L, -P);
					temp5=Wigner3jMy(L1, L, L2, -K1, K1+K2,-K2)*pot.GetOP(L, -P);
					if ((temp4==0)&&(temp5==0)) continue;

					for (int Q=-L; Q<=L; Q++)
					{
						temp2=2.0*cos(Q*MagneticParameter->JumpAngle[i])*
											real(ComplexWigner(L, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L, Q, K2-K1, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												);
						if((P+K1+K2)%2) temp2=-temp2;

						temp3=2.0*cos(Q*MagneticParameter->JumpAngle[i])*
											real(ComplexWigner(L, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L, Q,-K2-K1, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												);
						if((P+K1+L2)%2) temp3=-temp3;
						temp3 *= jK2;

						sum+=temp2*temp4 + temp3*temp5;
					}
				}
				for (int P=-L2; P<=L2; P++)
				{
					temp4=Wigner3jMy(L1, L, L2, -K1, K1-P, P);
					if(temp4==0) continue;
					for (int Q=-L2; Q<=L2; Q++)
					{
						temp2=real((exp(-I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1-P, i)+exp(I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1-P, i))
												*ComplexWigner(L2, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L2, Q, K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												);

						temp3=real((exp(-I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1-P, i)+exp(I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1-P, i))
												*ComplexWigner(L2, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L2, Q,-K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												);
						if((L2+K2)%2) temp3=-temp3;
						temp3 *= jK2;

						sum-=temp2*temp4 + temp3*temp4;
					}
				}
				}
				sum*=MagneticParameter->JumpFrequencyR[i];
				sum2+=sum*temp1;

			}

			}
			else // if(jK1 != jK2)
			{
			for(int i=0; i<MagneticParameter->NJumpAxes; i++)
			{
				temp4=Wigner3jMy(L1, L, L2, -K1, K1-K2, K2);
				temp5=Wigner3jMy(L1, L, L2, -K1, K1+K2,-K2);

				sum=0; 

				if ((MagneticParameter->JumpAxisFi[i]==0)&&(MagneticParameter->JumpAxisTeta[i]==0))
				{
					temp4=Wigner3jMy(L1, L, L2, -K1, K1-K2, K2);
					temp5=Wigner3jMy(L1, L, L2, -K1, K1+K2,-K2);

					temp2=-imag(exp(-I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1-K2, i)+exp( I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1-K2, i));

					temp3=-imag(exp( I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1+K2, i)+exp(-I*(double)K2*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1+K2, i));

					if((L2+K2)%2) temp3=-temp3;
					temp3 *= jK2;

					sum+=temp2*temp4 + temp3*temp5;
				}
				else
				{
				for (int P=-L; P<=L; P++)
					for (int Q=-L; Q<=L; Q++)
					{
						temp2=2.0*cos(Q*MagneticParameter->JumpAngle[i])*
											imag(ComplexWigner(L, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L, Q, K2-K1, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												*pot.GetOP(L, -P));
						if((P+K1+K2)%2) temp2=-temp2;

						temp3=2.0*cos(Q*MagneticParameter->JumpAngle[i])*
											imag(ComplexWigner(L, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L, Q,-K2-K1, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												*pot.GetOP(L, -P));
						if((P+K1+L2)%2) temp3=-temp3;
						temp3 *= jK2;

						sum+=temp2*temp4 + temp3*temp5;
					}
				for (int P=-L2; P<=L2; P++)
				{
					temp4=Wigner3jMy(L1, L, L2, -K1, K1-P, P);
				  if(temp4==0) continue;
					for (int Q=-L2; Q<=L2; Q++)
					{
						temp2=imag((exp(-I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1-P, i)+exp(I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1-P, i))
												*ComplexWigner(L2, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L2, Q, K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												);

						temp3=imag((exp(-I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_P(L, K1-P, i)+exp(I*(double)Q*MagneticParameter->JumpAngle[i])*pot.GetOP_JUMP_M(L, K1-P, i))
												*ComplexWigner(L2, P, Q,     MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)
											    *ComplexWigner(L2, Q,-K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i])
												);

						if((L2+K2)%2) temp3=-temp3;
						temp3 *= jK2;

						sum-=temp2*temp4 + temp3*temp4;

					}
				}

				}

				sum*=MagneticParameter->JumpFrequencyR[i];
				sum2+=sum*temp1;

			}



			}
//			sum2+=sum*temp1;
		}
		sum2*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum2=-sum2;
		if(jK1!=jK2) sum2=jK2*sum2;
		element=sum2;
	}

//	time_potential+=(clock()-time_start);
	return element;


}


std::complex<double> MatrixElementJumpsBSM_Iso(BasisVectorB basisr, BasisVectorB basisc, MagneticParameters *MagneticParameter)
{
	std::complex<double> element(0,0);
// definition of jK, eK=sqrt(jK);
	double jK1, jK2;
	double jM1, jM2;
//	clock_t time_start=clock();

	if(basisr.K>0)		jK1=1;
//definition of jK1, jK2
	else if(basisr.K==0)
	{
		if (basisr.L%2) jK1=-1;
		else			jK1=1;
	}
	else				jK1=-1;

	if(basisc.K>0)		jK2=1;
	else if (basisc.K==0)
	{
		if (basisc.L%2) jK2=-1;
		else			jK2=1;
	}
    else				jK2=-1;
//definition of jM1, jM2
	if(basisr.M>0)		jM1=1;
	else if(basisr.M==0)
	{
		if (basisr.L%2) jM1=-1;
		else			jM1=1;
	}
	else				jM1=-1;

	if(basisc.M>0)		jM2=1;
	else if (basisc.M==0)
	{
		if (basisc.L%2) jM2=-1;
		else			jM2=1;
	}
    else				jM2=-1;

	if((basisr.L  ==basisc.L)&&
	   (basisr.M  ==basisc.M)&&
	   (basisr.N  ==basisc.N)&&
	   (jM1==jM2)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		double sum=0;
		int L1, L2, K1, K2, M1, M2;
		L1=basisr.L;
		L2=basisc.L;
		K1=basisr.K;
		K2=basisc.K;
		M1=basisr.M;
		M2=basisc.M;
		double temp1, temp2, temp3;
		

		if(jK1==jK2)
		{

			for(int i=0; i<MagneticParameter->NJumpAxes; i++)
			{
				temp1=0;
				temp2=0;
				temp3=0;

				if(K1==K2)  temp1+=1;
				if(K1==-K2) 
				{
					if((L2+K2)%2)	temp1+=-jK2;
					else			temp1+= jK2;
				}
				for (int Q=-L2; Q<=L2; Q++)
				{
					temp2+=cos((double)Q*MagneticParameter->JumpAngle[i])*
											real(ComplexWigner(L2, K1, Q, MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)*ComplexWigner(L2, Q, K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i]));
					temp3+=cos((double)Q*MagneticParameter->JumpAngle[i])*
											real(ComplexWigner(L2, K1, Q, MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)*ComplexWigner(L2, Q,-K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i]));
				}
				if((L2+K2)%2)	temp3*=-jK2;
				else			temp3*= jK2;
				sum+=(temp1-temp2-temp3)*MagneticParameter->JumpFrequencyR[i];

				}

		}
		else
		{
			for(int i=0; i<MagneticParameter->NJumpAxes; i++)
			{
				temp1=0;
				temp2=0;
				temp3=0;

				for (int Q=-L2; Q<=L2; Q++)
				{
					temp2+=cos((double)Q*MagneticParameter->JumpAngle[i])*
											imag(ComplexWigner(L2, K1, Q, MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)*ComplexWigner(L2, Q, K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i]));
					temp3+=cos((double)Q*MagneticParameter->JumpAngle[i])*
											imag(ComplexWigner(L2, K1, Q, MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)*ComplexWigner(L2, Q,-K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i]));
				}
				if((L2+K2)%2)	temp3*=-jK2;
				else			temp3*= jK2;
				sum+=(-temp2-temp3)*MagneticParameter->JumpFrequencyR[i];

			}
		}


			
		sum*=Nk((double)K1, (double)K2);
		if(jK1!=jK2) sum=jK2*sum;
		element=2.0*sum;
	}

//	time_potential+=(clock()-time_start);
	return element;


}


std::complex<double> MatrixElementConformation(BasisVectorB basisr, BasisVectorB basisc, double* ConversionRate)
{
	std::complex<double> element(0,0);
	if((basisr.L  ==basisc.L)&&
	   (basisr.M  ==basisc.M)&&
	   (basisr.K  ==basisc.K)&&
	   (basisr.N  ==basisc.N)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		element=ConversionRate[MAX_PARTICLES*basisr.Q + basisc.Q];
	}

//	time_potential+=(clock()-time_start);
	return element;

}

std::complex<double> MatrixElementJDiffusion(BasisVectorB basisr, BasisVectorB basisc, double JDiffusionRate)
{
	std::complex<double> element(0,0);
	if((basisr.L  ==basisc.L)&&
	   (basisr.M  ==basisc.M)&&
	   (basisr.K  ==basisc.K)&&
	   (basisr.N  ==basisc.N)&&
	   (basisr.pS1==basisc.pS1)&&
	   (basisr.qS1==basisc.qS1)&&
	   (basisr.pI1==basisc.pI1)&&
	   (basisr.qI1==basisc.qI1)&&
	   (basisr.pS2==basisc.pS2)&&
	   (basisr.qS2==basisc.qS2)&&
	   (basisr.pI2==basisc.pI2)&&
	   (basisr.qI2==basisc.qI2))
	{
		element=JDiffusionRate*basisc.N*basisc.N;
	}

//	time_potential+=(clock()-time_start);
	return element;

}

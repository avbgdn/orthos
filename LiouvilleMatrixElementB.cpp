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
// Overview: every function below returns one element <row|operator|column> between two basis vectors of the Liouville space
// used by the stochastic Liouville equation (SLE). A basis vector (BasisVectorB, see BasisB.h) is |L,K,M,N,Q,pS1,qS1,pI1,qI1,pS2,qS2,pI2,qI2>:
// L,K,M label the orientational (Wigner function) part, N is an extra index coupled by the exchange term J2 and damped by J-diffusion,
// Q is the conformer index, and (p,q) pairs label the spin part for electrons S1,S2 and nuclei I1,I2.
// The spin-Hamiltonian elements are returned without the factor -i; the caller (CalculateLiouvillianB) multiplies by -i and adds the
// relaxation terms (rotational diffusion, diffusion in the ordering potential, jumps, conformer exchange, J-diffusion).
// Convention: jK=+1 or -1 labels which of the two K-symmetrised combinations of the +K and -K orientational functions a vector is.
// It is stored in the sign of K (K<0 means jK=-1; the functions below work with |K|) and, for K=0, in the parity of L (odd L means jK=-1).
// In the M-symmetrised monoradical basis (M>=0 only) jM is defined in the same way from M.
// Functions with BSM in the name are for the M-symmetrised monoradical basis; functions without it are for the unsymmetrised-in-M basis.
// Used by: LiouvillianB.cpp (CalculateLiouvillianB, which assembles the sparse Liouvillian from these elements) and Potential.cpp (XLK_ND).
// The declarations are in LiouvilleMatrixElementB.h.

#include "SphericalTensor.h"
#include "LiouvilleMatrixElementB.h"

//clock_t time_spin;
//clock_t time_diffusion;
//clock_t time_potential;

extern double Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3);

// based on the spin part of the basis vector, calculates the index of a row/column, corresponding to this quantum numbers in the spin matrix.
// Maps the spin part of a basis vector to one running index of the spin Liouville space (0...1295). The p,q labels give the quantum numbers
// of the two sides (primed = ket, double primed = bra) of the spin operator: n'=(q+p)/2 and n''=(q-p)/2 for each of S1, I1, S2, I2.
// The packing is hard-coded for electron spin 1/2 (2 states) and nuclear spin 1 (3 states); only the spin indices of the basis vector are used.
// The index addresses the spin-interaction matrices (AuL.GetElement) and the starting vector.
// Used by: MatrixElementBSK1 (this file); CalculateLiouvillianB (LiouvillianB.cpp, construction of the starting vector).
int SpinIndexB(BasisVectorB &basis)
{
// TODO!!!
	double ns1, nss1, ms1, mss1;
	double ns2, nss2, ms2, mss2;

	// ket and bra quantum numbers of electron and nuclear spins from the p,q labels
	ns1 =(basis.qS1+basis.pS1)/2.0;
    nss1=(basis.qS1-basis.pS1)/2.0;

    ms1 =(basis.qI1+basis.pI1)/2.0;
    mss1=(basis.qI1-basis.pI1)/2.0;

	ns2 =(basis.qS2+basis.pS2)/2.0;
    nss2=(basis.qS2-basis.pS2)/2.0;

    ms2 =(basis.qI2+basis.pI2)/2.0;
    mss2=(basis.qI2-basis.pI2)/2.0;

// mixed-radix running index: the table gives the number of states and the weight of each quantum number
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
// M-symmetrised spin-Hamiltonian element (without the factor -i) for the monoradical basis: the M-symmetrised combination
// of four MatrixElementBSK1 elements. The conjugation operator S reverses M and pI1 (qS1 is zero for every vector of the
// monoradical basis). The symmetrised vector is a*(|r>+eta*|Sr>) with eta=(-1)^(L+M); a vector with S r=r (M=0 and pI1=0)
// has a=1/2, every other vector a=1/sqrt(2). Arguments as for MatrixElementBSK1.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) for monoradical conformers.
std::complex<double> MatrixElementBSM(BasisVectorB basisr, BasisVectorB basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2)
{
	std::complex<double> a, b, c, d, temp1, temp2, temp3, temp4;
	bool invertM1=false, invertM2=false;

//	clock_t time_start=clock();
	// the spin interactions have rank 0 and 2, so L can change by at most 2
	if (abs(basisc.L-basisr.L)>2) 
	{
//		time_spin+=(clock()-time_start);
		return 0;
	}

	// unsymmetrised element between the row and column vectors
	temp1=MatrixElementBSK1(basisr, basisc, Interactions,NInteractions, psi, J2);
	// a vector is self-conjugate (S r=r) when M and pI1 are both zero
	bool rselfconj=((basisr.M==0)&&(basisr.pI1==0));
	bool cselfconj=((basisc.M==0)&&(basisc.pI1==0));

	// element with the row vector replaced by its conjugate partner S r, shared with temp1 if the row vector is self-conjugate
	if (rselfconj) temp2=temp1;
	else
	{
		basisr.M=-basisr.M;
		basisr.pI1=-basisr.pI1;   
		temp2=MatrixElementBSK1(basisr, basisc, Interactions,NInteractions, psi, J2);
		basisr.M=-basisr.M;
		basisr.pI1=-basisr.pI1;
	}

	// element with the column vector replaced by its conjugate partner S c, shared with temp1 if the column vector is self-conjugate
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

	// element with both vectors conjugated; coincides with temp3 or temp2 if the row or the column vector is self-conjugate
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

	// amplitude a (c) of the row (column) vector, 1/2 if it is self-conjugate and 1/sqrt(2) otherwise, and coefficient b (d)
	// of its conjugate partner, with the (-1)^(L+M) phase rule
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


	// bilinear combination of the four elements
	return a*c*temp1+b*c*temp2+a*d*temp3+b*d*temp4;


}

// calculation of small Wigner d-functions.
// Table of the rank-2 reduced Wigner functions d(2,m,ms;psi), m,ms=-2..2, stored as [m+2][ms+2] for the current orientation angle psi.
// A single shared global: filled by precalculate_wigner, read by wigner_d_stored, valid for the psi of the latest precalculate_wigner call.
double wigner_precalculated_array[5][5];

// Fills wigner_precalculated_array for the angle psi (the orientation angle passed to CalculateLiouvillianB; psi=0 is the DirectorParallel case),
// so that MatrixElementBSK1 does not recompute the Wigner functions for every matrix element. hi and ksi are not used.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp), once per orientation before the matrix is built.
void precalculate_wigner(double hi, double psi, double ksi)
{
	// all rank-2 components (m,ms)
	for(int m=-2; m<=2; m++)
		for(int ms=-2; ms<=2; ms++)
			wigner_precalculated_array[m+2][ms+2]=wigner_d(2, m, ms, psi);
}

// Returns d(l,m,ms;psi) for the psi of the last precalculate_wigner call: 1 for l=0, the stored rank-2 table otherwise (so valid only for l=0,2).
// Used by: MatrixElementBSK1 (this file).
double wigner_d_stored(int l, int m, int ms)
{
	// rank 0 is trivially 1, rank 2 is read from the table
	if (l==0) return 1.0;
	else      return wigner_precalculated_array[m+2][ms+2];

}

// K-symmetrized matrix elements
// Spin part of the Liouvillian (the caller multiplies it by -i) between row vector basisr and column vector basisc in the K-symmetrised basis.
// It is a sum over the spin interactions Interactions[0..NInteractions-1] (g1, g2, A1, A2, DJ, Lorentz linewidth) of products of: the
// orientational overlap of the two basis functions with the rank-l (l=0,2) interaction tensor (Wigner 3j symbols in L,M,K and the rotation
// d(l,m,ms;psi)), the molecular-frame tensor components Fdu, and the spin-space matrix element AuL, addressed through SpinIndexB.
// Vectors with different N are coupled only by the J2 term (Interactions[4]). psi is the orientation angle, J2 the exchange coupling.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) for biradical conformers; MatrixElementBSM (this file).
std::complex<double> MatrixElementBSK1(BasisVectorB basisr, BasisVectorB basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2)
// Calculates matrix elements of spin part of Liouvillian, in K-symetrized space, with analytic formulae
{
	double jK1, jK2;
	std::complex<double> eK1, eK2;
	
//	double psi=0;
	// indices of the row and column vectors in the spin Liouville space
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

	// eK=sqrt(jK) is 1 or i (computed here but not used further below)
	if(jK1==1)		eK1=1;
	else			eK1=I;

	if(jK2==1)		eK2=1;
	else			eK2=I;

// same N: spin-Hamiltonian terms of all interactions
if(basisr.N==basisc.N)
{
	// loop over the rank l=0,2 of the interaction tensors and their components m
	for (int l=0; l<=2; l+=2)
		for (int m=-l; m<=l; m++)
		{
			// orientational factors: 3j symbol in L,M (temp3) and rotation d(l,m,ms;psi) with the L normalisation (temp2); skip if either vanishes
			int ms=basisr.M-basisc.M;
			int mss=basisr.K-basisc.K;
			std::complex<double> temp3=Wigner3jMy(basisr.L, l, basisc.L, -basisr.M, ms, basisc.M);
			if(real(temp3)==0) continue;

			std::complex<double> temp1=0;
//			std::complex<double> temp2=wigner_d(l, m, ms, psi)*sqrt((double)(2*basisc.L+1)*(2*basisr.L+1));
			std::complex<double> temp2=wigner_d_stored(l, m, ms)*sqrt((double)(2*basisc.L+1)*(2*basisr.L+1));
			if(abs(temp2)==0) continue;

			// 3j symbols in L,K: temp4 for the K difference, temp5 for the K sum (with jK2 and the (L+K) parity sign)
			std::complex<double> temp4=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K-basisc.K, basisc.K);
			std::complex<double> temp5=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K+basisc.K, -basisc.K)*jK2;
			if( (basisc.L+basisc.K)%2 ) temp5=-temp5;

			// sum over the interactions: Fdu in the molecular frame (real part for equal jK, imaginary part otherwise) times the spin matrix element;
			// the last interaction (Lorentz linewidth) has imaginary matrix elements and enters with a factor i
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

			// accumulate the contribution of this (l,m)
			temp=temp+temp1*temp2*temp3;
		}

	// sign and normalisation of the K-symmetrised functions (1/2 for K=0, 1/sqrt(2) otherwise, for each of the two vectors)
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
	// different N: only the J2 term couples vectors with identical L, M, K and jK (|N difference|<=2)
	if ((basisr.L==basisc.L)&&
		(basisr.M==basisc.M)&&
		(basisr.K==basisc.K)&&
		(jK1==jK2))
	{
		// only the isotropic (rank 0) part of the J2 interaction contributes
		int l=0, m=0;
		int ms=basisr.M-basisc.M;
		int mss=basisr.K-basisc.K;
		std::complex<double> temp3=Wigner3jMy(basisr.L, l, basisc.L, -basisr.M, ms, basisc.M);
		std::complex<double> temp1=0;
		std::complex<double> temp2=wigner_d_stored(l, m, ms)*sqrt((double)(2*basisc.L+1)*(2*basisr.L+1));
		std::complex<double> temp4=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K-basisc.K, basisc.K);
		std::complex<double> temp5=Wigner3jMy(basisr.L, l, basisc.L, -basisr.K, basisr.K+basisc.K, -basisc.K)*jK2;
		if( (basisc.L+basisc.K)%2 ) temp5=-temp5;

		// J-only interaction Interactions[4]: weight 2 between N=0 and N=2, 2*sqrt(2) for other |N difference|<=2, zero otherwise
		int i=4;   // (only J)
		if((basisr.N==0)&&(basisc.N==2) || (basisr.N==2)&&(basisc.N==0)) temp1=( temp4*2.0000000000000000000000000000000*J2 + temp5*2.0000000000000000000000000000000*J2 )*conj(Interactions[i].AuL.GetElement(l, m, j, k));
		else if (abs(basisr.N-basisc.N)<=2)								 temp1=( temp4*2.8284271247461900976033774484194*J2 + temp5*2.8284271247461900976033774484194*J2 )*conj(Interactions[i].AuL.GetElement(l, m, j, k));
		else															 temp1=0;

		// combine the factors, then apply the same K-symmetrisation sign and normalisation as above
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
// Normalisation of the K-symmetrised pair (K1,K2): a factor 1/sqrt(2) for each of K1, K2 that is zero, otherwise 1.
// Used by: MatrixElementIsotropicBrownianB_ND, MatrixElementIsotropicBrownianB, MatrixElementIsotropicBrownianBSM_ND,
// MatrixElementPotentialB_ND, MatrixElementPotentialB, MatrixElementPotentialBSM_ND, MatrixElementJumpsBSM_Pot, MatrixElementJumpsBSM_Iso (this file).
double Nk(double K1, double K2)
{
	// each zero K contributes a factor 1/sqrt(2)
	double n=1.0;
	if (K1==0) n *= 0.7071067811865475;
	if (K2==0) n *= 0.7071067811865475;
	return n;
}

// Coefficient sqrt((L-K-1)(L-K)(L+K+1)(L+K+2)) of the step K -> K+2 of the rotational-diffusion operator within rank L.
// Used by: MatrixElementIsotropicBrownianB_ND, MatrixElementIsotropicBrownianBSM_ND, MatrixElementIsotropicBrownianB, XLK_ND, XLK (this file).
double Np(double L, double K)
{
	return sqrt((L-K-1)*(L-K)*(L+K+1)*(L+K+2));
}
// Coefficient sqrt((L+K-1)(L+K)(L-K+1)(L-K+2)) of the step K -> K-2 of the rotational-diffusion operator within rank L.
// Used by: MatrixElementIsotropicBrownianB_ND, MatrixElementIsotropicBrownianBSM_ND, MatrixElementIsotropicBrownianB, XLK_ND, XLK (this file).
double Nm(double L, double K)
{
	return sqrt((L+K-1)*(L+K)*(L-K+1)*(L-K+2));
}
// Coefficient sqrt((L+K+1)(L-K)) of the step K -> K+1 (K raising operator) within rank L.
// Used by: MatrixElementIsotropicBrownianB_ND, MatrixElementIsotropicBrownianBSM_ND, XLK_ND, XLK (this file).
double Mp(double L, double K)
{
	return sqrt((L+K+1)*(L-K));
}
// Coefficient sqrt((L-K+1)(L+K)) of the step K -> K-1 (K lowering operator) within rank L.
// Used by: MatrixElementIsotropicBrownianB_ND, MatrixElementIsotropicBrownianBSM_ND, XLK_ND, XLK (this file).
double Mm(double L, double K)
{
	return sqrt((L-K+1)*(L+K));
}
// Disabled: the definition below is inside a comment block. MatrixElementIsotropicBrownianB_ND_nosym was the rotational-diffusion element
// without K-symmetrisation (R is the 3x3 diffusion tensor).
// Used by: nothing in this code base (only referenced from commented-out code in MatrixElementIsotropicBrownianB_ND and
// MatrixElementIsotropicBrownianBSM_ND).
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
// Rotational-diffusion part of the Liouvillian in the K-symmetrised basis, for free diffusion in an isotropic medium (no ordering potential).
// It does not change L, M, N or any spin index. R is the 3x3 diffusion tensor (R(0,0), R(1,1), R(2,2) are the xx, yy, zz components); the diagonal
// element is ((Rx+Ry)/2)*(L(L+1)-K^2)+Rz*K^2 and the couplings |dK|=1,2 come from the non-axial components of R.
// If any entry of Diffusion_Type is nonzero, the L(L+1), K^2 and Rz*K^2 terms are damped by (1+R_Tau*x)^(-Diffusion_Type) (non-Brownian
// diffusion, diagonal elements only); R_Tau and Diffusion_Type have three entries each.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) for biradical conformers.
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
	// diffusion does not change L, M, N or any spin index
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
		// work with |K|; the sign information is carried by jK1, jK2
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
		// non-Brownian: only the diagonal element, with each rate term damped by (1+R_Tau*x)^(-Diffusion_Type)
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
		// Brownian, equal jK: diagonal element and the couplings dK=+-1 (via R(0,2)) and dK=+-2 (via R(0,0)-R(1,1))
		if(jK1==jK2)
		{
			if(K1==K2)
			{
				element=(R(0,0)+R(1,1))*0.5*(L1*(L1+1)-K1*K1)+R(2,2)*K1*K1;
				// special case K1=1: the K=+1 and K=-1 functions are coupled by the dK=2 term, which adds to the diagonal element
				if(K1==1)    // special case K1=-1, K2=1!
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
			// different jK: couplings that mix the two combinations, via the tensor components R(0,1) and R(1,2)
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
// Same as MatrixElementIsotropicBrownianB_ND for a diagonal diffusion tensor with components Rx, Ry, Rz, without the non-Brownian option:
// diagonal elements and the dK=+-2 couplings between vectors with equal jK.
// Used by: nothing in this code base.
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


	// diffusion does not change L, M, N, jK or any spin index
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
		// diagonal element and the dK=+-2 couplings of the axial-tensor diffusion operator
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
// Coefficient X(L,K) (M=0 component) of the diffusion-in-potential superoperator for the ordering potential pot and a 3x3 rotational diffusion
// tensor R that need not be diagonal in the potential frame. The first part is linear in the potential coefficients pot.Get(L,K) (terms with
// dK=0,+-1,+-2 weighted by components of R); the second part is a double sum over pairs of potential components (L1,K1), (L2,K2) coupled to rank L
// by 3j symbols. Zero for odd L or |K|>L. Real and imaginary parts serve the jK-equal and jK-different elements of the potential term.
// Used by: Potential::CalculateXLK (Potential.cpp), which stores the values in XLK_Matrix; MatrixElementPotentialB_ND and MatrixElementPotentialBSM_ND
// read the stored values through Potential::GetXLK (calls to XLK_ND there are commented out).
std::complex<double> XLK_ND(int L, int K, Potential &pot, Matrix<std::complex<double>, 3, 3> &R)
// changed R matrix from complex<double> to double!
// for non-diagonal diffusion tensor!
{
	std::complex<double> result;
	std::complex<double> sum(0, 0);
	std::complex<double> I(0,1);
	// only even L and |K|<=L contribute
	if(L%2) return 0;
	if (abs(K)>L) return 0;
//	if(K%2) return 0;
//	if(K<0) K=-K;
	// part linear in the potential: terms with dK=+2,-2,+1,-1,0 weighted by the diffusion-tensor components
	result=-0.5*pot.Get(L, K+2)*Nm((double)L, (double)(K+2))*0.25*(R(0,0)-R(1,1)+2.0*I*R(0,1));
	result-=0.5*pot.Get(L, K-2)*Np((double)L, (double)(K-2))*0.25*(R(0,0)-R(1,1)-2.0*I*R(0,1));
	result-=0.5*pot.Get(L, K+1)*Mm((double)L, (double)(K+1))*(2.0*(double)K+1)*0.5*(R(0,2)+I*R(1,2));
	result-=0.5*pot.Get(L, K-1)*Mp((double)L, (double)(K-1))*(2.0*(double)K-1)*0.5*(R(0,2)-I*R(1,2));
	result-=0.5*(pot.Get(L, K)*(L*(L+1)-K*K)*(R(0,0)+R(1,1))*0.5+pot.Get(L, K)*K*K*R(2,2));

//	printf("L %d K %d el=%lf\n", L, K, result);


	// part quadratic in the potential: sum over pairs of components (L1,K1), (L2,K2) coupled to rank L
	for(int L1=0; L1<=pot.Lmax; L1+=2)
		for(int L2=0; L2<=pot.Lmax; L2+=2)
			for(int K1=-L1; K1<=L1; K1+=2)
				for(int K2=-L2; K2<=L2; K2+=2)
				{
					// product of the two potential coefficients and the 3j symbol over L1,L,L2; skip if it vanishes
					std::complex<double> el;
					el=pot.Get(L1, K1)*pot.Get(L2, K2);
					el*=Wigner3jMy(L1, L, L2, 0, 0, 0);
					if(abs(el)==0) continue;
					// 3j symbols with K shifted by 0 or +-1, times the corresponding ladder coefficients
					double temp1, temp2, temp3, temp4, temp5, temp6;
					temp1=Mp((double)L1, (double)K1)*Mp((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2+1);
					temp2=Mm((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1-1, -K, K2-1);
					temp3=Mp((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2-1);
					temp4=												  K1*K2*Wigner3jMy(L1, L, L2, K1,   -K, K2);
					temp5=Mp((double)L1, (double)K1)*						 K2*Wigner3jMy(L1, L, L2, K1+1, -K, K2);
					temp6=Mm((double)L1, (double)K1)*						 K2*Wigner3jMy(L1, L, L2, K1-1, -K, K2);
					// weight each term by its diffusion-tensor component
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
	// combine the linear and quadratic parts
	if(L%2) return result+sum*0.25*(2.0*L+1);
	else	return result-sum*0.25*(2.0*L+1);


}

// XLK decomposition components in the case of diagonal rotation diffusion tensor.
// Real version of XLK_ND for a diagonal rotational diffusion tensor with components Rx, Ry, Rz in the potential frame; zero for odd L or odd K.
// Used by: MatrixElementPotentialB (this file), which is itself never called.
double XLK(int L, int K, Potential &pot, double Rx, double Ry, double Rz)
{
	double result, sum=0;
	// only even L and even K contribute
	if(L%2) return 0;
	if(K%2) return 0;
//	if(K<0) K=-K;
	// part linear in the potential
	result=-0.5*(pot.Get(L, K+2)*Nm((double)L, (double)(K+2))+pot.Get(L, K-2)*Np((double) L, double(K-2)))*0.25*(Rx-Ry);
	result-=0.5*(pot.Get(L, K)*(L*(L+1)-K*K)*(Rx+Ry)*0.5+pot.Get(L, K)*K*K*Rz);

//	printf("L %d K %d el=%lf\n", L, K, result);


	// part quadratic in the potential: sum over pairs of components (L1,K1), (L2,K2) coupled to rank L
	for(int L1=0; L1<=pot.Lmax; L1+=2)
		for(int L2=0; L2<=pot.Lmax; L2+=2)
			for(int K1=-L1; K1<=L1; K1+=2)
				for(int K2=-L2; K2<=L2; K2+=2)
				{
					// product of the two potential coefficients and the 3j symbol over L1,L,L2; skip if it vanishes
					double el=pot.Get(L1, K1)*pot.Get(L2, K2);
					el*=Wigner3jMy(L1, L, L2, 0, 0, 0);
					if(el==0) continue;
					// 3j symbols with K shifted by 0 or +-1, times the corresponding ladder coefficients
					double temp1, temp2, temp3, temp4;
					temp1=Mp((double)L1, (double)K1)*Mp((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2+1);
					temp2=Mm((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1-1, -K, K2-1);
					temp3=Mp((double)L1, (double)K1)*Mm((double)L2, (double)K2)*Wigner3jMy(L1, L, L2, K1+1, -K, K2-1);
					temp4=												  K1*K2*Wigner3jMy(L1, L, L2, K1, -K, K2);

					// weight each term by its diffusion-tensor component
					el*=0.25*(Rx-Ry)*(temp1+temp2)+0.5*(Rx+Ry)*temp3+Rz*temp4;
//					cout << Wigner3jMy(2, 2, 2, 1, -2, -1) << endl;
//					printf("L %d K %d L1 %d L2 %d K1 %d K2 %d el=%lf\n", L, K, L1, L2, K1, K2, el);
					sum+=el;
				}
	// combine the linear and quadratic parts
	return result-sum*0.25*(2.0*L+1);
}	

// Prefactor sqrt((2*L1+1)*(2*L2+1)) of the orientational overlap between functions of ranks L1 and L2.
// Used by: MatrixElementPotentialB_ND, MatrixElementPotentialB, MatrixElementPotentialBSM_ND, MatrixElementJumpsBSM_Pot (this file).
double Nl(double L1, double L2)
{
	return sqrt((2*L1+1)*(2*L2+1));
}

// Matrix element of Brownian diffusion in anisotropic medium, in the case of non-coaxial orientation and diffusion:
// Rotational diffusion in the ordering (mean-field) potential, K-symmetrised basis, for a diffusion tensor that need not be coaxial with the
// potential (_ND). The element couples vectors with different L and K through the stored coefficients X(L,K) (pot.GetXLK) summed over the even ranks
// L=0..2*pot.Lmax with Wigner 3j symbols; M, N and all spin indices are unchanged. The R argument is not used in this function: the diffusion
// tensor enters through the XLK table filled by Potential::CalculateXLK.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) when MFPotential.Lmax is nonzero, for biradical conformers.
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


	// the potential term does not change M, N or any spin index; L and K may change
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

		// sum over the even ranks L of the potential-diffusion coefficients
		for(int L=0; L<=2*pot.Lmax; L+=2)
		{
			// 3j symbol in L1,L,L2 and M; skip if forbidden
			temp1=Wigner3jMy(L1, L, L2, M1, 0, -M2);
			if (temp1==0) continue;
			// XLK for K1-K2 (temp2) and K1+K2 (temp3) times 3j symbols in K: real part for equal jK, imaginary part otherwise
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
		// normalisation and K-symmetrisation sign (and the factor jK2 for unequal jK)
		sum*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum=-sum;
		if(jK1!=jK2) sum=jK2*sum;
		element=sum;
	}

	return element;


}

// The case of coaxial diffusion and orientation (according to Freed at al.)
// Same as MatrixElementPotentialB_ND for coaxial diagonal diffusion tensor (Rx,Ry,Rz) and potential, with the coefficients XLK computed directly
// by XLK() rather than read from the stored table; only vectors with equal jK are coupled.
// Used by: nothing in this code base.
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


	// the potential term does not change M, N, jK or any spin index; L and K may change
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

		// sum over the even ranks L of the potential-diffusion coefficients
		for(int L=0; L<=2*pot.Lmax; L+=2)
		{
			double temp1, temp2, temp3;
			// 3j symbol in L1,L,L2 and M; skip if forbidden
			temp1=Wigner3jMy(L1, L, L2, M1, 0, -M1);
			if (temp1==0) continue;
			// XLK for K1-K2 (temp2) and K1+K2 (temp3) times 3j symbols in K
			temp2=XLK(L, K1-K2, pot, Rx, Ry, Rz)*Wigner3jMy(L1, L, L2, K1, K2-K1, -K2);
			temp3=XLK(L, K1+K2, pot, Rx, Ry, Rz)*Wigner3jMy(L1, L, L2, K1, -K2-K1, K2)*jK2;
			if((L2+K2)%2) temp3=-temp3;

			sum+=temp1*(temp2+temp3);
		}
		// normalisation and K-symmetrisation sign
		sum*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum=-sum;
		element=sum;
	}

	return element;
}



// M-Symmetrized, for High-field approximation
// M-symmetrised monoradical version of MatrixElementPotentialB_ND: the same stored XLK coefficients, 3j symbols and real/imaginary choice by jK,
// but vectors with different M-symmetrisation label jM are not coupled. R is not used (see MatrixElementPotentialB_ND).
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) when MFPotential.Lmax is nonzero, for monoradical conformers.
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

	// M, N, jM and all spin indices must match; L and K may change
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

		// sum over the even ranks L of the potential-diffusion coefficients
		for(int L=0; L<=2*pot.Lmax; L+=2)
		{
			// 3j symbol in L1,L,L2 and M; skip if forbidden
			temp1=Wigner3jMy(L1, L, L2, M1, 0, -M2);
			if (temp1==0) continue;
			// XLK for K1-K2 (temp2) and K1+K2 (temp3) times 3j symbols in K: real part for equal jK, imaginary part otherwise
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
		// normalisation and K-symmetrisation sign (and the factor jK2 for unequal jK)
		sum*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum=-sum;
		if(jK1!=jK2) sum=jK2*sum;
		element=sum;
	}

//	time_potential+=(clock()-time_start);
	return element;


}




// M-symmetrised monoradical version of MatrixElementIsotropicBrownianB_ND: the same K couplings of rotational diffusion, but vectors
// with different M-symmetrisation label jM are not coupled. Arguments as for MatrixElementIsotropicBrownianB_ND.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) for monoradical conformers.
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


	// diffusion does not change L, M, N, jM or any spin index
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
		// work with |K|; the sign information is carried by jK1, jK2
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
		// non-Brownian: only the diagonal element, with each rate term damped by (1+R_Tau*x)^(-Diffusion_Type)
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
		// Brownian, equal jK: diagonal element and the couplings dK=+-1 (via R(0,2)) and dK=+-2 (via R(0,0)-R(1,1))
		if(jK1==jK2)
		{
			if(K1==K2)
			{
				element=(R(0,0)+R(1,1))*0.5*(L1*(L1+1)-K1*K1)+R(2,2)*K1*K1;
				// special case K1=1: the K=+1 and K=-1 functions are coupled by the dK=2 term, which adds to the diagonal element
				if(K1==1)    // special case K1=-1, K2=1!
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
			// different jK: couplings that mix the two combinations, via the tensor components R(0,1) and R(1,2)
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


// Disabled: the definition below is inside a comment block, although the declaration remains in LiouvilleMatrixElementB.h.
// MatrixElementDiscreteJumpsB was a model of discrete jumps between NSites sites: JumpFrequency on the diagonal for K not divisible by NSites.
// Used by: nothing in this code base (calling it would fail at link time).
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



// Discrete orientational jumps (NJumpAxes jump axes with polar angles JumpAxisFi, JumpAxisTeta, rotation angle JumpAngle and rate JumpFrequencyR)
// for a monoradical in the M-symmetrised basis, in the presence of an ordering potential, which enters through the order parameters pot.GetOP and
// the per-axis jump order parameters pot.GetOP_JUMP_P / pot.GetOP_JUMP_M. Real parts are used for equal jK, imaginary parts for unequal jK.
// The remark 'for non-diagonal diffusion tensors!' below the signature is a leftover: this function takes no diffusion tensor.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) when the conformer has jump axes (NJumpAxes) and MFPotential.Lmax is nonzero (monoradicals only).
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

	// jumps do not change M, N, jM or any spin index
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

		// sum over the ranks L allowed by the triangle rule between L1 and L2
		for(int L=abs(L1-L2); L<=L1+L2; L++)
		{
			// 3j symbol in M with weight (2L+1); skip if forbidden
			temp1=Wigner3jMy(L1, L, L2, -M1, 0, M1);
			if (temp1==0) continue;
			temp1 *= (2*L+1);
			// equal jK: real parts (the unequal-jK branch further below uses imaginary parts)
			if(jK1==jK2)
			{

			// loop over the jump axes, each weighted by its jump rate JumpFrequencyR[i]
			for(int i=0; i<MagneticParameter->NJumpAxes; i++)
			{
				sum=0; 
				// jump axis along z (JumpAxisFi=JumpAxisTeta=0): the rotations into the axis frame are trivial, so the terms are written out directly
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

				// general axis: sum over components P,Q, rotating to the jump-axis frame and back with ComplexWigner
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
				// remaining terms, which involve the jump order parameters pot.GetOP_JUMP_P / pot.GetOP_JUMP_M
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
				// weight by the jump rate of this axis and accumulate with the M-dependent 3j factor
				sum*=MagneticParameter->JumpFrequencyR[i];
				sum2+=sum*temp1;

			}

			}
			else // if(jK1 != jK2)
			{
			// unequal jK: same structure as above, taking imaginary parts
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
		// normalisation and K-symmetrisation sign (and the factor jK2 for unequal jK)
		sum2*=Nl((double)L1, (double)L2)*Nk((double)K1, (double)K2);
		if ((M1+K1)%2) sum2=-sum2;
		if(jK1!=jK2) sum2=jK2*sum2;
		element=sum2;
	}

//	time_potential+=(clock()-time_start);
	return element;


}


// Discrete orientational jumps as in MatrixElementJumpsBSM_Pot, but without an ordering potential (isotropic distribution): diagonal in L and M.
// For each jump axis the element is rate*(unit term - rotation term), the rotation term being a sum over Q of cos(Q*JumpAngle) times two
// ComplexWigner rotations through the jump-axis frame.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp) when the conformer has jump axes (NJumpAxes) and MFPotential.Lmax is zero (monoradicals only).
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

	// jumps do not change L, M, N, jM or any spin index
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
		

		// equal jK: real parts, with the unit (no-jump) term for K1==K2 and K1==-K2
		if(jK1==jK2)
		{

			// loop over the jump axes
			for(int i=0; i<MagneticParameter->NJumpAxes; i++)
			{
				temp1=0;
				temp2=0;
				temp3=0;

				// unit term: K1==K2 and, for the mirror component, K1==-K2 with the (L2+K2) parity sign and jK2
				if(K1==K2)  temp1+=1;
				if(K1==-K2) 
				{
					if((L2+K2)%2)	temp1+=-jK2;
					else			temp1+= jK2;
				}
				// rotation terms for K2 (temp2) and -K2 (temp3): sum over Q of cos(Q*JumpAngle) times rotations to and from the jump-axis frame
				for (int Q=-L2; Q<=L2; Q++)
				{
					temp2+=cos((double)Q*MagneticParameter->JumpAngle[i])*
											real(ComplexWigner(L2, K1, Q, MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)*ComplexWigner(L2, Q, K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i]));
					temp3+=cos((double)Q*MagneticParameter->JumpAngle[i])*
											real(ComplexWigner(L2, K1, Q, MagneticParameter->JumpAxisFi[i], MagneticParameter->JumpAxisTeta[i], 0)*ComplexWigner(L2, Q,-K2, 0, -MagneticParameter->JumpAxisTeta[i], -MagneticParameter->JumpAxisFi[i]));
				}
				// parity sign and jK2 for the mirror term, then add this axis weighted by its jump rate
				if((L2+K2)%2)	temp3*=-jK2;
				else			temp3*= jK2;
				sum+=(temp1-temp2-temp3)*MagneticParameter->JumpFrequencyR[i];

				}

		}
		else
		{
			// unequal jK: imaginary parts only, no unit term
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


			
		// K-symmetrisation normalisation, factor jK2 for unequal jK, and the overall factor 2
		sum*=Nk((double)K1, (double)K2);
		if(jK1!=jK2) sum=jK2*sum;
		element=2.0*sum;
	}

//	time_potential+=(clock()-time_start);
	return element;


}


// Conformer-exchange part of the Liouvillian: between vectors that are identical in every index except the conformer index Q, the element is
// the table entry ConversionRate[Q of row][Q of column] (row-major with stride MAX_PARTICLES); for equal Q it is the diagonal entry of the table.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp): for vectors of different conformers, and added to every element of equal conformers.
std::complex<double> MatrixElementConformation(BasisVectorB basisr, BasisVectorB basisc, double* ConversionRate)
{
	std::complex<double> element(0,0);
	// only vectors identical in L,K,M,N and all spin indices are coupled
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
		// exchange rate for the pair of conformer indices (Q of the row, Q of the column)
		element=ConversionRate[MAX_PARTICLES*basisr.Q + basisc.Q];
	}

//	time_potential+=(clock()-time_start);
	return element;

}

// J-diffusion part of the Liouvillian: a diagonal damping JDiffusionRate*N*N of the extra index N (the index coupled by the J2 term of
// MatrixElementBSK1), identical for all vectors; zero for N=0.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp), added to every element.
std::complex<double> MatrixElementJDiffusion(BasisVectorB basisr, BasisVectorB basisc, double JDiffusionRate)
{
	std::complex<double> element(0,0);
	// diagonal only: all indices of the two vectors must match
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
		// damping grows as N squared
		element=JDiffusionRate*basisc.N*basisc.N;
	}

//	time_potential+=(clock()-time_start);
	return element;

}

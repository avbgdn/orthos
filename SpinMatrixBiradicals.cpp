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
//
// The calcluation of spin matrices corresponding to different Spin Hamiltonian terms for biradicals
//

// SpinMatrixBiradicals.cpp: the spin-operator matrices of the biradical Spin Hamiltonian in the spherical-tensor form. Each
// interaction is the contraction of a spherical tensor F(l,m) (class SphericalTensor) with the 36x36 spin matrices A(l,m)
// built here. The spin space is two electron spins S=1/2 and two nuclear spins I=1, in the Kronecker-product order electron
// 1, nucleus 1, electron 2, nucleus 2 (2*3*2*3=36 states). The file defines the class InteractionSphericalTensorB (stores the
// A(l,m) of one interaction and returns elements of the matching Liouville-space superoperator), the six builders Ag1L,
// AA1L, Ag2L, AA2L, ADL and AWLB, StartingVectorSpinB, and the helpers delta and PSMB.
// Used by: RWMagnetic.cpp (MagneticParameters::RotateTensors passes the builders to InteractionSphericalTensorB::Set),
// LiouvilleMatrixElementB.cpp (MatrixElementBSK1 reads the matrices through GetElement), LiouvillianB.cpp
// (CalculateLiouvillianB calls StartingVectorSpinB).

#include "BiradicalODF.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "SphericalTensor.h"
#include "SpinMatrixB.h"
// kron(A,B,C,D) is the Kronecker product of four matrices: electron 1, nucleus 1, electron 2, nucleus 2
#define kron(A, B, C, D) kroneckerProduct(A, kroneckerProduct(B, kroneckerProduct(C, D)))
//#define kron(A, B, C, D) kroneckerProduct(kroneckerProduct(kroneckerProduct(A, B), C), D)
// dimension of the spin Hilbert space: 2*3*2*3 states
#define MATSIZE 36

//extern double B0;

// global spin matrices; defined in SphericalTensor.cpp and filled by InitializeSpinMatrices
extern Matrix<std::complex<double>, 2, 2> Sx, Sy, Sz, Si, Sp, Sm;
extern Matrix<std::complex<double>, 3, 3> Ix, Iy, Iz, Ii, Ip, Im;


// Builds and stores the spin matrices of one interaction: calls ALFunction(0,0,B0) for rank 0 and ALFunction(2,m,B0) for
// m=-2..2. ALFunction is one of the builders Ag1L, AA1L, Ag2L, AA2L, ADL, AWLB, passed as a function pointer; B0 is the
// magnetic field (RWMagnetic.cpp computes it from CenterField and the trace of g1R), and only Ag1L and Ag2L use it.
// Used by: MagneticParameters::RotateTensors (RWMagnetic.cpp), once for each of the six interactions Interactions[0..5].
void InteractionSphericalTensorB::Set(Matrix<std::complex<double>, MATSIZE, MATSIZE> ALFunction (int, int, double), double B0)
{
// l=0
		AL0[0]=ALFunction(0, 0, B0);
// l=2
		int l=2;
		for(int m=-l; m<=l; m++)
			{
				AL2[m+2]=ALFunction(l, m, B0);
			}
	}
//	Matrix<std::complex<double>, MATSIZE, MATSIZE> * InteractionSphericalTensorB::GetMatrix(int l, int m)
//	{
//		if((l<0)||(l>2)) return (Matrix<std::complex<double>, MATSIZE, MATSIZE>*)NULL;
//		if(abs(m)>l) return (Matrix<std::complex<double>, MATSIZE, MATSIZE> *)NULL;
//		if (l==0) return AL0;
//		else if (l==2) return AL2[m+2];
//		else return (Matrix<std::complex<double>, MATSIZE, MATSIZE> *)NULL;
//
//	}
	// Prints the 36x36 matrix A(l,m) for l=0 or l=2 (and, for l=0, also its first Liouville-space element) to the screen.
	// Used by: nothing in this code base (a call in MagneticParameters::RotateTensors in RWMagnetic.cpp is commented out).
	void InteractionSphericalTensorB::PrintSmallMatrix(int l, int m)
	{
		// rank 0: print the matrix and the first element of its superoperator
		if(l==0) 
		{
//				for (int i=0; i<MATSIZE; i++)
//				{
//					for (int j=0; j<MATSIZE; j++)
//						printf("% 8.5lf   ", real(this->GetElement(l, m, i, j)));
//						printf("(% 8.2lf   % 8.2lf )  ", real(AL0[0](i, j)), imag(AL0[0](i, j)));
//					printf("\n");
//				}
			cout << l << m << endl << AL0[0] << endl << endl;
			cout <<	PSMB(AL0[0], 0, 0);

		}
		// rank 2: print the matrix of the requested m
		if(l==2) cout << l << m << endl << AL2[m+2] << endl << endl;
	}
	// Returns element (i,j) of the Liouville-space superoperator (commutator with the stored matrix A(l,m), see PSMB). i and j
	// run over the 36*36=1296 spin Liouville states; MatrixElementBSK1 obtains them from the spin parts of the two basis vectors
	// through SpinIndexB. Returns 0 for l outside 0..2, |m|>l or an index out of range.
	// Used by: MatrixElementBSK1 (LiouvilleMatrixElementB.cpp).
	std::complex<double> InteractionSphericalTensorB::GetElement(int l, int m, int i, int j)
	{
		// invalid rank or projection: no such matrix
		if((l<0)||(l>2)) return 0;
		if(abs(m)>l) return 0;
		if((i>=MATSIZE*MATSIZE)||(j>=MATSIZE*MATSIZE)) return 0;
	
		// superoperator element of the stored matrix of the requested rank
		if (l==0) return PSMB(AL0[0], i, j);
		else if(l==2) return PSMB(AL2[m+2], i, j);
		else return 0;
	}
	// Prints all stored matrices A(l,m) (l=0,2) as 1296x1296 Liouville-space matrices, element by element, to the screen.
	// Used by: nothing in this code base.
	void InteractionSphericalTensorB::PrintOut()
	{
		// loop over the ranks and projections and print every superoperator element
		for (int l=0; l<=2; l+=2)
			for(int m=-l; m<=l; m++)
			{
				printf("\n\nA(%d, %d):\n", l, m);
				for (int i=0; i<MATSIZE*MATSIZE; i++)
				{
					for (int j=0; j<MATSIZE*MATSIZE; j++)
//						printf("% 8.5lf   ", real(this->GetElement(l, m, i, j)));
						printf("(% 8.2lf   % 8.2lf )", real(this->GetElement(l, m, i, j)), imag(this->GetElement(l, m, i, j)));
					printf("\n");
				}

			}
		

	}

	// Allocates (with malloc) storage for one rank-0 matrix and five rank-2 matrices (m=-2..2), each 36x36.
	// Used by: implicit construction of SpinInteractionB objects, for example MagneticParameters::Interactions.
	InteractionSphericalTensorB::InteractionSphericalTensorB()
	{
		AL0=(Matrix<std::complex<double>, 36, 36> *)malloc(sizeof(Matrix<std::complex<double>, 36, 36>)*1);
		AL2=(Matrix<std::complex<double>, 36, 36> *)malloc(sizeof(Matrix<std::complex<double>, 36, 36>)*5);
	
	}
	// Releases the storage allocated by the constructor.
	// Used by: implicit destruction of SpinInteractionB objects.
	InteractionSphericalTensorB::~InteractionSphericalTensorB()
	{
		free(AL0);
		free(AL2);
	}




// Spin matrix A(l,m) of the electron Zeeman interaction of electron 1 with the field B0 along z, as a 36x36 matrix. The
// contraction of the field vector (0,0,B0) with the spin vector of electron 1 leaves l=0: -B0*Sz/sqrt(3); l=2, m=0:
// sqrt(2/3)*B0*Sz; l=2, m=+1: -0.5*B0*S+; l=2, m=-1: 0.5*B0*S-. Every other (l,m) gives the zero matrix.
// Used by: passed as a function pointer to InteractionSphericalTensorB::Set for Interactions[0] in
// MagneticParameters::RotateTensors (RWMagnetic.cpp).
Matrix<std::complex<double>, 36, 36> Ag1L( int l, int m , double B0)
{
	// make sure the global spin matrices are filled
	InitializeSpinMatrices();

//	Matrix<std::complex<double>, ColMajor, int> Result(1296, 1296);
	// the zero matrix is returned for (l,m) that do not occur
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);
	// each (l,m) component is a combination of Kronecker products of the spin matrices
	if((l==0)&&(m==0))
		result=-sqrt(0.33333333333333333333)*B0*kron(Sz, Ii, Si, Ii);		// sqrt(0.33333333333333333333)
	if((l==2)&&(m==-1))
        result=0.5*B0*kron(Sm, Ii, Si, Ii);
	if((l==2)&&(m==0))
        result=sqrt(0.66666666666666666667)*B0*kron(Sz, Ii, Si, Ii);		// sqrt(0.66666666666666666667)
	if((l==2)&&(m==1))
        result=-0.5*B0*kron(Sp, Ii, Si, Ii);
    
//	cout << B0 << endl;
//	cout << Sz << endl;
//	Result=kroneckerProduct(result, Ident)-kroneckerProduct(Ident, result.adjoint() );

   return result;
}


// Spin matrix A(l,m) of the hyperfine coupling of electron 1 with nucleus 1 (spin 1): the spherical components of the
// product of the spin vectors S and I. l=0: -(Sz*Iz+(S+*I-+S-*I+)/2)/sqrt(3); l=2: m=+-2: 0.5*S(+-)*I(+-); m=+-1:
// -+0.5*(S(+-)*Iz+Sz*I(+-)); m=0: sqrt(2/3)*(Sz*Iz-(S+*I-+S-*I+)/4). The argument B0 is not used.
// Used by: passed as a function pointer to InteractionSphericalTensorB::Set for Interactions[2] in
// MagneticParameters::RotateTensors (RWMagnetic.cpp).
Matrix<std::complex<double>, 36, 36> AA1L( int l, int m , double B0)
{
	// make sure the global spin matrices are filled
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	// the zero matrix is returned for (l,m) that do not occur
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);

	// each (l,m) component is a combination of Kronecker products of the spin matrices
	if((l==0)&&(m==0))
		result=-sqrt(0.33333333333333333333)*(kron(Sz, Iz, Si, Ii) + 0.5*kron(Sp, Im, Si, Ii) + 0.5*kron(Sm, Ip, Si, Ii));
	if((l==2)&&(m==-2))
        result=0.5*kron(Sm, Im, Si, Ii);
	if((l==2)&&(m==-1))
        result=0.5*(kron(Sm, Iz, Si, Ii) + kron(Sz, Im, Si, Ii));
	if((l==2)&&(m==0))
        result=sqrt(0.66666666666666666667)*(kron(Sz, Iz, Si, Ii) - 0.25*kron(Sp, Im, Si, Ii) - 0.25*kron(Sm, Ip, Si, Ii));
	if((l==2)&&(m==1))
        result=-0.5*(kron(Sp, Iz, Si, Ii) + kron(Sz, Ip, Si, Ii));
	if((l==2)&&(m==2))
        result=0.5*kron(Sp, Ip, Si, Ii);
//   Result=kroneckerProduct(result, Ident)-kroneckerProduct(Ident, result.adjoint() );

   return result;
}

// Spin matrix A(l,m) of the electron Zeeman interaction of electron 2 with the field B0 along z: the same components as
// Ag1L (l=0: -B0*Sz/sqrt(3); l=2: sqrt(2/3)*B0*Sz for m=0, -0.5*B0*S+ for m=+1, 0.5*B0*S- for m=-1) but for electron 2.
// Used by: passed as a function pointer to InteractionSphericalTensorB::Set for Interactions[1] in
// MagneticParameters::RotateTensors (RWMagnetic.cpp).
Matrix<std::complex<double>, 36, 36> Ag2L( int l, int m , double B0)
{
	// make sure the global spin matrices are filled
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	// the zero matrix is returned for (l,m) that do not occur
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);
	// each (l,m) component is a combination of Kronecker products of the spin matrices
	if((l==0)&&(m==0))
		result=-sqrt(0.33333333333333333333)*B0*kron(Si, Ii, Sz, Ii);
	if((l==2)&&(m==-1))
        result=0.5*B0*kron(Si, Ii, Sm, Ii);
	if((l==2)&&(m==0))
        result=sqrt(0.66666666666666666667)*B0*kron(Si, Ii, Sz, Ii);
	if((l==2)&&(m==1))
        result=-0.5*B0*kron(Si, Ii, Sp, Ii);
//   Result=kroneckerProduct(result, Ident)-kroneckerProduct(Ident, result.adjoint() );

   return result;
}

// Spin matrix A(l,m) of the hyperfine coupling of electron 2 with nucleus 2: the same components as AA1L but built from
// the operators of electron 2 and nucleus 2. The argument B0 is not used.
// Used by: passed as a function pointer to InteractionSphericalTensorB::Set for Interactions[3] in
// MagneticParameters::RotateTensors (RWMagnetic.cpp).
Matrix<std::complex<double>, 36, 36> AA2L( int l, int m , double B0)
{
	// make sure the global spin matrices are filled
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	// the zero matrix is returned for (l,m) that do not occur
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);

	// each (l,m) component is a combination of Kronecker products of the spin matrices
	if((l==0)&&(m==0))
		result=-sqrt(0.33333333333333333333)*(kron(Si, Ii, Sz, Iz) + 0.5*kron(Si, Ii, Sp, Im) + 0.5*kron(Si, Ii, Sm, Ip));
	if((l==2)&&(m==-2))
        result=0.5*kron(Si, Ii, Sm, Im);
	if((l==2)&&(m==-1))
        result=0.5*(kron(Si, Ii, Sm, Iz) + kron(Si, Ii, Sz, Im));
	if((l==2)&&(m==0))
        result=sqrt(0.66666666666666666667)*(kron(Si, Ii, Sz, Iz) - 0.25*kron(Si, Ii, Sp, Im) - 0.25*kron(Si, Ii, Sm, Ip));
	if((l==2)&&(m==1))
        result=-0.5*(kron(Si, Ii, Sp, Iz) + kron(Si, Ii, Sz, Ip));
	if((l==2)&&(m==2))
        result=0.5*kron(Si, Ii, Sp, Ip);
//   Result=kroneckerProduct(result, Ident)-kroneckerProduct(Ident, result.adjoint() );

   return result;
}

// Spin matrix A(l,m) of the coupling between electron 1 and electron 2 (dipolar and exchange tensor DJ): the same
// components as AA1L with the spin vector of electron 1 and that of electron 2 in place of S and I. The argument B0 is
// not used.
// Used by: passed as a function pointer to InteractionSphericalTensorB::Set for Interactions[4] in
// MagneticParameters::RotateTensors (RWMagnetic.cpp).
Matrix<std::complex<double>, 36, 36> ADL( int l, int m , double B0)
{
	// make sure the global spin matrices are filled
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	// the zero matrix is returned for (l,m) that do not occur
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);
	// each (l,m) component is a combination of Kronecker products of the spin matrices
	if((l==0)&&(m==0))
		result=-sqrt(0.33333333333333333333)*(kron(Sz, Ii, Sz, Ii) + 0.5*kron(Sp, Ii, Sm, Ii) + 0.5*kron(Sm, Ii, Sp, Ii));
	if((l==2)&&(m==-2))
        result=0.5*kron(Sm, Ii, Sm, Ii);
	if((l==2)&&(m==-1))
        result=0.5*(kron(Sm, Ii, Sz, Ii) + kron(Sz, Ii, Sm, Ii));
	if((l==2)&&(m==0))
        result=sqrt(0.66666666666666666667)*(kron(Sz, Ii, Sz, Ii) - 0.25*kron(Sp, Ii, Sm, Ii) - 0.25*kron(Sm, Ii, Sp, Ii));
	if((l==2)&&(m==1))
        result=-0.5*(kron(Sp, Ii, Sz, Ii) + kron(Sz, Ii, Sp, Ii));
	if((l==2)&&(m==2))
        result=0.5*kron(Sp, Ii, Sp, Ii);
//   Result=kroneckerProduct(result, Ident)-kroneckerProduct(Ident, result.adjoint() );

   return result;
}




// Spin matrix A(l,m) for the Lorentzian line width tensor: the components of the field-like contraction used in Ag1L but
// built from the total electron spin of both electrons (Sz, S+, S- of electron 1 plus the same of electron 2) and without
// the factor B0 (l=0: -(S1z+S2z)/sqrt(3); l=2: sqrt(2/3)*(S1z+S2z) for m=0, -0.5*(S1+ + S2+) for m=+1, 0.5*(S1- + S2-) for
// m=-1). MatrixElementBSK1 multiplies this last interaction by i, i.e. it acts as damping.
// Used by: passed as a function pointer to InteractionSphericalTensorB::Set for Interactions[5] in
// MagneticParameters::RotateTensors (RWMagnetic.cpp).
Matrix<std::complex<double>, 36, 36> AWLB( int l, int m , double B0)
{
	// make sure the global spin matrices are filled
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	// the zero matrix is returned for (l,m) that do not occur
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);

	// each (l,m) component is a combination of Kronecker products of the spin matrices
	if((l==0)&&(m==0))
		result=-sqrt(0.33333333333333333333)*(kron(Sz, Ii, Si, Ii)+kron(Si, Ii, Sz, Ii));
	if((l==2)&&(m==-1))
        result=0.5*(kron(Sm, Ii, Si, Ii)+kron(Si, Ii, Sm, Ii));
	if((l==2)&&(m==0))
        result=sqrt(0.66666666666666666667)*(kron(Sz, Ii, Si, Ii)+kron(Si, Ii, Sz, Ii));
	if((l==2)&&(m==1))
        result=-0.5*(kron(Sp, Ii, Si, Ii)+kron(Si, Ii, Sp, Ii));
//   Result=kroneckerProduct(result, Ident)-kroneckerProduct(Ident, result.adjoint() );

   return result;
}


// Spin part of the starting vector of the spectrum calculation, a vector of 1296=36*36 elements in the spin Liouville space:
// for each combined spin index (i,j,k,l) of electron 1, nucleus 1, electron 2, nucleus 2 one block of 36 elements
// kron(Si.col(i),Ii.col(j),Sx.col(k),Ii.col(l))+kron(Sx.col(i),Ii.col(j),Si.col(k),Ii.col(l)), i.e. the elements of the
// operator Sx of electron 1 plus Sx of electron 2.
// Used by: CalculateLiouvillianB (LiouvillianB.cpp), which multiplies it by the orientational equilibrium component of
// each basis vector (GetSOP).
Matrix<std::complex<double>, 1296, 1> StartingVectorSpinB()
{
	// make sure the global spin matrices are filled
	InitializeSpinMatrices();
	Matrix<std::complex<double>, 1296, 1> sv;
//	cout << kroneckerProduct(Sx.col(0), Ii.col(0)).cols() << "x" << kroneckerProduct(Sx.col(0), Ii.col(0)).rows();

    // loop over the four spin indices of a state; each combination fills one block of 36 elements
    for(int i=0; i<2; i++)
		for(int j=0; j<3; j++)
		    for(int k=0; k<2; k++)
				for(int l=0; l<3; l++)
//					cout << 36*(l+3*k+6*j+18*i) << endl;
					sv.block<MATSIZE, 1>(36*(l+3*k+6*j+18*i), 0)=kron(Si.col(i), Ii.col(j), Sx.col(k), Ii.col(l))+kron(Sx.col(i), Ii.col(j), Si.col(k), Ii.col(l));
//					sv.block<MATSIZE, 1>(36*(l+3*k+6*j+18*i), 0)=kron(Sx.col(i), Ii.col(j), Sx.col(k), Ii.col(l));

	return sv;
}




// Kronecker delta: 1.0 if a==b, otherwise 0.0.
// Used by: PSMB (this file).
double delta(int a, int b)
// kroneckerProductecker delta-symbol
{
	if (a==b) return 1.0;
	else return 0.0;
}

// Element (i,j) of the commutator superoperator of the 36x36 matrix sm in the 1296-dimensional spin Liouville space; the
// index is composite, i=36*a+b with Hilbert-space indices a and b. It is the sum of sm acting on the second index (when
// the first indices agree) and minus the transpose of sm acting on the first index (when the second indices agree), i.e.
// I x sm - sm^T x I; zero if an index is out of range or neither index pair agrees.
// Used by: InteractionSphericalTensorB::GetElement and PrintSmallMatrix (this file).
std::complex<double> PSMB ( const Matrix<std::complex<double>, 36, 36 > &sm, int i, int j )
{
	// reject indices outside the 1296-dimensional space
	if((i>=MATSIZE*MATSIZE) || (j>=MATSIZE*MATSIZE)) return 0;
	// the two Kronecker deltas tell which of the two Hilbert-space factors is untouched
	double d1=delta(i%MATSIZE, j%MATSIZE), d2=delta(i/MATSIZE, j/MATSIZE);
	if((d1==0) && (d2==0)) return 0;
//	cout << "Here!\n";
//	return sm(i/MATSIZE, j/MATSIZE)*d1 - d2*conj(sm(j%MATSIZE, i%MATSIZE));
// builds I⊗sm - sm'⊗I, the standard commutator-superoperator vec(sm.rho-rho.sm)=(I⊗sm-sm'⊗I)*vec(rho),
// per the column-major Kronecker identity vec(AXB)=(B'⊗A)vec(X); sm' denotes the plain transpose
// (index swap only, no conjugation), which is the correct operator for this identity regardless of
// whether sm is real or complex, so conj() is never applied here; every sm actually passed in (Ag1L,
// AA1L, Ag2L, AA2L, ADL, AWLB) happens to be real-valued (built from Sz, Sp, Sm, Iz, Ip, Im times a
// real B0), so this choice also matches the previous conj()-based formula numerically in all current
// usage, including an arbitrary non-orthorhombic tensor tilt (the tilt-dependent complex Wigner
// coefficients live in Fdu, split into real()/imag() before combining with this result, never in sm)
	return -sm(j/MATSIZE, i/MATSIZE)*d1 + d2*sm(i%MATSIZE, j%MATSIZE);
}

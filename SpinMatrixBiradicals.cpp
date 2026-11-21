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

#include "BiradicalODF.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "SphericalTensor.h"
#include "SpinMatrixB.h"
#define kron(A, B, C, D) kroneckerProduct(A, kroneckerProduct(B, kroneckerProduct(C, D)))
//#define kron(A, B, C, D) kroneckerProduct(kroneckerProduct(kroneckerProduct(A, B), C), D)
#define MATSIZE 36

//extern double B0;

extern Matrix<std::complex<double>, 2, 2> Sx, Sy, Sz, Si, Sp, Sm;
extern Matrix<std::complex<double>, 3, 3> Ix, Iy, Iz, Ii, Ip, Im;


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
	void InteractionSphericalTensorB::PrintSmallMatrix(int l, int m)
	{
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
		if(l==2) cout << l << m << endl << AL2[m+2] << endl << endl;
	}
	std::complex<double> InteractionSphericalTensorB::GetElement(int l, int m, int i, int j)
	{
		if((l<0)||(l>2)) return 0;
		if(abs(m)>l) return 0;
		if((i>=MATSIZE*MATSIZE)||(j>=MATSIZE*MATSIZE)) return 0;
	
		if (l==0) return PSMB(AL0[0], i, j);
		else if(l==2) return PSMB(AL2[m+2], i, j);
		else return 0;
	}
	void InteractionSphericalTensorB::PrintOut()
	{
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

	InteractionSphericalTensorB::InteractionSphericalTensorB()
	{
		AL0=(Matrix<std::complex<double>, 36, 36> *)malloc(sizeof(Matrix<std::complex<double>, 36, 36>)*1);
		AL2=(Matrix<std::complex<double>, 36, 36> *)malloc(sizeof(Matrix<std::complex<double>, 36, 36>)*5);
	
	}
	InteractionSphericalTensorB::~InteractionSphericalTensorB()
	{
		free(AL0);
		free(AL2);
	}




Matrix<std::complex<double>, 36, 36> Ag1L( int l, int m , double B0)
{
	InitializeSpinMatrices();

//	Matrix<std::complex<double>, ColMajor, int> Result(1296, 1296);
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);
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


Matrix<std::complex<double>, 36, 36> AA1L( int l, int m , double B0)
{
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);

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

Matrix<std::complex<double>, 36, 36> Ag2L( int l, int m , double B0)
{
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);
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

Matrix<std::complex<double>, 36, 36> AA2L( int l, int m , double B0)
{
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);

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

Matrix<std::complex<double>, 36, 36> ADL( int l, int m , double B0)
{
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);
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




Matrix<std::complex<double>, 36, 36> AWLB( int l, int m , double B0)
{
	InitializeSpinMatrices();

//	SparseMatrix<std::complex<double>, 1296, 1296> Result=Matrix<std::complex<double>, 1296, 1296>::Zero();
	Matrix<std::complex<double>, 36, 36> result=Matrix<std::complex<double>, 36, 36>::Zero();
	Matrix<std::complex<double>, 36, 36> Ident=kron(Si, Ii, Si, Ii);

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


Matrix<std::complex<double>, 1296, 1> StartingVectorSpinB()
{
	InitializeSpinMatrices();
	Matrix<std::complex<double>, 1296, 1> sv;
//	cout << kroneckerProduct(Sx.col(0), Ii.col(0)).cols() << "x" << kroneckerProduct(Sx.col(0), Ii.col(0)).rows();

    for(int i=0; i<2; i++)
		for(int j=0; j<3; j++)
		    for(int k=0; k<2; k++)
				for(int l=0; l<3; l++)
//					cout << 36*(l+3*k+6*j+18*i) << endl;
					sv.block<MATSIZE, 1>(36*(l+3*k+6*j+18*i), 0)=kron(Si.col(i), Ii.col(j), Sx.col(k), Ii.col(l))+kron(Sx.col(i), Ii.col(j), Si.col(k), Ii.col(l));
//					sv.block<MATSIZE, 1>(36*(l+3*k+6*j+18*i), 0)=kron(Sx.col(i), Ii.col(j), Sx.col(k), Ii.col(l));

	return sv;
}




double delta(int a, int b)
// kroneckerProductecker delta-symbol
{
	if (a==b) return 1.0;
	else return 0.0;
}

std::complex<double> PSMB ( const Matrix<std::complex<double>, 36, 36 > &sm, int i, int j )
{
	if((i>=MATSIZE*MATSIZE) || (j>=MATSIZE*MATSIZE)) return 0;
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

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

// SpinMatrixB.h: declares the spin-operator part of the biradical Spin Hamiltonian in the spherical-tensor form: the class
// InteractionSphericalTensorB (the 36x36 spin matrices A(l,m) of one interaction), the struct SpinInteractionB (one
// interaction: its spin matrices together with its spherical tensor Fdu) and the builders of the A(l,m) matrices. The spin
// space is two electron spins S=1/2 and two nuclear spins I=1 (2*3*2*3=36 states). All definitions are in
// SpinMatrixBiradicals.cpp.
// Used by: BiradicalODF.h (member Interactions of MagneticParameters), LiouvilleMatrixElementB.h and, through it,
// LiouvilleMatrixElementB.cpp and LiouvillianB.cpp; SpinMatrixBiradicals.cpp.
#ifndef SPINMARIXB_H

#define SPINMARIXB_H

#include "SphericalTensor.h"

// InteractionSphericalTensorB: the spin-operator matrices of one interaction, AL0 (one matrix, rank 0) and AL2 (five
// matrices, rank 2, m=-2..2), each 36x36; Set fills them by calling one of the builders declared below, and GetElement
// returns elements of the corresponding Liouville-space superoperator.
class InteractionSphericalTensorB
{
private:
	Matrix<std::complex<double>, 36, 36> *AL0;
	Matrix<std::complex<double>, 36, 36> *AL2;

public:
	void Set(Matrix<std::complex<double>, 36,36> ALFunction (int, int, double), double B0);
//	Matrix<std::complex<double>, 36, 36> *GetMatrix(int l, int m);
	std::complex<double> GetElement(int l, int m, int i, int j);
	void PrintOut();
	void PrintSmallMatrix(int l, int m);

	InteractionSphericalTensorB();
	~InteractionSphericalTensorB();
};


// SpinInteractionB: one term of the Spin Hamiltonian: AuL holds its spin matrices and Fdu its spherical tensor (the
// orientation-dependent part). MagneticParameters::Interactions holds six of them: g1, g2, A1, A2, DJ and the Lorentzian
// line width.
struct SpinInteractionB
{
	InteractionSphericalTensorB AuL;
	SphericalTensor Fdu;
};

// Builders of the spin matrices A(l,m) of the six interactions (definitions in SpinMatrixBiradicals.cpp): electron Zeeman
// term of electron 1 (Ag1L) and 2 (Ag2L), hyperfine coupling of electron 1 with nucleus 1 (AA1L) and electron 2 with
// nucleus 2 (AA2L), electron-electron coupling (ADL) and the Lorentzian line width (AWLB); they are passed to
// InteractionSphericalTensorB::Set as function pointers. StartingVectorSpinB gives the spin part of the starting vector,
// delta is the Kronecker delta and PSMB gives an element of the commutator superoperator of a 36x36 matrix.
Matrix<std::complex<double>, 36, 36> Ag1L( int l, int m , double B0);
Matrix<std::complex<double>, 36, 36> AA1L( int l, int m , double B0);
Matrix<std::complex<double>, 36, 36> AWLB( int l, int m , double B0);
Matrix<std::complex<double>, 36, 36> Ag2L( int l, int m , double B0);
Matrix<std::complex<double>, 36, 36> AA2L( int l, int m , double B0);
Matrix<std::complex<double>, 36, 36> ADL( int l, int m , double B0);
Matrix<std::complex<double>, 1296, 1> StartingVectorSpinB();

double delta(int a, int b);

std::complex<double> PSMB ( const Matrix<std::complex<double>, 36, 36 > &sm, int i, int j );


#endif
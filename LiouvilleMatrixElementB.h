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

// Declarations of the Liouvillian matrix-element functions defined in LiouvilleMatrixElementB.cpp (the 3j symbol Wigner3jMy is defined in
// Wigner3jMy.cpp). Each function gives one element between two basis vectors (BasisVectorB) of the Liouville space of the stochastic Liouville
// equation; see the comments in LiouvilleMatrixElementB.cpp for the meaning of each term.
// Included by: LiouvilleMatrixElementB.cpp, LiouvillianB.cpp (CalculateLiouvillianB uses the matrix elements, SpinIndexB and precalculate_wigner),
// Potential.cpp (XLK_ND) and Wigner3jMy.cpp.
#ifndef LIOUVILLEMATRIXELEMENT_H

#define LIOUVILLEMATRIXELEMENT_H

#include "BiradicalODF.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "SpinMatrixB.h"



// Fills the table of rank-2 reduced Wigner functions for the angle psi (hi, ksi unused). Used by: CalculateLiouvillianB (LiouvillianB.cpp).
void precalculate_wigner(double hi, double psi, double ksi);


// Wigner 3j symbol (L1 L2 L3; M1 M2 M3), defined in Wigner3jMy.cpp. Used by: the matrix-element functions, XLK and XLK_ND in LiouvilleMatrixElementB.cpp.
double Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3);
//extern double Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3);

// Running index of the spin part of a basis vector in the spin Liouville space. Used by: MatrixElementBSK1, CalculateLiouvillianB (starting vector).
int SpinIndexB(BasisVectorB &basis);
//std::complex<double> MatrixElementB   (BasisVectorB &basisr, BasisVectorB &basisc, SpinInteractionB *Interactions, int NInteractions);

// Spin-Hamiltonian element in the K-symmetrised basis. Used by: CalculateLiouvillianB (LiouvillianB.cpp), MatrixElementBSM.
std::complex<double> MatrixElementBSK1(BasisVectorB  basisr, BasisVectorB  basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2);
//std::complex<double> MatrixElementBSM (BasisVectorB  basisr, BasisVectorB  basisc, SpinInteractionB *Interactions, int NInteractions);

// Rotational-diffusion element, M-symmetrised monoradical basis. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementIsotropicBrownianBSM_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, double* Diffusion_Type, double* R_Tau);
// Rotational-diffusion element, K-symmetrised basis, non-diagonal diffusion tensor. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementIsotropicBrownianB_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, double* Diffusion_Type, double* R_Tau);
// Rotational-diffusion element for a diagonal diffusion tensor. Used by: nothing in this code base.
std::complex<double> MatrixElementIsotropicBrownianB(BasisVectorB basisr, BasisVectorB basisc, double Rx, double Ry, double Rz);
//std::complex<double> MatrixElementIsotropicBrownianB(BasisVectorB basisr, BasisVectorB basisc, double Rpar, double Rper);

// Diffusion-in-potential element, M-symmetrised monoradical basis. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementPotentialBSM_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, Potential &pot);
// Diffusion-in-potential element, K-symmetrised basis. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementPotentialB_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, Potential &pot);
// Diffusion-in-potential element for a diagonal diffusion tensor. Used by: nothing in this code base.
std::complex<double> MatrixElementPotentialB(BasisVectorB basisr, BasisVectorB basisc, double Rx, double Ry, double Rz, Potential &pot);

// Coefficient X(L,K) of the diffusion-in-potential superoperator, non-diagonal diffusion tensor. Used by: Potential::CalculateXLK (Potential.cpp).
std::complex<double> XLK_ND(int L, int K, Potential &pot, Matrix<std::complex<double>, 3, 3> &R);
//std::complex<double> XLK_ND(int L, int K, Potential &pot, Matrix<std::complex<double>, 3, 3> &R);
// Same for a diagonal diffusion tensor. Used by: MatrixElementPotentialB, which is itself never called.
double XLK(int L, int K, Potential &pot, double Rx, double Ry, double Rz);

// Spin-Hamiltonian element, M-symmetrised monoradical basis. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementBSM(BasisVectorB basisr, BasisVectorB basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2);

// Declared only: its definition in LiouvilleMatrixElementB.cpp is commented out. Used by: nothing in this code base.
std::complex<double> MatrixElementDiscreteJumpsB(BasisVectorB basisr, BasisVectorB basisc, int NSites, double JumpFrequency);

// Conformer-exchange element. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementConformation(BasisVectorB basisr, BasisVectorB basisc, double* ConversionRate);

// J-diffusion damping element. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementJDiffusion(BasisVectorB basisr, BasisVectorB basisc, double JDiffusionRate);


// Discrete-jump element without ordering potential, monoradical. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementJumpsBSM_Iso(BasisVectorB basisr, BasisVectorB basisc, MagneticParameters *MagneticParameter);
// Discrete-jump element with ordering potential, monoradical. Used by: CalculateLiouvillianB (LiouvillianB.cpp).
std::complex<double> MatrixElementJumpsBSM_Pot(BasisVectorB basisr, BasisVectorB basisc, MagneticParameters* MagneticParameter, Potential &pot);



#endif
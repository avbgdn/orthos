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

#ifndef LIOUVILLEMATRIXELEMENT_H

#define LIOUVILLEMATRIXELEMENT_H

#include "BiradicalODF.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "SpinMatrixB.h"



void precalculate_wigner(double hi, double psi, double ksi);


double Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3);
//extern double Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3);

int SpinIndexB(BasisVectorB &basis);
//std::complex<double> MatrixElementB   (BasisVectorB &basisr, BasisVectorB &basisc, SpinInteractionB *Interactions, int NInteractions);

std::complex<double> MatrixElementBSK1(BasisVectorB  basisr, BasisVectorB  basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2);
//std::complex<double> MatrixElementBSM (BasisVectorB  basisr, BasisVectorB  basisc, SpinInteractionB *Interactions, int NInteractions);

std::complex<double> MatrixElementIsotropicBrownianBSM_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, double* Diffusion_Type, double* R_Tau);
std::complex<double> MatrixElementIsotropicBrownianB_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, double* Diffusion_Type, double* R_Tau);
std::complex<double> MatrixElementIsotropicBrownianB(BasisVectorB basisr, BasisVectorB basisc, double Rx, double Ry, double Rz);
//std::complex<double> MatrixElementIsotropicBrownianB(BasisVectorB basisr, BasisVectorB basisc, double Rpar, double Rper);

std::complex<double> MatrixElementPotentialBSM_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, Potential &pot);
std::complex<double> MatrixElementPotentialB_ND(BasisVectorB basisr, BasisVectorB basisc, Matrix<std::complex<double>, 3, 3> &R, Potential &pot);
std::complex<double> MatrixElementPotentialB(BasisVectorB basisr, BasisVectorB basisc, double Rx, double Ry, double Rz, Potential &pot);

std::complex<double> XLK_ND(int L, int K, Potential &pot, Matrix<std::complex<double>, 3, 3> &R);
//std::complex<double> XLK_ND(int L, int K, Potential &pot, Matrix<std::complex<double>, 3, 3> &R);
double XLK(int L, int K, Potential &pot, double Rx, double Ry, double Rz);

std::complex<double> MatrixElementBSM(BasisVectorB basisr, BasisVectorB basisc, SpinInteractionB *Interactions, int NInteractions, double psi, double J2);

std::complex<double> MatrixElementDiscreteJumpsB(BasisVectorB basisr, BasisVectorB basisc, int NSites, double JumpFrequency);

std::complex<double> MatrixElementConformation(BasisVectorB basisr, BasisVectorB basisc, double* ConversionRate);

std::complex<double> MatrixElementJDiffusion(BasisVectorB basisr, BasisVectorB basisc, double JDiffusionRate);


std::complex<double> MatrixElementJumpsBSM_Iso(BasisVectorB basisr, BasisVectorB basisc, MagneticParameters *MagneticParameter);
std::complex<double> MatrixElementJumpsBSM_Pot(BasisVectorB basisr, BasisVectorB basisc, MagneticParameters* MagneticParameter, Potential &pot);



#endif
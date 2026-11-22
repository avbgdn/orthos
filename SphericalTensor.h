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

// SphericalTensor.h: declares class SphericalTensor, which stores a symmetric 3x3 Cartesian tensor (g, hyperfine, dipolar,
// line width, diffusion) as its spherical components T(l,m) with l=0,2 and rotates it by Euler angles, and declares
// InitializeSpinMatrices, which fills the global spin-operator matrices. Both are defined in SphericalTensor.cpp.
// Used by: BiradicalODF.h (members Fmg1, Fmg2, FmA1, FmA2, FmD, FmW, FmR of MagneticParameters), SpinMatrixB.h (member Fdu
// of SpinInteractionB), and directly by LiouvilleMatrixElementB.cpp, LiouvillianB.cpp, main.cpp, SphericalTensor.cpp and
// SpinMatrixBiradicals.cpp.

#ifndef SPHERICALTENSOR_H

#define SPHERICALTENSOR_H

//#include "BiradicalODF.h"

#include "Eigen/Dense"
#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <time.h>

#ifndef MAX_STRLEN
#define MAX_STRLEN 1000
#endif

using namespace Eigen;
using namespace std;


// fills the global spin matrices Sx..Im; defined in SphericalTensor.cpp, called by the A-matrix builders in SpinMatrixBiradicals.cpp
void InitializeSpinMatrices();

// SphericalTensor: the components T00, T1_1..T11 and T2_2..T22 of a tensor in the spherical basis (T1_1 means T(1,-1)). The
// constructors fill only rank 0 and rank 2 from a symmetric Cartesian matrix; the rank-1 components stay zero. Name is a
// label used when the tensor is printed.
class SphericalTensor
{
// TODO: 
private:
		std::complex<double> T00, T1_1, T10, T11, T2_2, T2_1, T20, T21, T22;
		std::complex<double> I;
public:
	char Name[MAX_STRLEN];
	std::complex<double> Get(int l, int m);
	void Set(Matrix<std::complex<double>, 3, 3> &F);
	void Set(int l, int m, std::complex<double> Value);
	SphericalTensor *Rotate(double EulerAngles[3]);
	void RotateMe(double EulerAngles[3]);
	void PrintOut(std::fstream & tensorfile);
	SphericalTensor(Matrix<std::complex<double>, 3, 3> &F);
	SphericalTensor(Matrix<double, 3, 3> &F);

	SphericalTensor();
	~SphericalTensor();
};



#endif

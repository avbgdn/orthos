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


void InitializeSpinMatrices();

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

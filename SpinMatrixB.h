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

#ifndef SPINMARIXB_H

#define SPINMARIXB_H

#include "SphericalTensor.h"

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


struct SpinInteractionB
{
	InteractionSphericalTensorB AuL;
	SphericalTensor Fdu;
};

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
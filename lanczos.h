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

#ifndef LANCZOS_H

#define LANCZOS_H

//#include "unsupported/Eigen/SparseExtra"

Matrix<std::complex<double>, Dynamic, 2> Lanczos(SparseMatrix<std::complex<double> > &Z,
												 const Matrix<std::complex<double>, Dynamic, 1> &x,
												 int nstep, int ndim,
												 MagneticParameters* MagneticParameter,
												 int &ActualSteps);

Matrix<double, Dynamic, 2> CFVD(Matrix<std::complex<double>, Dynamic, 2> ab, int ndim, double z0, int npoints, double w0, MagneticParameters* MagneticParameter);



#endif

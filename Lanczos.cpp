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

// Lanczos tridiagonalization and Continued fraction decomposition calculation (according to Freed et al.) 
// 

#include "BiradicalODF.h"
#include "Eigen/Sparse"
#include "lanczos.h"
#include <cfloat>
using namespace Eigen;

Matrix<std::complex<double>, Dynamic, 2> Lanczos(SparseMatrix<std::complex<double> > &Z,
												 const Matrix<std::complex<double>, Dynamic, 1> &x,
												 int nstep, int ndim,
												 MagneticParameters* MagneticParameter,
												 int &ActualSteps)
{
	Matrix<std::complex<double>, Dynamic, 1> y, X, temp;
	Matrix<std::complex<double>, Dynamic, 2> ab;

	y.resize(ndim, 1);
	X.resize(ndim, 1);
	temp.resize(ndim, 1);
	ab.resize(ndim, 2);

	ab=Matrix<std::complex<double>, Dynamic, Dynamic>::Zero(ndim, 2);
	y=Matrix<std::complex<double>, Dynamic, Dynamic>::Zero(ndim, 1);

	X=x;
	if(nstep>ndim)
		nstep=ndim;
	ActualSteps=nstep;
	for(int nl=0; nl<nstep; nl++)
	{
		y=Z*X+y;
		ab(nl, 0)=(X.adjoint().conjugate()*y)(0,0);
		y=y-ab(nl, 0)*X;
		std::complex<double> psnorm;
		psnorm=0;
		for(int k=0; k<ndim; k++)
			psnorm+=(y(k,0)*y(k,0));
		ab(nl, 1)=sqrt(psnorm);

		// Residual vanished to machine precision: the Krylov subspace has closed on an exact
		// invariant subspace (the starting vector decomposes into a finite sum of eigenvectors
		// already spanned), so stop before dividing the residual by a near-zero pivot
		double betatol=DBL_EPSILON*std::max(std::abs(ab(nl, 0)), 1.0);
		if(std::abs(ab(nl, 1))<=betatol)
		{
			ab(nl, 1)=0;
			ActualSteps=nl+1;
			break;
		}

		temp=y/ab(nl, 1);
	    y=-X*ab(nl, 1);
		X=temp;
// for debugging purposes, each 10th step of Lanczos we will output the spectrum:
/*
		if(!(nl%10) && (nl>0))
		{
			char buffer[MAX_STRLEN];
			sprintf(buffer, "spec%03d.out", nl/10);
			Matrix<double, Dynamic, 2> Spectrum;
			int npoints=1000;
			Spectrum.resize(npoints, 2);
			Spectrum=CFVD(ab, nl, MagneticParameter->CenterField, MagneticParameter->NPoints, 0, MagneticParameter);

			FILE *SpecOut;	
			SpecOut=fopen(buffer, "w");
			for(int i=0; i<npoints; i++)
				fprintf(SpecOut, "%lf    %e \n", Spectrum(i, 0), Spectrum(i, 1)*2/PI);	
			fclose(SpecOut);
		}
*/
	}
	ab(ActualSteps-1, 1)=0;


	return ab;

}

Matrix<double, Dynamic, 2> CFVD(Matrix<std::complex<double>, Dynamic, 2> ab, int ndim, double z0, int npoints, double w0, MagneticParameters* MagneticParameter)
{
	Matrix<double, Dynamic, 2> Spectrum;
	std::complex<double> z, x, s, tv, td;
	std::complex<double> I(0,1);
	Spectrum.resize(npoints, 2);
//	Spectrum=Matrix<std::complex<double>, Dynamic, 2>::Zero(2*npoints+1, 2);
	for(int iz=0; iz<npoints; iz++)
	{
// z0 = 2*B0 - Field;
	  z=2*z0-MagneticParameter->FieldNet[iz];
//      z=z0-dz*(std::complex<double>)iz/(2.0*(std::complex<double>)npoints);
	  x=w0+I*z;
      s=1.0/(ab(ndim-1, 0)+x);
      tv=s*ab(ndim-2, 1)*ab(ndim-2, 1);
      td=-tv*s;
	  for(int k=ndim-2; k>=1; k--)
	  {
          s=1.0/(ab(k, 0)+x-tv);
          tv=s*ab(k-1, 1)*ab(k-1, 1);
          td=-tv*(1.0-td)*s;
	  }
      s=1.0/(ab(0, 0)+x-tv);
      tv=s;
//	  Spectrum(iz+npoints, 0)=(double)(z0+dz*(std::complex<double>)iz/(2.0*(std::complex<double>)npoints)).real();
	  Spectrum(iz, 0)=MagneticParameter->FieldNet[iz];
	  Spectrum(iz, 1)=(double)(I*tv*(1.0-td)*s).real();
	}

	return Spectrum;
}

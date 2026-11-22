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

// File purpose: Lanczos tridiagonalisation of the Liouvillian of the slow-motion (stochastic Liouville equation)
// EPR calculation, and evaluation of the EPR lineshape from the tridiagonal coefficients as a continued fraction.
// Used by: LiouvillianB.cpp (CalculateLiouvillianB calls Lanczos() and CFVD() through lanczos.h);
// Lanczos.cpp itself only includes lanczos.h and BiradicalODF.h.

#include "BiradicalODF.h"
#include "Eigen/Sparse"
#include "lanczos.h"
#include <cfloat>
using namespace Eigen;

// Lanczos: reduces the complex-symmetric sparse Liouvillian Z to tridiagonal form, starting from the vector x.
// The three-term recursion uses plain (unconjugated) bilinear products, as is appropriate for a complex-symmetric
// matrix, and the recursion assumes that x has unit norm.
// Z: Liouvillian matrix (ndim x ndim); x: normalised starting vector of length ndim;
// nstep: requested number of steps (cut to ndim); MagneticParameter: only used by the commented-out debug output;
// ActualSteps: output, number of steps really done (less than nstep if the Krylov subspace closes early).
// Returns ab (ndim x 2): column 0 holds the diagonal coefficients, column 1 the off-diagonal coefficients
// of the tridiagonal matrix; only the first ActualSteps rows are filled.
// Used by: CalculateLiouvillianB() in LiouvillianB.cpp.
Matrix<std::complex<double>, Dynamic, 2> Lanczos(SparseMatrix<std::complex<double> > &Z,
												 const Matrix<std::complex<double>, Dynamic, 1> &x,
												 int nstep, int ndim,
												 MagneticParameters* MagneticParameter,
												 int &ActualSteps)
{
	Matrix<std::complex<double>, Dynamic, 1> y, X, temp;
	Matrix<std::complex<double>, Dynamic, 2> ab;

	// allocate the work vectors and the output array, and set the output array and the carry-over vector to zero
	y.resize(ndim, 1);
	X.resize(ndim, 1);
	temp.resize(ndim, 1);
	ab.resize(ndim, 2);

	ab=Matrix<std::complex<double>, Dynamic, Dynamic>::Zero(ndim, 2);
	y=Matrix<std::complex<double>, Dynamic, Dynamic>::Zero(ndim, 1);

	// start from the supplied vector and limit the number of steps to the matrix dimension
	X=x;
	if(nstep>ndim)
		nstep=ndim;
	ActualSteps=nstep;
	// Lanczos recursion: each pass gives one diagonal and one off-diagonal coefficient
	for(int nl=0; nl<nstep; nl++)
	{
		// apply the Liouvillian to the current Lanczos vector and add the carry-over from the previous step
		y=Z*X+y;
		ab(nl, 0)=(X.adjoint().conjugate()*y)(0,0);
		y=y-ab(nl, 0)*X;
		// norm of the residual vector, taken as the square root of the bilinear (unconjugated) square sum
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

		// normalise the residual to get the next Lanczos vector and store minus beta times the old one as carry-over
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
	// the last off-diagonal coefficient is not used by the continued fraction, so it is set to zero
	ab(ActualSteps-1, 1)=0;


	return ab;

}

// CFVD: evaluates the EPR lineshape on the field net from the Lanczos coefficients ab, using the continued fraction
// of the resolvent built from the tridiagonal matrix, evaluated bottom-up (from the last step to the first).
// The recursion also carries the derivative of the continued fraction with respect to x (variable td),
// so that the returned second column is the field derivative of the lineshape.
// ab: tridiagonal coefficients from Lanczos(); ndim: number of Lanczos steps to use (the rows ab(0..ndim-1));
// z0: centre field (CenterField); npoints: number of field points; w0: real part of the complex frequency x;
// MagneticParameter: supplies the field net FieldNet[0..npoints-1].
// Returns an npoints x 2 array: column 0 is the field, column 1 the lineshape value.
// Used by: CalculateLiouvillianB() in LiouvillianB.cpp (called with w0=0); also by the commented-out debug
// block inside Lanczos() above.
Matrix<double, Dynamic, 2> CFVD(Matrix<std::complex<double>, Dynamic, 2> ab, int ndim, double z0, int npoints, double w0, MagneticParameters* MagneticParameter)
{
	Matrix<double, Dynamic, 2> Spectrum;
	std::complex<double> z, x, s, tv, td;
	std::complex<double> I(0,1);
	Spectrum.resize(npoints, 2);
//	Spectrum=Matrix<std::complex<double>, Dynamic, 2>::Zero(2*npoints+1, 2);
	// loop over the field points of the spectrum
	for(int iz=0; iz<npoints; iz++)
	{
// z0 = 2*B0 - Field;
	  z=2*z0-MagneticParameter->FieldNet[iz];
//      z=z0-dz*(std::complex<double>)iz/(2.0*(std::complex<double>)npoints);
	  x=w0+I*z;
      s=1.0/(ab(ndim-1, 0)+x);
      tv=s*ab(ndim-2, 1)*ab(ndim-2, 1);
      td=-tv*s;
	  // continued fraction from the second-last step down to step 1, with the derivative recursion alongside
	  for(int k=ndim-2; k>=1; k--)
	  {
          s=1.0/(ab(k, 0)+x-tv);
          tv=s*ab(k-1, 1)*ab(k-1, 1);
          td=-tv*(1.0-td)*s;
	  }
      // top level of the continued fraction, then the lineshape value from the fraction and its derivative
      s=1.0/(ab(0, 0)+x-tv);
      tv=s;
//	  Spectrum(iz+npoints, 0)=(double)(z0+dz*(std::complex<double>)iz/(2.0*(std::complex<double>)npoints)).real();
	  Spectrum(iz, 0)=MagneticParameter->FieldNet[iz];
	  Spectrum(iz, 1)=(double)(I*tv*(1.0-td)*s).real();
	}

	return Spectrum;
}

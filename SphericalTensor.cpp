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

// SphericalTensor.cpp: the class SphericalTensor, which converts a symmetric 3x3 Cartesian interaction tensor (g,
// hyperfine, dipolar, line width, rotational diffusion) into its spherical components T(l,m), gives access to them, rotates
// them by Euler angles and prints them; and the definition of the global electron and nuclear spin matrices (Sx...Im),
// which InitializeSpinMatrices fills.
// Used by: RWMagnetic.cpp (MagneticParameters::RotateTensors builds one SphericalTensor per interaction),
// LiouvilleMatrixElementB.cpp and LiouvillianB.cpp (read the components with Get; LiouvillianB.cpp also prints them),
// SpinMatrixBiradicals.cpp (uses the global spin matrices through extern and calls InitializeSpinMatrices).

#include "BiradicalODF.h"
#include "SphericalTensor.h"

//
// The class for calculation, storage and manipulation of tensors in spherical representation
//

// Global spin-operator matrices in the Zeeman basis: electron spin 1/2 (2x2: Sx, Sy, Sz, identity Si, raising Sp, lowering
// Sm) and nuclear spin 1 (3x3: Ix, Iy, Iz, identity Ii, raising Ip, lowering Im). They are filled by InitializeSpinMatrices.
// Used by: SpinMatrixBiradicals.cpp (declared extern there).
Matrix<std::complex<double>, 2, 2> Sx, Sy, Sz, Si, Sp, Sm;
Matrix<std::complex<double>, 3, 3> Ix, Iy, Iz, Ii, Ip, Im;

// Returns the spherical component T(l,m) (l=0,1,2, |m|<=l); returns 0 for any other (l,m). The rank-1 components are not
// set by the constructors, so they are zero unless set with Set(l,m,Value).
// Used by: Rotate, RotateMe and PrintOut (this file); MatrixElementBSK1 (LiouvilleMatrixElementB.cpp);
// CalculateLiouvillianB (LiouvillianB.cpp, symmetry tests of the tensors and the isotropic hyperfine value).
std::complex<double> SphericalTensor::Get(int l, int m)
	{
		// look up the member that holds the requested component
		if ((l==0)&&(m== 0)) return T00;
		if ((l==1)&&(m==-1)) return T1_1;
		if ((l==1)&&(m== 0)) return T10;
		if ((l==1)&&(m== 1)) return T11;
		if ((l==2)&&(m==-2)) return T2_2;
		if ((l==2)&&(m==-1)) return T2_1;
		if ((l==2)&&(m== 0)) return T20;
		if ((l==2)&&(m== 1)) return T21;
		if ((l==2)&&(m== 2)) return T22;
		return 0.0;
	}
	// Sets the rank-0 and rank-2 components from a 3x3 complex Cartesian tensor F with the same formulae as the constructors
	// (the rank-1 components are left unchanged).
	// Used by: nothing in this code base.
	void SphericalTensor::Set(Matrix<std::complex<double>, 3, 3> &F)
	{
		// rank-0 and rank-2 components from the Cartesian elements of F
		T00 =-(F(0,0)+F(1,1)+F(2,2))/sqrt(3.0);
		T20 =sqrt(2.0/3.0)*(F(2,2)-0.5*(F(0,0)+F(1,1)));
		T21 =-(F(2,0)+I*F(2,1));
		T2_1= (F(2,0)-I*F(2,1));
        T22 =0.5*(F(0,0)-F(1,1)+2.0*I*F(0,1));
        T2_2=0.5*(F(0,0)-F(1,1)-2.0*I*F(0,1));
	}
	// Sets the single component T(l,m) to Value; (l,m) outside l=0,1,2, |m|<=l is ignored.
	// Used by: Rotate and RotateMe (this file).
	void SphericalTensor::Set(int l, int m, std::complex<double> Value)
	{
		// store the value in the member that holds the (l,m) component
		if ((l==0)&&(m== 0)) T00=Value;
		if ((l==1)&&(m==-1)) T1_1=Value;
		if ((l==1)&&(m== 0)) T10=Value;
		if ((l==1)&&(m== 1)) T11=Value;
		if ((l==2)&&(m==-2)) T2_2=Value;
		if ((l==2)&&(m==-1)) T2_1=Value;
		if ((l==2)&&(m== 0)) T20=Value;
		if ((l==2)&&(m== 1)) T21=Value;
		if ((l==2)&&(m== 2)) T22=Value;
	}

	// Rotates this tensor in place by the Euler angles (alpha, beta, gamma in radians): the new component T'(l,m) is the sum
	// over k of T(l,k)*D^l_km(alpha,beta,gamma), for l=0,1,2.
	// Used by: nothing in this code base.
	void SphericalTensor::RotateMe(double EulerAngles[3])
	{
// Euler angles, which convert initial coordinate system to final;
// Formulae are taken from Schneider, Freed, BMR.
		SphericalTensor *Result=new SphericalTensor();
//		std::complex<double> T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
        // rotate each rank separately: every new component is a sum of old components times a Wigner function
        for(int l=0; l<=2; l++)
			for(int m=-l; m<=l; m++)
			{
				std::complex<double> result=0;
				for(int k=-l; k<=l; k++)
				{
				    result=result+this->Get(l, k)*ComplexWigner(l, k, m, EulerAngles[0], EulerAngles[1], EulerAngles[2]);
					// NOTE that this result is COMPLEX conjugate with respect to MatLab Fdg, FdA, FdW formulae
					// But, this is correct rotation, cf Freed, Schneider (1989), eq. (14).
				}
				Result->Set(l, m, result);
			}
        // copy the rotated components back into this object
        for(int l=0; l<=2; l++)
			for(int m=-l; m<=l; m++)
			{
				this->Set(l, m, Result->Get(l, m));
			}
	}


	// Returns a pointer to a new SphericalTensor that holds this tensor rotated by the Euler angles (alpha, beta, gamma in
	// radians), T'(l,m)=sum over k of T(l,k)*D^l_km; this tensor is unchanged and the caller owns the returned object.
	// Used by: nothing in this code base.
	SphericalTensor * SphericalTensor::Rotate(double EulerAngles[3])
	{
// Euler angles, which convert initial coordinate system to final;
// Formulae are taken from Schneider, Freed, BMR.
		SphericalTensor *Result=new SphericalTensor();
//		std::complex<double> T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
        // rotate each rank separately: every new component is a sum of old components times a Wigner function
        for(int l=0; l<=2; l++)
			for(int m=-l; m<=l; m++)
			{
				std::complex<double> result=0;
				for(int k=-l; k<=l; k++)
				{
				    result=result+this->Get(l, k)*ComplexWigner(l, k, m, EulerAngles[0], EulerAngles[1], EulerAngles[2]);
					// NOTE that this result is COMPLEX conjugate with respect to MatLab Fdg, FdA, FdW formulae
					// But, this is correct rotation, cf Freed, Schneider (1989), eq. (14).
				}
				Result->Set(l, m, result);
			}
		return Result;

	}
	// Writes the Name of the tensor and its components T(l,m), l=0..2, to the stream tensorfile, one line per rank (diagnostics).
	// Used by: CalculateLiouvillianB (LiouvillianB.cpp), which prints the six interaction tensors to its tensors<N>_<rank>.out file.
	void SphericalTensor::PrintOut(std::fstream &tensorfile)
	{
		tensorfile << Name << " spherical tensorial components: \n";
		// one line per rank, with all its components
		for (int l=0; l<=2; l++)
		{
			for(int m=-l; m<=l; m++)
			{
				char buffer[MAX_STRLEN];
				sprintf(buffer, "T(%2d, %2d)=", l, m);
				tensorfile << buffer;
				tensorfile << this->Get(l, m) << "\t";
			}
			tensorfile<<endl;
		}
		tensorfile<<endl;
	}

	// Constructs the tensor from a symmetric 3x3 complex Cartesian matrix F: T(0,0)=-trace(F)/sqrt(3) and the rank-2
	// components T(2,m) from the elements of F (formulae of Varshalovich, see below); the rank-1 components are zero, which
	// is why F must be symmetric.
	// Used by: MagneticParameters::RotateTensors (RWMagnetic.cpp), for the rotational diffusion tensor DiffusionR.
	SphericalTensor::SphericalTensor(Matrix<std::complex<double>, 3, 3> &F)
	{
// Only works correctly for symmetric tensors!
		// imaginary unit and zero initial components
		I=std::complex<double>(0, 1);
		T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
// Formulae from Varshalovich		
		T00 =-(F(0,0)+F(1,1)+F(2,2))/sqrt(3.0);
		T20 =sqrt(2.0/3.0)*(F(2,2)-0.5*(F(0,0)+F(1,1)));
		T21 =-(F(2,0)+I*F(2,1));
		T2_1= (F(2,0)-I*F(2,1));
        T22 =0.5*(F(0,0)-F(1,1)+2.0*I*F(0,1));
        T2_2=0.5*(F(0,0)-F(1,1)-2.0*I*F(0,1));

	}

	// Constructs the tensor from a symmetric 3x3 real Cartesian matrix F, with the same formulae as the complex version.
	// Used by: MagneticParameters::RotateTensors (RWMagnetic.cpp), for the rotated g1, g2, A1, A2 tensors, the dipolar and
	// exchange tensor DJ and the Lorentzian line width tensor.
	SphericalTensor::SphericalTensor(Matrix<double, 3, 3> &F)
	{
// Only works correctly for symmetric tensors!
		// imaginary unit and zero initial components
		I=std::complex<double>(0, 1);
		T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
// Formulae from Varshalovich		
		T00 =-(F(0,0)+F(1,1)+F(2,2))/sqrt(3.0);
		T20 =sqrt(2.0/3.0)*(F(2,2)-0.5*(F(0,0)+F(1,1)));
		T21 =-(F(2,0)+I*F(2,1));
		T2_1= (F(2,0)-I*F(2,1));
        T22 =0.5*(F(0,0)-F(1,1)+2.0*I*F(0,1));
        T2_2=0.5*(F(0,0)-F(1,1)-2.0*I*F(0,1));

	}


	// Constructs a tensor with all components zero.
	// Used by: Rotate and RotateMe (this file) for the result object, and implicitly for the member Fdu of every
	// SpinInteractionB (SpinMatrixB.h).
	SphericalTensor::SphericalTensor()
	{
		// imaginary unit and zero components
		I=std::complex<double>(0, 1);
		T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
	}
	// Destructor: nothing to release.
	// Used by: implicit destruction of SphericalTensor objects.
	SphericalTensor::~SphericalTensor()
	{}



// Fills the global spin-operator matrices in the Zeeman basis: for the electron spin S=1/2 the 2x2 matrices Sx, Sy, Sz,
// the identity Si and the raising and lowering matrices Sp, Sm; for the nuclear spin I=1 the 3x3 matrices Ix, Iy, Iz, the
// identity Ii and Ip, Im (Ip and Im carry the spin-1 factor sqrt(2), the states are ordered m=+1, 0, -1).
// Used by: Ag1L, AA1L, Ag2L, AA2L, ADL, AWLB and StartingVectorSpinB (SpinMatrixBiradicals.cpp), each of which calls it
// first.
// Initializing Unitary matrices for Spin Hamiltonian
void InitializeSpinMatrices()
{
	// imaginary unit
	std::complex<double> I(0, 1);

	// electron spin 1/2 operators
	Sx << 0, 0.5,
		0.5, 0;

	Sy << 0, 0.5*I,
		-0.5*I, 0;

	Sz << 0.5, 0,
		0, -0.5;

	Sp << 0, 1,
		  0, 0;
	
	Sm << 0, 0,
		  1, 0;

	Si << 1, 0,
		0, 1;

	// nuclear spin 1 operators
	Ix << 0, 1, 0,
		  1, 0, 1,
		  0, 1, 0;
	Ix=Ix*sqrt(2.0);

	Iy << 0, I, 0,
		 -I, 0, I,
		  0,-I, 0;
	Iy=Iy*sqrt(2.0);

	Iz << 1, 0, 0,
		  0, 0, 0,
		  0, 0,-1;

	Ii << 1, 0, 0,
		  0, 1, 0,
		  0, 0, 1;

	Ip << 0, 1, 0,
		  0, 0, 1,
		  0, 0, 0;
	Ip=Ip*sqrt(2.0);

	Im << 0, 0, 0,
		  1, 0, 0,
		  0, 1, 0;
	Im=Im*sqrt(2.0);

}
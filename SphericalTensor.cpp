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

#include "BiradicalODF.h"
#include "SphericalTensor.h"

//
// The class for calculation, storage and manipulation of tensors in spherical representation
//

Matrix<std::complex<double>, 2, 2> Sx, Sy, Sz, Si, Sp, Sm;
Matrix<std::complex<double>, 3, 3> Ix, Iy, Iz, Ii, Ip, Im;

std::complex<double> SphericalTensor::Get(int l, int m)
	{
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
	void SphericalTensor::Set(Matrix<std::complex<double>, 3, 3> &F)
	{
		T00 =-(F(0,0)+F(1,1)+F(2,2))/sqrt(3.0);
		T20 =sqrt(2.0/3.0)*(F(2,2)-0.5*(F(0,0)+F(1,1)));
		T21 =-(F(2,0)+I*F(2,1));
		T2_1= (F(2,0)-I*F(2,1));
        T22 =0.5*(F(0,0)-F(1,1)+2.0*I*F(0,1));
        T2_2=0.5*(F(0,0)-F(1,1)-2.0*I*F(0,1));
	}
	void SphericalTensor::Set(int l, int m, std::complex<double> Value)
	{
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

	void SphericalTensor::RotateMe(double EulerAngles[3])
	{
// Euler angles, which convert initial coordinate system to final;
// Formulae are taken from Schneider, Freed, BMR.
		SphericalTensor *Result=new SphericalTensor();
//		std::complex<double> T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
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
        for(int l=0; l<=2; l++)
			for(int m=-l; m<=l; m++)
			{
				this->Set(l, m, Result->Get(l, m));
			}
	}


	SphericalTensor * SphericalTensor::Rotate(double EulerAngles[3])
	{
// Euler angles, which convert initial coordinate system to final;
// Formulae are taken from Schneider, Freed, BMR.
		SphericalTensor *Result=new SphericalTensor();
//		std::complex<double> T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
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
	void SphericalTensor::PrintOut(std::fstream &tensorfile)
	{
		tensorfile << Name << " spherical tensorial components: \n";
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

	SphericalTensor::SphericalTensor(Matrix<std::complex<double>, 3, 3> &F)
	{
// Only works correctly for symmetric tensors!
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

	SphericalTensor::SphericalTensor(Matrix<double, 3, 3> &F)
	{
// Only works correctly for symmetric tensors!
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


	SphericalTensor::SphericalTensor()
	{
		I=std::complex<double>(0, 1);
		T00=0, T1_1=0, T10=0, T11=0, T2_2=0, T2_1=0, T20=0, T21=0, T22=0;
	}
	SphericalTensor::~SphericalTensor()
	{}



// Initializing Unitary matrices for Spin Hamiltonian
void InitializeSpinMatrices()
{
	std::complex<double> I(0, 1);

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
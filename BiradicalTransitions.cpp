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
#include "unsupported/Eigen/KroneckerProduct"

extern int MinimizationFlag;

// The function generates the array of transitions (pair Field-Intensity)
// for rigid biradical with magnetic parameters given in the object MagneticParameter
// with orientation defined by rotation matrix Field_to_Molecule_RotMat.
// Also the values of effective Gaussian and Lorentzian linewidths corresponding to this orientation are calculated.
Array<double, 54, 2> CalculateTransitions(MagneticParameters *MagneticParameter,
										  Matrix<double, 3, 3> Field_to_Molecule_RotMat,
										  double *HGaussEff, double *HLorentzEff
										  )
{
// Initializing Magnetic parameters
	Matrix<double, 3, 3> g1M, g2M, A1M, A2M, HGaussM, HLorentzM;
	Matrix<double, 3, 3>  DJ;
/*
	double dx, dy, dz, dr;
	dx=MagneticParameter->dx;
	dy=MagneticParameter->dy;
	dz=MagneticParameter->dz;
*/
	FILE *transitions;
	if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) transitions=fopen("Transitions.out", "a+");

/*
// Computing fine splitting tensor in point dipole approximation;
    dr=sqrt(dx*dx+dy*dy+dz*dz);
    D << dr*dr-3*dx*dx, -3*dx*dy, -3*dx*dz,
		-3*dx*dy, dr*dr-3*dy*dy, -3*dy*dz,
		-3*dx*dz, -3*dy*dz, dr*dr-3*dz*dz;

	D=D*18561.0*pow(dr, -5);
    J << 1, 0, 0,
		 0, 1, 0,
		 0, 0, 1;
	J=J*MagneticParameter->J0;
	DJ=D+J;
// Rotating DJ tensor as it is coaxial with g1 tensor
	Matrix<double, 3, 3> RotMat;
    RotMat=RotationMatrix(MagneticParameter->g1Euler);
	DJ=RotMat*DJ*RotMat.transpose();
*/
	DJ=MagneticParameter->DJ;

// Initializing Array for Transition intensities
	Array<double, 54, 2> Transitions=Array<double, 54, 2>::Zero();
	int r=0;
// Initializing Magnetic Field related values;
	Matrix<double, 1, 3> B0, B1, nv;

	B1 << 0, 1, 0;
	nv << 0, 0, 1;
	B0 = MagneticParameter->MagneticParameterArray[0].CenterField/sqrt( (MagneticParameter->MagneticParameterArray[0].g1R.transpose() * MagneticParameter->MagneticParameterArray[0].g1R).trace() / 3.0 ) * nv;
//	B0 = MagneticParameter->CenterField/sqrt( (MagneticParameter->g1R.transpose() * MagneticParameter->g1R).trace() / 3.0 ) * nv;

// Rotation of all tensors to Resonator frame.
    Matrix<double, 3, 3> Field_to_Molecule_RotMat_T;
    Field_to_Molecule_RotMat_T=Field_to_Molecule_RotMat.transpose();
    
	g1M=Field_to_Molecule_RotMat_T*MagneticParameter->g1R*Field_to_Molecule_RotMat;
	g2M=Field_to_Molecule_RotMat_T*MagneticParameter->g2R*Field_to_Molecule_RotMat;
	A1M=Field_to_Molecule_RotMat_T*MagneticParameter->A1R*Field_to_Molecule_RotMat;
	A2M=Field_to_Molecule_RotMat_T*MagneticParameter->A2R*Field_to_Molecule_RotMat;
	DJ=Field_to_Molecule_RotMat_T*DJ*Field_to_Molecule_RotMat;
	HGaussM=Field_to_Molecule_RotMat_T*MagneticParameter->HGaussR*Field_to_Molecule_RotMat;
	HLorentzM=Field_to_Molecule_RotMat_T*MagneticParameter->HLorentzR*Field_to_Molecule_RotMat;


// Constructing effective A-values;
	double Aeff1, Aeff2, mz1, mz2;
			
	Aeff1=sqrt(nv*g1M*(A1M*A1M.transpose())*g1M.transpose()*nv.transpose())/sqrt((g1M.transpose()*g1M).trace()/3);
	Aeff2=sqrt(nv*g2M*(A2M*A2M.transpose())*g2M.transpose()*nv.transpose())/sqrt((g2M.transpose()*g2M).trace()/3);

	(*HGaussEff)=sqrt(nv*g1M*(HGaussM*HGaussM.transpose())*g1M.transpose()*nv.transpose())/sqrt((g1M.transpose()*g1M).trace()/3);
	(*HLorentzEff)=sqrt(nv*g1M*(HLorentzM*HLorentzM.transpose())*g1M.transpose()*nv.transpose())/sqrt((g1M.transpose()*g1M).trace()/3);

// Initializing Unitary matrices for Spin Hamiltonian
	Matrix<std::complex<double>, 2, 2> Sx, Sy, Sz, Si;
	std::complex<double> I(0, 1);

	Sx << 0, 0.5,
		0.5, 0;

	Sy << 0, 0.5*I,
		-0.5*I, 0;

	Sz << 0.5, 0,
		0, -0.5;

	Si << 1, 0,
		0, 1;

// Construction of SpinHamiltonian:
	Matrix<std::complex<double>, 4, 4> SpinHamiltonian, H11, H12, H21, H22, H3, Hp, Intensity;
	H11=Matrix<std::complex<double>, 4, 4>::Zero();
	H12=Matrix<std::complex<double>, 4, 4>::Zero();
	H21=Matrix<std::complex<double>, 4, 4>::Zero();
	H22=Matrix<std::complex<double>, 4, 4>::Zero();
	H3 =Matrix<std::complex<double>, 4, 4>::Zero();

    H11=    (B0(0)*g1M(0, 0)+B0(1)*g1M(1, 0)+B0(2)*g1M(2, 0))*kroneckerProduct(Sx, Si);
    H11=H11+(B0(0)*g1M(0, 1)+B0(1)*g1M(1, 1)+B0(2)*g1M(2, 1))*kroneckerProduct(Sy, Si);
    H11=H11+(B0(0)*g1M(0, 2)+B0(1)*g1M(1, 2)+B0(2)*g1M(2, 2))*kroneckerProduct(Sz, Si);

    H12=    (B0(0)*g2M(0, 0)+B0(1)*g2M(1, 0)+B0(2)*g2M(2, 0))*kroneckerProduct(Si, Sx);
    H12=H12+(B0(0)*g2M(0, 1)+B0(1)*g2M(1, 1)+B0(2)*g2M(2, 1))*kroneckerProduct(Si, Sy);
    H12=H12+(B0(0)*g2M(0, 2)+B0(1)*g2M(1, 2)+B0(2)*g2M(2, 2))*kroneckerProduct(Si, Sz);

    H3=   DJ(0, 0)*kroneckerProduct(Sx, Sx);
    H3=H3+DJ(0, 1)*kroneckerProduct(Sx, Sy);
    H3=H3+DJ(0, 2)*kroneckerProduct(Sx, Sz);

    H3=H3+DJ(1, 0)*kroneckerProduct(Sy, Sx);
    H3=H3+DJ(1, 1)*kroneckerProduct(Sy, Sy);
    H3=H3+DJ(1, 2)*kroneckerProduct(Sy, Sz);

    H3=H3+DJ(2, 0)*kroneckerProduct(Sz, Sx);
    H3=H3+DJ(2, 1)*kroneckerProduct(Sz, Sy);
    H3=H3+DJ(2, 2)*kroneckerProduct(Sz, Sz);

    Hp= (B1(0)*g1M(0, 0)+B1(1)*g1M(1, 0)+B1(2)*g1M(2, 0))*kroneckerProduct(Sx, Si)+
		(B1(0)*g1M(0, 1)+B1(1)*g1M(1, 1)+B1(2)*g1M(2, 1))*kroneckerProduct(Sy, Si)+
		(B1(0)*g1M(0, 2)+B1(1)*g1M(1, 2)+B1(2)*g1M(2, 2))*kroneckerProduct(Sz, Si)+
		(B1(0)*g2M(0, 0)+B1(1)*g2M(1, 0)+B1(2)*g2M(2, 0))*kroneckerProduct(Si, Sx)+
		(B1(0)*g2M(0, 1)+B1(1)*g2M(1, 1)+B1(2)*g2M(2, 1))*kroneckerProduct(Si, Sy)+
		(B1(0)*g2M(0, 2)+B1(1)*g2M(1, 2)+B1(2)*g2M(2, 2))*kroneckerProduct(Si, Sz);

	if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) 
		fprintf(transitions, "   mz1  mz2  ms1i  ms1f  ms2i  ms2f    Field    Intensity \n");

	for(mz1=-1; mz1<=1; mz1++)
		for(mz2=-1; mz2<=1; mz2++)
		{
 			H21=(Aeff1*mz1)*kroneckerProduct(Sz, Si);
			H22=(Aeff2*mz2)*kroneckerProduct(Si, Sz);
			SpinHamiltonian=H11+H12+H21+H22+H3;

// Computing eigenvalues:
			SelfAdjointEigenSolver<Matrix<std::complex<double>, 4, 4> > es;
			es.compute(SpinHamiltonian);
			Intensity=es.eigenvectors().adjoint() * Hp * es.eigenvectors();
    
			for (int k=0; k<4; k++)
			{
				for (int j=k+1; j<4; j++)
				{
					double ms1i, ms2i, ms1f, ms2f;
					ms1i=k/2;			ms1i=0.5-ms1i;
					ms2i=k%2;			ms2i=0.5-ms2i;

					ms1f=j/2;			ms1f=0.5-ms1f;
					ms2f=j%2;			ms2f=0.5-ms2f;

					Transitions(r, 0)=2*MagneticParameter->CenterField-(es.eigenvalues()(j)-es.eigenvalues()(k));
					Transitions(r, 1)=pow(abs(Intensity(k, j)), 2);
					if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) 
						fprintf(transitions, "  %2d  %2d  %+4.1lf   %+4.1lf  %+4.1lf  %+4.1lf  %8.3lf   %lf\n", (int)mz1, (int)mz2, ms1i, ms1f, ms2i, ms2f, Transitions(r, 0), Transitions(r, 1));
					r++;
				}
			}
//	std::cout << "SpinHamiltonian:\n" << SpinHamiltonian << std::endl << std::endl;
//	std::cout << "Perturbation Hamiltonian:\n" << Hp << std::endl << std::endl;
//	std::cout << "Eigenvalues:\n" << es.eigenvalues() << std::endl << std::endl;
//	std::cout << "Eigenvectors:\n" << es.eigenvectors() << std::endl << std::endl;
//	std::cout << "Intensities:\n" << Intensity << std::endl << std::endl;
		}
		if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) fclose(transitions);
  return Transitions;

}



//Calculation of transitions for Monoradical
Array<double, 3, 2> CalculateTransitionsM(MagneticParameters *MagneticParameter,
										   Matrix<double, 3, 3> Field_to_Molecule_RotMat,
										   double *HGaussEff, double *HLorentzEff
										   )
{
// Initializing Magnetic parameters
	Matrix<double, 3, 3> g1M, g2M, A1M, A2M, HGaussM, HLorentzM;
	Matrix<double, 3, 3>  DJ;
/*
	double dx, dy, dz, dr;
	dx=MagneticParameter->dx;
	dy=MagneticParameter->dy;
	dz=MagneticParameter->dz;
*/
	FILE *transitions;
	if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) transitions=fopen("Transitions.out", "a+");

/*
// Computing fine splitting tensor in point dipole approximation;
    dr=sqrt(dx*dx+dy*dy+dz*dz);
    D << dr*dr-3*dx*dx, -3*dx*dy, -3*dx*dz,
		-3*dx*dy, dr*dr-3*dy*dy, -3*dy*dz,
		-3*dx*dz, -3*dy*dz, dr*dr-3*dz*dz;

	D=D*18561.0*pow(dr, -5);
    J << 1, 0, 0,
		 0, 1, 0,
		 0, 0, 1;
	J=J*MagneticParameter->J0;
	DJ=D+J;
// Rotating DJ tensor as it is coaxial with g1 tensor
	Matrix<double, 3, 3> RotMat;
    RotMat=RotationMatrix(MagneticParameter->g1Euler);
	DJ=RotMat*DJ*RotMat.transpose();
*/
	DJ=MagneticParameter->DJ;

// Initializing Array for Transition intensities
	Array<double, 3, 2> Transitions=Array<double, 3, 2>::Zero();
	int r=0;
// Initializing Magnetic Field related values;
	Matrix<double, 1, 3> B0, B1, nv;

	B1 << 0, 1, 0;
	nv << 0, 0, 1;

	B0 = MagneticParameter->MagneticParameterArray[0].CenterField/sqrt( (MagneticParameter->MagneticParameterArray[0].g1R.transpose() * MagneticParameter->MagneticParameterArray[0].g1R).trace() / 3.0 ) * nv;
//    cout << B0 << endl;
	//	B0 = MagneticParameter->CenterField/sqrt( (MagneticParameter->g1R.transpose() * MagneticParameter->g1R).trace() / 3.0 ) * nv;

// Rotation of all tensors to Resonator frame.
    Matrix<double, 3, 3> Field_to_Molecule_RotMat_T;
    Field_to_Molecule_RotMat_T=Field_to_Molecule_RotMat.transpose();
    
	g1M=Field_to_Molecule_RotMat_T*MagneticParameter->g1R*Field_to_Molecule_RotMat;
//	g2M=Field_to_Molecule_RotMat_T*MagneticParameter->g2R*Field_to_Molecule_RotMat;
	A1M=Field_to_Molecule_RotMat_T*MagneticParameter->A1R*Field_to_Molecule_RotMat;
//	A2M=Field_to_Molecule_RotMat_T*MagneticParameter->A2R*Field_to_Molecule_RotMat;
//	DJ=Field_to_Molecule_RotMat_T*DJ*Field_to_Molecule_RotMat;
	HGaussM=Field_to_Molecule_RotMat_T*MagneticParameter->HGaussR*Field_to_Molecule_RotMat;
	HLorentzM=Field_to_Molecule_RotMat_T*MagneticParameter->HLorentzR*Field_to_Molecule_RotMat;


// Constructing effective A-values;
	double Aeff1, Aeff2, mz1, mz2;
			
	Aeff1=sqrt(nv*g1M*(A1M*A1M.transpose())*g1M.transpose()*nv.transpose())/sqrt((g1M.transpose()*g1M).trace()/3);
//	Aeff2=sqrt(nv*g2M*(A2M*A2M.transpose())*g2M.transpose()*nv.transpose())/sqrt((g2M.transpose()*g2M).trace()/3);

	(*HGaussEff)=sqrt(nv*g1M*(HGaussM*HGaussM.transpose())*g1M.transpose()*nv.transpose())/sqrt((g1M.transpose()*g1M).trace()/3);
	(*HLorentzEff)=sqrt(nv*g1M*(HLorentzM*HLorentzM.transpose())*g1M.transpose()*nv.transpose())/sqrt((g1M.transpose()*g1M).trace()/3);

// Initializing Unitary matrices for Spin Hamiltonian
	Matrix<std::complex<double>, 2, 2> Sx, Sy, Sz, Si;
	std::complex<double> I(0, 1);

	Sx << 0, 0.5,
		0.5, 0;

	Sy << 0, 0.5*I,
		-0.5*I, 0;

	Sz << 0.5, 0,
		0, -0.5;

	Si << 1, 0,
		0, 1;

// Construction of SpinHamiltonian:
	Matrix<std::complex<double>, 2, 2> SpinHamiltonian, H11, H21, Hp, Intensity;
	H11=Matrix<std::complex<double>, 2, 2>::Zero();
	H21=Matrix<std::complex<double>, 2, 2>::Zero();

    H11=    (B0(0)*g1M(0, 0)+B0(1)*g1M(1, 0)+B0(2)*g1M(2, 0))*Sx;
    H11=H11+(B0(0)*g1M(0, 1)+B0(1)*g1M(1, 1)+B0(2)*g1M(2, 1))*Sy;
    H11=H11+(B0(0)*g1M(0, 2)+B0(1)*g1M(1, 2)+B0(2)*g1M(2, 2))*Sz;
/*
    H12=    (B0(0)*g2M(0, 0)+B0(1)*g2M(1, 0)+B0(2)*g2M(2, 0))*kroneckerProduct(Si, Sx);
    H12=H12+(B0(0)*g2M(0, 1)+B0(1)*g2M(1, 1)+B0(2)*g2M(2, 1))*kroneckerProduct(Si, Sy);
    H12=H12+(B0(0)*g2M(0, 2)+B0(1)*g2M(1, 2)+B0(2)*g2M(2, 2))*kroneckerProduct(Si, Sz);

    H3=   DJ(0, 0)*kroneckerProduct(Sx, Sx);
    H3=H3+DJ(0, 1)*kroneckerProduct(Sx, Sy);
    H3=H3+DJ(0, 2)*kroneckerProduct(Sx, Sz);

    H3=H3+DJ(1, 0)*kroneckerProduct(Sy, Sx);
    H3=H3+DJ(1, 1)*kroneckerProduct(Sy, Sy);
    H3=H3+DJ(1, 2)*kroneckerProduct(Sy, Sz);

    H3=H3+DJ(2, 0)*kroneckerProduct(Sz, Sx);
    H3=H3+DJ(2, 1)*kroneckerProduct(Sz, Sy);
    H3=H3+DJ(2, 2)*kroneckerProduct(Sz, Sz);
*/
    Hp= (B1(0)*g1M(0, 0)+B1(1)*g1M(1, 0)+B1(2)*g1M(2, 0))*Sx+
		(B1(0)*g1M(0, 1)+B1(1)*g1M(1, 1)+B1(2)*g1M(2, 1))*Sy+
		(B1(0)*g1M(0, 2)+B1(1)*g1M(1, 2)+B1(2)*g1M(2, 2))*Sz;

	if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) 
		fprintf(transitions, "   mz1  mz2  ms1i  ms1f  ms2i  ms2f    Field    Intensity \n");

	for(mz1=-1; mz1<=1; mz1++)
//		for(mz2=-1; mz2<=1; mz2++)
		{
 			H21=(Aeff1*mz1)*Sz;
//			H22=(Aeff2*mz2)*kroneckerProduct(Si, Sz);
			SpinHamiltonian=H11+H21;

// Computing eigenvalues:
			SelfAdjointEigenSolver<Matrix<std::complex<double>, 2, 2> > es;
			es.compute(SpinHamiltonian);
			Intensity=es.eigenvectors().adjoint() * Hp * es.eigenvectors();
    
			for (int k=0; k<2; k++)
			{
				for (int j=k+1; j<2; j++)
				{
					double ms1i, ms2i, ms1f, ms2f;
					ms1i=k/2;			ms1i=0.5-ms1i;
					ms2i=k%2;			ms2i=0.5-ms2i;

					ms1f=j/2;			ms1f=0.5-ms1f;
					ms2f=j%2;			ms2f=0.5-ms2f;

					Transitions(r, 0)=2*MagneticParameter->CenterField-(es.eigenvalues()(j)-es.eigenvalues()(k));
					Transitions(r, 1)=pow(abs(Intensity(k, j)), 2)*12.0;
					if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) 
						fprintf(transitions, "  %2d  %2d  %+4.1lf   %+4.1lf  %+4.1lf  %+4.1lf  %8.3lf   %lf\n", (int)mz1, (int)mz2, ms1i, ms1f, ms2i, ms2f, Transitions(r, 0), Transitions(r, 1));
					r++;
				}
			}
//	std::cout << "SpinHamiltonian:\n" << SpinHamiltonian << std::endl << std::endl;
//	std::cout << "Perturbation Hamiltonian:\n" << Hp << std::endl << std::endl;
//	std::cout << "Eigenvalues:\n" << es.eigenvalues() << std::endl << std::endl;
//	std::cout << "Eigenvectors:\n" << es.eigenvectors() << std::endl << std::endl;
//	std::cout << "Intensities:\n" << Intensity << std::endl << std::endl;
		}
		if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) fclose(transitions);
  return Transitions;

}

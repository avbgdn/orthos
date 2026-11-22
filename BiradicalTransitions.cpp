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

// File purpose: rigid-limit EPR transitions of one fixed orientation: diagonalises the spin Hamiltonian of a rigid
// biradical (CalculateTransitions) or monoradical (CalculateTransitionsM) and returns the resonance fields and
// intensities of all allowed transitions together with the orientation-projected Gaussian and Lorentzian linewidths.
// Used by: Discrepancy.cpp, IntegratedSpectrum.cpp and VoitlanderTriangleIntegrals.cpp; both functions are declared in
// BiradicalODF.h. The results are turned into spectra by IndividualSpectrum.cpp.

#include "BiradicalODF.h"
#include "unsupported/Eigen/KroneckerProduct"

extern int MinimizationFlag;

// The function generates the array of transitions (pair Field-Intensity)
// for rigid biradical with magnetic parameters given in the object MagneticParameter
// with orientation defined by rotation matrix Field_to_Molecule_RotMat.
// Also the values of effective Gaussian and Lorentzian linewidths corresponding to this orientation are calculated.
// Spin system: two electron spins 1/2 (4x4 Hamiltonian: Zeeman terms of both electrons, spin-spin tensor DJ, and
// effective hyperfine terms Aeff*mz*Sz with mz=-1,0,1 for each electron). The 9 combinations of (mz1,mz2) times
// the 6 pairs of the 4 eigenstates give the 54 rows of the output: column 0 is the transition field (Gauss),
// column 1 the intensity |<k|Hp|j>|^2 of the microwave perturbation Hp. HGaussEff and HLorentzEff (outputs) are
// full widths at half maximum in Gauss.
// Used by: Discrepancy() in Discrepancy.cpp; CalculateIntegratedSpectrumLiouville() in IntegratedSpectrum.cpp;
// VoitlanderVertexEvalRigid() in VoitlanderTriangleIntegrals.cpp.
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
	// the list of transitions is logged to Transitions.out only for a single-orientation run without fitting
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
	// spin-spin tensor DJ (dipolar plus exchange) is taken ready from MagneticParameter, where RotateTensors() builds it
	DJ=MagneticParameter->DJ;

// Initializing Array for Transition intensities
	Array<double, 54, 2> Transitions=Array<double, 54, 2>::Zero();
	int r=0;
// Initializing Magnetic Field related values;
	Matrix<double, 1, 3> B0, B1, nv;

	// static field B0 along z of the resonator frame, scaled by the rms g value of the first particle, and microwave field B1 along y
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

	// linewidth tensors projected on the field direction through the g tensor and divided by the rms g value, as FWHM in Gauss
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

    // electron Zeeman terms of both spins in the resonator frame (field B0 contracted with the g tensors)
    H11=    (B0(0)*g1M(0, 0)+B0(1)*g1M(1, 0)+B0(2)*g1M(2, 0))*kroneckerProduct(Sx, Si);
    H11=H11+(B0(0)*g1M(0, 1)+B0(1)*g1M(1, 1)+B0(2)*g1M(2, 1))*kroneckerProduct(Sy, Si);
    H11=H11+(B0(0)*g1M(0, 2)+B0(1)*g1M(1, 2)+B0(2)*g1M(2, 2))*kroneckerProduct(Sz, Si);

    H12=    (B0(0)*g2M(0, 0)+B0(1)*g2M(1, 0)+B0(2)*g2M(2, 0))*kroneckerProduct(Si, Sx);
    H12=H12+(B0(0)*g2M(0, 1)+B0(1)*g2M(1, 1)+B0(2)*g2M(2, 1))*kroneckerProduct(Si, Sy);
    H12=H12+(B0(0)*g2M(0, 2)+B0(1)*g2M(1, 2)+B0(2)*g2M(2, 2))*kroneckerProduct(Si, Sz);

    // spin-spin coupling term: S1.DJ.S2
    H3=   DJ(0, 0)*kroneckerProduct(Sx, Sx);
    H3=H3+DJ(0, 1)*kroneckerProduct(Sx, Sy);
    H3=H3+DJ(0, 2)*kroneckerProduct(Sx, Sz);

    H3=H3+DJ(1, 0)*kroneckerProduct(Sy, Sx);
    H3=H3+DJ(1, 1)*kroneckerProduct(Sy, Sy);
    H3=H3+DJ(1, 2)*kroneckerProduct(Sy, Sz);

    H3=H3+DJ(2, 0)*kroneckerProduct(Sz, Sx);
    H3=H3+DJ(2, 1)*kroneckerProduct(Sz, Sy);
    H3=H3+DJ(2, 2)*kroneckerProduct(Sz, Sz);

    // microwave perturbation operator Hp: Zeeman form with B1 in place of B0, for both electrons
    Hp= (B1(0)*g1M(0, 0)+B1(1)*g1M(1, 0)+B1(2)*g1M(2, 0))*kroneckerProduct(Sx, Si)+
		(B1(0)*g1M(0, 1)+B1(1)*g1M(1, 1)+B1(2)*g1M(2, 1))*kroneckerProduct(Sy, Si)+
		(B1(0)*g1M(0, 2)+B1(1)*g1M(1, 2)+B1(2)*g1M(2, 2))*kroneckerProduct(Sz, Si)+
		(B1(0)*g2M(0, 0)+B1(1)*g2M(1, 0)+B1(2)*g2M(2, 0))*kroneckerProduct(Si, Sx)+
		(B1(0)*g2M(0, 1)+B1(1)*g2M(1, 1)+B1(2)*g2M(2, 1))*kroneckerProduct(Si, Sy)+
		(B1(0)*g2M(0, 2)+B1(1)*g2M(1, 2)+B1(2)*g2M(2, 2))*kroneckerProduct(Si, Sz);

	// header line of the optional transition log
	if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) 
		fprintf(transitions, "   mz1  mz2  ms1i  ms1f  ms2i  ms2f    Field    Intensity \n");

	// loop over the nuclear projections mz1 and mz2 of the two hyperfine couplings, one 4x4 diagonalisation each
	for(mz1=-1; mz1<=1; mz1++)
		for(mz2=-1; mz2<=1; mz2++)
		{
 			H21=(Aeff1*mz1)*kroneckerProduct(Sz, Si);
			H22=(Aeff2*mz2)*kroneckerProduct(Si, Sz);
			SpinHamiltonian=H11+H12+H21+H22+H3;

// Computing eigenvalues:
			SelfAdjointEigenSolver<Matrix<std::complex<double>, 4, 4> > es;
			es.compute(SpinHamiltonian);
			// microwave perturbation in the eigenbasis of the spin Hamiltonian, its off-diagonal elements give the intensities
			Intensity=es.eigenvectors().adjoint() * Hp * es.eigenvectors();
    
			// all 6 pairs of the 4 eigenstates: transition field and intensity, ms labels (from the level index) are for the log only
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
		// close the optional transition log and return the 54 transitions
		if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) fclose(transitions);
  return Transitions;

}



//Calculation of transitions for Monoradical
// Same as CalculateTransitions() for a rigid monoradical: one electron spin 1/2 (2x2 Hamiltonian: Zeeman term and
// effective hyperfine term Aeff*mz*Sz), 3 values mz=-1,0,1 and one transition each, so the output has 3 rows:
// column 0 the transition field (Gauss), column 1 the intensity, scaled by a constant factor 12. HGaussEff and
// HLorentzEff (outputs) are full widths at half maximum in Gauss. The mz2 and ms2 columns of the log are not
// meaningful here.
// Used by: Discrepancy() in Discrepancy.cpp; CalculateIntegratedSpectrumLiouville() in IntegratedSpectrum.cpp;
// VoitlanderVertexEvalRigid() in VoitlanderTriangleIntegrals.cpp.
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
	// the list of transitions is logged to Transitions.out only for a single-orientation run without fitting
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
	// spin-spin tensor DJ is read for symmetry with the biradical routine but not used by the monoradical Hamiltonian below
	DJ=MagneticParameter->DJ;

// Initializing Array for Transition intensities
	Array<double, 3, 2> Transitions=Array<double, 3, 2>::Zero();
	int r=0;
// Initializing Magnetic Field related values;
	Matrix<double, 1, 3> B0, B1, nv;

	// static field B0 along z of the resonator frame, scaled by the rms g value of the first particle, and microwave field B1 along y
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

	// linewidth tensors projected on the field direction through the g tensor and divided by the rms g value, as FWHM in Gauss
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

    // electron Zeeman term in the resonator frame (field B0 contracted with the g tensor)
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
    // microwave perturbation operator Hp: Zeeman form with B1 in place of B0
    Hp= (B1(0)*g1M(0, 0)+B1(1)*g1M(1, 0)+B1(2)*g1M(2, 0))*Sx+
		(B1(0)*g1M(0, 1)+B1(1)*g1M(1, 1)+B1(2)*g1M(2, 1))*Sy+
		(B1(0)*g1M(0, 2)+B1(1)*g1M(1, 2)+B1(2)*g1M(2, 2))*Sz;

	// header line of the optional transition log
	if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) 
		fprintf(transitions, "   mz1  mz2  ms1i  ms1f  ms2i  ms2f    Field    Intensity \n");

	// loop over the nuclear projection mz1 of the hyperfine coupling, one 2x2 diagonalisation each
	for(mz1=-1; mz1<=1; mz1++)
//		for(mz2=-1; mz2<=1; mz2++)
		{
 			H21=(Aeff1*mz1)*Sz;
//			H22=(Aeff2*mz2)*kroneckerProduct(Si, Sz);
			SpinHamiltonian=H11+H21;

// Computing eigenvalues:
			SelfAdjointEigenSolver<Matrix<std::complex<double>, 2, 2> > es;
			es.compute(SpinHamiltonian);
			// microwave perturbation in the eigenbasis of the spin Hamiltonian, its off-diagonal element gives the intensity
			Intensity=es.eigenvectors().adjoint() * Hp * es.eigenvectors();
    
			// the single pair of the 2 eigenstates: transition field and intensity, ms labels (from the level index) are for the log only
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
		// close the optional transition log and return the 3 transitions
		if((MinimizationFlag==0)&&(MagneticParameter->IntegrationNet==1)) fclose(transitions);
  return Transitions;

}

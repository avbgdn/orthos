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
#include "Eigen/Sparse"
// #include "unsupported/Eigen/kroneckerProduct"
#include "SphericalTensor.h"
#include "SpinMatrixB.h"
#include "lanczos.h"
#include "BasisB.h"
#include "LiouvilleMatrixElementB.h"
#include "fftw3.h"


// extern clock_t time_spin, time_diffusion, time_potential;

#define MAX_BASIS 200000
#define VERBOSE 1
//extern double B0;

bool FduReal;
bool FduEven;
bool FduAxial;
bool DirectorParallel;
bool IsGroupMonoradical;
std::fstream timing;

extern int MpiRank, MpiSize;
int ParticleCount;

// Важные замечания: 
// 1. оператор AgL работает с величиной B0, изменённой в g_ср раз, чтобы корректно предсказывать центр поля.
// 3. При расчёте матричных элементов (без симметрии) беру комплексное сопряжение от всех Fmu


// Pruning routine.
// The Stochastic Liouville equation is solved by conjugate gradient method in [nprunes] points of spectra separated by [dz] Gauss.
// Only basis vectors, which contribute to the resulting spectrum, are selected and written to "pruneN.pru" file, where N is the number of spectrum.
void Pruning(SparseMatrix<std::complex<double> > &Z, const Matrix<std::complex<double>, Dynamic, 1> &b, int ndim, BasisVectorB *basis, double z0, int nprune, double dz, double PruningTolerance)
{
	std::fstream pruningfile;
        char buff[MAX_STRLEN];
        sprintf(buff, "basis%d.pru", MpiRank);
	pruningfile.open(buff, std::fstream::out);
	std::complex<double> I(0,1);
//	int nprune=1;
        if(!MpiRank)
	timing << "Started pruning routine\n\n";
	Matrix<std::complex<double>, Dynamic, 1> x;
	x.resize(ndim, 1);

	Matrix<double, Dynamic, Dynamic> res;
	res.resize(ndim,nprune);

	SparseMatrix<std::complex<double> > U(ndim,ndim);
	SparseMatrix<std::complex<double> > ZZ(ndim,ndim);
	U.setIdentity();
	BiCGSTAB<SparseMatrix<std::complex<double> > > solver;

//	for(double z=-52.0; z<=52.0; z+=104.0/(nprune-1))
	for(int i=0; i<nprune; i++)
	{
		double z;
		if(nprune) z=(double)(2*i-nprune+1)*dz/(double)(2*nprune-2);
		else	   z=0;
//		double z=0;
		ZZ=Z+U*(z0+z)*I;
		solver.compute(ZZ);
		x=solver.solve(b);
		res.col(i)=x.cwiseAbs();
		timing << i << ", rel.field=" << z << ", #iterations:     " << solver.iterations() << ", estimated error: " << solver.error() << std::endl;
	}
	pruningfile << "Pruning result:\n";
	pruningfile << res << endl << endl;
	pruningfile << "Vectors:\n";
	int npruned=0;
	for(int i=0; i<ndim; i++)
	{
		if (res.row(i).maxCoeff()>=PruningTolerance)
		{
			char buffer[MAX_STRLEN];
			sprintf(buffer, " % d % d % d % d % d  % d  % d  % d  % d  % d  % d  % d  % d\n", basis[i].L, basis[i].K, basis[i].M, basis[i].N, basis[i].Q, basis[i].pS1, basis[i].qS1, basis[i].pS2, basis[i].qS2, basis[i].pI1, basis[i].qI1 , basis[i].pI2, basis[i].qI2 );
			pruningfile << buffer;
			npruned++;
		}


	}
	pruningfile << "Dimension of minimized basis is " << npruned << endl;

	pruningfile.close();
}

//
// The routine of spectrum modelling according Stochastic Liouville equatinon.
//
Matrix<double, Dynamic, 1> CalculateLiouvillianB(MagneticParameters *MagneticParameter, double psi, char* BasisSource )
{
	if(MagneticParameter->NConformers<1)
	{
			int npoints=MagneticParameter->NPoints;
			Matrix<double, Dynamic, 1> Spectrum;
			Spectrum.resize(npoints, 1);
			for(int i=0; i<npoints; i++)
				Spectrum(i,0)=0.0;
			return Spectrum;
	}
//	Matrix<std::complex<double>, 3, 3> g1, A1, g2, A2, W, DJ, R;
//	double Diffusion_Type[3], R_Tau[3];
	clock_t time;

	time=clock();
	MagneticParameters *CurrentParameter;
	SpinInteractionB *CurrentInteraction;


/*
	time_spin=(clock_t)0;
	time_diffusion=(clock_t)0;
	time_potential=(clock_t)0;
	clock_t time_inside=(clock_t)0;
	clock_t time_inside2=(clock_t)0;
*/

	std::complex<double> I(0, 1);

	ParticleCount=MagneticParameter->Count;
// Output files are prepared
	char buff[MAX_STRLEN];
	std::fstream lioufile;
        sprintf(buff, "liou%d_%d.out", MagneticParameter->Count, MpiRank);
	lioufile.open(buff, std::fstream::out);

	std::fstream tensorfile;
        sprintf(buff, "tensors%d_%d.out", MagneticParameter->Count, MpiRank);
	tensorfile.open(buff, std::fstream::out);

	std::fstream tridiagonalfile;
        sprintf(buff, "tridiag%d_%d.out", MagneticParameter->Count, MpiRank);
	tridiagonalfile.open(buff, std::fstream::out);

        if(!MpiRank)
	timing.open("timing.out", std::fstream::app | std::fstream::out);

        if(!MpiRank)
	timing << "Entered Liouvillian calculation routine\n";

// Magnetic and rotation parameters are prepared - all of this moved to MagneticParameter->RotateTensors() routine

/*
	SpinInteractionB Interactions[6];
	for(int j=0; j<3; j++)
	{
		for(int k=0; k<3; k++)
		{
			g1(j, k)=MagneticParameter->g1R(j, k);
			g2(j, k)=MagneticParameter->g2R(j, k);
			A1(j, k)=MagneticParameter->A1R(j, k);
			A2(j, k)=MagneticParameter->A2R(j, k);
			W(j, k) =MagneticParameter->HLorentzR(j, k);
			R(j, k) =MagneticParameter->DiffusionR(j, k);
			DJ(j, k)=MagneticParameter->DJ(j, k);
		}
		Diffusion_Type[j]=MagneticParameter->DiffusionType[j];
		R_Tau[j]=MagneticParameter->RTau[j];
	}
*/

//	double JumpFrequency=MagneticParameter->JumpFrequencyR;

//	B0=MagneticParameter->CenterField*3/(MagneticParameter->g1R.trace());
	
//	double w0=real(MagneticParameter->HLorentzR.trace())/3.0;   
//	double w0=real(W.trace())/3.0;   

/*
// Spherical tensors are constructed
	SphericalTensor *Fmg1, *Fmg2, *FmA1, *FmA2, *FmD, *FmW, *FmR;
	Fmg1=new SphericalTensor(g1);
	FmA1=new SphericalTensor(A1);
	Fmg2=new SphericalTensor(g2);
	FmA2=new SphericalTensor(A2);
	FmD=new SphericalTensor(DJ);
	FmW=new SphericalTensor(W);
	FmR=new SphericalTensor(R);

// Object, corresponding to each interaction in the spherical representation, is constructed for each term in the Spin Hamiltonian
	Interactions[0].Fdu=*Fmg1;
	Interactions[0].AuL.Set(Ag1L);
	strcpy(Interactions[0].Fdu.Name, "g1");

	Interactions[1].Fdu=*Fmg2;
	Interactions[1].AuL.Set(Ag2L);
	strcpy(Interactions[1].Fdu.Name, "g2");

	Interactions[2].Fdu=*FmA1;
	Interactions[2].AuL.Set(AA1L);
	strcpy(Interactions[2].Fdu.Name, "A1");

	Interactions[3].Fdu=*FmA2;
	Interactions[3].AuL.Set(AA2L);
	strcpy(Interactions[3].Fdu.Name, "A2");

	Interactions[4].Fdu=*FmD;
	Interactions[4].AuL.Set(ADL);
	strcpy(Interactions[4].Fdu.Name, "DJ");

	Interactions[5].Fdu=*FmW;
	Interactions[5].AuL.Set(AWLB);
	strcpy(Interactions[5].Fdu.Name, "Lorentz");
*/

// Analyzing spherical tensors for symmetry
//	FduReal, FduEven, FduAxial;

	if(psi==0)	DirectorParallel=true;
	else		DirectorParallel=false;
	FduReal=true;
	FduEven=true;
	FduAxial=true;
	IsGroupMonoradical=true;
	for(int q=0; q<MagneticParameter->NConformers; q++)  // Conformation analysis:
	{
//		printf("q=%d, groupmember=%d\n", q, MagneticParameter->ConformerGroupMembers[q]);
		CurrentParameter=&(MagneticParameter->MagneticParameterArray[MagneticParameter->ConformerGroupMembers[q]]);
		CurrentInteraction=(CurrentParameter->Interactions);

		if (!(CurrentParameter->isMonoradical)) IsGroupMonoradical=false;

		for(int i=0; i<6; i++)
		{
			CurrentInteraction[i].Fdu.PrintOut(tensorfile);
			for(int m=-2; m<=2; m++)
			{
				if((m%2)&& (abs( CurrentInteraction[i].Fdu.Get(2, m)) >= 1.0e-10 ))	FduEven =false;
				if((m  )&& (abs( CurrentInteraction[i].Fdu.Get(2, m)) >= 1.0e-10 ))	FduAxial=false;
				if((imag( CurrentInteraction[i].Fdu.Get(2, m)) >= 1.0e-10 ))			FduReal =false;
			}
		}
//        cout << FduReal << FduEven << FduAxial <<endl;
// Checking whether FmR is even, and/or axial, and/or real
		for(int m=-2; m<=2; m++)
		{
			if((m%2)&& (abs( CurrentParameter->FmR->Get(2, m)) >= 1.0e-10 ))	FduEven =false;
			if((m  )&& (abs( CurrentParameter->FmR->Get(2, m)) >= 1.0e-10 ))	FduAxial=false;
			if((imag( CurrentParameter->FmR->Get(2, m)) >= 1.0e-10 ))			FduReal =false;
		}
	}

// Basis generation / reading:	
	BasisVectorB* Basis, *Basis1;
        if(!MpiRank)
	timing << "Started basis generation!\n";
        if(!MpiRank)
        timing << "Calculation time: " << (float)(clock()-time)/CLOCKS_PER_SEC << " sec." << std::endl;
	Basis=(BasisVectorB*)malloc(MAX_BASIS*sizeof(BasisVectorB));
	Basis1=(BasisVectorB*)malloc(MAX_BASIS*sizeof(BasisVectorB));

	int NBasis;

// If the filename of the file containing basis vector is specified in name.nam file, then read basis from file
	if (strcmp(BasisSource, "truncation") && strcmp(BasisSource, ""))		
												NBasis=LoadBasisB(BasisSource, Basis);
// If the filename is not specified, generate basis.
	else
	{
		if(!IsGroupMonoradical)					NBasis=GenerateBasisB(MagneticParameter->LEmax, 
																	  MagneticParameter->LOmax, 
																	  MagneticParameter->Kmax, 
																	  MagneticParameter->Mmax, 
																	  MagneticParameter->Nmax, 2,
																	  MagneticParameter->NConformers, Basis, Basis1);
		else									NBasis=GenerateBasisM(MagneticParameter->LEmax, 
																	  MagneticParameter->LOmax, 
																	  MagneticParameter->Kmax, 
																	  MagneticParameter->Mmax, 
																	  MagneticParameter->Nmax, 2,
																	  MagneticParameter->NConformers, Basis, Basis1);
	}
//	Matrix<int, Dynamic, Dynamic> Conjugates;
//	Conjugates.resize(NBasis, NBasis);
//	Conjugates=Matrix<int, Dynamic, Dynamic>::Zero(NBasis, NBasis);

        if(!MpiRank)
	timing << "Starting Matrix Element Calculation.\n";
// Precalculate wigner d-functions for the current orientation:
	precalculate_wigner(0, psi, 0);

        if(!MpiRank)
        timing << "Calculation time: " << (float)(clock()-time)/CLOCKS_PER_SEC << " sec." << std::endl;



	std::vector< Eigen::Triplet<std::complex<double> > > tripletList;

//	tripletList.reserve(NBasis*NBasis/5);

	double LiouvilleTruncation=1.0e-5;

	for(int j=0; j<NBasis; j++)
		for(int i=j; i<NBasis; i++)
		{
			std::complex<double> liou, liou1, liou2;
			std::complex<double> pot1, pot2;
            if(Basis[j].Q!=Basis[i].Q)
			{
				liou=MatrixElementConformation(Basis[j], Basis[i], MagneticParameter->ConversionRateR);
			}
			else
			{
//			clock_t time_start=clock();
				CurrentParameter=&(MagneticParameter->MagneticParameterArray[MagneticParameter->ConformerGroupMembers[Basis[j].Q]]);

				if(!(CurrentParameter->isMonoradical))
					liou=MatrixElementIsotropicBrownianB_ND(Basis[j], Basis[i], CurrentParameter->DiffusionR, CurrentParameter->DiffusionType, CurrentParameter->RTau)			-I*MatrixElementBSK1(Basis[j], Basis[i], CurrentParameter->Interactions, 6, psi, CurrentParameter->J2);
				else
					liou=MatrixElementIsotropicBrownianBSM_ND(Basis[j], Basis[i], CurrentParameter->DiffusionR, CurrentParameter->DiffusionType, CurrentParameter->RTau)		-I*MatrixElementBSM(Basis[j], Basis[i], CurrentParameter->Interactions, 6, psi, CurrentParameter->J2);

				liou+=MatrixElementConformation(Basis[j], Basis[i], MagneticParameter->ConversionRateR);
				liou+=MatrixElementJDiffusion(Basis[j], Basis[i], MagneticParameter->JDiffusionR);


// adding term of diffusion in anisotropic potential:
				if(MagneticParameter->MFPotential.Lmax)
				  if(!(CurrentParameter->isMonoradical))
					liou+=MatrixElementPotentialB_ND(Basis[j], Basis[i], CurrentParameter->DiffusionR, CurrentParameter->MFPotential);
				  else
					liou+=MatrixElementPotentialBSM_ND(Basis[j], Basis[i], CurrentParameter->DiffusionR, CurrentParameter->MFPotential);

// adding term of Jumps:
				if(CurrentParameter->NJumpAxes)
				{
					if(!(CurrentParameter->isMonoradical))
					{
						// biradical (!Not working before DEBUG of Jumps!
					}
					else
					{
					// monoradical	

						if(CurrentParameter->MFPotential.Lmax)
							liou+=MatrixElementJumpsBSM_Pot(Basis[j], Basis[i], CurrentParameter, CurrentParameter->MFPotential);
						else
							liou+=MatrixElementJumpsBSM_Iso(Basis[j], Basis[i], CurrentParameter);

/*
					if(MagneticParameter->MFPotential.Lmax)
						liou1=MatrixElementJumpsBSM_Pot(Basis[j], Basis[i], MagneticParameter, MagneticParameter->MFPotential);
					else
						liou1=MatrixElementJumpsBSM_Iso(Basis[j], Basis[i], MagneticParameter);

					if(MagneticParameter->MFPotential.Lmax)
						liou2=MatrixElementJumpsBSM_Pot(Basis[i], Basis[j], MagneticParameter, MagneticParameter->MFPotential);
					else
						liou2=MatrixElementJumpsBSM_Iso(Basis[i], Basis[j], MagneticParameter);
*/
					}
				}

			}

//			time_inside+=(clock()-time_start);


			if(abs(liou)>=LiouvilleTruncation)
			{
				tripletList.push_back(Eigen::Triplet<std::complex<double> >(j,i, liou));
				if(i!=j) 
					tripletList.push_back(Eigen::Triplet<std::complex<double> >(i,j, liou));
//				Liou.insert(j, i)=liou;
//				if(i!=j) Liou.insert(i, j)=liou;			

				lioufile << "<"<<j+1<<"|L|"<<i+1<<"> = "<< liou << endl;
			}
//				lioufile << "<"<<j+1<<"|L|"<<i+1<<"> = "<< liou << "    " << liou1 << "    " <<liou2 << endl;
//			if( abs(liou1)>=LiouvilleTruncation ) lioufile << "<"<<j+1<<"|L|"<<i+1<<"> = " << liou1 << " " << liou2 << endl;

//			if( real(liou1-liou)>=LiouvilleTruncation ) lioufile << "!!!!!!!!DIFFUS!!!!!!!<"<<j+1<<"|L|"<<i+1<<"> = " << liou << " " << liou1 << " " << liou2 << endl;
//			if( abs(liou1-liou2)>=LiouvilleTruncation ) lioufile << "!!!!!!!!!!!!!!!!!!!!!<"<<j+1<<"|L|"<<i+1<<"> = " << liou << " " << liou1 << " " << liou2 << endl;
		}
//		clock_t time_start=clock();

		Eigen::SparseMatrix<std::complex<double> > Liou(NBasis, NBasis);

		Liou.setFromTriplets(tripletList.begin(), tripletList.end());


		tripletList.clear();
		tripletList=std::vector< Eigen::Triplet<std::complex<double> > >();
//		tripletList.shrink_to_fit();
//		time_inside2+=(clock()-time_start);


        if(!MpiRank)
	timing << "Calculated Liouvillian matrix. \n";
        if(!MpiRank)
	timing << "Calculation time: " << (float)(clock()-time)/CLOCKS_PER_SEC << " sec." << std::endl;

	Matrix<std::complex<double>, Dynamic, 1> StartingVector;
    Matrix<std::complex<double>, Dynamic, 2> ab;
	StartingVector.resize(NBasis, 1);
	ab.resize(NBasis, 2);

// constructing starting vector
    Matrix<std::complex<double>, 1296, 1> sv;
	sv=StartingVectorSpinB();
	for(int i=0; i<NBasis; i++)
	{
		CurrentParameter=&(MagneticParameter->MagneticParameterArray[MagneticParameter->ConformerGroupMembers[Basis[i].Q]]);

		std::complex<double> OP=(CurrentParameter->MFPotential.GetSOP(Basis[i].L, abs(Basis[i].K), Basis[i].M));
		if((Basis[i].L+Basis[i].K)%2) OP-=(CurrentParameter->MFPotential.GetSOP(Basis[i].L, -abs(Basis[i].K), Basis[i].M));
		else 						  OP+=(CurrentParameter->MFPotential.GetSOP(Basis[i].L, -abs(Basis[i].K), Basis[i].M));

		if(Basis[i].K==0) OP/=sqrt(2.0);

		StartingVector(i, 0)=sqrt(2.0*Basis[i].L+1)
							*OP
							*sqrt(fabs(CurrentParameter->AmplitudeFactor))
							*sv(SpinIndexB(Basis[i]));

		if(Basis[i].N!=0) StartingVector(i,0)=0;
//		cout << OP << "\t" << sqrt(fabs(CurrentParameter->AmplitudeFactor)) << "\t" << sv(SpinIndexB(Basis[i])) << "\t" << SpinIndexB(Basis[i]) << "\n";

	}
// Normalization of Starting Vector
	StartingVector=StartingVector/sqrt((StartingVector.adjoint()*StartingVector)(0,0));

	lioufile << "Starting Vector: \n" << StartingVector << endl << endl;

        if(!MpiRank)
	timing << "Calculated Starting vector. \n";
        if(!MpiRank)
	timing << "Calculation time: " << (float)(clock()-time)/CLOCKS_PER_SEC << " sec." << std::endl;

	if(MagneticParameter->NPrune)	Pruning(Liou, 
											StartingVector, 
											NBasis, 
											Basis, 
											MagneticParameter->CenterField, 
											MagneticParameter->NPrune,
											MagneticParameter->dHPrune,
											MagneticParameter->PruningTolerance);
        if(!MpiRank)
	timing << "Starting Lanczos routine\n";

	// The Lanczos truncation depth cannot exceed the dimension of the basis
	int LanczosSteps=MagneticParameter->LanczosSteps;
	if(LanczosSteps>NBasis)
	{
		printf("LanczosSteps=%d exceeds the basis dimension %d for Particle #%d, the Lanczos depth is reduced to the basis dimension\n", MagneticParameter->LanczosSteps, NBasis, MagneticParameter->Count);
		LanczosSteps=NBasis;
	}

	int ActualSteps;
	ab=Lanczos(Liou, StartingVector, LanczosSteps, NBasis, MagneticParameter, ActualSteps);
	if(ActualSteps<LanczosSteps)
		tridiagonalfile << "Krylov subspace closed early at step " << ActualSteps << " of " << LanczosSteps
						<< " requested: the residual vanished to machine precision (exact invariant" << endl
						<< "subspace/eigenvector found), so the recursion was stopped instead of dividing by a near-zero pivot." << endl << endl;
	tridiagonalfile << "a, b vector: \n" << ab << endl << endl;

	int npoints=MagneticParameter->NPoints;
	Matrix<double, Dynamic, 2> Spectrum;
	Spectrum.resize(npoints, 2);
	Spectrum=CFVD(ab, LanczosSteps, MagneticParameter->CenterField, npoints, 0, MagneticParameter);

//	FILE *SpecOut;
    fstream SpecOut;
	SpecOut.open("spec", std::fstream::out);
	SpecOut << Spectrum;
	SpecOut.close();

        if(!MpiRank)
	timing << "End Liouvillian calculation\n";
        if(!MpiRank)
        timing << "Calculation time: " << (float)(clock()-time)/CLOCKS_PER_SEC << " sec." << std::endl;

// time for Fourier transform
	double LWGauss;
// LWGauss=A+(B-A)*sin^2(psi)
//     where A=HG_zz, B=HG_xx
// so, at psi=0    LWGauss=HG_zz
//     at psi=pi/2 LWGauss=HG_xx
	LWGauss=fabs(MagneticParameter->HGauss(2,2)+sin(psi)*sin(psi)*(MagneticParameter->HGauss(0,0)-MagneticParameter->HGauss(2,2)));
	if(LWGauss)
	{
		double *fftwin;
		fftw_plan fftwp;

		fftwin = (double*) fftw_malloc(sizeof(fftw_complex) * npoints);

		fftwp=fftw_plan_r2r_1d(npoints, fftwin, fftwin,
						 FFTW_REDFT00, FFTW_ESTIMATE);
		for(int i=0; i<npoints; i++)
			fftwin[i]=Spectrum(i, 1);
		fftw_execute(fftwp); 
// A DCT-I of npoints samples spaced by fs is the transform of an even extension of period
//    2*(npoints-1) samples, so frequency index i is i/(2*(npoints-1)*fs) cycles per Gauss, and the
//    transform of a unit-area Gaussian of full width at half maximum LWGauss is
//    exp(-(PI*LWGauss*i/(4*sqrt(ln2)*(npoints-1)*fs))^2). The leading factor undoes the
//    2*(npoints-1) scaling that REDFT00 introduces when it is applied twice.
		for(int i=0; i<npoints; i++)
			fftwin[i]*=0.5*exp(-pow(LWGauss*PI*(double)i/(4*sqrt(log(2.0))*(npoints-1)*MagneticParameter->fs), 2))/(double)(npoints-1);

		fftw_execute(fftwp); 

		for(int i=0; i<npoints; i++)
			Spectrum(i, 1)=fftwin[i];

		fftw_destroy_plan(fftwp);
		fftw_free(fftwin); 
	}

/*
	cout << "Spin time: " << (float)time_spin/CLOCKS_PER_SEC << " sec." << endl;
	cout << "Diff time: " << (float)time_diffusion/CLOCKS_PER_SEC << " sec." << endl;
	cout << "Potl time: " << (float)time_potential/CLOCKS_PER_SEC << " sec." << endl;
	cout << "Insd time: " << (float)time_inside/CLOCKS_PER_SEC << " sec." << endl;
	cout << "Ins2 time: " << (float)time_inside2/CLOCKS_PER_SEC << " sec." << endl;
*/

	free(Basis);
	free(Basis1);
	lioufile.close();
	tensorfile.close();
	tridiagonalfile.close();
	timing.close();

//	if(psi&&(!(MagneticParameter->isMonoradical)))  return Spectrum.col(1)*2*MagneticParameter->AmplitudeFactor/PI;
	return Spectrum.col(1)*MagneticParameter->AmplitudeFactorAllConformers/PI;

}
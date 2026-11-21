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

#include <iostream>
#include <complex>
#include <cmath>
#include <math.h>
#include "BiradicalODF.h"
#include "SphericalTensor.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "Eigen/Eigenvalues"
#include <stdlib.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include "nlopt.h"



#ifdef __MPI
#include <mpi.h>
#endif

//The library, which contains minimizer NL2SOL
// #pragma comment(lib, "FortranLibrary64.lib")

int MpiRank, MpiSize;
int MinimizationFlag;
DataForCalculation *p_additionaldata;

// The file for writing time consumption of the program
extern std::fstream timing;

double cov[MAX_ARG][MAX_ARG];		// covariance matrix

// the objective function for minimization with NLOpt library optimization algorithms
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data);

// the objective function for minimization with NL2SOL algorithm
void *madr(int* N, int* P, double* X, int* NF, double* R, int* UIPARM, double* URPARM, void* UFPARM, int* IV);

// DECLARATION OF FUNCTIONS WITHIN NL2SOL library:
#ifdef __linux
  extern "C" void nl2sno_(int *n, int *p, double *x, void *CALCR(int*, int*, double*, int*, double*, int*, double*, void*, int*), int *iv, double *v, int *uiparm, double *urparm, void* ufparm);
  extern "C" void dfault_(int *iv, double *v);
#else
  extern "C" void NL2SNO(int *n, int *p, double *x, void *CALCR(int*, int*, double*, int*, double*, int*, double*, void*, int*), int *iv, double *v, int *uiparm, double *urparm, void* ufparm);
  extern "C" void DFAULT(int *iv, double *v);
#endif

 // COMMON BLOCKS FOR NL2SOL:
#ifdef __linux
  extern "C"
  {  
	extern struct
	{
		double STEPJAC;
		double STEPBOUND;
	} stepj_;	
  }

  extern "C"
  {
	extern struct
	{
		double PARFCTOL;
		double PARXCTOL;
	} converg_;	
  }

  extern "C"
  {
	extern struct
	{
		int keyErr;
		double ANoiseLev;
	} keyerr_;	
  }

  extern "C"
  {
	extern struct
	{
		int KCONTROL;
	} kctrl_;
  } 
#else
  extern "C"
  {
	extern struct
	{
		double STEPJAC;
		double STEPBOUND;
	} STEPJ;	
  }

  extern "C"
  {
	extern struct
	{
		double PARFCTOL;
		double PARXCTOL;
	} CONVERG;	
  }

  extern "C"
  {
	extern struct
	{
		int keyErr;
		double ANoiseLev;
	} KEYERR;	
  }

  extern "C"
  {
	extern struct
	{
		int KCONTROL;
	} KCTRL;
  } 

#endif
//     common / STEPJ  / stepJac
//     &,                 stepBound
//	 DOUBLE PRECISION PARFCTOL, PARXCTOL
//	 COMMON /CONVERG/ PARFCTOL, PARXCTOL !ax!
//	common/KeyErr/keyErr,ANoiseLev		 !ax!

//   INTEGER KCONTROL
//   COMMON /KCTRL/ KCONTROL     !AB!


// #include <direct.h>

extern double ***SpectrumSum;
double ***SpectrumSum;

extern void CalculateLiouvillian();

// The main function of the program
int main(int argc, char** argv)
{
	
#ifdef __MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &MpiSize);
#else
    MpiRank=0;
    MpiSize=1;
#endif
// Print program version:
	Version();
// filenames:
	char magnetic_in_filename[MAX_STRLEN], ODF_in_filename[MAX_STRLEN];
	char magnetic_out_filename[MAX_STRLEN], ODF_out_filename[MAX_STRLEN];

	int NParticles;
	double *R;
	bool CovarianceCalculated=false;
//	char buffer[MAX_STRLEN];
//	getcwd(buffer, MAX_STRLEN);
//	cout << buffer << endl;
//
	clock_t ct_time;
	ct_time=clock();

	time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    if(!MpiRank)
    {
	printf ( "Start program. Time: %s\n", asctime (timeinfo) );


// Prepare file for timing;
	timing.open("timing.out", std::fstream::out);
	timing.close();

// Prepare file for discrepancy;
	std::fstream discr;
	discr.open("discrepancy", ios::out);
	discr.close();

// Prepare file for transitions;
	FILE *transitions;
	transitions=fopen("Transitions.out", "w");
	fclose(transitions);
    }

// Check Number of Particles:
	ifstream namenam;
	namenam.open ("name.nam", ios::in);
	if(!namenam.is_open())
	{
		std::cout << "Could not open file name.nam\n";
		exit(0);
	}
	while (!namenam.eof())
	{
		char buffer[MAX_STRLEN], buffer2[MAX_STRLEN], buffer3[MAX_STRLEN];
		namenam.getline(buffer, MAX_STRLEN);
		if( sscanf(buffer, "%s", buffer2)==-1) continue;
		if (!strcmp(buffer2, "NParticles:"))
		{
			sscanf(buffer, "%s %d", buffer3, &NParticles);
			break;
		}
	}
	namenam.close();

// Each object MagneticParameter[i] contains the set of Magnetic (and some other parameters) of i-th paramagnetic particle
// note that biradical molecule is treated as a single paramagnetic particle,
// different particles are moieties wich do not interact magnetically.
	MagneticParameters *MagneticParameter=new MagneticParameters[NParticles];

// Each object ODF[i] contains the set of Order parameters of i-th paramagnetic particle
	OrderParameters *ODF=new OrderParameters[NParticles];

// Read file "name.nam" to find out filenames of parameter and spectra files:
	ReadName(&(MagneticParameter[0]), &(ODF[0]), magnetic_in_filename, ODF_in_filename, magnetic_out_filename, ODF_out_filename);

// For the simulation of angular dependencies of EPR spectra, the set of angles (hi_SF) an rotation axes indications (axis_SF), which
//     define the orientation of the sample in the magnetic field of the EPR spectrometer.
	for(int i=1; i<NParticles; i++)
	{
		ODF[i].NOrientations=ODF[0].NOrientations;
		for(int j=0; j<ODF[0].NOrientations; j++)
		{
			ODF[i].hi_SF[j]=ODF[0].hi_SF[j];
			ODF[i].axis_SF[j]=ODF[0].axis_SF[j];
		}
	}
// Read Magnetic parameters into the object MagneticParameter for each paramagnetic particle:
	for(int i=0; i<NParticles; i++)
	{
		MagneticParameter[i].Count=i;
		MagneticParameter[i].MagneticParameterArray=MagneticParameter;
		MagneticParameter[i].ReadFromFile(magnetic_in_filename);
		if(i>0)
		{
// All particles should have the same H0:
			MagneticParameter[i].CenterField0=MagneticParameter[0].CenterField0;
			MagneticParameter[i].Shift=MagneticParameter[0].Shift;
			MagneticParameter[i].ShiftVF=10;
			MagneticParameter[i].ShiftNArg=-1;
// Also all particles will have the same field-net related parameters:
			MagneticParameter[i].NPoints=MagneticParameter[0].NPoints;
			MagneticParameter[i].fs     =MagneticParameter[0].fs;


			if(MagneticParameter[i].AmplitudeFactorVF>=10)
				MagneticParameter[i].AmplitudeCoefficient=MagneticParameter[i].AmplitudeFactor/MagneticParameter[MagneticParameter[i].AmplitudeFactorVF-10].AmplitudeFactor;
			for(int k=0; k<MagneticParameter[i].HyperfineCount; k++)
			{
				if(MagneticParameter[i].HFS[k].HFConstantVF>=10)
				{
					if(MagneticParameter[MagneticParameter[i].HFS[k].HFConstantVF-10].HyperfineCount>k)
						MagneticParameter[i].HFS[k].HFConstantCoefficient=MagneticParameter[i].HFS[k].HFConstant/MagneticParameter[MagneticParameter[i].HFS[k].HFConstantVF-10].HFS[k].HFConstant;
					else
					{
						printf("Warning. Particle %d, Isotropic nucleus %d has improper equivalence. \n", i, k);
						MagneticParameter[i].HFS[k].HFConstantVF=0;
					}
				}
			}
		}
		MagneticParameter[i].GenerateFieldNet();
	}
// Conformation analysis:
	for(int i=0; i<NParticles; i++)
	{
		if(MagneticParameter[i].ConformerGroup != MagneticParameter[i].Count)
		{
			if((MagneticParameter[i].ConformerGroup<0)||(MagneticParameter[i].ConformerGroup>=NParticles))
			{
				printf("Warning! Bad number of conformer group particle: %d. Directive Ignored!\n", MagneticParameter[i].ConformerGroup);
					MagneticParameter[i].ConformerGroup=MagneticParameter[i].Count;
					MagneticParameter[i].NConformers=1;				
			}
			else
			{
				if(MagneticParameter[i].Liouville)
				{
					if(MagneticParameter[MagneticParameter[i].ConformerGroup].Liouville)
					{
						MagneticParameter[i].NConformers=0;
						MagneticParameter[MagneticParameter[i].ConformerGroup].NConformers++; 
						MagneticParameter[MagneticParameter[i].ConformerGroup].ConformerGroupMembers[MagneticParameter[MagneticParameter[i].ConformerGroup].NConformers-1]=i;
						MagneticParameter[i].LOmax=MagneticParameter[MagneticParameter[i].ConformerGroup].LOmax;
						MagneticParameter[i].LEmax=MagneticParameter[MagneticParameter[i].ConformerGroup].LEmax;
						MagneticParameter[i].Mmax=MagneticParameter[MagneticParameter[i].ConformerGroup].Mmax;
						MagneticParameter[i].Kmax=MagneticParameter[MagneticParameter[i].ConformerGroup].Kmax;
						MagneticParameter[i].HGaussVF(0,0)=10+MagneticParameter[i].ConformerGroup;
						MagneticParameter[i].HGaussVF(1,1)=10+MagneticParameter[i].ConformerGroup;
						MagneticParameter[i].HGaussVF(2,2)=10+MagneticParameter[i].ConformerGroup;
						MagneticParameter[i].HGaussEulerVF[0]=10+MagneticParameter[i].ConformerGroup;
						MagneticParameter[i].HGaussEulerVF[1]=10+MagneticParameter[i].ConformerGroup;
						MagneticParameter[i].HGaussEulerVF[2]=10+MagneticParameter[i].ConformerGroup;
					}
					else
					{
						printf("Warning! Non-liouvillian particle %d cannot host a conformer group! Directive ignored!\n", MagneticParameter[i].ConformerGroup);
						MagneticParameter[i].ConformerGroup=MagneticParameter[i].Count;
						MagneticParameter[i].NConformers=1;
					}
				}
				else		// Non-liouvillian particle cannot be one of the conformers!
				{
					printf("Warning! Non-liouvillian particle %d cannot be in a conformer group! Directive ignored!\n", i);
					MagneticParameter[i].ConformerGroup=MagneticParameter[i].Count;
					MagneticParameter[i].NConformers=1;
				}
			}
		}
	}
//	for(int i=0; i<NParticles; i++)
//	{
//		if(MagneticParameter[i].NConformers>1)
//		{
// TODO: Conformation analysis. What else do I need to initialize?
//		}
//	}
//

// Array SpectrumSum is 3-dimensional: 
//  SpectrumSum[i][j][k] is the spectrum intensity for
//      i-th paramagnetic particle
//      j-th orientation of the angular dependency
//      k-th point of the spectrum
//	Thus, the dimensionality of SpetrumSum[][][] is NParticles*NOrientations*NPoints
	SpectrumSum=new double**[NParticles];
	R=new double[ODF[0].NOrientations * MagneticParameter[0].NPoints];

	for(int j=0; j<NParticles; j++)
	{
		SpectrumSum[j]=new double* [ODF[0].NOrientations];
		for (int i=0; i<ODF[0].NOrientations; i++)
			SpectrumSum[j][i]=new double[MagneticParameter[0].NPoints];
	}

// Reading Experimental EPR spectrum and interpolating on the Field net used for simulation.
	Experimental *ExpSpectra=new Experimental[ ODF[0].NOrientations ];
	for(int i=0; i<ODF[0].NOrientations; i++)
	{
		ExpSpectra[i].LoadFromFile(ODF[0].ExperimentalFiles[i]);
		ExpSpectra[i].InterpolateWhole(MagneticParameter[0].FieldNet, MagneticParameter[0].NPoints);
	}

// Angular Dependence and ODF related parameters
	for(int i=0; i<NParticles; i++)
	{
		ODF[i].Count=i;
		ODF[i].ReadFromFile(ODF_in_filename);
		ODF[i].ODFArray=ODF;
	}
// After input all angles will be in radians:
	for(int k=0; k<NParticles; k++)
	{    
		ODF[k].ConvertAnglesToRadians();
		MagneticParameter[k].ConvertAnglesToRadians();

// Generate integration net for particles (domains) orientations:
		if(MagneticParameter[k].Liouville)	MagneticParameter[k].IntegrationNetKsi=1;
		else								MagneticParameter[k].IntegrationNetKsi=MagneticParameter[k].IntegrationNet;

		MagneticParameter[k].GenerateIntegrationNet();
		MagneticParameter[k].GenerateIntegrationNetKsi();

		/*
		double *w1, *w2, *nodes;
		nodes=(double*)malloc(8*sizeof(double));
		w1   =(double*)malloc(8*sizeof(double));
		w2   =(double*)malloc(8*sizeof(double));
		kronrod ( 7, 1.0e-6, nodes, w1, w2 );
  for ( int i = 0; i < 8; i++ )
  {
    printf ( "  %4d  %14f  %14f  %14f\n", i, nodes[i], w1[i], w2[i] );
  }
*/
	}


// Debug 20-01-2015
//	for(double psi=0; psi<=2*PI; psi+=PI/100)
//		printf("%lf   %lf \n", psi, ODF[0].ODF_value(0, PI/2, psi));
//	return 0;

	double Arguments[MAX_ARG];
	double lb[MAX_ARG], ub[MAX_ARG];

	int NArg=0;

// if at least one of the simulated particles is slow-motional (i.e. needs modelling by solving Stochastic Liouville Equation)
// the set of Wigner 3J-symbols should be calculated. In this program Wigner 3J-symbols are pre-calculated and stored for speed increase
//
	bool isAnyLiouville=false;
	int L_max=0;
	for(int i=0; i<NParticles; i++)
	{
		if (MagneticParameter[i].Liouville)
		{
			isAnyLiouville=true;
			if (MagneticParameter[i].NJumpAxes==0)
			{
				if (L_max<max(MagneticParameter[i].LEmax, MagneticParameter[i].LOmax))
					L_max=max(MagneticParameter[i].LEmax, MagneticParameter[i].LOmax);
			}
			else
			{
				if (L_max<2*max(MagneticParameter[i].LEmax, MagneticParameter[i].LOmax))
					L_max=2*max(MagneticParameter[i].LEmax, MagneticParameter[i].LOmax);
			}
			if (L_max<2*MagneticParameter[i].MFPotential.Lmax)
				L_max=2*MagneticParameter[i].MFPotential.Lmax;
			MagneticParameter[i].MFPotential.MagneticParameterParent=&(MagneticParameter[i]);
		}
	}
	Calculate_Wigner3jm(L_max);

// If minimization via nonlinear least-squares is needed, the array of the variated parameters should be generated.
// These variated parameters are referred to as Arguments. The array Arguments[NArg] is double-precision floating point type.
// Method Par2Arg, defined within classes MagneticParameters and OrderParameters, finds which parameters should be variated indepenently
// and creates corresponding array of Arguments. The inverse procedure is programmed in Arg2Par method (see below):
	for(int i=0; i<NParticles; i++)
	{
		MagneticParameter[i].Par2Arg(Arguments, ub, lb, &NArg);
		ODF[i].Par2Arg(Arguments, ub, lb, &NArg);
	}

// preparing data for nonlinear optimization subroutines
	DataForCalculation additionaldata;
	additionaldata.ODF=ODF;
	additionaldata.ExpSpectra=ExpSpectra;
	additionaldata.MagneticParameter=MagneticParameter;
	additionaldata.NParticles=NParticles;

        p_additionaldata=&additionaldata;
	if(MinimizationFlag == 0) NArg=0;

   if(MpiRank==0)
   {

    for (int i=0; i<MpiSize; i++)
	{
		cout << "Process " << i << " calculating:\n";
		for (int k=0; k<NParticles; k++)
			if(MagneticParameter[k].IntegrationNet <= 1)
                     {
                        if(StartMPI(i, ODF[0].NOrientations, MpiSize)<EndMPI(i, ODF[0].NOrientations, MpiSize))
				cout << "\tParticle " << k << ",  spectra  from " << StartMPI(i, ODF[0].NOrientations, MpiSize) << " to " << EndMPI(i, ODF[0].NOrientations, MpiSize)-1 << "\n";
                     }
			else
                     {
                        if(StartMPI(i, MagneticParameter[k].IntegrationNet, MpiSize)<EndMPI(i, MagneticParameter[k].IntegrationNet, MpiSize))
				cout << "\tParticle " << k << ", int.nodes from " << StartMPI(i, MagneticParameter[k].IntegrationNet, MpiSize) << " to " << EndMPI(i, MagneticParameter[k].IntegrationNet, MpiSize)-1 << "\n";
                     }
	}
	cout << endl;
	std::cout << "Argument count: " << NArg << std::endl;
    for(int i=0; i<NArg; i++)
		std::cout << "Argument " << i << ": " << lb[i] << " < " << Arguments[i] << " < " << ub[i] << std::endl;


// START OF MINIMIZATION 

	nlopt_opt opt;
	if(MinimizationFlag == 0) NArg=0;

// Nonlinear Least-Squares minimization with the NL2SOL algorithm:
	if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "NL2SOL"))
	{
		std::cout << "Algorithm NL2SOL\n";
	if(MinimizationFlag != 1) NArg=0;

	int N=ODF[0].NOrientations * MagneticParameter[0].NPoints, P=NArg;
	int* iv;
	double* v;
	iv=new int[60+P];
	v=new double[93+P*N+3*N+P*(3*P+33)];
#ifdef __linux
	stepj_.STEPJAC=MagneticParameter[0].InitialStep;
	stepj_.STEPBOUND=1.0;

	converg_.PARFCTOL = MagneticParameter[0].FTol;
	converg_.PARXCTOL = MagneticParameter[0].XTol;

	keyerr_.keyErr=3;    // Calculate errors!
	kctrl_.KCONTROL=0;
#else
	STEPJ.STEPJAC=MagneticParameter[0].InitialStep;
	STEPJ.STEPBOUND=1.0;

	CONVERG.PARFCTOL = MagneticParameter[0].FTol;
	CONVERG.PARXCTOL = MagneticParameter[0].XTol;

	KEYERR.keyErr=3;    // Calculate errors!
	KCTRL.KCONTROL=0;
#endif
	for(int i=0; i<60+P; i++) iv[i]=0;
	for(int i=0; i<93+P*N+3*N+P*(3*P+33); i++) v[i]=0.0;

	int uiparm[1];
	double urparm[1];
	void *ufparm=(void*)(&additionaldata); 

#ifdef __linux
	dfault_(iv, v);
#else
	DFAULT(iv,v); 
#endif
   	iv[0]=12; 
   	iv[20]=2;
#ifdef __linux
	nl2sno_(&N,&P,Arguments,&madr,iv,v,uiparm,urparm,ufparm);
#else
	NL2SNO(&N,&P,Arguments,&madr,iv,v,uiparm,urparm,ufparm);
#endif

	if (iv[25]<=0)
		printf("Hessian is indefined, iv[25]=%d .\n", iv[25]);
	else
	{
	CovarianceCalculated=true;
	for(int i=0; i<P; i++)
		for(int j=0; j<=i; j++)
			cov[i][j]=v[iv[25]+(i*(i+1))/2+j-1];

	FILE *covariancefile;
	covariancefile=fopen("covariance.out", "w");
	for(int i=0; i<P; i++)
		fprintf(covariancefile, "Arg[%d]=%lf +- %lf \n", i, Arguments[i], sqrt(cov[i][i]));

	fprintf(covariancefile, "Covariance matrix, iv[25]=%d:\n", iv[25]);
	for(int i=0; i<P; i++)
	{
		for(int j=0; j<=i; j++)
//			printf("%d \t", iv[25]+(i*(i+1))/2+j-1);
			fprintf(covariancefile, "%lf \t", cov[i][j]/sqrt(cov[i][i]*cov[j][j]));
		fprintf (covariancefile, "\n");
	}
	fclose(covariancefile);

	}


	}
// Nonlinear Least-Squares minimization with NLOPT library algorithms:
	else
	{

	if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "DIRECT"))
	{
		std::cout << "Algorithm DIRECT-L\n";
		opt = nlopt_create(NLOPT_GN_DIRECT_L , NArg); /* algorithm and dimensionality */
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
	}
	else if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "CRS"))
	{
		std::cout << "Algorithm CRS-2\n";
		opt = nlopt_create(NLOPT_GN_CRS2_LM, NArg); /* algorithm and dimensionality */
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
	}
	else if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "BOBYQA"))
	{
		std::cout << "Algorithm BOBYQA\n";
		opt = nlopt_create(NLOPT_LN_BOBYQA , NArg); /* algorithm and dimensionality */
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
	}
	else if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "COBYLA"))
	{
		std::cout << "Algorithm COBYLA\n";
		opt = nlopt_create(NLOPT_LN_COBYLA , NArg); /* algorithm and dimensionality */
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
	}
	else if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "PRAXIS"))
	{
		std::cout << "Algorithm PRAXIS\n";
		opt = nlopt_create(NLOPT_LN_PRAXIS , NArg); /* algorithm and dimensionality */
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
	}
	else if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "NELDERMEAD"))
	{
		std::cout << "Algorithm NELDERMEAD\n";
		opt = nlopt_create(NLOPT_LN_NELDERMEAD , NArg); /* algorithm and dimensionality */
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
	}
	else if(!strcmp(MagneticParameter[0].NLOptAlgorithm, "SBPLX"))
	{
		std::cout << "Algorithm SBPLX\n";
		opt = nlopt_create(NLOPT_LN_SBPLX , NArg); /* algorithm and dimensionality */
		nlopt_set_lower_bounds(opt, lb);
		nlopt_set_upper_bounds(opt, ub);
	}
	else
	{
		std::cout << "No detected minimization algorithm\n";
		opt = nlopt_create(NLOPT_LN_PRAXIS , NArg); /* algorithm and dimensionality */
		std::cout << "Discrepancy:" << Discrepancy(additionaldata.NParticles, additionaldata.ExpSpectra, additionaldata.MagneticParameter, additionaldata.ODF, Arguments, R) << std::endl;
		goto end_calc;
	}
//	opt = nlopt_create(NLOPT_LN_BOBYQA , NArg); /* algorithm and dimensionality */
	nlopt_set_min_objective(opt, myfunc, &additionaldata);
	nlopt_set_xtol_rel(opt, MagneticParameter[0].XTol);
	nlopt_set_ftol_rel(opt, MagneticParameter[0].FTol);
	nlopt_set_initial_step1(opt, MagneticParameter[0].InitialStep);
	
   

	double minf;

	if (nlopt_optimize(opt, Arguments, &minf) < 0) 
	{
		printf("nlopt failed!\n");
	}
	else {
		printf("found minimum at fmin = %e\n", minf);
	}
	nlopt_destroy(opt);

}

// END OF NLOPT MINIMIZATION
//	std::cout << "Discrepancy:" << Discrepancy(additionaldata.ExpSpectra, additionaldata.MagneticParameter, additionaldata.ODF, Arguments) << std::endl;

end_calc:

   std::cout << "Calculation time: " << (float)(clock()-ct_time)/CLOCKS_PER_SEC << " sec." << std::endl;

/*   for(int i=0; i<MagneticParameter->NPoints; i++)
		{
			std::cout << MagneticParameter->FieldNet[i] << "\t   ";
			for (int k=0; k<ODF->NOrientations; k++)
				std::cout << ExpSpectra[k].SpectrumInterpolated[i] << "\t  " << SpectrumSum[k][i] << "\t   ";
			std::cout << std::endl;
		}
	
			std::cout << std::endl;
*/

// For output all angles will be in degrees:
	for(int k=0; k<NParticles; k++)
	{    
		ODF[k].ConvertAnglesToDegrees();
		MagneticParameter[k].ConvertAnglesToDegrees();
	}
    
    if(!CovarianceCalculated)	WriteOutput(NParticles, ExpSpectra, SpectrumSum, MagneticParameter, ODF, 0,  Arguments);
	else						WriteOutput(NParticles, ExpSpectra, SpectrumSum, MagneticParameter, ODF, NArg, Arguments);
	FILE *outfile;
	outfile=fopen(magnetic_out_filename, "w");
	fclose(outfile);
	outfile=fopen(ODF_out_filename, "w");
	fclose(outfile);
	for(int i=0; i<NParticles; i++)
	{
		MagneticParameter[i].WriteToFile(magnetic_out_filename);
		ODF[i].WriteToFile(ODF_out_filename);
	}

// After output all angles will be in radians:
	for(int k=0; k<NParticles; k++)
	{    
		ODF[k].ConvertAnglesToRadians();
		MagneticParameter[k].ConvertAnglesToRadians();
	}

	outfile=fopen(magnetic_out_filename, "a+");
		fprintf(outfile, "\nVariateFlags:");
		fprintf(outfile, "   1 - variate independently\n");
		fprintf(outfile, "   2 - variate transformed value \n");
		fprintf(outfile, "   3 - calculate one g-tensor value based on constant average and 2 other values.\n");
		fprintf(outfile, "   4 - tensor value is equal to y-value\n");
		fprintf(outfile, "   5 - tensor value is equal to x-value \n");
		fprintf(outfile, "   6 - g2(A2) tensor values are equal to g1(A1)\n");
		fprintf(outfile, "   7 - angles of A-tensor (or linewidth tensor) are inverse to those of g1-tensor (this means, this tensor is in the molecular frame)\n");
		fprintf(outfile, "   7 - angles of Diffusion tensor are equal to those of g1-tensor (this means, diffusion tensor is in g1 frame)\n");
		fprintf(outfile, "   8 - angles of A-tensor (or linewidth tensor) are equal to those of g2-tensor\n\n");
		fprintf(outfile, "Possible Algorithm keywords: COBYLA, BOBYQA, PRAXIS, NELDERMEAD, SBPLX, details see http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms \n\n");

		fprintf(outfile, "Euler angles of all magnetic tensors, except g1 and Diffusion, convert magnetic tensors into g1 tensor.\n");
		fprintf(outfile, "DiffusionEuler angles convert Diffusion frame into molecular orientation frame.\n");
		fprintf(outfile, "g1Euler angles convert g1 frame into molecular orientation frame. D-tensor is coaxial with g1 tensor\n");
		fprintf(outfile, "#IntegrateGauss, IntegrateGaussKronrod and IntegrateVoitlander flags select the orientation integration method: a fixed Gauss-Legendre grid (IntegrateGauss), the same grid plus Kronrod-only nodes for a convergence diff (IntegrateGaussKronrod), or an adaptive hemisphere triangulation (IntegrateVoitlander)\n");
		fprintf(outfile, "\nTo safely stop the calculation on the current iteration, create the file named \"stop\" in the working directory\n");
	fclose(outfile);

        int flag_end=1;
#ifdef __MPI
        MPI_Bcast(&flag_end, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    }
    else     // if (MpiRank!=0)
    {
#ifdef __MPI
      while(1)
      {
// receive next data?
        int flag_end;
        MPI_Bcast(&flag_end, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if(flag_end) break;

// broadcast-receive Arguments
        MPI_Bcast(Arguments, NArg, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// calculate spectrum
        int npts=p_additionaldata->MagneticParameter[0].NPoints*p_additionaldata->ODF[0].NOrientations;
        double* R1=new double[npts];
	double disc=Discrepancy(p_additionaldata->NParticles, p_additionaldata->ExpSpectra, p_additionaldata->MagneticParameter, p_additionaldata->ODF, Arguments, R1);

// the residual buffer belongs to this pass of the slave loop only
        delete[] R1;
// receive data
//        MPI_Gather(R1, npts, MPI_DOUBLE, R1, npts, MPI_DOUBLE, 0, MPI_COMM_WORLD);        
      }
#endif
    }
#ifdef __MPI
    MPI_Finalize();
#endif

	remove ("stop");
	return 0;
}

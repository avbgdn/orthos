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

// File purpose: the objective function of the fit. Discrepancy() builds the theoretical spectra of all particles for the
// current optimiser arguments, compares them with the experimental spectra and returns half the sum of squared
// deviations; madr() and myfunc() are the adapters called by the NL2SOL and NLopt optimisers. The file also holds the
// helper functions that split a loop between MPI ranks and interpolate a spectrum.
// Used by: main.cpp (calls Discrepancy() directly and passes madr and myfunc to the optimisers);
// IntegratedSpectrum.cpp and VoitlanderIntegration.cpp (StartMPI and EndMPI).

#include "BiradicalODF.h"

#ifdef __MPI
#include <mpi.h>
#endif

 // COMMON BLOCKS FOR NL2SOL:
#ifdef __linux
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
		int KCONTROL;
	} KCTRL;
  } 

#endif
//   INTEGER KCONTROL
//   COMMON /KCTRL/ KCONTROL     !AB!

extern double ***SpectrumSum;		// 3-dimensional: NParticles*NOrientations*NPoints
extern DataForCalculation *p_additionaldata;
extern int MpiRank, MpiSize;

// The function calculates the value of discrepancy (0.5 the sum of squared deviations between experimental and theoretical EPR spectra)
// for the given set of magnetic and order parameters.

// (The two comment lines above describe Discrepancy(), defined further down in this file.)
// StartMPI: first index of the slice of N items (orientations, integration nodes or seed triangles) that the MPI
// process Rank handles when the items are divided as evenly as possible between P processes; returns 0 for a rank
// outside 0..P-1.
// Used by: Discrepancy() in this file and CalculateIntegratedSpectrumLiouville() in IntegratedSpectrum.cpp (both only in
// MPI builds); IntegrateVoitlanderMesh() in VoitlanderIntegration.cpp (MPI builds); main() in main.cpp (prints the
// work distribution).
int StartMPI(int Rank, int N, int P)
{
	if (Rank>=P) return 0;
	if (Rank<0) return 0;
	return (Rank*N)/P;
}

// EndMPI: one past the last index of the slice that process Rank handles, so process Rank loops from StartMPI() up to
// but not including EndMPI(); returns 0 for a rank outside 0..P-1.
// Used by: the same callers as StartMPI() (Discrepancy(), CalculateIntegratedSpectrumLiouville(),
// IntegrateVoitlanderMesh(), main()).
int EndMPI(int Rank, int N, int P)
{
	if (Rank>=P) return 0;
	if (Rank<0) return 0;
	return ((Rank+1)*N)/P;
}



// SpectrumInterpolate: linear interpolation of the spectrum Intensity, tabulated on the increasing field net FieldNet
// with NPoints points, at the field x; zero outside the net.
// Used by: Discrepancy() in this file (shifting the spectrum by the isotropic hyperfine multiplet lines).
double SpectrumInterpolate(double *FieldNet, double *Intensity, double x, int NPoints)
{
	if (x<FieldNet[0]) return 0;
	if (x>FieldNet[NPoints-1]) return 0;
	for(int i=0; i<NPoints-1; i++)
	{
		if(x==FieldNet[i]) return Intensity[i];
		if((x>FieldNet[i])&&(x<FieldNet[i+1]))
			return Intensity[i]+(Intensity[i+1]-Intensity[i])*(x-FieldNet[i])/(FieldNet[i+1]-FieldNet[i]);
	}
	if(x==FieldNet[NPoints-1]) return Intensity[NPoints-1];
	return 0;
}



// Discrepancy: objective function of the fit. For the optimiser arguments, the parameters of all particles are set
// up, the spectra of every particle for every sample orientation are computed (Liouville slow-motion or rigid limit),
// broadened by isotropic hyperfine splitting, summed over particles, and compared with the experimental spectra.
// NParticles: number of particles; ExpSpectra: experimental spectra; MagneticParameter: parameters of the particles;
// ODF: orientation distribution of the particles; Arguments: values of the fitted parameters;
// R: output residual vector of length NOrientations*NPoints, R[i*NPoints+j]=(calculated-experimental)*penalty for
// spectrum i and field point j. Returns half the sum of squares of R. The spectra are also left in SpectrumSum and
// the output files are written by rank 0.
// Used by: madr() and myfunc() in this file (through the NL2SOL and NLopt optimisers); main() in main.cpp (directly,
// when no minimisation algorithm is detected, and in the loop of non-root MPI ranks).
double Discrepancy(int NParticles,
				   Experimental *ExpSpectra, 
				   MagneticParameters *MagneticParameter, 
				   OrderParameters *ODF, const double *Arguments,
				   double *R
				  )
{
//
// PUTTING ZEROS TO RESULTING SPECTRUM:
//
	for(int k=0; k<NParticles; k++)
		for(int i=0; i<ODF[k].NOrientations; i++)
			for(int j=0; j<MagneticParameter[0].NPoints; j++)
				SpectrumSum[k][i][j]=0;
//
// FOR EACH PARTICLE:
//
	for(int k=0; k<NParticles; k++)
	{
// before each calculation:
// check for null dipolar vector:
		if(   (fabs(MagneticParameter[k].dx)<=0.1)
			&&(fabs(MagneticParameter[k].dy)<=0.1)
			&&(fabs(MagneticParameter[k].dz)<=0.1))
		{
			std::cout << "\nEncounter of zero dipolar vector" << std::endl;
			MagneticParameter[k].dx=0.11;
			MagneticParameter[k].dy=0.11;
			MagneticParameter[k].dz=0.11;
		}
// 0. Convert arguments back to parameters.
		MagneticParameter[k].Arg2Par(Arguments);
		ODF[k].Arg2Par(Arguments);
// 1. Realize ODF parameters
		ODF[k].Realize();
// 3. recalculate Euler_SF=Euler_SF0 .+ delta_Euler_SF.
		ODF[k].SetEuler_SF();
		ODF[k].ConvertAnglesToRadians();
// 4. Make all Euler Angles radians

// All angles are already Radians
//		ODF[k].ConvertAnglesToRadians();
//		MagneticParameter[k].ConvertAnglesToRadians();
// 2. Recalculate Tjpqi
		ODF[k].AssignTjpqi();
// 6. Field Shift
		MagneticParameter[k].CenterField=MagneticParameter[k].CenterField0+MagneticParameter[k].Shift;
// 5. Rotate all tensors.
// also, rotation diffusion tensor will be converted into gaussian units
		MagneticParameter[k].RotateTensors();
// Apply quasi-librations:
		MagneticParameter[k].QuasiLibrations();
// 5a. Recalculate XLK coefficients:
		if (MagneticParameter[k].Liouville)
		{
			MagneticParameter[k].MFPotential.CalculateOP_Jumps_GK(max(MagneticParameter[k].LOmax, MagneticParameter[k].LEmax));
			MagneticParameter[k].MFPotential.CalculateXLK(2*MagneticParameter[k].MFPotential.Lmax);
			MagneticParameter[k].MFPotential.PrintOut(max(MagneticParameter[k].LOmax, MagneticParameter[k].LEmax));
		}

	}

	// after all particles are converted according to Arg2Par etc. etc. We calculate transition rate and occupances for conformationally mobile radicals
	for(int k=0; k<NParticles; k++)
		MagneticParameter[k].CalculateTranstionRateMatrix();

	// compute the spectra of each particle: orientation integration if IntegrationNet>1, otherwise one spectrum per sample orientation
	for(int k=0; k<NParticles; k++)
	{
		if(MagneticParameter[k].NConformers < 1) continue;   // if Particle is a part of conformer group - its spectrum outside of the group is zero!
// IF INTEGRATION NET GREATER THAN 1:
		if(MagneticParameter[k].IntegrationNet > 1) CalculateIntegratedSpectrumLiouville(&(MagneticParameter[k]), &(ODF[k]), SpectrumSum[k]);
// IF INTEGRATION NET == 1 AND THE SPECTRUM IS IN RIGID LIMIT
//		else if(! MagneticParameter[k].Liouville)   CalculateIntegratedSpectrumLiouville(&(MagneticParameter[k]), &(ODF[k]), SpectrumSum[k]);
// IF MONODOMAIN OR MONOCRYSTAL:
		else 	
		{
// loop over the sample orientations (spectra); under MPI each rank takes its own slice given by StartMPI() and EndMPI()
#ifdef __MPI
	for(int i=StartMPI(MpiRank, ODF[0].NOrientations, MpiSize); i<EndMPI(MpiRank, ODF[0].NOrientations, MpiSize); i++)
#else
	for(int i=0; i<ODF[0].NOrientations; i++)
#endif	
			{
// IF SLOW MOTIONAL FREED AND MONODOMAIN:
				if(MagneticParameter[k].Liouville)
				{

//					cout << ODF->beta_SF[i] << endl;
					// slow-motion (stochastic Liouville) spectrum of this particle for this sample orientation
					Matrix<double, Dynamic, 1> Spectrum;
					Spectrum.resize(MagneticParameter[0].NPoints , 1);
			// for 0-th particle the basis vectors can be read from file
					if(k==0)
						Spectrum=CalculateLiouvillianB(& (MagneticParameter[k]), ODF[k].beta_SF[i], ODF[k].BasisSource[i] );
			// otherwise, the basis vectors are generated from L,M,K truncation values and selection rules.
					else
						Spectrum=CalculateLiouvillianB(& (MagneticParameter[k]), ODF[k].beta_SF[i], "truncation" );

					// store the spectrum in the accumulator of this particle and orientation
					for(int j=0; j<MagneticParameter[0].NPoints; j++)
						SpectrumSum[k][i][j]=Spectrum(j, 0);
				}
// IF RIGID LIMIT AND MONOCRYSTAL:
				else
				{
				// Initializing magnetic field net related parameters
					double *Spectrum;
					Spectrum=(double*)malloc(MagneticParameter[0].NPoints*sizeof(double));

					double Field_to_Molecule_Euler[3];
					Matrix<double, 3, 3> RotMat;
					double HGaussEff, HLorentzEff;

					// Euler angles field->molecule: the negated sample-frame angles gamma, beta, alpha (in this order), and the rotation matrix for them
					Field_to_Molecule_Euler[0]=-(ODF[k].gamma_SF[i]);
					Field_to_Molecule_Euler[1]=-(ODF[k].beta_SF [i]);
					Field_to_Molecule_Euler[2]=-(ODF[k].alpha_SF[i]);

					RotMat=RotationMatrix(Field_to_Molecule_Euler);
					// biradical: 54 transitions, monoradical: 3 transitions; the monoradical flag is read from particle 0 (MagneticParameter->, not [k])
					if(!(MagneticParameter->isMonoradical))
					{
					// Solving Spin-Hamiltonian and calculating ESR transitions for given orientation
						Array<double, 54, 2> Transitions=CalculateTransitions(&(MagneticParameter[k]), RotMat, &HGaussEff, &HLorentzEff);
					// Calculating spectrum lineshape for given orientation
						CalculateIndividualSpectrum (Transitions, HGaussEff, HLorentzEff, MagneticParameter[0].NPoints, MagneticParameter[0].FieldNet, Spectrum, MagneticParameter[k].TransitionIntensityCutoff);
					}
					else
					{
					// Solving Spin-Hamiltonian and calculating ESR transitions for given orientation
						Array<double, 3, 2> Transitions=CalculateTransitionsM(&(MagneticParameter[k]), RotMat, &HGaussEff, &HLorentzEff);
					// Calculating spectrum lineshape for given orientation
						CalculateIndividualSpectrumM (Transitions, HGaussEff, HLorentzEff, MagneticParameter[0].NPoints, MagneticParameter[0].FieldNet, Spectrum, MagneticParameter[k].TransitionIntensityCutoff);
					}
				


					// scale by AmplitudeFactor/36 and store in the accumulator of this particle and orientation
					for(int j=0; j<MagneticParameter[0].NPoints; j++)
						SpectrumSum[k][i][j]=MagneticParameter[k].AmplitudeFactor*Spectrum[j]/36.0;
				}
			}
		}



// All angles will be stored as radians!
//		MagneticParameter[k].ConvertAnglesToDegrees();

	}  // end for-cycle by k - particle number

	double discrep=0;
// Calculating penalty coefficient for negative volume of orientation distribution function (which should be non-negative).
	double PenaltyCoefficient=1.0;
	for (int k=0; k<NParticles; k++)
		PenaltyCoefficient+=ODF[k].Penalty*ODF[k].NegativeVolume;

#ifdef __MPI

	// collect the partial spectra of all ranks on rank 0 (each rank has filled only its own share, the rest is zero)
	if(MpiRank==0)
	{
		double RecvBuffer[NParticles*ODF[0].NOrientations*MagneticParameter[0].NPoints*MpiSize];
		double SendBuffer[NParticles*ODF[0].NOrientations*MagneticParameter[0].NPoints];

		for(int k=0; k<NParticles; k++)
			for (int i=0; i<ODF[k].NOrientations; i++)
				for(int j=0; j<MagneticParameter[0].NPoints; j++)
				  SendBuffer[j+i*MagneticParameter[0].NPoints+k*MagneticParameter[0].NPoints*ODF[0].NOrientations]=SpectrumSum[k][i][j];
		MPI_Gather(SendBuffer, NParticles*ODF[0].NOrientations*MagneticParameter[0].NPoints, MPI_DOUBLE, RecvBuffer, NParticles*ODF[0].NOrientations*MagneticParameter[0].NPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);        

		// sum the contributions of all ranks for every particle, orientation and field point
		for(int k=0; k<NParticles; k++)
			for (int i=0; i<ODF[k].NOrientations; i++)
				for(int j=0; j<MagneticParameter[0].NPoints; j++)
				{
					double summ=0;
					for (int n=0; n<MpiSize; n++)
						summ+=RecvBuffer[j+i*MagneticParameter[0].NPoints+k*MagneticParameter[0].NPoints*ODF[0].NOrientations+n*MagneticParameter[0].NPoints*ODF[0].NOrientations*NParticles];

					SpectrumSum[k][i][j]=summ;
				}
	
	}
	else
	{
		// the other ranks pack their partial spectra and send them to rank 0
		double SendBuffer[NParticles*ODF[0].NOrientations*MagneticParameter[0].NPoints];

		for(int k=0; k<NParticles; k++)
			for (int i=0; i<ODF[k].NOrientations; i++)
				for(int j=0; j<MagneticParameter[0].NPoints; j++)
				  SendBuffer[j+i*MagneticParameter[0].NPoints+k*MagneticParameter[0].NPoints*ODF[0].NOrientations]=SpectrumSum[k][i][j];
		MPI_Gather(SendBuffer, NParticles*ODF[0].NOrientations*MagneticParameter[0].NPoints, MPI_DOUBLE, SendBuffer, NParticles*ODF[0].NOrientations*MagneticParameter[0].NPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);        
	}

#endif

// TIME TO APPLY ISOTROPIC HFS!
// i - Spectrum No
// k - Particle No
// j - Point No
	for(int k=0; k<NParticles; k++)
	{
	// stick spectrum of the isotropic hyperfine splitting of this particle: line shifts relative to zero and relative intensities
	MagneticParameter[k].Mult.SetFromHFS(0, MagneticParameter[k].HFS, MagneticParameter[k].HyperfineCount);
	if (MagneticParameter[k].Mult.NLines>1)
		for (int i=0; i<ODF[0].NOrientations; i++)
		{
			// convolve each spectrum with the stick spectrum: sum of copies shifted by each line position (read by interpolation) and weighted by its intensity
			double *tempspec=new double[MagneticParameter[0].NPoints];
			for(int j=0; j<MagneticParameter[0].NPoints; j++)
				tempspec[j]=0;
			
			for (int m=0; m<MagneticParameter[k].Mult.NLines; m++)
			{
				for(int j=0; j<MagneticParameter[0].NPoints; j++)
					tempspec[j]+=SpectrumInterpolate(MagneticParameter[0].FieldNet, SpectrumSum[k][i], MagneticParameter[0].FieldNet[j]-MagneticParameter[k].Mult.Transitions(m, 0), MagneticParameter[0].NPoints)*MagneticParameter[k].Mult.Transitions(m, 1);
			}

			// replace the spectrum by the broadened one
			for(int j=0; j<MagneticParameter[0].NPoints; j++)
				SpectrumSum[k][i][j]=tempspec[j];
		}
	}

	if(MagneticParameter[0].AmplitudeNormalized)
		// NORMALIZE CALCULATED SPECTRUM ON AMPLITUDE
	{
		for (int i=0; i<ODF[0].NOrientations; i++)
		{
			// largest and smallest value of the summed spectrum of all particles
			double maxpt=0, minpt=0;
			for(int j=0; j<MagneticParameter[0].NPoints; j++)
			{
				double sum=0;
				for(int k=0; k<NParticles; k++)
					sum+=SpectrumSum[k][i][j]; 
				if (maxpt<sum) maxpt=sum;
				if (minpt>sum) minpt=sum;
			}

			// divide the spectra of all particles by the peak-to-peak amplitude of their sum
			for(int j=0; j<MagneticParameter[0].NPoints; j++)
				for(int k=0; k<NParticles; k++)
					SpectrumSum[k][i][j] /= (maxpt-minpt); 

		}
	}

	// residual vector: summed theoretical minus experimental spectrum, times the penalty coefficient, and the sum of its squares
	for (int i=0; i<ODF[0].NOrientations; i++)
		for(int j=0; j<MagneticParameter[0].NPoints; j++)
		{
			double sum=0;
			for(int k=0; k<NParticles; k++)
			{
			     sum+=SpectrumSum[k][i][j]; 
			}
			R[i*MagneticParameter[0].NPoints+j]=(sum-ExpSpectra[i].SpectrumInterpolated[j])*PenaltyCoefficient;
			discrep+=pow( (sum-ExpSpectra[i].SpectrumInterpolated[j])*PenaltyCoefficient, 2);
		}
	// discrepancy is half the sum of squared residuals, overall and per point
	MagneticParameter[0].CurrentDiscrepancyOverall=discrep/2.0;
	MagneticParameter[0].CurrentDiscrepancyPerPoint=discrep/(2.0*MagneticParameter[0].NPoints*ODF[0].NOrientations);
   // only the root rank writes the output files (angles are changed to degrees for output and back to radians afterwards)
   if(MpiRank==0)
   {
// For output all angles will be in degrees:
	for(int k=0; k<NParticles; k++)
	{    
		ODF[k].ConvertAnglesToDegrees();
		MagneticParameter[k].ConvertAnglesToDegrees();
	}
//	cout << "A200 = " << ODF[0].A_jmn(2, 0, 0) << endl;

// Generating output files:x	
	WriteOutput(NParticles, ExpSpectra, SpectrumSum, MagneticParameter, ODF, 0, Arguments);
	// empty the temporary parameter files and write the current parameters of every particle into them
	FILE *outfile;
	outfile=fopen("Magnetic.temp", "w");
	fclose(outfile);
	outfile=fopen("ODF.temp", "w");
	fclose(outfile);
	for(int k=0; k<NParticles; k++)
	{
		MagneticParameter[k].WriteToFile("Magnetic.temp");
		ODF[k].WriteToFile("ODF.temp");
	}
// After output all angles will be in radians:
	for(int k=0; k<NParticles; k++)
	{    
		ODF[k].ConvertAnglesToRadians();
		MagneticParameter[k].ConvertAnglesToRadians();
	}
   }
	return discrep/2.0;
}

// interface for NL2SOL:
// madr: adapter between NL2SOL and Discrepancy(). NL2SOL calls it with the argument vector X (P values) and expects the
// residual vector in R; the problem data come through the pointer UFPARM (a DataForCalculation structure).
// The other arguments (N, NF, UIPARM, URPARM, IV) are part of the NL2SOL calling convention and are not used.
// A file named stop in the working directory makes NL2SOL stop at its next test. After the evaluation X and half the
// sum of squared residuals are appended to the file discrepancy. Under MPI, rank 0 first broadcasts the arguments to the
// other ranks, which are waiting in the loop in main().
// Used by: main() in main.cpp, as a function pointer passed to nl2sno_ (or NL2SNO) of NL2SOL1.for, which calls it back.
void *madr(int* N, int* P, double* X, int* NF, double* R, int* UIPARM, double* URPARM, void* UFPARM, int* IV)
{
	// recover the problem data from the user pointer
	DataForCalculation *additionaldata;
	additionaldata=(DataForCalculation *)UFPARM;
	FILE *KControlFile;

	// if the file named stop exists, set the NL2SOL control flag KCONTROL to 1, which makes the NL2SOL stop test true
	KControlFile=fopen("stop", "r");
	if (KControlFile != NULL) 
	{
#ifdef __linux
	  kctrl_.KCONTROL=1;
#else
	  KCTRL.KCONTROL=1;
#endif
	fclose(KControlFile);
	}

//	std::cout << "Arguments: \n";
//	for(int i=0; i<(int)n; i++)
//		std::cout << x[i] << "   ";
#ifdef __MPI
// receive next data
        int flag_end=0;
        MPI_Bcast(&flag_end, 1, MPI_INT, 0, MPI_COMM_WORLD);
// broadcast Arguments
        MPI_Bcast(X, *P, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// calculate spectrum at root
        int npts=additionaldata->MagneticParameter[0].NPoints;
//        double* R1=new double[npts];
	double disc=Discrepancy(additionaldata->NParticles, additionaldata->ExpSpectra, additionaldata->MagneticParameter, additionaldata->ODF, X, R);
// receive data
//        MPI_Gather(R1, npts, MPI_DOUBLE, R, npts, MPI_DOUBLE, 0, MPI_COMM_WORLD);        
#else
	double disc=Discrepancy(additionaldata->NParticles, additionaldata->ExpSpectra, additionaldata->MagneticParameter, additionaldata->ODF, X, R);
#endif
//	std::cout << "\nDiscrepancy:" << disc << std::endl;
//    printf("!");
	// progress mark, then log the arguments and half the sum of squared residuals to the file discrepancy
	cout << "!";
//	cout << IV[3];
	std::fstream discr;
	discr.open("discrepancy", ios::app);
//	discr << 0.5*disc << endl;
	for (int f=0; f<*P; f++)
		discr << X[f] << "\t";
	discr << endl;
	double sum=0;
         for(int f=0; f<additionaldata->ODF[0].NOrientations*additionaldata->MagneticParameter[0].NPoints; f++)
	{
		sum+=R[f]*R[f];
//		discr << R[f] << endl;
	}
	discr << 0.5*sum << endl << endl;
	discr.close();
	return 0;
}

// interface for NLOPT:
// myfunc: adapter between NLopt and Discrepancy(). NLopt calls it with n arguments x and expects the discrepancy as the
// return value; the problem data come through my_func_data (a DataForCalculation structure). grad is never filled in.
// Under MPI, rank 0 first broadcasts the arguments to the other ranks.
// Used by: main() in main.cpp, as a function pointer given to nlopt_set_min_objective, which calls it back.
double myfunc(unsigned n, const double *x, double *grad, void *my_func_data)
{
		// recover the problem data and allocate the residual vector R and a copy X of the arguments
		DataForCalculation *additionaldata;
		additionaldata=(DataForCalculation *)my_func_data;
		double *R, *X;

		R=new double[additionaldata->ODF[0].NOrientations * additionaldata->MagneticParameter[0].NPoints];
		X=new double[n];

		std::cout << "Arguments: \n";
		// print the arguments and copy them to X
		for(int i=0; i<(int)n; i++)
              {
			std::cout << x[i] << "   ";
                     X[i]=x[i];
              }
#ifdef __MPI
// receive next data
        int flag_end=0;
        MPI_Bcast(&flag_end, 1, MPI_INT, 0, MPI_COMM_WORLD);
// broadcast Arguments
        MPI_Bcast(X, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// calculate spectrum at root
	double disc=Discrepancy(additionaldata->NParticles, additionaldata->ExpSpectra, additionaldata->MagneticParameter, additionaldata->ODF, x, R);
#else
	double disc=Discrepancy(additionaldata->NParticles, additionaldata->ExpSpectra, additionaldata->MagneticParameter, additionaldata->ODF, x, R);
#endif
		std::cout << "\nDiscrepancy:" << disc << std::endl;

// the residual and argument buffers belong to this objective evaluation only
		delete[] R;
		delete[] X;

		return disc;

}

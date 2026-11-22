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

// Reading, writing and preparing the magnetic parameters of one paramagnetic particle. This file implements the
// methods of class MagneticParameters (declared in BiradicalODF.h): parsing of one Particle: section of the magnetic
// input file, output of the parameters, conversion of angles between degrees and radians, conversion between the
// physical parameters and the argument vector of the optimiser (Par2Arg, Arg2Par), rotation of the interaction
// tensors into the molecular frame, quasi-librations, the conformer rate matrix and the field net.
// Used by: main.cpp (ReadFromFile, GenerateFieldNet, Par2Arg, angle conversion, WriteToFile), Discrepancy.cpp
// (Arg2Par, RotateTensors, QuasiLibrations, CalculateTranstionRateMatrix, angle conversion, WriteToFile),
// ReadExperimental.cpp (WriteToFileBrief) and ODF.cpp (NextArg). The members filled here are read by the
// spectrum calculation (LiouvillianB.cpp, IntegratedSpectrum.cpp, Discrepancy.cpp and others) through BiradicalODF.h.

#include "BiradicalODF.h"
#include <math.h>

// Reading/writing magnetic parameters file

// Reads the parameters of this particle (number Count) from the magnetic input file. The file holds one section per
// particle, started by a line Particle: N. The first loop skips to the section with N equal to Count, the second loop
// parses the lines Keyword: values of that section until the section of another particle starts. Most keywords carry
// pairs of a value and a variate flag (VF): the flag tells Par2Arg whether the value is fixed (0), fitted (1 or 2)
// or tied to another parameter (3 and above, see Arg2Par). Angles are read in degrees; main() converts them to
// radians after reading. Liouville is set to true (slow-motion calculation with the stochastic Liouville equation)
// when any rotational diffusion rate Rx, Ry, Rz is positive, otherwise the rigid-limit calculation is used.
// Input: filename, the magnetic input file. Outputs: the members of this object; the program exits if the file
// cannot be opened or NPoints is negative.
// Used by: main() in main.cpp, once for each particle.
void MagneticParameters::ReadFromFile(char* filename)
{
// Default values of rotation parameters: for
// backwards-compatibility with BiracicalODF2
		Rx=0.0;
		Ry=0.0;
		Rz=0.0;

		RxVF=0;
		RyVF=0;
		RzVF=0;

		NJumpAxes=0;

		DiffusionEuler[0]=0.0;
		DiffusionEuler[1]=0.0;
		DiffusionEuler[2]=0.0;

		DiffusionEulerVF[0]=0;
		DiffusionEulerVF[1]=0;
		DiffusionEulerVF[2]=0;

		DiffusionType[0]=0.0;
		DiffusionType[1]=0.0;
		DiffusionType[2]=0.0;

		RTau[0]=0.0;
		RTau[0]=0.0;
		RTau[0]=0.0;

		// Allocate the conformer conversion-rate matrix, set all rates to zero and mark them as not varied
		ConversionRateR=(double*)malloc(MAX_PARTICLES*MAX_PARTICLES*sizeof(double));

		for(int i=0; i<MAX_PARTICLES; i++)
			for(int j=0; j<MAX_PARTICLES; j++)
			{
				ConversionRate[i][j]=0.0;
				ConversionRateR[i*MAX_PARTICLES+j]=0.0;
				ConversionRateVF[i][j]=0;
				ConversionRateNArg[i][j]=-1;
			}
		// By default a particle is a conformer group of its own
		ConformerGroup=Count;
		NConformers=1;
		ConformerGroupMembers[0]=Count;

		// Defaults for the basis truncation limits, pruning and the number of Lanczos steps
		LOmax=0; LEmax=0; Mmax=0; Kmax=0; Nmax=0;
		NPrune=0; dHPrune=64.0; PruningTolerance=1.0e-3;
		LanczosSteps=500;
//
	// Defaults: no hyperfine nuclei, no zero-field splitting, no exchange and no mean-field potential
	HyperfineCount=0;
	 
	ZFS_E=0; ZFS_Angle=0;
	ZFS_EVF=0; ZFS_AngleVF=0;
	J0=0; J2=0; JDiffusion=0;


    MFPotential.Lmax=0;

	// Open the magnetic input file; stop if it cannot be opened
	ifstream magnetic_par;
	magnetic_par.open (filename, ios::in);
	if(!magnetic_par.is_open())
	{
		std::cout << "Could not open file " << filename << std::endl;
		exit(0);
	}
	// First pass: skip to the section Particle: Count, setting the default optimiser parameters on the way
	while (!magnetic_par.eof())
	{
		char buffer[MAX_STRLEN], buffer2[MAX_STRLEN], buffer3[MAX_STRLEN];

		XTol=1.0e-7;
		FTol=1.0e-7;
		InitialStep=1.0;
		strcpy(NLOptAlgorithm, "BOBYQA");

		magnetic_par.getline(buffer, MAX_STRLEN);
		if( sscanf(buffer, "%s", buffer2)==-1) continue;
		int npart;
		if (!strcmp(buffer2, "Particle:"))
		{
			sscanf(buffer, "%s %d", buffer3, &npart);
			if(npart==Count) break;
		}
	}
	// Defaults for the switches and the integration method that the section may change
	isMonoradical=false;
	IntegrationMethod=INTEGRATE_GAUSSKRONROD;
	AmplitudeNormalized=false;

	// Second pass: read the keyword lines of this particle until the section of another particle starts
	while (!magnetic_par.eof())
	{
		char buffer[MAX_STRLEN], buffer2[MAX_STRLEN], buffer3[MAX_STRLEN];

		magnetic_par.getline(buffer, MAX_STRLEN);
		if( sscanf(buffer, "%s", buffer2)==-1) continue;
		// Keyword Monoradical: a single electron spin, without second spin, dipolar, exchange and ZFS terms
		if (!strcmp(buffer2, "Monoradical"))
		{
			isMonoradical=true;
		}
		// Principal values of the g1, g2, A1 and A2 tensors and of the Gaussian and Lorentzian line-width tensors
		// (FWHM in Gauss), each followed by its variate flag
		if (!strcmp(buffer2, "g1Tensor:"))
		{
			double a, b, c;
			int i, j, k;
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, &a, &i, &b, &j, &c, &k);
			g1 << a, 0, 0,
				  0, b, 0,
				  0, 0, c;
			g1VF << i, 0, 0,
				    0, j, 0,
				    0, 0, k;
		}
		if (!strcmp(buffer2, "g2Tensor:"))
		{
			double a, b, c;
			int i, j, k;
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, &a, &i, &b, &j, &c, &k);
			g2 << a, 0, 0,
				  0, b, 0,
				  0, 0, c;
			g2VF << i, 0, 0,
				    0, j, 0,
				    0, 0, k;
		}
		if (!strcmp(buffer2, "A1Tensor:"))
		{
			double a, b, c;
			int i, j, k;
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, &a, &i, &b, &j, &c, &k);
			A1 << a, 0, 0,
				  0, b, 0,
				  0, 0, c;
			A1VF << i, 0, 0,
				    0, j, 0,
				    0, 0, k;
		}
		if (!strcmp(buffer2, "A2Tensor:"))
		{
			double a, b, c;
			int i, j, k;
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, &a, &i, &b, &j, &c, &k);
			A2 << a, 0, 0,
				  0, b, 0,
				  0, 0, c;
			A2VF << i, 0, 0,
				    0, j, 0,
				    0, 0, k;
		}
		if (!strcmp(buffer2, "GaussTensor:"))
		{
			double a, b, c;
			int i, j, k;
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, &a, &i, &b, &j, &c, &k);
			HGauss << a, 0, 0,
					  0, b, 0,
					  0, 0, c;
			HGaussVF << i, 0, 0,
			    	    0, j, 0,
						0, 0, k;
		}
		if (!strcmp(buffer2, "LorentzTensor:"))
		{
			double a, b, c;
			int i, j, k;
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, &a, &i, &b, &j, &c, &k);
			HLorentz << a, 0, 0,
					  0, b, 0,
					  0, 0, c;
			HLorentzVF << i, 0, 0,
			    	    0, j, 0,
						0, 0, k;
		}
		// Euler angles (degrees in the file) that orient the principal axes of each tensor, as pairs of value and flag
		if (!strcmp(buffer2, "g1Euler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(g1Euler[0]), &(g1EulerVF[0]), 
													&(g1Euler[1]), &(g1EulerVF[1]), 
													&(g1Euler[2]), &(g1EulerVF[2]));
		}
		if (!strcmp(buffer2, "g2Euler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(g2Euler[0]), &(g2EulerVF[0]), 
													&(g2Euler[1]), &(g2EulerVF[1]), 
													&(g2Euler[2]), &(g2EulerVF[2]));
		}
		if (!strcmp(buffer2, "A1Euler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(A1Euler[0]), &(A1EulerVF[0]), 
													&(A1Euler[1]), &(A1EulerVF[1]), 
													&(A1Euler[2]), &(A1EulerVF[2]));
		}
		if (!strcmp(buffer2, "A2Euler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(A2Euler[0]), &(A2EulerVF[0]), 
													&(A2Euler[1]), &(A2EulerVF[1]), 
													&(A2Euler[2]), &(A2EulerVF[2]));
		}
		if (!strcmp(buffer2, "GaussEuler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(HGaussEuler[0]), &(HGaussEulerVF[0]), 
													&(HGaussEuler[1]), &(HGaussEulerVF[1]), 
													&(HGaussEuler[2]), &(HGaussEulerVF[2]));
		}
		if (!strcmp(buffer2, "LorentzEuler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(HLorentzEuler[0]), &(HLorentzEulerVF[0]), 
													&(HLorentzEuler[1]), &(HLorentzEulerVF[1]), 
													&(HLorentzEuler[2]), &(HLorentzEulerVF[2]));
		}
		// Dipolar vector between the two electron spins, which defines the dipole-dipole tensor
		if (!strcmp(buffer2, "Dipolar_vector:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(dx), &(dxVF), 
													&(dy), &(dyVF), 
													&(dz), &(dzVF));
		}

		// Zero-field splitting: parameter E of its non-axial part and the angle ZFS_Angle
		if (!strcmp(buffer2, "ZFS:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d", buffer3, 
													&(ZFS_E), &(ZFS_EVF), 
													&(ZFS_Angle), &(ZFS_AngleVF));
		}

		// Quasi-libration amplitudes about three axes and Euler angles of the libration frame (used by QuasiLibrations);
		// LibDistrWidth and LibDistrNSites are read and written here but are not used by QuasiLibrations
		if (!strcmp(buffer2, "Libration:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(LibrationAmplitude[0]), &(LibrationAmplitudeVF[0]), 
													&(LibrationAmplitude[1]), &(LibrationAmplitudeVF[1]), 
													&(LibrationAmplitude[2]), &(LibrationAmplitudeVF[2]));
		}
		if (!strcmp(buffer2, "LibrationEuler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(LibrationEuler[0]), &(LibrationEulerVF[0]), 
													&(LibrationEuler[1]), &(LibrationEulerVF[1]), 
													&(LibrationEuler[2]), &(LibrationEulerVF[2]));
		}
		if (!strcmp(buffer2, "LibDistrWidth:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(LibDistrWidth[0]), &(LibDistrWidthVF[0]), 
													&(LibDistrWidth[1]), &(LibDistrWidthVF[1]), 
													&(LibDistrWidth[2]), &(LibDistrWidthVF[2]));
		}
		if (!strcmp(buffer2, "LibDistrNSites:"))
		{
			sscanf(buffer, "%s %d", buffer3, 
													&(LibDistrNSites));
		}


		// Rotational diffusion rates Rx, Ry, Rz (1/s, converted to Gauss units in RotateTensors) and their axes
		if (!strcmp(buffer2, "Diffusion:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(Rx), &(RxVF), 
													&(Ry), &(RyVF), 
													&(Rz), &(RzVF));
		}
		if (!strcmp(buffer2, "DiffusionEuler:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(DiffusionEuler[0]), &(DiffusionEulerVF[0]), 
													&(DiffusionEuler[1]), &(DiffusionEulerVF[1]), 
													&(DiffusionEuler[2]), &(DiffusionEulerVF[2]));
		}

		// Diffusion type of each principal axis; a value other than 0, 0.5 or 1 is replaced by 0 (Brownian diffusion)
		if (!strcmp(buffer2, "DiffusionType:"))
		{
			sscanf(buffer, "%s %lf %lf %lf", buffer3, 
													&(DiffusionType[0]), 
													&(DiffusionType[1]), 
													&(DiffusionType[2]));
          for(int i=0; i<3; i++)		// check for correct Diffusion Type
		  {
			if ((DiffusionType[i]!=0.0)&&
				(DiffusionType[i]!=0.5)&&
				(DiffusionType[i]!=1.0))
			{
				printf("Wrong DiffusionType[%d], setting to 0 (Brownian diffusion)\n", i);
				DiffusionType[i]=0.0;
			}
		  }
		}
		// Residence times RTau of the jump-diffusion model, one for each axis
		if (!strcmp(buffer2, "RTau:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
													&(RTau[0]), &(RTauVF[0]), 
													&(RTau[1]), &(RTauVF[1]), 
													&(RTau[2]), &(RTauVF[2]));
		}


		// Mean-field orientational potential: lines of L, K, value and flag up to a blank line; Potential: no means none
		if (!strcmp(buffer2, "Potential:"))
		{
			char buffer4[MAX_STRLEN];
			if (sscanf(buffer, "%s %s", buffer3, buffer4) < 2) strcpy(buffer4, "yes");
			if(strcmp(buffer4, "no"))
			{
				while(magnetic_par.getline(buffer3, MAX_STRLEN))
				{
					if((strlen(buffer3)==0) || (!strcmp(buffer3, "\r"))) break;
					int L, K, VF;
					double Value;
					sscanf(buffer3, "%d %d %lf %d", &L, &K, &Value, &VF);
					MFPotential.Add(L, K, Value, VF);
				}
			}
		sscanf(buffer, "%s", buffer2);

		}

		// Integration net (theta and phi points) for the potential
		if (!strcmp(buffer2, "NetPotential:"))
		{
			sscanf(buffer, "%s %d %d", buffer3, 
										&(MFPotential.INT_NET), &(MFPotential.INT_NET_FI) );
		}


		// Basis truncation limits LE, LO, K, M, N (used by GenerateBasisB and GenerateBasisM)
		if (!strcmp(buffer2, "Truncation:"))
		{
			sscanf(buffer, "%s %d %d %d %d %d", buffer3, 
													&(LEmax), &(LOmax), 
													&(Kmax), &(Mmax), &(Nmax));
		}
		// Basis pruning settings: NPrune, dHPrune and PruningTolerance
		if (!strcmp(buffer2, "Pruning:"))
		{
			sscanf(buffer, "%s %d %lf %lf ", buffer3, 
													&(NPrune), &(dHPrune), 
													&(PruningTolerance));
		}
		// One discrete-jump axis per line: axis angles fi and teta, jump angle and jump frequency, each with a flag
		if (!strcmp(buffer2, "Jumps:"))
		{
			if(NJumpAxes==MAX_JUMPS)
			{
				printf("Maximum number of discrete jump axes for Particle #%d exceeded, further jumps are ignored!\n", Count);
			}
			else
			{
				sscanf(buffer, "%s %lf %d %lf %d %lf %d %lf %d", buffer3, 
													&(JumpAxisFi   [NJumpAxes]), &(JumpAxisFiVF   [NJumpAxes]), 
													&(JumpAxisTeta [NJumpAxes]), &(JumpAxisTetaVF [NJumpAxes]), 
													&(JumpAngle    [NJumpAxes]), &(JumpAngleVF    [NJumpAxes]),
													&(JumpFrequency[NJumpAxes]), &(JumpFrequencyVF[NJumpAxes]));
				NJumpAxes++;
			}
		}
		// Number of Lanczos steps
		if (!strcmp(buffer2, "LanczosSteps:"))
		{
			sscanf(buffer, "%s %d ", buffer3, &LanczosSteps);
		}


		// Exchange parameters J0, J2 and JDiffusion; a negative J0 is stored as its modulus with the sign kept in J0sign
		if (!strcmp(buffer2, "Exchange:"))
		{
			sscanf(buffer, "%s %lf %d %lf %d %lf %d", buffer3, 
										&(J0), &(J0VF), &(J2), &(J2VF), &(JDiffusion), &(JDiffusionVF));
			if (J0>=0)  J0sign=1.0;
			else		
			{
				J0=-J0;
				J0sign=-1.0;
			}

		}	


		// Isotropic hyperfine nuclei: lines of nucleus name, number of equivalent nuclei, constant and flag up to a blank line
		if (!strcmp(buffer2, "Hyperfine:"))
		{
			char buffer4[MAX_STRLEN];
			if (sscanf(buffer, "%s %s", buffer3, buffer4) < 2) strcpy(buffer4, "yes");
			
//			printf("\n\n_%s_ _%s_\n\n", buffer3, buffer4);
			if(strcmp(buffer4, "no"))
			{
//				printf("HERE!\n");
				while(magnetic_par.getline(buffer3, MAX_STRLEN))
				{
					if((strlen(buffer3)==0) || (!strcmp(buffer3, "\r"))) break;
//					char atomtype[MAX_STRLEN];
					double Value;
					sscanf(buffer3, "%s %d %lf %d", &(HFS[HyperfineCount].AtomType), &(HFS[HyperfineCount].Degeneracy), &(HFS[HyperfineCount].HFConstant), &(HFS[HyperfineCount].HFConstantVF));
					HyperfineCount++;
					

				}
			}
		sscanf(buffer, "%s", buffer2);

		}


		// Centre field of the spectrum, and the field shift Shift with its variate flag
		if (!strcmp(buffer2, "FCenter:"))
		{
			sscanf(buffer, "%s %lf %lf %d ", buffer3, 
											&CenterField0, &(Shift), &(ShiftVF));
		}	
		// Overall amplitude factor of the spectrum of this particle
		if (!strcmp(buffer2, "Amplitude:"))
		{
			sscanf(buffer, "%s %lf %d ", buffer3, 
											&(AmplitudeFactor), &(AmplitudeFactorVF));
		}

		// Field net: number of points and field step DeltaH in Gauss
		if (!strcmp(buffer2, "NPoints:"))
		{
			sscanf(buffer, "%s %d ", buffer3, &NPoints);
			if(NPoints<0)
			{
				std::cout << "Error reading " << filename << ": NPoints=" << NPoints << " is negative." << std::endl;
				exit(0);
			}
		}
		if (!strcmp(buffer2, "DeltaH:"))
		{
			sscanf(buffer, "%s %lf ", buffer3, &fs);
		}
		// Orientation integration: size of the grid and the method (Gauss, Gauss-Kronrod or Voitlander)
		if (!strcmp(buffer2, "Integration_Net:"))
		{
			sscanf(buffer, "%s %d ", buffer3, &IntegrationNet);
		}

		if (!strcmp(buffer2, "IntegrateGauss"))
		{
//			sscanf(buffer, "%s %d ", buffer3, &IntegrationNet);
			IntegrationMethod=INTEGRATE_GAUSS;
		}
		if (!strcmp(buffer2, "IntegrateGaussKronrod"))
		{
			IntegrationMethod=INTEGRATE_GAUSSKRONROD;
//			sscanf(buffer, "%s %d ", buffer3, &IntegrationNet);
		}
		if (!strcmp(buffer2, "IntegrateVoitlander"))
		{
			IntegrationMethod=INTEGRATE_VOITLANDER;
		}
		// Normalise the calculated spectrum by its peak-to-peak amplitude (applied in Discrepancy)
		if (!strcmp(buffer2, "AmplitudeNormalized"))
		{
			AmplitudeNormalized=true;
//			sscanf(buffer, "%s %d ", buffer3, &IntegrationNet);
		}


		// Intensity cutoff handed to the rigid-limit spectrum routines (CalculateIndividualSpectrum)
		if (!strcmp(buffer2, "TransitionCutoff:"))
		{
			sscanf(buffer, "%s %lf ", buffer3, &TransitionIntensityCutoff);
		}
		// Optimiser settings: algorithm name, tolerances XTol and FTol, and initial step
		if (!strcmp(buffer2, "Algorithm:"))
		{
			sscanf(buffer, "%s %s ", buffer3, NLOptAlgorithm);
		}
		if (!strcmp(buffer2, "XTol:"))
		{
			sscanf(buffer, "%s %lf ", buffer3, &XTol);
		}
		if (!strcmp(buffer2, "FTol:"))
		{
			sscanf(buffer, "%s %lf ", buffer3, &FTol);
		}	
		if (!strcmp(buffer2, "InitialStep:"))
		{
			sscanf(buffer, "%s %lf ", buffer3, &InitialStep);
		}		

		// Number of the conformer group to which this particle belongs
		if (!strcmp(buffer2, "ConformerGroup:"))
		{
			sscanf(buffer, "%s %d ", buffer3, &ConformerGroup);
		}

		// Conversion rates between conformers: lines of i, j, value and flag (i<j) up to a blank line
		if (!strcmp(buffer2, "ConversionRate:"))
		{
			while(magnetic_par.getline(buffer3, MAX_STRLEN))
			{
				if((strlen(buffer3)==0) || (!strcmp(buffer3, "\r"))) break;
				int I, J, VF;
				double Value;
				sscanf(buffer3, "%d %d %lf %d", &I, &J, &Value, &VF);
				if(I<J) 
				{
					ConversionRate[I][J]=Value;
					ConversionRateVF[I][J]=VF;
				}
				else
				{
					printf("Warning: for conformer transition probability w_ij, it should be i<j. Value ignored!\n");
					printf("Line: %s\n", buffer3);
				}
			}
			sscanf(buffer, "%s", buffer2);

		}


		// The header of another particle ends the section of this particle
		if (!strcmp(buffer2, "Particle:"))
		{
			int npart;
			sscanf(buffer, "%s %d", buffer3, &npart);
			if(npart!=Count) break;
		}

		// The slow-motion (stochastic Liouville) calculation is used if any rotational diffusion rate is positive
		if((Rx>0)||(Ry>0)||(Rz>0))	Liouville=true;
		else						Liouville=false;


	}

}

// Converts all angles stored in this object from radians to degrees: Euler angles of the g, A, Gaussian, Lorentzian
// and diffusion tensors, libration amplitudes, frame angles and distribution widths, discrete-jump angles and the ZFS
// angle. The parameter files use degrees, the calculation itself uses radians.
// Used by: Discrepancy() in Discrepancy.cpp and main() in main.cpp, immediately before the parameters are written.
void MagneticParameters::ConvertAnglesToDegrees(void)
{
		// Euler angles of the g, A, line-width and diffusion tensors
		g1Euler[0]/=PI/180;
		g1Euler[1]/=PI/180;
		g1Euler[2]/=PI/180;

		g2Euler[0]/=PI/180;
		g2Euler[1]/=PI/180;
		g2Euler[2]/=PI/180;

		A1Euler[0]/=PI/180;
		A1Euler[1]/=PI/180;
		A1Euler[2]/=PI/180;

		A2Euler[0]/=PI/180;
		A2Euler[1]/=PI/180;
		A2Euler[2]/=PI/180;

		HGaussEuler[0]/=PI/180;
		HGaussEuler[1]/=PI/180;
		HGaussEuler[2]/=PI/180;

		HLorentzEuler[0]/=PI/180;
		HLorentzEuler[1]/=PI/180;
		HLorentzEuler[2]/=PI/180;

		DiffusionEuler[0]/=PI/180;
		DiffusionEuler[1]/=PI/180;
		DiffusionEuler[2]/=PI/180;

		// Libration amplitudes, Euler angles of the libration frame and distribution widths
		LibrationAmplitude[0]/=PI/180;
		LibrationAmplitude[1]/=PI/180;
		LibrationAmplitude[2]/=PI/180;

		LibrationEuler[0]/=PI/180;
		LibrationEuler[1]/=PI/180;
		LibrationEuler[2]/=PI/180;

		LibDistrWidth[0]/=PI/180;
		LibDistrWidth[1]/=PI/180;
		LibDistrWidth[2]/=PI/180;

		// Axis angles and jump angles of the discrete jumps
		for(int i=0; i<NJumpAxes; i++)
		{
			JumpAxisFi[i]/=PI/180;
			JumpAxisTeta[i]/=PI/180;
			JumpAngle[i]/=PI/180;
		}

		// Angle of the zero-field splitting tensor
		ZFS_Angle/=PI/180;

}

// Converts all angles stored in this object from degrees to radians (the inverse of ConvertAnglesToDegrees).
// Used by: main() in main.cpp (after reading the input files and again after writing the final output) and
// Discrepancy() in Discrepancy.cpp (after writing the temporary parameter files).
void MagneticParameters::ConvertAnglesToRadians(void)
{
		// Euler angles of the g, A, line-width and diffusion tensors
		g1Euler[0]*=PI/180;
		g1Euler[1]*=PI/180;
		g1Euler[2]*=PI/180;

		g2Euler[0]*=PI/180;
		g2Euler[1]*=PI/180;
		g2Euler[2]*=PI/180;

		A1Euler[0]*=PI/180;
		A1Euler[1]*=PI/180;
		A1Euler[2]*=PI/180;

		A2Euler[0]*=PI/180;
		A2Euler[1]*=PI/180;
		A2Euler[2]*=PI/180;

		HGaussEuler[0]*=PI/180;
		HGaussEuler[1]*=PI/180;
		HGaussEuler[2]*=PI/180;

		HLorentzEuler[0]*=PI/180;
		HLorentzEuler[1]*=PI/180;
		HLorentzEuler[2]*=PI/180;

		DiffusionEuler[0]*=PI/180;
		DiffusionEuler[1]*=PI/180;
		DiffusionEuler[2]*=PI/180;

		// Libration amplitudes, Euler angles of the libration frame and distribution widths
		LibrationAmplitude[0]*=PI/180;
		LibrationAmplitude[1]*=PI/180;
		LibrationAmplitude[2]*=PI/180;

		LibrationEuler[0]*=PI/180;
		LibrationEuler[1]*=PI/180;
		LibrationEuler[2]*=PI/180;

		LibDistrWidth[0]*=PI/180;
		LibDistrWidth[1]*=PI/180;
		LibDistrWidth[2]*=PI/180;

		// Axis angles and jump angles of the discrete jumps
		for(int i=0; i<NJumpAxes; i++)
		{
			JumpAxisFi[i]*=PI/180;
			JumpAxisTeta[i]*=PI/180;
			JumpAngle[i]*=PI/180;
		}
		
		// Angle of the zero-field splitting tensor
		ZFS_Angle*=PI/180;

}
// Appends the full parameter listing of this particle to a file, as Keyword: value lines in the style of the magnetic
// input file. Angles are written as currently stored (degrees when called after ConvertAnglesToDegrees). The file is opened in append mode, so the sections of successive particles follow
// each other. Input: filename, the file to append to.
// Used by: main() in main.cpp (final output, to the output magnetic file) and Discrepancy() in Discrepancy.cpp
// (to the file Magnetic.temp after a calculation).
void MagneticParameters::WriteToFile(char* filename)
{
		FILE *file;
		// Open in append mode and write the particle number, the monoradical flag and the conformer group
		file=fopen(filename, "a+");
		fprintf(file, "\nParticle: %d\n", Count);
		if (isMonoradical)
			fprintf(file, "Monoradical\n");
		fprintf(file, "ConformerGroup: %d\n", ConformerGroup);
		// Conformer conversion rates, only for more than one conformer
		if(NConformers>1)
		{
			fprintf(file, "ConversionRate:\n");
			for(int i=0; i<NConformers-1; i++)
				for(int j=i+1; j<NConformers; j++)
					if((ConversionRate[i][j]>=1.0e-10)||(ConversionRateVF[i][j]!=0))
						fprintf(file, " %d  %d  %10e    %d \n", i, j, ConversionRate[i][j], ConversionRateVF[i][j]);
			fprintf(file, "\n");
		}
		

		// Principal values, variate flags and Euler angles of the g and hyperfine tensors
		fprintf(file, "g1Tensor:         %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", g1(0, 0), g1VF(0, 0), g1(1, 1), g1VF(1, 1), g1(2, 2), g1VF(2, 2));
		fprintf(file, "g1Euler:          %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", g1Euler[0], g1EulerVF[0], g1Euler[1], g1EulerVF[1], g1Euler[2], g1EulerVF[2]);

		fprintf(file, "g2Tensor:         %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", g2(0, 0), g2VF(0, 0), g2(1, 1), g2VF(1, 1), g2(2, 2), g2VF(2, 2));
		fprintf(file, "g2Euler:          %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n\n", g2Euler[0], g2EulerVF[0], g2Euler[1], g2EulerVF[1], g2Euler[2], g2EulerVF[2]);

		fprintf(file, "A1Tensor:         %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", A1(0, 0), A1VF(0, 0), A1(1, 1), A1VF(1, 1), A1(2, 2), A1VF(2, 2));
		fprintf(file, "A1Euler:          %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", A1Euler[0], A1EulerVF[0], A1Euler[1], A1EulerVF[1], A1Euler[2], A1EulerVF[2]);
 
		fprintf(file, "A2Tensor:         %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", A2(0, 0), A2VF(0, 0), A2(1, 1), A2VF(1, 1), A2(2, 2), A2VF(2, 2));
		fprintf(file, "A2Euler:          %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n\n", A2Euler[0], A2EulerVF[0], A2Euler[1], A2EulerVF[1], A2Euler[2], A2EulerVF[2]);

		// Gaussian and Lorentzian line widths (FWHM in Gauss) with their Euler angles
		fprintf(file, "GaussTensor:      %10.6lf    %d    %10.6lf    %d    %10.6lf    %d     // FWHM in Gauss \n", HGauss(0, 0), HGaussVF(0, 0), HGauss(1, 1), HGaussVF(1, 1), HGauss(2, 2), HGaussVF(2, 2));
		fprintf(file, "GaussEuler:       %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", HGaussEuler[0], HGaussEulerVF[0], HGaussEuler[1], HGaussEulerVF[1], HGaussEuler[2], HGaussEulerVF[2]);

		fprintf(file, "LorentzTensor:    %10.6lf    %d    %10.6lf    %d    %10.6lf    %d     // FWHM in Gauss \n", HLorentz(0, 0), HLorentzVF(0, 0), HLorentz(1, 1), HLorentzVF(1, 1), HLorentz(2, 2), HLorentzVF(2, 2));
		fprintf(file, "LorentzEuler:     %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n\n", HLorentzEuler[0], HLorentzEulerVF[0], HLorentzEuler[1], HLorentzEulerVF[1], HLorentzEuler[2], HLorentzEulerVF[2]);

		// Dipolar vector, zero-field splitting and exchange parameters
		fprintf(file, "Dipolar_vector:   %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", dx, dxVF, dy, dyVF, dz, dzVF);
		fprintf(file, "ZFS:              %10.6lf    %d    %+10.6lf   %d    \n", ZFS_E, ZFS_EVF, ZFS_Angle, ZFS_AngleVF);
		fprintf(file, "Exchange:         %10.5lf    %d    %10.5lf    %d    %10e    %d    \n\n", J0, J0VF, J2, J2VF, JDiffusion, JDiffusionVF);

		// Quasi-libration parameters
		fprintf(file, "Libration:        %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", LibrationAmplitude[0], LibrationAmplitudeVF[0], LibrationAmplitude[1], LibrationAmplitudeVF[1], LibrationAmplitude[2], LibrationAmplitudeVF[2]);
		fprintf(file, "LibrationEuler:   %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", LibrationEuler[0], LibrationEulerVF[0], LibrationEuler[1], LibrationEulerVF[1], LibrationEuler[2], LibrationEulerVF[2]);
		fprintf(file, "LibDistrWidth:    %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", LibDistrWidth[0], LibDistrWidthVF[0], LibDistrWidth[1], LibDistrWidthVF[1], LibDistrWidth[2], LibDistrWidthVF[2]);
		fprintf(file, "LibDistrNSites:   %d    \n\n", LibDistrNSites);

		// Rotational diffusion rates, diffusion types and jump-diffusion residence times
		fprintf(file, "Diffusion:        %10e    %d    %10e    %d    %10e    %d    \n", Rx, RxVF, Ry, RyVF, Rz, RzVF);
		fprintf(file, "DiffusionEuler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", DiffusionEuler[0], DiffusionEulerVF[0], DiffusionEuler[1], DiffusionEulerVF[1], DiffusionEuler[2], DiffusionEulerVF[2]);

		fprintf(file, "DiffusionType:      %4.2lf    %4.2lf    %4.2lf        // Eperp, Eperp*, Epar \n",  DiffusionType[0], DiffusionType[1], DiffusionType[2]);
		fprintf(file, "RTau:             %10.6lf    %d    %10.6lf    %d    %10.6lf    %d          \n",  RTau[0], RTauVF[0],  RTau[1], RTauVF[1],  RTau[2], RTauVF[2]);

		// Discrete-jump axes
		for(int i=0; i<NJumpAxes; i++)
		{
			fprintf(file, "Jumps:            %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    %10.6lf    %d      // fi, teta, angle, frequency\n", JumpAxisFi[i], JumpAxisFiVF[i], JumpAxisTeta[i], JumpAxisTetaVF[i], JumpAngle[i], JumpAngleVF[i], JumpFrequency[i], JumpFrequencyVF[i]);
		}
		fprintf(file, "\n");

		// Mean-field potential coefficients, or Potential: no
		if(MFPotential.Lmax)
		{
			fprintf(file, "Potential:\n");
			for(int L=2; L<=MFPotential.Lmax; L+=2)
				for(int K=0; K<=L; K+=2)
					fprintf(file, "%d %d %10.5lf     %d \n", L, K, MFPotential.Get(L, K), MFPotential.GetVF(L, K));
			fprintf(file, "\n");

//			fprintf(file, "NetPotential:         %d   %d          // Integration net for teta, fi\n", MFPotential.INT_NET, MFPotential.INT_NET_FI);
		}

		else
			fprintf(file, "Potential:        no\n\n");


		// Isotropic hyperfine nuclei, or Hyperfine: no
		if(HyperfineCount)
		{
			fprintf(file, "Hyperfine:\n");
			for(int i=0; i<HyperfineCount; i++)
					fprintf(file, "%4s %2d %10.5lf     %d \n", HFS[i].AtomType, HFS[i].Degeneracy, HFS[i].HFConstant, HFS[i].HFConstantVF);
			fprintf(file, "\n");

//			fprintf(file, "NetPotential:         %d   %d          // Integration net for teta, fi\n", MFPotential.INT_NET, MFPotential.INT_NET_FI);
		}

		else
			fprintf(file, "Hyperfine:        no\n\n");


		
		// Basis truncation, pruning and Lanczos settings, then field net, amplitude and
		// orientation integration settings
		fprintf(file, "Truncation:       %d      %d      %d      %d      %d        // LE, LO, K, M, N \n\n", LEmax, LOmax, Kmax, Mmax, Nmax);
		fprintf(file, "Pruning:          %d      %10.5lf    %10.5lf    // NPrune dHPrune PruningTolerance\n\n", NPrune, dHPrune, PruningTolerance);
		fprintf(file, "LanczosSteps:     %d\n", LanczosSteps);
		fprintf(file, "FCenter:         %8.3lf            %10.6lf    %d   \n", CenterField0, Shift, ShiftVF);
		fprintf(file, "Amplitude:        %10.6lf    %d    \n", AmplitudeFactor, AmplitudeFactorVF);
		fprintf(file, "NPoints:             %d    \n", NPoints);
		fprintf(file, "DeltaH:           %lf  \n", fs);
/*
		if(IntegrateGaussKronrod)
			fprintf(file, "Integration_Net:      %d    \n", KronrodNodes);
		else
			fprintf(file, "Integration_Net:      %d    \n", 2*GaussNodes+1);
*/
		fprintf(file, "Integration_Net:      %d    \n", IntegrationNet);

		// Keyword of the orientation integration method
		if(IntegrationMethod==INTEGRATE_GAUSSKRONROD)
			fprintf(file, "IntegrateGaussKronrod\n");
		else if(IntegrationMethod==INTEGRATE_VOITLANDER)
			fprintf(file, "IntegrateVoitlander\n");
		else
			fprintf(file, "IntegrateGauss\n");

		if(AmplitudeNormalized)
			fprintf(file, "AmplitudeNormalized\n");
		
		fprintf(file, "TransitionCutoff: %lf  \n", TransitionIntensityCutoff);
		// The optimiser settings are stored once, in the section of particle 0
		if(Count==0)
		{
		fprintf(file, "\nAlgorithm:   %s  \n", NLOptAlgorithm);
		fprintf(file, "XTol: %e  \n", XTol);
		fprintf(file, "FTol: %e  \n", FTol);
		fprintf(file, "InitialStep: %e  \n", InitialStep);
		}

		fclose(file);
}


// Appends a compact, human-readable summary of the parameters of this particle to a file: one tensor with its Euler
// angles per line, leaving out items that do not apply (the second electron and the spin-spin terms of a monoradical,
// diffusion, jumps and potential of a rigid-limit calculation). The file is opened in append mode.
// Input: filename, the file to append to.
// Used by: WriteOutput (ReadExperimental.cpp), for every particle after the header of each spectrum file.
void MagneticParameters::WriteToFileBrief(char* filename)
{
		FILE *file;
		// Open in append mode and write the particle number and the conformer data
		file=fopen(filename, "a+");
		fprintf(file, "Particle: %d\n", Count);
		fprintf(file, "ConformerGroup: %d\n", ConformerGroup);
		if(NConformers>1)
		{
			fprintf(file, "ConversionRate:\n");
			for(int i=0; i<NConformers-1; i++)
				for(int j=i+1; j<NConformers; j++)
					if((ConversionRate[i][j]>=1.0e-10)||(ConversionRateVF[i][j]!=0))
						fprintf(file, "c %d  %d  %10e    %d \n", i, j, ConversionRate[i][j], ConversionRateVF[i][j]);
			fprintf(file, "\n");

		}
		// Monoradical flag, then the tensors with variate flags and Euler angles; the items of the second electron and the
		// spin-spin terms are written only for biradicals
		if (isMonoradical)
			fprintf(file, "Monoradical\n");

		fprintf(file, "g1Tensor:       %10.6lf    %d    %10.6lf    %d    %10.6lf    %d", g1(0, 0), g1VF(0, 0), g1(1, 1), g1VF(1, 1), g1(2, 2), g1VF(2, 2));
		fprintf(file, ", Euler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", g1Euler[0], g1EulerVF[0], g1Euler[1], g1EulerVF[1], g1Euler[2], g1EulerVF[2]);

		if(!isMonoradical) fprintf(file, "g2Tensor:       %10.6lf    %d    %10.6lf    %d    %10.6lf    %d", g2(0, 0), g2VF(0, 0), g2(1, 1), g2VF(1, 1), g2(2, 2), g2VF(2, 2));
		if(!isMonoradical) fprintf(file, ", Euler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", g2Euler[0], g2EulerVF[0], g2Euler[1], g2EulerVF[1], g2Euler[2], g2EulerVF[2]);

		fprintf(file, "A1Tensor:       %10.6lf    %d    %10.6lf    %d    %10.6lf    %d", A1(0, 0), A1VF(0, 0), A1(1, 1), A1VF(1, 1), A1(2, 2), A1VF(2, 2));
		fprintf(file, ", Euler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", A1Euler[0], A1EulerVF[0], A1Euler[1], A1EulerVF[1], A1Euler[2], A1EulerVF[2]);
 
		if(!isMonoradical) fprintf(file, "A2Tensor:       %10.6lf    %d    %10.6lf    %d    %10.6lf    %d", A2(0, 0), A2VF(0, 0), A2(1, 1), A2VF(1, 1), A2(2, 2), A2VF(2, 2));
		if(!isMonoradical) fprintf(file, ", Euler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", A2Euler[0], A2EulerVF[0], A2Euler[1], A2EulerVF[1], A2Euler[2], A2EulerVF[2]);

		fprintf(file, "GaussTensor:    %10.6lf    %d    %10.6lf    %d    %10.6lf    %d", HGauss(0, 0), HGaussVF(0, 0), HGauss(1, 1), HGaussVF(1, 1), HGauss(2, 2), HGaussVF(2, 2));
		fprintf(file, ", Euler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", HGaussEuler[0], HGaussEulerVF[0], HGaussEuler[1], HGaussEulerVF[1], HGaussEuler[2], HGaussEulerVF[2]);

		fprintf(file, "LorentzTensor:  %10.6lf    %d    %10.6lf    %d    %10.6lf    %d", HLorentz(0, 0), HLorentzVF(0, 0), HLorentz(1, 1), HLorentzVF(1, 1), HLorentz(2, 2), HLorentzVF(2, 2));
		fprintf(file, ", Euler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d    \n", HLorentzEuler[0], HLorentzEulerVF[0], HLorentzEuler[1], HLorentzEulerVF[1], HLorentzEuler[2], HLorentzEulerVF[2]);

		if(!isMonoradical) fprintf(file, "Dipolar_vector: %10.6lf    %d    %10.6lf    %d    %10.6lf    %d", dx, dxVF, dy, dyVF, dz, dzVF);
		if(!isMonoradical) fprintf(file, ", ZFS: %10.6lf    %d   %+10.6lf   %d", ZFS_E, ZFS_EVF, ZFS_Angle, ZFS_AngleVF);
		if(!isMonoradical) fprintf(file, ", Exchange:%10.5lf    %d    %10.5lf    %d    %10e    %d    \n\n", J0, J0VF, J2, J2VF, JDiffusion, JDiffusionVF);

	// Libration parameters, written only if a libration amplitude is non-zero
	if (LibrationAmplitude[0] || LibrationAmplitude[1] || LibrationAmplitude[0])
	{
		fprintf(file, "Libration:        %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", LibrationAmplitude[0], LibrationAmplitudeVF[0], LibrationAmplitude[1], LibrationAmplitudeVF[1], LibrationAmplitude[2], LibrationAmplitudeVF[2]);
		fprintf(file, "LibrationEuler:   %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", LibrationEuler[0], LibrationEulerVF[0], LibrationEuler[1], LibrationEulerVF[1], LibrationEuler[2], LibrationEulerVF[2]);
		if(LibDistrNSites) fprintf(file, "LibDistrWidth:    %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    \n", LibDistrWidth[0], LibDistrWidthVF[0], LibDistrWidth[1], LibDistrWidthVF[1], LibDistrWidth[2], LibDistrWidthVF[2]);
		if(LibDistrNSites) fprintf(file, "LibDistrNSites:   %d    \n", LibDistrNSites);
	}
    // Slow-motion parameters: diffusion, jump axes and potential, for the Liouville calculation only
    if(Liouville)
	{
		fprintf(file, "Diffusion:   %10e    %d    %10e    %d    %10e    %d ", Rx, RxVF, Ry, RyVF, Rz, RzVF);

		fprintf(file, ", Euler:   %+10.6lf    %d    %+10.6lf    %d    %+10.6lf    %d, ", DiffusionEuler[0], DiffusionEulerVF[0], DiffusionEuler[1], DiffusionEulerVF[1], DiffusionEuler[2], DiffusionEulerVF[2]);
		fprintf(file, "DiffusionType: %4.2lf  %4.2lf  %4.2lf , ",  DiffusionType[0], DiffusionType[1], DiffusionType[2]);
		fprintf(file, "RTau: %10.6lf %d    %10.6lf %d    %10.6lf %d  \n",  RTau[0], RTauVF[0],  RTau[1], RTauVF[1],  RTau[2], RTauVF[2]);
	
		for(int i=0; i<NJumpAxes; i++)
		{
			fprintf(file, "Jumps:   %10.6lf    %d    %10.6lf    %d    %10.6lf    %d    %10.6lf    %d      // fi, teta, angle, frequency\n", JumpAxisFi[i], JumpAxisFiVF[i], JumpAxisTeta[i], JumpAxisTetaVF[i], JumpAngle[i], JumpAngleVF[i], JumpFrequency[i], JumpFrequencyVF[i]);
		}

		if(MFPotential.Lmax)
		{
			fprintf(file, "Potential:\n");
			for(int L=2; L<=MFPotential.Lmax; L+=2)
			{
				for(int K=0; K<=L; K+=2)
					fprintf(file, "c %d %d %10.5lf     %d,\t", L, K, MFPotential.Get(L, K), MFPotential.GetVF(L, K));
				fprintf(file, "\n");
			}
//			fprintf(file, "NetPotential:         %d   %d          // Integration net for teta, fi\n", MFPotential.INT_NET, MFPotential.INT_NET_FI);
		}

	  }
		// Isotropic hyperfine nuclei
		if(HyperfineCount)
		{
			fprintf(file, "Hyperfine:\n");
			for(int i=0; i<HyperfineCount; i++)
					fprintf(file, "%4s %2d %10.5lf     %d \n", HFS[i].AtomType, HFS[i].Degeneracy, HFS[i].HFConstant, HFS[i].HFConstantVF);
			fprintf(file, "\n");

//			fprintf(file, "NetPotential:         %d   %d          // Integration net for teta, fi\n", MFPotential.INT_NET, MFPotential.INT_NET_FI);
		}

		
		// Truncation limits and Lanczos steps for the Liouville calculation, pruning only if it is switched on
		if(Liouville) fprintf(file, "Truncation:       %d      %d      %d      %d     %d         // LE, LO, K, M, N \n", LEmax, LOmax, Kmax, Mmax, Nmax);
		if(Liouville) fprintf(file, "LanczosSteps:     %d\n", LanczosSteps);
		if(NPrune) fprintf(file, "Pruning:          %d      %10.5lf    %10.5lf    // NPrune dHPrune PruningTolerance\n", NPrune, dHPrune, PruningTolerance);
		// The field net settings are written for particle 0 only
		if(Count==0)
		{
			fprintf(file, "FCenter:   %8.3lf   %10.6lf   %d, ", CenterField0, Shift, ShiftVF);
			fprintf(file, "NPoints: %d, ", NPoints);
			fprintf(file, "DeltaH: %lf \n", fs);
		}
		// Amplitude and orientation integration settings
		fprintf(file, "Amplitude:   %10.6lf    %d, ", AmplitudeFactor, AmplitudeFactorVF);
		fprintf(file, "Integration_Net:      %d, ", IntegrationNet);

		if(AmplitudeNormalized)
			fprintf(file, "AmplitudeNormalized\n");
		if(IntegrationMethod==INTEGRATE_GAUSSKRONROD)
			fprintf(file, "IntegrateGaussKronrod\n");
		else if(IntegrationMethod==INTEGRATE_VOITLANDER)
			fprintf(file, "IntegrateVoitlander\n");
		else
			fprintf(file, "IntegrateGauss\n");

		fprintf(file, "TransitionCutoff: %lf  \n", TransitionIntensityCutoff);
		// The optimiser settings are written for particle 0 only
		if(Count==0)
		{
		fprintf(file, "Algorithm:   %s, ", NLOptAlgorithm);
		fprintf(file, "XTol: %e, ", XTol);
		fprintf(file, "FTol: %e, ", FTol);
		fprintf(file, "InitialStep: %e  \n\n", InitialStep);
		}

		fclose(file);
}


// Prepares the interaction tensors for the current parameter values. For a biradical it builds the 3x3 spin-spin
// tensor DJ from the dipolar vector, the zero-field splitting and the isotropic exchange J0. It then rotates the
// tensors from their own principal-axes frames into the g1 frame and from there into the molecular frame, using the
// Euler angles and RotationMatrix, and rotates the diffusion tensor. The diffusion, jump, JDiffusion and R1 rates
// are converted from 1/s to Gauss units. Finally it builds the spherical-tensor form (SphericalTensor) of every
// interaction and stores it in Interactions[0..5] for the construction of the Liouvillian.
// Inputs: the parameters of this object (angles in radians) and CenterField. Outputs: g1R, g2R, A1R, A2R,
// HGaussR, HLorentzR, DJ, DiffusionR, the rates in Gauss units, Fm* and Interactions of this object.
// Used by: Discrepancy() in Discrepancy.cpp, once for each particle at every evaluation of the discrepancy.
void MagneticParameters::RotateTensors()
{
    Matrix<double, 3, 3> RotMat, RotMat_T;
// Creating DJ matrix for spin-spin dipole and exchange coupling:

	if(!(isMonoradical))
	{
		
		// Length of the dipolar vector, and the Euler angles of the ZFS axis system derived from its direction (plus ZFS_Angle)
		double dr=sqrt(dx*dx+dy*dy+dz*dz);
		Matrix<double, 3, 3> D, E;
		double ksi, tet;
		double DEuler[3];
		tet=acos(dz/dr);
		if(dx==0) ksi=PI/2;
		else     ksi=PI-atan(dy/dx);

		DEuler[0]=ZFS_Angle;
		DEuler[1]=tet;
		DEuler[2]=ksi;

		// Dipolar tensor D=18561*(dr^2*I-3*d*d^T)/dr^5 for the vector d=(dx,dy,dz)
		D << dr*dr-3*dx*dx, -3*dx*dy, -3*dx*dz,
		    -3*dx*dy, dr*dr-3*dy*dy, -3*dy*dz,
		    -3*dx*dz, -3*dy*dz, dr*dr-3*dz*dz;
		D=D*18561.0*pow(dr, -5);

		// Non-axial zero-field splitting tensor diag(-ZFS_E, ZFS_E, 0), rotated by the Euler angles DEuler below
		E << -ZFS_E, 0, 0,
			  0, ZFS_E, 0,
			  0, 0, 0;

/*
		E << 18561.0*pow(dr, -3), 0, 0,
		     0, 18561.0*pow(dr, -3), 0,
			 0, 0, -18561.0*pow(dr, -3)*2;
*/


		RotMat=RotationMatrix(DEuler);
		E=RotMat*E*RotMat.transpose();

/*
		cout << D << endl << endl;
		cout << E << endl << endl;
*/

		// Isotropic exchange tensor J0*I, and the total spin-spin tensor DJ=D+J+E
		Matrix<double, 3, 3> J;
		J << 1, 0, 0,
			 0, 1, 0,
			 0, 0, 1;
		J=J*J0;
		DJ=D+J+E;
	}
	else   // if isMonoradical
	{
		// A monoradical has no spin-spin tensor
		DJ=Matrix<double, 3, 3>::Zero(3,3);
	}

// Euler angles in the program are built in the following way:
//   Euler angles of all magnetic tensors, except g1 are convert magnetic tensors into g1 tensor.
//   g1Euler angles convert g1 frame into molecular orientation frame.
// Rotation of tensors into g1-frame
	RotMat=RotationMatrix(g2Euler);
	g2R=RotMat*g2*RotMat.transpose();

    RotMat=RotationMatrix(A1Euler);
	A1R=RotMat*A1*RotMat.transpose();

    RotMat=RotationMatrix(A2Euler);
	A2R=RotMat*A2*RotMat.transpose();

    RotMat=RotationMatrix(HGaussEuler);
	HGaussR=RotMat*HGauss*RotMat.transpose();

    RotMat=RotationMatrix(HLorentzEuler);
	HLorentzR=RotMat*HLorentz*RotMat.transpose();
// Rotation of tensors into molecular frame:
    RotMat=RotationMatrix(g1Euler);
	g1R=RotMat*g1*RotMat.transpose();
	g2R=RotMat*g2R*RotMat.transpose();
	A1R=RotMat*A1R*RotMat.transpose();
	A2R=RotMat*A2R*RotMat.transpose();
	DJ =RotMat*DJ*RotMat.transpose();
	HGaussR=RotMat*HGaussR*RotMat.transpose();
	HLorentzR=RotMat*HLorentzR*RotMat.transpose();

    // Diffusion tensor diag(Rx,Ry,Rz) rotated by DiffusionEuler
    RotMat=RotationMatrix(DiffusionEuler);

	Matrix<double, 3, 3> Diffusion, DiffR;
	Diffusion << Rx, 0, 0,
				  0, Ry, 0,
				  0, 0, Rz;

	DiffR=RotMat*Diffusion*RotMat.transpose();
	// Store the rotated diffusion tensor as a complex matrix
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			DiffusionR(i, j)=DiffR(i, j);

// also, diffusion tensor will be in Gauss units
	double normR=2.002322*927.31e-30/1.05443e-34; // g0*beta_e/h_bar;

//	Rx/=normR;
//	Ry/=normR;
//	Rz/=normR;

//	JumpFrequencyR=JumpFrequency/normR;
	// Jump frequencies in Gauss units
	for(int i=0; i<NJumpAxes; i++)
		JumpFrequencyR[i]=JumpFrequency[i]/normR;

	// Diffusion tensor and JDiffusion in Gauss units
	for(int j=0; j<3; j++)
		for(int k=0; k<3; k++)
		{
			DiffusionR(j, k)/=normR;
		}
	JDiffusionR=JDiffusion/normR;


// Calculation of Spherical tensors:
// the spherical tensors built by the previous call are released before reallocation
	delete Fmg1;
	delete Fmg2;
	delete FmA1;
	delete FmA2;
	delete FmD;
	delete FmW;
	delete FmR;

	Fmg1=new SphericalTensor(g1R);
	FmA1=new SphericalTensor(A1R);
	Fmg2=new SphericalTensor(g2R);
	FmA2=new SphericalTensor(A2R);
	FmD=new SphericalTensor(DJ);

// the SLE relaxation path turns the isotropic trace of this tensor directly into a diagonal
// decay rate G of the resolvent (CFVD, LiouvilleMatrixElementB.cpp), and an exponential
// autocorrelation exp(-G|t|) Fourier-transforms to a Lorentzian of FWHM=2G, so the tensor fed
// to the SLE path must be halved to keep HLorentz a FWHM, matching the rigid-limit convention;
// HLorentzR itself is left untouched because CalculateTransitions uses it directly as FWHM
	Matrix<double, 3, 3> HLorentzSLE=0.5*HLorentzR;
	FmW=new SphericalTensor(HLorentzSLE);
	FmR=new SphericalTensor(DiffusionR);

// Calculation of Interactions:
	// Object, corresponding to each interaction in the spherical representation, is constructed for each term in the Spin Hamiltonian

	// Field scale B0=CenterField/g_iso, with g_iso the isotropic g value (trace/3) of particle 0
	double B0=CenterField*3/(MagneticParameterArray[0].g1R.trace());

//	cout << B0 << endl;
	// Fill the six interaction records for g1, g2, A1, A2, DJ and the Lorentzian line width: the spherical tensor (Fdu)
	// and the spin-space matrices of the interaction for the field scale B0 (AuL)
	Interactions[0].Fdu=*Fmg1;
	Interactions[0].AuL.Set(Ag1L, B0);
	strcpy(Interactions[0].Fdu.Name, "g1");

//	Interactions[0].AuL.PrintSmallMatrix(0, 0);

	Interactions[1].Fdu=*Fmg2;
	Interactions[1].AuL.Set(Ag2L, B0);
	strcpy(Interactions[1].Fdu.Name, "g2");

	Interactions[2].Fdu=*FmA1;
	Interactions[2].AuL.Set(AA1L, B0);
	strcpy(Interactions[2].Fdu.Name, "A1");

	Interactions[3].Fdu=*FmA2;
	Interactions[3].AuL.Set(AA2L, B0);
	strcpy(Interactions[3].Fdu.Name, "A2");

	Interactions[4].Fdu=*FmD;
	Interactions[4].AuL.Set(ADL, B0);
	strcpy(Interactions[4].Fdu.Name, "DJ");

	Interactions[5].Fdu=*FmW;
	Interactions[5].AuL.Set(AWLB, B0);
	strcpy(Interactions[5].Fdu.Name, "Lorentz");

}

// For conformationally mobile biradicals, builds the matrix ConversionRateR of the exchange between conformers in
// Gauss units: off-diagonal elements are minus the conversion rates, each diagonal element is minus the sum of the
// other elements of its row weighted by sqrt(|AmplitudeFactor of the other conformer / own AmplitudeFactor|).
// It also sums the absolute amplitude factors of the members of the conformer group into AmplitudeFactorAllConformers.
// Inputs: ConversionRate, NConformers and the AmplitudeFactor of the particles (through MagneticParameterArray).
// Outputs: ConversionRateR and AmplitudeFactorAllConformers of this object.
// Used by: Discrepancy() in Discrepancy.cpp, for every particle after all particles have been updated.
void MagneticParameters::CalculateTranstionRateMatrix()
{
	// Factor g0*beta_e/h_bar that converts rates in 1/s into Gauss units
	double normR=2.002322*927.31e-30/1.05443e-34; // g0*beta_e/h_bar;
	// Rate matrix of the exchange between conformers: off-diagonal elements first, then the diagonal
	if(NConformers>1)
	{
		for(int i=0; i<NConformers; i++)
		{
			for(int j=i+1; j<NConformers; j++)
			{
				ConversionRateR[i*MAX_PARTICLES+j]=-ConversionRate[i][j]/normR;
				ConversionRateR[j*MAX_PARTICLES+i]=-ConversionRate[i][j]/normR;
			}
			// diagonal element calculation:
			double sum=0.0;
			for(int j=0; j<NConformers; j++)
			{
				if(j!=i) sum-=ConversionRateR[i*MAX_PARTICLES+j]*sqrt(fabs(MagneticParameterArray[j].AmplitudeFactor/MagneticParameterArray[i].AmplitudeFactor));
			}
			ConversionRateR[i*MAX_PARTICLES+i]=sum;
		}
	}
	// Sum of the absolute amplitude factors over the members of the conformer group
	AmplitudeFactorAllConformers=0;
	for (int i=0; i<NConformers; i++)
	{
		AmplitudeFactorAllConformers+=fabs(MagneticParameterArray[ConformerGroupMembers[i]].AmplitudeFactor);
	}
/*
	std::fstream discr;
	discr.open("discrepancy", ios::app);
	for(int i=0; i<NConformers; i++)
	{
		for(int j=0; j<NConformers; j++)
		{
			discr << ConversionRateR[i*MAX_PARTICLES+j] << " ";
		}
		discr << endl;
	}
	discr << endl;
	for(int i=0; i<NConformers; i++)
		discr << fabs(MagneticParameterArray[ConformerGroupMembers[i]].AmplitudeFactor);
	discr << endl;
	discr << AmplitudeFactorAllConformers << endl << endl;

	discr.close();
	for(int i=0; i<NConformers; i++)
	{
		for(int j=0; j<NConformers; j++)
		{
			printf("%10e  ", ConversionRateR[i*MAX_PARTICLES+j]);
		}
		printf("\n");
	}
	printf("Amplitude=%lf \n", AmplitudeFactorAllConformers);
*/
}


// Averages the interaction tensors g1R, g2R, A1R, A2R and DJ over small-angle librations of the molecule. The tensors
// are rotated into the libration frame (LibrationEuler) and averaged over a set of orientations that covers the
// libration amplitudes LibrationAmplitude[0..2]: a tilt with weight sin(teta) when the first two amplitudes are
// non-zero, and a twist by an angle psi from minus to plus the third amplitude. The averaged tensors are rotated back
// to the molecular frame and replace the stored ones; if all amplitudes are below LIB_MIN they are left unchanged.
// The sampling uses NS steps for each angle.
// Used by: Discrepancy() in Discrepancy.cpp, after RotateTensors, for every particle.
void MagneticParameters::QuasiLibrations()
{
    Matrix<double, 3, 3> g1Lib, g2Lib, A1Lib, A2Lib, DLib;

// now all tensors are in Molecular frame. We rotate them to Libration Frame
    Matrix<double, 3, 3> RotMat;
    RotMat=RotationMatrix(LibrationEuler);

	g1R=RotMat*g1R*RotMat.transpose();
	g2R=RotMat*g2R*RotMat.transpose();
	A1R=RotMat*A1R*RotMat.transpose();
	A2R=RotMat*A2R*RotMat.transpose();
	DJ =RotMat*DJ*RotMat.transpose();
// now all tensors are in Librations frame.
// cout << "Before:" << endl << g1R << endl << endl;

	// Libration amplitudes about the three axes of the libration frame (radians)
	double LX=LibrationAmplitude[0];
	double LY=LibrationAmplitude[1];
	double LZ=LibrationAmplitude[2];
    Matrix<double, 3,3> RLib, g1L, g2L, A1L, A2L, DJL, g1sum, g2sum, A1sum, A2sum, DJsum;
    double EulerLib[3];
	double sum=0;
   g1sum << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   A1sum << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   g2sum << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   A2sum << 0, 0, 0, 0, 0, 0, 0, 0, 0;
   DJsum << 0, 0, 0, 0, 0, 0, 0, 0, 0;
 //cout << "Amplitudes:" << LX << " " << LY << " " << LZ << endl;

   // The twist amplitude LZ is limited to pi
   if (LZ>=PI) LZ=PI;
   double fi, teta, psi;
 // Case 1: only the twist amplitude LZ is non-zero
 if((fabs(LX)<=LIB_MIN)&&(fabs(LY)<=LIB_MIN))		// ONLY PSI!
 {
	 if(fabs(LZ)<=LIB_MIN)
	 {
//		cout << LZ;
// no need to change tensors
		g1sum=g1sum+g1R;
		A1sum=A1sum+A1R;
		g2sum=g2sum+g2R;
		A2sum=A2sum+A2R;
		DJsum=DJsum+DJ;	 
	 }
	 else
	 {
		// Average over twist angles psi from -LZ to LZ in steps of LZ/NS
		fi=0; teta=0;
		for(psi=-fabs(LZ); psi<=fabs(LZ);psi+=fabs(LZ)/NS)		
		   {
				// Euler angles of this libration step; rotate every tensor and add it to the running sums
				EulerLib[0]=fi;
				EulerLib[1]=teta;
				EulerLib[2]=-fi+psi;

				RLib=RotationMatrix(EulerLib);
				g1L=RLib*g1R*RLib.transpose();
				A1L=RLib*A1R*RLib.transpose();
				g2L=RLib*g2R*RLib.transpose();
				A2L=RLib*A2R*RLib.transpose();
				DJL=RLib* DJ*RLib.transpose();
	
				g1sum=g1sum+g1L;
				A1sum=A1sum+A1L;
				g2sum=g2sum+g2L;
				A2sum=A2sum+A2L;
				DJsum=DJsum+DJL;
				sum+=1;

//				cout << g1L << endl << endl;
		   }
		   // Mean over the sampled twist angles
		   g1sum=g1sum/sum;
		   A1sum=A1sum/sum;
		   g2sum=g2sum/sum;
		   A2sum=A2sum/sum;
		   DJsum=DJsum/sum;

	 }
 }
 else	//  also Lx, Ly
 {
   // Case 2: tilt libration. Average over the azimuth fi and the tilt angle teta, weighting each tilt by sin(teta)
   for(fi=0; fi<2*PI; fi+=2*PI/NS)
   {
	   // Largest tilt angle at this azimuth: LY for cos(fi)^2=1 and LX for cos(fi)^2=0, at most pi/2
	   double cos2fi=cos(fi)*cos(fi);
//	   double sin2fi=sin(fi)*sin(fi);
	   double tetafinal;
	   if((LX>=PI/2)&&(LY>=PI/2)) tetafinal=PI/2;
	   else 
	   {
		   tetafinal=LX+(LY-LX)*cos2fi;        // too many trouble with cone: atan(sqrt(1/(sin2fi/b1+cos2fi/a1)));
		   if(tetafinal>PI/2) tetafinal=PI/2;				// check of maximum teta!
	   }
	   // Tilt angles from tetafinal/NS up to, but not including, tetafinal
	   for (teta=tetafinal/NS; teta<tetafinal; teta+=tetafinal/NS)
	   {
		   // Without a twist amplitude only the tilt is sampled, otherwise twist angles psi are added to every tilt
		   if(fabs(LZ)<=LIB_MIN)
		   {
				// Tilt only: rotate every tensor and add it with weight sin(teta)
				EulerLib[0]=fi;
				EulerLib[1]=teta;
				EulerLib[2]=-fi;

				RLib=RotationMatrix(EulerLib);
				g1L=RLib*g1R*RLib.transpose();
				A1L=RLib*A1R*RLib.transpose();
				g2L=RLib*g2R*RLib.transpose();
				A2L=RLib*A2R*RLib.transpose();
				DJL=RLib* DJ*RLib.transpose();
	
				g1sum=g1sum+g1L*sin(teta);
				A1sum=A1sum+A1L*sin(teta);
				g2sum=g2sum+g2L*sin(teta);
				A2sum=A2sum+A2L*sin(teta);
				DJsum=DJsum+DJL*sin(teta);
				sum+=sin(teta);				

		   }
		   else
		   for(psi=-fabs(LZ); psi<=fabs(LZ);psi+=fabs(LZ)/NS)		// NOT UNIFORM!!! TODO: SYmmetric relative to psi=0;
		   {
				// Tilt and twist: rotate every tensor and add it with weight sin(teta)
				EulerLib[0]=fi;
				EulerLib[1]=teta;
				EulerLib[2]=-fi+psi;

				RLib=RotationMatrix(EulerLib);
				g1L=RLib*g1R*RLib.transpose();
				A1L=RLib*A1R*RLib.transpose();
				g2L=RLib*g2R*RLib.transpose();
				A2L=RLib*A2R*RLib.transpose();
				DJL=RLib* DJ*RLib.transpose();
	
				g1sum=g1sum+g1L*sin(teta);
				A1sum=A1sum+A1L*sin(teta);
				g2sum=g2sum+g2L*sin(teta);
				A2sum=A2sum+A2L*sin(teta);
				DJsum=DJsum+DJL*sin(teta);
				sum+=sin(teta);			
		   }
	   }
   }
   // Weighted mean over all sampled orientations
   g1sum=g1sum/sum;
   A1sum=A1sum/sum;
   g2sum=g2sum/sum;
   A2sum=A2sum/sum;
   DJsum=DJsum/sum;

 }

 //cout << "After:" << endl << g1sum << endl << endl;
// converting tensors back to Molecular Frame
	g1R=RotMat.transpose()*g1sum*RotMat;
	g2R=RotMat.transpose()*g2sum*RotMat;
	A1R=RotMat.transpose()*A1sum*RotMat;
	A2R=RotMat.transpose()*A2sum*RotMat;
	DJ =RotMat.transpose()*DJsum*RotMat;
// Correction: after the back-rotation above the tensors are in the molecular frame again, not in the libration frame
// now all tensors are in Librations frame.

}

// Allocates FieldNet and Spectrum (NPoints values each) and fills FieldNet with NPoints equidistant fields of step
// fs centred on CenterField0: FieldNet[i]=CenterField0+(i-(NPoints-1)/2)*fs, where (NPoints-1)/2 is an integer division.
// Inputs: CenterField0, NPoints and fs of this object. Outputs: FieldNet and Spectrum of this object.
// Used by: main() in main.cpp, once for each particle after its input has been read.
void MagneticParameters::GenerateFieldNet()
{
// Generating FieldNet:
	FieldNet=(double*)malloc(NPoints*sizeof(double));
	Spectrum=(double*)malloc(NPoints*sizeof(double));

	// Equidistant field values with step fs around the centre field
	for(int i=0; i<NPoints; i++)
		FieldNet[i]=CenterField0+(i-(NPoints-1)/2)*fs;


}


	// Constructor: sets the pointers to the spherical tensors (Fmg1, Fmg2, FmA1, FmA2, FmD, FmW, FmR) to NULL, so that
	// RotateTensors can delete them safely at its first call. All other members are filled by ReadFromFile, Par2Arg
	// and RotateTensors.
	// Used by: implicitly, when main() creates the array new MagneticParameters[NParticles] (main.cpp).
	MagneticParameters::MagneticParameters()
	{

// the spherical tensors are first allocated by RotateTensors, which deletes whatever it finds here
		Fmg1=NULL; Fmg2=NULL; FmA1=NULL; FmA2=NULL; FmD=NULL; FmW=NULL; FmR=NULL;
	}

	// Destructor: does nothing, the arrays and spherical tensors allocated by other methods are not released here.
	// Used by: nothing in this code base (the MagneticParameters array created in main.cpp is never deleted).
	MagneticParameters::~MagneticParameters()
	{
	}

// Advances the index into the flat optimiser argument arrays, with a bound check against MAX_ARG
// Used by: MagneticParameters::Par2Arg (this file, after each parameter that becomes an argument) and
// OrderParameters::Par2Arg (ODF.cpp).
void NextArg(int *NArg)
{
	(*NArg)++;

// The next parameter would be written past the end of Arguments, lb, ub and cov
	if(*NArg>=MAX_ARG)
	{
		printf("Number of independently varied parameters reached the limit MAX_ARG=%d. Reduce the number of VariateFlags set to 1 or 2, or raise MAX_ARG in BiradicalODF.h and rebuild.\n", MAX_ARG);
		exit(0);
	}
}


		// Collects the parameters that the optimiser may vary into the flat arrays Arguments (starting values), lb and ub
		// (lower and upper bounds, from the MIN_* and MAX_* constants in BiradicalODF.h), advancing the index *NArg. A
		// parameter is collected if its variate flag VF is 1 or 2: 1 varies the value itself, 2 varies a transformed value
		// (1000*(g-2) for g values, log10 of the value for the exchange, diffusion, conversion and jump rates; if such a
		// rate is exactly zero the variation is switched to direct with a warning); for the other parameters 2 acts like 1.
		// The position of each collected parameter in Arguments is stored in the matching ...NArg member, -1 if the
		// parameter is not varied, so that Arg2Par can write the values back. Potential coefficients are collected only for
		// flag 1. The function also sets the nuclear spin and Sx2=2*spin of every hyperfine nucleus and stores the trace
		// sums used by variate flag 3. Inputs: the parameters of this object. Outputs: Arguments, ub, lb and *NArg are
		// extended, ...NArg members are set.
		// Used by: main() in main.cpp, once for each particle, before the optimisation starts.
		void MagneticParameters::Par2Arg(double* Arguments, double *ub, double *lb, int *NArg)
	{
	// Tensor components (g, A, Gaussian and Lorentzian widths): flag 1 or 2 stores the value and its bounds in Arguments
	/////////// g1 Tensor
		if ((g1VF(0,0)==1) || (g1VF(0,0)==2))
		{
			if (g1VF(0,0)==1) 
			{
				Arguments[*NArg]=g1(0,0);
				lb[*NArg]=MIN_GX;
				ub[*NArg]=MAX_GX;
			}
			else              
			{
				Arguments[*NArg]=1000.0*(g1(0,0)-2.0);
				lb[*NArg]=1000.0*(MIN_GX-2.0);
				ub[*NArg]=1000.0*(MAX_GX-2.0);
			}
			g1NArg(0,0)=*NArg;
			NextArg(NArg);
		}
		else g1NArg(0,0)=-1;
		if ((g1VF(1,1)==1) || (g1VF(1,1)==2))
		{
			if (g1VF(1,1)==1) 
			{
				Arguments[*NArg]=g1(1,1);
				lb[*NArg]=MIN_GY;
				ub[*NArg]=MAX_GY;
			}
			else              
			{
				Arguments[*NArg]=1000.0*(g1(1,1)-2.0);
				lb[*NArg]=1000.0*(MIN_GY-2.0);
				ub[*NArg]=1000.0*(MAX_GY-2.0);
			}
			g1NArg(1,1)=*NArg;
			NextArg(NArg);
		}		
		else g1NArg(1,1)=-1;
		if ((g1VF(2,2)==1) || (g1VF(2,2)==2))
		{
			if (g1VF(2,2)==1) 
			{
				Arguments[*NArg]=g1(2,2);
				lb[*NArg]=MIN_GZ;
				ub[*NArg]=MAX_GZ;
			}
			else              
			{
				Arguments[*NArg]=1000.0*(g1(2,2)-2.0);
				lb[*NArg]=1000.0*(MIN_GZ-2.0);
				ub[*NArg]=1000.0*(MAX_GZ-2.0);
			}
			g1NArg(2,2)=*NArg;
			NextArg(NArg);
		}
		else g1NArg(2,2)=-1;
	/////////// g2 Tensor
		if ((g2VF(0,0)==1) || (g2VF(0,0)==2))
		{
			if (g2VF(0,0)==1) 
			{
				Arguments[*NArg]=g2(0,0);
				lb[*NArg]=MIN_GX;
				ub[*NArg]=MAX_GX;
			}
			else              
			{
				Arguments[*NArg]=1000.0*(g2(0,0)-2.0);
				lb[*NArg]=1000.0*(MIN_GX-2.0);
				ub[*NArg]=1000.0*(MAX_GX-2.0);
			}
			g2NArg(0,0)=*NArg;
			NextArg(NArg);
		}
		else g2NArg(0,0)=-1;
		if ((g2VF(1,1)==1) || (g2VF(1,1)==2))
		{
			if (g2VF(1,1)==1) 
			{
				Arguments[*NArg]=g2(1,1);
				lb[*NArg]=MIN_GY;
				ub[*NArg]=MAX_GY;
			}
			else              
			{
				Arguments[*NArg]=1000.0*(g2(1,1)-2.0);
				lb[*NArg]=1000.0*(MIN_GY-2.0);
				ub[*NArg]=1000.0*(MAX_GY-2.0);
			}
			g2NArg(1,1)=*NArg;
			NextArg(NArg);
		}		
		else g2NArg(1,1)=-1;
		if ((g2VF(2,2)==1) || (g2VF(2,2)==2))
		{
			if (g2VF(2,2)==1) 
			{
				Arguments[*NArg]=g2(2,2);
				lb[*NArg]=MIN_GZ;
				ub[*NArg]=MAX_GZ;
			}
			else             
			{
				Arguments[*NArg]=1000.0*(g2(2,2)-2.0);
				lb[*NArg]=1000.0*(MIN_GZ-2.0);
				ub[*NArg]=1000.0*(MAX_GZ-2.0);
			}
			g2NArg(2,2)=*NArg;
			NextArg(NArg);
		}	
		else g2NArg(2,2)=-1;
	/////////// A1 Tensor
		if ((A1VF(0,0)==1) || (A1VF(0,0)==2))
		{
			Arguments[*NArg]=A1(0,0);
			lb[*NArg]=MIN_AX;
			ub[*NArg]=MAX_AX;
			A1NArg(0,0)=*NArg;
			NextArg(NArg);
		}
		else A1NArg(0,0)=-1;
		if ((A1VF(1,1)==1) || (A1VF(1,1)==2))
		{
			Arguments[*NArg]=A1(1,1);
			lb[*NArg]=MIN_AY;
			ub[*NArg]=MAX_AY;
			A1NArg(1,1)=*NArg;
			NextArg(NArg);
		}		
		else A1NArg(1,1)=-1;
		if ((A1VF(2,2)==1) || (A1VF(2,2)==2))
		{
			Arguments[*NArg]=A1(2,2);
			lb[*NArg]=MIN_AZ;
			ub[*NArg]=MAX_AZ;
			A1NArg(2,2)=*NArg;
			NextArg(NArg);
		}
		else A1NArg(2,2)=-1;
	/////////// A2 Tensor
		if ((A2VF(0,0)==1) || (A2VF(0,0)==2))
		{
			Arguments[*NArg]=A2(0,0);
			lb[*NArg]=MIN_AX;
			ub[*NArg]=MAX_AX;
			A2NArg(0,0)=*NArg;
			NextArg(NArg);
		}
		else A2NArg(0,0)=-1;
		if ((A2VF(1,1)==1) || (A2VF(1,1)==2))
		{
			Arguments[*NArg]=A2(1,1);
			lb[*NArg]=MIN_AY;
			ub[*NArg]=MAX_AY;
			A2NArg(1,1)=*NArg;
			NextArg(NArg);
		}		
		else A2NArg(1,1)=-1;
		if ((A2VF(2,2)==1) || (A2VF(2,2)==2))
		{
			Arguments[*NArg]=A2(2,2);
			lb[*NArg]=MIN_AZ;
			ub[*NArg]=MAX_AZ;
			A2NArg(2,2)=*NArg;
			NextArg(NArg);
		}
		else A2NArg(2,2)=-1;
	/////////// HGauss Tensor
		if ((HGaussVF(0,0)==1) || (HGaussVF(0,0)==2))
		{
			Arguments[*NArg]=HGauss(0,0);
			lb[*NArg]=MIN_LW;
			ub[*NArg]=MAX_LW;
			HGaussNArg(0,0)=*NArg;
			NextArg(NArg);
		}
		else HGaussNArg(0,0)=-1;
		if ((HGaussVF(1,1)==1) || (HGaussVF(1,1)==2))
		{
			Arguments[*NArg]=HGauss(1,1);
			lb[*NArg]=MIN_LW;
			ub[*NArg]=MAX_LW;
			HGaussNArg(1,1)=*NArg;
			NextArg(NArg);
		}		
		else HGaussNArg(1,1)=-1;
		if ((HGaussVF(2,2)==1) || (HGaussVF(2,2)==2))
		{
			Arguments[*NArg]=HGauss(2,2);
			lb[*NArg]=MIN_LW;
			ub[*NArg]=MAX_LW;
			HGaussNArg(2,2)=*NArg;
			NextArg(NArg);
		}
		else HGaussNArg(2,2)=-1;
	/////////// HLorentz Tensor
		if ((HLorentzVF(0,0)==1) || (HLorentzVF(0,0)==2))
		{
			lb[*NArg]=MIN_LW;
			ub[*NArg]=MAX_LW;
			Arguments[*NArg]=HLorentz(0,0);
			HLorentzNArg(0,0)=*NArg;
			NextArg(NArg);
		}
		else HLorentzNArg(0,0)=-1;
		if ((HLorentzVF(1,1)==1) || (HLorentzVF(1,1)==2))
		{
			lb[*NArg]=MIN_LW;
			ub[*NArg]=MAX_LW;
			Arguments[*NArg]=HLorentz(1,1);
			HLorentzNArg(1,1)=*NArg;
			NextArg(NArg);
		}		
		else HLorentzNArg(1,1)=-1;
		if ((HLorentzVF(2,2)==1) || (HLorentzVF(2,2)==2))
		{
			lb[*NArg]=MIN_LW;
			ub[*NArg]=MAX_LW;
			Arguments[*NArg]=HLorentz(2,2);
			HLorentzNArg(2,2)=*NArg;
			NextArg(NArg);
		}
		else HLorentzNArg(2,2)=-1;
// Euler angles of the tensor axes, diffusion axes and libration frame, then libration amplitudes and widths
///      EULER ANGLES
	/////////// g1Euler Euler
		if ((g1EulerVF[0]==1) || (g1EulerVF[0]==2))
		{
			Arguments[*NArg]=g1Euler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			g1EulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else g1EulerNArg[0]=-1;
		if ((g1EulerVF[1]==1) || (g1EulerVF[1]==2))
		{
			Arguments[*NArg]=g1Euler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			g1EulerNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else g1EulerNArg[1]=-1;
		if ((g1EulerVF[2]==1) || (g1EulerVF[2]==2))
		{
			Arguments[*NArg]=g1Euler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			g1EulerNArg[2]=*NArg;
			NextArg(NArg);
		}
		else g1EulerNArg[2]=-1;
	/////////// g2Euler Euler
		if ((g2EulerVF[0]==1) || (g2EulerVF[0]==2))
		{
			Arguments[*NArg]=g2Euler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			g2EulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else g2EulerNArg[0]=-1;
		if ((g2EulerVF[1]==1) || (g2EulerVF[1]==2))
		{
			Arguments[*NArg]=g2Euler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			g2EulerNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else g2EulerNArg[1]=-1;
		if ((g2EulerVF[2]==1) || (g2EulerVF[2]==2))
		{
			Arguments[*NArg]=g2Euler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			g2EulerNArg[2]=*NArg;
			NextArg(NArg);
		}	
		else g2EulerNArg[2]=-1;
	/////////// A1Euler Euler
		if ((A1EulerVF[0]==1) || (A1EulerVF[0]==2))
		{
			Arguments[*NArg]=A1Euler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			A1EulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else A1EulerNArg[0]=-1;
		if ((A1EulerVF[1]==1) || (A1EulerVF[1]==2))
		{
			Arguments[*NArg]=A1Euler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			A1EulerNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else A1EulerNArg[1]=-1;
		if ((A1EulerVF[2]==1) || (A1EulerVF[2]==2))
		{
			Arguments[*NArg]=A1Euler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			A1EulerNArg[2]=*NArg;
			NextArg(NArg);
		}
		else A1EulerNArg[2]=-1;
	/////////// A2Euler Euler
		if ((A2EulerVF[0]==1) || (A2EulerVF[0]==2))
		{
			Arguments[*NArg]=A2Euler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			A2EulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else A2EulerNArg[0]=-1;
		if ((A2EulerVF[1]==1) || (A2EulerVF[1]==2))
		{
			Arguments[*NArg]=A2Euler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			A2EulerNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else A2EulerNArg[1]=-1;
		if ((A2EulerVF[2]==1) || (A2EulerVF[2]==2))
		{
			Arguments[*NArg]=A2Euler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			A2EulerNArg[2]=*NArg;
			NextArg(NArg);
		}
		else A2EulerNArg[2]=-1;
	/////////// HGaussEuler Euler
		if ((HGaussEulerVF[0]==1) || (HGaussEulerVF[0]==2))
		{
			Arguments[*NArg]=HGaussEuler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			HGaussEulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else HGaussEulerNArg[0]=-1;
		if ((HGaussEulerVF[1]==1) || (HGaussEulerVF[1]==2))
		{
			Arguments[*NArg]=HGaussEuler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			HGaussEulerNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else HGaussEulerNArg[1]=-1;
		if ((HGaussEulerVF[2]==1) || (HGaussEulerVF[2]==2))
		{
			Arguments[*NArg]=HGaussEuler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			HGaussEulerNArg[2]=*NArg;
			NextArg(NArg);
		}
		else HGaussEulerNArg[2]=-1;
	/////////// HLorentzEuler Euler
		if ((HLorentzEulerVF[0]==1) || (HLorentzEulerVF[0]==2))
		{
			Arguments[*NArg]=HLorentzEuler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			HLorentzEulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else HLorentzEulerNArg[0]=-1;
		if ((HLorentzEulerVF[1]==1) || (HLorentzEulerVF[1]==2))
		{
			Arguments[*NArg]=HLorentzEuler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			HLorentzEulerNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else HLorentzEulerNArg[1]=-1;
		if ((HLorentzEulerVF[2]==1) || (HLorentzEulerVF[2]==2))
		{
			Arguments[*NArg]=HLorentzEuler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			HLorentzEulerNArg[2]=*NArg;
			NextArg(NArg);
		}
		else HLorentzEulerNArg[2]=-1;
	///// DiffusionEuler
		if ((DiffusionEulerVF[0]==1) || (DiffusionEulerVF[0]==2))
		{
			Arguments[*NArg]=DiffusionEuler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			DiffusionEulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else DiffusionEulerNArg[0]=-1;
		if ((DiffusionEulerVF[1]==1) || (DiffusionEulerVF[1]==2))
		{
			Arguments[*NArg]=DiffusionEuler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			DiffusionEulerNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else DiffusionEulerNArg[1]=-1;
		if ((DiffusionEulerVF[2]==1) || (DiffusionEulerVF[2]==2))
		{
			Arguments[*NArg]=DiffusionEuler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			DiffusionEulerNArg[2]=*NArg;
			NextArg(NArg);
		}
		else DiffusionEulerNArg[2]=-1;
	///// LibrationEuler
		if ((LibrationEulerVF[0]==1) || (LibrationEulerVF[0]==2))
		{
			Arguments[*NArg]=LibrationEuler[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibrationEulerNArg[0]=*NArg;
			NextArg(NArg);
		}
		else LibrationEulerNArg[0]=-1;
		if ((LibrationEulerVF[1]==1) || (LibrationEulerVF[1]==2))
		{
			Arguments[*NArg]=LibrationEuler[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibrationEulerNArg[1]=*NArg;
			NextArg(NArg);
		}
		else LibrationEulerNArg[1]=-1;
		if ((LibrationEulerVF[2]==1) || (LibrationEulerVF[2]==2))
		{
			Arguments[*NArg]=LibrationEuler[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibrationEulerNArg[2]=*NArg;
			NextArg(NArg);
		}
		else LibrationEulerNArg[2]=-1;

	///// LibrationAmplitude
		if ((LibrationAmplitudeVF[0]==1) || (LibrationAmplitudeVF[0]==2))
		{
			Arguments[*NArg]=LibrationAmplitude[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibrationAmplitudeNArg[0]=*NArg;
			NextArg(NArg);
		}
		else LibrationAmplitudeNArg[0]=-1;
		if ((LibrationAmplitudeVF[1]==1) || (LibrationAmplitudeVF[1]==2))
		{
			Arguments[*NArg]=LibrationAmplitude[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibrationAmplitudeNArg[1]=*NArg;
			NextArg(NArg);
		}
		else LibrationAmplitudeNArg[1]=-1;
		if ((LibrationAmplitudeVF[2]==1) || (LibrationAmplitudeVF[2]==2))
		{
			Arguments[*NArg]=LibrationAmplitude[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibrationAmplitudeNArg[2]=*NArg;
			NextArg(NArg);
		}
		else LibrationAmplitudeNArg[2]=-1;
	///// LibDistrWidth
		if ((LibDistrWidthVF[0]==1) || (LibDistrWidthVF[0]==2))
		{
			Arguments[*NArg]=LibDistrWidth[0];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibDistrWidthNArg[0]=*NArg;
			NextArg(NArg);
		}
		else LibDistrWidthNArg[0]=-1;
		if ((LibDistrWidthVF[1]==1) || (LibDistrWidthVF[1]==2))
		{
			Arguments[*NArg]=LibDistrWidth[1];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibDistrWidthNArg[1]=*NArg;
			NextArg(NArg);
		}
		else LibDistrWidthNArg[1]=-1;
		if ((LibDistrWidthVF[2]==1) || (LibDistrWidthVF[2]==2))
		{
			Arguments[*NArg]=LibDistrWidth[2];
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			LibDistrWidthNArg[2]=*NArg;
			NextArg(NArg);
		}
		else LibDistrWidthNArg[2]=-1;

// Non-tensors:
		// Components of the dipolar vector
		if((dxVF==1)||(dxVF==2))
		{
			Arguments[*NArg]=dx;
			lb[*NArg]=MIN_DX;
			ub[*NArg]=MAX_DX;
			dxNArg=*NArg;
			NextArg(NArg);
		}
		else dxNArg=-1;
		if((dyVF==1)||(dyVF==2))
		{
			Arguments[*NArg]=dy;
			lb[*NArg]=MIN_DY;
			ub[*NArg]=MAX_DY;
			dyNArg=*NArg;
			NextArg(NArg);
		}
		else dyNArg=-1;
		if((dzVF==1)||(dzVF==2))
		{
			Arguments[*NArg]=dz;
			lb[*NArg]=MIN_DZ;
			ub[*NArg]=MAX_DZ;
			dzNArg=*NArg;
			NextArg(NArg);
		}
		else dzNArg=-1;
		// Exchange parameters: flag 1 varies the value, flag 2 its log10
		if((J0VF==1)||(J0VF==2))
		{
			if((J0VF==2)&&(J0==0.0))
			{
				printf("Warning: J0 is varied as log10 but its starting value is exactly 0, switching to direct variation\n");
				J0VF=1;
			}
			if (J0VF==1) 
			{
				Arguments[*NArg]=J0;
				lb[*NArg]=MIN_EXCHANGE;
				ub[*NArg]=MAX_EXCHANGE;
			}
			else
			{
				Arguments[*NArg]=log10(J0);
				lb[*NArg]=-4;
				ub[*NArg]=log10((double)MAX_EXCHANGE);
			}

			J0NArg=*NArg;
			NextArg(NArg);
		}
		else J0NArg=-1;

		if((J2VF==1)||(J2VF==2))
		{
			if((J2VF==2)&&(J2==0.0))
			{
				printf("Warning: J2 is varied as log10 but its starting value is exactly 0, switching to direct variation\n");
				J2VF=1;
			}
			if (J2VF==1) 
			{
				Arguments[*NArg]=J2;
				lb[*NArg]=MIN_EXCHANGE;
				ub[*NArg]=MAX_EXCHANGE;
			}
			else
			{
				Arguments[*NArg]=log10(J2);
				lb[*NArg]=-4;
				ub[*NArg]=log10((double)MAX_EXCHANGE);
			}

			J2NArg=*NArg;
			NextArg(NArg);
		}
		else J2NArg=-1;

		if((JDiffusionVF==1)||(JDiffusionVF==2))
		{
			if((JDiffusionVF==2)&&(JDiffusion==0.0))
			{
				printf("Warning: JDiffusion is varied as log10 but its starting value is exactly 0, switching to direct variation\n");
				JDiffusionVF=1;
			}
			if (JDiffusionVF==1) 
			{
				Arguments[*NArg]=JDiffusion;
				lb[*NArg]=MIN_DIFFUSION;
				ub[*NArg]=MAX_DIFFUSION;
			}
			else
			{
				Arguments[*NArg]=log10(JDiffusion);
				lb[*NArg]=-4;
				ub[*NArg]=log10((double)MAX_DIFFUSION);
			}

			JDiffusionNArg=*NArg;
			NextArg(NArg);
		}
		else JDiffusionNArg=-1;


		// Zero-field splitting parameter E and angle
		if((ZFS_EVF==1)||(ZFS_EVF==2))
		{
			Arguments[*NArg]=ZFS_E;
			lb[*NArg]=-10000;
			ub[*NArg]=10000;
			ZFS_ENArg=*NArg;
			NextArg(NArg);
		}
		else ZFS_ENArg=-1;
		if((ZFS_AngleVF==1)||(ZFS_AngleVF==2))
		{
			Arguments[*NArg]=ZFS_Angle;
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			ZFS_AngleNArg=*NArg;
			NextArg(NArg);
		}
		else ZFS_AngleNArg=-1;


		// Rotational diffusion rates Rx, Ry, Rz: the value for flag 1, its log10 for flag 2
		if((RxVF==1)||(RxVF==2))
		{
			if((RxVF==2)&&(Rx==0.0))
			{
				printf("Warning: Rx is varied as log10 but its starting value is exactly 0, switching to direct variation\n");
				RxVF=1;
			}
			if (RxVF==1) 
			{
				Arguments[*NArg]=Rx;
				lb[*NArg]=MIN_DIFFUSION;
				ub[*NArg]=MAX_DIFFUSION;
			}
			else
			{
				Arguments[*NArg]=log10(Rx);
				lb[*NArg]=-4;
				ub[*NArg]=log10((double)MAX_DIFFUSION);
			}

			RxNArg=*NArg;
			NextArg(NArg);
		}
		else RxNArg=-1;
		if((RyVF==1)||(RyVF==2))
		{
			if((RyVF==2)&&(Ry==0.0))
			{
				printf("Warning: Ry is varied as log10 but its starting value is exactly 0, switching to direct variation\n");
				RyVF=1;
			}
			if (RyVF==1) 
			{
				Arguments[*NArg]=Ry;
				lb[*NArg]=MIN_DIFFUSION;
				ub[*NArg]=MAX_DIFFUSION;
			}
			else
			{
				Arguments[*NArg]=log10(Ry);
				lb[*NArg]=-4;
				ub[*NArg]=log10((double)MAX_DIFFUSION);
			}

			RyNArg=*NArg;
			NextArg(NArg);
		}
		else RyNArg=-1;
		if((RzVF==1)||(RzVF==2))
		{
			if((RzVF==2)&&(Rz==0.0))
			{
				printf("Warning: Rz is varied as log10 but its starting value is exactly 0, switching to direct variation\n");
				RzVF=1;
			}
			if (RzVF==1) 
			{
				Arguments[*NArg]=Rz;
				lb[*NArg]=MIN_DIFFUSION;
				ub[*NArg]=MAX_DIFFUSION;
			}
			else
			{
				Arguments[*NArg]=log10(Rz);
				lb[*NArg]=-4;
				ub[*NArg]=log10((double)MAX_DIFFUSION);
			}

			RzNArg=*NArg;
			NextArg(NArg);
		}
		else RzNArg=-1;
// Jump diffusion residence times RTau:
		if ((RTauVF[0]==1) || (RTauVF[0]==2))
		{
			Arguments[*NArg]=RTau[0];
			lb[*NArg]=0;
			ub[*NArg]=10000;
			RTauNArg[0]=*NArg;
			NextArg(NArg);
		}
		else RTauNArg[0]=-1;
		if ((RTauVF[1]==1) || (RTauVF[1]==2))
		{
			Arguments[*NArg]=RTau[1];
			lb[*NArg]=0;
			ub[*NArg]=10000;
			RTauNArg[1]=*NArg;
			NextArg(NArg);
		}		
		else RTauNArg[1]=-1;
		if ((RTauVF[2]==1) || (RTauVF[2]==2))
		{
			Arguments[*NArg]=RTau[2];
			lb[*NArg]=0;
			ub[*NArg]=10000;
			RTauNArg[2]=*NArg;
			NextArg(NArg);
		}
		else RTauNArg[2]=-1;
// Conformation analysis:
	if(NConformers>1)
	{
		for(int i=0; i<NConformers; i++)
			for(int j=i+1; j<NConformers; j++)
				if((ConversionRateVF[i][j]==1)||(ConversionRateVF[i][j]==2))
				{
					if((ConversionRateVF[i][j]==2)&&(ConversionRate[i][j]==0.0))
					{
						printf("Warning: ConversionRate[%d][%d] is varied as log10 but its starting value is exactly 0, switching to direct variation\n", i, j);
						ConversionRateVF[i][j]=1;
					}
					if (ConversionRateVF[i][j]==1) 
					{
						Arguments[*NArg]=ConversionRate[i][j];
						lb[*NArg]=MIN_DIFFUSION;
						ub[*NArg]=MAX_DIFFUSION;
					}
					else
					{
						Arguments[*NArg]=log10(ConversionRate[i][j]);
						lb[*NArg]=-4;
						ub[*NArg]=log10((double)MAX_DIFFUSION);
					}

					ConversionRateNArg[i][j]=*NArg;
					NextArg(NArg);
				}
				else ConversionRateNArg[i][j]=-1;
	
	}
// Hyperfine:
		for (int i=0; i<HyperfineCount; i++)
		{
			if((HFS[i].HFConstantVF==1)||(HFS[i].HFConstantVF==2))
			{
				Arguments[*NArg]=HFS[i].HFConstant;
				lb[*NArg]=0.0;
				ub[*NArg]=1.0e11;
				HFS[i].HFConstantNArg=*NArg;
				NextArg(NArg);
			}
			else HFS[i].HFConstantNArg=-1;
		}

		// Nuclear spin of each hyperfine nucleus from its name, and Sx2=2*spin
		for (int i=0; i<HyperfineCount; i++)
		{
			HFS[i].SetMagneticFromName();
			HFS[i].Sx2=(int)(2*HFS[i].NuclearSpin);
		}

// DiscreteJumps:
	for(int i=0; i<NJumpAxes; i++)
	{
		if((JumpFrequencyVF[i]==1)||(JumpFrequencyVF[i]==2))
		{
			if((JumpFrequencyVF[i]==2)&&(JumpFrequency[i]==0.0))
			{
				printf("Warning: JumpFrequency[%d] is varied as log10 but its starting value is exactly 0, switching to direct variation\n", i);
				JumpFrequencyVF[i]=1;
			}
			if (JumpFrequencyVF[i]==1) 
			{
				Arguments[*NArg]=JumpFrequency[i];
				lb[*NArg]=MIN_DIFFUSION;
				ub[*NArg]=MAX_DIFFUSION;
			}
			else
			{
				Arguments[*NArg]=log10(JumpFrequency[i]);
				lb[*NArg]=-4;
				ub[*NArg]=log10((double)MAX_DIFFUSION);
			}

			JumpFrequencyNArg[i]=*NArg;
			NextArg(NArg);
		}
		else JumpFrequencyNArg[i]=-1;

		if((JumpAxisFiVF[i]==1)||(JumpAxisFiVF[i]==2))
		{
				Arguments[*NArg]=JumpAxisFi[i];
				lb[*NArg]=MIN_ANGLE;
				ub[*NArg]=MAX_ANGLE;

			JumpAxisFiNArg[i]=*NArg;
			NextArg(NArg);
		}
		else JumpAxisFiNArg[i]=-1;

		if((JumpAxisTetaVF[i]==1)||(JumpAxisTetaVF[i]==2))
		{
				Arguments[*NArg]=JumpAxisTeta[i];
				lb[*NArg]=MIN_ANGLE;
				ub[*NArg]=MAX_ANGLE;

			JumpAxisTetaNArg[i]=*NArg;
			NextArg(NArg);
		}
		else JumpAxisTetaNArg[i]=-1;

		if((JumpAngleVF[i]==1)||(JumpAngleVF[i]==2))
		{
				Arguments[*NArg]=JumpAngle[i];
				lb[*NArg]=MIN_ANGLE;
				ub[*NArg]=MAX_ANGLE;

			JumpAngleNArg[i]=*NArg;
			NextArg(NArg);
		}
		else JumpAngleNArg[i]=-1;

	}




		// Mean-field potential coefficients: only flag 1 makes a coefficient an argument, the others are marked -1
		for(int L=2; L<=MFPotential.Lmax; L+=2)
			for(int K=0; K<=L; K+=2)
			{
				if(MFPotential.GetVF(L, K)==1)
				{
					Arguments[*NArg]=MFPotential.Get(L, K);
					lb[*NArg]=MIN_POTENTIAL;
					ub[*NArg]=MAX_POTENTIAL;

					MFPotential.SetNArg(L, K, *NArg);
					NextArg(NArg);
				}
				else
					MFPotential.SetNArg(L, K, -1);

			}

			


		// Field shift and amplitude factor
		if((ShiftVF==1)||(ShiftVF==2))
		{
			Arguments[*NArg]=Shift;
			lb[*NArg]=MIN_SDVIG;
			ub[*NArg]=MAX_SDVIG;
			ShiftNArg=*NArg;
			NextArg(NArg);
		}
		else ShiftNArg=-1;
		if((AmplitudeFactorVF==1)||(AmplitudeFactorVF==2))
		{
			Arguments[*NArg]=AmplitudeFactor;
			lb[*NArg]=MIN_CONT;
			ub[*NArg]=MAX_CONT;
			AmplitudeFactorNArg=*NArg;
			NextArg(NArg);
		}
		else AmplitudeFactorNArg=-1;
		// Sums of the three principal values, which variate flag 3 keeps constant in Arg2Par
		g1Average3=g1(0,0)+g1(1,1)+g1(2,2);
		g2Average3=g2(0,0)+g2(1,1)+g2(2,2);
		A1Average3=A1(0,0)+A1(1,1)+A1(2,2);
		A2Average3=A2(0,0)+A2(1,1)+A2(2,2);

	}

	// The inverse of Par2Arg: writes the values of the optimiser arguments back into the parameters, then applies the
	// dependent variate flags. Flag 2 values are converted back (g=2+0.001*argument, 10^argument for the exchange,
	// diffusion, conversion and jump rates). Flag 3 sets the z component of g1, g2, A1 or A2 so that the trace is kept;
	// flag 4 copies the y component to the x and z components (also for Rx, Rz and RTau); flag 5 copies the x component
	// to y and z; flag 6 sets the g2 and A2 components equal to those of g1 and A1; flag 7 sets the Euler angles of A1,
	// A2, the Gaussian and the Lorentzian tensor to the inverse of those of g1 (the diffusion angles are set equal to
	// those of g1); flag 8 sets them equal to those of g2; flags 10 and above copy the parameter from particle VF-10.
	// Afterwards the line widths are made positive, the sign of J0 is restored and, for the slow-motion calculation,
	// the order parameters of the mean-field potential are recalculated.
	// Input: Arguments, the vector of the optimiser. Outputs: the parameters of this object.
	// Used by: Discrepancy() in Discrepancy.cpp, at the start of every evaluation of the discrepancy, for each particle.
	void MagneticParameters::Arg2Par(const double* Arguments)
	{
// Independently varied values:
		for (int j=0; j<3; j++)
		{
					if(g1NArg(j, j)!=-1) 
					{
						if(g1VF(j,j)==1) g1(j, j)=Arguments[g1NArg(j, j)];
						if(g1VF(j,j)==2) g1(j, j)=2.0+0.001*Arguments[g1NArg(j, j)];
					}
					if(g2NArg(j, j)!=-1) 
					{
						if(g2VF(j,j)==1) g2(j, j)=Arguments[g2NArg(j, j)];
						if(g2VF(j,j)==2) g2(j, j)=2.0+0.001*Arguments[g2NArg(j, j)];
					}
					if(A1NArg(j, j)!=-1) A1(j, j)=Arguments[A1NArg(j, j)];
					if(A2NArg(j, j)!=-1) A2(j, j)=Arguments[A2NArg(j, j)];
					if(HGaussNArg(j, j)!=-1) HGauss(j, j)=Arguments[HGaussNArg(j, j)];
					if(HLorentzNArg(j, j)!=-1) HLorentz(j, j)=Arguments[HLorentzNArg(j, j)];

					if(g1EulerNArg[j]!=-1) g1Euler[j]=Arguments[g1EulerNArg[j]];
					if(g2EulerNArg[j]!=-1) g2Euler[j]=Arguments[g2EulerNArg[j]];
					if(A1EulerNArg[j]!=-1) A1Euler[j]=Arguments[A1EulerNArg[j]];
					if(A2EulerNArg[j]!=-1) A2Euler[j]=Arguments[A2EulerNArg[j]];
					if(HGaussEulerNArg[j]!=-1) HGaussEuler[j]=Arguments[HGaussEulerNArg[j]];
					if(HLorentzEulerNArg[j]!=-1) HLorentzEuler[j]=Arguments[HLorentzEulerNArg[j]];
					if(DiffusionEulerNArg[j]!=-1) DiffusionEuler[j]=Arguments[DiffusionEulerNArg[j]];
					if(LibrationEulerNArg[j]!=-1) LibrationEuler[j]=Arguments[LibrationEulerNArg[j]];
					if(LibrationAmplitudeNArg[j]!=-1) LibrationAmplitude[j]=Arguments[LibrationAmplitudeNArg[j]];
					if(LibDistrWidthNArg[j]!=-1) LibDistrWidth[j]=Arguments[LibDistrWidthNArg[j]];

		}
		// Dipolar vector and zero-field splitting parameters
		if(dxNArg!=-1) dx=Arguments[dxNArg];
		if(dyNArg!=-1) dy=Arguments[dyNArg];
		if(dzNArg!=-1) dz=Arguments[dzNArg];
		if(ZFS_ENArg!=-1) ZFS_E=Arguments[ZFS_ENArg];
		if(ZFS_AngleNArg!=-1) ZFS_Angle=Arguments[ZFS_AngleNArg];


		// Exchange and diffusion rates: the value for flag 1, 10^argument for flag 2
		if(J0NArg!=-1) 
		{
			if(J0VF==1) J0=Arguments[J0NArg];
			if(J0VF==2) J0=pow(10.0, Arguments[J0NArg]);
		}
		if(J2NArg!=-1) 
		{
			if(J2VF==1) J2=Arguments[J2NArg];
			if(J2VF==2) J2=pow(10.0, Arguments[J2NArg]);
		}
		if(JDiffusionNArg!=-1) 
		{
			if(JDiffusionVF==1) JDiffusion=Arguments[JDiffusionNArg];
			if(JDiffusionVF==2) JDiffusion=pow(10.0, Arguments[JDiffusionNArg]);
		}
		if(RxNArg!=-1) 
		{
			if(RxVF==1) Rx=Arguments[RxNArg];
			if(RxVF==2) Rx=pow(10.0, Arguments[RxNArg]);
		}
		if(RyNArg!=-1) 
		{
			if(RyVF==1) Ry=Arguments[RyNArg];
			if(RyVF==2) Ry=pow(10.0, Arguments[RyNArg]);
		}
		if(RzNArg!=-1) 
		{
			if(RzVF==1) Rz=Arguments[RzNArg];
			if(RzVF==2) Rz=pow(10.0, Arguments[RzNArg]);
		}
		// Conformer conversion rates
		if(NConformers>1)
		{
			for(int i=0; i<NConformers; i++)
				for(int j=i+1; j<NConformers; j++)
					if(ConversionRateNArg[i][j]!=-1) 
					{
						if(ConversionRateVF[i][j]==1) ConversionRate[i][j]=Arguments[ConversionRateNArg[i][j]];
						if(ConversionRateVF[i][j]==2) ConversionRate[i][j]=pow(10.0, Arguments[ConversionRateNArg[i][j]]);
					}
		}


		// Jump-diffusion residence times
		for(int i=0; i<3; i++)
			if(RTauNArg[i]!=-1)
			{
				RTau[i]=Arguments[RTauNArg[i]];
			}
// Hyperfine
		for (int i=0; i<HyperfineCount; i++)
		{
			if(HFS[i].HFConstantNArg!=-1)  HFS[i].HFConstant =Arguments[HFS[i].HFConstantNArg ];
		}

// DiscreteJumps:
		for(int i=0; i<NJumpAxes; i++)
		{
			if(JumpAxisFiNArg  [i]!=-1) JumpAxisFi  [i]=Arguments[JumpAxisFiNArg  [i]];
			if(JumpAxisTetaNArg[i]!=-1) JumpAxisTeta[i]=Arguments[JumpAxisTetaNArg[i]];
			if(JumpAngleNArg   [i]!=-1) JumpAngle   [i]=Arguments[JumpAngleNArg   [i]];
            if(JumpFrequencyNArg[i]!=-1)
			{
				if(JumpFrequencyVF[i]==1) JumpFrequency[i]=Arguments[JumpFrequencyNArg[i]];
				if(JumpFrequencyVF[i]==2) JumpFrequency[i]=pow(10.0, Arguments[JumpFrequencyNArg[i]]);
			}
		}

		// Mean-field potential coefficients are written back with Add
		for(int L=2; L<=MFPotential.Lmax; L+=2)
			for(int K=0; K<=L; K+=2)
			{
				if(MFPotential.GetNArg(L, K)!=-1)
				{
					MFPotential.Add(L, K, Arguments[MFPotential.GetNArg(L, K)], MFPotential.GetVF(L, K));
				}

			}


		// Field shift and amplitude factor
		if(ShiftNArg!=-1) Shift=Arguments[ShiftNArg];
		if(AmplitudeFactorNArg!=-1) AmplitudeFactor=Arguments[AmplitudeFactorNArg];
// Dependencies
// Variate flag 3
		if (g1VF(2,2)==3)
			g1(2,2)=g1Average3-g1(0,0)-g1(1,1);
		if (g2VF(2,2)==3)
			g2(2,2)=g2Average3-g2(0,0)-g2(1,1);
		if (A1VF(2,2)==3)
			A1(2,2)=A1Average3-A1(0,0)-A1(1,1);
		if (A2VF(2,2)==3)
			A2(2,2)=A2Average3-A2(0,0)-A2(1,1);
// Variate flag 4, value equal to Tensor_Y
	 for (int j=0; j<=2; j+=2)		// j=0, 2
	 {
		if (g1VF(j,j)==4)
			g1(j,j)=g1(1,1);
		if (g2VF(j,j)==4)
			g2(j,j)=g2(1,1);
		if (A1VF(j,j)==4)
			A1(j,j)=A1(1,1);
		if (A2VF(j,j)==4)
			A2(j,j)=A2(1,1);
		if (HGaussVF(j,j)==4)
			HGauss(j,j)=HGauss(1,1);
		if (HLorentzVF(j,j)==4)
			HLorentz(j,j)=HLorentz(1,1);
	    if (LibrationAmplitudeVF[j]==4)
		    LibrationAmplitude[j]=LibrationAmplitude[1];
	 }
	if (RxVF==4)
		Rx=Ry;
	if (RzVF==4)
		Rz=Ry;

	if (RTauVF[0]==4)
		RTau[0]=RTau[1];
	if (RTauVF[2]==4)
		RTau[2]=RTau[1];

// Variate flag 5, value equal to Tensor_X
	 for (int j=1; j<=2; j++)		// j=1, 2
	 {
		if (g1VF(j,j)==5)
			g1(j,j)=g1(0,0);
		if (g2VF(j,j)==5)
			g2(j,j)=g2(0,0);
		if (A1VF(j,j)==5)
			A1(j,j)=A1(0,0);
		if (A2VF(j,j)==5)
			A2(j,j)=A2(0,0);
		if (HGaussVF(j,j)==5)
			HGauss(j,j)=HGauss(0,0);
		if (HLorentzVF(j,j)==5)
			HLorentz(j,j)=HLorentz(0,0);
	    if (LibrationAmplitudeVF[j]==5)
		    LibrationAmplitude[j]=LibrationAmplitude[0];
	 }
	if (RyVF==5)
		Ry=Rx;
	if (RzVF==5)
		Rz=Rx;

	if (RTauVF[1]==5)
		RTau[1]=RTau[0];
	if (RTauVF[2]==5)
		RTau[2]=RTau[0];

// Variate flag 6, value equal to Tensor1
	 for (int j=0; j<=2; j++)		// j=0, 1, 2
	 {
		if (g2VF(j,j)==6)
			g2(j,j)=g1(j,j);
		if (A2VF(j,j)==6)
			A2(j,j)=A1(j,j);
	 }
// Variate flag 7, EulerAngle inverse to those of g1-tensor
	 for (int j=0; j<=2; j++)		// j=0, 1, 2
	 {
		if (A1EulerVF[j]==7)
			A1Euler[j]=-g1Euler[2-j];
		if (A2EulerVF[j]==7)
			A2Euler[j]=-g1Euler[2-j];
		if (HGaussEulerVF[j]==7)
			HGaussEuler[j]=-g1Euler[2-j];
		if (HLorentzEulerVF[j]==7)
			HLorentzEuler[j]=-g1Euler[2-j];
		if (DiffusionEulerVF[j]==7)
			DiffusionEuler[j]=g1Euler[j];	 }
// Variate flag 8, EulerAngle equal to those of g2-tensor
	 for (int j=0; j<=2; j++)		// j=0, 1, 2
	 {
		if (A1EulerVF[j]==8)
			A1Euler[j]=g2Euler[j];
		if (A2EulerVF[j]==8)
			A2Euler[j]=g2Euler[j];
		if (HGaussEulerVF[j]==8)
			HGaussEuler[j]=g2Euler[j];
		if (HLorentzEulerVF[j]==8)
			HLorentzEuler[j]=g2Euler[j];
		if (DiffusionEulerVF[j]==8)
			DiffusionEuler[j]=g2Euler[j];
	 }
 // Variate flags 10+
	 // Parameters tied to particle VF-10 are copied from MagneticParameterArray[VF-10]
	 for(int j=0;j<=2; j++)
	 {
		 if (g1VF(j, j)>=10)
			 g1(j,j)=MagneticParameterArray[g1VF(j,j)-10].g1(j,j);
		 if (g2VF(j, j)>=10)
			 g2(j,j)=MagneticParameterArray[g2VF(j,j)-10].g2(j,j);
		 if (A1VF(j, j)>=10)
			 A1(j,j)=MagneticParameterArray[A1VF(j,j)-10].A1(j,j);
		 if (A2VF(j, j)>=10)
			 A2(j,j)=MagneticParameterArray[A2VF(j,j)-10].A2(j,j);
		 if (HGaussVF(j, j)>=10)
			 HGauss(j,j)=MagneticParameterArray[HGaussVF(j,j)-10].HGauss(j,j);
		 if (HLorentzVF(j, j)>=10)
			 HLorentz(j,j)=MagneticParameterArray[HLorentzVF(j,j)-10].HLorentz(j,j);

		 if (g1EulerVF[j]>=10)
			 g1Euler[j]=MagneticParameterArray[g1EulerVF[j]-10].g1Euler[j];
		 if (g2EulerVF[j]>=10)
			 g2Euler[j]=MagneticParameterArray[g2EulerVF[j]-10].g2Euler[j];
		 if (A1EulerVF[j]>=10)
			 A1Euler[j]=MagneticParameterArray[A1EulerVF[j]-10].A1Euler[j];
		 if (A2EulerVF[j]>=10)
			 A2Euler[j]=MagneticParameterArray[A2EulerVF[j]-10].A2Euler[j];
		 if (HGaussEulerVF[j]>=10)
			 HGaussEuler[j]=MagneticParameterArray[HGaussEulerVF[j]-10].HGaussEuler[j];
		 if (HLorentzEulerVF[j]>=10)
			 HLorentzEuler[j]=MagneticParameterArray[HLorentzEulerVF[j]-10].HLorentzEuler[j];
		 if (DiffusionEulerVF[j]>=10)
			 DiffusionEuler[j]=MagneticParameterArray[DiffusionEulerVF[j]-10].DiffusionEuler[j];
	 }
	 if (dxVF>=10)
		 dx=MagneticParameterArray[dxVF-10].dx;
	 if (dyVF>=10)
		 dy=MagneticParameterArray[dyVF-10].dy;
	 if (dzVF>=10)
		 dz=MagneticParameterArray[dzVF-10].dz;
	 
	 // Conversion rates tied to another particle
	 if(NConformers>1)
	 {
		 for(int i=0; i<NConformers; i++)
			 for(int j=i+1; j<NConformers; j++)
			 {
				 if (ConversionRateVF[i][j]>=10)
					 ConversionRate[i][j]=MagneticParameterArray[ConversionRateVF[i][j]-10].ConversionRate[i][j];
			 }
	 }

	 if (J0VF>=10)
		 J0=MagneticParameterArray[J0VF-10].J0;

	 if (RxVF>=10)
		 Rx=MagneticParameterArray[RxVF-10].Rx;
	 if (RyVF>=10)
		 Ry=MagneticParameterArray[RyVF-10].Ry;
	 if (RzVF>=10)
		 Rz=MagneticParameterArray[RzVF-10].Rz;
	 for(int i=0; i<3; i++)
	 {
	   if (RTauVF[i]>=10)
		 RTau[i]=MagneticParameterArray[RTauVF[i]-10].RTau[i];
	 }
	// Hyperfine constants tied to another particle, scaled with HFConstantCoefficient
	for (int i=0; i<HyperfineCount; i++)
	{
		if(HFS[i].HFConstantVF>=10)
			HFS[i].HFConstant=MagneticParameterArray[HFS[i].HFConstantVF-10].HFS[i].HFConstant*HFS[i].HFConstantCoefficient;
	}
  // Discrete Jumps:
	for(int i=0; i<NJumpAxes; i++)
	{
		if(JumpAxisFiVF   [i]>=10) JumpAxisFi   [i]=MagneticParameterArray[JumpAxisFiVF   [i]-10].JumpAxisFi   [i];
		if(JumpAxisTetaVF [i]>=10) JumpAxisTeta [i]=MagneticParameterArray[JumpAxisTetaVF [i]-10].JumpAxisTeta [i];
		if(JumpAngleVF    [i]>=10) JumpAngle    [i]=MagneticParameterArray[JumpAngleVF    [i]-10].JumpAngle    [i];
		if(JumpFrequencyVF[i]>=10) JumpFrequency[i]=MagneticParameterArray[JumpFrequencyVF[i]-10].JumpFrequency[i];

	}

	 if (AmplitudeFactorVF>=10)
		 AmplitudeFactor=MagneticParameterArray[AmplitudeFactorVF-10].AmplitudeFactor*AmplitudeCoefficient;
	 // All particles use the field shift of particle 0
	 Shift=MagneticParameterArray[0].Shift;

	 // For the slow-motion calculation recalculate the order parameters of the mean-field potential
	 if(Liouville)
	 {
		 MFPotential.MyCount=Count;
//		 MFPotential.CalculateOP(max(LOmax, LEmax));
		 MFPotential.CalculateOP_GK(max(LOmax, LEmax));
	 }
//	Check for negative linewidths
	 HLorentz(0,0)=fabs(HLorentz(0,0));
	 HLorentz(1,1)=fabs(HLorentz(1,1));
	 HLorentz(2,2)=fabs(HLorentz(2,2));

	 HGauss(0,0)=fabs(HGauss(0,0));
	 HGauss(1,1)=fabs(HGauss(1,1));
	 HGauss(2,2)=fabs(HGauss(2,2));

// Checking the sign of J0;
	 J0=fabs(J0)*J0sign;

	}




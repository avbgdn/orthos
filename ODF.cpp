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


// ODF.cpp: orientation distribution function (ODF) of the paramagnetic particles in the sample, and the Wigner
// rotation functions it needs. Contents: factorial helpers, the reduced Wigner d-function wigner_d, the complex Wigner
// D-function ComplexWigner, and the member functions of the class OrderParameters (declared in BiradicalODF.h). That class
// holds the ODF expansion coefficients A_jmn over Wigner D-functions, the sample orientation angles of each recorded
// spectrum, the variate flags and optimiser argument numbers of the fitted ODF parameters, the ODF input and output
// file format, and the negative-volume penalty of the ODF.
// Used by: main.cpp (creates and reads the ODF array, registers the fitted parameters, writes the result),
// Discrepancy.cpp (updates the ODF parameters and sample angles at every evaluation), IntegratedSpectrum.cpp and
// VoitlanderIntegration.cpp (ODF_reduced, the ODF weight of each particle orientation), ReadExperimental.cpp
// (WriteToFileBrief); ComplexWigner and wigner_d are also used by Potential.cpp, SphericalTensor.cpp and
// LiouvilleMatrixElementB.cpp.

#include "BiradicalODF.h"
#include <math.h>
#include <cmath>
#include <complex>
#include <stdlib.h>

// Factorial n! in integer arithmetic; returns -1 for negative n.
// Used by: nothing in this code base (wigner_d uses dfact instead).
// factorial function
long int fact(int n)
{
	// recursive definition: 0!=1 and n!=n*(n-1)!
	if(n==0) return 1;
	else if (n>=0) return n*fact(n-1);
	else return -1;
}

// Factorial n! in long double arithmetic (wider range than fact); returns -1 for negative n.
// Used by: wigner_d (this file).
long double dfact(int n)
{
	// recursive definition: 0!=1 and n!=n*(n-1)!
	if(n==0) return 1;
	else if (n>=0) return n*dfact(n-1);
	else return -1;
}

// Reduced Wigner d-function d^j_mn(beta), the polar-angle part of the Wigner rotation matrix D^j_mn: an explicit
// alternating sum over k of powers of cos(beta/2) and sin(beta/2) divided by factorials, times a square-root prefactor.
// Inputs: integer rank j and indices m, n, polar angle beta in radians. Returns 0 if |m| or |n| exceeds j.
// Used by: ComplexWigner (this file); precalculate_wigner (LiouvilleMatrixElementB.cpp, rank j=2 only).
// small d-Wigner function
double wigner_d(int j, int m, int n, double beta)
{
	int NMIN, NMAX, k;
	double sum=0;
	// outside the allowed index range the function is zero
	if (abs(m)>j) return 0;
	if (abs(n)>j) return 0;

	// The choice of summation bounds
	if (m+n>0) NMIN=m+n;
	else       NMIN=0;

	if (j+m<j+n) NMAX=j+m;
	else         NMAX=j+n;

	// alternating sum over k (the sign is (-1)^k, chosen by the parity of k)
	for(k=NMIN; k<=NMAX; k++)
	{
        if (k%2)				// k nechet
			sum-=pow(cos(beta/2), 2*k-m-n)*pow(sin(beta/2), 2*j+m+n-2*k)/(double)(dfact(k)*dfact(j+m-k)*dfact(j+n-k)*dfact(k-m-n));
		else					// k chet
			sum+=pow(cos(beta/2), 2*k-m-n)*pow(sin(beta/2), 2*j+m+n-2*k)/(double)(dfact(k)*dfact(j+m-k)*dfact(j+n-k)*dfact(k-m-n));
	}

	// prefactor (-1)^(j+m)*sqrt((j+m)!(j-m)!(j+n)!(j-n)!)
	if((j+m)%2)                   // j+m nechet
		sum*=-sqrt((double)dfact(j+m)*dfact(j-m)*dfact(j+n)*dfact(j-n));
	else                        // j+m chet
		sum*= sqrt((double)dfact(j+m)*dfact(j-m)*dfact(j+n)*dfact(j-n));

	return sum;
}

// Complex Wigner D-function D^j_mn(alpha,beta,gamma)=d^j_mn(beta)*exp(-i*m*alpha-i*n*gamma) for the Euler angles
// alpha, beta, gamma in radians (the phase convention is the one in the return statement below).
// Used by: OrderParameters::AssignTjpqi, ODF_value and ODF_reduced (this file); Potential::CalculateOP,
// CalculateOP_GK and CalculateOP_Jumps_GK (Potential.cpp); SphericalTensor::Rotate and RotateMe (SphericalTensor.cpp);
// MatrixElementJumpsBSM_Pot and MatrixElementJumpsBSM_Iso (LiouvilleMatrixElementB.cpp).
// Complex Wigner D-function
std::complex<double> ComplexWigner(int j, int m, int n, double alpha, double beta, double gamma)
{
	std::complex<double> I(0, 1);
	std::complex<double> Value;
	// reduced d-function times the two phase factors
	return wigner_d(j, m, n, beta)*exp(-gamma*I*(double)n-alpha*I*(double)m);
}

//////////////////////////////////////////////////////////////////////////
////////  CLASS FOR ODF AND ORDER PARAMETER MANIPULATION /////////////////
//////////////////////////////////////////////////////////////////////////


// OrderParameters holds the orientation distribution function (ODF) of one kind of paramagnetic particle: the expansion
// coefficients A_jmn over Wigner D-functions D^j_mn (j even, up to Jmax<=MAX_RANK), stored in
// Ajmn_container[j][m+Jmax][n+Jmax] so that negative m and n map to non-negative array indices. For each of the
// NOrientations recorded spectra it also keeps the Euler angles of the sample orientation (alpha_SF, beta_SF, gamma_SF),
// the products T_jpqi (see AssignTjpqi) and the ODF weights CurrentODF. Fitted parameters carry a variate flag (0: fixed;
// 1 or 2: fitted, placed in the optimiser argument array by Par2Arg; 10 or more: copied from particle flag-10 by Arg2Par)
// and an argument number (-1 if not fitted).
// Used by: main.cpp, Discrepancy.cpp, IntegratedSpectrum.cpp, VoitlanderIntegration.cpp, ReadExperimental.cpp,
// ReadNamenam.cpp.
//////////////////////////////////////////////////////////////////////////
	// Constructor: sets all coefficients A_jmn, products T_jpqi and sample angles to zero, the rank Jmax to 0 and the number
	// of sample orientations to 1 (isotropic ODF, single spectrum).
	// Used by: main.cpp, which allocates the ODF array with new OrderParameters[NParticles].
	OrderParameters::OrderParameters()
	{
		// clear the coefficient and product containers and the sample angles
		for (int j=0; j<=MAX_RANK; j++)
			for(int m=0; m<=2*MAX_RANK; m++)
				for(int n=0; n<=2*MAX_RANK; n++)
				{
					Ajmn_container[j][m][n]=0;
					for(int i=0; i<MAX_ORIENTATIONS; i++)
					{
						Tjpqi_container[j][m][n][i]=0;
						alpha_SF[i]=0;
						beta_SF[i]=0;
						gamma_SF[i]=0;
					}
				}
		// by default the ODF is rank 0 with a single sample orientation
		Jmax=0;
		NOrientations=1;
	}
//////////////////////////////////////////////////////////////////////////
// Destructor: nothing to release, all storage is in fixed-size members.
// Used by: compiler-generated destruction of OrderParameters objects (main.cpp allocates the ODF array with new[]).
OrderParameters::~OrderParameters()
	{}

//////////////////////////////////////////////////////////////////////////
	// Precomputes T_jpqi=A_jpq*D^j_p0(alpha_SF[i],beta_SF[i],gamma_SF[i]) for every even j<=Jmax, p and q in -j..j and each
	// sample orientation i, so that ODF_reduced only has to multiply by the remaining Wigner function. Angles in radians.
	// Used by: Discrepancy (Discrepancy.cpp), once per particle at every evaluation, after Arg2Par, Realize, SetEuler_SF
	// and ConvertAnglesToRadians.
	void OrderParameters::AssignTjpqi()
	{
		// loop over even rank, both indices and all sample orientations
		for (int j=0; j<=Jmax; j+=2)
			for(int p=-j; p<=j; p++)
				for(int q=-j; q<=j; q++)
					for(int i=0; i<NOrientations; i++)
						Tjpqi_container[j][p+Jmax][q+Jmax][i]=Ajmn_container[j][p+Jmax][q+Jmax]*ComplexWigner(j, p, 0, alpha_SF[i], beta_SF[i], gamma_SF[i]);
	}
//////////////////////////////////////////////////////////////////////////
	// Returns the ODF coefficient A_jmn; m and n may be negative (they are shifted by Jmax inside the array).
	// Used by: ODF_value, WriteToFile, WriteToFileBrief, Par2Arg and Arg2Par (this file); Arg2Par also calls it on another
	// particle through ODFArray.
	std::complex<double> OrderParameters::A_jmn(int j, int m, int n)
	{
		return Ajmn_container[j][m+Jmax][n+Jmax];
	}
	// Returns the variate flag of A_jmn (0 fixed, 1 or 2 fitted, 10 or more copied from particle flag-10).
	// Used by: WriteToFile, WriteToFileBrief, Par2Arg and Arg2Par (this file).
	int OrderParameters::A_jmnVF(int j, int m, int n)
	{
		return Ajmn_containerVF[j][m+Jmax][n+Jmax];
	}
	// Returns the position of A_jmn in the optimiser argument array, as stored by Par2Arg (-1 if not fitted).
	// Used by: Arg2Par (this file).
	int OrderParameters::A_jmnNArg(int j, int m, int n)
	{
		return Ajmn_containerNArg[j][m+Jmax][n+Jmax];
	}
//////////////////////////////////////////////////////////////////////////
	// Returns the product T_jpqi stored by AssignTjpqi for sample orientation i.
	// Used by: ODF_reduced (this file).
	std::complex<double> OrderParameters::T_jpqi(int j, int p, int q, int i)
	{
		return Tjpqi_container[j][p+Jmax][q+Jmax][i];
	}
//////////////////////////////////////////////////////////////////////////
	// Evaluates the ODF itself at the Euler angles ksi, teta, psi (radians): the real part of the sum of
	// A_jmn*D^j_mn(ksi,teta,psi) over even j, m and n.
	// Used by: Arg2Par (this file), to integrate the ODF on a grid for the negative-volume estimate (a call in main.cpp is
	// commented out).
	double OrderParameters::ODF_value(double ksi, double teta, double psi)
	{
		std::complex<double> sum(0,0);
		// sum of coefficient times Wigner function over even j, m and n
		for(int j=0; j<=Jmax; j+=2)
			for(int m=-j; m<=j; m+=2)
				for(int n=-j; n<=j; n+=2)
				{
					sum+=A_jmn(j, m, n)*ComplexWigner(j, m, n, ksi, teta, psi);
				}
		return real(sum);

	}
	// Fills CurrentODF[i] for every sample orientation i: the sum over even j and over p, q of
	// Re(T_jpqi*D^j_0q(0,teta,ksi)), i.e. the number density of particles at orientation angles teta and ksi (radians).
	// teta and ksi are the Euler angles of the Field-to-Molecule rotation at the current integration point; the caller
	// multiplies the spectrum of that orientation by CurrentODF[i] and by the quadrature weight.
	// Used by: CalculateIntegratedSpectrumLiouville (IntegratedSpectrum.cpp) and CreateVertex (VoitlanderIntegration.cpp).
	void OrderParameters::ODF_reduced(double teta, double ksi)
// Given the orientation distribution function A_jmn,
// and the set of angles _SF for all sample orientations,
// this subroutine calculates the number density of
// particles which are oriented in resonator by angles ksi, teta.
// The calculated values are stored in array CurrentODF[0..NOrientations],
// to be used in the block of spectra calculation.
	{
		// clear the accumulators
		for (int i=0; i<NOrientations; i++)
			CurrentODF[i]=0;

		// add the contribution of every coefficient for every sample orientation
		for(int j=0; j<=Jmax; j+=2)
			for(int p=-j; p<=j; p++)
				for(int q=-j; q<=j; q++)
					for (int i=0; i<NOrientations; i++)
						CurrentODF[i]+=real(T_jpqi(j, p, q, i)*ComplexWigner(j, 0, q, 0, teta, ksi));
		return;
	}
//////////////////////////////////////////////////////////////////////////
	// Stores the ODF coefficient A_jmn (m and n are shifted by Jmax inside the array).
	// Used by: ReadFromFile and Arg2Par (this file).
	void OrderParameters::Set_A_jmn(int j, int m, int n, std::complex<double> value)
// Sets the value of specified order parameter
	{
		Ajmn_container[j][m+Jmax][n+Jmax]=value;
	}
	// Stores the variate flag of A_jmn.
	// Used by: ReadFromFile (this file).
	void OrderParameters::Set_A_jmnVF(int j, int m, int n, int value)
// Sets the value of specified order parameter variate flag
	{
		Ajmn_containerVF[j][m+Jmax][n+Jmax]=value;
	}
	// Stores the position of A_jmn in the optimiser argument array (-1 if not fitted).
	// Used by: Par2Arg (this file).
	void OrderParameters::Set_A_jmnNArg(int j, int m, int n, int value)
// Sets the value of specified order parameter number in Arg array
	{
		Ajmn_containerNArg[j][m+Jmax][n+Jmax]=value;
	}
//////////////////////////////////////////////////////////////////////////
	// Imposes reality of the ODF: for m<=0 sets A_j,m,n to the complex conjugate of A_j,-m,-n (only the m,n>=0 half is
	// read from the input file or fitted), and makes A_j,0,0 real.
	// Used by: Discrepancy (Discrepancy.cpp), once per particle at every evaluation, after Arg2Par.
	void OrderParameters::Realize()
// This function makes the ODF A_jmn real by formula (3.5-6) in dissertation
	{
		for (int j=0; j<=Jmax; j++)
		{
			// coefficients with negative m and n are the complex conjugates of the positive ones
			for(int m=-j; m<=0; m++)
// note that in fact we do not have non-zero coefficients A_jmn with m and n of different sign, 
// so we can make the next line: for(int n=-j; n<=0; n++)
				for(int n=-j; n<=j; n++)
					Ajmn_container[j][m+Jmax][n+Jmax]=conj(Ajmn_container[j][-m+Jmax][-n+Jmax]);
			// A_j00 is its own conjugate, so only its real part is kept
			Ajmn_container[j][0+Jmax][0+Jmax]=real(Ajmn_container[j][0+Jmax][0+Jmax]);
		}
	}

	// Reads the ODF parameters of particle number Count from the ODF input file: finds the block 'Particle: Count', then
	// reads ODFRank (Jmax), Penalty, DelHi0, DelHi1 and DelHi2 (value and variate flag) and the OrderParameters table
	// (value and variate flag of A_jmn for even j, m and n with 0<=m,n<=j; the coefficients with negative m and n are not
	// read, see Realize). Terminates the program if the file cannot be opened or the rank is outside 0..MAX_RANK; if the
	// block of this particle is missing it sets an isotropic ODF.
	// Used by: main (main.cpp), once per particle.
	void OrderParameters::ReadFromFile(char* filename)
	{
	// open the ODF input file
	ifstream ODF_par;
	ODF_par.open (filename, ios::in);
	if(!ODF_par.is_open())
	{
		std::cout << "Could not open file " << filename << std::endl;
		exit(0);
	}
	// skip lines until the block of this particle is found
	bool particlefound=false;
	while (!ODF_par.eof())
	{
		char buffer[MAX_STRLEN], buffer2[MAX_STRLEN], buffer3[MAX_STRLEN];
		ODF_par.getline(buffer, MAX_STRLEN);
		sscanf(buffer, "%s", buffer2);
		int npart;
		if (!strcmp(buffer2, "Particle:"))
		{
			sscanf(buffer, "%s %d", buffer3, &npart);
			if(npart==Count)
			{
				particlefound=true;
				break;
			}
		}
	}
// the ODF parameter file has no "Particle:" block for this particle: force an isotropic ODF (rank 0) instead of silently keeping whatever was last read
	if(!particlefound)
	{
		std::cout << "Warning: file " << filename << " has no ODF block for particle " << Count << ". Forcing isotropic ODF (rank 0)." << std::endl;
		Jmax=0;
		Penalty=0.0;
		Set_A_jmn(0, 0, 0, 1.0);
		Set_A_jmnVF(0, 0, 0, 0);
		return;
	}
	// read the parameters of this particle line by line until the order-parameter table has been read
	while (!ODF_par.eof())
	{
		char buffer[MAX_STRLEN], buffer2[MAX_STRLEN], buffer3[MAX_STRLEN];
		ODF_par.getline(buffer, MAX_STRLEN);
		sscanf(buffer, "%s", buffer2);
		if (!strcmp(buffer2, "ODFRank:"))
		{
			sscanf(buffer, "%s %d", buffer3, &Jmax);
// reject a rank outside the fixed-size Ajmn_container/Tjpqi_container bound instead of silently corrupting adjacent heap data
			if((Jmax<0)||(Jmax>MAX_RANK))
			{
				std::cout << "Error reading " << filename << ": particle " << Count << " requests ODFRank " << Jmax << ", outside the valid range [0, " << MAX_RANK << "]." << std::endl;
				exit(0);
			}
		}
		// scalar parameters: penalty weight and the three angle errors with their variate flags
		if (!strcmp(buffer2, "Penalty:"))
			sscanf(buffer, "%s %lf", buffer3, &Penalty);
		if (!strcmp(buffer2, "DelHi0:"))
			sscanf(buffer, "%s %lf %d", buffer3, &DelHi0, &DelHi0VF);
		if (!strcmp(buffer2, "DelHi1:"))
			sscanf(buffer, "%s %lf %d", buffer3, &DelHi1, &DelHi1VF);
		if (!strcmp(buffer2, "DelHi2:"))
			sscanf(buffer, "%s %lf %d", buffer3, &DelHi2, &DelHi2VF);
		// order-parameter table: one block per even j, one row per even m, pairs (value, variate flag) for even n
		if (!strcmp(buffer2, "OrderParameters:"))
		{
			ODF_par.getline(buffer, MAX_STRLEN);    //       m      A_jm0        A_jm2        A_jm4        A_jm6        A_jm8
			for (int j=0; j<=Jmax; j+=2)
			{
				ODF_par.getline(buffer, MAX_STRLEN);    //        j=0  
				for(int m=0; m<=j; m+=2)
				{
					int mtemp;
					ODF_par >> mtemp;       // reading dummy value of m
					for(int n=0; n<=j; n+=2)
					{
						double dtemp;
						ODF_par >> dtemp;    // Order Parameter
						ODF_par >> mtemp;    // Variate Flag
						Set_A_jmn(j, m, n, dtemp);
						Set_A_jmnVF(j, m, n, mtemp);
//						std::cout << A_jmn(j, m, n) << "\t";
					}	
				    ODF_par.getline(buffer, MAX_STRLEN);    //        end line  
				}
			}
			return;
		}
		// the block of another particle starts: stop reading
		if (!strcmp(buffer2, "Particle:"))
		{
			int npart;
			sscanf(buffer, "%s %d", buffer3, &npart);
			if(npart!=Count) break;
		}

	}


	}
	// Converts the stored sample Euler angles (alpha_SF, beta_SF, gamma_SF) from radians to degrees, for printing.
	// Used by: Discrepancy (Discrepancy.cpp, before the temporary output files are written) and main (main.cpp, before the
	// final output is written).
	void OrderParameters::ConvertAnglesToDegrees(void)
	{
		// convert the angles of all orientation slots
		for(int i=0; i<MAX_ORIENTATIONS; i++)
		{
				alpha_SF[i]*=180/PI;
				beta_SF[i]*=180/PI;
				gamma_SF[i]*=180/PI;
		}
	}
	// Converts the stored sample Euler angles (alpha_SF, beta_SF, gamma_SF) from degrees to radians.
	// Used by: main (main.cpp, after reading the input and after writing the output) and Discrepancy (Discrepancy.cpp, after
	// SetEuler_SF and after writing the temporary output files).
	void OrderParameters::ConvertAnglesToRadians(void)
	{
		// convert the angles of all orientation slots
		for(int i=0; i<MAX_ORIENTATIONS; i++)
		{
				alpha_SF[i]/=180/PI;
				beta_SF[i]/=180/PI;
				gamma_SF[i]/=180/PI;
		}
	}
	// Builds the sample Euler angles of every recorded spectrum from its nominal rotation angle hi_SF[i], its rotation axis
	// axis_SF[i] (0: x, 1: y, 2: z) and the fitted angle errors DelHi0, DelHi1 and DelHi2. The angles are produced in degrees;
	// Discrepancy.cpp converts them to radians right after. For the slow-motion (Liouville) case beta_SF[i] is the angle psi
	// passed to CalculateLiouvillianB (Discrepancy.cpp).
	// Used by: Discrepancy (Discrepancy.cpp), once per particle at every evaluation.
	void OrderParameters::SetEuler_SF(void)
	{
		// choose the Euler angles according to the rotation axis of each orientation
		for(int i=0; i<MAX_ORIENTATIONS; i++)
		{
			if(axis_SF[i]==0)
			{
				// x-rotation
				alpha_SF[i]=90;
				beta_SF[i]=-(hi_SF[i]+DelHi0);
				gamma_SF[i]=0;
			}
			else if(axis_SF[i]==1)
			{
				// y-rotation
				alpha_SF[i]=0;
				beta_SF[i]=(hi_SF[i]+DelHi1);
				gamma_SF[i]=0;
			}
			else if(axis_SF[i]==2)
			{
				// z-rotation
				alpha_SF[i]=(hi_SF[i]+DelHi2);
				beta_SF[i]=90;
				gamma_SF[i]=90;
			}
		}


	}
	// Appends the ODF block of this particle to the file, in the format that ReadFromFile reads: Particle number, ODFRank,
	// Penalty with the NegativeVolume values, DelHi0..2 with variate flags, and the table of A_jmn (real parts, even j,
	// m and n>=0, each with its variate flag).
	// Used by: main (main.cpp, final output) and Discrepancy (Discrepancy.cpp, temporary file ODF.temp).
	void OrderParameters::WriteToFile(char* filename)
	{
		FILE *file;
		// open the file for appending
		file=fopen(filename, "a+");
		// header lines: particle number, rank, penalty and angle errors
		fprintf(file, "\nParticle: %d\n", Count);
		fprintf(file, "ODFRank:       %d\n", Jmax);
		fprintf(file, "Penalty:       %lf      NegativeVolume: %e  (ODF3: %e)\n", Penalty, NegativeVolume, NegativeVolumeODF3);
		fprintf(file, "DelHi0:        %3.5lf  %d\n", DelHi0, DelHi0VF);
		fprintf(file, "DelHi1:        %3.5lf  %d\n", DelHi1, DelHi1VF);
		fprintf(file, "DelHi2:        %3.5lf  %d\n", DelHi2, DelHi2VF);
		// table of coefficients: one block per even j, one row per even m
		fprintf(file, "\nOrderParameters:\n");
		fprintf(file, "      m      A_jm0        A_jm2        A_jm4        A_jm6        A_jm8\n");
		for(int j=0; j<=Jmax; j+=2)
		{
		fprintf(file, " j=%d\n", j);
		for (int m=0; m<=j; m+=2)
		{
			fprintf(file, "      %d  ", m);
			for (int n=0; n<=j; n+=2)
				fprintf(file, " %2.6lf  %d  ", real(A_jmn(j, m, n)), A_jmnVF(j, m, n));
			fprintf(file, "\n");
		}

		}
		fclose(file);

	}

	// Appends a shorter report of the ODF to the file: the same quantities as WriteToFile, without the Particle header, with
	// the DelHi values on one line and each row of the coefficient table starting with the letter c.
	// Used by: WriteOutput (ReadExperimental.cpp).
	void OrderParameters::WriteToFileBrief(char* filename)
	{
		FILE *file;
		// open the file for appending
		file=fopen(filename, "a+");
		// header lines: rank, penalty and angle errors
		fprintf(file, "ODFRank: %d, ", Jmax);
		fprintf(file, "Penalty:  %lf,  NegativeVolume: %e  (ODF3: %e)\n", Penalty, NegativeVolume, NegativeVolumeODF3);
		fprintf(file, "DelHi0,1,2:        %3.5lf  %d  %3.5lf  %d  %3.5lf  %d\n", DelHi0, DelHi0VF, DelHi1, DelHi1VF, DelHi2, DelHi2VF);
		// table of coefficients: one block per even j, one row per even m
		fprintf(file, "OrderParameters:\n");
		fprintf(file, "      m      A_jm0        A_jm2        A_jm4        A_jm6        A_jm8\n");
		for(int j=0; j<=Jmax; j+=2)
		{
		fprintf(file, " j=%d\n", j);
		for (int m=0; m<=j; m+=2)
		{
			fprintf(file, "c      %d  ", m);
			for (int n=0; n<=j; n+=2)
				fprintf(file, " %2.6lf  %d  ", real(A_jmn(j, m, n)), A_jmnVF(j, m, n));
			fprintf(file, "\n");
		}

		}
		fprintf(file, "\n");
		fclose(file);

	}
	
	// Registers the fitted ODF parameters in the optimiser argument array. For DelHi0, DelHi1, DelHi2 and every A_jmn (even
	// j, m and n with 0<=m,n<=j) whose variate flag is 1 or 2 it stores the current value in Arguments[*NArg], sets the
	// bounds lb and ub (angles: MIN_ANGLE..MAX_ANGLE; A_jmn: -(2j+1)..2j+1), remembers the argument number and advances
	// *NArg with NextArg. Parameters that are not fitted get argument number -1.
	// Used by: main (main.cpp), once per particle after the magnetic parameters have been registered.
	void OrderParameters::Par2Arg(double* Arguments, double *ub, double *lb, int *NArg)
	{
		// angle errors: the fitted ones are placed in the argument array
		if((DelHi0VF==1)||(DelHi0VF==2))
		{
			Arguments[*NArg]=DelHi0;
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			DelHi0NArg=*NArg;
			NextArg(NArg);
		}
		else DelHi0NArg=-1;
		if((DelHi1VF==1)||(DelHi1VF==2))
		{
			Arguments[*NArg]=DelHi1;
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			DelHi1NArg=*NArg;
			NextArg(NArg);
		}
		else DelHi1NArg=-1;
		if((DelHi2VF==1)||(DelHi2VF==2))
		{
			Arguments[*NArg]=DelHi2;
			lb[*NArg]=MIN_ANGLE;
			ub[*NArg]=MAX_ANGLE;
			DelHi2NArg=*NArg;
			NextArg(NArg);
		}
		else DelHi2NArg=-1;
		// order parameters, registered in the same way
		for (int j=0; j<=Jmax; j+=2)
			for(int m=0; m<=j; m+=2)
				for(int n=0; n<=j; n+=2)
				{
					if((A_jmnVF(j, m, n)==1)||(A_jmnVF(j, m, n)==2))
					{
						Arguments[*NArg]=real(A_jmn(j, m, n));
						lb[*NArg]=-2*j-1;
						ub[*NArg]=2*j+1;
						Set_A_jmnNArg(j, m, n, (*NArg));
						NextArg(NArg);
					}
					else Set_A_jmnNArg(j, m, n, -1);;
				}
}
	// Inverse of Par2Arg: copies the optimiser arguments back into DelHi0..2 and A_jmn, copies the parameters whose variate
	// flag is 10 or more from particle flag-10 (through ODFArray), and recomputes NegativeVolume and NegativeVolumeODF3,
	// which measure how much of the ODF is negative; Discrepancy uses NegativeVolume in its penalty coefficient.
	// Used by: Discrepancy (Discrepancy.cpp), once per particle at every evaluation.
	void OrderParameters::Arg2Par(const double* Arguments)
	{
		// angle errors from the argument array
		if(DelHi0NArg!=-1) DelHi0=Arguments[DelHi0NArg];
		if(DelHi1NArg!=-1) DelHi1=Arguments[DelHi1NArg];
		if(DelHi2NArg!=-1) DelHi2=Arguments[DelHi2NArg];

		// angle errors tied to another particle are copied from it
		if (DelHi0VF>=10)
	 		DelHi0=ODFArray[DelHi0VF-10].DelHi0;
		if (DelHi1VF>=10)
	 		DelHi1=ODFArray[DelHi1VF-10].DelHi1;
		if (DelHi2VF>=10)
	 		DelHi2=ODFArray[DelHi2VF-10].DelHi2;

		// order parameters: from the argument array, or copied from another particle
		for (int j=0; j<=Jmax; j+=2)
			for(int m=0; m<=j; m+=2)
				for(int n=0; n<=j; n+=2)
				{
					if(A_jmnNArg(j,m,n)!=-1) Set_A_jmn(j, m, n, Arguments[A_jmnNArg(j,m,n)]);
					if (A_jmnVF(j, m, n)>=10)
	 					Set_A_jmn(j, m, n, ODFArray[A_jmnVF(j, m, n)-10].A_jmn(j, m, n));
				}

// Calculation of Negative Volume
  // set up the integration grid and accumulators for the negative-volume estimate
  int NNET=50;		// integration net for ODF
  double integral=0, integral_abs=0;
  // integrate the ODF and its absolute value over a uniform grid in ksi, cos(teta) and psi
  for(int iksi=0; iksi<NNET; iksi++)
	  for(int iteta=0; iteta<NNET; iteta++)
		  for(int ipsi=0; ipsi<NNET; ipsi++)
	  {
	  // Define integration angles and Euler angles Field->Molecule:
		  double ksi=(double)iksi*2.0*PI/(double)NNET;
		  double costeta=-1.0+(double)iteta*2.0/(double)NNET;
		  double psi=(double)ipsi*2.0*PI/(double)NNET;
		  double temp=ODF_value(ksi, acos(costeta), psi);
		  integral+=temp;
		  integral_abs+=fabs(temp);
	  }
		  // (integral_abs-integral)/2 is the integral of the negative part of the ODF; it is divided by the integral of the ODF
		  // for NegativeVolume and by the integral of |ODF| for NegativeVolumeODF3
		  NegativeVolume=(integral_abs-integral)/(2.0*integral);		// Note that this is NOT as in program ODF3!
		  NegativeVolumeODF3=(integral_abs-integral)/(2.0*integral_abs);		// For comparison with ODF3

		  // the negative volume is undefined if an integral of the ODF over orientations vanishes: reset both values
		  // (the trailing comment below says spectrum, but these are integrals of the ODF)
		  if ((integral==0)||(integral_abs==0))
		  {							// Negative volume not applicable: integral over the spectrum is zero.
			  NegativeVolume=0;		// Note that this is NOT as in program ODF3!
			  NegativeVolumeODF3=0;		// For comparison with ODF3
		  }
	}



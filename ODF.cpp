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
#include <math.h>
#include <cmath>
#include <complex>
#include <stdlib.h>

// factorial function
long int fact(int n)
{
	if(n==0) return 1;
	else if (n>=0) return n*fact(n-1);
	else return -1;
}

long double dfact(int n)
{
	if(n==0) return 1;
	else if (n>=0) return n*dfact(n-1);
	else return -1;
}

// small d-Wigner function
double wigner_d(int j, int m, int n, double beta)
{
	int NMIN, NMAX, k;
	double sum=0;
	if (abs(m)>j) return 0;
	if (abs(n)>j) return 0;

	// The choice of summation bounds
	if (m+n>0) NMIN=m+n;
	else       NMIN=0;

	if (j+m<j+n) NMAX=j+m;
	else         NMAX=j+n;

	for(k=NMIN; k<=NMAX; k++)
	{
        if (k%2)				// k nechet
			sum-=pow(cos(beta/2), 2*k-m-n)*pow(sin(beta/2), 2*j+m+n-2*k)/(double)(dfact(k)*dfact(j+m-k)*dfact(j+n-k)*dfact(k-m-n));
		else					// k chet
			sum+=pow(cos(beta/2), 2*k-m-n)*pow(sin(beta/2), 2*j+m+n-2*k)/(double)(dfact(k)*dfact(j+m-k)*dfact(j+n-k)*dfact(k-m-n));
	}

	if((j+m)%2)                   // j+m nechet
		sum*=-sqrt((double)dfact(j+m)*dfact(j-m)*dfact(j+n)*dfact(j-n));
	else                        // j+m chet
		sum*= sqrt((double)dfact(j+m)*dfact(j-m)*dfact(j+n)*dfact(j-n));

	return sum;
}

// Complex Wigner D-function
std::complex<double> ComplexWigner(int j, int m, int n, double alpha, double beta, double gamma)
{
	std::complex<double> I(0, 1);
	std::complex<double> Value;
	return wigner_d(j, m, n, beta)*exp(-gamma*I*(double)n-alpha*I*(double)m);
}

//////////////////////////////////////////////////////////////////////////
////////  CLASS FOR ODF AND ORDER PARAMETER MANIPULATION /////////////////
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
	OrderParameters::OrderParameters()
	{
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
		Jmax=0;
		NOrientations=1;
	}
//////////////////////////////////////////////////////////////////////////
OrderParameters::~OrderParameters()
	{}

//////////////////////////////////////////////////////////////////////////
	void OrderParameters::AssignTjpqi()
	{
		for (int j=0; j<=Jmax; j+=2)
			for(int p=-j; p<=j; p++)
				for(int q=-j; q<=j; q++)
					for(int i=0; i<NOrientations; i++)
						Tjpqi_container[j][p+Jmax][q+Jmax][i]=Ajmn_container[j][p+Jmax][q+Jmax]*ComplexWigner(j, p, 0, alpha_SF[i], beta_SF[i], gamma_SF[i]);
	}
//////////////////////////////////////////////////////////////////////////
	std::complex<double> OrderParameters::A_jmn(int j, int m, int n)
	{
		return Ajmn_container[j][m+Jmax][n+Jmax];
	}
	int OrderParameters::A_jmnVF(int j, int m, int n)
	{
		return Ajmn_containerVF[j][m+Jmax][n+Jmax];
	}
	int OrderParameters::A_jmnNArg(int j, int m, int n)
	{
		return Ajmn_containerNArg[j][m+Jmax][n+Jmax];
	}
//////////////////////////////////////////////////////////////////////////
	std::complex<double> OrderParameters::T_jpqi(int j, int p, int q, int i)
	{
		return Tjpqi_container[j][p+Jmax][q+Jmax][i];
	}
//////////////////////////////////////////////////////////////////////////
	double OrderParameters::ODF_value(double ksi, double teta, double psi)
	{
		std::complex<double> sum(0,0);
		for(int j=0; j<=Jmax; j+=2)
			for(int m=-j; m<=j; m+=2)
				for(int n=-j; n<=j; n+=2)
				{
					sum+=A_jmn(j, m, n)*ComplexWigner(j, m, n, ksi, teta, psi);
				}
		return real(sum);

	}
	void OrderParameters::ODF_reduced(double teta, double ksi)
// Given the orientation distribution function A_jmn,
// and the set of angles _SF for all sample orientations,
// this subroutine calculates the number density of
// particles which are oriented in resonator by angles ksi, teta.
// The calculated values are stored in array CurrentODF[0..NOrientations],
// to be used in the block of spectra calculation.
	{
		for (int i=0; i<NOrientations; i++)
			CurrentODF[i]=0;

		for(int j=0; j<=Jmax; j+=2)
			for(int p=-j; p<=j; p++)
				for(int q=-j; q<=j; q++)
					for (int i=0; i<NOrientations; i++)
						CurrentODF[i]+=real(T_jpqi(j, p, q, i)*ComplexWigner(j, 0, q, 0, teta, ksi));
		return;
	}
//////////////////////////////////////////////////////////////////////////
	void OrderParameters::Set_A_jmn(int j, int m, int n, std::complex<double> value)
// Sets the value of specified order parameter
	{
		Ajmn_container[j][m+Jmax][n+Jmax]=value;
	}
	void OrderParameters::Set_A_jmnVF(int j, int m, int n, int value)
// Sets the value of specified order parameter variate flag
	{
		Ajmn_containerVF[j][m+Jmax][n+Jmax]=value;
	}
	void OrderParameters::Set_A_jmnNArg(int j, int m, int n, int value)
// Sets the value of specified order parameter number in Arg array
	{
		Ajmn_containerNArg[j][m+Jmax][n+Jmax]=value;
	}
//////////////////////////////////////////////////////////////////////////
	void OrderParameters::Realize()
// This function makes the ODF A_jmn real by formula (3.5-6) in dissertation
	{
		for (int j=0; j<=Jmax; j++)
		{
			for(int m=-j; m<=0; m++)
// note that in fact we do not have non-zero coefficients A_jmn with m and n of different sign, 
// so we can make the next line: for(int n=-j; n<=0; n++)
				for(int n=-j; n<=j; n++)
					Ajmn_container[j][m+Jmax][n+Jmax]=conj(Ajmn_container[j][-m+Jmax][-n+Jmax]);
			Ajmn_container[j][0+Jmax][0+Jmax]=real(Ajmn_container[j][0+Jmax][0+Jmax]);
		}
	}

	void OrderParameters::ReadFromFile(char* filename)
	{
	ifstream ODF_par;
	ODF_par.open (filename, ios::in);
	if(!ODF_par.is_open())
	{
		std::cout << "Could not open file " << filename << std::endl;
		exit(0);
	}
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
		if (!strcmp(buffer2, "Penalty:"))
			sscanf(buffer, "%s %lf", buffer3, &Penalty);
		if (!strcmp(buffer2, "DelHi0:"))
			sscanf(buffer, "%s %lf %d", buffer3, &DelHi0, &DelHi0VF);
		if (!strcmp(buffer2, "DelHi1:"))
			sscanf(buffer, "%s %lf %d", buffer3, &DelHi1, &DelHi1VF);
		if (!strcmp(buffer2, "DelHi2:"))
			sscanf(buffer, "%s %lf %d", buffer3, &DelHi2, &DelHi2VF);
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
		if (!strcmp(buffer2, "Particle:"))
		{
			int npart;
			sscanf(buffer, "%s %d", buffer3, &npart);
			if(npart!=Count) break;
		}

	}


	}
	void OrderParameters::ConvertAnglesToDegrees(void)
	{
		for(int i=0; i<MAX_ORIENTATIONS; i++)
		{
				alpha_SF[i]*=180/PI;
				beta_SF[i]*=180/PI;
				gamma_SF[i]*=180/PI;
		}
	}
	void OrderParameters::ConvertAnglesToRadians(void)
	{
		for(int i=0; i<MAX_ORIENTATIONS; i++)
		{
				alpha_SF[i]/=180/PI;
				beta_SF[i]/=180/PI;
				gamma_SF[i]/=180/PI;
		}
	}
	void OrderParameters::SetEuler_SF(void)
	{
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
	void OrderParameters::WriteToFile(char* filename)
	{
		FILE *file;
		file=fopen(filename, "a+");
		fprintf(file, "\nParticle: %d\n", Count);
		fprintf(file, "ODFRank:       %d\n", Jmax);
		fprintf(file, "Penalty:       %lf      NegativeVolume: %e  (ODF3: %e)\n", Penalty, NegativeVolume, NegativeVolumeODF3);
		fprintf(file, "DelHi0:        %3.5lf  %d\n", DelHi0, DelHi0VF);
		fprintf(file, "DelHi1:        %3.5lf  %d\n", DelHi1, DelHi1VF);
		fprintf(file, "DelHi2:        %3.5lf  %d\n", DelHi2, DelHi2VF);
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

	void OrderParameters::WriteToFileBrief(char* filename)
	{
		FILE *file;
		file=fopen(filename, "a+");
		fprintf(file, "ODFRank: %d, ", Jmax);
		fprintf(file, "Penalty:  %lf,  NegativeVolume: %e  (ODF3: %e)\n", Penalty, NegativeVolume, NegativeVolumeODF3);
		fprintf(file, "DelHi0,1,2:        %3.5lf  %d  %3.5lf  %d  %3.5lf  %d\n", DelHi0, DelHi0VF, DelHi1, DelHi1VF, DelHi2, DelHi2VF);
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
	
	void OrderParameters::Par2Arg(double* Arguments, double *ub, double *lb, int *NArg)
	{
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
	void OrderParameters::Arg2Par(const double* Arguments)
	{
		if(DelHi0NArg!=-1) DelHi0=Arguments[DelHi0NArg];
		if(DelHi1NArg!=-1) DelHi1=Arguments[DelHi1NArg];
		if(DelHi2NArg!=-1) DelHi2=Arguments[DelHi2NArg];

		if (DelHi0VF>=10)
	 		DelHi0=ODFArray[DelHi0VF-10].DelHi0;
		if (DelHi1VF>=10)
	 		DelHi1=ODFArray[DelHi1VF-10].DelHi1;
		if (DelHi2VF>=10)
	 		DelHi2=ODFArray[DelHi2VF-10].DelHi2;

		for (int j=0; j<=Jmax; j+=2)
			for(int m=0; m<=j; m+=2)
				for(int n=0; n<=j; n+=2)
				{
					if(A_jmnNArg(j,m,n)!=-1) Set_A_jmn(j, m, n, Arguments[A_jmnNArg(j,m,n)]);
					if (A_jmnVF(j, m, n)>=10)
	 					Set_A_jmn(j, m, n, ODFArray[A_jmnVF(j, m, n)-10].A_jmn(j, m, n));
				}

// Calculation of Negative Volume
  int NNET=50;		// integration net for ODF
  double integral=0, integral_abs=0;
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
		  NegativeVolume=(integral_abs-integral)/(2.0*integral);		// Note that this is NOT as in program ODF3!
		  NegativeVolumeODF3=(integral_abs-integral)/(2.0*integral_abs);		// For comparison with ODF3

		  if ((integral==0)||(integral_abs==0))
		  {							// Negative volume not applicable: integral over the spectrum is zero.
			  NegativeVolume=0;		// Note that this is NOT as in program ODF3!
			  NegativeVolumeODF3=0;		// For comparison with ODF3
		  }
	}



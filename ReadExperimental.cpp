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
#include <fstream>
#include <iostream>
#include <stdlib.h>

char Output_Modifier[MAX_STRLEN];


// Reading of experimental EPR spectrum files.
// The files may contain any information. At some place in the file, there should be the line "Spectr:"
// followed by the spectrum data, as two columns (Field Intensity) separated with spaces
extern int MpiRank, MpiSize;
extern double cov[MAX_ARG][MAX_ARG];		// covariance matrix

void Experimental::LoadFromFile(char * filename)
{
	char buffer[MAX_STRLEN];
	ifstream file, file2;
	file.open (filename, ios::in);
	if(!file.is_open())
	{
		std::cout << "Could not open file " << filename << std::endl;
		exit(0);
	}
	NPoints=0;
	int Nlines=0, addlines=1;
// How Many Points are there in spectrum?
	file.getline(buffer, MAX_STRLEN);
	while (!file.eof())
	{
		Nlines+=addlines;
		NPoints++;
		if(!strcmp(buffer, "Spectr:") || !strcmp(buffer, "Spectr:\r")) 
		{
			NPoints=0;
			addlines=0;
		}
		file.getline(buffer, MAX_STRLEN);
	}
	file.close();

	NPoints--;
	SpectrumX=new double[NPoints];
	SpectrumY=new double[NPoints];

	file2.open(filename, ios::in);
	if(!file2.is_open())
	{
		std::cout << "Could not open file " << filename << std::endl;
		exit(0);
	}
	for(int i=0; i<Nlines; i++)
	{
		file2.getline(buffer, MAX_STRLEN);
	}
	for(int i=0; i<NPoints; i++)
	{
		file2.getline(buffer, MAX_STRLEN);
		file2 >> SpectrumX[i];
		file2 >> SpectrumY[i];
//		std::cout << buffer << std::endl;
//		std::cout << SpectrumX[i] << "   "<< SpectrumY[i] << std::endl;
	}
	file.close();

}

double Experimental::Interpolate(double X)
{
	int k=0;
	if (X<SpectrumX[0]) return 0;
	else if (X>SpectrumX[NPoints-1]) return 0;
	else
		for (k=0; SpectrumX[k+1] < X; k++);
	if (SpectrumX[k+1]==X) return SpectrumY[k+1];
	return SpectrumY[k]+(X-SpectrumX[k])*(SpectrumY[k+1]-SpectrumY[k])/(SpectrumX[k+1]-SpectrumX[k]);
}

void Experimental::InterpolateWhole(double *FieldNet, int Npts)
{
	SpectrumInterpolated=new double[Npts];
	for(int i=0; i<Npts; i++)
		SpectrumInterpolated[i]=Interpolate(FieldNet[i]);


}

Experimental::Experimental()
{}
Experimental::~Experimental()
{}


void Experimental::CalculateStats(double Hmax)
{
	MeanL=0;
	VarianceL=0;
	int Nstat=0;

	for(int i=0; i<NPoints; i++)
		if(SpectrumX[i]<=Hmax)
		{
			Nstat++;
			MeanL+=SpectrumY[i];
		}
	MeanL/=Nstat;
	for(int i=0; i<NPoints; i++)
		if(SpectrumX[i]<=Hmax)
		{
			VarianceL+=pow(SpectrumY[i]-MeanL, 2);
		}
	VarianceL/=(Nstat-1);

	return;



}

void WriteOutput(int NParticles,
				 Experimental *ExpSpectra, 
				 double ***SpectrumSum, 
				 MagneticParameters *MagneticParameter, 
				 OrderParameters *ODF,
				 int NCovarianceArguments,
				 const double * Arguments)
{
	FILE *output;
	char buffer[MAX_STRLEN];

	for(int i=0; i<ODF->NOrientations; i++)
	{
		sprintf(buffer, "%s%s.dat", Output_Modifier, ODF->ExperimentalFiles[i]);
		output=fopen( buffer, "w");
		fprintf(output, "Discrepancy: %e  , Per Point:  %e ,", MagneticParameter[0].CurrentDiscrepancyOverall, MagneticParameter[0].CurrentDiscrepancyPerPoint);
// Calculating heuristic measure of experimental noise.
// First find the left wing of the spectrum, where the spectrum intensity is less than 10^-5 x spectrum amplitude
// find the spec amplitude
		double SMax=0;
		double SMin=0;
		double Hmax;
		for (int j=0; j<MagneticParameter[0].NPoints; j++)
		{
			double summ=0;
			for (int k=0; k<NParticles; k++)
				summ+=SpectrumSum[k][i][j];
			if(SMax<summ) SMax=summ;
			if(SMin>summ) SMin=summ;
		}

		for (int j=0; j<MagneticParameter[0].NPoints; j++)
		{
			double summ=0;
			for (int k=0; k<NParticles; k++)
				summ+=SpectrumSum[k][i][j];

			if (fabs(summ)>1.0e-4*(SMax-SMin))
			{
				Hmax=MagneticParameter[0].FieldNet[j];
				break;
			} 
		}
		ExpSpectra[i].CalculateStats(Hmax);
		fprintf(output, "0.5*Exp.sp.variance: %e (mean %e) , calc. at fields < %7.2lf Gauss \n", 0.5*ExpSpectra[i].VarianceL, ExpSpectra[i].MeanL, Hmax);
	

		fclose(output);


		for (int j=0; j<NParticles; j++)
		{
			MagneticParameter[j].WriteToFileBrief(buffer);
			ODF[j].WriteToFileBrief(buffer);
		}
		output=fopen( buffer, "a+");
		if (NCovarianceArguments)
		{
			fprintf(output, "\nArgument errors:\n");
			for(int j=0; j<NCovarianceArguments; j++)
				fprintf(output, "Arg[%d]=%lf +- %lf \n", j, Arguments[j], sqrt(cov[j][j]));
			fprintf(output, "Covariance matrix:\n");
			for(int k=0; k<NCovarianceArguments; k++)
			{
				fprintf (output, "c");
				for(int j=0; j<=k; j++)
					fprintf(output, "%lf \t", cov[k][j]/sqrt(cov[k][k]*cov[j][j]));
				fprintf (output, "\n");
	}

		}
		fprintf(output, "\n   H_exp       S_exp             H        Exp_interp       ");
		for(int k=0; k<NParticles; k++)
			fprintf(output,"Part%d          ", k);
		fprintf(output, "Sum\n");
		for (int j=0; j<max(MagneticParameter->NPoints, ExpSpectra[i].NPoints); j++)
		{
			if (j<ExpSpectra[i].NPoints)
				fprintf(output, "  %5.3lf   %+3.6lf      ", ExpSpectra[i].SpectrumX[j], 
															ExpSpectra[i].SpectrumY[j]);
			else
				fprintf(output, "      ,        ,        ");

			if(j<MagneticParameter[0].NPoints)
					fprintf(output, "  %5.3lf    %+5.6lf   ", MagneticParameter[0].FieldNet[j], ExpSpectra[i].SpectrumInterpolated[j]);
			else
						fprintf(output, "      ,      ,     ");
			double sum=0;
			for(int k=0; k<NParticles; k++)
			{
				if(j<MagneticParameter[k].NPoints)
				{
						sum+=SpectrumSum[k][i][j];
						fprintf(output, "   %+3.6lf   ", SpectrumSum[k][i][j]);
				}
				else
						fprintf(output, "      ,        ");
			}
				if(j<MagneticParameter[0].NPoints)
						fprintf(output, "   %+3.6lf   ", sum);
				else
						fprintf(output, "      ,         ");

			fprintf(output, "\n");
																			 
		}
	fclose(output);
	}
}

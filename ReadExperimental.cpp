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

// Experimental EPR spectra. Class Experimental (declared in BiradicalODF.h) holds one measured spectrum as columns
// of field and intensity, loads it from a file, interpolates it onto the field net of the simulation and estimates
// its noise level. WriteOutput writes the result files, with the experimental, interpolated and simulated spectra
// and the parameters of all particles. The file also defines the global string Output_Modifier (prefix of the
// output files), which is set by ReadName in ReadNamenam.cpp.
// Used by: main.cpp (loading and interpolation of the spectra, final call of WriteOutput) and Discrepancy.cpp
// (Discrepancy compares the interpolated spectra with the simulation and calls WriteOutput).

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

// Reads an experimental spectrum from a text file. Lines up to and including the line Spectr: are header, after it
// come two columns (field, intensity); fills SpectrumX, SpectrumY and NPoints. The line directly below Spectr: is
// consumed and discarded by the first getline of the reading loop, so it acts as a column heading.
// Input: filename, the experimental file. Output: the arrays SpectrumX and SpectrumY of this object.
// The program exits if the file cannot be opened.
// Used by: main() in main.cpp, once for each spectrum of the set.
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
	// First pass: count the header lines (up to and including Spectr:) and the lines after it
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

	// The number of points is one less than the number of lines counted after Spectr:; allocate the arrays
	NPoints--;
	SpectrumX=new double[NPoints];
	SpectrumY=new double[NPoints];

	// Second pass: skip the header, then read the (field, intensity) pairs
	file2.open(filename, ios::in);
	if(!file2.is_open())
	{
		std::cout << "Could not open file " << filename << std::endl;
		exit(0);
	}
	// Skip the header lines, including the line Spectr:
	for(int i=0; i<Nlines; i++)
	{
		file2.getline(buffer, MAX_STRLEN);
	}
	// Read one pair per point; the getline at the start of each pass discards the rest of the previous line
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

// Linear interpolation of the experimental spectrum at the field X; zero outside the measured field range.
// SpectrumX must be in increasing order. Input: X, a field value in the units of SpectrumX.
// Output: the interpolated intensity.
// Used by: Experimental::InterpolateWhole (this file).
double Experimental::Interpolate(double X)
{
	int k=0;
	// Outside the measured range the spectrum is zero; inside, search upwards for the interval that contains X
	if (X<SpectrumX[0]) return 0;
	else if (X>SpectrumX[NPoints-1]) return 0;
	else
		for (k=0; SpectrumX[k+1] < X; k++);
	if (SpectrumX[k+1]==X) return SpectrumY[k+1];
	// Linear interpolation between the two neighbouring experimental points
	return SpectrumY[k]+(X-SpectrumX[k])*(SpectrumY[k+1]-SpectrumY[k])/(SpectrumX[k+1]-SpectrumX[k]);
}

// Evaluates the experimental spectrum at every field of the simulation net and stores the result in
// SpectrumInterpolated, so that experimental and simulated spectra can be compared point by point.
// Inputs: FieldNet, the Npts field values of the simulation. Output: the array SpectrumInterpolated.
// Used by: main() in main.cpp, once for each spectrum of the set (the result is read by Discrepancy).
void Experimental::InterpolateWhole(double *FieldNet, int Npts)
{
	// Allocate the array and interpolate the experimental spectrum at each field of the net
	SpectrumInterpolated=new double[Npts];
	for(int i=0; i<Npts; i++)
		SpectrumInterpolated[i]=Interpolate(FieldNet[i]);


}

// Constructor: does nothing, the spectrum is filled by LoadFromFile.
// Used by: implicitly, when main() creates the array new Experimental[NOrientations] (main.cpp).
Experimental::Experimental()
{}
// Destructor: does nothing (the arrays allocated by LoadFromFile and InterpolateWhole are not freed).
// Used by: nothing in this code base (the array of Experimental objects in main.cpp is never deleted).
Experimental::~Experimental()
{}


// Calculates the mean MeanL and the variance VarianceL (sum of squared deviations divided by the number of points
// minus one) of the experimental intensity over the points with field not above Hmax. WriteOutput applies it to the
// low-field wing, where the simulated signal is negligible, as a heuristic measure of the experimental noise.
// Input: Hmax, the highest field included, in the units of SpectrumX. Outputs: MeanL and VarianceL of this object.
// Used by: WriteOutput (this file).
void Experimental::CalculateStats(double Hmax)
{
	// Mean intensity over the points with field not above Hmax
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
	// Sum of squared deviations from the mean over the same points, divided by Nstat-1 below
	for(int i=0; i<NPoints; i++)
		if(SpectrumX[i]<=Hmax)
		{
			VarianceL+=pow(SpectrumY[i]-MeanL, 2);
		}
	VarianceL/=(Nstat-1);

	return;



}

// Writes one result file for each spectrum of the set. The file is named from Output_Modifier and the name of the
// experimental file, with the extension .dat, and holds the discrepancy, the noise estimate of the experimental spectrum,
// the brief parameters of every particle and of its ODF, optionally the errors of the fitted arguments with their
// correlation matrix, and a table with the experimental spectrum, its interpolation onto the simulation field net,
// the simulated spectrum of each particle and their sum.
// Inputs: NParticles; ExpSpectra, the experimental spectra; SpectrumSum[particle][spectrum][field point], the simulated
// spectra; MagneticParameter and ODF, the parameters of all particles; NCovarianceArguments, the number of fitted
// arguments with errors available (0 for none); Arguments, the values of the fitted arguments.
// Used by: Discrepancy() in Discrepancy.cpp (after a calculation) and main() in main.cpp (at the end of the run,
// with the covariance information of the fit when it has been calculated).
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

	// One output file per spectrum of the set
	for(int i=0; i<ODF->NOrientations; i++)
	{
		// Name the file from the output prefix and the experimental file name, and write the discrepancy first
		sprintf(buffer, "%s%s.dat", Output_Modifier, ODF->ExperimentalFiles[i]);
		output=fopen( buffer, "w");
		fprintf(output, "Discrepancy: %e  , Per Point:  %e ,", MagneticParameter[0].CurrentDiscrepancyOverall, MagneticParameter[0].CurrentDiscrepancyPerPoint);
// Calculating heuristic measure of experimental noise.
// First find the left wing of the spectrum, where the spectrum intensity is less than 10^-5 x spectrum amplitude
// find the spec amplitude
		double SMax=0;
		double SMin=0;
		double Hmax;
		// Largest and smallest value of the summed simulated spectrum
		for (int j=0; j<MagneticParameter[0].NPoints; j++)
		{
			double summ=0;
			for (int k=0; k<NParticles; k++)
				summ+=SpectrumSum[k][i][j];
			if(SMax<summ) SMax=summ;
			if(SMin>summ) SMin=summ;
		}

		// Hmax is the first field where the summed spectrum exceeds 1e-4 of its range
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
		// Noise estimate from the experimental points at fields below Hmax
		ExpSpectra[i].CalculateStats(Hmax);
		fprintf(output, "0.5*Exp.sp.variance: %e (mean %e) , calc. at fields < %7.2lf Gauss \n", 0.5*ExpSpectra[i].VarianceL, ExpSpectra[i].MeanL, Hmax);
	

		fclose(output);


		// Append the brief parameter listing of every particle and of its ODF (these methods open the file by name)
		for (int j=0; j<NParticles; j++)
		{
			MagneticParameter[j].WriteToFileBrief(buffer);
			ODF[j].WriteToFileBrief(buffer);
		}
		// Reopen the file for appending: fit errors and correlations, then the table of spectra
		output=fopen( buffer, "a+");
		// Standard errors of the fitted arguments (square roots of the covariance diagonal) and their correlation matrix
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
		// Table header, then one row per field point: experimental data, interpolated experimental spectrum, simulated
		// spectrum of each particle and their sum; columns that have ended are filled with commas
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

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

// Reading of the control file name.nam, which gives the names of the input and output parameter files,
// the set of experimental spectra with their sample orientations, and whether the
// parameters are to be fitted (Minimization: yes) or a single simulation is wanted.
// Used by: main.cpp (main calls ReadName once, before any parameter file is read). ReadName also sets the global
// flag MinimizationFlag (defined in main.cpp, read by main.cpp and BiradicalTransitions.cpp) and the global
// string Output_Modifier (defined in ReadExperimental.cpp, read by WriteOutput there).

#include "BiradicalODF.h"
#include <string.h>

// Reading file "name.nam" containing basic information about simulation and simulation files.
extern int MinimizationFlag;
extern char Output_Modifier[MAX_STRLEN];

// Reads name.nam from the working directory and fills in the file names and the list of spectra.
// Keywords (first word of a line): NSpectra (number of spectra), Minimization (yes or no), InputODF, OutputODF,
// InputMagnetic and OutputMagnetic (file names), OutputPrefix (prefix of the output spectrum files, default out),
// and Spectra, followed by one line per spectrum with four fields: experimental file name, sample rotation angle
// (hi_SF), rotation axis (axis_SF: 0=x, 1=y, 2=z) and BasisSource. BasisSource is the word truncation or an empty
// name for a basis built from the truncation
// limits, otherwise it is handed to LoadBasisB to load a stored basis. NSpectra must precede Spectra.
// Inputs: ODF receives NOrientations and the per-spectrum arrays; the four file-name buffers are filled in;
// MagneticParameter is not used. Outputs: sets MinimizationFlag and Output_Modifier, always returns 0;
// the program exits if name.nam cannot be opened; the fields of a Spectra line are not checked.
// Used by: main() in main.cpp, called once at start-up (on ODF[0] and MagneticParameter[0]).
int ReadName(MagneticParameters *MagneticParameter, OrderParameters *ODF, char* Magnetic_in_filename, char *ODF_in_filename, char* Magnetic_out_filename, char *ODF_out_filename)
{
	// Defaults: no minimisation, output prefix out
	MinimizationFlag=0;
	strcpy(Output_Modifier, "out");
	// Open name.nam in the working directory; stop if it is missing
	ifstream magnetic_par;
	magnetic_par.open ("name.nam", ios::in);
	if(!magnetic_par.is_open())
	{
		std::cout << "Could not open file name.nam\n";
		exit(0);
	}
	// Read name.nam line by line: the first word of each line selects the setting it carries
	while (!magnetic_par.eof())
	{
		char buffer[MAX_STRLEN], buffer2[MAX_STRLEN], buffer3[MAX_STRLEN];
		magnetic_par.getline(buffer, MAX_STRLEN);
		if( sscanf(buffer, "%s", buffer2)==-1) continue;
		// Number of spectra (sample orientations) in the set
		if (!strcmp(buffer2, "NSpectra:"))
			sscanf(buffer, "%s %d", buffer3, &(ODF->NOrientations));
		// Minimization: yes switches the least-squares fit on
		if (!strcmp(buffer2, "Minimization:"))
		{
			char buffer4[MAX_STRLEN];
			sscanf(buffer, "%s %s", buffer3, buffer4);
			if(!strcmp(buffer4, "yes"))
				MinimizationFlag=1;
//			else if(!strcmp(buffer, "transitions"))
//				MinimizationFlag=2;
		}
		// Names of the input and output ODF and magnetic parameter files, and the prefix of the output files
		if (!strcmp(buffer2, "InputODF:"))
			sscanf(buffer, "%s %s", buffer3, ODF_in_filename);
		if (!strcmp(buffer2, "OutputPrefix:"))
			sscanf(buffer, "%s %s", buffer3, Output_Modifier);
		if (!strcmp(buffer2, "OutputODF:"))
			sscanf(buffer, "%s %s", buffer3, ODF_out_filename);
		if (!strcmp(buffer2, "InputMagnetic:"))
			sscanf(buffer, "%s %s", buffer3, Magnetic_in_filename);
		if (!strcmp(buffer2, "OutputMagnetic:"))
			sscanf(buffer, "%s %s", buffer3, Magnetic_out_filename);
		// Spectra table: one line per spectrum, NOrientations lines in all
		if (!strcmp(buffer2, "Spectra:"))
		{
			for(int i=0; i<ODF->NOrientations; i++)
			{
//				char buffer1[30];
				magnetic_par.getline(buffer, MAX_STRLEN);
				sscanf(buffer, "%s %lf %d %s", ODF->ExperimentalFiles[i], &(ODF->hi_SF[i]), &(ODF->axis_SF[i]), (ODF->BasisSource[i]) );
//				ODF->alpha_SF0[i] *= PI/180;
//				ODF->beta_SF0[i]  *= PI/180;
//				ODF->gamma_SF0[i] *= PI/180;
			}
		}
		
	}
	return 0;

}
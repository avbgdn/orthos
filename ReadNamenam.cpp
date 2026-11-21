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
#include <string.h>

// Reading file "name.nam" containing basic information about simulation and simulation files.
extern int MinimizationFlag;
extern char Output_Modifier[MAX_STRLEN];

int ReadName(MagneticParameters *MagneticParameter, OrderParameters *ODF, char* Magnetic_in_filename, char *ODF_in_filename, char* Magnetic_out_filename, char *ODF_out_filename)
{
	MinimizationFlag=0;
	strcpy(Output_Modifier, "out");
	ifstream magnetic_par;
	magnetic_par.open ("name.nam", ios::in);
	if(!magnetic_par.is_open())
	{
		std::cout << "Could not open file name.nam\n";
		exit(0);
	}
	while (!magnetic_par.eof())
	{
		char buffer[MAX_STRLEN], buffer2[MAX_STRLEN], buffer3[MAX_STRLEN];
		magnetic_par.getline(buffer, MAX_STRLEN);
		if( sscanf(buffer, "%s", buffer2)==-1) continue;
		if (!strcmp(buffer2, "NSpectra:"))
			sscanf(buffer, "%s %d", buffer3, &(ODF->NOrientations));
		if (!strcmp(buffer2, "Minimization:"))
		{
			char buffer4[MAX_STRLEN];
			sscanf(buffer, "%s %s", buffer3, buffer4);
			if(!strcmp(buffer4, "yes"))
				MinimizationFlag=1;
//			else if(!strcmp(buffer, "transitions"))
//				MinimizationFlag=2;
		}
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
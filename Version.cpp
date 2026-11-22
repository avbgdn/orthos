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

// Program version of Orthos. Version() records the current version string in a small text file named
// version in the working directory; the history of changes is kept as comments inside the function.
// Used by: main.cpp (main calls Version once at start-up).

#include "BiradicalODF.h"

// Writes the program version string (Orthos version and its date) to the file version in the working
// directory; the comment block inside the function is the change log of the program.
// No inputs or outputs.
// Used by: main() in main.cpp, called once at start-up (by every process when MPI is used).
void Version(void)
{
////////////////////////////	Program version		/////////////////////////////
//  13/08/2015 - basic version by AB
//  14/08/2015 - forced equal fs and NPoints for all particles
//				 discrepancy file, which is updated with current discrepancy
//				 added sign of J0, J0 is possible to variate as logarithm.
//  17/08/2015 - fixed the bug with interpreting void lines in input files.
//  11/09/2015 - fixed output of variated order parameters (real-complex error stuff)
//               made the axes of D-tensor coaxial with g1 tensor
//  14/09/2015 - angular dependence of Gaussian linewidth for liouvillian calculation
//  15/09/2015 - for zeroth particle independent pruned basis sets can be used for different orientations
//  21/09/2015 - Linux-compatible version
//  16/10/2015 - MPI implemented
//  02/11/2015 - Fixed bug with the starting vector.
//  06/11/2015 - Separate output for different particle numbers
//  30/12/2015 - Updated NL2SOL minimizer to 2015 version
//  16/01/2016 - Included Discrete Jump model (according to Freed), Monoradicals are calculated in M-symmetrized basis
//  19/01/2016 - Debugged function ODF::Realize()
//  15/04/2016 - Added and debugged penalty for negative volume
//  31/10/2016 - Fixed memory leak in NL2SOL invocation
//				 Added interface for my Jump Diffusion model
//				 Enabled Gauss-Kronrod integration
//
//				 Added librations
//				 Started to store XLK coefficients for speed.
//  19/12/2016 - Model of jump reorientation
//               improved parallel version
//
//  08/04/2017 - Orthos 1.3 version, improvements:
//               - Gauss-Kronrod integration on arbitrary net size
//               - Isotropic hyperfine splitting
//               - More detailed output information in output files (not only spectra)
//               - Output Prefix in name.nam file
//  15/10/2018   - Non-axial ZFS tensor
//
//  31/01/2019   - non-Brownian diffusion types
//
//  23/04/2020   Orthos 2.0 version, improvements
//               - added conformational mobility of biradicals.
//  03/05/2023   - fixed bugs with isotropic hyperfine
//

	// Version string of the current release
	char ProgramVersion[MAX_STRLEN];
	strcpy(ProgramVersion, "Orthos v2.0, 03/05/2023");
	// Write the version string to the file named version
	std::fstream prog_vers;
	prog_vers.open("version", std::fstream::out);
	prog_vers << ProgramVersion << endl;
	prog_vers.close();

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

}
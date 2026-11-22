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

// File purpose: generates the nodes and weights of the quadrature nets over the two orientation angles of the
// molecule relative to the field: cos(teta) in [-1,1] (Node, KronWeight, GausWeight) and the azimuth ksi
// (NodeKsi, KronWeightKsi, GausWeightKsi, mapped by the user to ksi=PI*(NodeKsi+1) in [0,2*PI]). Both are Gauss-Kronrod
// rules, whose nodes and weights are computed by kronrod() in kronrod.cpp.
// Used by: main.cpp calls GenerateIntegrationNet() and GenerateIntegrationNetKsi() for each particle;
// the generated arrays are read by CalculateIntegratedSpectrumLiouville() in IntegratedSpectrum.cpp.

#include "BiradicalODF.h"
#define USE_GAUSS_KRONROD 2
#include <stdio.h>
#include "kronrod.h"

// GenerateIntegrationNet: fills the nodes Node[] on [-1,1] (cos(teta)) and the Kronrod and Gauss weights KronWeight[]
// and GausWeight[] for IntegrationNet nodes. The weights are normalised to add up to 1, so the quadrature gives an
// average over the interval. IntegrationNet itself is adjusted here: cut to MAX_GKNODES, made odd and at least 3
// (unless it is 1, which is a single node at 1 with unit weight; GaussNodes and KronrodNodes are set accordingly).
// Used by: main() in main.cpp, once for each particle.
void MagneticParameters::GenerateIntegrationNet()
{
// a single orientation: one node at cos(teta)=1 with unit weights
if (IntegrationNet==1)
{
	GaussNodes=1; KronrodNodes=1;
	Node[ 0]= 1.0000e+00; KronWeight[ 0]=1.0000; GausWeight[ 0]=1.0000;
	return;
}
// Gauss-Kronrod net: adjust the number of nodes, then get the nodes and weights from kronrod()
if (USE_GAUSS_KRONROD)
{
	if (IntegrationNet>MAX_GKNODES) 
		IntegrationNet=MAX_GKNODES;
	if (!(IntegrationNet%2)) IntegrationNet++;		// if integration net is even, make it odd, so the mirrored net is symmetric about a central node
	if (IntegrationNet<3) IntegrationNet=3;			// floor to the smallest net for which kronrod(IntegrationNet/2,...) receives a meaningful order (>=1)
	GaussNodes=IntegrationNet; KronrodNodes=IntegrationNet;
	// nodes and weights of one half of the symmetric rule (Gauss weights are zero at the Kronrod-only nodes)
	kronrod ( IntegrationNet/2, 1.0e-9, Node, KronWeight, GausWeight );

	// fill the remaining entries by mirroring the generated nodes with a changed sign
	for(int i=IntegrationNet/2+1; i<IntegrationNet; i++)
	{
		Node[i]=-Node[i-IntegrationNet/2-1];
		KronWeight[i]=KronWeight[i-IntegrationNet/2-1];
		GausWeight[i]=GausWeight[i-IntegrationNet/2-1];
	}
	// halve the weights so that they add up to 1 instead of 2 over [-1,1]
	for(int i=0; i<IntegrationNet; i++)				// NORMALIZATION
	{
		KronWeight[i]/=2.0;
		GausWeight[i]/=2.0;
	}
/*
	for(int i=0; i<IntegrationNet; i++)
		cout << Node[i] << " \t" << KronWeight[i] << " \t" << GausWeight[i] << endl;
	cout << endl;
*/
}
else  // if (USE_GAUSS_KRONROD==0)
{
	// equidistant nodes with equal weights; not used as long as USE_GAUSS_KRONROD is nonzero (it is defined as 2 at the top of this file)
	GaussNodes=IntegrationNet; KronrodNodes=IntegrationNet;
	for(int i=0; i<IntegrationNet; i++)
	{
		Node[i]=(double)(i*2.0)/(double)IntegrationNet - 1.0;
		KronWeight[i]=1/(double)IntegrationNet;
		GausWeight[i]=1/(double)IntegrationNet;
	}
}

}





// GenerateIntegrationNetKsi: the same as GenerateIntegrationNet() for the azimuth ksi: fills NodeKsi[] on [-1,1],
// KronWeightKsi[] and GausWeightKsi[] for IntegrationNetKsi nodes, with weights adding up to 1, and adjusts
// IntegrationNetKsi in the same way (cut to MAX_GKNODES, made odd and at least 3, unless it is 1).
// Used by: main() in main.cpp, once for each particle (IntegrationNetKsi is 1 for Liouville particles, equal to
// IntegrationNet otherwise).
void MagneticParameters::GenerateIntegrationNetKsi()
{
// a single azimuth: one node at 1 with unit weights
if (IntegrationNetKsi==1)
{
	GaussNodesKsi=1; KronrodNodesKsi=1;
	NodeKsi[ 0]= 1.0000e+00; KronWeightKsi[ 0]=1.0000; GausWeightKsi[ 0]=1.0000;
	return;
}

// Gauss-Kronrod net: adjust the number of nodes, then get the nodes and weights from kronrod()
if (USE_GAUSS_KRONROD)
{
	if (IntegrationNetKsi>MAX_GKNODES) 
		IntegrationNetKsi=MAX_GKNODES;
	if (!(IntegrationNetKsi%2)) IntegrationNetKsi++;		// if integration net is even, make it odd, so the mirrored net is symmetric about a central node
	if (IntegrationNetKsi<3) IntegrationNetKsi=3;			// floor to the smallest net for which kronrod(IntegrationNetKsi/2,...) receives a meaningful order (>=1)
	GaussNodesKsi=IntegrationNetKsi; KronrodNodesKsi=IntegrationNetKsi;
	// nodes and weights of one half of the symmetric rule (Gauss weights are zero at the Kronrod-only nodes)
	kronrod ( IntegrationNetKsi/2, 1.0e-9, NodeKsi, KronWeightKsi, GausWeightKsi );

	// fill the remaining entries by mirroring the generated nodes with a changed sign
	for(int i=IntegrationNetKsi/2+1; i<IntegrationNetKsi; i++)
	{
		NodeKsi[i]=-NodeKsi[i-IntegrationNetKsi/2-1];
		KronWeightKsi[i]=KronWeightKsi[i-IntegrationNetKsi/2-1];
		GausWeightKsi[i]=GausWeightKsi[i-IntegrationNetKsi/2-1];
	}
	// halve the weights so that they add up to 1 instead of 2 over [-1,1]
	for(int i=0; i<IntegrationNetKsi; i++)				// NORMALIZATION
	{
		KronWeightKsi[i]/=2.0;
		GausWeightKsi[i]/=2.0;
	}

}
else  // if (USE_GAUSS_KRONROD==0)
{
	// equidistant nodes with equal weights; not used as long as USE_GAUSS_KRONROD is nonzero (it is defined as 2 at the top of this file)
	GaussNodesKsi=IntegrationNetKsi; KronrodNodesKsi=IntegrationNetKsi;
	for(int i=0; i<IntegrationNetKsi; i++)
	{
		NodeKsi[i]=(double)(i*2.0)/(double)IntegrationNetKsi - 1.0;
		KronWeightKsi[i]=1/(double)IntegrationNetKsi;
		GausWeightKsi[i]=1/(double)IntegrationNetKsi;
	}
}

}
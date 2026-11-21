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
#define USE_GAUSS_KRONROD 2
#include <stdio.h>
#include "kronrod.h"

void MagneticParameters::GenerateIntegrationNet()
{
if (IntegrationNet==1)
{
	GaussNodes=1; KronrodNodes=1;
	Node[ 0]= 1.0000e+00; KronWeight[ 0]=1.0000; GausWeight[ 0]=1.0000;
	return;
}
if (USE_GAUSS_KRONROD)
{
	if (IntegrationNet>MAX_GKNODES) 
		IntegrationNet=MAX_GKNODES;
	if (!(IntegrationNet%2)) IntegrationNet++;		// if integration net is even, make it odd, so the mirrored net is symmetric about a central node
	if (IntegrationNet<3) IntegrationNet=3;			// floor to the smallest net for which kronrod(IntegrationNet/2,...) receives a meaningful order (>=1)
	GaussNodes=IntegrationNet; KronrodNodes=IntegrationNet;
	kronrod ( IntegrationNet/2, 1.0e-9, Node, KronWeight, GausWeight );

	for(int i=IntegrationNet/2+1; i<IntegrationNet; i++)
	{
		Node[i]=-Node[i-IntegrationNet/2-1];
		KronWeight[i]=KronWeight[i-IntegrationNet/2-1];
		GausWeight[i]=GausWeight[i-IntegrationNet/2-1];
	}
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
	GaussNodes=IntegrationNet; KronrodNodes=IntegrationNet;
	for(int i=0; i<IntegrationNet; i++)
	{
		Node[i]=(double)(i*2.0)/(double)IntegrationNet - 1.0;
		KronWeight[i]=1/(double)IntegrationNet;
		GausWeight[i]=1/(double)IntegrationNet;
	}
}

}





void MagneticParameters::GenerateIntegrationNetKsi()
{
if (IntegrationNetKsi==1)
{
	GaussNodesKsi=1; KronrodNodesKsi=1;
	NodeKsi[ 0]= 1.0000e+00; KronWeightKsi[ 0]=1.0000; GausWeightKsi[ 0]=1.0000;
	return;
}

if (USE_GAUSS_KRONROD)
{
	if (IntegrationNetKsi>MAX_GKNODES) 
		IntegrationNetKsi=MAX_GKNODES;
	if (!(IntegrationNetKsi%2)) IntegrationNetKsi++;		// if integration net is even, make it odd, so the mirrored net is symmetric about a central node
	if (IntegrationNetKsi<3) IntegrationNetKsi=3;			// floor to the smallest net for which kronrod(IntegrationNetKsi/2,...) receives a meaningful order (>=1)
	GaussNodesKsi=IntegrationNetKsi; KronrodNodesKsi=IntegrationNetKsi;
	kronrod ( IntegrationNetKsi/2, 1.0e-9, NodeKsi, KronWeightKsi, GausWeightKsi );

	for(int i=IntegrationNetKsi/2+1; i<IntegrationNetKsi; i++)
	{
		NodeKsi[i]=-NodeKsi[i-IntegrationNetKsi/2-1];
		KronWeightKsi[i]=KronWeightKsi[i-IntegrationNetKsi/2-1];
		GausWeightKsi[i]=GausWeightKsi[i-IntegrationNetKsi/2-1];
	}
	for(int i=0; i<IntegrationNetKsi; i++)				// NORMALIZATION
	{
		KronWeightKsi[i]/=2.0;
		GausWeightKsi[i]/=2.0;
	}

}
else  // if (USE_GAUSS_KRONROD==0)
{
	GaussNodesKsi=IntegrationNetKsi; KronrodNodesKsi=IntegrationNetKsi;
	for(int i=0; i<IntegrationNetKsi; i++)
	{
		NodeKsi[i]=(double)(i*2.0)/(double)IntegrationNetKsi - 1.0;
		KronWeightKsi[i]=1/(double)IntegrationNetKsi;
		GausWeightKsi[i]=1/(double)IntegrationNetKsi;
	}
}

}
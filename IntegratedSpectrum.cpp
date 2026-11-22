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

// 
// Part of EPR spectrum calculation in the rigid limit.
// Individual spectra for all biradical orientations are integrated with orientation distribution function,
// to give the "powder" spectrum.

// File purpose: orientational (powder or ODF-weighted) average of the spectrum of one particle over all molecular
// orientations, for the rigid limit and for the slow-motion (Liouville) case; the existing comment above says
// rigid limit only, but CalculateIntegratedSpectrumLiouville() below also handles Liouville particles.
// Used by: Discrepancy.cpp (Discrepancy() calls CalculateIntegratedSpectrumLiouville() when IntegrationNet>1);
// IntegratedSpectrum.cpp itself uses VoitlanderIntegration.cpp and VoitlanderTriangleIntegrals.cpp for the Voitlander mesh.
#include "BiradicalODF.h"
#include "VoitlanderIntegration.h"
#include <math.h>
#include <iostream>

extern int MpiRank, MpiSize;

// forward declarations of the three rigid-limit callbacks defined in VoitlanderTriangleIntegrals.cpp,
// they are passed to IntegrateVoitlanderMesh() below
void VoitlanderVertexEvalRigid(MagneticParameters *MagneticParameter, VoitlanderVertex *Vertex);
void VoitlanderFreePayloadRigid(void *Payload);
void VoitlanderTriangleIntegralRigid(VoitlanderVertex *V1, VoitlanderVertex *V2, VoitlanderVertex *V3,
                                      double SolidAngle, int NPoints, double *FieldNet, double *Contribution);

// rigid-limit orientation average via the adaptive hemisphere triangulation
//    (VoitlanderIntegration.h) instead of the fixed Gauss/Gauss-Kronrod grid;
//    SLE/Liouville particles are not handled here (the per-triangle numerical
//    cubature callback a non-analytical SLE lineshape needs is explicitly
//    deferred to the future MOMD task, see
//    knowledge/epr_agents/orthos_sle_fwhm_voitlander_2026-11-02.md), so
//    CalculateIntegratedSpectrumLiouville only calls this for !Liouville particles
//    and falls through to its own existing grid otherwise
// CalculateIntegratedSpectrumVoitlander: orientation average for a rigid-limit particle by the adaptive hemisphere
// mesh; fills SpectrumSum[k][i] (k: spectrum number, i: field point) for the ODF components of ODF.
// Used by: CalculateIntegratedSpectrumLiouville() in this file (when IntegrationMethod is INTEGRATE_VOITLANDER and the
// particle is not in the Liouville regime).
void CalculateIntegratedSpectrumVoitlander(MagneticParameters *MagneticParameter,
                                            OrderParameters *ODF,
                                            double **SpectrumSum)
{
    // clear the spectrum accumulator
    for(int i=0; i<MagneticParameter->NPoints; i++)
        for(int k=0; k<ODF->NOrientations; k++)
            SpectrumSum[k][i]=0;

    // integrate over the hemisphere mesh, adding the ODF-weighted contribution of every triangle to SpectrumSum
    IntegrateVoitlanderMesh(MagneticParameter, ODF, VoitlanderVertexEvalRigid, VoitlanderFreePayloadRigid, VoitlanderTriangleIntegralRigid, SpectrumSum);

// the mesh engine accumulates a raw solid-angle integral in steradians (every
//    leaf triangle's true measure, summing to exactly 2*PI over the whole
//    hemisphere by construction); the existing Gauss/Gauss-Kronrod grid instead
//    accumulates a NORMALISED powder average (its quadrature weights are
//    pre-scaled to sum to 1 over the whole sphere, see GausWeight/GausWeightKsi
//    in GaussKronrodIntegration.cpp), so dividing by 2*PI here reduces this
//    mesh's raw integral to the same normalised-average convention; restricting
//    the mesh to one hemisphere needs no further factor of two, because the EPR
//    resonance condition is itself inversion-symmetric (Bres(n)=Bres(-n)), so a
//    hemisphere average already equals the full-sphere average
    for(int i=0; i<MagneticParameter->NPoints; i++)
        for(int k=0; k<ODF->NOrientations; k++)
            SpectrumSum[k][i]*=MagneticParameter->AmplitudeFactor/(36.0*2.0*PI);
}

// CalculateIntegratedSpectrumLiouville: orientation average of the spectrum of one particle over the molecular
// orientations (polar angle teta, azimuth ksi) by Gauss or Gauss-Kronrod quadrature, with each spectrum k weighted
// by the orientation distribution function of that sample orientation (ODF->CurrentODF[k]). The orientation spectrum
// is the rigid-limit lineshape sum or, for Liouville particles, the slow-motion spectrum of CalculateLiouvillianB().
// If the Voitlander method is selected for a rigid-limit particle the work is passed on to
// CalculateIntegratedSpectrumVoitlander() instead.
// MagneticParameter: parameters of one particle; ODF: orientation distribution of that particle;
// SpectrumSum: output [ODF->NOrientations][NPoints]. Under MPI each rank adds only its own share of the teta nodes,
// the sum over ranks is done afterwards in Discrepancy(). A log file integr<particle>_<rank>.out is written.
// Used by: Discrepancy() in Discrepancy.cpp, for every particle with IntegrationNet>1 (nothing else calls it).
void CalculateIntegratedSpectrumLiouville(MagneticParameters *MagneticParameter,
										  OrderParameters *ODF,
										  double **SpectrumSum)
{
// the Voitlander mesh is only derived for the rigid limit so far; Liouville
//    particles fall through unchanged to the grid below, selected exactly as
//    if IntegrateGaussKronrod had been chosen (see the comment above
//    CalculateIntegratedSpectrumVoitlander)
    if((MagneticParameter->IntegrationMethod==INTEGRATE_VOITLANDER)&&(!(MagneticParameter->Liouville)))
    {
        CalculateIntegratedSpectrumVoitlander(MagneticParameter, ODF, SpectrumSum);
        return;
    }

// Initializing Magnetic parameters
    double Field_to_Molecule_Euler[3];
    Matrix<double, 3, 3> RotMat, RotMat_T;
	double HGaussEff, HLorentzEff;
// Initializing magnetic field net related parameters
	double *Spectrum;
	Spectrum=(double*)malloc(MagneticParameter->NPoints*sizeof(double));

// make array for Gaussian integration
	double **SpectrumSumGauss;
	SpectrumSumGauss=new double* [ODF->NOrientations];
	for (int i=0; i<ODF->NOrientations; i++)
		SpectrumSumGauss[i]=new double[MagneticParameter->NPoints];
	

	// clear both spectrum accumulators
	for(int i=0; i<MagneticParameter->NPoints; i++)
		for(int k=0; k<ODF->NOrientations; k++)
		{
			SpectrumSum[k][i]=0;
			SpectrumSumGauss[k][i]=0;
		}

	// open the log file of this rank, with the integration net sizes in its header
	std::fstream intfile;
	char buff[MAX_STRLEN];

	sprintf(buff, "integr%d_%d.out", MagneticParameter->Count, MpiRank);

	intfile.open(buff, std::fstream::out);
	intfile << "Integration net: " << MagneticParameter->IntegrationNet << " over teta, " << MagneticParameter->IntegrationNetKsi << " over ksi." << endl;
	intfile << "Gauss nodes: " << MagneticParameter->GaussNodes << ", Kronrod Nodes: " << MagneticParameter->KronrodNodes << endl << endl;


// Starting integration over orientations:	
//	ODF->WriteToFile("temp_ODF");
#ifdef __MPI
	for(int iteta=StartMPI(MpiRank, MagneticParameter->IntegrationNet, MpiSize); iteta<EndMPI(MpiRank, MagneticParameter->IntegrationNet, MpiSize); iteta++)
#else
	for(int iteta=0; iteta<MagneticParameter->IntegrationNet; iteta++)
#endif	
	{
	  // loop over the integration nodes in ksi (main.cpp sets IntegrationNetKsi=1 for Liouville particles)
	  for(int iksi=0; iksi<MagneticParameter->IntegrationNetKsi; iksi++)
	  {
	  // Define integration angles and Euler angles Field->Molecule:
//		  double ksi=(double)iksi*2.0*PI/(double)MagneticParameter->IntegrationNetKsi;
//		  double costeta=-1.0+(double)iteta*2.0/(double)MagneticParameter->IntegrationNet;
		if(!(MagneticParameter->IntegrationMethod==INTEGRATE_GAUSSKRONROD))			// if we are integrating only on Gauss nodes, skip nodes, which are Kronrod-only.
			if((MagneticParameter->GausWeightKsi[iksi]==0)||(MagneticParameter->GausWeight[iteta]==0)) continue;



		  // ksi from the node mapped from [-1,1] to [0,2*PI], cos(teta) taken directly from the node
		  double ksi=PI*(MagneticParameter->NodeKsi[iksi]+1.0);
		  double costeta=MagneticParameter->Node[iteta];
		  double teta=acos(costeta);

		  Field_to_Molecule_Euler[0]=0;
		  Field_to_Molecule_Euler[1]=teta;
          Field_to_Molecule_Euler[2]=ksi;

		if(!(MagneticParameter->Liouville))		// rigid limit
		{
		  // rotation matrix from the Euler angles field->molecule, then transitions and lineshape sum for this orientation
		  RotMat=RotationMatrix(Field_to_Molecule_Euler);
		  if(!(MagneticParameter->isMonoradical))
		  {
		  // Solving Spin-Hamiltonian and calculating ESR transitions for given orientation
			  Array<double, 54, 2> Transitions=CalculateTransitions(MagneticParameter, RotMat, &HGaussEff, &HLorentzEff);
		  // Calculating spectrum lineshape for given orientation
			  CalculateIndividualSpectrum (Transitions, HGaussEff, HLorentzEff, MagneticParameter->NPoints, MagneticParameter->FieldNet, Spectrum, MagneticParameter->TransitionIntensityCutoff);
		  }
		  else
		  {
		  // Solving Spin-Hamiltonian and calculating ESR transitions for given orientation
			  Array<double, 3, 2> Transitions=CalculateTransitionsM(MagneticParameter, RotMat, &HGaussEff, &HLorentzEff);
		  // Calculating spectrum lineshape for given orientation
			  CalculateIndividualSpectrumM (Transitions, HGaussEff, HLorentzEff, MagneticParameter->NPoints, MagneticParameter->FieldNet, Spectrum, MagneticParameter->TransitionIntensityCutoff);
		  }


		}
		else		 // Liouvillian calculation
		{
			// slow-motion spectrum for this teta, by solving the Liouville-space problem over the whole field net
			Matrix<double, Dynamic, 1> SpectrumL;
			SpectrumL.resize(MagneticParameter->NPoints , 1);		
// the basis vectors are generated from L,M,K truncation values and selection rules.
			SpectrumL=CalculateLiouvillianB((MagneticParameter), teta, "truncation" );
			for(int j=0; j<MagneticParameter->NPoints; j++)
				Spectrum[j]=SpectrumL(j, 0);					// Translate the spectrum into double* format
		}


	  // Number density of particles oriented at angles teta, ksi in resonator are pre-calculated and stored in ODF->CurrentODF[k]:
		  ODF->ODF_reduced(Field_to_Molecule_Euler[1], ksi);

	  // Adding contribution of given orientation to spectrum:
		  for(int k=0; k<ODF->NOrientations; k++)
			  for(int i=0; i<MagneticParameter->NPoints; i++)
			  {
			    // Gauss-Kronrod: add with Kronrod weights (result) and with Gauss weights (only written to the log file); otherwise add with Gauss weights
			    if(MagneticParameter->IntegrationMethod==INTEGRATE_GAUSSKRONROD)
				{
				  SpectrumSum[k][i]     +=(ODF->CurrentODF[k])*Spectrum[i]*MagneticParameter->KronWeight[iteta]*MagneticParameter->KronWeightKsi[iksi];
				  SpectrumSumGauss[k][i]+=(ODF->CurrentODF[k])*Spectrum[i]*MagneticParameter->GausWeight[iteta]*MagneticParameter->GausWeightKsi[iksi];
				}
				else
				{
				  SpectrumSum[k][i]+=(ODF->CurrentODF[k])*Spectrum[i]*MagneticParameter->GausWeight[iteta]*MagneticParameter->GausWeightKsi[iksi];
				}
			  }
	  }
	}


	  // Normalization:
   // rigid-limit spectra are scaled by AmplitudeFactor/36 here, Liouville spectra are not rescaled in this function
   if(!(MagneticParameter->Liouville))
	  for(int i=0; i<MagneticParameter->NPoints; i++)
		   for (int k=0; k<ODF->NOrientations; k++)
		   {
			     SpectrumSum[k][i]      *= MagneticParameter->AmplitudeFactor/36.0;
			     SpectrumSumGauss[k][i] *= MagneticParameter->AmplitudeFactor/36.0;
		   }

// Output:


	// write the spectra of all sample orientations (both weighted sums) to the log file
	for (int k=0; k<ODF->NOrientations; k++)
   		for(int i=0; i<MagneticParameter->NPoints; i++)
		{
			intfile << k << ", " << i << ", " << i+k*MagneticParameter->NPoints << ", " << SpectrumSum[k][i] << ", " << SpectrumSumGauss[k][i] << endl;
		}

		intfile.close();

// the Gauss-quadrature accumulator and the scratch spectrum are local to this call
		for (int i=0; i<ODF->NOrientations; i++)
			delete[] SpectrumSumGauss[i];
		delete[] SpectrumSumGauss;
		free(Spectrum);

}

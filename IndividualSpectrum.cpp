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

// Part of EPR spectrum calculation in the rigid limit.
// For given orientation of biradical molecule, EPR spectrum with convolution of Lorentzizan and Gaussian lineshape is calculated, 
//    based on the transition fields and intensities calculated in function BiradicalTransitions().
//

#include "BiradicalODF.h"
#include <math.h>
#include <iostream>

#include "Faddeeva.hh"

// Linewidth convention used throughout this program: every HG and HL quantity, from the
//    GaussTensor and LorentzTensor entries of Magnetic.par through the fitted parameters and the
//    orientation-projected effective widths HGaussEff and HLorentzEff, is a FULL WIDTH AT HALF
//    MAXIMUM in Gauss. For the Gaussian component HG=2*sqrt(2*ln2)*sigma, for the Lorentzian
//    component HL=2*gamma.
// Each of the three functions below returns the first derivative with respect to the field of an
//    absorption lineshape of UNIT AREA. The absorption therefore integrates to 1 and the
//    derivative returned here double-integrates to 1, independently of HG and HL.
// The Voigt function uses the standard Faddeeva parameterisation of the Voigt profile,
//    z=sqrt(ln2)*(2*Hrel+i*HL)/HG and V(Hrel)=2*sqrt(ln2/PI)*real(w(z))/HG, so that it reduces
//    exactly to the Gaussian as HL->0 and exactly to the Lorentzian as HG->0, with no step at the
//    width ratio where the dispatch below switches between them.

double GaussDerivativeLineshape(double Hrel, double HG, double HL)
{
  return -16*pow(log(2.0), 1.5)*Hrel*exp(-4*log(2.0)*Hrel*Hrel/(HG*HG))/(sqrt(PI)*HG*HG*HG);
}

double LorentzDerivativeLineshape(double Hrel, double HG, double HL)
{
  return -16*Hrel*HL*pow(4*Hrel*Hrel+HL*HL, -2)/(PI);
}

double VoigtDerivativeLineshape(double Hrel, double HG, double HL)
{
	std::complex<double> z(2*sqrt(log(2.0))*Hrel/HG, sqrt(log(2.0))*HL/HG);
	double result=-real(z*Faddeeva::w(z))*8*log(2.0)/(HG*HG*sqrt(PI));
	return result;
}

// TODO: Convolution!

void CalculateIndividualSpectrum (Array<double, 54, 2> &Transitions,
								  double HGaussEff,
								  double HLorentzEff,
								  int NPoints,
								  double *FieldNet,
								  double *SpectrumAmplitude,
								  double IntensityCutoff)
{
	for(int i=0; i<NPoints; i++)
	{
		SpectrumAmplitude[i]=0;
//		HGaussEff=fabs(HGaussEff);
//		HLorentzEff=fabs(HLorentzEff);
		double CurrentField=FieldNet[i];
		for (int j=0; j<54; j++)
		{
			if(Transitions(j, 1)<=IntensityCutoff) continue;
			double amp, Hrel=CurrentField-Transitions(j, 0);
			if(HGaussEff<=HLorentzEff*0.001) amp=LorentzDerivativeLineshape(Hrel, HGaussEff, HLorentzEff);
			else if (HLorentzEff<=HGaussEff*0.001) amp=GaussDerivativeLineshape(Hrel, HGaussEff, HLorentzEff);
			else amp=VoigtDerivativeLineshape(Hrel, HGaussEff, HLorentzEff);
            SpectrumAmplitude[i]+=Transitions(j, 1)*amp;
		}
	}


}

// Same for monoradical
void CalculateIndividualSpectrumM(Array<double, 3, 2> &Transitions,
								  double HGaussEff,
								  double HLorentzEff,
								  int NPoints,
								  double *FieldNet,
								  double *SpectrumAmplitude,
								  double IntensityCutoff)
{
	for(int i=0; i<NPoints; i++)
	{
		SpectrumAmplitude[i]=0;
//		HGaussEff=fabs(HGaussEff);
//		HLorentzEff=fabs(HLorentzEff);
		double CurrentField=FieldNet[i];
		for (int j=0; j<3; j++)
		{
			if(Transitions(j, 1)<=IntensityCutoff) continue;
			double amp, Hrel=CurrentField-Transitions(j, 0);
			if(HGaussEff<=HLorentzEff*0.001) amp=LorentzDerivativeLineshape(Hrel, HGaussEff, HLorentzEff);
			else if (HLorentzEff<=HGaussEff*0.001) amp=GaussDerivativeLineshape(Hrel, HGaussEff, HLorentzEff);
			else amp=VoigtDerivativeLineshape(Hrel, HGaussEff, HLorentzEff);
            SpectrumAmplitude[i]+=Transitions(j, 1)*amp;
		}
	}


}

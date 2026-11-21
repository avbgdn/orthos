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

// Rigid-limit vertex evaluation and the six per-triangle integral forms plugged
//    into the generic engine of VoitlanderIntegration.cpp as one combined
//    VoitlanderTriangleIntegral callback (VoitlanderTriangleIntegralRigid). Each
//    of the six forms below is derived in full in
//    knowledge/epr_agents/orthos_sle_fwhm_voitlander_2026-11-02.md; only the
//    conclusions are restated here as code comments.
// Convention: within one triangle, the j-th transition's resonance position and
//    intensity are taken to vary affinely in the triangle's barycentric
//    coordinates, matched across the three vertices by transition INDEX (the
//    same ordinal slot of the Array<double,54,2>/Array<double,3,2> that
//    CalculateTransitions/CalculateTransitionsM already return); this is only
//    accurate once the refine test has made the triangle small enough that
//    transition ordering has not changed between its three vertices.

#include "VoitlanderIntegration.h"
#include <cmath>
#include <algorithm>
#include <iostream>

#define MAX_TRANSITIONS 54

// one vertex's cached rigid-limit physics: up to 54 (position, intensity)
//    transition pairs plus the orientation-projected effective Gaussian and
//    Lorentzian widths, exactly as CalculateTransitions/CalculateTransitionsM
//    already compute them for the existing fixed-grid integration
struct VoitlanderRigidPayload
{
    double Position[MAX_TRANSITIONS];
    double Intensity[MAX_TRANSITIONS];
    double HGaussEff, HLorentzEff;
    double IntensityCutoff;
};

// evaluates one orientation via the existing CalculateTransitions (biradical) or
//    CalculateTransitionsM (monoradical) routine, exactly as the fixed Gauss and
//    Gauss-Kronrod grids already do; the generic refine-test proxy (Position,
//    Width) is reduced from this to the intensity-weighted mean transition
//    position and the narrowest of the two orientation-projected widths
void VoitlanderVertexEvalRigid(MagneticParameters *MagneticParameter, VoitlanderVertex *Vertex)
{
    double eulerangles[3]={0, Vertex->Teta, Vertex->Ksi};
    Matrix<double, 3, 3> rotmat=RotationMatrix(eulerangles);

    VoitlanderRigidPayload *payload=new VoitlanderRigidPayload;
    double hgausseff, hlorentzeff;
    if(!MagneticParameter->isMonoradical)
    {
        Array<double, 54, 2> transitions=CalculateTransitions(MagneticParameter, rotmat, &hgausseff, &hlorentzeff);
        for(int j=0; j<54; j++)
        {
            payload->Position[j]=transitions(j, 0);
            payload->Intensity[j]=transitions(j, 1);
        }
    }
    else
    {
        Array<double, 3, 2> transitions=CalculateTransitionsM(MagneticParameter, rotmat, &hgausseff, &hlorentzeff);
        for(int j=0; j<3; j++)
        {
            payload->Position[j]=transitions(j, 0);
            payload->Intensity[j]=transitions(j, 1);
        }
        for(int j=3; j<MAX_TRANSITIONS; j++)
        {
            payload->Position[j]=0;
            payload->Intensity[j]=0;
        }
    }
    payload->HGaussEff=hgausseff;
    payload->HLorentzEff=hlorentzeff;
    payload->IntensityCutoff=MagneticParameter->TransitionIntensityCutoff;
    Vertex->Payload=payload;

    // position of the single STRONGEST transition: a single representative
    // resonance-position proxy for the refine test, which cannot track all 54
    // transitions individually without knowing the engine's generic contract.
    // An intensity-weighted MEAN over all active transitions was tried first,
    // but a weak transition flickering in and out of the IntensityCutoff
    // between two almost-identical orientations moves that mean by a large,
    // spurious amount and defeats the refine test's convergence; the single
    // strongest transition's identity and position are far more stable under
    // small orientation changes, at the cost of not tracking weaker,
    // well-separated lines explicitly, a known limitation (confession section)
    double strongestintensity=-1, strongestposition=0;
    for(int j=0; j<MAX_TRANSITIONS; j++)
        if(payload->Intensity[j]>strongestintensity)
        {
            strongestintensity=payload->Intensity[j];
            strongestposition=payload->Position[j];
        }
    Vertex->Position=strongestposition;

    // narrowest of the two orientation-projected widths, with two numerical
    // floors: (1) a tensor input of EXACTLY zero (e.g. a dataset with no
    // Gaussian component at all, GaussTensor=0) means that channel is simply
    // absent, not narrower than the other one, so it must not be allowed to
    // pull the minimum down to zero; (2) CalculateTransitions' effective-width
    // projection can also legitimately pass through (near) zero at a genuine
    // orientation-space turning point of Bres even when both channels are
    // nominally present, which the fixed Gauss/Gauss-Kronrod grid never
    // notices (it never divides by the width) but this refine test does. Both
    // cases are floored the same way: treat anything below 1e-3 Gauss (this
    // dataset's empirically observed numerical noise floor) as absent, and
    // take the narrowest of whatever remains
    double hgausswidth=fabs(payload->HGaussEff), hlorentzwidth=fabs(payload->HLorentzEff);
    if(hgausswidth<1e-3) hgausswidth=1e300;
    if(hlorentzwidth<1e-3) hlorentzwidth=1e300;
    Vertex->Width=std::max(1e-3, std::min(hgausswidth, hlorentzwidth));
}

void VoitlanderFreePayloadRigid(void *Payload)
{
    delete (VoitlanderRigidPayload*)Payload;
}

// unit-area Gaussian absorption profile (not its derivative); FWHM convention
//    HG=2*sqrt(2*ln2)*sigma, matching IndividualSpectrum.cpp
static double GaussAbsorption(double Hrel, double HG, double HL)
{ return 2.0*sqrt(log(2.0)/PI)/HG*exp(-4.0*log(2.0)*Hrel*Hrel/(HG*HG)); }

// antiderivative of GaussAbsorption: a rescaled error function, since
//    d/dx[(1/2)(1+erf(ax))]=(a/sqrt(pi))*exp(-a^2 x^2) matches GaussAbsorption
//    for a=2*sqrt(ln2)/HG
static double GaussCDF(double Hrel, double HG, double HL)
{ return 0.5*(1.0+erf(2.0*sqrt(log(2.0))*Hrel/HG)); }

// unit-area Lorentzian absorption profile; FWHM convention HL=2*gamma
static double LorentzAbsorption(double Hrel, double HG, double HL)
{ return (1.0/PI)*(HL/2.0)/(Hrel*Hrel+HL*HL/4.0); }

// antiderivative of LorentzAbsorption: the pole decomposition
//    1/(x^2+a^2)=(1/2ia)[1/(x-ia)-1/(x+ia)] integrates to an arctan
static double LorentzCDF(double Hrel, double HG, double HL)
{ return (1.0/PI)*atan(2.0*Hrel/HL)+0.5; }

typedef double (*VoitlanderLineshapeFn)(double Hrel, double HG, double HL);

// isotropic (zeroth-order) triangle contribution for one transition: the
//    triangle's resonance position and widths are treated as constant at their
//    three-vertex average, so the flat-triangle area integral of the derivative
//    lineshape is just the lineshape evaluated once, times the triangle's true
//    solid angle; this is the literal unified Task A/commit 39013f8 derivative
//    function, reused unchanged (Gauss, Lorentz and Voigt alike)
static double IsotropicContribution(double B, double Bbar, double HGbar, double HLbar, double Intensity, double SolidAngle, VoitlanderLineshapeFn Deriv)
{ return Intensity*SolidAngle*Deriv(B-Bbar, HGbar, HLbar); }

// anisotropic (first-order) triangle contribution for one transition, exact for
//    an affine resonance position Bres(u,v)=B1+(B2-B1)u+(B3-B1)v over the flat
//    reference triangle (0,0),(1,0),(0,1): writing f'(B-Bres)=
//    -(1/|g|^2) g.grad[f(B-Bres)] (g=(B2-B1,B3-B1) is the constant gradient of
//    Bres) and applying the divergence theorem turns the area integral of the
//    derivative lineshape f' into a sum of three boundary line integrals of the
//    absorption lineshape f itself, each evaluated via its own antiderivative
//    (CDF): erf-based for Gauss (the "Green's theorem" reduction), arctan-based
//    for Lorentz (the "pole decomposition" of 1/(x^2+a^2)). Degenerate edges
//    (equal endpoint position) and a near-zero gradient both reduce smoothly to
//    the isotropic case above, which this function returns directly once the
//    gradient is numerically negligible compared to the width.
static double AnisotropicTriangleIntegral(double B, double B1, double B2, double B3, double HG, double HL, double WidthScale,
                                           VoitlanderLineshapeFn Deriv, VoitlanderLineshapeFn Absorp, VoitlanderLineshapeFn CDF)
{
    double gx=B2-B1, gy=B3-B1;
    double gsq=gx*gx+gy*gy;
    double Bbar=(B1+B2+B3)/3.0;
    if(gsq<1e-12*WidthScale*WidthScale) return 0.5*Deriv(B-Bbar, HG, HL);

    double Bstart[3]={B1, B2, B3}, Bend[3]={B2, B3, B1};
// g.n_edge for the three edges of the (0,0),(1,0),(0,1) reference triangle,
// traversed counter-clockwise, with outward normal*length (dv,-du) per edge
    double gdotn[3]={B1-B3, B2+B3-2.0*B1, B1-B2};

    double sum=0;
    for(int e=0; e<3; e++)
    {
        double deltaB=Bend[e]-Bstart[e];
        double edgeintegral=(fabs(deltaB)<1e-12*WidthScale) ? Absorp(B-Bstart[e], HG, HL) : (CDF(B-Bstart[e], HG, HL)-CDF(B-Bend[e], HG, HL))/deltaB;
        sum+=gdotn[e]*edgeintegral;
    }
    return -sum/gsq;
}

// anisotropic Voigt: a fixed 3-point Hammer-node cubature, degree-two accurate,
//    weight 1/3 per node, matching Spinach's vwcubature.m cub_order=3 rule
//    exactly (kernel/line_shapes/vwcubature.m, case 3); the three barycentric
//    nodes are (2/3,1/6,1/6) and its two cyclic permutations, placed strictly
//    inside the triangle rather than on its edges, and the resonance position
//    and both linewidth parameters are barycentrically interpolated to each
//    node using that node's own weights, mirroring Spinach's nu/g/s construction
static double AnisotropicVoigtCubature(double B, double B1, double B2, double B3, double HG1, double HG2, double HG3, double HL1, double HL2, double HL3, double Intensity, double SolidAngle)
{
    double lam[3][3]={{2.0/3.0, 1.0/6.0, 1.0/6.0}, {1.0/6.0, 2.0/3.0, 1.0/6.0}, {1.0/6.0, 1.0/6.0, 2.0/3.0}};
    double sum=0;
    for(int node=0; node<3; node++)
    {

// barycentric interpolation of the resonance position and both linewidths
//    to this Hammer node, using the node's own weights
        double Bnode=lam[node][0]*B1+lam[node][1]*B2+lam[node][2]*B3;
        double HGnode=lam[node][0]*HG1+lam[node][1]*HG2+lam[node][2]*HG3;
        double HLnode=lam[node][0]*HL1+lam[node][1]*HL2+lam[node][2]*HL3;
        sum+=VoigtDerivativeLineshape(B-Bnode, HGnode, HLnode);
    }
    return Intensity*SolidAngle*sum/3.0;
}

// combined rigid-limit per-triangle integral, plugged into the generic engine
//    as the one VoitlanderTriangleIntegral callback for IntegrateVoitlander;
//    dispatches per transition between the six forms above using the same
//    HG<=0.001*HL / HL<=0.001*HG convention IndividualSpectrum.cpp already uses
//    to pick Gauss/Lorentz/Voigt, plus a resonance-position-gradient check to
//    pick isotropic vs anisotropic within the Gauss/Lorentz regimes
void VoitlanderTriangleIntegralRigid(VoitlanderVertex *V1, VoitlanderVertex *V2, VoitlanderVertex *V3,
                                      double SolidAngle, int NPoints, double *FieldNet, double *Contribution)
{
    VoitlanderRigidPayload *p1=(VoitlanderRigidPayload*)V1->Payload;
    VoitlanderRigidPayload *p2=(VoitlanderRigidPayload*)V2->Payload;
    VoitlanderRigidPayload *p3=(VoitlanderRigidPayload*)V3->Payload;

    for(int i=0; i<NPoints; i++) Contribution[i]=0;

    double HGbar=(p1->HGaussEff+p2->HGaussEff+p3->HGaussEff)/3.0;
    double HLbar=(p1->HLorentzEff+p2->HLorentzEff+p3->HLorentzEff)/3.0;

    for(int j=0; j<MAX_TRANSITIONS; j++)
    {
        if((p1->Intensity[j]<=p1->IntensityCutoff)||(p2->Intensity[j]<=p2->IntensityCutoff)||(p3->Intensity[j]<=p3->IntensityCutoff))
            continue;

        double B1=p1->Position[j], B2=p2->Position[j], B3=p3->Position[j];
        double intensitybar=(p1->Intensity[j]+p2->Intensity[j]+p3->Intensity[j])/3.0;
        double span=std::max(B1, std::max(B2, B3))-std::min(B1, std::min(B2, B3));

        for(int i=0; i<NPoints; i++)
        {
            double B=FieldNet[i], value;
            if(HGbar<=HLbar*0.001)
            {
// pure Lorentzian regime
                if(span<=1e-3*HLbar)
                    value=IsotropicContribution(B, (B1+B2+B3)/3.0, HGbar, HLbar, intensitybar, SolidAngle, LorentzDerivativeLineshape);
                else
                    value=intensitybar*SolidAngle*2.0*AnisotropicTriangleIntegral(B, B1, B2, B3, HGbar, HLbar, HLbar, LorentzDerivativeLineshape, LorentzAbsorption, LorentzCDF);
            }
            else if(HLbar<=HGbar*0.001)
            {
// pure Gaussian regime
                if(span<=1e-3*HGbar)
                    value=IsotropicContribution(B, (B1+B2+B3)/3.0, HGbar, HLbar, intensitybar, SolidAngle, GaussDerivativeLineshape);
                else
                    value=intensitybar*SolidAngle*2.0*AnisotropicTriangleIntegral(B, B1, B2, B3, HGbar, HLbar, HGbar, GaussDerivativeLineshape, GaussAbsorption, GaussCDF);
            }
            else
            {
// Voigt regime: comparable Gaussian and Lorentzian widths
                double narrowest=std::min(HGbar, HLbar);
                if(span<=1e-3*narrowest)
                    value=IsotropicContribution(B, (B1+B2+B3)/3.0, HGbar, HLbar, intensitybar, SolidAngle, VoigtDerivativeLineshape);
                else
                    value=AnisotropicVoigtCubature(B, B1, B2, B3, p1->HGaussEff, p2->HGaussEff, p3->HGaussEff, p1->HLorentzEff, p2->HLorentzEff, p3->HLorentzEff, intensitybar, SolidAngle);
            }
            Contribution[i]+=value;
        }
    }
}

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

// Generic adaptive-triangulation orientation-integration engine ("Voitlander"),
//    named after Spinach's kernels/voitlander.m, whose per-triangle analytical
//    forms inspired the ones in VoitlanderTriangleIntegrals.cpp. This header
//    declares only the engine's generic contract (mesh vertex record, the two
//    injected callbacks, and the driving entry point); the engine itself
//    (hemisphere mesh, recursive refinement, vertex caching, the refine test,
//    static per-seed-triangle MPI split) lives in VoitlanderIntegration.cpp and
//    is completely agnostic of which physical model filled the vertices.
// A vertex is evaluated once by VoitlanderVertexEval and is shared, together
//    with its cached Payload, by every triangle (at every subdivision depth)
//    that touches it. Position and Width are a generic resonance-position and
//    linewidth proxy, read only by the engine's refine test; a method-specific
//    callback may reduce an arbitrarily rich physics result (a parametric list
//    of transitions for the rigid limit, or a full computed SLE lineshape for a
//    future MOMD extension) down to this same pair of scalars without the
//    engine ever inspecting Payload itself.

#ifndef VOITLANDERINTEGRATION_H
#define VOITLANDERINTEGRATION_H

#include "BiradicalODF.h"

struct VoitlanderVertex
{
    double Teta, Ksi;                  // orientation Euler angles (polar, azimuthal)
    double Position, Width;            // generic resonance-position / linewidth proxy, refine test only
    double CurrentODF[MAX_ORIENTATIONS];   // ODF weight at this orientation, one per ODF->NOrientations component
    void *Payload;                     // method-specific cached physics data, opaque to the engine
};

// evaluates one orientation and fills Vertex->Position, Vertex->Width and Vertex->Payload;
//    Vertex->Teta, Vertex->Ksi and Vertex->CurrentODF are already set by the engine itself
typedef void (*VoitlanderVertexEval)(MagneticParameters *MagneticParameter, VoitlanderVertex *Vertex);

// releases a Payload allocated by a VoitlanderVertexEval callback
typedef void (*VoitlanderFreePayload)(void *Payload);

// overwrites Contribution[0..NPoints-1] with one leaf triangle's field-domain
//    spectral contribution (ODF weighting and orientation-index bookkeeping are
//    handled by the engine, not by this callback); SolidAngle is the triangle's
//    true measure on the unit hemisphere, in steradians
typedef void (*VoitlanderTriangleIntegral)(VoitlanderVertex *V1, VoitlanderVertex *V2, VoitlanderVertex *V3,
                                            double SolidAngle, int NPoints, double *FieldNet, double *Contribution);

// drives the generic hemisphere mesh/refine/vertex-cache engine to the requested
//    accuracy using the injected VertexEval/TriangleIntegral pair, and adds the
//    ODF-weighted result into SpectrumSum[0..ODF->NOrientations-1][0..NPoints-1]
void IntegrateVoitlanderMesh(MagneticParameters *MagneticParameter, OrderParameters *ODF,
                              VoitlanderVertexEval VertexEval, VoitlanderFreePayload FreePayload,
                              VoitlanderTriangleIntegral TriangleIntegral, double **SpectrumSum);

#endif

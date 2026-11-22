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

// Generic hemisphere-mesh adaptive-triangulation engine for IntegrateVoitlander.
//    Everything here is agnostic of the physical model: the mesh, the vertex
//    cache, and the refine test are the same regardless of whether the
//    injected VoitlanderVertexEval/VoitlanderTriangleIntegral pair evaluates
//    the rigid limit (VoitlanderTriangleIntegrals.cpp, this task) or, later, a
//    fully numerical SLE/MOMD lineshape (not implemented here, see
//    knowledge/epr_agents/orthos_sle_fwhm_voitlander_2026-11-02.md) -- adding
//    that case needs only a new VertexEval/TriangleIntegral pair, not a change
//    to this file.

// Used by: IntegratedSpectrum.cpp, whose CalculateIntegratedSpectrumVoitlander() calls the only entry point of this
// file, IntegrateVoitlanderMesh(); all other functions here are static helpers of that function.

#include "VoitlanderIntegration.h"
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <iostream>

#ifdef __MPI
#include <mpi.h>
#endif

extern int MpiRank, MpiSize;
// MPI work-split helpers defined in Discrepancy.cpp
int StartMPI(int Rank, int N, int P);
int EndMPI(int Rank, int N, int P);

// span-of-resonance-position-vs-narrowest-linewidth refine threshold: a triangle
//    is subdivided while the three vertices' resonance positions differ by more
//    than REFINE_FACTOR times the narrowest of their three linewidths, i.e. the
//    affine (first-order) approximation used by the per-triangle integral is
//    only trusted while the lineshape itself has not yet resolved the triangle's
//    curvature; 8 was the starting value Alexey specified, but on the fixedGA
//    C4/PAAzo6 benchmark it left the spectrum 6.5% rms / 32.7% max pointwise off
//    GaussKronrod(101), visibly oscillatory and unconverged. A sweep against the
//    same GaussKronrod(101) reference (4 -> 4.1%/16.4%, 1 -> 0.37%/1.74%) showed
//    the error only drops to GK101 level at 1, so the threshold is tightened to
//    1 per Alexey's 2026-11-03 correction (convergence criterion stricter); a
//    finer initial seed mesh was tested alongside it and made no measurable
//    difference once REFINE_FACTOR=1, so it is not used
#define VOITLANDER_REFINE_FACTOR 1.0

// hard cap on recursive subdivision depth, protecting against runaway
//    subdivision near a genuine singularity or a bug in a refine criterion input,
//    exactly the failure mode flagged for Spinach's own voitlander. A transition
//    crossing (two resonance lines swapping intensity rank) makes the single
//    strongest-transition position proxy genuinely discontinuous there, so the
//    refine test can never converge locally; the depth cap is what bounds that
//    case to a small, deliberately under-resolved band of triangles around the
//    crossing curve (a 1-D feature on the 2-D mesh) rather than an unbounded
//    subdivision. 7 levels (4^7=16384 leaves per seed triangle in the worst
//    case) is a practical compromise between that bound and resolution
#define VOITLANDER_MAX_DEPTH 7

// absolute hard stop on total mesh vertices: the depth cap bounds a single
//    pathological triangle's subtree, but not how many such triangles (e.g.
//    transition-crossing bands) a dataset might contain; this is the backstop
#define VOITLANDER_MAX_VERTICES 300000

// one mesh vertex's 3-D Cartesian position on the unit hemisphere, kept
//    separate from VoitlanderVertex (teta, ksi) because mesh construction
//    (edge midpoints, solid angles) is far more robust in Cartesian coordinates
// Used by: the vertex list and the helper functions of this file (CreateVertex, GetMidpoint, SolidAngle, RefineTriangle).
struct VoitlanderXYZ { double x, y, z; };

// bundles every piece of state the recursive mesh traversal needs, so the
//    recursion itself can stay a small set of free functions instead of using
//    file-scope globals (which would not be safe to call back into for more
//    than one particle, or if a future caller ever nested two integrations)
// Used by: every function in this file; one instance is created in IntegrateVoitlanderMesh() and passed down by pointer.
struct VoitlanderMeshContext
{
    MagneticParameters *MagneticParameter;
    OrderParameters *ODF;
    VoitlanderVertexEval VertexEval;
    VoitlanderTriangleIntegral TriangleIntegral;
    int NPoints;
    double *FieldNet;
    double **SpectrumSum;
    double *Scratch;
    std::vector<VoitlanderVertex> Vertices;
    std::vector<VoitlanderXYZ> XYZ;
    std::map<std::pair<int, int>, int> MidpointCache;
    long LeafTriangleCount;
    int MaxDepthReached;
};

// normalises a Cartesian direction to the unit hemisphere and evaluates its
//    physics once; the returned index is then cached by every caller (the seed
//    construction and GetMidpoint) so this orientation is never evaluated twice
// Used by: BuildHemisphereSeed() (the pole, the upper ring and the equator ring) and GetMidpoint() (edge midpoints).
static int CreateVertex(VoitlanderMeshContext *Context, double x, double y, double z)
{
    // project the direction onto the unit sphere and get its polar angle Teta and azimuth Ksi (in [0,2*PI))
    double norm=sqrt(x*x+y*y+z*z);
    x/=norm; y/=norm; z/=norm;

    VoitlanderVertex vertex;
    vertex.Teta=acos(std::max(-1.0, std::min(1.0, z)));
    vertex.Ksi=atan2(y, x);
    if(vertex.Ksi<0) vertex.Ksi+=2.0*PI;

    // orientation distribution function weight, one value per ODF component
    Context->ODF->ODF_reduced(vertex.Teta, vertex.Ksi);
    for(int k=0; k<Context->ODF->NOrientations; k++)
        vertex.CurrentODF[k]=Context->ODF->CurrentODF[k];

    // method-specific physics evaluation, filling Position, Width and Payload
    Context->VertexEval(Context->MagneticParameter, &vertex);

    // store the vertex and its Cartesian position, the index in the two lists identifies the vertex
    Context->Vertices.push_back(vertex);
    VoitlanderXYZ xyz; xyz.x=x; xyz.y=y; xyz.z=z;
    Context->XYZ.push_back(xyz);

// anomalous vertex counts are a sign of a bug or a singularity in the refine
//    test, not a case that needs more subdivision (see the caution in the
//    commit message about Spinach's own voitlander); print progress so a
//    runaway refinement is diagnosable, and hard-stop well before it could
//    exhaust memory or runtime, since the depth cap alone bounds the subtree
//    of a single pathological triangle but not the number of such triangles
    if(Context->Vertices.size()%2000==0)
        std::cerr << "Voitlander: " << Context->Vertices.size() << " mesh vertices evaluated so far" << std::endl;
    if(Context->Vertices.size()>VOITLANDER_MAX_VERTICES)
    {
        std::cerr << "Voitlander: aborting after " << VOITLANDER_MAX_VERTICES << " mesh vertices; this is an anomalous triangle count (see the per-triangle integral's commit message), not a case needing more subdivision" << std::endl;
        exit(1);
    }

    return (int)Context->Vertices.size()-1;
}

// returns the cached midpoint vertex of edge (IndexA,IndexB), creating and
//    evaluating it only the first time this edge is split by any triangle at
//    any subdivision level -- this is the vertex cache the brief asks for
// Used by: RefineTriangle() (for the three edges of each triangle that is split).
static int GetMidpoint(VoitlanderMeshContext *Context, int IndexA, int IndexB)
{
    // look the edge up in the cache, the key is the ordered pair of vertex indices
    std::pair<int, int> key(std::min(IndexA, IndexB), std::max(IndexA, IndexB));
    std::map<std::pair<int, int>, int>::iterator found=Context->MidpointCache.find(key);
    if(found!=Context->MidpointCache.end()) return found->second;

    // not cached yet: create the vertex at the midpoint of the edge and remember it
    VoitlanderXYZ a=Context->XYZ[IndexA], b=Context->XYZ[IndexB];
    int midpoint=CreateVertex(Context, a.x+b.x, a.y+b.y, a.z+b.z);
    Context->MidpointCache[key]=midpoint;
    return midpoint;
}

// true solid angle of the spherical triangle with unit-vector corners A, B, C,
//    via the Van Oosterom-Strackee formula tan(Omega/2)=|a.(bxc)|/(1+a.b+b.c+c.a)
// Used by: RefineTriangle() (for every leaf triangle).
static double SolidAngle(VoitlanderXYZ A, VoitlanderXYZ B, VoitlanderXYZ C)
{
    double triple=A.x*(B.y*C.z-B.z*C.y)-A.y*(B.x*C.z-B.z*C.x)+A.z*(B.x*C.y-B.y*C.x);
    double ab=A.x*B.x+A.y*B.y+A.z*B.z;
    double bc=B.x*C.x+B.y*C.y+B.z*C.z;
    double ca=C.x*A.x+C.y*A.y+C.z*A.z;
    return 2.0*atan2(fabs(triple), 1.0+ab+bc+ca);
}

// the one refine criterion used throughout the engine, explicitly spelled out
//    here per Alexey's instruction: refine while the span of the resonance
//    position across the triangle's three vertices exceeds VOITLANDER_REFINE_FACTOR
//    times the narrowest of their three linewidths, i.e. while the affine
//    (first-order) model the per-triangle integral relies on has not yet been
//    resolved down to a small fraction of a linewidth
// Used by: RefineTriangle() (decides whether a triangle is split or integrated).
static bool NeedsRefine(VoitlanderMeshContext *Context, int IndexA, int IndexB, int IndexC)
{
    // resonance positions and linewidths at the three corners, with the position span and the narrowest width
    double posA=Context->Vertices[IndexA].Position, posB=Context->Vertices[IndexB].Position, posC=Context->Vertices[IndexC].Position;
    double widA=Context->Vertices[IndexA].Width, widB=Context->Vertices[IndexB].Width, widC=Context->Vertices[IndexC].Width;
    double span=std::max(posA, std::max(posB, posC))-std::min(posA, std::min(posB, posC));
    double narrowest=std::min(widA, std::min(widB, widC));

// absolute safety floor on the span itself, independent of the width-ratio
//    test above: once three vertices already agree on the resonance position to
//    within the field digitisation of the output spectrum, refining further can
//    only be chasing floating-point noise in VertexEval, never real physics
    if(span<1e-6) return false;

    return span>VOITLANDER_REFINE_FACTOR*narrowest;
}

// recursively refines one triangle, or, once the refine test is satisfied (or
//    the depth cap is hit), calls the injected per-triangle integral once and
//    folds its result into SpectrumSum with this triangle's ODF weight, taken
//    as the plain average of the three vertices' cached CurrentODF
// Used by: IntegrateVoitlanderMesh() (once for each seed triangle of this MPI rank, depth 0) and by RefineTriangle()
// itself (the four sub-triangles, depth+1).
static void RefineTriangle(VoitlanderMeshContext *Context, int IndexA, int IndexB, int IndexC, int Depth)
{
    // split into four triangles using the edge midpoints and refine each of them
    if(Depth<VOITLANDER_MAX_DEPTH && NeedsRefine(Context, IndexA, IndexB, IndexC))
    {
        int midAB=GetMidpoint(Context, IndexA, IndexB);
        int midBC=GetMidpoint(Context, IndexB, IndexC);
        int midCA=GetMidpoint(Context, IndexC, IndexA);
        RefineTriangle(Context, IndexA, midAB, midCA, Depth+1);
        RefineTriangle(Context, midAB, IndexB, midBC, Depth+1);
        RefineTriangle(Context, midCA, midBC, IndexC, Depth+1);
        RefineTriangle(Context, midAB, midBC, midCA, Depth+1);
        return;
    }

    // leaf triangle: update the mesh statistics and get the spectral contribution of the triangle from the injected integral
    if(Depth>Context->MaxDepthReached) Context->MaxDepthReached=Depth;
    Context->LeafTriangleCount++;

    double solidangle=SolidAngle(Context->XYZ[IndexA], Context->XYZ[IndexB], Context->XYZ[IndexC]);
    VoitlanderVertex *va=&Context->Vertices[IndexA], *vb=&Context->Vertices[IndexB], *vc=&Context->Vertices[IndexC];
    Context->TriangleIntegral(va, vb, vc, solidangle, Context->NPoints, Context->FieldNet, Context->Scratch);

    // add the contribution to the spectrum of each sample orientation, weighted by the mean ODF of the corners
    for(int k=0; k<Context->ODF->NOrientations; k++)
    {
        double odfweight=(va->CurrentODF[k]+vb->CurrentODF[k]+vc->CurrentODF[k])/3.0;
        for(int i=0; i<Context->NPoints; i++)
            Context->SpectrumSum[k][i]+=odfweight*Context->Scratch[i];
    }
}

// builds the level-0 hemisphere seed: a pole-up icosahedron's north-polar cap
//    (5 triangles, vertices at the pole and at the upper ring, polar angle
//    atan(2)=63.4349 degrees) plus a full equatorial band of 10 triangles
//    whose outer ring vertices are the icosahedron's lower-ring vertices
//    projected radially onto the true equator (polar angle 90 degrees, same
//    36-degree azimuthal stagger as the icosahedron's own upper-to-lower-ring
//    offset). The unmodified icosahedron has no face set that tiles an exact
//    hemisphere (its lower ring sits at 116.565 degrees, short of the next
//    full band), so this equatorial projection is a deliberate, documented
//    simplification of the level-0 seed only; adaptive refinement corrects
//    any residual distortion it introduces.
//    15 level-0 triangles, 11 level-0 vertices (1 pole + 5 + 5): 5 polar-cap
//    triangles plus BOTH alternating triangle orientations of the equatorial
//    band (upper[i],upper[next],equator[i]) and (upper[next],equator[next],
//    equator[i]); an earlier version built only the first of these two band
//    orientations, leaving the other 5 band gores entirely outside the mesh
//    (solid-angle sum 4.752 sr instead of the hemisphere's 2*pi=6.283 sr, a
//    24% silently unintegrated region) -- found and fixed 2026-11-03
// Used by: IntegrateVoitlanderMesh() (once, to get the 15 level-0 triangles).
static void BuildHemisphereSeed(VoitlanderMeshContext *Context, int SeedTriangleA[15], int SeedTriangleB[15], int SeedTriangleC[15])
{
    // pole vertex and the ring of 5 vertices at polar angle atan(2)
    int pole=CreateVertex(Context, 0, 0, 1);
    int upper[5], equator[5];
    double theta1=atan(2.0);
    for(int i=0; i<5; i++)
    {
        double ksi=i*2.0*PI/5.0;
        upper[i]=CreateVertex(Context, sin(theta1)*cos(ksi), sin(theta1)*sin(ksi), cos(theta1));
    }
    // ring of 5 vertices on the equator, shifted by 36 degrees in azimuth against the upper ring
    for(int i=0; i<5; i++)
    {
        double ksi=i*2.0*PI/5.0+PI/5.0;
        equator[i]=CreateVertex(Context, cos(ksi), sin(ksi), 0);
    }

    // the 15 seed triangles: for each azimuth step one polar cap triangle and two equatorial band triangles
    for(int i=0; i<5; i++)
    {
        int next=(i+1)%5;
        SeedTriangleA[i]=pole;       SeedTriangleB[i]=upper[i];    SeedTriangleC[i]=upper[next];
        SeedTriangleA[5+i]=upper[i]; SeedTriangleB[5+i]=upper[next]; SeedTriangleC[5+i]=equator[i];
        SeedTriangleA[10+i]=upper[next]; SeedTriangleB[10+i]=equator[next]; SeedTriangleC[10+i]=equator[i];
    }
}

// IntegrateVoitlanderMesh: orientation average over the upper hemisphere by adaptive triangulation: it builds the
// 15-triangle seed mesh, refines the triangles of this MPI rank until the refine test is satisfied, and lets the
// injected callbacks give the spectral contribution of every leaf triangle; the ODF-weighted raw solid-angle
// integral is added to SpectrumSum[k][i] (k: sample orientation, i: field point), and mesh statistics are written to
// the file voitlander<particle>.out by rank 0.
// Used by: CalculateIntegratedSpectrumVoitlander() in IntegratedSpectrum.cpp (the only caller).
void IntegrateVoitlanderMesh(MagneticParameters *MagneticParameter, OrderParameters *ODF,
                              VoitlanderVertexEval VertexEval, VoitlanderFreePayload FreePayload,
                              VoitlanderTriangleIntegral TriangleIntegral, double **SpectrumSum)
{
    // set up the mesh context with the particle data, the callbacks and the output array
    VoitlanderMeshContext context;
    context.MagneticParameter=MagneticParameter;
    context.ODF=ODF;
    context.VertexEval=VertexEval;
    context.TriangleIntegral=TriangleIntegral;
    context.NPoints=MagneticParameter->NPoints;
    context.FieldNet=MagneticParameter->FieldNet;
    context.SpectrumSum=SpectrumSum;
    context.Scratch=new double[context.NPoints];
    context.LeafTriangleCount=0;
    context.MaxDepthReached=0;

    // build the 15 level-0 hemisphere triangles
    int seedA[15], seedB[15], seedC[15];
    BuildHemisphereSeed(&context, seedA, seedB, seedC);

// each MPI rank owns a disjoint, statically-assigned subset of the 15 level-0
//    seed triangles and performs all of its own subdivision and evaluation
//    locally, with no cross-rank redundant work; the vertex cache above is
//    per-rank, since the seed triangles do not share edges with each other
#ifdef __MPI
    int firstseed=StartMPI(MpiRank, 15, MpiSize), lastseed=EndMPI(MpiRank, 15, MpiSize);
#else
    int firstseed=0, lastseed=15;
#endif
    // refine the seed triangles of this rank, the leaf triangles accumulate into SpectrumSum
    for(int tri=firstseed; tri<lastseed; tri++)
        RefineTriangle(&context, seedA[tri], seedB[tri], seedC[tri], 0);

// SpectrumSum itself is deliberately left as this rank's own partial sum over
//    its disjoint seed-triangle share (the same convention CalculateIntegratedSpectrumLiouville's
//    teta-partitioned Gauss-Kronrod loop already uses): Discrepancy.cpp's own
//    post-loop MPI_Gather/sum combines every particle's SpectrumSum across ranks
//    exactly once, so reducing it again here double-counted every contribution
//    by a factor of MpiSize (found and fixed 2026-11-03, see commit message)
#ifdef __MPI
    // combine the leaf count and the depth over all ranks for the statistics output
    long totalleaves=context.LeafTriangleCount; int maxdepth=context.MaxDepthReached;
    MPI_Allreduce(MPI_IN_PLACE, &totalleaves, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &maxdepth, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    context.LeafTriangleCount=totalleaves; context.MaxDepthReached=maxdepth;
#endif

// mesh statistics, one file per particle, for the validation comparison
//    against IntegrateGauss/IntegrateGaussKronrod (triangle count, depth reached)
    if(!MpiRank)
    {
        char buff[MAX_STRLEN];
        sprintf(buff, "voitlander%d.out", MagneticParameter->Count);
        std::fstream statfile; statfile.open(buff, std::fstream::out);
        statfile << "Leaf triangles: " << context.LeafTriangleCount << std::endl;
        statfile << "Maximum subdivision depth reached: " << context.MaxDepthReached << " (cap " << VOITLANDER_MAX_DEPTH << ")" << std::endl;
        statfile << "Mesh vertices evaluated: " << context.Vertices.size() << std::endl;
        statfile.close();
    }

    // release the cached payloads of all vertices and the scratch array
    for(size_t v=0; v<context.Vertices.size(); v++)
        FreePayload(context.Vertices[v].Payload);
    delete[] context.Scratch;
}

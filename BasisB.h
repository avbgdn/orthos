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

// BasisB.h: declares BasisVectorB, the record that labels one vector of the basis in which the
// stochastic Liouville equation (SLE) is written, and the three routines that build or read the
// list of such vectors (implemented in BasisB.cpp).
// Used by: BasisB.cpp and LiouvillianB.cpp (which include it), and through BiradicalODF.h (which
// includes it) by every source file that includes BiradicalODF.h, in particular
// LiouvilleMatrixElementB.cpp/.h, where BasisVectorB is the argument type of the matrix element functions.

#ifndef BASISB_H
#define BASISB_H

// One basis vector: an orientational function times a product of spin operators.
// L, K, M: rank and the two projection indices of the orientational function (the symmetry flags
//   FduReal, FduEven and FduAxial restrict K in GenerateBasisB; for a director parallel to the
//   field M is tied to the nuclear coherence orders through pI1+pI2=M).
// N: index of the exchange-coupling (J) diffusion variable, its diagonal Liouvillian element is
//   JDiffusionRate*N*N (see MatrixElementJDiffusion); N is 0 for a monoradical.
// Q: conformer index.
// pS1,qS1 and pI1,qI1 (electron 1 and its nucleus), pS2,qS2 and pI2,qI2 (electron 2 and its nucleus):
//   each spin carries the operator |m_ket><m_bra| with p=m_ket-m_bra (coherence order) and
//   q=m_ket+m_bra (see SpinIndexB); for a monoradical the second electron and nucleus are placeholders
//   (pS2=0, qS2=1, pI2=0, qI2=0).
struct BasisVectorB
{
	int L, K, M, N, Q, pS1, qS1, pI1, qI1, pS2, qS2, pI2, qI2;
};

// Reads the basis vectors from a pruned-basis file and returns how many were read.
// Used by: CalculateLiouvillianB in LiouvillianB.cpp.
int LoadBasisB(char *filename, BasisVectorB *basis);
// Generates the basis for a biradical group into basis (and an identical copy in basis1) and returns its dimension.
// Used by: CalculateLiouvillianB in LiouvillianB.cpp.
int GenerateBasisB(int LEmax, int LOmax, int Kmax, int Mmax, int Nmax, int Ix2, int NConformers, BasisVectorB *basis, BasisVectorB *basis1);
// Generates the M-symmetrised monoradical basis (M>=0) into basis (and an identical copy in basis1) and returns its dimension.
// Used by: CalculateLiouvillianB in LiouvillianB.cpp.
int GenerateBasisM(int LEmax, int LOmax, int Kmax, int Mmax, int Nmax, int Ix2, int NConformers, BasisVectorB *basis, BasisVectorB *basis1);

#endif
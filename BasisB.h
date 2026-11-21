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

#ifndef BASISB_H
#define BASISB_H

struct BasisVectorB
{
	int L, K, M, N, Q, pS1, qS1, pI1, qI1, pS2, qS2, pI2, qI2;
};

int LoadBasisB(char *filename, BasisVectorB *basis);
int GenerateBasisB(int LEmax, int LOmax, int Kmax, int Mmax, int Nmax, int Ix2, int NConformers, BasisVectorB *basis, BasisVectorB *basis1);
int GenerateBasisM(int LEmax, int LOmax, int Kmax, int Mmax, int Nmax, int Ix2, int NConformers, BasisVectorB *basis, BasisVectorB *basis1);

#endif
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

// This file contains functions for generating Basis Sets for SLE calculations. 

extern bool FduReal;
extern bool FduEven;
extern bool FduAxial;
extern bool DirectorParallel;
extern std::fstream timing;

extern int MpiRank, MpiSize, ParticleCount;

// Load basis set from file filename, for example in case of calculated in pruned basis.
// Basis file structure should be the following.
// The file can contain any infromation, but in some place should contain the line "Vectors:", 
//  after which the basis vectors for SLE calculation are listed, line by line, in the form: 
//  L  K  M  N  Q pS1 qS1 pS2 qS2 pI1 qI1 pI2 qI2  (thirteen numbers separated with spaces),
//  until the end of the file is reached.
// The pruning subroutine produces basis files in the appropriate form.
int LoadBasisB(char *filename, BasisVectorB *basis)
{
	int L, K, M, N, Q, pS1, pS2, qS1, qS2, pI1, pI2, qI1, qI2;
	int nbasis=0;
	std::fstream inputfile;
	inputfile.open(filename, std::fstream::in);
	char buffer[MAX_STRLEN], buffer1[MAX_STRLEN];
        char buff[MAX_STRLEN];
        sprintf(buff, "basis%d_%d.out", ParticleCount, MpiRank);
	FILE *basisfile;
        basisfile=fopen(buff, "w");
	fprintf(basisfile, "Started basis readout from %s. \n      |  L  K  M  N  Q pS1 qS1 pS2 qS2 pI1 qI1 pI2 qI2 >\n" , filename);

	while(inputfile.getline(buffer, MAX_STRLEN))
		if(!strcmp(buffer, "Vectors:") || !strcmp(buffer, "Vectors:\r"))
		{
			// начать считывать базисные векторы
			inputfile.getline(buffer, MAX_STRLEN);
			while (!inputfile.eof())
			{
				sscanf(buffer, "%d %d %d %d %d %d %d %d %d %d %d %d %d", &L, &K, &M, &N, &Q, &pS1, &qS1, &pS2, &qS2, &pI1, &qI1, &pI2, &qI2);
				basis[nbasis].L=L;
				basis[nbasis].K=K;
				basis[nbasis].M=M;
				basis[nbasis].N=N;
				basis[nbasis].Q=Q;
				basis[nbasis].pS1=pS1;
				basis[nbasis].qS1=qS1;
				basis[nbasis].pI1=pI1;
				basis[nbasis].qI1=qI1;
				basis[nbasis].pS2=pS2;
				basis[nbasis].qS2=qS2;
				basis[nbasis].pI2=pI2;
				basis[nbasis].qI2=qI2;	


				fprintf(basisfile, "%5d | % d % d % d % d % d % d % d  % d  % d  % d  % d  % d  % d  >\n", nbasis+1, L, K, M, N, Q, pS1, qS1,pS2, qS2, pI1, qI1 , pI2, qI2 );
				nbasis++;

				inputfile.getline(buffer, MAX_STRLEN);
				sscanf(buffer, "%s", buffer1);
				if (!strcmp(buffer1, "Dimension")) break;
			}
		}
	fprintf(basisfile, "Dimension of basis is %d\n", nbasis);
        if(!MpiRank)
	timing << "Finished basis read.";
        if(!MpiRank)
	timing << "Dimension of basis is " << nbasis << endl;
	fclose(basisfile);
	return nbasis;

}

// The subroutine for generating basis set for SLE calculation for biradical.
// Some symmetries of magnetic tensors reduce the basis dimension by the following rules:
//
// Let Fd(l, m) (with l=0,2; m=-l..+l) be the spherical tensor representation of the magnetic interaction tensor in the frame of reference
// of particle orientation (ordering potential frame) for ordered systems, and of rotation diffusion tensor for isotropic systems.
// 1. If for any interaction, Fd(l, m) are real, the FduReal flag is set true, 
//     and only basis vectors with jK (see Schneider, Freed, 1989) are considered.
// 2. If for any interaction, Fd(l, m) are non-zero only for even m, the FduEven flag is set true
//     and only basis vectors with even K are considered.
// 3. If for any interaction, Fd(l, m) are non-zero only for m=0, the FduAxial flag is set true
//     and only basis vectors with even L and K=0 are considered.
int GenerateBasisB(int LEmax, int LOmax, int Kmax, int Mmax, int Nmax, int Ix2, int NConformers, BasisVectorB *basis, BasisVectorB *basis1)
{
	FILE *basisfile;
        char buff[MAX_STRLEN];
        sprintf(buff, "basis%d_%d.out", ParticleCount, MpiRank);
	basisfile=fopen(buff, "w");
	int nbasis=0;
        if(!MpiRank)
	timing << "Started basis generation. \n";
	
	if(FduReal)
		fprintf(basisfile, "FDU Real!\n");
	if(FduEven)
		fprintf(basisfile, "FDU Even!\n");
	if(FduAxial)
		fprintf(basisfile, "FDU Axial!\n");

	fprintf(basisfile, "Started basis generation. \n      |  L  K  M  N  Q  pS1 qS1 pS2 qS2 pI1 qI1 pI2 qI2 >\n" );


	int K0, M0;
	K0=0; M0=0;
for(int Q=0; Q<NConformers; Q+=1)
{
	for(int L=0; L<=max(LEmax, LOmax); L+=1)
	{
		if(L%2)
		{
			if(L>LOmax) continue;
		}
		else
			if(L>LEmax) continue;

		for(int K=-min(L, Kmax); K<=min(L, Kmax); K+=1)
		{
			if(FduReal)
			{
				if (K<0) continue;
				if ((K==0)&&(L%2)) continue;		// only jK=+1!
			}
			if(FduEven)
			{
				if (K%2) continue;					// only even K;
			}
			if(FduAxial)
			{
				if (L%2) continue;
				if (K!=0) continue;					// only even L and K=0;
			}

			for(int M=-min(L, Mmax); M<=min(L, Mmax); M+=1)
			  for(int N=0; N<=Nmax; N+=2)
				for(int pS1=0; pS1<=1; pS1++)
					for(int qS1=abs(pS1)-1; qS1<=1-abs(pS1); qS1+=2)
//						for(int pS2=0; pS2<=0; pS2+=2)   // only 0
						for(int pS2=1-abs(pS1); pS2<=1-abs(pS1); pS2+=2)   // if pS1=0 then pS2=+1; if pS1=+1 then pS2=0;
//							for(int qS2=1; qS2<=1; qS2+=2)
							for(int qS2=abs(pS2)-1; qS2<=1-abs(pS2); qS2+=2)
								for(int pI1=-2; pI1<=2; pI1++)
									for(int qI1=abs(pI1)-Ix2; qI1<=Ix2-abs(pI1); qI1+=2)
//									for(int qI1=0; qI1<=0; qI1+=2)
										for(int pI2=-2; pI2<=2; pI2++)
										{
											if (DirectorParallel&&(pI1+pI2 != M)) continue;
											for(int qI2=abs(pI2)-Ix2; qI2<=Ix2-abs(pI2); qI2+=2)
//											for(int qI2=0; qI2<=0; qI2+=2)	// only 0
							{
								basis[nbasis].L=L;
								basis[nbasis].K=K;
								basis[nbasis].M=M;
								basis[nbasis].N=N;
								basis[nbasis].pS1=pS1;
								basis[nbasis].qS1=qS1;
								basis[nbasis].pI1=pI1;
								basis[nbasis].qI1=qI1;
								basis[nbasis].pS2=pS2;
								basis[nbasis].qS2=qS2;
								basis[nbasis].pI2=pI2;
								basis[nbasis].qI2=qI2;
								basis[nbasis].Q=Q;

								basis1[nbasis].L=L;
								basis1[nbasis].K=K;
								basis1[nbasis].M=M;
								basis1[nbasis].N=N;
								basis1[nbasis].pS1=pS1;
								basis1[nbasis].qS1=qS1;
								basis1[nbasis].pI1=pI1;
								basis1[nbasis].qI1=qI1;
								basis1[nbasis].pS2=pS2;
								basis1[nbasis].qS2=qS2;
								basis1[nbasis].pI2=pI2;
								basis1[nbasis].qI2=qI2;
								basis1[nbasis].Q=Q;


								fprintf(basisfile, "%5d | % d % d % d % d % d % d  % d  % d  % d  % d  % d  % d  % d  >\n", nbasis+1, L, K, M, N, Q, pS1, qS1,pS2, qS2, pI1, qI1 , pI2, qI2 );
								nbasis++;
							}
										}
		}
	}
}
	fprintf(basisfile, "Dimension of basis is %d\n", nbasis);

	timing << "Finished basis generation.";
	timing << "Dimension of basis is " << nbasis << endl;
	fclose(basisfile);
	return nbasis;
}

// Basis generation subroutine for monoradical
// The same symmetry selection rules as for biradical molecule apply.
//
int GenerateBasisM(int LEmax, int LOmax, int Kmax, int Mmax, int Nmax, int Ix2, int NConformers, BasisVectorB *basis, BasisVectorB *basis1)
{
	FILE *basisfile;
        char buff[MAX_STRLEN];
        sprintf(buff, "basis%d_%d.out", ParticleCount, MpiRank);
	basisfile=fopen(buff, "w");
	int nbasis=0;
        if(!MpiRank)
	timing << "Started basis generation. \n";
	
	if(FduReal)
		fprintf(basisfile, "FDU Real!\n");
	if(FduEven)
		fprintf(basisfile, "FDU Even!\n");
	if(FduAxial)
		fprintf(basisfile, "FDU Axial!\n");

	fprintf(basisfile, "Started basis generation. \n      |  L  K  M  N  Q  pS1 qS1 pS2 qS2 pI1 qI1 pI2 qI2 >\n" );


	int K0, M0;
	int N=0;
	K0=0; M0=0;
for(int Q=0; Q<NConformers; Q++)
{
	for(int L=0; L<=max(LEmax, LOmax); L+=1)
	{
		if(L%2)
		{
			if(L>LOmax) continue;
		}
		else
			if(L>LEmax) continue;

		for(int K=-min(L, Kmax); K<=min(L, Kmax); K+=1)
		{
			if(FduReal)
			{
				if (K<0) continue;
				if ((K==0)&&(L%2)) continue;		// only jK=+1!
			}
			if(FduEven)
			{
				if (K%2) continue;					// only even K;
			}
			if(FduAxial)
			{
				if (L%2) continue;
				if (K!=0) continue;					// only even L and K=0;
			}

// AB: for monoradicals in M-symmetrized basis replaced this:
//			for(int M=-min(L, Mmax); M<=min(L, Mmax); M+=1)
// for this:
			for(int M=0; M<=min(L, Mmax); M+=1)
				for(int pS1=1; pS1<=1; pS1++)
					for(int qS1=abs(pS1)-1; qS1<=1-abs(pS1); qS1+=2)
						for(int pS2=0; pS2<=0; pS2+=2)   // only 0
//						for(int pS2=1-abs(pS1); pS2<=1-abs(pS1); pS2+=2)   // if pS1=0 then pS2=+1; if pS1=+1 then pS2=0;
							for(int qS2=1; qS2<=1; qS2+=2)
//							for(int qS2=abs(pS2)-1; qS2<=1-abs(pS2); qS2+=2)
								for(int pI1=-2; pI1<=2; pI1++)
									for(int qI1=abs(pI1)-Ix2; qI1<=Ix2-abs(pI1); qI1+=2)
//									for(int qI1=0; qI1<=0; qI1+=2)
										for(int pI2=0; pI2<=0; pI2++)
										{
											if (DirectorParallel&&(pI1+pI2 != M)) continue;

											// the M-symmetrised vector is a(|r>+(-1)^(L+M)|Sr>), so a self-conjugate vector with odd L
											// carries the phase -1 and the combination vanishes; such a null symmetrised vector would
											// put a zero row and column into the Liouvillian and must be kept out of the basis
											if((M==0)&&(pI1==0)&&(pI2==0)&&(L%2)) continue;

											// at M=0 the symmetrised vector built on -pI1 is (-1)^L times the one built on +pI1, so
											// admitting both signs would have the same vector entered twice and the basis would be
											// linearly dependent; lbasix.f of the Freed code restricts pI1 likewise at M=0
											if((M==0)&&(pI1<0)) continue;
//											for(int qI2=abs(pI2)-Ix2; qI2<=Ix2-abs(pI2); qI2+=2)
											for(int qI2=0; qI2<=0; qI2+=2)	// only 0
							{
								basis[nbasis].L=L;
								basis[nbasis].K=K;
								basis[nbasis].M=M;
								basis[nbasis].N=0;
								basis[nbasis].pS1=pS1;
								basis[nbasis].qS1=qS1;
								basis[nbasis].pI1=pI1;
								basis[nbasis].qI1=qI1;
								basis[nbasis].pS2=pS2;
								basis[nbasis].qS2=qS2;
								basis[nbasis].pI2=pI2;
								basis[nbasis].qI2=qI2;
								basis[nbasis].Q=Q;

								basis1[nbasis].L=L;
								basis1[nbasis].K=K;
								basis1[nbasis].M=M;
								basis1[nbasis].N=0;
								basis1[nbasis].pS1=pS1;
								basis1[nbasis].qS1=qS1;
								basis1[nbasis].pI1=pI1;
								basis1[nbasis].qI1=qI1;
								basis1[nbasis].pS2=pS2;
								basis1[nbasis].qS2=qS2;
								basis1[nbasis].pI2=pI2;
								basis1[nbasis].qI2=qI2;
								basis1[nbasis].Q=Q;


								fprintf(basisfile, "%5d | % d % d % d % d % d % d  % d  % d  % d  % d  % d  % d  % d  >\n", nbasis+1, L, K, M, N, Q, pS1, qS1,pS2, qS2, pI1, qI1 , pI2, qI2 );
								nbasis++;
							}
										}
		}
	}
}
	fprintf(basisfile, "Dimension of basis is %d\n", nbasis);

        if(!MpiRank)
	timing << "Finished basis generation.";
        if(!MpiRank)
	timing << "Dimension of basis is " << nbasis << endl;
	fclose(basisfile);
	return nbasis;
}
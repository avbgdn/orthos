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
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include "LiouvilleMatrixElementB.h"

double *Wigner3jSymbols;

// Coefficient A(j) of the three-term recursion over j1=j for the symbol (j j2 j3; m1 m2 m3)
static double recursion_A(double j, double j2, double j3, double m1)
{
	return sqrt((j*j-(j2-j3)*(j2-j3))*((j2+j3+1)*(j2+j3+1)-j*j)*(j*j-m1*m1));
}

// Coefficient B(j) of the three-term recursion over j1=j for the symbol (j j2 j3; m1 m2 m3)
static double recursion_B(double j, double j2, double j3, double m1, double m2, double m3)
{
	return -(2*j+1)*((j2*(j2+1)-j3*(j3+1))*m1-j*(j+1)*(m3-m2));
}

// Computes (j L2 L3; -M2-M3 M2 M3) for all j between jmin and L2+L3 into f[0...L2+L3-jmin], jmin=max(|L2-L3|,|M2+M3|).
// The method is the three-term recursion over j of K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975):
//   j A(j+1) f(j+1) + B(j) f(j) + (j+1) A(j) f(j-1) = 0
// The recursion is run upwards from jmin as long as |f| grows, and downwards from L2+L3 to the same point. The two
// halves are joined there, the normalisation is sum_j (2j+1) f(j)^2=1 and the sign is fixed by the sign of f(L2+L3),
// which is (-1)^(L2-L3-M1).
static void recursion_3j(int L2, int L3, int M2, int M3, int jmin, int jmax, double *f, double *g)
{
	double j2=(double)L2;
	double j3=(double)L3;
	double m1=(double)(-M2-M3);
	double m2=(double)M2;
	double m3=(double)M3;
	int n=jmax-jmin+1;
	int join=0;

	// Upward recursion from jmin, stopped at the first local maximum of the absolute value
	f[0]=1.0;
	if(jmin>0)
	{
		for(join=0; join<n-1; join++)
		{
			double j=(double)(jmin+join);
			double below=0.0;
			if(join>0) below=(j+1)*recursion_A(j, j2, j3, m1)*f[join-1];
			f[join+1]=-(recursion_B(j, j2, j3, m1, m2, m3)*f[join]+below)/(j*recursion_A(j+1, j2, j3, m1));
			if(fabs(f[join+1])<=fabs(f[join])) break;
			if(fabs(f[join+1])>1e50)
				for(int k=0; k<=join+1; k++) f[k]*=1e-50;
		}
	}

	// Downward recursion from the maximum j, down to the joining point
	if(join<n-1)
	{
		g[n-1]=1.0;
		for(int k=n-1; k>join; k--)
		{
			double j=(double)(jmin+k);
			double above=0.0;
			if(k<n-1) above=j*recursion_A(j+1, j2, j3, m1)*g[k+1];
			g[k-1]=-(above+recursion_B(j, j2, j3, m1, m2, m3)*g[k])/((j+1)*recursion_A(j, j2, j3, m1));
			if(fabs(g[k-1])>1e50)
				for(int p=k-1; p<n; p++) g[p]*=1e-50;
		}

		// Scale the downward solution to meet the upward one
		double scale=f[join]/g[join];
		for(int k=join+1; k<n; k++) f[k]=g[k]*scale;
	}

	// Normalise and fix the overall sign
	double norm=0.0;
	for(int k=0; k<n; k++) norm+=(2.0*(jmin+k)+1.0)*f[k]*f[k];
	norm=1.0/sqrt(norm);
	if(((L2-L3+M2+M3)%2!=0)!=(f[n-1]<0)) norm=-norm;
	for(int k=0; k<n; k++) f[k]*=norm;
}

unsigned long int index_3j(unsigned long int L1, unsigned long int L2, unsigned long int L3, unsigned long int M2, unsigned long int M3)
{
	return L2 + M2 + (L3 + (L3*(L3 - 1))/2)*(2*L2 + 1) + M3*(2*L2 + 1) + (L2*(3*L2 - 1)*(L2 + 1)*(L2 + 2))/12 + (L1*(6*L1 - 1)*(L1 + 1)*(L1 + 2)*(L1 + 3))/120;
}

double Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3)
{
	if (M1+M2+M3) return 0;
	if (abs(M1)>L1) return 0;
	if (abs(M2)>L2) return 0;
	if (abs(M3)>L3) return 0;

	if (L1<0) return 0;
	if (L2<0) return 0;
	if (L3<0) return 0;

	double result, sign;
	if ((L1+L2+L3)%2) sign=-1;
	else			  sign=1;
	if((L1>=L2)&&(L2>=L3))
	{
		// l1 l2 l3	even
		if (M3>=0) return      Wigner3jSymbols[index_3j(L1, L2, L3, M2, M3)];
		else	   return sign*Wigner3jSymbols[index_3j(L1, L2, L3,-M2,-M3)];
	}
	if((L1<L2)&&(L2<L3))
	{
		// l3 l2 l1	odd
		if (M1>=0) return sign*Wigner3jSymbols[index_3j(L3, L2, L1, M2, M1)];
		else	   return      Wigner3jSymbols[index_3j(L3, L2, L1,-M2,-M1)];
	}
	if((L1>=L2)&&(L2<L3)&&(L1>=L3))
	{
		// l1 l3 l2	odd
		if (M2>=0) return sign*Wigner3jSymbols[index_3j(L1, L3, L2, M3, M2)];
		else	   return      Wigner3jSymbols[index_3j(L1, L3, L2,-M3,-M2)];
	}
	if((L1>=L2)&&(L2<L3)&&(L1<L3))
	{
		// l3 l1 l2 even
		if (M2>=0) return      Wigner3jSymbols[index_3j(L3, L1, L2, M1, M2)];
		else	   return sign*Wigner3jSymbols[index_3j(L3, L1, L2,-M1,-M2)];
	}
	if((L1<L2)&&(L2>=L3)&&(L1>=L3))
	{
		// l2 l1 l3 odd
		if (M3>=0) return sign*Wigner3jSymbols[index_3j(L2, L1, L3, M1, M3)];
		else	   return      Wigner3jSymbols[index_3j(L2, L1, L3,-M1,-M3)];
	}
	if((L1<L2)&&(L2>=L3)&&(L1<L3))
	{
		// l2 l3 l1 even
//		printf("%d %d %d %d %d %d \n", L1, L2, L3, M1, M2, M3);
		if (M1>=0) return      Wigner3jSymbols[index_3j(L2, L3, L1, M3, M1)];
		else	   return sign*Wigner3jSymbols[index_3j(L2, L3, L1,-M3,-M1)];
	}

	return result;
}
void Calculate_Wigner3jm(int Lmax)
{
	Wigner3jSymbols=(double*)calloc(index_3j(Lmax+1, 0, 0, 0, 0), sizeof(double));

	if(Wigner3jSymbols==NULL)
	{
		printf("Allocation fault!, Lmax=%d, index_3j=%lu\n", Lmax, index_3j(Lmax+1, 0, 0, 0, 0));
		return;
	}

	// Work arrays for the recursion, long enough for any j between 0 and 2*Lmax
	double *forward=(double*)malloc((2*Lmax+1)*sizeof(double));
	double *backward=(double*)malloc((2*Lmax+1)*sizeof(double));

	// Symbols that violate the triangle or the projection rules stay at the zero set by calloc
	int L2, L3, M2, M3;
	for(L2=0; L2<=Lmax; L2++)
		for(L3=0; L3<=L2; L3++)
			for(M3=0; M3<=L3; M3++)
				for(M2=-L2; M2<=L2; M2++)
				{
					int jmin=std::max(L2-L3, abs(M2+M3));
					int jmax=L2+L3;
					if(jmin>jmax) continue;
					recursion_3j(L2, L3, M2, M3, jmin, jmax, forward, backward);
					for(int L1=std::max(jmin, L2); (L1<=jmax)&&(L1<=Lmax); L1++)
						Wigner3jSymbols[index_3j(L1, L2, L3, M2, M3)]=forward[L1-jmin];
				}

	free(forward);
	free(backward);
}


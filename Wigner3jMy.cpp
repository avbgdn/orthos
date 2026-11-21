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
#include "LiouvilleMatrixElementB.h"

double *Wigner3jSymbols;

double _Wigner3jMy(int L1, int L2, int L3, int M1, int M2, int M3)
{
	int errflag;
	double l1=(double)L1;
	double l2=(double)L2;
	double l3=(double)L3;
	double m1=(double)M1;
	double m3=(double)M3;
	double m2min, m2max;

	if(M1+M2+M3) return 0;
	double thrcof[MAX_WIGNER];
	
	Wigner3jm(l1, l2, l3, m1, m2min, m2max, thrcof, MAX_WIGNER, errflag);
	if(errflag==6) printf("Too small thrcof in Wigner3JM calculation. Consider recompiling\n");
	if(errflag) return 0;
	if((-m1-m2max<=m3)&&(-m1-m2min>=m3))
		return thrcof[(int)(-m1-m3-m2min)];
	else return 0;
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
	Wigner3jSymbols=(double*)malloc(index_3j(Lmax+1, 0, 0, 0, 0)*sizeof(double));

	if(Wigner3jSymbols==NULL)
	{
		printf("Allocation fault!, Lmax=%d, index_3j=%lu\n", Lmax, index_3j(Lmax+1, 0, 0, 0, 0));
		return;
	}
	int L1, L2, L3, M2, M3;
	int index=0;
	for(L1=0; L1<Lmax+1; L1++)
		for(L2=0; L2<=L1; L2++)
			for(L3=0; L3<=L2; L3++)
				for(M3=0; M3<=L3; M3++)
					for(M2=-L2; M2<=L2; M2++)
					{
//						int index=L2 + M2 + (L3 + (L3*(L3 - 1))/2)*(2*L2 + 1) + M3*(2*L2 + 1) + (L2*(3*L2 - 1)*(L2 + 1)*(L2 + 2))/12 + (L1*(6*L1 - 1)*(L1 + 1)*(L1 + 2)*(L1 + 3))/120;
						Wigner3jSymbols[index]=_Wigner3jMy(L1, L2, L3, -M2-M3, M2, M3);
						index++;
					}

}



/***********************************************************************
*                                                                      *
*     Program to calculate the first kind modified Bessel function     *
*  of integer order N, for any REAL X, using the function BESSI(N,X).  *
*                                                                      *
* -------------------------------------------------------------------- *
*    SAMPLE RUN:                                                       *
*                                                                      *
*    (Calculate Bessel function for N=2, X=0.75).                      *
*                                                                      *
*    Bessel function of order 2 for X =  0.7500:                       *
*                                                                      *
*         Y = 0.073667                                                 *
*                                                                      *
* -------------------------------------------------------------------- *
*    Reference: From Numath Library By Tuan Dang Trong in Fortran 77.  *
*                                                                      *
*                               C++ Release 1.1 By J-P Moreau, Paris.  *
*                                        (www.jpmoreau.fr)             *
*                                                                      * 
*    Version 1.1: corected value of P4 in BESSIO (P4=1.2067492 and not *
*                 1.2067429) Aug. 2011.                                *
***********************************************************************/
// File purpose: modified Bessel functions of the first kind, I0, I1 and I_N, adapted from the program described
// in the banner above.
// Used by: nothing in this code base; the file is compiled and linked (it is listed in the Makefile) but no other
// source file calls these functions or declares them.
#include <stdio.h>
#include <math.h>


  // forward declarations of the helper functions defined below
  double BESSI0(double X);
  double BESSI1(double X);

// ---------------------------------------------------------------------
  // BESSI: modified Bessel function I_N(X) of integer order N for real X. N=0 and N=1 are handed to BESSI0 and BESSI1;
  // for higher N the downward (Miller) recurrence is used, normalised with I0(X). Note: the code below always uses
  // the downward recurrence for N>=2; the forward recursion for X>N mentioned in the text below is not implemented.
  // Used by: nothing in this code base.
  double BESSI(int N, double X) {
/*----------------------------------------------------------------------
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
------------------------------------------------------------------------*/

      // accuracy parameter of the starting index of the Miller recurrence and the limits used to rescale against overflow
      int IACC = 40; 
	  double BIGNO = 1e10, BIGNI = 1e-10;
      double TOX, BIM, BI, BIP, BSI;
      int J, M;

      // orders 0 and 1 and the point X=0 are handled directly
      if (N==0)  return (BESSI0(X));
      if (N==1)  return (BESSI1(X));
      if (X==0.0) return 0.0;

      // start the downward recurrence at a high index M with arbitrary small starting values
      TOX = 2.0/X;
      BIP = 0.0;
      BI  = 1.0;
      BSI = 0.0;
      M = (int) (2*((N+floor(sqrt((double)IACC*N)))));
      // downward recurrence, rescaling all values when they grow too large, and keeping the value at index N
      for (J = M; J>0; J--) {
        BIM = BIP+J*TOX*BI;
        BIP = BI;
        BI  = BIM;
        if (fabs(BI) > BIGNO) {
          BI  = BI*BIGNI;
          BIP = BIP*BIGNI;
          BSI = BSI*BIGNI;
        }
        if (J==N)  BSI = BIP;
      }
      // normalise the kept value using the accurately known I0(X) against the recurrence value at index 0
      return (BSI*BESSI0(X)/BI);
  }

// ----------------------------------------------------------------------
//  Auxiliary Bessel functions for N=0, N=1
  // BESSI0: I0(X) from polynomial approximations: a polynomial in (X/3.75)^2 for |X|<3.75 and exp(|X|)/sqrt(|X|)
  // times a polynomial in 3.75/|X| otherwise.
  // Used by: BESSI() in this file (for N=0 and for normalising the recurrence).
  double BESSI0(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067492;
      P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
      Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
      Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
      Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
      // small argument: polynomial in (X/3.75)^2
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      }
      else {
        // large argument: asymptotic form exp(|X|)/sqrt(|X|) times a polynomial in 3.75/|X|
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

// ---------------------------------------------------------------------
  // BESSI1: I1(X) from polynomial approximations of the same kind as in BESSI0(). Note: for X<=-3.75 the large-argument
  // branch returns a positive value, i.e. the odd symmetry I1(-X)=-I1(X) is not restored there.
  // Used by: BESSI() in this file (for N=1).
  double BESSI1(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=0.5; P2=0.87890594; P3=0.51498869; P4=0.15084934;
      P5=0.2658733e-1; P6=0.301532e-2; P7=0.32411e-3;
      Q1=0.39894228; Q2=-0.3988024e-1; Q3=-0.362018e-2;
      Q4=0.163801e-2; Q5=-0.1031555e-1; Q6=0.2282967e-1;
      Q7=-0.2895312e-1; Q8=0.1787654e-1; Q9=-0.420059e-2;
      // small argument: X times a polynomial in (X/3.75)^2
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
      }
      else {
        // large argument: asymptotic form exp(|X|)/sqrt(|X|) times a polynomial in 3.75/|X|
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

// ---------------------------------------------------------------------

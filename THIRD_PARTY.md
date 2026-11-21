# Third-party components

Orthos itself is GNU GPL v3 (see the headers of the orthos source files; the full text is in LICENSE).

| Component | Files | Source | Licence as stated in the files |
|---|---|---|---|
| Eigen 3.2.5 | `Eigen/`, `unsupported/Eigen/KroneckerProduct` | https://eigen.tuxfamily.org | Primarily MPL 2.0; some files carry BSD or LGPL terms. The upstream licence files of release 3.2.5 (COPYING.MPL2, COPYING.BSD, COPYING.LGPL, COPYING.GPL, COPYING.MINPACK, COPYING.README) are in `Eigen/`. They were not present in the tree this repository was made from and were added from the 3.2.5 tag of the Eigen repository. |
| Faddeeva package | `faddeeva.cc`, `Faddeeva.hh` | http://ab-initio.mit.edu/Faddeeva | MIT licence, Copyright (c) 2012 Massachusetts Institute of Technology (full text in the file headers) |
| Kronrod quadrature | `kronrod.cpp`, `kronrod.h` | C version by John Burkardt of the FORTRAN77 code of R. Piessens and M. Branders (Math. Comp. 28 (1974) 135-139) | "This code is distributed under the GNU LGPL license." (header of each routine) |
| NL2SOL | `NL2SOL1.for`, `Nl2def1.for` | NL2SOL version 2.2 (Dennis, Gay, Welsch), ACM TOMS Algorithm 573, from netlib | The file headers state no licence ("LATEST REVISION - 03/15/90 (JRD)" only). **Licence status unclear; needs review before publication.** |

System libraries used and not bundled: FFTW3 (GPL), NLopt (MIT/LGPL), libgfortran/libgomp (GPL with runtime exception).

The Wigner 3j routine in `Wigner3jMy.cpp` is original orthos code (GPL v3) implementing the three-term recursion of K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975), written from the published method. The earlier routine adapted from Matpack (`Wigner3jm.cpp`, `Wigner3jm.hpp`, `simError.cpp`, `simError.h`) was removed because the Matpack licence does not permit distribution of modified versions.

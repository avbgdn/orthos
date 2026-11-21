# orthos

Orthos simulates continuous-wave EPR spectra of nitroxide mono- and biradicals in slowly tumbling, partially ordered environments. It builds the stochastic Liouville equation (SLE) in a basis of orientational, electron-spin and nuclear-spin operator products, tridiagonalises the Liouvillian by the Lanczos algorithm, and returns the field-swept absorption derivative. Rotational diffusion (Brownian, free or jump), an orientational potential, discrete jumps, librations, g and hyperfine tensors, Gaussian and Lorentzian line broadening and dipolar/exchange coupling in biradicals are supported. Parameters can be fitted to experimental spectra with the NL2SOL least-squares minimiser (NLopt is also linked).

The original program and its documentation: http://orthos-epr.sourceforge.net/ . Copyright (C) 2015-2016 Alexey Bogdanov, GNU GPL version 3 (see the file headers).

This repository holds the **unsaturated** program (orthos 1.3). A version with microwave saturation is in the sibling repository `orthos_sat`.

## Change in this version: director tilt

For a tilted director (angle psi between the director and the field, ODF/orientation not along the field) the M-symmetrised monoradical basis had two defects:

* in `MatrixElementBSM` the symmetrised vector was treated as self-conjugate whenever M=0; it is self-conjugate only when M=0 and pI1=0 (the index pI1 is flipped together with M by the symmetrisation operator). The two conditions coincide at psi=0, where the rule pI1+pI2=M ties pI1 to M;
* `GenerateBasisM` admitted the odd-L self-conjugate vector, whose symmetrised combination vanishes, and admitted both signs of pI1 at M=0, where the symmetrised vector at -pI1 is (-1)^L times the one at +pI1, so one vector entered the basis twice (compare lbasix.f of the Freed code).

Both are fixed in `LiouvilleMatrixElementB.cpp` and `BasisB.cpp`. At psi=0 the spectra are unchanged (bit identical to the version without the fix). At psi=45 degrees with an ordering potential (c20=+2) the normalised RMS difference to the Freed code (frnmin5) is 0.0036 after the fix, against 0.140 before (shift-compensated, amplitude-normalised comparison); without a potential and with an isotropic ODF the spectrum must not depend on psi, and the RMS difference between psi=0 and psi=45 falls from 0.093 to 0.0028. These numbers were measured by the maintainer on one nitroxide parameter set and are not a general guarantee.

## Requirements

* g++ with C++17 (GCC 12 or newer, the Makefile uses `-ftrivial-auto-var-init=zero`)
* gfortran (the NL2SOL minimiser is compiled from the Fortran sources, `-std=legacy -ffixed-form`)
* FFTW3 (`libfftw3-dev`) and NLopt (`libnlopt-dev`): headers and libraries
* GNU make
* Eigen 3.2.5 is bundled (`Eigen/`, `unsupported/Eigen/KroneckerProduct`), no installation needed

On Debian/Ubuntu: `sudo apt install build-essential gfortran libfftw3-dev libnlopt-dev`.

## Build

    make -j4

The executable is `build/orthos`. Without root, `make TOOLCHAIN=/path/to/prefix` uses compilers and libraries from a self-contained sysroot-style prefix (`prefix/usr/bin/g++-15`, `prefix/usr/bin/gfortran-15`, `prefix/usr/lib/x86_64-linux-gnu`).

### Windows

A `Makefile.win` for MinGW-w64 (MSYS2) is provided but **has not been tested** by the maintainer of this repository. The flags in it are those of the MinGW batch builds used earlier for validation executables; `-D__linux` selects the POSIX code paths in `main.cpp` and `Discrepancy.cpp`. A Windows build needs MSYS2 with the packages `mingw-w64-x86_64-gcc`, `mingw-w64-x86_64-gcc-fortran`, `mingw-w64-x86_64-fftw`, `mingw-w64-x86_64-nlopt` and `make`; from an MSYS2 MinGW64 shell run `make -f Makefile.win`. No Visual Studio project is included.

## Run

Orthos reads its input from the current directory and writes its output there; run it in a copy of an input directory:

    mkdir run && cp examples/biradical_l616/* run/ && cd run
    ../build/orthos

Input files:

* `name.nam` - run control: number of spectra, whether to fit (`Minimization: yes/no`), number of particles (the distinct ODF/conformer entries), file names, and one line per spectrum: `experimental_file  director_angle_deg  rotation_axis  [basis_source]`.
* `Magnetic.par` - per particle: type (`Monoradical` or biradical by default), g and hyperfine tensors with Euler angles, line widths, diffusion tensor and type, potential (`Potential: yes` followed by lines `L K Value VF`), dipolar/exchange coupling, basis truncation (LE, LO, K, M), Lanczos steps, field centre, number of points, field step, integration net. Each parameter is followed by its variate flag for the fit (0 = fixed).
* `ODF.par` - orientation distribution function (order parameters) per particle.
* the experimental spectrum (`*.esr`).

Output: `spec` / `out*.dat` (simulated spectrum), `Magnetic.out` and `ODF.out` (parameters as used or fitted), `discrepancy`, and diagnostic files (`basis*.out`, `liou*.out`, `tridiag*.out`, `tensors*.out`, `potential*.out`, `timing.out`).

`examples/biradical_l616` is a small two-particle biradical example (`l616-000n.esr`, 810 points).

## Third-party code and licences

See THIRD_PARTY.md. Read it before redistributing: two components have licence terms that need review.

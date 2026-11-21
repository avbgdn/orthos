# Makefile for a native Linux build of orthos with g++ and gfortran.
#
# Usage:
#   make                              # system g++/gfortran, libraries from the default search paths
#   make -j4
#   make TOOLCHAIN=/path/to/prefix    # compilers and libraries from a self-contained prefix
#                                     # (sysroot style, e.g. unpacked with dpkg-deb -x)
#   make clean
#
# Requirements: g++ with C++17, gfortran, FFTW3 and NLopt (headers and libraries).
# Eigen 3.2.5 is bundled (Eigen/ and unsupported/Eigen/KroneckerProduct).
# The NL2SOL minimiser is compiled from NL2SOL1.for and Nl2def1.for.

CXX ?= g++
FC ?= gfortran

ifdef TOOLCHAIN
SYSROOT := --sysroot=$(TOOLCHAIN)
LIBDIR := $(TOOLCHAIN)/usr/lib/x86_64-linux-gnu
CXX := $(TOOLCHAIN)/usr/bin/g++-15
FC := $(TOOLCHAIN)/usr/bin/gfortran-15
# the compiler binaries themselves are linked against libraries in the prefix
export LD_LIBRARY_PATH := $(LIBDIR)
LDFLAGS := $(SYSROOT) -L$(LIBDIR) -Wl,-rpath,$(LIBDIR)
endif

# -ftrivial-auto-var-init=zero works around an uninitialised buffer in
# RWMagnetic.cpp (ReadFromFile): on a blank line of Magnetic.par the sscanf
# for the keyword fails and leaves stale stack content in the keyword buffer,
# which can make the per-particle field scan end early. Zero-initialising
# locals makes the read deterministic. Requires GCC 12 or newer.
CXXFLAGS := $(SYSROOT) -std=c++17 -O2 -ftrivial-auto-var-init=zero -I.
FFLAGS := $(SYSROOT) -std=legacy -ffixed-form -ffixed-line-length-none -O2
LDLIBS := -lfftw3 -lnlopt -lgfortran -lm

OUTDIR := build
OBJDIR := $(OUTDIR)/obj
TARGET := $(OUTDIR)/orthos

CPPSRCS := BasisB.cpp BesselFunctions.cpp BiradicalTransitions.cpp Discrepancy.cpp \
    GaussKronrodIntegration.cpp Hyperfine.cpp IndividualSpectrum.cpp IntegratedSpectrum.cpp \
    Lanczos.cpp LiouvilleMatrixElementB.cpp LiouvillianB.cpp Multiplet.cpp ODF.cpp \
    Potential.cpp RWMagnetic.cpp RWODF.cpp ReadExperimental.cpp ReadNamenam.cpp \
    RotationMatrix.cpp SphericalTensor.cpp SpinMatrixBiradicals.cpp Version.cpp \
    VoitlanderIntegration.cpp VoitlanderTriangleIntegrals.cpp \
    Wigner3jMy.cpp kronrod.cpp main.cpp
CCSRCS := faddeeva.cc
FSRCS := NL2SOL1.for Nl2def1.for

CPPOBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(CPPSRCS))
CCOBJS := $(patsubst %.cc,$(OBJDIR)/%.o,$(CCSRCS))
FOBJS := $(patsubst %.for,$(OBJDIR)/%.o,$(FSRCS))

.PHONY: all clean
all: $(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.cc | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.for | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(TARGET): $(CPPOBJS) $(CCOBJS) $(FOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(CPPOBJS) $(CCOBJS) $(FOBJS) $(LDLIBS)

clean:
	rm -rf $(OUTDIR)

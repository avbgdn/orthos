#ifndef BIRADICALODF_H

#define BIRADICALODF_H

// #define __MPI

#include "Eigen/Dense"
#include <fstream>
#include <iostream>
#include <complex>
#include <cmath>
#include <time.h>
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

// BiradicalODF.h: the main header of Orthos. It defines the global size limits and the default bounds of the fitted
// parameters (MAX_RANK, MAX_ORIENTATIONS, MAX_POT, MAX_TRUNC, MIN_GX, MAX_GX, ...), declares the classes shared by the whole
// program (Hyperfine, Multiplet, Potential, OrderParameters, MagneticParameters, Experimental, DataForCalculation) and
// declares many free functions, grouped under the name of the .cpp file that defines them. It contains no function
// definitions, only declarations.
// Used by: included directly by BasisB.cpp, BiradicalTransitions.cpp, Discrepancy.cpp, GaussKronrodIntegration.cpp,
// Hyperfine.cpp, IndividualSpectrum.cpp, IntegratedSpectrum.cpp, Lanczos.cpp, LiouvillianB.cpp, Multiplet.cpp, ODF.cpp,
// RWMagnetic.cpp, ReadExperimental.cpp, ReadNamenam.cpp, RotationMatrix.cpp, SphericalTensor.cpp,
// SpinMatrixBiradicals.cpp, Version.cpp and main.cpp, and through LiouvilleMatrixElementB.h and VoitlanderIntegration.h
// by the files that include those (for example Potential.cpp and Wigner3jMy.cpp).

#include "BasisB.h"
#include "SphericalTensor.h"
#include "SpinMatrixB.h"

using namespace Eigen;
using namespace std;

#define PI 3.14159265357989
// Array limits: MAX_RANK is the highest ODF rank Jmax (class OrderParameters), MAX_ORIENTATIONS the highest number of
// recorded sample orientations (spectra), and MAX_ARG the highest number of fitted parameters in the optimiser argument
// array (checked by NextArg).
#define MAX_RANK 20
#define MAX_ORIENTATIONS 50

#define MAX_ARG 50

#ifndef MAX_STRLEN
#define MAX_STRLEN 1000
#endif

// MAX_POT is the highest even rank L of the mean-field potential expansion (Potential::A); the order parameters of the
// potential are stored up to L=2*MAX_TRUNC (Potential::OP, SOP and the jump arrays); MAX_PARTICLES sizes the per-particle
// arrays such as the conformer conversion-rate matrix.
#define MAX_POT 6
#define MAX_TRUNC 60

#define MAX_PARTICLES 10

// Maximum number of axes for discrete jump model:
#define MAX_JUMPS 10

// Lower and upper bounds!

#define MIN_GX 2.0080
#define MAX_GX 2.0100
#define MIN_GY 2.0040
#define MAX_GY 2.0080
#define MIN_GZ 2.0000
#define MAX_GZ 2.0040

#define MIN_AX  3.0
#define MAX_AX  8.0
#define MIN_AY  3.0
#define MAX_AY  8.0
#define MIN_AZ 30.0
#define MAX_AZ 36.0

#define MIN_ANGLE 0.0
#define MAX_ANGLE 2*PI

#define MIN_LW 0.0
#define MAX_LW 20.0

#define MIN_DX 0.0
#define MAX_DX 30.0

#define MIN_DY 0.0
#define MAX_DY 30.0

#define MIN_DZ 0.0
#define MAX_DZ 30.0

#define MIN_POTENTIAL -10
#define MAX_POTENTIAL 100

#define MIN_CONT 0.0
#define MAX_CONT 2.0

#define MIN_SDVIG -5.0
#define MAX_SDVIG 5.0

#define MIN_EXCHANGE 0.0
#define MAX_EXCHANGE 5000

#define MIN_DIFFUSION 0.0
#define MAX_DIFFUSION 1.0e11

#define MAX_HFS 10
// minimum libration amplitude:
#define LIB_MIN 0.0001
// integration net for librations:
#define NS 30

// maximum integration net for Gauss-Kronrod procedure
#define MAX_GKNODES 2001






// Hyperfine Splitting:
// File Hyperfine.cpp
// Hyperfine: one group of equivalent nuclei with an isotropic hyperfine coupling: atom type, number of equivalent nuclei
// (Degeneracy), nuclear spin (set from the atom type by SetMagneticFromName), and the coupling constant HFConstant with its
// variate flag and optimiser argument number. MagneticParameters::HFS holds up to MAX_HFS such groups; member functions
// are in Hyperfine.cpp.
class Hyperfine
{
public:
	char AtomType[MAX_STRLEN];
	int Degeneracy;
	double NuclearSpin;
	int Sx2;
	double HFConstant, HFConstantCoefficient;
	int HFConstantNArg, HFConstantVF;
	int Count;

	void SetMagneticFromName();
	Hyperfine();
	~Hyperfine();
};

// File Multiplet.cpp
// Multiplet: the stick spectrum of the isotropic hyperfine splitting of all groups in HFS: NLines lines, each with a
// field shift (column 0 of Transitions) and a relative intensity (column 1); filled by SetFromHFS (Multiplet.cpp) and
// applied to the simulated spectra in Discrepancy.cpp.
class Multiplet
{
public:
	Array<double, Dynamic, 3> Transitions;
//	Array<double, Dynamic, 1> PrincipalSpinProjection;
	int NLines;

	Multiplet();
	~Multiplet();

	void SetFromHFS(double Field0, Hyperfine * HFS, int NNuclei);

};


void CalculateMultiplet(int Sx2, int N, int *mult);


// IndividualSpectrum.cpp
double LorentzDerivativeLineshape(double Hrel, double HG, double HL);
double   GaussDerivativeLineshape(double Hrel, double HG, double HL);
double   VoigtDerivativeLineshape(double Hrel, double HG, double HL);

void CalculateIndividualSpectrum (Array<double, 54, 2> &Transitions,
								  double HGaussEff,
								  double HLorentzEff,
								  int NPoints,
								  double *FieldNet,
								  double *SpectrumAmplitude,
								  double IntensityCutoff);
     
void CalculateIndividualSpectrumM(Array<double, 3, 2> &Transitions,
								  double HGaussEff,
								  double HLorentzEff,
								  int NPoints,
								  double *FieldNet,
								  double *SpectrumAmplitude,
								  double IntensityCutoff);

// RotationMatrix.cpp
Matrix<double, 3,3> RotationMatrix(double EulerAngles[3]);


// File Potential.cpp
// Potential: the mean-field orienting potential of a particle and the quantities derived from it (member functions in
// Potential.cpp). A holds the potential coefficients A(L,K) (even L and K, K>=0) of the expansion over Wigner functions
// D^L_0K, with their variate flags VF and optimiser argument numbers NArg. OP holds the order parameters of the
// equilibrium orientation distribution and SOP those of its square root; OP2 and SOP2 are the same integrals computed with
// the embedded Gauss rule (error estimate); OP_JUMP_P, OP_JUMP_M and the *2 copies are the order parameters of the
// discrete-jump model for each jump axis; XLK_Matrix holds the coefficients of the rotational diffusion operator in the
// potential. INT_NET and INT_NET_FI are the integration nets used by CalculateOP. Each MagneticParameters object owns one
// Potential (member MFPotential) and MagneticParameterParent points back to it.
class MagneticParameters;
class Potential
{
public:
	int MyCount;		//Particle count
	double A[(MAX_POT+2)*(MAX_POT+4)/8];
	double SOP[(2*MAX_TRUNC+2)*(2*MAX_TRUNC+4)/8];   // sqrt of order parameters 
	double OP[(2*MAX_TRUNC+2)*(2*MAX_TRUNC+4)/8];    // order parameters
	double SOP2[(2*MAX_TRUNC+2)*(2*MAX_TRUNC+4)/8];   // sqrt of order parameters - Gauss integration
	double OP2[(2*MAX_TRUNC+2)*(2*MAX_TRUNC+4)/8];    // order parameters - Gauss integration
	std::complex<double> OP_JUMP_P[MAX_JUMPS][MAX_TRUNC*MAX_TRUNC];
	std::complex<double> OP_JUMP_M[MAX_JUMPS][MAX_TRUNC*MAX_TRUNC];

	std::complex<double> OP_JUMP_P2[MAX_JUMPS][MAX_TRUNC*MAX_TRUNC];  // Gauss Integration
	std::complex<double> OP_JUMP_M2[MAX_JUMPS][MAX_TRUNC*MAX_TRUNC];  // Gauss Integration

	bool isAxial;
	int VF[(MAX_POT+2)*(MAX_POT+4)/8];
	int NArg[(MAX_POT+2)*(MAX_POT+4)/8];
	std::complex<double> XLK_Matrix[2*MAX_POT+1][4*MAX_POT+1][4*MAX_POT+1];
	int Lmax;
	double Normalization_OP;
	MagneticParameters* MagneticParameterParent;
	double Get(int L, int K);
	int GetVF(int L, int K);
	int GetNArg(int L, int K);
	double GetOP(int L, int K);
	double GetSOP(int L, int K, int M);
	std::complex<double> GetOP_JUMP_P(int L, int K, int NAxis);
	std::complex<double> GetOP_JUMP_M(int L, int K, int NAxis);

	void CalculateOP(int MaxL);
	void CalculateOP_GK(int MaxL);

	void CalculateOP_Jumps_GK(int MaxL);

	void CalculateXLK(int MaxL);
	std::complex<double> GetXLK(int L, int K, int M);
	void Add(int L, int K, double Value, int VariateFlag);
	void SetNArg(int L, int K, int NA);
	void PrintOut(int MaxL);
	int INT_NET, INT_NET_FI;
	Potential();
	~Potential();
};


//////////////////////////////////////////////////////////////////////////
////////  CLASS FOR ODF AND ORDER PARAMETER MANIPULATION /////////////////
//////////////////////////////////////////////////////////////////////////

double wigner_d(int j, int m, int n, double beta);
std::complex<double> ComplexWigner(int j, int m, int n, double alpha, double beta, double gamma);

// OrderParameters: the orientation distribution function (ODF) of the particles, the sample orientations of the recorded
// spectra and the fitted angle errors; member functions are in ODF.cpp, where the members are described in more detail.
class OrderParameters
{

	std::complex<double> Ajmn_container[MAX_RANK+1][2*MAX_RANK+1][2*MAX_RANK+1];
	std::complex<double> Tjpqi_container[MAX_RANK+1][2*MAX_RANK+1][2*MAX_RANK+1][MAX_ORIENTATIONS];
public:
// self:
	OrderParameters *ODFArray;
	int Count;
	int Jmax, NOrientations;
	double alpha_SF[MAX_ORIENTATIONS];		// for each of the orientations
	double beta_SF[MAX_ORIENTATIONS];
	double gamma_SF[MAX_ORIENTATIONS];
	char BasisSource[MAX_ORIENTATIONS][30];

	double hi_SF[MAX_ORIENTATIONS];		// for each of the orientations
	int axis_SF[MAX_ORIENTATIONS];      // Axis for rotation: 0 - x
										//                    1 - y
										//                    2 - z

	double DelHi0, DelHi1, DelHi2;      // Angle error
	double CurrentODF[MAX_ORIENTATIONS];
	double Penalty, NegativeVolume, NegativeVolumeODF3;
	char ExperimentalFiles[MAX_ORIENTATIONS][MAX_STRLEN];

// Variate Flags:
	int Ajmn_containerVF[MAX_RANK+1][2*MAX_RANK+1][2*MAX_RANK+1];
    int DelHi0VF, DelHi1VF, DelHi2VF;

// Number in Arguments Array:
	int Ajmn_containerNArg[MAX_RANK+1][2*MAX_RANK+1][2*MAX_RANK+1];
	int DelHi0NArg, DelHi1NArg, DelHi2NArg;
//////////////////////////////////////////////////////////////////////////
	OrderParameters();
//////////////////////////////////////////////////////////////////////////
	~OrderParameters();

//////////////////////////////////////////////////////////////////////////
	void AssignTjpqi();
//////////////////////////////////////////////////////////////////////////
	std::complex<double> A_jmn(int j, int m, int n);
	std::complex<double> T_jpqi(int j, int p, int q, int i);
	int A_jmnVF(int j, int m, int n);
	int A_jmnNArg(int j, int m, int n);
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
	void ODF_reduced(double teta, double ksi);
	double ODF_value(double ksi, double teta, double psi);
// Given the orientation distribution function A_jmn,
// and the set of angles _SF for all sample orientations,
// this subroutine calculates the number density of
// particles which are oriented in resonator by angles ksi, teta.
// The calculated values are stored in array CurrentODF[0..NOrientations],
// to be used in the block of spectra calculation.
//////////////////////////////////////////////////////////////////////////
	void Set_A_jmn(int j, int m, int n, std::complex<double> value);
	void Set_A_jmnVF(int j, int m, int n, int value);
	void Set_A_jmnNArg(int j, int m, int n, int value);

// Sets the value of specified order parameter and variate flag

//////////////////////////////////////////////////////////////////////////
	void Realize();
// This function makes the ODF A_jmn real by formula (3.5-6) in dissertation
	void ReadFromFile(char* filename);
	void WriteToFile(char* filename);
	void WriteToFileBrief(char* filename);
	void SetEuler_SF(void);
	void ConvertAnglesToDegrees(void);
	void ConvertAnglesToRadians(void);

	void Par2Arg(double* Arguments, double *ub, double *lb, int *NArg);
	void Arg2Par(const double* Arguments);

};


// RWMagnetic.cpp
// MagneticParameters: all spin Hamiltonian, relaxation, motion and numerical parameters of one particle (conformer): g
// and hyperfine tensors, dipolar and exchange coupling, line widths, rotational diffusion and the mean-field potential
// MFPotential, discrete-jump and libration models, the
// integration nets and the basis truncations (LOmax, LEmax, Mmax, Kmax, Nmax). Every fitted parameter has a variate flag
// (VF) and an optimiser argument number (NArg). Member functions are in RWMagnetic.cpp.
class MagneticParameters
{
public:
// self:
	MagneticParameters *MagneticParameterArray;


// Whether we should calculate spectrum via SLE procedure:
	bool Liouville;
	bool isMonoradical;
// Minimizator parameters:
	char NLOptAlgorithm[MAX_STRLEN];
	double XTol;
	double FTol;
	double InitialStep;
// Coefficient for bound optimization
	double AmplitudeCoefficient;
// Initializing Magnetic parameters
// Particle number
	int Count;

// Conformation analysis:
	int ConformerGroup, NConformers;
	double ConversionRate[MAX_PARTICLES][MAX_PARTICLES];
	double *ConversionRateR;
	int ConformerGroupMembers[MAX_PARTICLES];

	Matrix<double, 3, 3> g1, g2, A1, A2, HGauss, HLorentz;
	Matrix<double, 3, 3> DJ;
    double dx, dy, dz, J0, J0sign, J2, JDiffusion, JDiffusionR;
	double ZFS_E, ZFS_Angle;
	double Rx, Ry, Rz;
	double g1Average3, g2Average3, A1Average3, A2Average3;
// Mean field potential:
	Potential MFPotential;
// DiscreteJumps
//	int JumpNSites;
//	double JumpFrequency, JumpFrequencyR;
// Discrete Jumps by my model:
	int NJumpAxes;
	double JumpAxisFi[MAX_JUMPS], JumpAxisTeta[MAX_JUMPS], JumpAngle[MAX_JUMPS], JumpFrequency[MAX_JUMPS], JumpFrequencyR[MAX_JUMPS];
// Jump diffusion and free diffusion possibility:
	double DiffusionType[3];
	double RTau[3];
// Quasi-librations
	double LibrationAmplitude[3];
	double LibrationEuler[3];
	double LibDistrWidth[3];
	int LibDistrNSites;
// AfterRotation:
	Matrix<double, 3, 3> g1R, g2R, A1R, A2R, HGaussR, HLorentzR; 
	Matrix<std::complex<double>, 3, 3> DiffusionR;

// Elements for calculation of Liouvillian:
	SphericalTensor *Fmg1, *Fmg2, *FmA1, *FmA2, *FmD, *FmW, *FmR;
	SpinInteractionB Interactions[6];


	double CenterField0, Shift, CenterField;
	double g1Euler[3], g2Euler[3], A1Euler[3], A2Euler[3], HGaussEuler[3], HLorentzEuler[3], DiffusionEuler[3];

//	double HGaussEff, HLorentzEff;
// Initializing magnetic field net related parameters
	double *FieldNet, *Spectrum, fs;
	int NPoints;

	double AmplitudeFactor, AmplitudeFactorAllConformers;
	double TransitionIntensityCutoff;
// Integration net for orientations
	int IntegrationNet, IntegrationNetKsi;

// orientation-integration method selector: IntegrateGauss (fixed Gauss-Legendre
// grid only), IntegrateGaussKronrod (the same grid plus Kronrod-only nodes, for
// a convergence diff) or IntegrateVoitlander (adaptive hemisphere triangulation,
// VoitlanderIntegration.h); replaces the former bool IntegrateGaussKronrod
#define INTEGRATE_GAUSS 0
#define INTEGRATE_GAUSSKRONROD 1
#define INTEGRATE_VOITLANDER 2

	short IntegrationMethod;
	bool AmplitudeNormalized;
// Variate Flags:
	Matrix<int, 3, 3> g1VF, g2VF, A1VF, A2VF, HGaussVF, HLorentzVF;
	int dxVF, dyVF, dzVF, J0VF, J2VF, JDiffusionVF;
	int ZFS_EVF, ZFS_AngleVF;
	int RxVF, RyVF, RzVF;
	int ShiftVF;
	int g1EulerVF[3], g2EulerVF[3], A1EulerVF[3], A2EulerVF[3], HGaussEulerVF[3], HLorentzEulerVF[3], DiffusionEulerVF[3];
	int AmplitudeFactorVF;
//	int JumpFrequencyVF;
	int JumpAxisFiVF[MAX_JUMPS], JumpAxisTetaVF[MAX_JUMPS], JumpAngleVF[MAX_JUMPS], JumpFrequencyVF[MAX_JUMPS];  // Discrete jumps by my model
// Jump diffusion:
	int RTauVF[3];
// Conformer analysis:
	int ConversionRateVF[MAX_PARTICLES][MAX_PARTICLES];

	int LibrationAmplitudeVF[3], LibrationEulerVF[3], LibDistrWidthVF[3];
// Numbers in Arguments Array:
	Matrix<int, 3, 3> g1NArg, g2NArg, A1NArg, A2NArg, HGaussNArg, HLorentzNArg;
	int dxNArg, dyNArg, dzNArg, J0NArg, J2NArg, JDiffusionNArg;
	int ZFS_ENArg, ZFS_AngleNArg;
	int RxNArg, RyNArg, RzNArg;
	int ShiftNArg;
	int g1EulerNArg[3], g2EulerNArg[3], A1EulerNArg[3], A2EulerNArg[3], HGaussEulerNArg[3], HLorentzEulerNArg[3], DiffusionEulerNArg[3];
	int AmplitudeFactorNArg;
//	int JumpFrequencyNArg;
// Jump diffusion:
	int RTauNArg[3];
	int LibrationAmplitudeNArg[3], LibrationEulerNArg[3], LibDistrWidthNArg[3];
	int JumpAxisFiNArg[MAX_JUMPS], JumpAxisTetaNArg[MAX_JUMPS], JumpAngleNArg[MAX_JUMPS], JumpFrequencyNArg[MAX_JUMPS];  // Discrete jumps by my model

// Conformer analysis:
	int ConversionRateNArg[MAX_PARTICLES][MAX_PARTICLES];


// Truncation values for basis
	int LOmax, LEmax, Mmax, Kmax, Nmax;
// Pruning parameters
	int NPrune;
	double dHPrune, PruningTolerance;
// Maximum number of steps for Lanczos routine
	int LanczosSteps;

// Isotropic Hyperfine
	Hyperfine HFS[MAX_HFS];
	int HyperfineCount;

	Multiplet Mult;

	double CurrentDiscrepancyOverall;
	double CurrentDiscrepancyPerPoint;


// Gauss-Kronrod integration routine
	int GaussNodes, KronrodNodes;
	double Node[MAX_GKNODES], KronWeight[MAX_GKNODES], GausWeight[MAX_GKNODES];

	int GaussNodesKsi, KronrodNodesKsi;
	double NodeKsi[MAX_GKNODES], KronWeightKsi[MAX_GKNODES], GausWeightKsi[MAX_GKNODES];

//
	void ReadFromFile(char* filename);
	void  WriteToFile(char* filename);
	void  WriteToFileBrief(char* filename);
	void RotateTensors();
	void CalculateTranstionRateMatrix();
	void GenerateFieldNet();
	MagneticParameters();
	~MagneticParameters();
	void ConvertAnglesToDegrees(void);
	void ConvertAnglesToRadians(void);
	void QuasiLibrations(void);

	void Par2Arg(double* Arguments, double *ub, double *lb, int *NArg);
	void Arg2Par(const double* Arguments);

	void GenerateIntegrationNet();
	void GenerateIntegrationNetKsi();


};








// BiradicalTransitions.cpp
Array<double, 54, 2> CalculateTransitions(MagneticParameters *MagneticParameter,
										  Matrix<double, 3, 3> Field_to_Molecule_RotMat,
										  double *HGaussEff, double *HLorentzEff
										  );

Array<double, 3, 2> CalculateTransitionsM(MagneticParameters *MagneticParameter,
										  Matrix<double, 3, 3> Field_to_Molecule_RotMat,
										  double *HGaussEff, double *HLorentzEff
										  );

// file IntegratedSpectrum.cpp
void CalculateIntegratedSpectrumLiouville(MagneticParameters *MagneticParameter,
										  OrderParameters *ODF,
										  double **SpectrumSum);
// ReadNamenam.cpp
int ReadName(MagneticParameters *MagneticParameter, OrderParameters *ODF, char* Magnetic_in_filename, char *ODF_in_filename, char* Magnetic_out_filename, char *ODF_out_filename);

// ReadExperimental.cpp
// Experimental: one experimental spectrum read from file (SpectrumX, SpectrumY), its interpolation on the simulation field
// net (SpectrumInterpolated) and noise statistics; member functions are in ReadExperimental.cpp.
class Experimental
{
public:
	double *SpectrumX;
	double *SpectrumY;
	double *SpectrumInterpolated;
	int NPoints; 

// noise statistics
	double MeanL, VarianceL;

	void CalculateStats(double Hmax);

	void LoadFromFile(char * filename);
	double Interpolate(double X);
	void InterpolateWhole(double *FieldNet, int Npts);

	Experimental();
	~Experimental();


};


void WriteOutput(int NParticles, 
				 Experimental *ExpSpectra, 
				 double ***SpectrumSum, 
				 MagneticParameters *MagneticParameter, 
				 OrderParameters *ODF,
				 int NCovarianceArguments,
				 const double * Arguments);

// File RWMagnetic.cpp

// Advances the index into the flat optimiser argument arrays, with a bound check against MAX_ARG
void NextArg(int *NArg);


// File Discrepancy.cpp

// for multiprocessor usage: in integration over N orientations by P processes,
// Processor MpiRank integrates over orientations k,  StartMPI() <= k < EndMPI.
int StartMPI(int Rank, int N, int P);
int EndMPI(int Rank, int N, int P);


double Discrepancy(int NParticles,
				   Experimental *ExpSpectra, 
				   MagneticParameters *MagneticParameter, 
				   OrderParameters *ODF, const double *Arguments, double *R);


// DataForCalculation: the problem data (particles, experimental spectra, magnetic parameters and ODFs) that the optimiser
// hands back to the residual functions madr and myfunc (Discrepancy.cpp) through their user-data pointer.
struct DataForCalculation
{
	int NParticles;
	Experimental *ExpSpectra;
	MagneticParameters *MagneticParameter;
	OrderParameters *ODF;
};

// File LiouvillianB.cpp

Matrix<double, Dynamic, 1> CalculateLiouvillianB(MagneticParameters *MagneticParameter, double psi, char* BasisSource);
void Calculate_Wigner3jm(int Lmax);

// File Version.cpp

void Version(void);


#endif

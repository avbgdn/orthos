// Isotropic hyperfine coupling to nuclei. An object of class Hyperfine (declared in BiradicalODF.h) describes one
// group of Degeneracy equivalent nuclei of one type (AtomType) with the isotropic hyperfine constant HFConstant.
// This file only provides the constructor, the destructor and the lookup of the nuclear spin from the nucleus name.
// Used by: RWMagnetic.cpp (ReadFromFile fills the array HFS of MagneticParameters, Par2Arg calls
// SetMagneticFromName, Arg2Par updates HFConstant), Multiplet.cpp (Multiplet::SetFromHFS turns the list of nuclei
// into hyperfine line positions and intensities), Discrepancy.cpp (applies that multiplet to the spectra)
// and main.cpp (ties the constants of different particles together).

#include "BiradicalODF.h"


// Constructor: does nothing, the members are filled by MagneticParameters::ReadFromFile and Par2Arg (RWMagnetic.cpp).
// Used by: implicitly, when each MagneticParameters object contains its array HFS[MAX_HFS] (objects created
// by new MagneticParameters in main.cpp).
Hyperfine::Hyperfine()
{}

// Destructor: does nothing.
// Used by: implicitly, when a MagneticParameters object holding the array HFS is destroyed (never done explicitly).
Hyperfine::~Hyperfine()
{}

// Sets the nuclear spin NuclearSpin from the nucleus name AtomType: H, F, P, C13, N15 have spin 1/2; D, H2, N have
// spin 1; Cu, Cu63, Cu65 have spin 3/2. Any other name is read as a number giving the spin itself; a name that is
// neither known nor a non-zero number prints a warning.
// Input: AtomType of this object. Output: NuclearSpin of this object.
// Used by: MagneticParameters::Par2Arg (RWMagnetic.cpp), for every nucleus of the array HFS.
void Hyperfine::SetMagneticFromName()
{
	// Table of known nuclei: set the spin and return at the first name that matches
	if(!strcmp(AtomType, "H"))
	{
		NuclearSpin=0.5;
		return;
	}
	if(!strcmp(AtomType, "D"))
	{
		NuclearSpin=1;
		return;
	}
	if(!strcmp(AtomType, "H2"))
	{
		NuclearSpin=1;
		return;
	}
	if(!strcmp(AtomType, "C13"))
	{
		NuclearSpin=0.5;
		return;
	}
	if(!strcmp(AtomType, "N"))
	{
		NuclearSpin=1.0;
		return;
	}
	if(!strcmp(AtomType, "N15"))
	{
		NuclearSpin=0.5;
		return;
	}
	if(!strcmp(AtomType, "F"))
	{
		NuclearSpin=0.5;
		return;
	}
	if(!strcmp(AtomType, "P"))
	{
		NuclearSpin=0.5;
		return;
	}
	if(!strcmp(AtomType, "Cu"))
	{
		NuclearSpin=1.5;
		return;
	}
	if(!strcmp(AtomType, "Cu63"))
	{
		NuclearSpin=1.5;
		return;
	}
	if(!strcmp(AtomType, "Cu65"))
	{
		NuclearSpin=1.5;
		return;
	}
	
	// Unknown name: read AtomType as a number giving the nuclear spin, and warn if that gives zero
	if(!(NuclearSpin=atof(AtomType)))
		printf("Unrecognized nucleus No. %d : %s", Count, AtomType);

}

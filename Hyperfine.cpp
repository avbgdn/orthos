#include "BiradicalODF.h"


Hyperfine::Hyperfine()
{}

Hyperfine::~Hyperfine()
{}

void Hyperfine::SetMagneticFromName()
{
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
	
	if(!(NuclearSpin=atof(AtomType)))
		printf("Unrecognized nucleus No. %d : %s", Count, AtomType);

}

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

// Rotation matrix of the z-y-z Euler angles, used throughout the program to turn interaction tensors
// and field directions from one axis system into another.
// Used by: RWMagnetic.cpp (RotateTensors, QuasiLibrations), Discrepancy.cpp (Discrepancy),
// IntegratedSpectrum.cpp (CalculateIntegratedSpectrumLiouville) and VoitlanderTriangleIntegrals.cpp
// (VoitlanderVertexEvalRigid); BiradicalTransitions.cpp contains only commented-out calls.

#include "BiradicalODF.h"
#include <math.h>

// Builds the 3x3 rotation matrix R=R3*R2*R1 of the Euler angles EulerAngles[0..2]=(fi, teta, ksi), in radians.
// R1 and R3 are rotations about z by fi and ksi with matrix [cos sin 0; -sin cos 0; 0 0 1], R2 is a rotation
// about y by teta with matrix [cos 0 -sin; 0 1 0; sin 0 cos]. RotateTensors and QuasiLibrations rotate a tensor T
// as R*T*transpose(R).
// Input: EulerAngles, three angles in radians. Output: the 3x3 matrix R.
// Used by: RotateTensors and QuasiLibrations (RWMagnetic.cpp), Discrepancy (Discrepancy.cpp, rotation from
// the field frame to the molecular frame for the rigid-limit calculation), CalculateIntegratedSpectrumLiouville
// (IntegratedSpectrum.cpp, same purpose) and VoitlanderVertexEvalRigid (VoitlanderTriangleIntegrals.cpp).
Matrix<double, 3,3> RotationMatrix(double EulerAngles[3])
{
 // Unpack the three Euler angles
 Matrix<double, 3,3> rotationMatrix, R1, R2, R3;
 double fi, teta, ksi;
 fi=EulerAngles[0];
 teta=EulerAngles[1];
 ksi=EulerAngles[2];
 
 // Elementary rotations: R1 about z by fi, R2 about y by teta, R3 about z by ksi
 R1<< cos(fi), sin(fi), 0,
	 -sin(fi), cos(fi), 0,
	 0, 0, 1;
 R2<< cos(teta), 0, -sin(teta),
	 0, 1, 0,
	 sin(teta), 0, cos(teta);
 R3<< cos(ksi), sin(ksi), 0,
	 -sin(ksi), cos(ksi), 0,
	 0, 0, 1;

 // Combined rotation
 return R3*R2*R1;

}
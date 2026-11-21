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

#include "BiradicalODF.h"
#include <math.h>

Matrix<double, 3,3> RotationMatrix(double EulerAngles[3])
{
 Matrix<double, 3,3> rotationMatrix, R1, R2, R3;
 double fi, teta, ksi;
 fi=EulerAngles[0];
 teta=EulerAngles[1];
 ksi=EulerAngles[2];
 
 R1<< cos(fi), sin(fi), 0,
	 -sin(fi), cos(fi), 0,
	 0, 0, 1;
 R2<< cos(teta), 0, -sin(teta),
	 0, 1, 0,
	 sin(teta), 0, cos(teta);
 R3<< cos(ksi), sin(ksi), 0,
	 -sin(ksi), cos(ksi), 0,
	 0, 0, 1;

 return R3*R2*R1;

}
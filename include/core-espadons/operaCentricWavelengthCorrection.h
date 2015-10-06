#ifndef OPERACENTRICWAVELENGTHCORRECTION_H
#define OPERACENTRICWAVELENGTHCORRECTION_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaCentricWavelengthCorrection
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: May/2015
 
 Copyright (C) 2011  Opera Pipeline team, Canada France Hawaii Telescope
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see:
 http://software.cfht.hawaii.edu/licenses
 -or-
 http://www.gnu.org/licenses/gpl-3.0.html
 ********************************************************************/

#include <string>

/*! 
 * operaSpectrumCalibration
 * \author Eder Martioli / Christopher Usher
 * \brief Calculate and apply Heliocentric or Barycentric velocity wavelength correction.
 * \file operaCentricWavelengthCorrection.h
 * \ingroup libraries
 */

const double JDTime1 = 2400000.5;

int CentricWavelengthCorrection(int argc, char *argv[], const std::string modulename, const bool barycentric);

double CalculateSOFABarycentricCorrection(double latitudeInDeg, double longitudeInDeg, double elevation, double MJDTime, double raInHrs, double decInDeg, double haInHrs, bool debug, std::string modulename);

#endif

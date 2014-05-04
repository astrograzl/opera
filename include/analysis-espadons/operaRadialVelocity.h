#ifndef OPERARADIALVELOCITY_H
#define OPERARADIALVELOCITY_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaRadialVelocity
 Version: 1.0
 Description: Caluclate Radial Velocity
 to start up with an OPERA module.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: SEP/2012
 Contact: teeple@cfht.hawaii.edu
 
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

/*! \brief Apply wavelength correction based on telluric lines. */
/*! \file operaRadialVelocity.h */
/*! \ingroup core */

#define MAXNUMBEROFPOINTSINTELLURICSPECTRUM 2400000
#define MAXNUMBEROFLINESINTELLURICDATABASE 100000
#define MAXLENGTHOFLINEINTELLURICDATABASE 160
#define TYPICAL_PRESSURE_AT_MAUNAKEA 61000  // Pa
#define TYPICAL_TEMPERATURE_AT_MAUNAKEA 273 // K
#define k_BOLTZMANN_CONSTANT 1.3806503e23 // m2 kg s-2 K-1

#define TYPICAL_ATMOSPHERE_PATH_LENGTH  843500 // cm


/* prototypes */
static void printUsageSyntax(char *prgname);

unsigned readTelluricLines(string telluric_database_file, int *telluricMoleculeNumber, double *telluricLinesWavelength, double *telluricLinesIntensity);

unsigned getTelluricLinesRange(double wl0, double wlf, double **wl, double **intensity);

void generateSyntheticTelluricSpectrumUsingGaussianProfile(unsigned np, double *wavelengthVector, double *ouputSpectrum, double resolution);

unsigned readTelluricSpectrum(string telluric_spectrum, double *telluricSpectrumWavelength, double *telluricSpectrumIntensity);

unsigned getTelluricSpectrumRange(double wl0, double wlf, double **wl, double **intensity);

unsigned matchTelluricReferencewithObjectLines(double acceptableMismatch,double lineSigma, operaSpectralLines *telluricLines, operaSpectralLines *objectLines, Polynomial *wlcorrection);

void calculateWavelengthShiftByXCorr(operaSpectralElements *compSpectrum, double DWavelengthRange, double DWavelengthStep, double threshold, double *maxDWavelength, double *maxcorr);

#endif

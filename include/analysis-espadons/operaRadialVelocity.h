#ifndef OPERARADIALVELOCITY_H
#define OPERARADIALVELOCITY_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaRadialVelocity
 Version: 1.0
 Description: Measure radial velocity of source spectrum
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope
 Location: Hawaii USA
 Date: Jan/2015
 Contact: opera@cfht.hawaii.edu
 
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

/*! \brief This module measures the radial velocity of the source spectrum. */
/*! \file operaRadialVelocity.h */
/*! \ingroup core */

#define MAXNUMBEROFPOINTSINSTELLARSPECTRUM 2400000
#define MAXNUMBEROFLINESINTELLURICDATABASE 100000
#define MAXLENGTHOFLINEINTELLURICDATABASE 160
#define TYPICAL_PRESSURE_AT_MAUNAKEA 61000  // Pa
#define TYPICAL_TEMPERATURE_AT_MAUNAKEA 273 // K
#define k_BOLTZMANN_CONSTANT 1.3806503e23 // m2 kg s-2 K-1

#define TYPICAL_ATMOSPHERE_PATH_LENGTH  843500 // cm

unsigned readTelluricLines(string telluric_database_file, int *telluricMoleculeNumber, double *telluricLinesWavelength, double *telluricLinesIntensity);

unsigned getTelluricLinesRange(double wl0, double wlf, double **wl, double **intensity);

void generateSyntheticTelluricSpectrumUsingGaussianProfile(unsigned np, double *wavelengthVector, double *ouputSpectrum, double resolution);

unsigned readStellarSpectrum(string inputStellarSpectrum, double *stellarSpectrumWavelength, double *stellarSpectrumIntensity);

unsigned getStellarSpectrumRange(double wl0, double wlf, double **wl, double **intensity);

unsigned matchTelluricReferencewithObjectLines(double acceptableMismatch,double lineSigma, operaSpectralLines *telluricLines, operaSpectralLines *objectLines, Polynomial *wlcorrection, unsigned order, double wl_central, ofstream *flinesdata);

bool calculateSourceRVShiftByXCorr(unsigned nelem, double *wavelength, double *objectSpectrum, double radialVelocityRange, double radialVelocityStep, double threshold, double *maxRV, double *sigRV, double *maxcorr, ostream *fxcorrdata, ostream *fxcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double *chisqr, double centralRV);

void GenerateTelluricXCorrelationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string cleanDataFileName);

void GenerateTelluricSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename);

#endif

#ifndef OPERAWAVELENGTHCALIBRATION_H
#define OPERAWAVELENGTHCALIBRATION_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operWavelengthCalibration
 Version: 1.0
 Description: Perform wavelength calibration 
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
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

/*! \brief Perform wavelength calibration. */
/*! \file operWavelengthCalibration.h */
/*! \ingroup core */

#define MAXREFWAVELENGTHS 21000
#define MAXFULLREFWAVELENGTHS 220000
#define MAXORDEROFPOLYNOMIAL 6

class operaSpectralElements;
class operaSpectralOrder;

unsigned getRawLines(operaSpectralOrder *spectralOrder, double *rawlinecenter, double *rawlinecenterError, double *rawlineflux, double *rawlinesigma, double *rawlinewidth, double *rawlinewidth_err);

unsigned getRawLinesFromUncalSpectrum(operaSpectralOrder *spectralOrder, double *rawlinecenter, double *rawlinecenterError, double *rawlineflux, double *rawlinesigma, double *rawlinewidth, double *rawlinewidth_err, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip);

unsigned getAtlasLinesFromSpectrum(operaWavelength *wavelength, double rawlinewidth, double uncalibrated_linewidth, int order, double *atlasLineswl,double *atlasLineswlError,double *atlasLinesflux, ofstream& fatlasdata, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, bool generate3DPlot);

unsigned readThoriumArgonAtlas(string thorium_argon_atlas, double *thAtlasWavelength, double *thAtlasIntensity);

unsigned getThoriumArgonAtlasRange(double wl0, double wlf, double **thwl, double **thi);

unsigned readAtlasSpectrum(string atlas_spectrum, double *atlasWavelength, double *atlasIntensity, double *atlasVariance);

unsigned getAtlasSpectrumRange(double wl0, double wlf, double **thwl, double **thi, double **thvar);

//unsigned createSimulatedSpectrum(unsigned nLines, double *lineCenters, double *lineAmplitudes, double lineSigma, double minwl, double maxwl, double wlstep, double *outputwl, double *outputSpectrum);

void GenerateWavelengthOrdersPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display);

void GenerateWavelengthSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, string linesdatafilename, bool subtractCentralWavelength, bool display);

void GenerateWavelength3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, bool subtractCentralWavelength, bool display);

void GenerateWavelengthSolutionPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display);

void calculateInitialSolutionFromLineSet(string inputLineSetFilename, operaSpectralOrder *spectralOrder, int order, unsigned maxorderofpolynomial);

unsigned readLineSet(string inputLineSetFilename, int order, double *wavelengthData, double *wavelengthErrors, double *distanceData);

#endif

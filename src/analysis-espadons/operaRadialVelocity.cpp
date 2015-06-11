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


#include <stdio.h>
#include <stdarg.h>
#include <getopt.h>
#include <fstream>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaSpectralLines.h"
#include "libraries/Gaussian.h"
#include "libraries/Polynomial.h"						// for Polynomial

#include "analysis-espadons/operaRadialVelocity.h"
#include "core-espadons/operaWavelengthCalibration.h"

#include "libraries/operaSpectralTools.h"

#include "libraries/operaLibCommon.h"					// for doubleValue_t
#include "libraries/operaLib.h"							// for itos
#include "libraries/operaMath.h"						// for LengthofPolynomial
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFFT.h"							// for operaXCorrelation
#include "libraries/operaFit.h"							// for operaFitSplineDouble
#include "libraries/gzstream.h"							// for gzstream - read compressed reference spectra

#define NOTPROVIDED -999

/*! \file operaRadialVelocity.cpp */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*
 * the reference Telluric spectrum
 */
static unsigned nPointsInStellarSpectrum = 0;
static double stellarSpectrumWavelength[MAXNUMBEROFPOINTSINSTELLARSPECTRUM];
static double stellarSpectrumIntensity[MAXNUMBEROFPOINTSINSTELLARSPECTRUM];

/*
 * the reference Telluric spectral lines
 */
static unsigned ntelluriclines = 0;
static int telluricMoleculeNumber[MAXNUMBEROFLINESINTELLURICDATABASE];
static double telluricLinesWavelength[MAXNUMBEROFLINESINTELLURICDATABASE];
static double telluricLinesIntensity[MAXNUMBEROFLINESINTELLURICDATABASE];

/* prototypes */
static void printUsageSyntax(char *prgname);

/*!
 * operaRadialVelocity
 * \author Eder Martioli
 * \brief This module measures the radial velocity of the source spectrum.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \note --wave=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
    int opt;
    
    string inputWaveFile;
    string inputObjectSpectrum;
    string inputFlatFluxCalibration;
    
    string outputRVFile;
    string telluric_lines;      // HITRAN Library
    string inputStellarSpectrum;   // A telluric spectrum
    
    string inputBarycentricCorrection;
    
    string objectname;          // OBJECT
    
    int expnum=0;                 // EXPNUM
    double mjd=0;                 // MJDATE
    double exptime=0;             // EXPTIME
    double airmass=0;             // AIRMASS
    double outsideTemperature=0;  // TEMPERAT
    double windspeed=0;           // WINDSPED
    double relativeHumididy=0;    // RELHUMID
    double atmoPressure=0;        // PRESSURE
    
    /*
     * Parameters for telluric wavelength correction
     */
    
    int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;
    bool maxorderprovided = false;
    
    bool StarPlusSky = false;
    /*
     * The parameters below we don't know yet whether they would be useful if used as input
     */
    double spectralResolution = 80000; // Input spectral resolution as reference for line detection.
    
    double radialVelocitySearchRange = 200;  // km/s
    double radialVelocitySearchStep = 0.5;  // km/s
    
    double XCorrelationThreshold = 0.05;    //
    unsigned normalizationBinsize = 110;
    
    string inputWavelengthMask;
    
    bool useFitToFindMaximum = TRUE;
    
    bool plot=false;
    
    string xcorrsplotfilename;
    string specplotfilename;
    
    string xcorrscriptfilename;
    string specscriptfilename;
    
    string specdatafilename;
    
    string xcorrdatafilename;
    string xcorrfitdatafilename;
    
    struct option longopts[] = {
        {"inputWaveFile",1, NULL, 'w'},				// input wavelength calibration file (.wcal)
        {"inputObjectSpectrum",1, NULL, 'i'},		// input object spectrum file (.e or .p)
        {"inputFlatFluxCalibration",1, NULL, 'm'},  // flat field spectrum ff_
        {"outputRVFile",1, NULL, 'o'},			// output wavelength calibration file (.auto or .pauto)
        {"telluric_lines",1, NULL, 'L'},            // atlas of telluric lines (HITRAN)
        {"inputStellarSpectrum",1, NULL, 'T'},      // stellar spectrum
        {"inputBarycentricCorrection",1, NULL, 'B'},// barycentric correction
        {"spectralResolution",1, NULL, 'R'},        // Input spectral resolution (wl/dwl) as reference for line detection
        {"radialVelocitySearchRange",1, NULL, 'r'},       // radial Velocity Range (in km/s) to scan for first order correction
        {"radialVelocitySearchStep",1, NULL, 's'},        // radial Velocity Step step (in km/s) to scan for first orer correction
        {"XCorrelationThreshold",1, NULL, 'x'},     // X-correlation lower threshold to consider a match between telluric and object spectra
        {"normalizationBinsize",1, NULL, 'b'},      // normalization binsize
        {"StarPlusSky",1, NULL, 'k'},               // starPlusSky mode
        {"useFitToFindMaximum",1, NULL, 'f'},       // use fit to find maximum

        {"headerData",1, NULL, 'H'},       // a set of header data
        
        {"ordernumber",			1, NULL, 'O'},
        {"minorder",			1, NULL, 'M'},
        {"maxorder",			1, NULL, 'X'},
        {"inputWavelengthMask",			1, NULL, 'A'}, // Telluric wavelength mask
        
        {"xcorrsplotfilename",1, NULL, 'P'},
        {"specplotfilename",1, NULL, 'Q'},
        {"xcorrscriptfilename",1, NULL, 'S'},
        {"specscriptfilename",1, NULL, 'U'},
        
        {"xcorrdatafilename",1, NULL, 'D'},
        {"xcorrfitdatafilename",1, NULL, 'F'},
        {"specdatafilename",1, NULL, 'G'},
        
        {"plot",0, NULL, 'p'},
        {"verbose",0, NULL, 'v'},
        {"debug",0, NULL, 'd'},
        {"trace",0, NULL, 't'},
        {"help",0, NULL, 'h'},
        {0,0,0,0}};
    
    while((opt = getopt_long(argc, argv, "w:i:m:o:L:T:B:R:r:s:x:b:k:f:H:O:M:X:A:P:Q:S:U:D:F:G:p::v::d::t::h",
                             longopts, NULL))  != -1)
    {
        switch(opt)
        {
            case 'w':
                inputWaveFile = optarg;
                break;
            case 'i':
                inputObjectSpectrum = optarg;
                break;
            case 'm':
                inputFlatFluxCalibration = optarg;
                break;
            case 'o':
                outputRVFile = optarg;
                break;
            case 'L':       // atlas of telluric lines
                telluric_lines = optarg;
                break;
            case 'T':       // stellar spectrum
                inputStellarSpectrum = optarg;
                break;
            case 'B':       // barycentric correction
                inputBarycentricCorrection = optarg;
                break;
            case 'R':
                spectralResolution = atof(optarg);
                break;
            case 'r':
                radialVelocitySearchRange = atof(optarg);
                break;
            case 's':
                radialVelocitySearchStep = atof(optarg);
                break;
            case 'x':
                XCorrelationThreshold = atof(optarg);
                break;
            case 'b':
                normalizationBinsize = atoi(optarg);
                break;
            case 'k':
                StarPlusSky = (atoi(optarg)?true:false);
                break;
            case 'f':
                useFitToFindMaximum = (atoi(optarg)?true:false);
                break;
            case 'H':
                if (strlen(optarg))
                    sscanf(optarg, "%d %lf %lf %lf %lf %lf %lf %lf", &expnum,&mjd,&exptime,&airmass,&outsideTemperature,&windspeed,&relativeHumididy,&atmoPressure);
                break;
            case 'O':
                ordernumber = atoi(optarg);
                break;
            case 'M':
                minorder = atoi(optarg);
                minorderprovided = true;
                break;
            case 'X':
                maxorder = atoi(optarg);
                maxorderprovided = true;
                break;
            case 'A':
                inputWavelengthMask = optarg;
                break;
            case 'P':
                xcorrsplotfilename = optarg;
                plot = 1;
                break;
            case 'Q':
                specplotfilename = optarg;
                plot = 1;
                break;
            case 'S':
                xcorrscriptfilename = optarg;
                break;
            case 'U':
                specscriptfilename = optarg;
                break;
            case 'D':
                xcorrdatafilename = optarg;
                break;
            case 'F':
                xcorrfitdatafilename = optarg;
                break;
            case 'G':
                specdatafilename = optarg;
                break;
            case 'p':
                plot = true;
                break;
            case 'v':
                verbose = true;
                break;
            case 'd':
                debug = true;
                break;
            case 't':
                trace = true;
                break;
            case 'h':
                printUsageSyntax(argv[0]);
                exit(EXIT_SUCCESS);
                break;
            case '?':
                printUsageSyntax(argv[0]);
                exit(EXIT_SUCCESS);
                break;
        }
    }
    
    /*Start the module here*/
    
    try {
        // we need an input wavelength calibration file ...
        if (inputWaveFile.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        // we need an input object spectrum file ...
        if (inputObjectSpectrum.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        // we need an output wavelength calibration file ...
        if (outputRVFile.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        // we need an input atlas of telluric lines ...
        if (telluric_lines.empty() && inputStellarSpectrum.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        if (inputWavelengthMask.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        if (inputBarycentricCorrection.empty()) {
            throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        
        if (verbose) {
            cout << "operaRadialVelocity: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaRadialVelocity: inputObjectSpectrum = " << inputObjectSpectrum << endl;
            cout << "operaRadialVelocity: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
            cout << "operaRadialVelocity: outputRVFile = " << outputRVFile << endl;
            cout << "operaRadialVelocity: telluric_lines =" << telluric_lines << endl;
            cout << "operaRadialVelocity: inputStellarSpectrum =" << inputStellarSpectrum << endl;
            cout << "operaRadialVelocity: spectralResolution =" << spectralResolution << endl;
            cout << "operaRadialVelocity: radialVelocitySearchRange =" << radialVelocitySearchRange << endl;
            cout << "operaRadialVelocity: radialVelocitySearchStep =" << radialVelocitySearchStep << endl;
            cout << "operaRadialVelocity: XCorrelationThreshold =" << XCorrelationThreshold << endl;
            cout << "operaRadialVelocity: normalizationBinsize =" << normalizationBinsize << endl;
            cout << "operaRadialVelocity: StarPlusSky = " << StarPlusSky << endl;
            cout << "operaRadialVelocity: useFitToFindMaximum = " << useFitToFindMaximum << endl;
            cout << "operaRadialVelocity: inputWavelengthMask = " << inputWavelengthMask << endl;
            cout << "operaRadialVelocity: inputBarycentricCorrection = " << inputBarycentricCorrection << endl;
            
            if(ordernumber != NOTPROVIDED) {
                cout << "operaRadialVelocity: ordernumber = " << ordernumber << endl;
            }
            if(plot) {
                cout << "operaRadialVelocity: xcorrsplotfilename = " << xcorrsplotfilename << endl;
                cout << "operaRadialVelocity: specplotfilename = " << specplotfilename << endl;
                cout << "operaRadialVelocity: xcorrscriptfilename = " << xcorrscriptfilename << endl;
                cout << "operaRadialVelocity: specscriptfilename = " << specscriptfilename << endl;
                cout << "operaRadialVelocity: xcorrdatafilename = " << xcorrdatafilename << endl;
                cout << "operaRadialVelocity: xcorrfitdatafilename = " << xcorrfitdatafilename << endl;
                cout << "operaRadialVelocity: specdatafilename = " << specdatafilename << endl;
            }
            
        }
        ofstream *fxcorrdata = NULL;
        ofstream *fxcorrfitdata = NULL;
        ofstream *fspecdata = NULL;
        
        if (!xcorrdatafilename.empty()) {
            fxcorrdata = new ofstream();
            fxcorrdata->open(xcorrdatafilename.c_str());
        }
        
        if (!xcorrfitdatafilename.empty()) {
            fxcorrfitdata = new ofstream();
            fxcorrfitdata->open(xcorrfitdatafilename.c_str());
        }
        
        if (!specdatafilename.empty()) {
            fspecdata = new ofstream();
            fspecdata->open(specdatafilename.c_str());
        }
        
        ofstream *foutRVfile = NULL;
        foutRVfile = new ofstream();
        foutRVfile->open(outputRVFile.c_str());
        
        operaSpectralOrderVector spectralOrders(inputObjectSpectrum);
        
        spectralOrders.ReadSpectralOrders(inputBarycentricCorrection);
        
        double barycentricRV = spectralOrders.getRadialVelocityCorrection();
        
        if (verbose) {
            cout << "operaRadialVelocity: Barycentric Radial Velocity =" << barycentricRV << " km/s" << endl;
        }
        
        spectralOrders.ReadSpectralOrders(inputWaveFile); // This merges in the wavelength calibration information
        
        if(!minorderprovided) {
            minorder = spectralOrders.getMinorder();
        }
        if(!maxorderprovided) {
            maxorder = spectralOrders.getMaxorder();
        }
        
        if(ordernumber != NOTPROVIDED) {
            minorder = ordernumber;
            maxorder = ordernumber;
        }
		      
        if (verbose)
        cout << "operaRadialVelocity: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        
        double integratedFlux, meanFlux;
        double maxSNR, meanSNR;
        
        spectralOrders.calculateRawFluxQuantities(minorder, maxorder,&integratedFlux,&meanFlux,&maxSNR,&meanSNR);
        
        //---------------------------------
        // Correct for flat-field
        if (!inputFlatFluxCalibration.empty()) {
            bool starplusskyInvertSkyFiber = false;
            spectralOrders.correctFlatField(inputFlatFluxCalibration, minorder, maxorder, StarPlusSky, starplusskyInvertSkyFiber);
        }
        //---------------------------------
        
        //---------------------------------
        // Get object spectrum within wavelength mask input ranges defined in "inputWavelengthMask"
        double *objectSpectrum = new double[MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
        double *objectSpectrumVariance = new double[MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
        double *wavelength = new double [MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
        
        unsigned nelem = spectralOrders.getSpectrumWithinTelluricMask(inputWavelengthMask, minorder, maxorder, TRUE, normalizationBinsize, wavelength, objectSpectrum, objectSpectrumVariance);
        
#ifdef PRINT_DEBUG
        for(unsigned i=0; i<nelem; i++) {
            cout << wavelength[i] << " " << objectSpectrum[i] << " " << objectSpectrumVariance[i] << endl;
        }
#endif
        //---------------------------------
        
        //---------------------------------
        // Define wavelength full range
        double wl0 = wavelength[0] - radialVelocitySearchRange*wavelength[0]/SPEED_OF_LIGHT_KMS;
        double wlf = wavelength[nelem-1] + radialVelocitySearchRange*wavelength[nelem-1]/SPEED_OF_LIGHT_KMS;
        if(verbose) {
            cout << "operaRadialVelocity: wl0=" <<  wl0 << " wlf=" << wlf << endl;
        }
        //---------------------------------
        
        /*
         * Read telluric lines database
         *		lambda vs. intensity
         */
        double maxatlasflux = -BIG;
        if (!telluric_lines.empty()) {
            if (debug) {
                cout << "operaRadialVelocity: reading telluric lines database " << telluric_lines << endl;
            }
            ntelluriclines = readTelluricLines(telluric_lines,telluricMoleculeNumber,telluricLinesWavelength,telluricLinesIntensity);
#ifdef PRINT_DEBUG
            for(unsigned line=0;line<ntelluriclines;line++) {
                cout << line << " " << telluricMoleculeNumber[line] << " " << telluricLinesWavelength[line] << " " << telluricLinesIntensity[line] << endl;
            }
#endif
            //---------------------------------
            // Find line with maximum absorption for normalization
            for (unsigned i=0; i<ntelluriclines; i++) {
                if((1-telluricLinesIntensity[i]) > maxatlasflux)
                maxatlasflux = (1-telluricLinesIntensity[i]);
            }
            //---------------------------------
        } else {
            maxatlasflux = 1;
        }
        
        //---------------------------------
        /* Read telluric reference files
         *
         * Read Telluric reference spectrum
         *		lambda vs. intensity
         */
        if (!inputStellarSpectrum.empty()) {
            nPointsInStellarSpectrum = readStellarSpectrum(inputStellarSpectrum, stellarSpectrumWavelength, stellarSpectrumIntensity);
        }
        //---------------------------------
        
        //---------------------------------
        // Calculate radial velocity shift by cross-correlation:
        double rvshift = 0;
        double rvshifterror = 0;
        double maxcorr = 0;
        double chisqr = 0;
        
        /*
         * Below is the main function to calculate teh radial velocity correction by cross-correlation
         * between observed spectrum and telluric reference.
         */
        if (verbose) {
            cout << "operaRadialVelocity: calculating cross-correlation for radialVelocitySearchRange=" << radialVelocitySearchRange << " km/s and radialVelocitySearchStep=" << radialVelocitySearchStep << " km/s" << endl;
        }

        double centralRV = 0.0;    // km/s

        bool validXCorrelation = calculateSourceRVShiftByXCorr(nelem, wavelength, objectSpectrum, radialVelocitySearchRange, radialVelocitySearchStep, XCorrelationThreshold, &rvshift, &rvshifterror, &maxcorr, fxcorrdata, fxcorrfitdata, spectralResolution, useFitToFindMaximum, &chisqr, centralRV);
        
        if(validXCorrelation) {
            if (verbose) {
                cout << "\noperaRadialVelocity: Radial Velocity correction = " << rvshift << " +/- " << rvshifterror << " km/s, maxXCorr=" << maxcorr << ", chisqr=" << chisqr << "\n" << endl;
            }
        } else {
            rvshift = 0;
            rvshifterror = 0;
            maxcorr = 0;
            chisqr = 0;
            cout << "\noperaRadialVelocity: Invalid cross-correlation. Radial Velocity correction = " << rvshift << " +/- " << rvshifterror << " km/s, maxXCorr=" << maxcorr << ", chisqr=" << chisqr << "\n" << endl;
        }
        //---------------------------------
        
        if(fspecdata!=NULL) {
            /*
             * Generate telluric synthetic spectrum using both HITRAN lines and input telluric spectrum:
             */
            //---------------------------------
            // First use hitran lines to generate synthetic spectrum
            double *hitranTelluricSpectrum = new double[nelem];
            if (!telluric_lines.empty()) {
                generateSyntheticTelluricSpectrumUsingGaussianProfile(nelem,wavelength,hitranTelluricSpectrum,spectralResolution);
            }
            //---------------------------------
            
            //---------------------------------
            // Then read telluric spectrum from input spectrum file.
            double *syntheticStellarSpectrum = new double[nelem];
            
            if (!inputStellarSpectrum.empty()) {
                double *wl,*stellarSpectrum;
                unsigned npointsInShortStellarSpectrum = getStellarSpectrumRange(wl0,wlf,&wl,&stellarSpectrum);
                double *shiftedwavelength = new double [npointsInShortStellarSpectrum];
                for (unsigned i=0; i<npointsInShortStellarSpectrum; i++) {
                    double DWavelength = rvshift * wl[i] / SPEED_OF_LIGHT_KMS;
                    shiftedwavelength[i] = wl[i] + DWavelength;
                }
                double *convolvedStellarSpectrum = new double [npointsInShortStellarSpectrum];
                convolveSpectrumWithGaussianByResolution(npointsInShortStellarSpectrum,shiftedwavelength,stellarSpectrum,convolvedStellarSpectrum,spectralResolution);
                
                // Spline to match same sampling as in observed spectrum
                operaFitSplineDouble(npointsInShortStellarSpectrum,shiftedwavelength,convolvedStellarSpectrum,nelem,wavelength,syntheticStellarSpectrum);
            }
            //---------------------------------
            
            //---------------------------------
            // Print spectral data to file
            for(unsigned i=0; i<nelem; i++) {
                *fspecdata << wavelength[i] << " " << objectSpectrum[i] << " " << hitranTelluricSpectrum[i] << " " << syntheticStellarSpectrum[i] << endl;
            }
            //---------------------------------
        }
        
        *foutRVfile << "#EXPNUM MJD RVSHIFT RVSHIFT_ERR BARYRVCORR INTEGRATEDFLUX MEANFLUX MAXSNR MEANSNR EXPTIME AIRMASS OUTTEMP WINDSPEED RELHUMID ATMPRESSURE XCORR" << endl;
        *foutRVfile << expnum << ' ';
        *foutRVfile << fixed << setprecision(7) << mjd << ' ';
        *foutRVfile << fixed << setprecision(4) << rvshift << ' ' << rvshifterror << ' ' << barycentricRV << ' ';
        *foutRVfile << scientific << setprecision(6) << integratedFlux << ' ' << meanFlux << ' ';
        *foutRVfile << fixed << setprecision(3) << maxSNR << ' ' << meanSNR << ' '
        << exptime << ' '
        << airmass << ' '
        << outsideTemperature << ' '
        << windspeed << ' '
        << relativeHumididy << ' '
        << atmoPressure << ' ' << maxcorr << endl;
        
        foutRVfile->close();
        
        /*
         * Telluric wavelength correction orders info plot:
         */
        if (fxcorrdata != NULL && fxcorrfitdata != NULL) {
            fxcorrdata->close();
            fxcorrfitdata->close();
            if (!xcorrscriptfilename.empty()) {
                GenerateTelluricXCorrelationPlot(xcorrscriptfilename, xcorrsplotfilename, xcorrdatafilename, xcorrfitdatafilename);
            }
        }
        
        /*
         * Spectrum plot: plot observed and reference telluric spectra.
         */
        
        if (fspecdata != NULL) {
            fspecdata->close();
            
            if (!specscriptfilename.empty()) {
                GenerateTelluricSpecPlot(specscriptfilename, specplotfilename, specdatafilename);
            }
        }
        
        //delete[] objectSpectrum;
        //delete[] wavelength;
        //delete[] objectSpectrumVariance;
        
    }
    catch (operaException e) {
        cerr << "operaRadialVelocity: " << e.getFormattedMessage() << endl;
        return EXIT_FAILURE;
    }
    catch (...) {
        cerr << "operaRadialVelocity: " << operaStrError(errno) << endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/*
 * Generate multiple plot containing statistical info about telluric wavelength correction
 */
void GenerateTelluricXCorrelationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string cleanDataFileName)
{
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    *fgnu << "reset" << endl;
    
    *fgnu << "\nset xlabel \"Radial Velocity (km/s)\"" << endl;
    *fgnu << "set ylabel \"cross-correlation\"" << endl;
    *fgnu << "set pointsize 1.5" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        *fgnu << "plot \"" << dataFileName << "\" u 1:2 t \"XCorr data\" w p pt 6 ps 0.75, ";
        *fgnu << "\"" << cleanDataFileName << "\" u 1:2 t \"gaussian fit\" w l lt 3 lw 2, ";
        *fgnu << "\"" << cleanDataFileName << "\" u 1:3:4 t \"fit data\" w yerr pt 7 ps 0.75" << endl;
        
        *fgnu << "\n#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
        *fgnu << "#replot" << endl;
    } else {
        *fgnu << endl;
        
        *fgnu << "plot \"" << dataFileName << "\" u 1:2 t \"XCorr data\" w p pt 6 ps 0.75, ";
        *fgnu << "\"" << cleanDataFileName << "\" u 1:2 t \"gaussian fit\" w l lt 3 lw 2, ";
        *fgnu << "\"" << cleanDataFileName << "\" u 1:3:4 t \"fit data\" w yerr pt 7 ps 0.75" << endl;
        
        *fgnu << endl;
        
        *fgnu << "\n#set terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    
    if(!outputPlotEPSFileName.empty())
    systemf("gnuplot %s",gnuScriptFileName.c_str());
}


/*
 * Generate 2D plot for spectra of atlas + comparison + identified lines
 */
void GenerateTelluricSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string specdatafilename) {
    
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "#unset key" << endl;
    
    *fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
    
    *fgnu << "set ylabel \"norm flux\"" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        *fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 3, ";
        *fgnu << "\"" << specdatafilename << "\" u 1:($3*$4) t \"Synthetic Spectrum (Telluric*Stellar)\" w l lt 1 lw 2";

        *fgnu << endl;
        
        *fgnu << "\n#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
        *fgnu << "#replot" << endl;
    } else {
        *fgnu << endl;
        
        *fgnu << "plot \"" << specdatafilename << "\" u 1:2 t \"Object Spectrum\" w l lt 3, ";
        *fgnu << "\"" << specdatafilename << "\" u 1:($3*$4) t \"Synthetic Spectrum (Telluric*Stellar)\" w l lt 1 lw 2";
        *fgnu << endl;
        *fgnu << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if(!outputPlotEPSFileName.empty())
    systemf("gnuplot %s",gnuScriptFileName.c_str());
}


/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
    cout <<
    "\n"
    " Usage: "+string(modulename)+"  [-vdth]" +
    " --inputWaveFile=<WAVE_FILE>"
    " --inputObjectSpectrum=<SPECTRUM_FILE>"
    " --outputRVFile=<WAVE_FILE>"
    " --telluric_lines=<LINES_FILE>"
    " --inputStellarSpectrum=<SPECTRUM_FILE>"
    " --spectralResolution=<DBL_VALUE>"
    " --radialVelocitySearchRange=<DBL_VALUE>"
    " --radialVelocitySearchStep=<DBL_VALUE>"
    " --XCorrelationThreshold=<DBL_VALUE>"
    " --normalizationBinsize=<UNS_VALUE>"
    " --StarPlusSky=<BOOL>"
    " --useFitToFindMaximum=<BOOL>"
    " --ordernumber=<INT_VALUE"
    " --minorder=<INT_VALUE>"
    " --maxorder=<INT_VALUE>"
    " --xcorrsplotfilename=<EPS_FILE>"
    " --specplotfilename=<EPS_FILE>"
    " --xcorrscriptfilename=<GNUPLOT_FILE>"
    " --specscriptfilename=<GNUPLOT_FILE>"
    " --xcorrdatafilename=<DATA_FILE>"
    " --xcorrfitdatafilename=<DATA_FILE>"
    " --specdatafilename=<DATA_FILE>\n\n"
    " Example: "+string(modulename)+" --inputObjectSpectrum= ... -v \n\n"
    "  -h, --help  display help message\n"
    "  -v, --verbose,  Turn on message sending\n"
    "  -d, --debug,  Turn on debug messages\n"
    "  -t, --trace,  Turn on trace messages\n"
    "  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file \n"
    "  -s, --inputObjectSpectrum=<SPECTRUM_FILE>, Input object spectrum file (.s)\n"
    "  -o, --outputRVFile=<WAVE_FILE>, Output wavelength calibration file to store final solution\n"
    "  -L, --telluric_lines=<LINES_FILE>, Atlas of telluric lines \n"
    "  -T, --inputStellarSpectrum=<SPECTRUM_FILE>, Spectrum of telluric lines\n"
    "  -R, --spectralResolution=<DBL_VALUE>, Input spectral resolution (wl/dwl) as reference for line detection\n"
    "  -r, --radialVelocitySearchRange=<DBL_VALUE>, Radial velocity range (in nm) to scan for first order correction\n"
    "  -s, --radialVelocitySearchStep=<DBL_VALUE>, Radial velocity step (in nm) to scan for first orer correction\n"
    "  -x, --XCorrelationThreshold=<DBL_VALUE>, X-correlation lower threshold to consider a match between telluric and object spectra\n"
    "  -b, --normalizationBinsize=<UNS_VALUE>, Binsize to normalize input object spectrum \n"
    "  -k, --StarPlusSky=<BOOL>, Star plus sky mode \n"
    "  -f, --useFitToFindMaximum=<BOOL>, Use gaussian fit to find maximum xcorr and corresponding RV \n"
    "  -O, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
    "  -N, --minorder=<INT_VALUE>, Define minimum order number\n"
    "  -X, --maxorder=<INT_VALUE>, Define maximum order number\n"
    "  -S, --xcorrscriptfilename=<GNUPLOT_FILE>\n"
    "  -U, --specscriptfilename=<GNUPLOT_FILE>\n"
    "  -P, --xcorrsplotfilename=<EPS_FILE>\n"
    "  -Q, --specplotfilename=<EPS_FILE>\n"
    "  -D, --xcorrdatafilename=<DATA_FILE>\n"
    "  -F, --xcorrfitdatafilename=<DATA_FILE>\n"
    "  -G, --specdatafilename=<DATA_FILE>\n\n";
}

/*
 * Read the entire set of telluric lines in HITRAN database
 */

unsigned readTelluricLines(string telluric_database_file, int *telluricMoleculeNumber, double *telluricLinesWavelength, double *telluricLinesIntensity) {
   	igzstream astream;
    string dataline;
    int tmpnumber = -1;
    double tmpwn = -1.0;
    float tmpi = -1.0;
    unsigned line = 0;
    char *buff = new char[MAXLENGTHOFLINEINTELLURICDATABASE];
    
    int *tmp_telluricMoleculeNumber = new int[MAXNUMBEROFLINESINTELLURICDATABASE];
    double *tmp_telluricLinesWavelength = new double[MAXNUMBEROFLINESINTELLURICDATABASE];
    double *tmp_telluricLinesIntensity = new double[MAXNUMBEROFLINESINTELLURICDATABASE];
    
    astream.open(telluric_database_file.c_str());
    if (astream.is_open()) {
        while (astream.good()) {
            getline(astream, dataline);
            if (strlen(dataline.c_str())) {
                if (dataline.c_str()[0] == '#') {
                    // skip comments
                } else {
                    sscanf(dataline.c_str(), "%d %lf %G %[^\n]", &tmpnumber, &tmpwn, &tmpi, buff);
                    
                    //tmpi = tmpi*pow(10,tmpexp);
                    double wavelength_in_nm = 1e7/tmpwn;
                    
                    tmp_telluricMoleculeNumber[line] = tmpnumber;
                    
                    tmp_telluricLinesWavelength[line] = convertVacuumToAirWavelength(wavelength_in_nm*10)/10;
                    
                    double NoverV = TYPICAL_PRESSURE_AT_MAUNAKEA/(TYPICAL_TEMPERATURE_AT_MAUNAKEA*k_BOLTZMANN_CONSTANT);
                    tmp_telluricLinesIntensity[line] = ((double)tmpi/(NoverV*1e-6))/TYPICAL_ATMOSPHERE_PATH_LENGTH;
                    
                    //                    printf("%u %d %lf %lf\n",line,tmpnumber,wavelength_in_nm,(double)tmpi*ONEMOLE);
                    
                    line++;
                }	// skip comments
            }	// if strlen
        } // while (astream.good())
        line--;
        if (line > 0) {
            if (verbose) {
                printf("          [Telluric] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", line, tmp_telluricLinesWavelength[0], tmp_telluricLinesWavelength[line/2], tmp_telluricLinesWavelength[line-1]);
            }
        } else {
            printf("          [Telluric] no lines found in telluric database.\n");
        }
        astream.close();
    }	// if (astream.open()
    delete[] buff;
    
    unsigned np = line;
    for(unsigned i=0; i<line; i++) {
        np--;
        telluricMoleculeNumber[i] = tmp_telluricMoleculeNumber[np];
        telluricLinesWavelength[i] = tmp_telluricLinesWavelength[np];
        telluricLinesIntensity[i] = tmp_telluricLinesIntensity[np];
        //        cout << "operaRadialVelocity: " << telluricMoleculeNumber[i] << " " << telluricLinesWavelength[i] << " " << telluricLinesIntensity[i] << endl;
    }
    
    delete[] tmp_telluricMoleculeNumber;
    delete[] tmp_telluricLinesWavelength;
    delete[] tmp_telluricLinesIntensity;
    
    return line;
}

/*
 * get a subset of the telluric lines for this order only, between wl0 and wlf
 */
unsigned getTelluricLinesRange(double wl0, double wlf, double **wl, double **intensity) {
    
    unsigned firstline = 0;
    unsigned line = 0;
    
    for (line=0; line<ntelluriclines; line++) {
        if (telluricLinesWavelength[line] >= wl0) {
            if (firstline == 0) {
                *intensity = &telluricLinesIntensity[line];
                *wl = &telluricLinesWavelength[line];
                firstline = line;
            }
            if (telluricLinesWavelength[line] > wlf)
            break;
        }
    }
    if (line)
    line--;
    if (line > firstline) {
        return (line-firstline);
    } else {
        return 0;
    }
}

/*
 * Read the the full atmospheric transmission spectrum
 */
unsigned readStellarSpectrum(string inputStellarSpectrum, double *stellarSpectrumWavelength, double *stellarSpectrumIntensity) {
    igzstream astream;
    string dataline;
    
    double tmpwl = -1.0;
    double tmpi = -1.0;
    double tmpni = -1.0;
    unsigned np = 0;
    
    astream.open(inputStellarSpectrum.c_str());
    if (astream.is_open()) {
        while (astream.good()) {
            getline(astream, dataline);
            if (strlen(dataline.c_str())) {
                if (dataline.c_str()[0] == '#') {
                    // skip comments
                } else {
                    sscanf(dataline.c_str(), " %lf %lf %lf ", &tmpwl, &tmpi, &tmpni);
                    
                    //double DWavelength = barycentricRV * (tmpwl/10.0)  / SPEED_OF_LIGHT_KMS;
                    
                    stellarSpectrumWavelength[np] = tmpwl/10.0;
                    stellarSpectrumIntensity[np] = tmpni;
                    np++;
                }	// skip comments
            }
        } // while (astream.good())
        
        if (np > 0) {
            if (verbose) {
                printf("          [Stellar Synthetic Spetrum] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, stellarSpectrumWavelength[0], stellarSpectrumWavelength[np/2], stellarSpectrumWavelength[np-1]);
            }
        } else {
            printf("          [Stellar Synthetic Spetrum] no points found in telluric spectrum.\n");
        }
        astream.close();
    }	// if (astream.open()
    return np;
}

/*
 * get a subset of the telluric lines for this order only, between wl0 and wlf
 */
unsigned getStellarSpectrumRange(double wl0, double wlf, double **wl, double **intensity) {
    
    unsigned firstline = 0;
    unsigned line = 0;
    
    for (line=0; line<nPointsInStellarSpectrum; line++) {
        if (stellarSpectrumWavelength[line] >= wl0) {
            if (firstline == 0) {
                *intensity = &stellarSpectrumIntensity[line];
                *wl = &stellarSpectrumWavelength[line];
                firstline = line;
            }
            if (stellarSpectrumWavelength[line] > wlf)
            break;
        }
    }
    if (line)
    line--;
    if (line > firstline) {
        return (line-firstline);
    } else {
        return 0;
    }
}

unsigned matchTelluricReferencewithObjectLines(double acceptableMismatch,double lineSigma, operaSpectralLines *telluricLines, operaSpectralLines *objectLines, Polynomial *wlcorrection, unsigned order, double wl_central, ofstream *flinesdata) {
    /*
     * vectors for uncalibrated spectral lines
     */
    if (objectLines->getnLines() == 0) {
        throw operaException("operaRadialVelocity: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    if (telluricLines->getnLines() == 0) {
        throw operaException("operaRadialVelocity: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    double *objectlinecenter = new double[objectLines->getnLines()];
    double *objectlinecenterError = new double[objectLines->getnLines()];
    double *objectlineflux = new double[objectLines->getnLines()];
    double *objectlinesigma = new double[objectLines->getnLines()];
    
    double *objectMatchedData = new double[objectLines->getnLines()];
    unsigned *objectMatchedindex = new unsigned[objectLines->getnLines()];
    double *telluricMatchedData = new double[telluricLines->getnLines()];
    unsigned *telluricMatchedindex = new unsigned[telluricLines->getnLines()];
    
    /*
     * vectors for telluric spectral lines
     */
    double *telluricLineswl = new double[telluricLines->getnLines()];
    double *telluricLineswlError = new double[telluricLines->getnLines()];
    double *telluricLinesflux = new double[telluricLines->getnLines()];
    
    unsigned line = 0;
    
    for(unsigned feature=0;feature<objectLines->getNFeatures();feature++) {
        operaSpectralFeature *currentFeature = objectLines->getSpectralFeature(feature);
        double *center = currentFeature->getGaussianFit()->getCenterVector();
        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector();
        for(unsigned l=0; l<currentFeature->getnLines(); l++) {
            objectlinecenter[line] = center[l];
            objectlinecenterError[line] = centerError[l];
            objectlineflux[line] = amplitude[l];
            objectlinesigma[line] = sigma[l];
            if(debug)
            cout << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << endl;
            line++;
        }
    }
    line = 0;
    
    for(unsigned feature=0;feature<telluricLines->getNFeatures();feature++) {
        operaSpectralFeature *currentFeature = telluricLines->getSpectralFeature(feature);
        double *center = currentFeature->getGaussianFit()->getCenterVector();
        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector();
        for(unsigned l=0; l<currentFeature->getnLines(); l++) {
            telluricLineswl[line] = center[l];
            telluricLineswlError[line] = centerError[l];
            telluricLinesflux[line] = amplitude[l];
            if(debug)
            cout << center[l] <<  " " << amplitude[l] << " " << sigma[l] << endl;
            line++;
        }
    }
    
    /*
     * Below it identifies and select the set of lines that match both comparison and atlas.
     * The criteria for matching is that the difference between centers must be < acceptableMismatch x sigma
     */
    unsigned nmatch = 0;
    unsigned nextfirstline = 0;
    
    double acceptMismatchInwlUnits = acceptableMismatch*lineSigma;
    
    for (unsigned i=0; i<objectLines->getnLines(); i++) {
        
        unsigned bestAtlasMatchIndex = 0;
        double mindifference = acceptMismatchInwlUnits;
        
        for(unsigned l=nextfirstline;l<telluricLines->getnLines();l++) {
            
            double difference = fabs(objectlinecenter[i] - telluricLineswl[l]);
#ifdef PRINT_DEBUG
            cout << "operaRadialVelocity: " << "mindiff=" << mindifference << " diff=" << difference << " object[" << i << "]=" << objectlinecenter[i] << " telluric[" << l << "]=" << telluricLineswl[l] << endl;
#endif
            if(objectlinecenter[i] > telluricLineswl[l]  && difference < mindifference) {
                mindifference = difference;
                bestAtlasMatchIndex = l;
            } else if (objectlinecenter[i] <= telluricLineswl[l] && difference < mindifference) {
                objectMatchedData[nmatch] = objectlinecenter[i];
                telluricMatchedData[nmatch] = telluricLineswl[l];
                telluricMatchedindex[nmatch] = l;
                objectMatchedindex[nmatch] = i;
                nextfirstline = l+1;
                nmatch++;
                break;
            } else if (objectlinecenter[i] <= telluricLineswl[l]  && difference > mindifference) {
                if(bestAtlasMatchIndex) {
                    objectMatchedData[nmatch] = objectlinecenter[i];
                    telluricMatchedData[nmatch] = telluricLineswl[bestAtlasMatchIndex];
                    telluricMatchedindex[nmatch] = bestAtlasMatchIndex;
                    objectMatchedindex[nmatch] = i;
                    nextfirstline = bestAtlasMatchIndex+1;
                    nmatch++;
                }
                break;
            }
            if(l==telluricLines->getnLines()-1) {
                if(bestAtlasMatchIndex) {
                    objectMatchedData[nmatch] = objectlinecenter[i];
                    telluricMatchedData[nmatch] = telluricLineswl[bestAtlasMatchIndex];
                    telluricMatchedindex[nmatch] = bestAtlasMatchIndex;
                    objectMatchedindex[nmatch] = i;
                    nextfirstline = bestAtlasMatchIndex+1;
                    nmatch++;
                }
            }
        }
    }
    
    if (nmatch > 0) {
        if(flinesdata != NULL){
            for(unsigned index=0; index<nmatch;index++){
                *flinesdata << order << " " << objectMatchedData[index] << " " << telluricMatchedData[index] << " " << objectlineflux[objectMatchedindex[index]] << " " << objectlinesigma[objectMatchedindex[index]] << " " << wl_central << endl;
            }
            *flinesdata << endl;
        }
        if (wlcorrection) {
            int npar = wlcorrection->getOrderOfPolynomial();
            double *par = (double *)wlcorrection->getVector();
            double chisqr;
            
            operaLMFitPolynomial(nmatch, telluricMatchedData, objectMatchedData, npar, par, &chisqr);
            
            wlcorrection->setChisqr(chisqr);
        }
   	}
    
    delete[] objectMatchedData;
    delete[] objectMatchedindex;
    delete[] telluricMatchedData;
    delete[] telluricMatchedindex;
    
    delete[] objectlinecenter;
    delete[] objectlinecenterError;
    delete[] objectlineflux;
    delete[] objectlinesigma;
    delete[] telluricLineswl;
    delete[] telluricLineswlError;
    delete[] telluricLinesflux;
    
    return nmatch;
}

bool calculateSourceRVShiftByXCorr(unsigned nelem, double *wavelength, double *objectSpectrum, double radialVelocityRange, double radialVelocityStep, double threshold, double *maxRV, double *sigRV, double *maxcorr, ostream *fxcorrdata, ostream *fxcorrfitdata, double spectralResolution, bool useFitToFindMaximum, double *chisqr, double centralRV) {
    
    bool status = true;
    
    double wl0 = wavelength[0] * (1 - radialVelocityRange / SPEED_OF_LIGHT_KMS);
    double wlf = wavelength[nelem-1] * (1 + radialVelocityRange / SPEED_OF_LIGHT_KMS);
    
    double *telluricSpectrum = new double[nelem];
    generateSyntheticTelluricSpectrumUsingGaussianProfile(nelem,wavelength,telluricSpectrum,spectralResolution);

    double *inputStellarWavelength, *inputStellarSpectrum;
    unsigned npointsInInputStellarSpectrum = getStellarSpectrumRange(wl0,wlf,&inputStellarWavelength,&inputStellarSpectrum);
    double *stellarShiftedWavelength = new double[npointsInInputStellarSpectrum];
    double *convolvedStellarSpectrum = new double[npointsInInputStellarSpectrum];
    convolveSpectrumWithGaussianByResolution(npointsInInputStellarSpectrum,inputStellarWavelength,inputStellarSpectrum,convolvedStellarSpectrum,spectralResolution);
    
    double *syntheticStellarSpectrum = new double[nelem];

    double *syntheticSpectrum = new double[nelem];
    
    unsigned nDataPoints = (unsigned)ceil(radialVelocityRange/radialVelocityStep);
    double firstRV = centralRV - radialVelocityRange/2.0;
    double deltaRV = firstRV;
    unsigned jmax = 0;
    *maxcorr = -BIG;
    *maxRV = 0;
    
    double *crosscorrelation = new double [nDataPoints];
    double *crosscorrerror = new double [nDataPoints];
    double *dRV = new double [nDataPoints];
    double xcorrerror = 2e-04;

    for(unsigned j=0; j<nDataPoints;j++) {
        
        for (unsigned i=0; i<npointsInInputStellarSpectrum; i++) {
            double DWavelength = deltaRV * inputStellarWavelength[i] / SPEED_OF_LIGHT_KMS;
            stellarShiftedWavelength[i] = inputStellarWavelength[i] + DWavelength;
        }
        
        // Spline to match same sampling as in observed spectrum
        operaFitSplineDouble(npointsInInputStellarSpectrum,stellarShiftedWavelength,convolvedStellarSpectrum,nelem,wavelength,syntheticStellarSpectrum);
        
        for (unsigned i=0; i<nelem; i++) {
            syntheticSpectrum[i] = telluricSpectrum[i] * syntheticStellarSpectrum[i];
        }

        crosscorrelation[j] = operaCrossCorrelation(nelem,objectSpectrum,syntheticSpectrum);
        crosscorrerror[j] = xcorrerror;
        dRV[j] = deltaRV;
        
        if(debug) {
            cout << dRV[j] << " " << crosscorrelation[j] << endl;
        }
        // Test (I) :
        if(crosscorrelation[j] > *maxcorr && crosscorrelation[j] > threshold) {
            *maxcorr = crosscorrelation[j];
            *maxRV = deltaRV;
            *sigRV = radialVelocityStep;
            jmax = j;
        }
        
        deltaRV+=radialVelocityStep;
    }
    
    if (*maxRV == firstRV) {
        *maxcorr = 0;
        *maxRV = 0;
        *sigRV = 0;
        *chisqr = 0;
        status = false;
    }
    
    if(useFitToFindMaximum && status == TRUE) {
        
        double *peakRV = new double[nDataPoints];
        double *peakXCorr = new double[nDataPoints];
        double *peakXCorrErr = new double[nDataPoints];
        
        unsigned np = 0;
        
        double RVFitRange = (SPEED_OF_LIGHT_KMS/spectralResolution);
        double startRV = dRV[jmax] - 2*RVFitRange/3;
        double endRV = dRV[jmax] + 2*RVFitRange/3;
        
        for (unsigned j=0; j<nDataPoints; j++) {
            if(dRV[j] >= startRV && dRV[j] <= endRV) {
                peakRV[np] = dRV[j];
                peakXCorr[np] = crosscorrelation[j];
                np++;
            }
        }
        
        double a=(double)crosscorrelation[jmax];
        double x0=(double)dRV[jmax];
        double sig=(double)RVFitRange/4;
        double ea;
        double ex0;
        double esig;
        double firstchisqr;
        
        operaLMFitGaussian(np, peakRV, peakXCorr, &a, &x0, &sig, &firstchisqr);
        
        double rmsResiduals = 0;
        for (unsigned j=0; j<np; j++) {
            double x = (double)peakRV[j];
            double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
            rmsResiduals += ((peakXCorr[j] - gaussfunc)*(peakXCorr[j] - gaussfunc))/(double)np;
        }
        rmsResiduals = sqrt(rmsResiduals);
        for (unsigned j=0; j<np; j++) {
            peakXCorrErr[j] = rmsResiduals;
        }
        
        double fitchisqr = 0;
        operaMPFitGaussian(np, peakRV, peakXCorr, peakXCorrErr, &a, &ea, &x0, &ex0, &sig, &esig, &fitchisqr);

        if(debug) {
            cout << a << "+/-" << ea << endl;
            cout << x0 << "+/-" << ex0 << endl;
            cout << sig << "+/-" << esig <<  " fitchisqr=" << fitchisqr << endl;
        }
        
        // Below is for plotting
        if(fxcorrdata != NULL) {
            for(unsigned j=0; j<nDataPoints;j++) {
                *fxcorrdata << dRV[j] << " " <<   crosscorrelation[j] << " " <<  crosscorrerror[j] << endl;
            }
            *fxcorrdata << endl;
        }
        
        if(fxcorrdata != NULL) {
            for(unsigned j=0; j<np;j++) {
                double x = (double)peakRV[j];
                double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
                
                *fxcorrfitdata << peakRV[j] << " " <<  gaussfunc << " " <<  peakXCorr[j] << " " <<  peakXCorrErr[j] << " " << peakXCorr[j] - gaussfunc << endl;
            }
            *fxcorrfitdata << endl;
        }
        
        *maxcorr = a;
        *maxRV = x0;
        *sigRV = ex0;
        *chisqr = fitchisqr;
        
    } else if (!useFitToFindMaximum && status == TRUE) {
        
        // Below is for plotting
        if(fxcorrdata != NULL) {
            for(unsigned j=0; j<nDataPoints;j++) {
                *fxcorrdata  << dRV[j] << " " <<  crosscorrelation[j] << " " <<  crosscorrerror[j] << endl;
            }
            *fxcorrdata << endl;
        }
        
        // Below is for plotting
        if(fxcorrfitdata != NULL) {
            *fxcorrfitdata  << *maxRV << " " << *maxcorr << " " <<  *maxcorr  <<  " " <<  crosscorrerror[jmax] <<  " " << 0.0 << endl;
        }
        *chisqr = 0;
    }
    
    delete[] crosscorrelation;
    delete[] crosscorrerror;
    delete[] dRV;
    
    delete[] stellarShiftedWavelength;
    delete[] telluricSpectrum;
    delete[] syntheticStellarSpectrum;
    delete[] syntheticSpectrum;

    return status;
}

void generateSyntheticTelluricSpectrumUsingGaussianProfile(unsigned np, double *wavelengthVector, double *ouputSpectrum, double resolution) {
    
    for(unsigned i=0;i<np;i++) {
        ouputSpectrum[i] = 1.0;
    }
    
    double *wl, *intensity;
    
    for(unsigned i=0;i<np;i++) {
        double gaussianWidth = (wavelengthVector[i]/resolution);
        
        unsigned nlinesInRange = getTelluricLinesRange(wavelengthVector[i] - 5*gaussianWidth,wavelengthVector[i] + 5*gaussianWidth,&wl,&intensity);
        
        if(debug) {
            cout << "operaRadialVelocity: " << i <<
            " gaussianwidth=" << gaussianWidth <<
            " wl0=" << wavelengthVector[i] - 5*gaussianWidth <<
            " wlf=" << wavelengthVector[i] + 5*gaussianWidth <<
            " nlinesInRange=" << nlinesInRange << endl;
        }
        for(unsigned j=0; j<nlinesInRange; j++) {
            double opticaldepth = intensity[j]*exp(-((wl[j] - wavelengthVector[i])*(wl[j] - wavelengthVector[i])/(2*gaussianWidth*gaussianWidth)))/(sqrt(2*M_PI)*gaussianWidth);
            ouputSpectrum[i] *= exp(-opticaldepth);
        }
    }
}

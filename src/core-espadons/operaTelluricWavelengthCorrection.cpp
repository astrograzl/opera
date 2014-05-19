/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaTelluricWavelengthCorrection
 Version: 1.0
 Description: Apply wavelength correction based on telluric lines
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
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

#include "core-espadons/operaTelluricWavelengthCorrection.h"
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
#define MINELEMENTS 20

/*! \file operaTelluricWavelengthCorrection.cpp */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*
 * the reference Telluric spectrum
 */
static unsigned nPointsInTelluricSpectrum = 0;
static double telluricSpectrumWavelength[MAXNUMBEROFPOINTSINTELLURICSPECTRUM];
static double telluricSpectrumIntensity[MAXNUMBEROFPOINTSINTELLURICSPECTRUM];

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
 * operaTelluricWavelengthCorrection
 * \author Eder Martioli
 * \brief Calculate and apply wavelength correction based on telluric lines.
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
	string outputWaveFile;
    string telluric_lines;
    string telluric_spectrum;
    /*
     * Parameters for telluric wavelength correction
     */

	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;
    bool maxorderprovided = false;
    
    bool interactive = false;
    
	//bool debug=false, verbose=false, trace=false;
    
    bool plot=false;
    
    string xcorrsplotfilename;
    string specplotfilename;
	
    string xcorrscriptfilename;
    string specscriptfilename;
    
    string atlasdatafilename;
	string compdatafilename;
	string linesdatafilename;
	string xcorrdatafilename;
    string xcorrfitdatafilename;
    
    bool subtractCentralWavelength = TRUE;
    
    bool apply2ndOrderCorrection = FALSE;
    /*
     * The parameters below we don't know yet whether they would be useful if used as input
     */
    double DetectionThreshold = 0.03;    // threshold to regulate the sensitivity of line detection. Must be between 0 and 1.
    double LocalMaxFilterWidth = 4.0;    // parameter to set a window filter to guarantee a line is not detected twice. It's in units of line width
    double MinPeakDepth = 0.1;           // limit that also regulates the sensitity of line detection in units of noise.
    
    double spectralResolution = 80000; // Input spectral resolution as reference for line detection
    
    double initialWavelengthRange = 0.1;
    double initialWavelengthStep = 0.0001;
    double XCorrelationThreshold = 0.05;
    float sigmaThreshold = 1.0;
    
    unsigned normalizationBinsize = 110;
    
    
    /*
     * halfslitsize below sets a window of valid points.
     * A gaussian will be fit to points within the window in order to
     * obtain the maximum cross-correlation.
     */
    int halfslitsize = 5;
    
	struct option longopts[] = {
		{"inputWaveFile",1, NULL, 'w'},				// input wavelength calibration file (.wcal)
 		{"inputObjectSpectrum",1, NULL, 'i'},		// input object spectrum file (.e or .p)
        {"outputWaveFile",1, NULL, 'o'},			// output wavelength calibration file (.auto or .pauto)
		{"telluric_lines",1, NULL, 'L'},            // atlas of telluric lines
		{"telluric_spectrum",1, NULL, 'T'},         // spectrum of telluric lines
		{"spectralResolution",1, NULL, 'R'},        // Input spectral resolution (wl/dwl) as reference for line detection
		{"initialWavelengthRange",1, NULL, 'r'},    // Wavelength shift range (in nm) to scan for first order correction
		{"initialWavelengthStep",1, NULL, 's'},     // Wavelength step (in nm) to scan for first orer correction
		{"XCorrelationThreshold",1, NULL, 'x'},     // X-correlation lower threshold to consider a match between telluric and object spectra
		{"sigmaThreshold",1, NULL, 'g'},            // X-correlation lower threshold in sigma units
		{"normalizationBinsize",1, NULL, 'b'},      // normalization binsize
        {"subtractCentralWavelength",1, NULL, 'K'},
        {"ordernumber",			1, NULL, 'O'},
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'},

        {"xcorrsplotfilename",1, NULL, 'P'},
        {"specplotfilename",1, NULL, 'Q'},
        {"xcorrscriptfilename",1, NULL, 'S'},
        {"specscriptfilename",1, NULL, 'U'},
        {"xcorrdatafilename",1, NULL, 'D'},
        {"xcorrfitdatafilename",1, NULL, 'F'},      
        {"atlasdatafilename",1, NULL, 'G'},
        {"compdatafilename",1, NULL, 'H'},
        {"linesdatafilename",1, NULL, 'J'},
        
		{"interactive",0, NULL, 'I'},
		
		{"plot",0, NULL, 'p'},
		{"verbose",0, NULL, 'v'},
		{"debug",0, NULL, 'd'},
		{"trace",0, NULL, 't'},
		{"help",0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:w:s:o:O:L:T:R:r:s:x:g:b:M:X:P:Q:D:F:G:H:J:S:U:K:p::v::d::t::h",
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
			case 'o':
				outputWaveFile = optarg;
				break;
			case 'L':       // atlas of telluric lines
				telluric_lines = optarg;
				break;
			case 'T':       // spectrum of telluric lines
				telluric_spectrum = optarg;
				break;
			case 'R':
				spectralResolution = atof(optarg);
				break;
			case 'r':
				initialWavelengthRange = atof(optarg);
				break;
			case 's':
				initialWavelengthStep = atof(optarg);
				break;
			case 'x':
				XCorrelationThreshold = atof(optarg);
				break;
			case 'g':
				sigmaThreshold = atof(optarg);
				break;
			case 'b':
				normalizationBinsize = atoi(optarg);
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
            case 'P':
				xcorrsplotfilename = optarg;
				plot = 1;
				break;
            case 'Q':
				specplotfilename = optarg;
				plot = 1;
				break;
			case 'D':
				xcorrdatafilename = optarg;
				break;
			case 'F':
				xcorrfitdatafilename = optarg;
				break;
			case 'G':
				atlasdatafilename = optarg;
				break;
			case 'H':
				compdatafilename = optarg;
				break;
			case 'J':
				linesdatafilename = optarg;
				break;
			case 'S':
				xcorrscriptfilename = optarg;
				break;
			case 'U':
				specscriptfilename = optarg;
				break;
			case 'K':
                subtractCentralWavelength = (atoi(optarg)?true:false);
				break;                
			case 'I':		// for interactive plots
				interactive = true;
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
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input object spectrum file ...        
		if (inputObjectSpectrum.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output wavelength calibration file ...
		if (outputWaveFile.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input atlas of telluric lines ...
		if (telluric_lines.empty() && telluric_spectrum.empty()) {
			throw operaException("operaTelluricWavelengthCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (verbose) {
			cout << "operaTelluricWavelengthCorrection: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaTelluricWavelengthCorrection: inputObjectSpectrum = " << inputObjectSpectrum << endl;
			cout << "operaTelluricWavelengthCorrection: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaTelluricWavelengthCorrection: telluric_lines =" << telluric_lines << endl;
			cout << "operaTelluricWavelengthCorrection: telluric_spectrum =" << telluric_spectrum << endl;
			cout << "operaTelluricWavelengthCorrection: spectralResolution =" << spectralResolution << endl;
			cout << "operaTelluricWavelengthCorrection: initialWavelengthRange =" << initialWavelengthRange << endl;
			cout << "operaTelluricWavelengthCorrection: initialWavelengthStep =" << initialWavelengthStep << endl;
			cout << "operaTelluricWavelengthCorrection: XCorrelationThreshold =" << XCorrelationThreshold << endl;
			cout << "operaTelluricWavelengthCorrection: sigmaThreshold =" << sigmaThreshold << endl;
			cout << "operaTelluricWavelengthCorrection: normalizationBinsize =" << normalizationBinsize << endl;
            
            if(ordernumber != NOTPROVIDED) {
                cout << "operaTelluricWavelengthCorrection: ordernumber = " << ordernumber << endl;
            }
            if(plot) {
                cout << "operaTelluricWavelengthCorrection: xcorrsplotfilename = " << xcorrsplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specplotfilename = " << specplotfilename << endl;
                cout << "operaTelluricWavelengthCorrection: xcorrscriptfilename = " << xcorrscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: specscriptfilename = " << specscriptfilename << endl;
                cout << "operaTelluricWavelengthCorrection: xcorrdatafilename = " << xcorrdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: xcorrfitdatafilename = " << xcorrfitdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: atlasdatafilename = " << atlasdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: compdatafilename = " << compdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: linesdatafilename = " << linesdatafilename << endl;
                cout << "operaTelluricWavelengthCorrection: subtractCentralWavelength = " << subtractCentralWavelength << endl;
                if(interactive) {
                    cout << "operaTelluricWavelengthCorrection: interactive = YES" << endl;
                } else {
                    cout << "operaTelluricWavelengthCorrection: interactive = NO" << endl;
                }
            }
            
		}
		ofstream *fatlasdata = NULL;
		ofstream *fcompdata = NULL;
		ofstream *flinesdata = NULL;
		ofstream *fxcorrdata = NULL;
		ofstream *fxcorrfitdata = NULL;
        
        if (!atlasdatafilename.empty()) {
            fatlasdata = new ofstream();
            fatlasdata->open(atlasdatafilename.c_str());
        }
        if (!compdatafilename.empty()) {
            fcompdata = new ofstream();
            fcompdata->open(compdatafilename.c_str());
        }
        if (!linesdatafilename.empty()) {
            flinesdata = new ofstream();
            flinesdata->open(linesdatafilename.c_str());
        }
        
        if (!xcorrdatafilename.empty()) {
            fxcorrdata = new ofstream();
            fxcorrdata->open(xcorrdatafilename.c_str());
        }
        
        if (!xcorrfitdatafilename.empty()) {
            fxcorrfitdata = new ofstream();
            fxcorrfitdata->open(xcorrfitdatafilename.c_str());
        }
        
		operaSpectralOrderVector spectralOrders(inputObjectSpectrum);
        
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
			cout << "operaTelluricWavelengthCorrection: minorder ="<< minorder << " maxorder=" << maxorder << endl;
    
		/*
		 * Read telluric reference files
		 */
        
		/*
		 * Read Telluric reference spectrum
		 *		lambda vs. intensity
		 */
		if (!telluric_spectrum.empty()) {
			if (debug) {
				cout << "operaTelluricWavelengthCorrection: telluric reference spectrum " << telluric_spectrum << endl;
			}
            nPointsInTelluricSpectrum = readTelluricSpectrum(telluric_spectrum, telluricSpectrumWavelength, telluricSpectrumIntensity);
        }
		/*
		 * Read telluric lines database
		 *		lambda vs. intensity
		 */
		if (!telluric_lines.empty()) {
			if (debug) {
				cout << "operaTelluricWavelengthCorrection: reading telluric lines database " << telluric_lines << endl;
			}
            ntelluriclines = readTelluricLines(telluric_lines,telluricMoleculeNumber,telluricLinesWavelength,telluricLinesIntensity);
        }

/*
        for(unsigned line=0;line<ntelluriclines;line++) {
            cout << line << " " << telluricMoleculeNumber[line] << " " << telluricLinesWavelength[line] << " " << telluricLinesIntensity[line] << endl;
        }
*/
        
		for (unsigned order=(unsigned)minorder; order<=(unsigned)maxorder; order++) {
            if(debug)
                cout << "\noperaTelluricWavelengthCorrection: Processing order = " << order << endl;
            
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if (spectralOrder->gethasWavelength()) {
                
                operaWavelength *wavelength =  spectralOrder->getWavelength();
                Polynomial *wavelengthPolynomial =  wavelength->getWavelengthPolynomial();
                
                if (spectralOrder->gethasSpectralElements()) {
					
					// DT May 20 2014 don't try to process order if there are not enough spectral elements
					if (spectralOrder->getSpectralElements()->getnSpectralElements() < MINELEMENTS)
						continue;
					
                    spectralOrder->applyNormalization(normalizationBinsize,0,FALSE,NULL,NULL,TRUE,0);
                    
                    if (debug) {
						cout << "operaTelluricWavelengthCorrection: reading object spectrum " << inputObjectSpectrum << endl;
					}
                    
                    operaSpectralElements *compSpectrum = spectralOrder->getSpectralElements();
                    compSpectrum->setwavelengthsFromCalibration(wavelength);
                    operaFluxVector *compfluxvector = compSpectrum->getFluxVector();

#ifdef PRINT_DEBUG
                    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                        cout << compSpectrum->getdistd(i) << " "<< compSpectrum->getwavelength(i) << " " << compfluxvector->getflux(i) << " " << compSpectrum->getXCorrelation(i) << endl;
                    }
#endif
                    
                    double wl_central = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()/2);
                    double wl0 = compSpectrum->getwavelength(0);
                    double wlf = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()-1);
                    
                    if(debug)
                        cout << "operaTelluricWavelengthCorrection: wavelength range: wl0=" << wl0 << ", wl_central=" << wl_central << ", wlf=" << wlf << endl;
                    
                    double linewidth = wl_central/spectralResolution;
                    
                    if(debug)
                        cout << "operaTelluricWavelengthCorrection: linewidth=" << linewidth << " nm" << endl;
 
                    operaSpectralLines *compLines = new operaSpectralLines(compSpectrum,linewidth,wavelength_disp);
                    
                    double CompLocalMaxFilterWidth = LocalMaxFilterWidth*linewidth;
					
                    double meanVariance = 0;
                    unsigned nvarpoints = 0;
                    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                        if(!isnan(compfluxvector->getvariance(i))) {
                            meanVariance += compfluxvector->getvariance(i);
                            nvarpoints++;
                        }
                    }
                    double CompMinPeakDepth = MinPeakDepth*sqrt(meanVariance/(double)nvarpoints);
					
                    compLines->detectAbsorptionSpectralFeatures(DetectionThreshold,CompLocalMaxFilterWidth,CompMinPeakDepth);
                    
                    if (compLines->getnLines() < 3) {	// changed from == 0 Nov 15 2012, since quicksort is called on data & it crashes with only 1 line, for example
                        if (verbose)
							printf("operaTelluricWavelengthCorrection: Warning: order %d: [Object Spectrum] %d lines detected from input object. Skipping calibration.\n", order, compLines->getnLines());
                        delete(compLines);
						continue;
                    } else {
                        if (debug) {
                            printf("operaTelluricWavelengthCorrection: order %d: [Object Spectrum] %d lines in object between wl0 = %.2f and wlf = %.2f.\n", order,compLines->getnLines(), wl0, wlf);
                        }
                    }                   
                     /*
                     * Below it reads wavelength and flux from Telluric reference spectrum
                     */
                    if (!telluric_spectrum.empty()) {
                        if (debug) {
                            cout << "operaTelluricWavelengthCorrection: calculating first order correction (wlshift) using cross-correlation.." << endl;
                        }

                        double wlshift = 0;
                        double maxcorr = 0;
                        
                        bool validXCorrelation = calculateWavelengthShiftByXCorr(compSpectrum, initialWavelengthRange, initialWavelengthStep, XCorrelationThreshold, halfslitsize, sigmaThreshold, &wlshift, &maxcorr, fxcorrdata, fxcorrfitdata, order);
                        
                        if (!validXCorrelation) {
                            wlshift = 0;
                        }
                        
                        if(debug)
                            cout << "operaTelluricWavelengthCorrection: order=" << order << " wlshift=" << wlshift <<  " maxcorr=" << maxcorr << endl;
                        
                        double zeroOrderCoeff = wavelengthPolynomial->getCoefficient(0);
                        wavelengthPolynomial->setCoefficient(0,zeroOrderCoeff+wlshift);
                        compSpectrum->setwavelengthsFromCalibration(wavelength);
                        
                        double wl_central = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()/2);
                        double wl0 = compSpectrum->getwavelength(0);
                        double wlf = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()-1);
                        
                        if(debug)
                            cout << "operaTelluricWavelengthCorrection: new wavelength range: wl0=" << wl0 << ", wl_central=" << wl_central << ", wlf=" << wlf << endl;
                        
                        double *wl,*transmission;
                        unsigned npointsInShortTelluricSpectrum = getTelluricSpectrumRange(wl0-fabs(wlf-wl0)/100,wlf+fabs(wlf-wl0)/100,&wl,&transmission);
                        
                        if (npointsInShortTelluricSpectrum > 0) {
                            if(verbose)
                                cout << "operaTelluricWavelengthCorrection: order=" << order << " validXCorrelation=" << validXCorrelation << " wlshift=" << wlshift  <<  " maxcorr=" << maxcorr << endl;
                            
                            if (debug) {
                                cout << "operaTelluricWavelengthCorrection: calculating second order correction based on telluric lines for order " << order << endl;
                                cout << "operaTelluricWavelengthCorrection: reading telluric reference spectrum " << telluric_spectrum << endl;
                            }
                            
                            /*
                             * Below it degrades the resolution of the atlas to the resolution of object spectrum.
                             * The degradation is done by convolving the spectrum with a gaussian.
                             */
                            if (debug) {
                                cout << "operaTelluricWavelengthCorrection: convolving spectrum with Gaussian.." << endl;
                            }
                            double *convolvedTelluricFlux = new double[npointsInShortTelluricSpectrum];
                            
                            convolveSpectrumWithGaussian(npointsInShortTelluricSpectrum,wl,transmission,convolvedTelluricFlux,linewidth);
                            
#ifdef PRINT_DEBUG
                            for(unsigned i=0;i<npointsInShortTelluricSpectrum;i++)
                                cout << "operaTelluricWavelengthCorrection: " << wl[i] << ' ' << transmission[i] << ' ' << convolvedTelluricFlux[i] << endl;
#endif
                            
                            /*
                             * Below it resamples telluric spectrum.
                             */
                            unsigned nTelluricElements = compSpectrum->getnSpectralElements();
                            double *telluricWavelength = new double[nTelluricElements];
                            double *telluricFlux = new double[nTelluricElements];
                            
                            for (unsigned i=0; i<nTelluricElements; i++) {
                                telluricWavelength[i] = compSpectrum->getwavelength(i);
                            }
                            operaFitSplineDouble(npointsInShortTelluricSpectrum,wl,convolvedTelluricFlux,nTelluricElements,telluricWavelength,telluricFlux);
                            
        
                            if (fatlasdata != NULL) { // for plotting
                                double maxatlasflux = -BIG;
                                for (unsigned i=0; i<nTelluricElements; i++) {
                                    if((1-telluricFlux[i]) > maxatlasflux)
                                        maxatlasflux = (1-telluricFlux[i]);
                                }
                                for (unsigned i=0; i<nTelluricElements; i++) {
                                    *fatlasdata << order << " " << telluricWavelength[i] << " " << (1 - (1-telluricFlux[i])/maxatlasflux) << " " << wl_central << endl;
                                }
                                *fatlasdata << endl;
                            }
                            
                            if (debug) {
                                cout << "operaTelluricWavelengthCorrection: npointsInShortTelluricSpectrum=" << npointsInShortTelluricSpectrum << endl;
                            }
                            
                            if (apply2ndOrderCorrection) {
                                
                                /*
                                 * Below it calculates the cross-correlation between the telluric spectrum and a gaussian function.
                                 */
                                if (debug) {
                                    cout << "operaTelluricWavelengthCorrection: calculating cross correlation with Gaussian.." << endl;
                                }
                                for (unsigned i=0; i<nTelluricElements; i++) {
                                    telluricFlux[i] = 1 - telluricFlux[i];
                                }
                                double *telluricSpectrumXCorr = new double[nTelluricElements];
                                
                                calculateXCorrWithGaussian(nTelluricElements,telluricWavelength,telluricFlux,telluricSpectrumXCorr,linewidth);
                                
                                /*
                                 * Below it reads the telluric spectrum into an operaSpectralElements class
                                 */
                                if (debug) {
                                    cout << "operaTelluricWavelengthCorrection: reading telluric spectrum into operaSpectralElement.." << endl;
                                }
                                operaSpectralElements *telluricSpectrum = new operaSpectralElements(nTelluricElements);
                                
                                for (unsigned i=0; i<nTelluricElements; i++) {
                                    telluricSpectrum->setXCorrelation(telluricSpectrumXCorr[i], i);
                                    telluricSpectrum->setwavelength(telluricWavelength[i], i);
                                    telluricSpectrum->setFlux(telluricFlux[i],i);
                                    telluricSpectrum->setFluxVariance(0.0001,i);
                                }
                                
                                telluricSpectrum->setHasXCorrelation(true);
                                telluricSpectrum->setHasWavelength(true);
                                telluricSpectrum->setHasRawSpectrum(true);
                                
                                if(debug) {
                                    for (unsigned i=0; i<nTelluricElements; i++) {
                                        cout << telluricSpectrum->getwavelength(i)  << " " << telluricSpectrum->getFlux(i)  << " " << telluricSpectrum->getFluxVariance(i) << " " << telluricSpectrum->getXCorrelation(i)<< " " << compfluxvector->getflux(i) << " " << compSpectrum->getXCorrelation(i) << endl;
                                    }
                                }
                                
                                /*
                                 * Below it creates an operaSpectralLines class for the telluric lines
                                 */
                                if (debug) {
                                    cout << "operaTelluricWavelengthCorrection: creating an operaSpectralLines class for the telluric lines.." << endl;
                                }
                                operaSpectralLines *telluricLines = new operaSpectralLines(telluricSpectrum, linewidth, wavelength_disp);
                                
                                /*
                                 * Below it sets detection thresholds and run the algorithm to detect spectral lines in the reference telluric spectrum
                                 */
                                if (debug) {
                                    cout << "operaTelluricWavelengthCorrection: setting detection thresholds and running the algorithm to detect spectral lines in the reference telluric spectrum.." << endl;
                                }
                                double TelluricLocalMaxFilterWidth = LocalMaxFilterWidth*linewidth;
                                double TelluricMinPeakDepth = 0.000001;
                                telluricLines->detectSpectralFeatures(DetectionThreshold,TelluricLocalMaxFilterWidth,TelluricMinPeakDepth);
                                
                                /*
                                 * Below it reads the telluric lines information
                                 */
                                if (debug) {
                                    cout << "operaTelluricWavelengthCorrection: reading the telluric lines information.." << endl;
                                }
                                
                                if (telluricLines->getnLines() < 3) {	// changed from == 0 Nov 15 2012, since quicksort is called on data & it crashes with only 1 line, for example
									if (verbose)
										printf("operaTelluricWavelengthCorrection: Warning:  order %d: [Telluric] %d lines detected from input telluric reference. Skipping calibration.\n", order, telluricLines->getnLines());
                                    delete(telluricLines);
                                    continue;
                                } else {
                                    if (debug) {
                                        printf("operaTelluricWavelengthCorrection: order %d: [Telluric] %d lines detected in telluric reference within the range wl0=%.2f, wlf=%.2f.\n", order,  telluricLines->getnLines(), wl0, wlf);
                                    }
                                }
                                
                                /*
                                 *
                                 * Here it should be checked the quality of telluric reference.
                                 * Say we need a minimum number of lines, with good coverage and with good signal.
                                 * If any of the conditions above fail it must not be a good source for
                                 * calibration.
                                 *
                                 *
                                 * Mininum number of lines is arbitrary but say a minimum of 50
                                 *
                                 */
                                
                                unsigned Coeffs = 2;
                                
                                Polynomial *wlcorrection = new Polynomial(Coeffs);
                                
                                unsigned numberOfMatchedLines = matchTelluricReferencewithObjectLines(1.0,linewidth,telluricLines,compLines,wlcorrection,order,wl_central,flinesdata);
                                
                                double *par = (double *)wlcorrection->getVector();
                                double *parerr = (double *)wlcorrection->getErrorVector();
                                
                                double zeroOrderModelVariance = wavelengthPolynomial->getCoefficientError(0)*wavelengthPolynomial->getCoefficientError(0);
                                
                                if (verbose) {
                                    cout << "operaTelluricWavelengthCorrection: order=" << order << ", numberOfMatchedLines=" << numberOfMatchedLines << ", zeroOrderModelVariance=" << zeroOrderModelVariance << ", chisqr=" << wlcorrection->getChisqr() << endl;
                                }
                                
                                /*
                                 * Apply wavelength correction up to 2nd order
                                 *
                                 *   From previous calibration we have: wl = a + b*x + c*x*x + d*x*x*x
                                 *   The new telluric correction gives: new_wl = u + v*wl
                                 *
                                 *   Therefore: new_wl = (u + v*a) + v*b*x + v*c*x*x + v*d*x*x*x
                                 *
                                 *   The correction will only be applied if chisqr < variance(wl)
                                 */
                                
                                if(wlcorrection->getChisqr() < zeroOrderModelVariance || zeroOrderModelVariance == 0) {
                                    unsigned outputNpar = wavelengthPolynomial->getOrderOfPolynomial();
                                    
                                    double *newpar = new double[outputNpar];
                                    double *newparerr = new double[outputNpar];
                                    
                                    newpar[0] = par[0] + par[1]*wavelengthPolynomial->getCoefficient(0);
                                    newparerr[0] = sqrt(parerr[0]*parerr[0] + (wavelengthPolynomial->getCoefficient(0)*parerr[1] + wavelengthPolynomial->getCoefficientError(0)*par[1])*(wavelengthPolynomial->getCoefficient(0)*parerr[1] + wavelengthPolynomial->getCoefficientError(0)*par[1]));
                                    
                                    for(unsigned coeff=1;coeff<outputNpar;coeff++) {
                                        newpar[coeff] = par[1]*wavelengthPolynomial->getCoefficient(coeff);
                                        newparerr[coeff] = sqrt((wavelengthPolynomial->getCoefficient(coeff)*parerr[1] + wavelengthPolynomial->getCoefficientError(coeff)*par[1])*(wavelengthPolynomial->getCoefficient(coeff)*parerr[1] + wavelengthPolynomial->getCoefficientError(coeff)*par[1]));
                                    }
                                    
                                    // cout << wlcorrection->getChisqr() << " " << zeroOrderModelVariance << " " << zeroOrderCorrectedVariance << endl;
                                    for(unsigned coeff=0;coeff<outputNpar;coeff++) {
                                        wavelengthPolynomial->setCoefficient(coeff,newpar[coeff]);
                                        wavelengthPolynomial->setCoefficientError(coeff, newparerr[coeff]);
                                    }
                                    
                                    delete[] newpar;
                                    delete[] newparerr;
                                } else {
                                    if (verbose) {
                                        cout << "operaTelluricWavelengthCorrection: no wavelength correction has been applied using telluric lines." << endl;
                                    }
                                }
                                

                                delete[] telluricSpectrumXCorr;
                                delete(telluricLines);
                                delete(wlcorrection);
                            }
                            
                            delete[] convolvedTelluricFlux;
                            delete[] telluricWavelength;
                            delete[] telluricFlux;
                            
                        } else {
                            if(verbose)
                                cout << "operaTelluricWavelengthCorrection: order=" << order << " no points in the telluric spectrum." << endl;
                        }
                        
                        if (fcompdata != NULL){ // for plotting
                            for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                                double dist = compSpectrum->getdistd(i);
                                double lambda = compSpectrum->getwavelength(i);
                                double flux = compSpectrum->getFlux(i);
                                *fcompdata << order << " " << dist << " " << lambda << " " << flux << " " << wl_central << endl;
                            }
                            *fcompdata << endl;
                        }
                        
                    } else {
                        if(verbose)
                            cout << "operaTelluricWavelengthCorrection: order=" << order << " empty input telluric spectrum. Skipping calibration." << endl;
                        continue;
                    } // if (!telluric_spectrum.empty())
					//delete(compLines);
                } else {
                    if(verbose)
                        cout << "operaTelluricWavelengthCorrection: order=" << order << " no spectral elements in input spectrum. Skipping calibration." << endl;
                    continue;
                } // if (spectralOrder->gethasSpectralElements())
            
                if(debug) {
                    cout << "operaTelluricWavelengthCorrection: corrected wavelength solution for order " << order << endl;
                    cout << "operaTelluricWavelengthCorrection: ";
                    wavelengthPolynomial->printEquation(&cout);
                    cout << endl;
                }
            } else {
                if(verbose)
                    cout << "operaTelluricWavelengthCorrection: order=" << order << " no wavelength calibration for order " << order << endl;
                continue;
            } // if (spectralOrder->gethasWavelength())
 		}

		// output a new wavelength calibration file
		spectralOrders.WriteSpectralOrders(outputWaveFile, Wave);
        
        /*
         * Telluric wavelength correction orders info plot: 
         */
        if (fxcorrdata != NULL && fxcorrfitdata != NULL) {
            fxcorrdata->close();
            fxcorrfitdata->close();
            if (!xcorrscriptfilename.empty()) {
                GenerateTelluricXCorrelationPlot(xcorrscriptfilename, xcorrsplotfilename, xcorrdatafilename, xcorrfitdatafilename, interactive);
            }
        }
        
        /*
         * Telluric wavelength correction plot: plot atlas and comparison spectra and final set of matched lines.
         */

        if (fatlasdata != NULL && fcompdata != NULL && flinesdata != NULL) {
            fatlasdata->close();
            fcompdata->close();
            flinesdata->close();
            
            if (!specscriptfilename.empty()) {
                GenerateTelluricSpecPlot(specscriptfilename, specplotfilename, atlasdatafilename, compdatafilename, linesdatafilename, subtractCentralWavelength, interactive, apply2ndOrderCorrection);
            }
        }
        
	}
	catch (operaException e) {
		cerr << "operaTelluricWavelengthCorrection: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaTelluricWavelengthCorrection: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/*
 * Generate multiple plot containing statistical info about telluric wavelength correction
 */
void GenerateTelluricXCorrelationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, string cleanDataFileName, bool display)
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

    *fgnu << "\nset xlabel \"{/Symbol Dl} (nm)\"" << endl;
    *fgnu << "set ylabel \"order number + cross-correlation\"" << endl;
    *fgnu << "set pointsize 1.5" << endl;

    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        *fgnu << "plot \"" << dataFileName << "\" u 3:($1+$5) t \"gaussian fit\" w l lt 3 lw 2, ";
        *fgnu << "\"" << dataFileName << "\" u 3:($1+$6) t \"XCorr data\" w p pt 6, ";
        *fgnu << "\"" << dataFileName << "\" u 3:($1+$7) t \"SigmaThreshold\" w l, ";
        *fgnu << "\"" << dataFileName << "\" u 3:($1+$8) t \"XCorrThreshold\" w l, ";
        *fgnu << "\"" << cleanDataFileName << "\" u 3:($1+$6) t \"fit data\" w p pt 7" << endl;
        
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        *fgnu << endl;
        
        *fgnu << "plot \"" << dataFileName << "\" u 3:($1+$5) t \"gaussian fit\" w l lt 3 lw 2, ";
        *fgnu << "\"" << dataFileName << "\" u 3:($1+$6) t \"XCorr data\" w p pt 6, ";
        *fgnu << "\"" << dataFileName << "\" u 3:($1+$7) t \"SigmaThreshold\" w l, ";
        *fgnu << "\"" << dataFileName << "\" u 3:($1+$8) t \"XCorrThreshold\" w l, ";
        *fgnu << "\"" << cleanDataFileName << "\" u 3:($1+$6) t \"fit data\" w p pt 7" << endl;
        
        *fgnu << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }
}


/*
 * Generate 2D plot for spectra of atlas + comparison + identified lines
 */
void GenerateTelluricSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, string linesdatafilename, bool subtractCentralWavelength, bool display, bool apply2ndOrderCorrection) {
    
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    if(subtractCentralWavelength) {
        *fgnu << "\nset xlabel \"{/Symbol l} - {/Symbol l}_c (nm)\"" << endl;
    } else {
        *fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
    }
    *fgnu << "set ylabel \"order number + norm flux\"" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        if(subtractCentralWavelength) {
            *fgnu << "plot \"" << atlasdatafilename << "\" u ($2-$4):($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$5):($1 + $4) w l lt 3";
            if(apply2ndOrderCorrection) {
                *fgnu << ",\"" << linesdatafilename << "\" u ($3-$6):($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
            } else {
                *fgnu << endl;
            }
        } else {
            *fgnu << "plot \"" << atlasdatafilename << "\" u 2:($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u 3:($1 + $4) w l lt 3";
            if(apply2ndOrderCorrection) {
                *fgnu << ",\"" << linesdatafilename << "\" u 3:($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
            } else {
                *fgnu << endl;
            }
        }
        
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        *fgnu << endl;
        
        if(subtractCentralWavelength) {
            *fgnu << "plot \"" << atlasdatafilename << "\" u ($2-$4):($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$5):($1 + $4) w l lt 3";
            if(apply2ndOrderCorrection) {
                *fgnu << ",\"" << linesdatafilename << "\" u ($3-$6):($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
            } else {
                *fgnu << endl;
            }
        } else {
            *fgnu << "plot \"" << atlasdatafilename << "\" u 2:($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u 3:($1 + $4) w l lt 3";
            if(apply2ndOrderCorrection) {
                *fgnu << ",\"" << linesdatafilename << "\" u 3:($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
            } else {
                *fgnu << endl;
            }
        }
        *fgnu << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }    
}


/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --inputWaveFile=<WAVE_FILE>"
	" --inputObjectSpectrum=<SPECTRUM_FILE>"
	" --outputWaveFile=<WAVE_FILE>"
	" --telluric_lines=<LINES_FILE>"
	" --telluric_spectrum=<SPECTRUM_FILE>"
	" --spectralResolution=<DBL_VALUE>"
	" --initialWavelengthRange=<DBL_VALUE>"
	" --initialWavelengthStep=<DBL_VALUE>"
	" --XCorrelationThreshold=<DBL_VALUE>"
	" --sigmaThreshold=<FLT_VALUE>"
	" --normalizationBinsize=<UNS_VALUE>"
    " --ordernumber=<INT_VALUE"
	" --minorder=<INT_VALUE>"
	" --maxorder=<INT_VALUE>"
    " --xcorrsplotfilename=<EPS_FILE>"
    " --specplotfilename=<EPS_FILE>"
    " --xcorrscriptfilename=<GNUPLOT_FILE>"
    " --specscriptfilename=<GNUPLOT_FILE>"
    " --xcorrdatafilename=<DATA_FILE>"
    " --xcorrfitdatafilename=<DATA_FILE>"
    " --atlasdatafilename=<DATA_FILE>"
    " --compdatafilename=<DATA_FILE>"
    " --linesdatafilename=<DATA_FILE>"
    " --subtractCentralWavelength=<BOOL>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputObjectSpectrum=/Users/edermartioli//opera//spectra/51Peg-12BQ10-Dec01/1599045i.e.gz --inputWaveFile=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/OLAPAa_pol_Normal.wcal.gz --telluric_lines=/Users/edermartioli/opera-1.0//config/opera_HITRAN08-extracted.par.gz --telluric_spectrum=/Users/edermartioli/opera-1.0//config/KPNO_atmtrans.dat.gz --spectralResolution=80000 --initialWavelengthRange=0.1 --initialWavelengthStep=0.002 --XCorrelationThreshold=0.1 --sigmaThreshold=1.25 --normalizationBinsize=110 --xcorrdatafilename=1599045itellxcorr.pdat --xcorrfitdatafilename=1599045itellxcorrfit.pdat --atlasdatafilename=1599045itellatlas.pdat --compdatafilename=1599045itellcomp.pdat --linesdatafilename=1599045itelllines.pdat --specscriptfilename=1599045itell.gnu --xcorrscriptfilename=1599045itellxcorr.gnu --subtractCentralWavelength=1 --outputWaveFile=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/1599045i.tell.gz -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file \n"
	"  -s, --inputObjectSpectrum=<SPECTRUM_FILE>, Input object spectrum file (.s)\n"
	"  -o, --outputWaveFile=<WAVE_FILE>, Output wavelength calibration file to store final solution\n"
	"  -L, --telluric_lines=<LINES_FILE>, Atlas of telluric lines \n"
	"  -T, --telluric_spectrum=<SPECTRUM_FILE>, Spectrum of telluric lines\n"
	"  -R, --spectralResolution=<DBL_VALUE>, Input spectral resolution (wl/dwl) as reference for line detection\n"
	"  -r, --initialWavelengthRange=<DBL_VALUE>, Wavelength shift range (in nm) to scan for first order correction\n"
	"  -s, --initialWavelengthStep=<DBL_VALUE>, Wavelength step (in nm) to scan for first orer correction\n"
	"  -x, --XCorrelationThreshold=<DBL_VALUE>, X-correlation lower threshold to consider a match between telluric and object spectra\n"
	"  -g, --sigmaThreshold=<FLT_VALUE>, X-correlation lower threshold based on the dispersion of xcorrelation data\n"
	"  -b, --normalizationBinsize=<UNS_VALUE>, Binsize to normalize input object spectrum \n"
    "  -O, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
	"  -N, --minorder=<INT_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<INT_VALUE>, Define maximum order number\n"
	"  -P, --xcorrsplotfilename=<EPS_FILE>\n"
	"  -Q, --specplotfilename=<EPS_FILE>\n"
	"  -D, --xcorrdatafilename=<DATA_FILE>\n"
	"  -F, --xcorrfitdatafilename=<DATA_FILE>\n"
	"  -G, --atlasdatafilename=<DATA_FILE>\n"
	"  -H, --compdatafilename=<DATA_FILE>\n"
	"  -J, --linesdatafilename=<DATA_FILE>\n"
	"  -S, --xcorrscriptfilename=<GNUPLOT_FILE>\n"
	"  -U, --specscriptfilename=<GNUPLOT_FILE>\n"
	"  -K, --subtractCentralWavelength=<BOOL>,Choose to subtract order central wavelength for plot\n"
	"  -I, --interactive=<BOOL>\n\n";
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
//        cout << "operaTelluricWavelengthCorrection: " << telluricMoleculeNumber[i] << " " << telluricLinesWavelength[i] << " " << telluricLinesIntensity[i] << endl;
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


void generateSyntheticTelluricSpectrumUsingGaussianProfile(unsigned np, double *wavelengthVector, double *ouputSpectrum, double resolution) {
    
    for(unsigned i=0;i<np;i++) {    
        ouputSpectrum[i] = 1.0;
    }
    
    double *wl, *intensity;
    
    for(unsigned i=0;i<np;i++) {
        double gaussianWidth = wavelengthVector[i]/resolution;

        unsigned nlinesInRange = getTelluricLinesRange(wavelengthVector[i] - 5*gaussianWidth,wavelengthVector[i] + 5*gaussianWidth,&wl,&intensity);
        
        if(debug) {
            cout << "operaTelluricWavelengthCorrection: " << i <<
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

/*
 * Read the the full atmospheric transmission spectrum
 */
unsigned readTelluricSpectrum(string telluric_spectrum, double *telluricSpectrumWavelength, double *telluricSpectrumIntensity) {
	igzstream astream;
	string dataline;
    
	double tmpwl = -1.0;
	double tmpi = -1.0;
	unsigned np = 0;
	
	astream.open(telluric_spectrum.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi);
                    
                    telluricSpectrumWavelength[np] = tmpwl;
                    telluricSpectrumIntensity[np] = tmpi;
                    np++;
                }	// skip comments
            }
		} // while (astream.good())
		
		if (np > 0) {
			if (verbose) {
				printf("          [Telluric] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, telluricSpectrumWavelength[0], telluricSpectrumWavelength[np/2], telluricSpectrumWavelength[np-1]);
			}
		} else {
			printf("          [Telluric] no points found in telluric spectrum.\n");
		}
		astream.close();
	}	// if (astream.open()
	return np;
}

/*
 * get a subset of the telluric lines for this order only, between wl0 and wlf
 */
unsigned getTelluricSpectrumRange(double wl0, double wlf, double **wl, double **intensity) {
    
    unsigned firstline = 0;
	unsigned line = 0;
	
	for (line=0; line<nPointsInTelluricSpectrum; line++) {
		if (telluricSpectrumWavelength[line] >= wl0) {
			if (firstline == 0) {
				*intensity = &telluricSpectrumIntensity[line];
				*wl = &telluricSpectrumWavelength[line];
				firstline = line;
			}
			if (telluricSpectrumWavelength[line] > wlf)
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
		throw operaException("operaTelluricWavelengthCorrection: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (telluricLines->getnLines() == 0) {
		throw operaException("operaTelluricWavelengthCorrection: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
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
            cout << "operaTelluricWavelengthCorrection: " << "mindiff=" << mindifference << " diff=" << difference << " object[" << i << "]=" << objectlinecenter[i] << " telluric[" << l << "]=" << telluricLineswl[l] << endl;
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

bool calculateWavelengthShiftByXCorr(operaSpectralElements *compSpectrum, double DWavelengthRange, double DWavelengthStep, double threshold, int halfslitsize, float nsigcut, double *maxDWavelength, double *maxcorr, ostream *fxcorrdata, ostream *fxcorrfitdata, unsigned order) {
    
    bool status = true;
    
    operaFluxVector *compfluxvector = compSpectrum->getFluxVector();
    
    double *objectSpectrum = new double[compSpectrum->getnSpectralElements()];
    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
        objectSpectrum[i] = compfluxvector->getflux(i);
    }
    
    double *telluricWavelength = new double[compSpectrum->getnSpectralElements()];
    double *telluricSpectrum = new double[compSpectrum->getnSpectralElements()];
//    double *hitranTelluricSpectrum = new double[compSpectrum->getnSpectralElements()];
    
    double firstDWavelength = -DWavelengthRange/2;
    //double lastDWavelength = +DWavelengthRange/2;
    
    double wl0 = compSpectrum->getwavelength(0)*(1 - DWavelengthRange);
    double wlf = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()-1)*(1 + DWavelengthRange);
    double wl_central = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()/2);

    double *wl,*transmission_p;
    unsigned npointsInShortTelluricSpectrum = getTelluricSpectrumRange(wl0,wlf,&wl,&transmission_p);

    double *transmission = new double [npointsInShortTelluricSpectrum];

    /*
     * Normalize telluric reference
     */
    double maxatlasflux = -BIG;
    for (unsigned i=0; i<npointsInShortTelluricSpectrum; i++) {
        if((1-transmission_p[i]) > maxatlasflux)
            maxatlasflux = (1-transmission_p[i]);
    }
    for (unsigned i=0; i<npointsInShortTelluricSpectrum; i++) {
       transmission[i]= 1 - (1-transmission_p[i])/maxatlasflux;
    }
    
    unsigned nDataPoints = (unsigned)ceil(DWavelengthRange/DWavelengthStep);
    
    double DWavelength = firstDWavelength;
    unsigned jmax = 0;
    *maxcorr = -BIG;
    *maxDWavelength = 0;
	
	if (npointsInShortTelluricSpectrum == 0) {
		*maxcorr = NAN;
		*maxDWavelength = NAN;
		return false;
	}
    
    float *crosscorrelation = new float [nDataPoints];
    float *dwl = new float [nDataPoints];
    float *dRV = new float [nDataPoints];

    for(unsigned j=0; j<nDataPoints;j++) {
        for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
            telluricWavelength[i] = compSpectrum->getwavelength(i) + DWavelength;
        }
        
        operaFitSplineDouble(npointsInShortTelluricSpectrum,wl,transmission,compSpectrum->getnSpectralElements(),telluricWavelength,telluricSpectrum);
        
        crosscorrelation[j] = (float)(operaCrossCorrelation(compSpectrum->getnSpectralElements(),objectSpectrum,telluricSpectrum));
        dwl[j] = (float)DWavelength;
        dRV[j] = (float)(DWavelength/wl_central)*SPEED_OF_LIGHT_M*(1e-3);
        
        //     generateSyntheticTelluricSpectrumUsingGaussianProfile(compSpectrum->getnSpectralElements(),telluricWavelength,hitranTelluricSpectrum,80000);
        //     double crosscorrelation = operaCrossCorrelation(compSpectrum->getnSpectralElements(),objectSpectrum,hitranTelluricSpectrum);
        
        if((double)(crosscorrelation[j]) > *maxcorr && (double)(crosscorrelation[j]) > threshold) {
            *maxcorr = (double)(crosscorrelation[j]);
            *maxDWavelength = DWavelength;
            jmax = j;
        }
        
        DWavelength+=DWavelengthStep;
    }
    if (*maxDWavelength == firstDWavelength) {
		*maxcorr = NAN;
		*maxDWavelength = NAN;
        status = false;
	}
    if(debug)
        cout << "calculateWavelengthShiftByXCorr: status=" << status << " order=" << order << " wlshift=" << *maxDWavelength  <<  " maxcorr=" << *maxcorr << endl;

    /*
     * Two tests will be performed:
     *  I) Whether it is an isolated maximum point, i.e. three points before and three points after all must be lower than max
     *
     *  II) Whether it's a significant maximum, i.e. maxcorrelation > median + 3*sigma
     */
    
    halfslitsize = 5;
    
    double *peakXdata = new double[2*halfslitsize+2];
    double *peakYdata = new double[2*halfslitsize+2];
    unsigned np = 0;
    // Test (I)
    
    if((int)jmax - halfslitsize < 0 || jmax + halfslitsize >= nDataPoints) {
        *maxcorr = NAN;
		*maxDWavelength = NAN;
        status = false;
    } else {
        float goingup = -BIG;
        for(unsigned j=jmax-halfslitsize;j<jmax;j++) {
            peakXdata[np] = (double)dwl[j];
            peakYdata[np] = (double)crosscorrelation[j];
            np++;
            if(crosscorrelation[j] > goingup) {
                goingup = crosscorrelation[j];
            } else {
                *maxcorr = NAN;
                *maxDWavelength = NAN;
                np=0;
                status = false;
                break;
            }
        }
  
        peakXdata[np] = (double)dwl[jmax];
        peakYdata[np] = (double)crosscorrelation[jmax];
        np++;
        
        float goingdown = BIG;
        for(unsigned j=jmax+1;j<=jmax+halfslitsize;j++) {
            peakXdata[np] = (double)dwl[j];
            peakYdata[np] = (double)crosscorrelation[j];
            np++;
            if(crosscorrelation[j] < goingdown) {
                goingdown = crosscorrelation[j];
            } else {
                *maxcorr = NAN;
                *maxDWavelength = NAN;
                np=0;
                status = false;
                break;
            }
        }
    }
    if(debug)
        cout << "calculateWavelengthShiftByXCorr: Test I status=" << status << endl;

    // Test (II)
    
    float medianXcorr = operaArrayMedian(nDataPoints,crosscorrelation);
    float medsigXcorr = operaArrayMedianSigma(nDataPoints,crosscorrelation,medianXcorr);
    
    if(debug)
        cout << "calculateWavelengthShiftByXCorr:maxcorr=" << crosscorrelation[jmax] << "  medianXcorr=" << medianXcorr << "  medsigXcorr=" << medsigXcorr << "  median+n*sig=" << medianXcorr + nsigcut*medsigXcorr << endl;

    if(crosscorrelation[jmax] < medianXcorr + nsigcut*medsigXcorr) {
        *maxcorr = NAN;
        *maxDWavelength = NAN;
        status = false;
    }
    
    if(debug)
        cout << "calculateWavelengthShiftByXCorr: Test II status=" << status << endl;

    if(status==true) {
        double a=(double)crosscorrelation[jmax];
        double x0=(double)dwl[jmax];
        double sig=(double)(fabs(dwl[jmax+halfslitsize/2]-dwl[jmax-halfslitsize/2]))/2.0;
        double chisqr;
        
        operaLMFitGaussian(np, peakXdata, peakYdata, &a, &x0, &sig, &chisqr);
        *maxcorr = a;
        *maxDWavelength = x0;
        
        // Below is for plotting
        if(fxcorrfitdata != NULL) {
            for(unsigned j=0; j<np;j++) {
                double x = (double)peakXdata[j];
                double rv = (double)(peakXdata[j]/wl_central)*SPEED_OF_LIGHT_M*(1e-3);
                double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
                *fxcorrfitdata <<  order  << " " << wl_central << " " << x << " " << rv << " " <<  gaussfunc << " " <<  peakYdata[j] << " " << medianXcorr + nsigcut*medsigXcorr <<  " " << threshold << endl;
            }
            *fxcorrfitdata << endl;
        }
        
        // Below is for plotting
        if(fxcorrdata != NULL) {
            for(unsigned j=0; j<nDataPoints;j++) {
                double x = (double)dwl[j];
                double gaussfunc = a*exp(-(x-x0)*(x-x0)/(2*sig*sig));
                
                *fxcorrdata <<  order  << " " << wl_central << " " << dwl[j] << " " << dRV[j] << " " <<  gaussfunc << " " <<  crosscorrelation[j] << " " << medianXcorr + nsigcut*medsigXcorr << " " << threshold <<  endl;
            }
            *fxcorrdata << endl;
        }
    } else {
        // Below is for plotting
        if(fxcorrdata != NULL) {
            for(unsigned j=0; j<nDataPoints;j++) {
                *fxcorrdata <<  order  << " " << wl_central << " " << dwl[j] << " " << dRV[j] << " " <<  0.0 << " " <<  crosscorrelation[j] << " " << medianXcorr + nsigcut*medsigXcorr << " " << threshold << endl;
            }
            *fxcorrdata << endl;
        }
    }

    delete[] peakXdata;
    delete[] peakYdata;
    
    delete[] crosscorrelation;
    delete[] dwl;
    delete[] dRV;
    
    delete[] objectSpectrum;
    delete[] telluricWavelength;
    delete[] telluricSpectrum;
    //    delete[] hitranTelluricSpectrum;

    delete[] transmission;
        
	return status;
}



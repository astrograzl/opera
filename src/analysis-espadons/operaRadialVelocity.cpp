/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaRadialVelocity
 Version: 1.0
 Description: Find Radial Velocity
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

/*! 
 * operaRadialVelocity
 * \author Eder Martioli
 * \brief Find Radial Velocity.
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
     * Parameters for heliocentric wavelength correction
     */
	
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;
    bool maxorderprovided = false;
    
    bool interactive = false;
    
	//bool debug=false, verbose=false, trace=false;
    
    bool plot=false;
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
    
    /*
     * The parameters below we don't know yet whether they will be input
     */
    double DetectionThreshold = 0.05;    // threshold to regulate the sensitivity of line detection. Must be between 0 and 1.
    double LocalMaxFilterWidth = 3.0;    // parameter to set a window filter to guarantee a line is not detected twice. It's in units of line width
    double MinPeakDepth = 0.3;           // limit that also regulates the sensitity of line detection in units of noise.
    
    double spectralResolution = 80000; // Input spectral resolution as reference for line detection
    
    double initialWavelengthRange = 0.1;
    double initialWavelengthStep = 0.0001;
    double XCorrelationThreshold = 0.05;
    
	struct option longopts[] = {
		{"inputWaveFile",1, NULL, 'w'},				// input wavelength calibration file (.wcal)
 		{"inputObjectSpectrum",1, NULL, 'i'},		// input object spectrum file (.e or .p)
        {"outputWaveFile",1, NULL, 'o'},			// output wavelength calibration file (.tell or .ptell)
		{"telluric_lines",1, NULL, 'L'},            // atlas of telluric lines
		{"telluric_spectrum",1, NULL, 'T'},         // spectrum of telluric lines
		{"spectralResolution",1, NULL, 'R'},        // Input spectral resolution (wl/dwl) as reference for line detection
		{"initialWavelengthRange",1, NULL, 'r'},    // Wavelength shift range (in nm) to scan for first order correction
		{"initialWavelengthStep",1, NULL, 's'},     // Wavelength step (in nm) to scan for first orer correction
		{"XCorrelationThreshold",1, NULL, 'x'},     // X-correlation lower threshold to consider a match between telluric and object spectra
		{"ordernumber",			1, NULL, 'O'},
		{"minorder",			1, NULL, 'M'},
		{"maxorder",			1, NULL, 'X'},
        {"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},
		{"interactive",0, NULL, 'I'},
		
		{"plot",0, NULL, 'p'},
		{"verbose",0, NULL, 'v'},
		{"debug",0, NULL, 'd'},
		{"trace",0, NULL, 't'},
		{"help",0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:w:s:o:O:L:T:R:r:s:x:M:X:P:F:S:I:p::v::d::t::h",
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
				plotfilename = optarg;
				plot = 1;
				break;
			case 'F':
				datafilename = optarg;
				break;
			case 'S':
				scriptfilename = optarg;
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
			throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input object spectrum file ...        
		if (inputObjectSpectrum.empty()) {
			throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output wavelength calibration file ...
		if (outputWaveFile.empty()) {
			throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input atlas of telluric lines ...
		if (telluric_lines.empty() && telluric_spectrum.empty()) {
			throw operaException("operaRadialVelocity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (verbose) {
			cout << "operaRadialVelocity: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaRadialVelocity: inputObjectSpectrum = " << inputObjectSpectrum << endl;
			cout << "operaRadialVelocity: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaRadialVelocity: telluric_lines =" << telluric_lines << endl;
			cout << "operaRadialVelocity: telluric_spectrum =" << telluric_spectrum << endl;
			cout << "operaRadialVelocity: spectralResolution =" << spectralResolution << endl;
			cout << "operaRadialVelocity: initialWavelengthRange =" << initialWavelengthRange << endl;
			cout << "operaRadialVelocity: initialWavelengthStep =" << initialWavelengthStep << endl;
			cout << "operaRadialVelocity: XCorrelationThreshold =" << XCorrelationThreshold << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaRadialVelocity: ordernumber = " << ordernumber << endl;
            }
            if(plot) {
                cout << "operaRadialVelocity: plotfilename = " << plotfilename << endl;
                cout << "operaRadialVelocity: datafilename = " << datafilename << endl;
                cout << "operaRadialVelocity: scriptfilename = " << scriptfilename << endl;
                if(interactive) {
                    cout << "operaRadialVelocity: interactive = YES" << endl;
                } else {
                    cout << "operaRadialVelocity: interactive = NO" << endl;
                }
            }
            
		}
        
        ofstream *fdata = NULL;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());
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
			cout << "operaRadialVelocity: minorder ="<< minorder << " maxorder=" << maxorder << endl;
		
		/*
		 * Read telluric reference files
		 */
        
		/*
		 * Read Telluric reference spectrum
		 *		lambda vs. intensity
		 */
		if (!telluric_spectrum.empty()) {
			if (verbose) {
				cout << "operaRadialVelocity: telluric reference spectrum " << telluric_spectrum << endl;
			}
            nPointsInTelluricSpectrum = readTelluricSpectrum(telluric_spectrum, telluricSpectrumWavelength, telluricSpectrumIntensity);
        }
		/*
		 * Read telluric lines database
		 *		lambda vs. intensity
		 */
		if (!telluric_lines.empty()) {
			if (verbose) {
				cout << "operaRadialVelocity: reading telluric lines database " << telluric_lines << endl;
			}
            ntelluriclines = readTelluricLines(telluric_lines,telluricMoleculeNumber,telluricLinesWavelength,telluricLinesIntensity);
        }
		
		/*
		 for(unsigned line=0;line<ntelluriclines;line++) {
		 cout << line << " " << telluricMoleculeNumber[line] << " " << telluricLinesWavelength[line] << " " << telluricLinesIntensity[line] << endl;
		 }
		 */
        
		for (unsigned order=(unsigned)minorder; order<=(unsigned)maxorder; order++) {
            if(verbose)
                cout << "\noperaRadialVelocity: *** Processing order = " << order << endl;
            
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if (spectralOrder->gethasWavelength()) {
                
                operaWavelength *wavelength =  spectralOrder->getWavelength();
                Polynomial *wavelengthPolynomial =  wavelength->getWavelengthPolynomial();
                
                if (spectralOrder->gethasSpectralElements()) {
					
                    spectralOrder->applyNormalization(90,0,FALSE,NULL,NULL,TRUE,0);
                    
                    if (verbose) {
						cout << "operaRadialVelocity: reading object spectrum " << inputObjectSpectrum << endl;
					}
                    
                    operaSpectralElements *compSpectrum = spectralOrder->getSpectralElements();
                    compSpectrum->setwavelengthsFromCalibration(wavelength);
                    operaFluxVector *compfluxvector = compSpectrum->getFluxVector();
                    /*                    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                     cout << compSpectrum->getdistd(i) << " "<< compSpectrum->getwavelength(i) << " " << compfluxvector->getflux(i) << " " << compSpectrum->getXCorrelation(i) << endl;
                     }
                     */
                    double wl_central = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()/2);
                    double wl0 = compSpectrum->getwavelength(0);
                    double wlf = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()-1);
                    
                    if(verbose)
                        cout << "operaRadialVelocity: wavelength range: wl0=" << wl0 << ", wl_central=" << wl_central << ", wlf=" << wlf << endl;
                    
                    double linewidth = wl_central/spectralResolution;
                    
                    if(verbose)
                        cout << "operaRadialVelocity: linewidth=" << linewidth << " nm" << endl;
                    
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
                        printf("operaRadialVelocity: Warning: order %d: [Object Spectrum] %d lines detected from input object. Skipping calibration.\n", order, compLines->getnLines());
                        delete(compLines);
						continue;
                    } else {
                        if (verbose) {
                            printf("operaRadialVelocity: order %d: [Object Spectrum] %d lines in object between wl0 = %.2f and wlf = %.2f.\n", order,compLines->getnLines(), wl0, wlf);
                        }
                    }
                    
					/*
                     * Below it reads wavelength and flux from Telluric reference spectrum
                     */
                    if (!telluric_spectrum.empty()) {
                        if (verbose) {
                            cout << "operaRadialVelocity: calculating first order correction (wlshift) using cross-correlation.." << endl;
                        }
						
                        double wlshift = 0;
                        double maxcorr = 0;
                        
                        calculateWavelengthShiftByXCorr(compSpectrum, initialWavelengthRange, initialWavelengthStep, XCorrelationThreshold, &wlshift, &maxcorr);
                        
                        if(verbose)
                            cout << "operaRadialVelocity: order=" << order << " wlshift=" << wlshift <<  " maxcorr=" << maxcorr << endl;
                        
                        double zeroOrderCoeff = wavelengthPolynomial->getCoefficient(0);
                        
                        wavelengthPolynomial->setCoefficient(0,zeroOrderCoeff+wlshift);
                        
                        compSpectrum->setwavelengthsFromCalibration(wavelength);
                        
                        double wl_central = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()/2);
                        double wl0 = compSpectrum->getwavelength(0);
                        double wlf = compSpectrum->getwavelength(compSpectrum->getnSpectralElements()-1);
                        
                        if(verbose)
                            cout << "operaRadialVelocity: new wavelength range: wl0=" << wl0 << ", wl_central=" << wl_central << ", wlf=" << wlf << endl;
                        
                        if (verbose) {
                            cout << "operaRadialVelocity: reading telluric reference spectrum " << telluric_spectrum << endl;
                        }
                        double *wl,*transmission;
                        unsigned npointsInShortTelluricSpectrum = getTelluricSpectrumRange(wl0-fabs(wlf-wl0)/100,wlf+fabs(wlf-wl0)/100,&wl,&transmission);
                        
                        
                        if (verbose) {
                            cout << "operaRadialVelocity: npointsInShortTelluricSpectrum=" << npointsInShortTelluricSpectrum << endl;
                        }
                        
                        /*
                         * Below it degrades the resolution of the atlas to the resolution of object spectrum.
                         * The degradation is done by convolving the spectrum with a gaussian.
                         */
                        if (verbose) {
                            cout << "operaRadialVelocity: convolving spectrum with Gaussian.." << endl;
                        }
                        double *convolvedTelluricFlux = new double[npointsInShortTelluricSpectrum];
                        
                        convolveSpectrumWithGaussian(npointsInShortTelluricSpectrum,wl,transmission,convolvedTelluricFlux,linewidth);
                        
                        if(debug) {
                            for(unsigned i=0;i<npointsInShortTelluricSpectrum;i++){
                                cout << "operaRadialVelocity: " << wl[i] << ' ' << transmission[i] << ' ' << convolvedTelluricFlux[i] << endl;
                            }
                        }
						
                        /*
                         * Below it resamples telluric spectrum.
                         */
                        unsigned nTelluricElements = compSpectrum->getnSpectralElements();
                        double *telluricWavelength = new double[nTelluricElements];
                        double *telluricFlux = new double[nTelluricElements];
                        double *telluricVar = new double[nTelluricElements];
                        
                        for (unsigned i=0; i<nTelluricElements; i++) {
                            telluricWavelength[i] = compSpectrum->getwavelength(i);
                        }
                        operaFitSplineDouble(npointsInShortTelluricSpectrum,wl,convolvedTelluricFlux,nTelluricElements,telluricWavelength,telluricFlux);
						
                        /*
                         * Below it calculates the cross-correlation between the telluric spectrum and a gaussian function.
                         */
                        if (verbose) {
                            cout << "operaRadialVelocity: calculating cross correlation with Gaussian.." << endl;
                        }
                        for (unsigned i=0; i<nTelluricElements; i++) {
                            telluricFlux[i] = 1 - telluricFlux[i];
                        }
                        double *telluricSpectrumXCorr = new double[nTelluricElements];
                        
                        calculateXCorrWithGaussian(nTelluricElements,telluricWavelength,telluricFlux,telluricSpectrumXCorr,linewidth);
						
                        /*
                         * Below it reads the telluric spectrum into an operaSpectralElements class
                         */
                        if (verbose) {
                            cout << "operaRadialVelocity: reading telluric spectrum into operaSpectralElement.." << endl;
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
                                cout << "operaRadialVelocity: " << telluricSpectrum->getwavelength(i)  << " " << telluricSpectrum->getFlux(i)  << " " << telluricSpectrum->getFluxVariance(i) << " " << telluricSpectrum->getXCorrelation(i)<< " " << compfluxvector->getflux(i) << endl;
                            }
                        }
                        
                        /*
                         * Below it creates an operaSpectralLines class for the telluric lines
                         */
                        if (verbose) {
                            cout << "operaRadialVelocity: creating an operaSpectralLines class for the telluric lines.." << endl;
                        }
                        operaSpectralLines *telluricLines = new operaSpectralLines(telluricSpectrum, linewidth, wavelength_disp);
						
                        /*
                         * Below it sets detection thresholds and run the algorithm to detect spectral lines in the reference telluric spectrum
                         */
                        if (verbose) {
                            cout << "operaRadialVelocity: setting detection thresholds and running the algorithm to detect spectral lines in the reference telluric spectrum.." << endl;
                        }
                        double TelluricLocalMaxFilterWidth = LocalMaxFilterWidth*linewidth;
                        double TelluricMinPeakDepth = 0.000001;
                        telluricLines->detectSpectralFeatures(DetectionThreshold,TelluricLocalMaxFilterWidth,TelluricMinPeakDepth);
                        
                        /*
                         * Below it reads the telluric lines information
                         */
                        if (verbose) {
                            cout << "operaRadialVelocity: reading the telluric lines information.." << endl;
                        }
                        
                        if (telluricLines->getnLines() < 3) {	// changed from == 0 Nov 15 2012, since quicksort is called on data & it crashes with only 1 line, for example
                            printf("operaRadialVelocity: Warning:  order %d: [Telluric] %d lines detected from input telluric reference. Skipping calibration.\n", order, telluricLines->getnLines());
                            delete(telluricLines);
							continue;
                        } else {
                            if (verbose) {
                                printf("operaRadialVelocity: order %d: [Telluric] %d lines detected in telluric reference within the range wl0=%.2f, wlf=%.2f.\n", order,  telluricLines->getnLines(), wl0, wlf);
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
                        
                        unsigned numberOfMatchedLines = matchTelluricReferencewithObjectLines(1.0,linewidth,telluricLines,compLines,wlcorrection);
						
                        double *par = (double *)wlcorrection->getVector();
                        double *parerr = (double *)wlcorrection->getErrorVector();
						
                        double zeroOrderModelVariance = wavelengthPolynomial->getCoefficientError(0)*wavelengthPolynomial->getCoefficientError(0);
                        
                        if (verbose) {
                            cout << "operaRadialVelocity: order=" << order << ", numberOfMatchedLines=" << numberOfMatchedLines << ", zeroOrderModelVariance=" << zeroOrderModelVariance << ", chisqr=" << wlcorrection->getChisqr() << endl;
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
                                cout << "operaRadialVelocity: no wavelength correction has been applied using telluric lines." << endl;
                            }
                        }
						
                        delete[] convolvedTelluricFlux;
                        delete[] telluricSpectrumXCorr;
                        delete[] telluricWavelength;
                        delete[] telluricFlux;
                        delete[] telluricVar;
						
						delete(telluricLines);
                        delete(wlcorrection);
                    } // if (!telluric_spectrum.empty())
					delete(compLines);
                } //if (spectralOrder->gethasSpectralElements())
				
                if(debug) {
                    cout << "operaRadialVelocity: corrected wavelength solution for order " << order << endl;
                    cout << "operaRadialVelocity: ";
                    wavelengthPolynomial->printEquation(&cout);
                    cout << endl;
                }
            }
 		}
		
		// output a new wavelength calibration file
		spectralOrders.WriteSpectralOrders(outputWaveFile, Wave);
        
		if (fdata != NULL) {
			fdata->close();
        }
        
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
    " --ordernumber=<INT_VALUE"
	" --minorder=<INT_VALUE>"
	" --maxorder=<INT_VALUE>"
    " --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" -inputWaveFile=Vega-Polar.wcal.gz --outputWaveFile=telluric.wcal.gz --inputObjectSpectrum=/Users/edermartioli/opera/spectra/Vega-Polar/1316861.e.gz --telluric_lines=/Users/edermartioli/opera-1.0//config/opera_HITRAN08-extracted.par.gz --telluric_spectrum=/Users/edermartioli/opera-1.0//config/KPNO_atmtrans.dat.gz --spectralResolution=80000 --initialWavelengthRange=0.1 --initialWavelengthStep=0.0001 --XCorrelationThreshold=0.2\n\n"
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
    "  -O, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
	"  -N, --minorder=<INT_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<INT_VALUE>, Define maximum order number\n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
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


void generateSyntheticTelluricSpectrumUsingGaussianProfile(unsigned np, double *wavelengthVector, double *ouputSpectrum, double resolution) {
    
    for(unsigned i=0;i<np;i++) {    
        ouputSpectrum[i] = 1.0;
    }
    
    double *wl, *intensity;
    
    for(unsigned i=0;i<np;i++) {
        double gaussianWidth = wavelengthVector[i]/resolution;
		
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

unsigned matchTelluricReferencewithObjectLines(double acceptableMismatch,double lineSigma, operaSpectralLines *telluricLines, operaSpectralLines *objectLines, Polynomial *wlcorrection) {
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
                cout << "operaRadialVelocity: " << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << endl;
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
                cout << "operaRadialVelocity: " << center[l] <<  " " << amplitude[l] << " " << sigma[l] << endl;
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
        if(debug){
            for(unsigned index=0; index<nmatch;index++){
				cout << "operaRadialVelocity: "  << telluricMatchedData[index] << ' '  << objectMatchedData[index] << ' ' << objectlineflux[objectMatchedindex[index]] << ' ' << objectlinesigma[objectMatchedindex[index]] << endl;
            }
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

void calculateWavelengthShiftByXCorr(operaSpectralElements *compSpectrum, double DWavelengthRange, double DWavelengthStep, double threshold, double *maxDWavelength, double *maxcorr) {
	
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
    
    double *wl,*trasmission;
    unsigned npointsInShortTelluricSpectrum = getTelluricSpectrumRange(wl0,wlf,&wl,&trasmission);
    
    unsigned nDataPoints = (unsigned)ceil(DWavelengthRange/DWavelengthStep);
    
    double DWavelength = firstDWavelength;
    unsigned jmax = 0;
    *maxcorr = -BIG;    
    *maxDWavelength = 0;
    
    double *maxtelluricSpectrum = NULL;
    
    for(unsigned j=0; j<nDataPoints;j++) {
        for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
            telluricWavelength[i] = compSpectrum->getwavelength(i) + DWavelength;
        }
        
        operaFitSplineDouble(npointsInShortTelluricSpectrum,wl,trasmission,compSpectrum->getnSpectralElements(),telluricWavelength,telluricSpectrum);
        
        double crosscorrelation = operaCrossCorrelation(compSpectrum->getnSpectralElements(),objectSpectrum,telluricSpectrum);
        
        //     generateSyntheticTelluricSpectrumUsingGaussianProfile(compSpectrum->getnSpectralElements(),telluricWavelength,hitranTelluricSpectrum,80000);
        //     double crosscorrelation = operaCrossCorrelation(compSpectrum->getnSpectralElements(),objectSpectrum,hitranTelluricSpectrum);
        
        if(crosscorrelation > *maxcorr && crosscorrelation > threshold) {
            *maxcorr = crosscorrelation;
            *maxDWavelength = DWavelength;
            maxtelluricSpectrum = telluricSpectrum;
            jmax = j;
        }
        
        // Below it prints the wavelength shift, the cross correlation and the threshold
        if(debug)
            cout << DWavelength << " " << crosscorrelation << " " << threshold << endl;
        
        DWavelength+=DWavelengthStep;
    }
    
    // Below it prints out the original wavelengths, the shifted wavelengths for maximum correlation, the object and telluric spectra.
    if(debug) {
        for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
            telluricWavelength[i] = compSpectrum->getwavelength(i) + DWavelength;
            cout << compSpectrum->getwavelength(i) << " " << telluricWavelength[i] << " " << compSpectrum->getFlux(i) << " " << maxtelluricSpectrum[i] << endl;
        }
    }
    
    delete[] objectSpectrum;
    delete[] telluricWavelength;
    delete[] telluricSpectrum;
    //    delete[] hitranTelluricSpectrum;
    
    if(jmax == 0 || jmax == nDataPoints-1) {
		*maxDWavelength = 0;
    }
}

/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaStitchOrders
 Version: 1.0
 Description: This module stitches orders together
 to start up with an OPERA module.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope
 Location: Hawaii USA
 Date: Jan/2014
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

#include "core-espadons/operaStitchOrders.h"
#include "core-espadons/operaWavelengthCalibration.h"

#include "libraries/operaSpectralTools.h"

#include "libraries/operaLibCommon.h"					// for doubleValue_t
#include "libraries/operaLib.h"							// for itos
#include "libraries/operaMath.h"						// for LengthofPolynomial
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFFT.h"							// for operaXCorrelation
#include "libraries/operaFit.h"							// for operaFitSplineDouble
#include "libraries/gzstream.h"							// for gzstream - read compressed reference spectra

bool calculateWavelengthShiftByXCorrInRange(operaSpectralElements *RefSpectrum, operaSpectralElements *compSpectrum, double wl0, double wlf, double DWavelengthStep, double DWavelengthRange,  float nsigcut, double threshold, double &maxDWavelength, double &maxcorr);
    
#define NOTPROVIDED -999

/*! \brief stitch orders together. */
/*! \file operaStitchOrders.cpp */
/*! \package operaStitchOrders */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*!
 * operaStitchOrders
 * \author Eder Martioli
 * \brief Stitch orders.
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
    
    string inputSpectrum;
	string inputWaveFile;
	string outputWaveFile;
	
    int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;
    bool maxorderprovided = false;
    
    int orderOfReference = 51;
    double DWavelengthRange = 0.1;
    double DWavelengthStep = 0.0001;
    double XCorrelationThreshold = 0.05;
    float sigmaThreshold = 1.0;
    
    bool interactive = false;
    
	struct option longopts[] = {
		{"inputSpectrum",       1, NULL, 's'},        
		{"inputWaveFile",       1, NULL, 'w'},
		{"outputWaveFile",      1, NULL, 'o'},

		{"orderOfReference",        1, NULL, 'r'},
		{"DWavelengthRange",        1, NULL, 'R'},
		{"DWavelengthStep",         1, NULL, 'S'},
		{"XCorrelationThreshold",   1, NULL, 'X'},
		{"sigmaThreshold",          1, NULL, 'T'},
        
        
		{"interactive",         0, NULL, 'I'},
		{"plot",                0, NULL, 'p'},
		{"verbose",             0, NULL, 'v'},
		{"debug",               0, NULL, 'd'},
		{"trace",               0, NULL, 't'},
		{"help",                0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "s:w:o:r:R:S:X:T:I::p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt)
		{
			case 's':
				inputSpectrum = optarg;
				break;
			case 'w':
                inputWaveFile = optarg;
				break;
			case 'o':		// output
				outputWaveFile = optarg;
				break;

			case 'r':		
				orderOfReference = atoi(optarg);
				break;
			case 'R':		
				DWavelengthRange = atof(optarg);
				break;
			case 'S':		
				DWavelengthStep = atof(optarg);
				break;
			case 'X':		
				XCorrelationThreshold = atof(optarg);
				break;
			case 'T':		
				sigmaThreshold = atof(optarg);
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
	
	try {
		// we need an input spectrum...
		if (inputSpectrum.empty()) {
			throw operaException("operaStitchOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        if (outputWaveFile.empty()) {
			throw operaException("operaStitchOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

        if (verbose) {
			cout << "operaStitchOrders: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaStitchOrders: inputSpectrum = " << inputSpectrum << endl;            
			cout << "operaStitchOrders: outputWaveFile = " << outputWaveFile << endl;
			cout << "operaStitchOrders: orderOfReference = " << orderOfReference << endl;
			cout << "operaStitchOrders: DWavelengthRange = " << DWavelengthRange << endl;
			cout << "operaStitchOrders: DWavelengthStep = " << DWavelengthStep << endl;
			cout << "operaStitchOrders: XCorrelationThreshold = " << XCorrelationThreshold << endl;
			cout << "operaStitchOrders: sigmaThreshold = " << sigmaThreshold << endl;            

            if(ordernumber != NOTPROVIDED) {
                cout << "operaStitchOrders: ordernumber = " << ordernumber << endl;
            }
		}
        
 		operaSpectralOrderVector spectralOrders(inputSpectrum);
        
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
			cout << "operaStitchOrders: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        float *wlShifts[2];
        wlShifts[0] = new float[MAXORDERS];
        wlShifts[1] = new float[MAXORDERS];
        
        for (unsigned i=0; i<MAXORDERS; i++) {
            wlShifts[0][i] = 0;
            wlShifts[1][i] = 0;
        }
        
        
        unsigned nord = 0;
        int orderIndex[MAXORDERS];
        
        for (unsigned order=(unsigned)minorder; order<=(unsigned)maxorder; order++) {
            if (verbose) {
                cout << "operaStitchOrders: processing order " << order << endl;
            }
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();

                orderIndex[nord] = (int)order;
                
                operaWavelength *wavelength = spectralOrder->getWavelength();
                spectralElements->setwavelengthsFromCalibration(wavelength);
                
                
                if((int)order - 1 > minorder) {
                    operaSpectralOrder *neighborSpectralOrder = spectralOrders.GetSpectralOrder(order-1);
                    if (neighborSpectralOrder->gethasSpectralElements() && neighborSpectralOrder->gethasWavelength()) {
                        operaSpectralElements *neighborElements = neighborSpectralOrder->getSpectralElements();
                        operaWavelength *neighborWavelength = neighborSpectralOrder->getWavelength();
                        neighborElements->setwavelengthsFromCalibration(neighborWavelength);
                        
                        double wl0=0, wlf=0;
                        bool overlap = getOverlappingWLRange(neighborElements,spectralElements,wl0,wlf);
                        if(overlap) {
                            double delta_wl = 0;
                            double maxcorr = 0;
                            bool xcorrstatus = calculateWavelengthShiftByXCorrInRange(neighborElements,spectralElements,wl0,wlf,DWavelengthStep,DWavelengthRange,sigmaThreshold,XCorrelationThreshold,delta_wl,maxcorr);
                            if(xcorrstatus) {
                                wlShifts[0][nord] = (float)delta_wl;
                            }
                        }
                    }
                }
                if((int)order + 1 <= maxorder) {
                    operaSpectralOrder *neighborSpectralOrder = spectralOrders.GetSpectralOrder(order+1);
                    if (neighborSpectralOrder->gethasSpectralElements() && neighborSpectralOrder->gethasWavelength()) {
                        operaSpectralElements *neighborElements = neighborSpectralOrder->getSpectralElements();
                        operaWavelength *neighborWavelength = neighborSpectralOrder->getWavelength();
                        neighborElements->setwavelengthsFromCalibration(neighborWavelength);
                        
                        double wl0=0, wlf=0;
                        bool overlap = getOverlappingWLRange(neighborElements,spectralElements,wl0,wlf);
                        if(overlap) {
                            double delta_wl = 0;
                            double maxcorr = 0;
                            bool xcorrstatus = calculateWavelengthShiftByXCorrInRange(neighborElements,spectralElements,wl0,wlf,DWavelengthStep,DWavelengthRange,sigmaThreshold,XCorrelationThreshold,delta_wl,maxcorr);
                            if(xcorrstatus) {
                               wlShifts[1][nord] = (float)delta_wl;
                            }
                        }
                    }
                }
                nord++;
            }
        }

        double *orderWlShift = new double[MAXORDERS];
        float *shift = new float[MAXORDERS];
        int referenceOrderIndex = NOTPROVIDED;
        unsigned np = 0;
                
        for(unsigned ord=0; ord<nord; ord++) {
            orderWlShift[ord] = 0;
            if(orderIndex[ord] == orderOfReference) {
                referenceOrderIndex = (int)ord;
           }
            if(wlShifts[0][ord] != 0) {
                shift[np++] = wlShifts[0][ord];
            }
            if(wlShifts[1][ord] != 0) {
                shift[np++] = wlShifts[1][ord];
            }
            if(debug) {
                cout << orderIndex[ord] << " " << wlShifts[0][ord] << " " << wlShifts[1][ord] << endl;
            }
        }

        float medianShift = operaArrayMedian(np,shift);
        float medianShiftSigma = operaArrayMedianSigma(np,shift,medianShift);
        if (debug) {
            cout << medianShift << " +/- " << medianShiftSigma << endl;
        }
        
        if(debug)
            cout << "referenceOrderIndex = " << referenceOrderIndex << " refOrder = " << orderIndex[referenceOrderIndex] << endl;

        
        for(int ord=referenceOrderIndex+1; ord<(int)nord; ord++) {
            if(wlShifts[0][ord] > medianShift - medianShiftSigma &&
               wlShifts[0][ord] < medianShift + medianShiftSigma) {
                for(int o = ord; o<(int)nord; o++) {
                    orderWlShift[o] += wlShifts[0][ord];
                }
            }
        }

        for(int ord=referenceOrderIndex-1; ord>=0; ord--) {
            if(wlShifts[1][ord] > medianShift - medianShiftSigma &&
               wlShifts[1][ord] < medianShift + medianShiftSigma) {
                for(int o=ord; o>=0; o--) {
                   orderWlShift[o] += wlShifts[1][ord];
                }
            }
        }

        for(unsigned ord=0; ord<nord; ord++) {
            if(debug) {
                cout << orderIndex[ord] << " " << orderWlShift[ord] << endl;
            }
            
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(orderIndex[ord]);
            
            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                operaWavelength *wavelength = spectralOrder->getWavelength();
                Polynomial *wavelengthPolynomial =  wavelength->getWavelengthPolynomial();
                double newWavelengthPolynomialCoeff = wavelengthPolynomial->getCoefficient(0) + (double)orderWlShift[ord];
                wavelengthPolynomial->setCoefficient(0, newWavelengthPolynomialCoeff);
                spectralOrder->sethasWavelength(true);
            }

        }
        
        spectralOrders.WriteSpectralOrders(outputWaveFile, Wave);
	}
	catch (operaException e) {
		cerr << "operaStitchOrders: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaStitchOrders: " << operaStrError(errno) << endl;
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
	" --inputSpectrum=<SPECTRUM_FILE>"
	" --outputWaveFile=<WAVE_FILE>"
	" --orderOfReference=<INT_VALUE>"
	" --DWavelengthRange=<DBL_VALUE>"
	" --DWavelengthStep=<DBL_VALUE>"
	" --XCorrelationThreshold=<DBL_VALUE>"
	" --sigmaThreshold=<FLT_VALUE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputSpectrum=/Users/edermartioli/opera//calibrations/PolarTest-08BQ00-Oct17/th_EEV1_pol_Normal.e.gz --inputWaveFile=/Users/edermartioli/opera//calibrations/PolarTest-08BQ00-Oct17/EEV1_pol_Normal.wcar.gz --outputWaveFile=/Users/edermartioli/opera//calibrations/PolarTest-08BQ00-Oct17/EEV1_pol_Normal.wcal.gz --orderOfReference=50 --DWavelengthRange=0.1 --DWavelengthStep=0.0001 --XCorrelationThreshold=0.05 --sigmaThreshold=1.0  -v -t -p \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file \n"
	"  -s, --inputSpectrum=<SPECTRUM_FILE>, Input object spectrum file (.s)\n"
	"  -o, --outputWaveFile=<WAVE_FILE>, Output wavelength calibration file to store final solution\n"
	"  -r, --orderOfReference=<INT_VALUE>, Order taken as reference (fixed) for stitching \n"
	"  -R, --DWavelengthRange=<DBL_VALUE>, Wavelength range to search for shift \n"
	"  -S, --DWavelengthStep=<DBL_VALUE>, Wavelength precision to search for shift \n"
	"  -X, --XCorrelationThreshold=<DBL_VALUE>, XCorrelation minimum treshold to accept shift \n"
	"  -T, --sigmaThreshold=<FLT_VALUE>, Threshold in units of sigma to find peak correlation \n"
	"  -I, --interactive=<BOOL>\n\n";
}


bool calculateWavelengthShiftByXCorrInRange(operaSpectralElements *RefSpectrum, operaSpectralElements *compSpectrum, double wl0, double wlf, double DWavelengthStep, double DWavelengthRange,  float nsigcut, double threshold, double &maxDWavelength, double &maxcorr) {
    
    bool status = false;
    
    double *refSpectrumFlux = new double[RefSpectrum->getnSpectralElements()];
    double *refSpectrumWavelength = new double[RefSpectrum->getnSpectralElements()];
    unsigned nref = getSpectrumWithinWLRange(RefSpectrum,wl0,wlf,refSpectrumFlux,refSpectrumWavelength);
    for (unsigned i=0; i<nref; i++) {
        if(isnan(refSpectrumFlux[i])) {
            refSpectrumFlux[i] = 0;
        }
    }
    
    double *compSpectrumFlux = new double[compSpectrum->getnSpectralElements()];
    double *compSpectrumWavelength = new double[compSpectrum->getnSpectralElements()];
    unsigned ncomp = getSpectrumWithinWLRange(compSpectrum,wl0-DWavelengthRange/2.0,wlf+DWavelengthRange/2.0,compSpectrumFlux,compSpectrumWavelength);
    for (unsigned i=0; i<ncomp; i++) {
        if(isnan(compSpectrumFlux[i])) {
            compSpectrumFlux[i] = 0;
        }
    }
    double *compSpectrumWavelength_mod = new double[compSpectrum->getnSpectralElements()];
    double *compSpectrumFlux_mod = new double[RefSpectrum->getnSpectralElements()];

    unsigned nDataPoints = (unsigned)ceil( (DWavelengthRange / DWavelengthStep) + 1.0);
    double DWavelength = -DWavelengthRange/2.0;
    
    unsigned jmax = 0;
    maxcorr = -BIG;
    maxDWavelength = 0;
	
	if (nref == 0 || ncomp ==0) {
		maxcorr = NAN;
		maxDWavelength = NAN;
		return false;
	}

    float *crosscorrelation = new float [nDataPoints];
    float *dwl = new float [nDataPoints];

    for(unsigned j=0; j<nDataPoints;j++) {
        for (unsigned i=0; i<ncomp; i++) {
            compSpectrumWavelength_mod[i] = compSpectrumWavelength[i] + DWavelength;
        }
        
        operaFitSplineDouble(ncomp,compSpectrumWavelength_mod,compSpectrumFlux,nref,refSpectrumWavelength,compSpectrumFlux_mod);
        
        crosscorrelation[j] = (float)(operaCrossCorrelation(nref,refSpectrumFlux,compSpectrumFlux_mod));
        dwl[j] = (float)DWavelength;
                
        if((double)(crosscorrelation[j]) > maxcorr && (double)(crosscorrelation[j]) > threshold) {
            maxcorr = (double)(crosscorrelation[j]);
            maxDWavelength = DWavelength;
            jmax = j;
            status = true;
        }

        if(debug) {
            cout << dwl[j] << " "  << crosscorrelation[j] << endl;
        }
        
        DWavelength+=DWavelengthStep;
    }
    
    double firstMaxDWavelength = maxDWavelength;
    double firstMaxcorr = maxDWavelength;    
  
    if (maxDWavelength == -DWavelengthRange/2.0) {
		maxcorr = NAN;
		maxDWavelength = NAN;
        status = false;
	}
    
    if(debug)
        cout << "calculateWavelengthShiftByXCorrInRange: status=" << status << " wlshift=" << maxDWavelength  <<  " maxcorr=" << maxcorr << endl;

     /*
     * Two tests will be performed:
     *  I) Whether it is an isolated maximum point, i.e. three points before and three points after all must be lower than max
     *
     *  II) Whether it's a significant maximum, i.e. maxcorrelation > median + 3*sigma
     */
    
    unsigned halfslitsize = 5;
    
    double *peakXdata = new double[2*halfslitsize+2];
    double *peakYdata = new double[2*halfslitsize+2];
    unsigned np = 0;
    // Test (I)
    
    if((int)jmax - halfslitsize < 0 || jmax + halfslitsize >= nDataPoints) {
        maxcorr = NAN;
		maxDWavelength = NAN;
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
                maxcorr = NAN;
                maxDWavelength = NAN;
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
                maxcorr = NAN;
                maxDWavelength = NAN;
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
        maxcorr = NAN;
        maxDWavelength = NAN;
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
        maxcorr = a;
        maxDWavelength = x0;
    } else {
        if(firstMaxcorr > threshold) {
            maxDWavelength = firstMaxDWavelength;
            maxcorr = firstMaxcorr;
            status = true;
        }
    }
    
    if(debug)
        cout << "calculateWavelengthShiftByXCorr: maxDWavelength=" << maxDWavelength << endl;
    
    delete[] peakXdata;
    delete[] peakYdata;

    delete[] crosscorrelation;
    delete[] dwl;
    
    delete[] refSpectrumFlux;
    delete[] refSpectrumWavelength;
    delete[] compSpectrumFlux;
    delete[] compSpectrumWavelength_mod;
    delete[] compSpectrumFlux_mod;

	return status;
}




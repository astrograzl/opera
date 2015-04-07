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
#include "libraries/operaArgumentHandler.h"

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

operaArgumentHandler args;

int main(int argc, char *argv[])
{
    string inputSpectrum;
    string inputWaveFile;
    string outputWaveFile;
    int orderOfReference;
    double DWavelengthRange;
    double DWavelengthStep;
    double XCorrelationThreshold;
    double sigmaThreshold;
    bool interactive;
    args.AddRequiredArgument("inputSpectrum", inputSpectrum, "Input object spectrum file (.s)");
	args.AddRequiredArgument("outputWaveFile", outputWaveFile, "Output wavelength calibration file to store final solution");
	args.AddOptionalArgument("inputWaveFile", inputWaveFile, "", "Input wavelength calibration file");
    args.AddOptionalArgument("orderOfReference", orderOfReference, 51, "Order taken as reference (fixed) for stitching");
    args.AddOptionalArgument("DWavelengthRange", DWavelengthRange, 0.1, "Wavelength range to search for shift");
    args.AddOptionalArgument("DWavelengthStep", DWavelengthStep, 0.0001, "Wavelength precision to search for shift");
    args.AddOptionalArgument("XCorrelationThreshold", XCorrelationThreshold, 0.05, "XCorrelation minimum treshold to accept shift");
    args.AddOptionalArgument("sigmaThreshold", sigmaThreshold, 1.0, "Threshold in units of sigma to find peak correlation");
    args.AddSwitch("interactive", interactive, "for interactive plots");
    //" Example: "+string(modulename)+" --inputSpectrum=/Users/edermartioli/opera//calibrations/PolarTest-08BQ00-Oct17/th_EEV1_pol_Normal.e.gz --inputWaveFile=/Users/edermartioli/opera//calibrations/PolarTest-08BQ00-Oct17/EEV1_pol_Normal.wcar.gz --outputWaveFile=/Users/edermartioli/opera//calibrations/PolarTest-08BQ00-Oct17/EEV1_pol_Normal.wcal.gz --orderOfReference=50 --DWavelengthRange=0.1 --DWavelengthStep=0.0001 --XCorrelationThreshold=0.05 --sigmaThreshold=1.0  -v -t -p \n\n"
    
    int ordernumber = NOTPROVIDED;
    int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
	
	try {
		args.Parse(argc, argv);
		// we need an input spectrum...
		if (inputSpectrum.empty()) {
			throw operaException("operaStitchOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        if (outputWaveFile.empty()) {
			throw operaException("operaStitchOrders: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

        if (args.verbose) {
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
        
        if(minorder == NOTPROVIDED) minorder = spectralOrders.getMinorder();
        if(maxorder == NOTPROVIDED) maxorder = spectralOrders.getMaxorder();        
        if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}
        if (args.verbose) cout << "operaStitchOrders: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        float wlShifts[2][MAXORDERS];
        for (unsigned i=0; i<MAXORDERS; i++) {
            wlShifts[0][i] = 0;
            wlShifts[1][i] = 0;
        }
        
        unsigned nord = 0;
        int orderIndex[MAXORDERS];
        
        for (unsigned order=(unsigned)minorder; order<=(unsigned)maxorder; order++) {
            if (args.verbose) cout << "operaStitchOrders: processing order " << order << endl;
            
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

        double orderWlShift[MAXORDERS];
        float shift[MAXORDERS];
        int referenceOrderIndex = NOTPROVIDED;
        unsigned np = 0;
                
        for(unsigned ord=0; ord<nord; ord++) {
            orderWlShift[ord] = 0;
            if(wlShifts[0][ord] != 0) shift[np++] = wlShifts[0][ord];
            if(wlShifts[1][ord] != 0) shift[np++] = wlShifts[1][ord];
            if(orderIndex[ord] == orderOfReference) referenceOrderIndex = (int)ord;
            if(args.debug) {
                cout << orderIndex[ord] << " " << wlShifts[0][ord] << " " << wlShifts[1][ord] << endl;
            }
        }

        float medianShift = operaArrayMedian(np,shift);
        float medianShiftSigma = operaArrayMedianSigma(np,shift,medianShift);
        if (args.debug) {
            cout << medianShift << " +/- " << medianShiftSigma << endl;
			cout << "referenceOrderIndex = " << referenceOrderIndex << " refOrder = " << orderIndex[referenceOrderIndex] << endl;
		}

        
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
            if(args.debug) cout << orderIndex[ord] << " " << orderWlShift[ord] << endl;
            
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

        if(args.debug) {
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
    
    if(args.debug)
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
    if(args.debug)
        cout << "calculateWavelengthShiftByXCorr: Test I status=" << status << endl;
    
    // Test (II)
    
    float medianXcorr = operaArrayMedian(nDataPoints,crosscorrelation);
    float medsigXcorr = operaArrayMedianSigma(nDataPoints,crosscorrelation,medianXcorr);
    
    if(args.debug)
        cout << "calculateWavelengthShiftByXCorr:maxcorr=" << crosscorrelation[jmax] << "  medianXcorr=" << medianXcorr << "  medsigXcorr=" << medsigXcorr << "  median+n*sig=" << medianXcorr + nsigcut*medsigXcorr << endl;
    
    if(crosscorrelation[jmax] < medianXcorr + nsigcut*medsigXcorr) {
        maxcorr = NAN;
        maxDWavelength = NAN;
        status = false;
    }
    
    if(args.debug)
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
    
    if(args.debug)
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




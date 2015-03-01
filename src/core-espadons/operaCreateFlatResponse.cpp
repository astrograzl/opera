/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateFlatResponse
 Version: 1.0
 Description: Flat Response Flux Calibration with Standard or Moon spectrum
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Feb/2015
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
#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/Polynomial.h"
#include "libraries/operaFFT.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaSpectralTools.h"			// void calculateUniformSample, getFluxAtWavelength
#include "libraries/operaArgumentHandler.h"

#define MAXFLUXREFERENCELENGTH 20000
#define MAXNUMBEROFREFWLRANGES 1000

/*! \file operaCreateFlatResponse.cpp */

using namespace std;

operaArgumentHandler args;

/*
 * the reference spectrum
 */
static unsigned nPointsInReferenceSpectrum = 0;
static double referenceWavelength[MAXFLUXREFERENCELENGTH];
static double referenceIntensity[MAXFLUXREFERENCELENGTH];
static double referenceNormIntensity[MAXFLUXREFERENCELENGTH];
static double referenceVariance[MAXFLUXREFERENCELENGTH];  

unsigned getReferenceSpectrumRange(double wl0, double wlf, double **wl, double **flux, double **normflux, double **fluxvar);
unsigned readReferenceSpectrum(string reference_spectrum, double *referenceWavelength, double *referenceIntensity, double *referenceVariance);

void normalizeIntensityByMaximum(unsigned np, double *intensity, double *variance);
void normalizeIntensityByReferenceWavelength(unsigned np, double *intensity, double *wavelength, double *outputNormIntensity, double refWavelength);
double getReferenceFlux(unsigned np, double *intensity, double *wavelength, double refWavelength);

double operaArrayMaxValue_d(unsigned np, const double *xarray, const double *yarray, double *maxx);
unsigned getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, float *refContinuumwl,float *refContinuumflux,float *refContinuumNormflux);
unsigned getReferenceSpectrumRange(unsigned nRefContinuum,double *refContinuumwl,double *refContinuumflux,double *refContinuumNormflux,double wl0,double wlf, double **wl, double **flux, double **normflux);

/*! 
 * operaCreateFlatResponse
 * \author Eder Martioli
 * \brief Flux Calibration with Standard source.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{	
	const double DELTA_WL = 1.0; // wavelength (in nm) range for stiching non-overlapping orders
    const int NOT_PROVIDED = -999;
    
	string inputUncalibratedSpectrum;
	string inputCalibratedSpectrum;
    string inputFlatFluxCalibration;
	string inputWaveFile;
    string inputWavelengthMaskForRefContinuum;
    string inputWavelengthMaskForUncalContinuum;
	string outputFlatResponseFile;
    double wavelengthForNormalization = 548;
    int ordernumber = NOT_PROVIDED;
    int minorder = 22;
    int maxorder = 62;    
	unsigned numberOfPointsInUniformSample = 200;
    unsigned numberOfPointsInUniformRefSample = 70;
    unsigned binsize = 100;
	
	args.AddRequiredArgument("inputUncalibratedSpectrum", inputUncalibratedSpectrum, "Spectrophotometric standard extracted uncalibrated spectrum");
	args.AddRequiredArgument("inputCalibratedSpectrum", inputCalibratedSpectrum, "Spectrophotometric standard template calibrated spectrum");
	args.AddRequiredArgument("inputFlatFluxCalibration", inputFlatFluxCalibration, "");
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "Input wavelength calibration file");
	args.AddRequiredArgument("inputWavelengthMaskForRefContinuum", inputWavelengthMaskForRefContinuum, "");
	args.AddRequiredArgument("inputWavelengthMaskForUncalContinuum", inputWavelengthMaskForUncalContinuum, "");
	args.AddRequiredArgument("outputFlatResponseFile", outputFlatResponseFile, "Output flux calibration conversion file");
	
	args.AddOptionalArgument("wavelengthForNormalization", wavelengthForNormalization, 548, "Wavelength (nm) for normalization of reference spectrum");
	args.AddOptionalArgument("ordernumber", ordernumber, NOT_PROVIDED, "Absolute order number to extract (default=all)");
	args.AddOptionalArgument("minorder", minorder, NOT_PROVIDED, "Define minimum order number");
	args.AddOptionalArgument("maxorder", maxorder, NOT_PROVIDED, "Define maximum order number");
	args.AddOptionalArgument("numberOfPointsInUniformSample", numberOfPointsInUniformSample, 200, "Define lowest order to consider in the fit across orders");
	args.AddOptionalArgument("numberOfPointsInUniformRefSample", numberOfPointsInUniformRefSample, 70, "Define highest order to consider in the fit across orders");
	args.AddOptionalArgument("binsize", binsize, 100, "Number of points to bin for continuum estimate");
	
	//"Example: "+string(modulename)+" --inputCalibratedSpectrum=HR1544_operaFluxCal.dat --inputUncalibratedSpectrum=1515004.e.gz --inputWaveFile=/Users/edermartioli/opera/calibrations/GalileanMoons/OLAPAa_pol_Normal.wcar.gz --outputFlatResponseFile=1515004.fcal.gz --normalizeCalibratedSpectrum=1 --binsize=210 --spectrumDataFilename=1515004.spec --continuumDataFilename=1515004.cont --scriptfilename=1515004fcal.gnu -v"
	
	try {
		// we need an input uncalibrated spectrum...
		if (inputUncalibratedSpectrum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input template spectrum...
		if (inputCalibratedSpectrum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input flat flux calibration spectrum...
		if (inputFlatFluxCalibration.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for ref...
		if (inputWavelengthMaskForRefContinuum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for uncal...
		if (inputWavelengthMaskForUncalContinuum.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an output...
		if (outputFlatResponseFile.empty()) {
			throw operaException("operaCreateFlatResponse: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
        }
        
		if (args.verbose) {
			cout << "operaCreateFlatResponse: input uncalibrated spectrum file = " << inputUncalibratedSpectrum << endl;
			cout << "operaCreateFlatResponse: input calibrated spectrum file = " << inputCalibratedSpectrum << endl;
			cout << "operaCreateFlatResponse: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
			cout << "operaCreateFlatResponse: inputWavelengthMaskForRefContinuum = " << inputWavelengthMaskForRefContinuum << endl;
			cout << "operaCreateFlatResponse: inputWavelengthMaskForUncalContinuum = " << inputWavelengthMaskForUncalContinuum << endl;
			cout << "operaCreateFlatResponse: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaCreateFlatResponse: output flux calibration file = " << outputFlatResponseFile << endl;
            cout << "operaCreateFlatResponse: wavelengthForNormalization= " << wavelengthForNormalization << " nm" << endl;
			cout << "operaCreateFlatResponse: numberOfPointsInUniformSample = " << numberOfPointsInUniformSample << endl;
			cout << "operaCreateFlatResponse: numberOfPointsInUniformRefSample = " << numberOfPointsInUniformRefSample << endl;
            if(ordernumber != NOT_PROVIDED) cout << "operaCreateFlatResponse: ordernumber = " << ordernumber << endl;
            cout << "operaCreateFlatResponse: binsize = " << binsize << endl;
		}
        
		operaSpectralOrderVector spectralOrders(inputUncalibratedSpectrum);
        spectralOrders.ReadSpectralOrders(inputWaveFile);

        if(minorder == NOT_PROVIDED) minorder = spectralOrders.getMinorder();
        if(maxorder == NOT_PROVIDED) maxorder = spectralOrders.getMaxorder();
        if(ordernumber != NOT_PROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}
        
		if (args.verbose) cout << "operaCreateFlatResponse: minorder ="<< minorder << " maxorder=" << maxorder << endl;        

		/*
		 * Flux calibration reference file:
         *
		 * Read reference calibrated spectrum
		 *		lambda vs. intensity, intensityVariance (optional)
		 */
        nPointsInReferenceSpectrum = readReferenceSpectrum(inputCalibratedSpectrum, referenceWavelength, referenceIntensity, referenceVariance);        
        normalizeIntensityByReferenceWavelength(nPointsInReferenceSpectrum,referenceIntensity,referenceWavelength,referenceNormIntensity,wavelengthForNormalization);
        
        //---------------------------------
        // Loop over orders to set maximum number of elements, set wavelength and the number of beams
        // --> maxNElements & NumberofBeams
        unsigned NumberofBeams = spectralOrders.getNumberofBeams(minorder, maxorder);
        
        for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                operaWavelength *wavelength = spectralOrder->getWavelength();
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();                
                SpectralElements->setwavelengthsFromCalibration(wavelength);
            }
        }
        if (args.verbose) cout << "operaCreateFlatResponse: NumberofBeams = " << NumberofBeams << endl;
        if(NumberofBeams == 0) throw operaException("operaCreateFlatResponse: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);

        
        //---------------------------------
        // Correct flat-field
        if (!inputFlatFluxCalibration.empty()) {
            spectralOrders.correctFlatField(inputFlatFluxCalibration, minorder, maxorder, false);
        }
        
        //---------------------------------
        // Calculate a clean sample of the continuum from the ref spectrum
        float *refContinuumwl = new float[MAXNUMBEROFREFWLRANGES];
        float *refContinuumflux = new float[MAXNUMBEROFREFWLRANGES];
        float *refContinuumNormflux = new float[MAXNUMBEROFREFWLRANGES];

        unsigned nRefContinuum = getContinuumFromInputReferenceSpectrum(inputWavelengthMaskForRefContinuum,refContinuumwl,refContinuumflux,refContinuumNormflux);

        float *uniformRef_wl = new float[numberOfPointsInUniformRefSample];
        float *uniformRef_flux = new float[numberOfPointsInUniformRefSample];
        
        calculateUniformSample(nRefContinuum,refContinuumwl,refContinuumflux,numberOfPointsInUniformRefSample,uniformRef_wl,uniformRef_flux);

        float *uniform_wl = new float[numberOfPointsInUniformSample];
        float *uniform_flux = new float[numberOfPointsInUniformSample];
        float *uniform_Beamflux[MAXNUMBEROFBEAMS];
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            uniform_Beamflux[beam] = new float[numberOfPointsInUniformSample];
        }
        
        spectralOrders.calculateCleanUniformSampleOfContinuum(minorder,maxorder,binsize,DELTA_WL,inputWavelengthMaskForUncalContinuum,numberOfPointsInUniformSample,uniform_wl,uniform_flux,uniform_Beamflux,TRUE);

        float *calibratedModelFlux = new float[numberOfPointsInUniformSample];
        operaFitSpline(numberOfPointsInUniformRefSample,uniformRef_wl,uniformRef_flux,numberOfPointsInUniformSample,uniform_wl,calibratedModelFlux);

        float *flatResp = new float[numberOfPointsInUniformSample];

        for(unsigned i=0;i<numberOfPointsInUniformSample;i++) {
            flatResp[i] = uniform_flux[i]/calibratedModelFlux[i];
        }
        double flatRespForNormalization = getFluxAtWavelength(numberOfPointsInUniformSample,uniform_wl,flatResp,wavelengthForNormalization);

        /*
         * and write out flatresponse (LE *.s)
         */
		ofstream frespoutput(outputFlatResponseFile.c_str());
        frespoutput << "***" << endl;
        frespoutput << numberOfPointsInUniformSample << " 1" << endl;
        for(unsigned i=0;i<numberOfPointsInUniformSample;i++) {
            flatResp[i] /= flatRespForNormalization;
            frespoutput << uniform_wl[i] << ' ' << flatResp[i] << endl;
        }
        frespoutput.close();
        
	}
	catch (operaException e) {
		cerr << "operaCreateFlatResponse: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateFlatResponse: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/*
 * Read the the full reference spectrum
 */
unsigned readReferenceSpectrum(string reference_spectrum, double *referenceWavelength, double *referenceIntensity, double *referenceVariance) {
	ifstream astream;
	string dataline;
    
	double tmpwl = -1.0; 
	double tmpi = -1.0; 
	unsigned np = 0;
	
	astream.open(reference_spectrum.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi);
                    
                    referenceWavelength[np] = tmpwl;
                    referenceIntensity[np] = tmpi;
                    referenceVariance[np] = tmpi;
                    np++;  
                }	// skip comments
            }
		} // while (astream.good())
        
		if (np > 0) {
			if (args.verbose) {
				printf("          [Reference] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, referenceWavelength[0], referenceWavelength[np/2], referenceWavelength[np-1]);
			}
		} else {
			printf("          [Reference] no points found in flux reference file.\n");
		}
		astream.close();
	}	// if (astream.open())
	return np;
}

/*
 * get a subset of the reference spectrum for this order only, between wl0 and wlf
 */
unsigned getReferenceSpectrumRange(double wl0, double wlf, double **wl, double **flux, double **normflux, double **fluxvar) {
	unsigned firstline = 0;
	unsigned np = 0;
    
	for (np=0; np<nPointsInReferenceSpectrum; np++) {
		if (referenceWavelength[np] >= wl0) {
			if (firstline == 0 && np) {
                    *flux = &referenceIntensity[np-1];
                    *normflux = &referenceNormIntensity[np-1];
                    *wl = &referenceWavelength[np-1];
                    *fluxvar = &referenceVariance[np-1];
                    firstline = np-1;
			}
			if (referenceWavelength[np] > wlf) {
                np++;
				break;
            }
		}
	}
	if (np == nPointsInReferenceSpectrum) np--;
	if (np > firstline) return (np-firstline);
	return 0;
}

void normalizeIntensityByReferenceWavelength(unsigned np, double *intensity, double *wavelength, double *outputNormIntensity, double refWavelength) {
    double referenceFlux = getReferenceFlux(np,intensity,wavelength,refWavelength);
    
	for(unsigned i=0;i<np;i++) {
        outputNormIntensity[i] = intensity[i]/referenceFlux;
	}
}

double getReferenceFlux(unsigned np, double *intensity, double *wavelength, double refWavelength) {
    
    float *wavelengthData_f = new float[np];
    float *fluxData_f = new float[np];
    
    for(unsigned i=0;i<np;i++) {
        wavelengthData_f[i] = (float)wavelength[i];
        fluxData_f[i] =  (float)intensity[i];
    }
    
    unsigned nElements = 1;
    
    float *referenceFlux = new float[nElements];
    float *referencewl = new float[nElements];
    
    for(unsigned i=0;i<nElements;i++) {
        referencewl[i] = (float)refWavelength;
    }
    operaFitSpline(np,wavelengthData_f,fluxData_f,nElements,referencewl,referenceFlux);
    
    double outputflux = (double)(referenceFlux[0]);
    
    delete[] wavelengthData_f;
    delete[] fluxData_f;
    delete[] referenceFlux;
    delete[] referencewl;
    
	return outputflux;
}

void normalizeIntensityByMaximum(unsigned np, double *intensity, double *variance) {
	double maxIntensity = -3.4e+38;
	
	for(unsigned i=0;i<np;i++) {
		if(intensity[i] > maxIntensity)
			maxIntensity = intensity[i];
	}
	
	for(unsigned i=0;i<np;i++) {
        intensity[i] /= maxIntensity;
        variance[i] /= maxIntensity;
	}
}


double operaArrayMaxValue_d(unsigned np, const double *xarray, const double *yarray, double *maxx) {
	double ymax = -3.4e+38;
	double xmax = 0;
    
	while (np--) {
		if(*yarray > ymax) {
			ymax = *yarray;
            xmax = *xarray;
        }
		yarray++;
        xarray++;
	}
    *maxx = xmax;
	return ymax;
}


unsigned getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, float *refContinuumwl,float *refContinuumflux,float *refContinuumNormflux) {
    
    double *wl0_vector = new double[MAXNUMBEROFREFWLRANGES];
    double *wlf_vector = new double[MAXNUMBEROFREFWLRANGES];
    
    unsigned nRangesInWLMask = readContinuumWavelengthMask(inputWavelengthMaskForRefContinuum,wl0_vector,wlf_vector);
    
    unsigned nTotalPoints = 0;
    for(unsigned k=0;k<nRangesInWLMask; k++){
        
        double *refwl = NULL, *refflux = NULL, *refnormflux = NULL, *reffluxvar = NULL;
        unsigned nPointsInReference = getReferenceSpectrumRange(wl0_vector[k],wlf_vector[k],&refwl,&refflux,&refnormflux,&reffluxvar);
        
        double ref_wl=0;
        double ref_maxFlux = 0;
        double ref_maxNormFlux = 0;
        
        ref_maxFlux = operaArrayMaxValue_d(nPointsInReference,refwl,refflux,&ref_wl);
        ref_maxNormFlux = operaArrayMaxValue_d(nPointsInReference,refwl,refnormflux, &ref_wl);
        
        if(args.debug) {
            cout << k << " "
            << nPointsInReference << " "
            << wl0_vector[k] << " "
            << wlf_vector[k] << " "
            << ref_wl << " "
            << ref_maxFlux << " "
            << ref_maxNormFlux << endl;
        }
        
        if(ref_wl && ref_maxFlux && ref_maxNormFlux) {
            refContinuumwl[nTotalPoints] = (float)ref_wl;
            refContinuumflux[nTotalPoints] = (float)ref_maxFlux;
            refContinuumNormflux[nTotalPoints] = (float)ref_maxNormFlux;
            nTotalPoints++;
        }
    }
    
    delete[] wl0_vector;
    delete[] wlf_vector;
    
    return nTotalPoints;
}

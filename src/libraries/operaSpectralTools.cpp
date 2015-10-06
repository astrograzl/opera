/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Library name: operaSpectralTools - common C++ library functions
Version: 1.0
Author(s): CFHT OPERA team
Affiliation: Canada France Hawaii Telescope 
Location: Hawaii USA
Date: Aug/2011

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

/*!
 * operaSpectralTools
 * \author Eder Martioli
 * \brief operaSpectralTools - common C++ library functions.
 * \file operaSpectralTools.cpp
 * \ingroup libraries
 */


#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaSpectralTools.h"

#include "libraries/operaSpectralLines.h"       // for operaSpectralLines
#include "libraries/operaSpectralFeature.h"    // for operaSpectralFeature
#include "libraries/operaLibCommon.h"           // for SPEED_OF_LIGHT_M
#include "libraries/operaStats.h"               // for operaCrossCorrelation
#include "libraries/operaFFT.h"                 // for operaXCorrelation

using namespace std;


/*
 * calculateXCorrWithGaussian(unsigned np, double *wavelength, double *flux, double *outputXcorr, double sigma)
 * \brief This function calculates the cross-correlation between an input spectrum and gaussian function
 * \param unsigned np
 * \param double *wavelength
 * \param double *flux
 * \param return double *outputXcorr
 * \param double sigma
 * \return void
 */
void calculateXCorrWithGaussian(unsigned np, double *wavelength, double *flux, double *outputXcorr, double sigma) {
    
    double wlstep;
    
    //set up window function
    double *windowFunc = (double *) malloc (np * sizeof(double));
    double *mainFunc = (double *) malloc (np * sizeof(double));
    
    for(unsigned i=0;i<np;i++) {
        if(i==0) {
            wlstep = fabs(wavelength[i+1] - wavelength[0]);
        } else if (i==np-1) {
            wlstep = fabs(wavelength[np-1] - wavelength[i-1]);
        } else {
            wlstep = fabs(wavelength[i+1] - wavelength[i-1])/2.0;
        }
        
        unsigned window = (unsigned)ceil(4*sigma/(2*wlstep));
        
        if (window > np/2) {
            window = np/2;
        }
        //calculate window function
        for(unsigned j=0;j<2*window;j++) {
            windowFunc[j] = exp(-double((j - window)*(j - window))/(2.0*double(2*window*2*window)/(4.0*4.0)))/(sqrt(2.0*M_PI)*(double(window)/2));
        }
    
        //calculate xcorrelation function
        for(unsigned j=0;j<2*window;j++) {
            if((i-window+j) >= 0 && (i-window+j) < np) {
                mainFunc[j] = flux[i-window+j];
            } else {
                mainFunc[j] = 0; // zeropad function
            }
        }
        outputXcorr[i] = operaCrossCorrelation(2*window,windowFunc,mainFunc);
    }
    
    free(windowFunc);
    free(mainFunc);
}

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux) {
	double maxflux = 0;
    
    for (unsigned i=0; i<nLines; i++) {
		if(lineflux[i] > maxflux) {
            maxflux = lineflux[i];
        }
	}
	for (unsigned i=0; i<nLines; i++) {
		lineflux[i] /= maxflux/100;
	}
}

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes and respective variances by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \param double *linevariance
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux, double *linevariance) {
	double maxflux = 0;
    
    for (unsigned i=0; i<nLines; i++) {
		if(lineflux[i] > maxflux) {
            maxflux = lineflux[i];
        }
	}
	for (unsigned i=0; i<nLines; i++) {
		lineflux[i] /= maxflux/100;
        linevariance[i] /= (maxflux/100)*(maxflux/100);
	}
}

/*
 * convolveSpectrumWithGaussian(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double sigma)
 * \brief This function calculates the convolution between an input spectrum and a gaussian function
 * \param unsigned np
 * \param double *wavelength
 * \param double *flux
 * \param return double *convolvedSpectrum
 * \param double sigma
 * \return void
 */
void convolveSpectrumWithGaussian(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double sigma) {
    double wlstep;
	
    memset(convolvedSpectrum, 0, sizeof(double)*np);

    for(unsigned i=0;i<np;i++) {
        if(i==0) {
            wlstep = fabs(wavelength[i+1] - wavelength[0]);
        } else if (i==np-1) {
            wlstep = fabs(wavelength[np-1] - wavelength[i-1]);
        } else {
            wlstep = fabs(wavelength[i+1] - wavelength[i-1])/2.0;
        }
        
        unsigned window = (unsigned)ceil(4*sigma/(2*wlstep));
        if (window > np/2) {
            window = np/2;
        }
        double weighSum = 0;
        for(unsigned j=0;j<2*window;j++) {
            if((i-window+j) >= 0 && (i-window+j) < np) {
                convolvedSpectrum[i] += flux[i-window+j]*exp(-((wavelength[i-window+j] - wavelength[i])*(wavelength[i-window+j] - wavelength[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma);
                weighSum += exp(-((wavelength[i-window+j] - wavelength[i])*(wavelength[i-window+j] - wavelength[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma);
            }
        }
        if (weighSum) {
            convolvedSpectrum[i] /= weighSum;
        }
    }
}


/*
 * convolveSpectrumWithGaussianByResolution(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double spectralResolution)
 * \brief This function calculates the convolution between an input spectrum and a gaussian function using the spectral Resolution to calculate line width
 * \param unsigned np
 * \param double *wavelength
 * \param double *flux
 * \param return double *convolvedSpectrum
 * \param double spectralResolution
 * \return void
 */
void convolveSpectrumWithGaussianByResolution(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double spectralResolution) {
    double wlstep;
    
    memset(convolvedSpectrum, 0, sizeof(double)*np);
    
    for(unsigned i=0;i<np;i++) {
        
        double sigma = wavelength[i] / spectralResolution;
        
        if(i==0) {
            wlstep = fabs(wavelength[i+1] - wavelength[0]);
        } else if (i==np-1) {
            wlstep = fabs(wavelength[np-1] - wavelength[i-1]);
        } else {
            wlstep = fabs(wavelength[i+1] - wavelength[i-1])/2.0;
        }
        
        unsigned window = (unsigned)ceil(4*sigma/(2*wlstep));
        if (window > np/2) {
            window = np/2;
        }
        double weighSum = 0;
        for(unsigned j=0;j<2*window;j++) {
            if((i-window+j) >= 0 && (i-window+j) < np) {
                convolvedSpectrum[i] += flux[i-window+j]*exp(-((wavelength[i-window+j] - wavelength[i])*(wavelength[i-window+j] - wavelength[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma);
                weighSum += exp(-((wavelength[i-window+j] - wavelength[i])*(wavelength[i-window+j] - wavelength[i])/(2*sigma*sigma)))/(sqrt(2*M_PI)*sigma);
            }
        }
        if (weighSum) {
            convolvedSpectrum[i] /= weighSum;
        }
    }
}


/*
 * double convertVacuumToAirWavelength(double vac_wl)
 * \brief This function converts from vacuum to air wavelength using the IAU standard for conversion
 * \brief  from air to vacuum wavelengths as given in Morton (1991, ApJS, 77, 119)
 * \param double vac_wl (in Angstrom)
 * \return double air_wl (in Angstrom)
 */
double convertVacuumToAirWavelength(double vac_wl) {
    return vac_wl/(1.0 + 0.0002735182 + (131.4182/(vac_wl*vac_wl)) + (276249000.0/(vac_wl*vac_wl*vac_wl*vac_wl)));
}


/*
 * double calculateBlackBodyVFlux(double Temperature)
 * \brief This function calculates the flux (in ph/(m^2 s m)) for a Black Body
 * \param double Temperature (in Kelvin)
 * \return double flux (in ph/(m^2 s m))
 */
double calculateBlackBodyVFlux(double Temperature) {
    
    double lamb[24], Vresp[24];
    double dlambda = 1e-8;
    
    for(unsigned i=0;i<24;i++) {
        lamb[i] = 470e-9 + (double)i*dlambda;
    }
    // Below is the spectral transmission for the V-band Johnson filter
    Vresp[0] = 0.000;
    Vresp[1] = 0.030;
    Vresp[2] = 0.163;
    Vresp[3] = 0.458;
    Vresp[4] = 0.780;
    Vresp[5] = 0.967;
    Vresp[6] = 1.000;
    Vresp[7] = 0.973;
    Vresp[8] = 0.898;
    Vresp[9] = 0.792;
    Vresp[10] = 0.684;
    Vresp[11] = 0.574;
    Vresp[12] = 0.461;
    Vresp[13] = 0.359;
    Vresp[14] = 0.270;
    Vresp[15] = 0.197;
    Vresp[16] = 0.135;
    Vresp[17] = 0.081;
    Vresp[18] = 0.045;
    Vresp[19] = 0.025;
    Vresp[20] = 0.017;
    Vresp[21] = 0.013;
    Vresp[22] = 0.009;
    Vresp[23] = 0.000;
    
    double fluxV = 0;
    
    for(unsigned i=0;i<24;i++) {
        fluxV += PlanckFunction(Temperature,lamb[i])*Vresp[i]*dlambda;
    }
    return fluxV;
}

/*
 * double planck(double T, double wl)
 * \brief This function returns the spectral radiance of a black body in photons/(s m^2 dlambda).
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double PlanckFunction(double T, double wl) {
    // http://spiff.rit.edu/classes/phys317/lectures/planck.html
    double K_B = 1.380658e-23;
    double H_PLANCK = 6.6260755e-34;
    
    // in units of W/(m^2 dlambda)
    //  flux = (2*M_PI*H_PLANCK*SPEED_OF_LIGHT_M*SPEED_OF_LIGHT_M)/(pow(wl,5)*(exp(H_PLANCK*SPEED_OF_LIGHT_M/(wl*K_B*T)) - 1.));
    
    // in units of photons/(s . m^2  dlambda)
    double flux = (2*M_PI*SPEED_OF_LIGHT_M)/(pow(wl,4)*(exp(H_PLANCK*SPEED_OF_LIGHT_M/(wl*K_B*T)) - 1.));
    
    return flux;
}


/*
 * double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T)
 * \brief This function returns the integrated spectral flux of a black body in photons/(s m^2) for a given
 * \brief wavelength range using the simpson method.
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T) {
    double f;
    unsigned N = 1000;
    double h = fabs(wlf - wl0)/((double)N - 1.);
    double xi = wl0;
    
    double sum = 0;
    for(unsigned i=0;i<=N;i++)
    {
        if(i==0) {
            f = PlanckFunction(T, wl0);
        } else if(i==N) {
            f = PlanckFunction(T, wlf);
        } else {
            xi += h;
            double aux = fmodf((double)(i),2.0);
            if(aux != 0) {
                f = 4.*PlanckFunction(T,xi);
            } else {
                f = 2.*PlanckFunction(T,xi);
            }
        }
        sum += (h/3.)*f;
    }
    
    return sum;
}


double getFactorToMatchFluxesBetweenElements(operaSpectralElements *refElements,operaSpectralElements *elementsToMatch,double delta_wl) {
    bool debug = false;
    double ref_wl0 = refElements->getwavelength(0);
    double ref_wlf = refElements->getwavelength(refElements->getnSpectralElements()-1);
    double elem2match_wl0 = elementsToMatch->getwavelength(0);
    double elem2match_wlf = elementsToMatch->getwavelength(elementsToMatch->getnSpectralElements()-1);
    
    if(ref_wl0 >= ref_wlf || elem2match_wl0 >= elem2match_wlf) {
        throw operaException("getOverlappingWLRange: initial wl must not be greater than final wl. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    // find intersecting range:
    bool elementsIntersect = false;
    double intersect_wl0 = 0;
    double intersect_wlf = 0;
    
    if(elem2match_wl0 >= ref_wl0 && elem2match_wl0 <= ref_wlf) {
        intersect_wl0 = elem2match_wl0;
        intersect_wlf = ref_wlf;
        elementsIntersect = true;
    }
    
    if (elem2match_wlf >= ref_wl0 && elem2match_wlf <= ref_wlf) {
        intersect_wlf = elem2match_wlf;
        if(intersect_wl0 == 0) {
            intersect_wl0 = ref_wl0;
        }
        elementsIntersect = true;
    }
    
    if(debug) {
        cout << "refElements     = " << ref_wl0 << " -> " << ref_wlf << endl;
        cout << "elementsToMatch = " << elem2match_wl0 << " -> " << elem2match_wlf << endl;
        if (elementsIntersect == true) {
            cout << "Intersection    = " << intersect_wl0 << " -> " << intersect_wlf << endl;
        } else {
            cout << " Elements Do Not Intersect" << endl;
        }
    }
    
    
    // first collect data:
    unsigned refNElements = refElements->getnSpectralElements();
    float *refFlux = new float[refNElements];
    unsigned nref = 0;
    
    unsigned nElements = elementsToMatch->getnSpectralElements();
    float *elemFlux = new float[nElements];
    unsigned nelem = 0;
    double factorToMatch = 1.0;
    
    if (elementsIntersect == true) {
        
        for(unsigned elemIndex=0; elemIndex<refNElements;elemIndex++) {
            if(refElements->getwavelength(elemIndex) >= intersect_wl0 &&
               refElements->getwavelength(elemIndex) <= intersect_wlf) {
                refFlux[nref] = refElements->getFlux(elemIndex);
                nref++;
            }
        }
        
        double medianRefFlux = (double)operaArrayMedian(nref,refFlux);
        
        for(unsigned elemIndex=0; elemIndex<nElements;elemIndex++) {
            if(elementsToMatch->getwavelength(elemIndex) >= intersect_wl0 &&
               elementsToMatch->getwavelength(elemIndex) <= intersect_wlf) {
                elemFlux[nelem] = elementsToMatch->getFlux(elemIndex);
                nelem++;
            }
        }
        double medianElemFlux = (double)operaArrayMedian(nelem,elemFlux);
        
        factorToMatch = medianRefFlux/medianElemFlux;
        if(debug) {
            cout << medianRefFlux << " " << medianElemFlux << " " << factorToMatch << endl;
        }
    } else {
        // find out which element set comes first
        // to find the wavelength between the two orders

        double refinflimit_wl, refsuplimit_wl;
        double eleminflimit_wl, elemsuplimit_wl;
        
        if(ref_wlf < elem2match_wl0) { // ref comes first
            refinflimit_wl = ref_wlf - delta_wl;
            refsuplimit_wl = ref_wlf;
            eleminflimit_wl = elem2match_wl0;
            elemsuplimit_wl = elem2match_wl0 + delta_wl;
            
        } else {            
            refinflimit_wl = ref_wl0;
            refsuplimit_wl = ref_wl0 + delta_wl;
            eleminflimit_wl = elem2match_wlf - delta_wl;
            elemsuplimit_wl = elem2match_wlf;
        }
        if(debug) {
            cout << "refinflimit_wl=" << refinflimit_wl << " refsuplimit_wl=" << refsuplimit_wl << endl;
            cout << "eleminflimit_wl=" << eleminflimit_wl << " elemsuplimit_wl=" << elemsuplimit_wl << endl;
        }
        
        for(unsigned elemIndex=0; elemIndex<refNElements;elemIndex++) {
            if(refElements->getwavelength(elemIndex) >= refinflimit_wl &&
               refElements->getwavelength(elemIndex) <= refsuplimit_wl) {
                refFlux[nref] = refElements->getFlux(elemIndex);
                nref++;
            }
        }
        
        double medianRefFlux = (double)operaArrayMedian(nref,refFlux);
        
        for(unsigned elemIndex=0; elemIndex<nElements;elemIndex++) {
            if(elementsToMatch->getwavelength(elemIndex) >= eleminflimit_wl &&
               elementsToMatch->getwavelength(elemIndex) <= elemsuplimit_wl) {
                elemFlux[nelem] = elementsToMatch->getFlux(elemIndex);
                nelem++;
            }
        }
        double medianElemFlux = (double)operaArrayMedian(nelem,elemFlux);
        
        factorToMatch = medianRefFlux/medianElemFlux;
        if(debug) {
            cout << medianRefFlux << " " << medianElemFlux << " " << factorToMatch << endl;
        }
    }
    delete[] refFlux;
    delete[] elemFlux;
    
    return factorToMatch;
}

void calculateUniformSample(unsigned np,float *wl,float *flux, unsigned npout, float *uniform_wl, float *uniform_flux) {
    float wl0 = wl[0];
    float wlf = wl[np-1];
    
    float wlstep = fabs(wlf - wl0)/(float)(npout-1);
    unsigned lastk=0;
    
    for(unsigned i=0;i<npout;i++) {
        uniform_wl[i] = wl0 + (float)i*wlstep;
        
        for (unsigned k=lastk; k<np;k++) {
            if(wl[k] >= uniform_wl[i] && k>0) {
                float slope = (flux[k] - flux[k - 1])/ (wl[k] - wl[k - 1]);
                float intercept = flux[k] - slope*wl[k];
                uniform_flux[i] = intercept + slope*uniform_wl[i];
                lastk = k-1;
                break;
            }
        }
    }
    uniform_flux[npout-1] = flux[np-1]; //otherwise we sometimes miss the last point due to rounding errors
}

float getFluxAtWavelength(unsigned np,float *wl,float *flux,float wavelengthForNormalization) {
    float wl0 = wl[0];
    float wlf = wl[np-1];
    
    if(wavelengthForNormalization < wl0 || wavelengthForNormalization > wlf) {
        cout << wl0 << " " << wavelengthForNormalization << " " << wlf << endl;

        throw operaException("operaSpectralTools: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    float fluxAtwavelengthForNormalization = 0;
    
    for (unsigned k=0; k<np;k++) {
        if(wl[k] >= wavelengthForNormalization && k>0) {
            float slope = (flux[k] - flux[k - 1])/ (wl[k] - wl[k - 1]);
            float intercept = flux[k] - slope*wl[k];
            fluxAtwavelengthForNormalization = intercept + slope*wavelengthForNormalization;
            break;
        }
    }
    return  fluxAtwavelengthForNormalization;
}

/*
 * Read mask
 */
unsigned readContinuumWavelengthMask(string wavelength_mask, double *wl0, double *wlf) {
	ifstream astream;
	string dataline;
    
	double tmpwl0 = -1.0;
	double tmpwlf = -1.0;
	unsigned np = 0;
	
	astream.open(wavelength_mask.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl0, &tmpwlf);
                    wl0[np] = tmpwl0;
                    wlf[np] = tmpwlf;
                    np++;
                }	// skip comments
            }
		} // while (astream.good())
		astream.close();
	}	// if (astream.open())
	return np;
}

unsigned getSpectrumWithinWLRange(operaSpectralElements *inputSpectrum, double wl0, double wlf, double *outputFlux, double *outputWavelength) {
    if(!inputSpectrum->getHasWavelength()) {
        throw operaException("getSpectrumWithinWLRange: no wavelength in input SpectralElements. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    unsigned np = 0;
    for (unsigned i=0; i<inputSpectrum->getnSpectralElements(); i++) {
        if(inputSpectrum->getwavelength(i) >= wl0 &&
           inputSpectrum->getwavelength(i) <= wlf ) {
            
            outputFlux[np] = inputSpectrum->getFlux(i);
            outputWavelength[np] = inputSpectrum->getwavelength(i);
            np++;
        }
    }
    return np;
}

bool getOverlappingWLRange(operaSpectralElements *refElements, operaSpectralElements *elementsToMatch, double &wl0, double &wlf) {
    bool debug = false;
    double ref_wl0 = refElements->getwavelength(0);
    double ref_wlf = refElements->getwavelength(refElements->getnSpectralElements()-1);
    double elem2match_wl0 = elementsToMatch->getwavelength(0);
    double elem2match_wlf = elementsToMatch->getwavelength(elementsToMatch->getnSpectralElements()-1);
    
    if(ref_wl0 > ref_wlf || elem2match_wl0 > elem2match_wlf) {
        throw operaException("getOverlappingWLRange: initial wl must not be greater than final wl. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    // find overlapping range:
    bool overlap = false;
    double intersect_wl0 = 0;
    double intersect_wlf = 0;
    
    if(elem2match_wl0 >= ref_wl0 && elem2match_wl0 <= ref_wlf) {
        intersect_wl0 = elem2match_wl0;
        intersect_wlf = ref_wlf;
        overlap = true;
    }
    
    if (elem2match_wlf >= ref_wl0 && elem2match_wlf <= ref_wlf) {
        intersect_wlf = elem2match_wlf;
        if(intersect_wl0 == 0) {
            intersect_wl0 = ref_wl0;
        }
        overlap = true;
    }
    
    wl0 = intersect_wl0;
    wlf = intersect_wlf;
    
    if(debug) {
        cout << "refElements     = " << ref_wl0 << " -> " << ref_wlf << endl;
        cout << "elementsToMatch = " << elem2match_wl0 << " -> " << elem2match_wlf << endl;
        if (overlap == true) {
            cout << "Intersection    = " << intersect_wl0 << " -> " << intersect_wlf << endl;
        } else {
            cout << " Elements do not overlap" << endl;
        }
    }
    
    return overlap;
}


unsigned detectSpectralLinesInSpectralOrder(operaSpectralOrder *spectralOrder, double *linecenter, double *linecenterError, double *lineflux, double *linesigma, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip, double spectralResolution,bool emissionSpectrum) {
    
    double linewidth = 0;
    double linewidth_err = 0;

    if (!spectralOrder->gethasSpectralElements() || !spectralOrder->gethasWavelength()) {
        throw operaException("detectSpectralLinesInSpectralOrder: order has no elements/wavelength. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned nlinesinorder = 0;
    
    operaSpectralElements *compSpectrum = spectralOrder->getSpectralElements();
    
    double *saveSpectrum = NULL;
    
    operaWavelength *wavelength =  spectralOrder->getWavelength();
    compSpectrum->setwavelengthsFromCalibration(wavelength);
    
    if(!emissionSpectrum) {
        saveSpectrum = new double[compSpectrum->getnSpectralElements()];
        // transform spectrum to emission, so we can use the line detection algorithm
        for(unsigned indexElem=0;indexElem < compSpectrum->getnSpectralElements(); indexElem++) {
            saveSpectrum[indexElem] = compSpectrum->getFlux(indexElem);
            
            double invertedFlux = 1.0 - compSpectrum->getFlux(indexElem) + DetectionThreshold;
            if(invertedFlux < 0 || isnan(invertedFlux)) {
                invertedFlux = 0.0;
            }
            compSpectrum->setFlux(invertedFlux,indexElem);
            
            //cout << indexElem << " " << compSpectrum->getwavelength(indexElem) << " " << invertedFlux << " " << saveSpectrum[indexElem] << endl;
        }
    }
    
    linewidth = wavelength->getcentralWavelength()/spectralResolution;
    
    operaSpectralLines compLines(compSpectrum,linewidth,wavelength_disp);
    
    operaFluxVector *compfluxvector = compSpectrum->getFluxVector();
    
#ifdef PRINT_DEBUG
    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
        cout << compSpectrum->getdistd(i) << " " << compSpectrum->getwavelength(i) << " " << compfluxvector->getflux(i) << " " << compSpectrum->getXCorrelation(i) << endl;
    }
#endif
    if(!compSpectrum->getHasXCorrelation()){
        double compSpectrumwl[MAXPOINTSINSIMULATEDSPECTRUM];
        double compSpectrumflux[MAXPOINTSINSIMULATEDSPECTRUM];
        double compXcorr[MAXPOINTSINSIMULATEDSPECTRUM];
        
        for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
            compSpectrumwl[i] = compSpectrum->getwavelength(i);
            compSpectrumflux[i] = compfluxvector->getflux(i);
        }
        
        calculateXCorrWithGaussian(compSpectrum->getnSpectralElements(), compSpectrumwl, compSpectrumflux, compXcorr, linewidth);
        for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
            compSpectrum->setXCorrelation(compXcorr[i], i);
        }
        compSpectrum->setHasXCorrelation(true);
    }
    double CompLocalMaxFilterWidth = LocalMaxFilterWidth*(linewidth);
    
    double meanVariance = 0;
    unsigned nvarpoints = 0;
    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
        if(!isnan(compfluxvector->getvariance(i))) {
            meanVariance += compfluxvector->getvariance(i);
            nvarpoints++;
        }
    }
    double CompMinPeakDepth = MinPeakDepth*sqrt(meanVariance/(double)nvarpoints);
    
    compLines.detectSpectralFeatures(DetectionThreshold,CompLocalMaxFilterWidth,CompMinPeakDepth);
    
    unsigned line = 0;
    float flinesigma[MAXREFWAVELENGTHSPERORDER];
    
    for(unsigned feature=0;feature<compLines.getNFeatures();feature++) {
        operaSpectralFeature *currentFeature = compLines.getSpectralFeature(feature);
        double *center = currentFeature->getGaussianFit()->getCenterVector();
        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector();
        for(unsigned l=0; l<currentFeature->getnLines(); l++) {
            linecenter[line] = center[l];
            linecenterError[line] = centerError[l];
            lineflux[line] = amplitude[l];
            linesigma[line] = sigma[l];
            flinesigma[line] = (float)linesigma[line];
#ifdef PRINT_DEBUG
            //    cout << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << " " << " " << currentFeature->getnLines() << " " << currentFeature->getGaussianFit()->getGaussianChisqr() << endl;
#endif
            line++;
        }
    }
    
    nlinesinorder = line;
    if (nlinesinorder > 0) {
        linewidth = (double)operaArrayMedian(nlinesinorder,flinesigma);
        linewidth_err = (double)operaArrayMedianSigma(nlinesinorder,flinesigma,(float)(linewidth));
        
        line = 0;
        for(unsigned feature=0;feature<compLines.getNFeatures();feature++) {
            operaSpectralFeature *currentFeature = compLines.getSpectralFeature(feature);
            double *center = currentFeature->getGaussianFit()->getCenterVector();
            double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
            double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
            double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector();
            for(unsigned l=0; l<currentFeature->getnLines(); l++) {
                if(sigma[l] > linewidth - (double)nsigclip*(linewidth_err) && sigma[l] < linewidth + (double)nsigclip*(linewidth_err)) {
                    linecenter[line] = center[l];
                    linecenterError[line] = centerError[l];
                    lineflux[line] = amplitude[l];
                    linesigma[line] = sigma[l];
                    flinesigma[line] = (float)linesigma[line];
#ifdef PRINT_DEBUG
                    //cout << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << " " << " " << currentFeature->getnLines() << " " << currentFeature->getGaussianFit()->getGaussianChisqr() << endl;
#endif
                    line++;
                }
            }
        }
        
        nlinesinorder = line;
        if (nlinesinorder > 0) {
            linewidth = (double)operaArrayMedian(nlinesinorder,flinesigma);
            linewidth_err = (double)operaArrayMedianSigma(nlinesinorder,flinesigma,(float)(linewidth));
            compSpectrum->setHasWavelength(true);
        }
    }

    if(!emissionSpectrum) {
        // transform spectrum back to absorption
        for(unsigned indexElem=0;indexElem < compSpectrum->getnSpectralElements(); indexElem++) {
            compSpectrum->setFlux(saveSpectrum[indexElem],indexElem);
        }

        for(unsigned i=0; i<nlinesinorder; i++) {
            double emissionFlux = lineflux[i];
            lineflux[i] = 1.0 - emissionFlux;
        }
        delete[] saveSpectrum;
    }

    doubleValue_t Resolution;
    
    Resolution.value = wavelength->getcentralWavelength()/linewidth;
    Resolution.error = linewidth_err * wavelength->getcentralWavelength()/(linewidth*linewidth);
    
    wavelength->setSpectralResolution(Resolution);

    return nlinesinorder;
}

double calculateDeltaRadialVelocityInKPS(double telluricWL, double observedWL, double spectralResolution) {
    double linewidth = telluricWL/spectralResolution;
    double radialVelocity = (telluricWL - observedWL)*SPEED_OF_LIGHT_KMS/telluricWL;
    return radialVelocity;
}


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

#include "libraries/operaLibCommon.h"       // for SPEED_OF_LIGHT_M
#include "libraries/operaStats.h"       // for operaCrossCorrelation
#include "libraries/operaFFT.h"         // for operaXCorrelation

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

//convolveSpectrumWithGaussian(npointsInShortTelluricSpectrum,wl,trasmission,convolvedTelluricSpectrum,actuallinewidth);


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


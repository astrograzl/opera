#ifndef OPERASPECTRALTOOLS_H
#define OPERASPECTRALTOOLS_H

/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Library name: operaSpectralTools
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

#include <stdarg.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include "libraries/operaSpectralElements.h"    // for operaSpectralElements
#include "libraries/ladfit.h"						// for ladfit

#define MAXNUMBEROFWLRANGES 1000
#define MINNUMBEROFPOINTSINSIDEBIN 5


/*! 
 * \brief general library routines.
 * \file operaSpectralTools.h
 * \ingroup libraries
 */

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
void calculateXCorrWithGaussian(unsigned np, double *wavelength, double *flux, double *outputXcorr, double sigma);

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux);

/*
 * normalizeSpectrum(unsigned nLines, double *lineflux)
 * \brief This function normalize an input vector of fluxes and respective variances by the maximum value in the array
 * \param unsigned nLines
 * \param double *lineflux
 * \param double *linevariance 
 * \return void
 */
void normalizeSpectrum(unsigned nLines, double *lineflux, double *linevariance);

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
void convolveSpectrumWithGaussian(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double sigma);


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
void convolveSpectrumWithGaussianByResolution(unsigned np, double *wavelength, double *flux, double *convolvedSpectrum, double spectralResolution);

/*
 * double convertVacuumToAirWavelength(double vac_wl)
 * \brief This function converts from vacuum to air wavelength using the IAU standard for conversion
 * \brief  from air to vacuum wavelengths as given in Morton (1991, ApJS, 77, 119)
 * \param double vac_wl (in Angstrom)
 * \return double air_wl (in Angstrom)
 */
double convertVacuumToAirWavelength(double vac_wl);

/*
 * double planck(double T, double wl)
 * \brief This function returns the spectral radiance of a black body in photons/(s m^2 dlambda).
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double PlanckFunction(double T, double wl);

/*
 * double calculateBlackBodyVFlux(double Temperature)
 * \brief This function calculates the flux (in ph/(m^2 s m)) for a Black Body
 * \param double Temperature (in Kelvin)
 * \return double flux (in ph/(m^2 s m))
 */
double calculateBlackBodyVFlux(double Temperature);

/*
 * double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T)
 * \brief This function returns the integrated spectral flux of a black body in photons/(s m^2) for a given
 * \brief wavelength range using the simpson method.
 * \param Temperature is a double input that represents the temperature of the black body in Kelvins
 * \param Wavelength is a double input that represents the wavelength at which the black body is observed in nanometers
 * \return double value for the spectral radiance of the black body in ph/(s m^2 dlambda)
 */
double IntegrateSpectralElementOfBlackBody(double wl0, double wlf, double T);

double getFactorToMatchFluxesBetweenElements(operaSpectralElements *refElements,operaSpectralElements *elementsToMatch, double delta_wl);

void calculateUniformSample(unsigned np,float *wl,float *flux, unsigned npout, float *uniform_wl, float *uniform_flux);

float getFluxAtWavelength(unsigned np,float *wl,float *flux,float wavelengthForNormalization);

unsigned readContinuumWavelengthMask(string wavelength_mask, double *wl0, double *wlf);

unsigned getSpectrumWithinWLRange(operaSpectralElements *inputSpectrum, double wl0, double wlf, double *outputFlux, double *outputWavelength);

bool getOverlappingWLRange(operaSpectralElements *refElements, operaSpectralElements *elementsToMatch, double &wl0, double &wlf);

unsigned detectSpectralLinesInSpectralOrder(operaSpectralOrder *spectralOrder, double *linecenter, double *linecenterError, double *lineflux, double *linesigma, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip, double spectralResolution,bool emissionSpectrum);

double calculateDeltaRadialVelocityInKPS(double telluricWL, double observedWL, double spectralResolution);

#endif

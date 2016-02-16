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
#include "libraries/operaSpectralOrder.h"    // for operaSpectralElements
#include "libraries/ladfit.h"						// for ladfit

#define MAXNUMBEROFWLRANGES 1000
#define MINNUMBEROFPOINTSINSIDEBIN 5


/*! 
 * \brief general library routines.
 * \file operaSpectralTools.h
 * \ingroup libraries
 */

using namespace std;

class operaSpectrum {
private:
	operaVector wavelength;
	operaFluxVector intensity;
public:
	operaSpectrum() { }
	operaSpectrum(unsigned presize) : wavelength(presize), intensity(presize) { }
	unsigned int size() const { return wavelength.size(); }
	bool empty() const { return wavelength.empty(); }
	void insert(double wl, double flux) { wavelength.insert(wl); intensity.insert(flux, 0); }
	void insert(double wl, double flux, double variance) { wavelength.insert(wl); intensity.insert(flux, variance); }
	void reverse() { wavelength.reverse(); intensity.reverse(); }
	void resize(unsigned newsize) { wavelength.resize(newsize); intensity.resize(newsize); }
	double* wavelength_ptr(unsigned offset = 0) { return wavelength.datapointer()+offset; }
	double* flux_ptr(unsigned offset = 0) { return intensity.getfluxpointer()+offset; }
	double* variance_ptr(unsigned offset = 0) { return intensity.getvariancepointer()+offset; }
	const double* wavelength_ptr(unsigned offset = 0) const { return wavelength.datapointer()+offset; }
	const double* flux_ptr(unsigned offset = 0) const { return intensity.getfluxpointer()+offset; }
	const double* variance_ptr(unsigned offset = 0) const { return intensity.getvariancepointer()+offset; }
	double firstwl() const { return wavelength.first(); }
	double midwl() const { return wavelength[wavelength.size()/2]; }
	double lastwl() const { return wavelength.last(); }
	double getwavelength(unsigned i) const { return wavelength[i]; }
	double getflux(unsigned i) const { return intensity.getflux(i); }
	double getvariance(unsigned i) const { return intensity.getvariance(i); }
	void Sort() { operaIndexMap indexmap = wavelength.indexsort(); wavelength.reorder(indexmap); intensity.reorder(indexmap); }
    const operaVector& wavelengthvector() const { return wavelength; }
    const operaFluxVector& getintensity() const { return intensity; }
};

class operaWavelengthRange {
private:
	double wl0;
	double wlf;
public:
	operaWavelengthRange(double start, double end) : wl0(start), wlf(end) { }
    bool contains(double wavelength) { return wl0 <= wavelength && wavelength <= wlf; }
    double getwl0() const { return wl0; }
    double getwlf() const { return wlf; }
};

class operaWavelengthRanges {
private:
	vector <operaWavelengthRange> ranges;
public:
    void addWavelengthRange(double start, double end) { ranges.push_back(operaWavelengthRange(start, end)); }
    double wl0(unsigned i) const { return ranges[i].getwl0(); }
    double wlf(unsigned i) const { return ranges[i].getwlf(); }
    bool contains(double wavelength) { for(unsigned i = 0; i < ranges.size(); i++) if(ranges[i].contains(wavelength)) return true; return false; }
    bool contains(double wavelength, unsigned i) { if(ranges[i].contains(wavelength)) return true; return false; }
    operaWavelengthRange getrange(unsigned i) { return ranges[i];}
    unsigned int size() const { return ranges.size(); }
};

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
operaWavelengthRanges readContinuumWavelengthMask(string wavelength_mask);
operaWavelengthRanges readIndexedWavelengthMask(string wavelength_mask);


unsigned getSpectrumWithinWLRange(operaSpectralElements *inputSpectrum, double wl0, double wlf, double *outputFlux, double *outputWavelength);

bool getOverlappingWLRange(operaSpectralElements *refElements, operaSpectralElements *elementsToMatch, double &wl0, double &wlf);

unsigned detectSpectralLinesInSpectralOrder(operaSpectralOrder *spectralOrder, double *linecenter, double *linecenterError, double *lineflux, double *linesigma, double LocalMaxFilterWidth, double MinPeakDepth, double DetectionThreshold, double nsigclip, double spectralResolution,bool emissionSpectrum);

operaWavelengthRanges getWavelengthMaskAroundLines(const operaSpectrum sourceLines, double spectralResolution, double nsig);

double calculateDeltaRadialVelocityInKPS(double telluricWL, double observedWL);

/*!
 * \brief Resamples a spectrum to a new set of wavelengths.
 * \param inputWavelength The wavelengths of the input spectrum
 * \param inputFlux The flux of the input spectrum
 * \param outputWavelength The wavelengths to resample to
 * \return The flux of the fit spectrum at each point along outputWavelength
 */
operaVector fitSpectrum(const operaVector& inputWavelength, const operaVector& inputFlux, const operaVector& outputWavelength);

/*!
 * \brief Resamples a spectrum to a new set of wavelengths.
 * \param inputSpectrum Input operaSpectrum
 * \param outputWavelength The wavelengths to resample to
 * \return The fit spectrum at each point along outputWavelength
 */
operaSpectrum fitSpectrum(operaSpectrum inputSpectrum, const operaVector& outputWavelength);

/*!
 * \brief Returns the subset of an input spectrum, which lies within the mask ranges
 * \param string The wavelength ranges to return output spectrum
 * \param inputSpectrum Input operaSpectrum
 * \return Output operaSpectrum
 */
operaSpectrum getSpectrumWithinMask(string wavelength_mask, operaSpectrum inputSpectrum);

/*!
 * \brief Returns the subset of an input spectrum, which lies within a given opera wavelength range
 * \param operaWavelengthRange The wavelength range to return output spectrum
 * \param inputSpectrum Input operaSpectrum
 * \return Output operaSpectrum
 */
operaSpectrum getSpectrumWithinRange(operaWavelengthRange wlrange, operaSpectrum inputSpectrum);

/*!
 * \brief Returns a masked spectrum, where it removes points at +/- nsig x resolution element away from line centers, from list of lines provided
 * \param operaSpectrum The input spectrum
 * \param operaSpectrum telluricLines
 * \param double spectral resolution
 * \param double nsig -- define size of mask around each line in units of resolution element.
 * \return Output operaSpectrum
 */
operaSpectrum maskSpectrumAroundLines(const operaSpectrum inputSpectrum, const operaSpectrum telluricLines, double spectralResolution, double nsig);

operaVector convolveSpectrum(operaSpectrum inputSpectrum, double spectralResolution);

operaFluxVector fitSpectrumToPolynomial(const operaVector& inputWavelength, const operaFluxVector& inputFlux, unsigned order, double portionToUseForFit);

unsigned findClosestInSortedRange(const operaVector& vector, double target, unsigned startindex, unsigned endindex);

#endif

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralElements
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"		// for operaSpectralOrder_t
#include "libraries/operaSpectralElements.h"
#include "libraries/operaFluxVector.h"          // for operaFluxVectors

#include "libraries/operaMath.h"
#include "libraries/operaFit.h"

/*!
 * operaSpectralElements
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates a vector of spectral elements.
 * \file operaSpectralElements.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectralElements
 * \brief A spectral element (SE) consists of a spectrally resolved 
 * \brief subdivision of a spectral order taken in the dispersion direction. 
 * \brief The minimal width of the spectral element is defined by the resolution 
 * \brief of the spectrograph (RE). The actual width could be defined as any size 
 * \brief smaller than the order size and as long as RE < RS, where RE is the 
 * \brief resolution of the element.
 * \return none
 */
/*
 * Constructor
 */
operaSpectralElements::operaSpectralElements() :
nSpectralElements(0),
maxnSpectralElements(0),
elementHeight(0),
SpectrumType(None),
fluxvector(NULL),
XCorrelation(NULL),
photoCenterX(NULL),
photoCenterY(NULL),
distd(NULL),
wavelength(NULL),
fluxSNR(NULL),
tell(NULL),
rvel(NULL),
rawFlux(NULL),
normalizedFlux(NULL),
fcalFlux(NULL),
hasRawFlux(false),
hasStandardFlux(false),
hasOptimalFlux(false),
hasOperaOptimalFlux(false),
hasExtendedBeamFlux(false),
hasXCorrelation(false),
hasWavelength(false),    
hasDistance(false),   
hasFluxSNR(false)    
{
	
}

operaSpectralElements::operaSpectralElements(unsigned nElements):
nSpectralElements(0),
maxnSpectralElements(0),
elementHeight(0),
SpectrumType(None),
fluxvector(NULL),
XCorrelation(NULL),
photoCenterX(NULL),
photoCenterY(NULL),
distd(NULL),
wavelength(NULL),
fluxSNR(NULL),
tell(NULL),
rvel(NULL),
rawFlux(NULL),
normalizedFlux(NULL),
fcalFlux(NULL),
hasRawFlux(false),
hasStandardFlux(false),
hasOptimalFlux(false),
hasOperaOptimalFlux(false),
hasExtendedBeamFlux(false),
hasXCorrelation(false),
hasWavelength(false),    
hasDistance(false),   
hasFluxSNR(false)    
{
	Createvectors(nElements);
}

operaSpectralElements::operaSpectralElements(unsigned nElements, operaSpectralOrder_t format, bool extended):
nSpectralElements(0),
maxnSpectralElements(0),
elementHeight(0),
SpectrumType(None),
fluxvector(NULL),
XCorrelation(NULL),
photoCenterX(NULL),
photoCenterY(NULL),
distd(NULL),
wavelength(NULL),
fluxSNR(NULL),
tell(NULL),
rvel(NULL),
rawFlux(NULL),
normalizedFlux(NULL),
fcalFlux(NULL),
hasRawFlux(false),
hasStandardFlux(false),
hasOptimalFlux(false),
hasOperaOptimalFlux(false),
hasExtendedBeamFlux(false),
hasXCorrelation(false),
hasWavelength(false),    
hasDistance(false),   
hasFluxSNR(false)    
{
	Createvectors(nElements, format, extended);
}

/*
 * Destructor
 */
operaSpectralElements::~operaSpectralElements() {
	Deletevectors();
}

/*
 * getters/setters
 */

void operaSpectralElements::Createvectors(unsigned nElements, bool extended) {
	if (nElements == 0) {
		throw operaException("operaSpectralElements: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (maxnSpectralElements != 0) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	maxnSpectralElements = nElements;
	nSpectralElements = nElements;
	SpectrumType = None;
	fluxvector = new operaFluxVector(nElements);   
	photoCenterX = (double *)malloc(nElements*sizeof(double));
	photoCenterY = (double *)malloc(nElements*sizeof(double));
	distd = (double *)malloc(nElements*sizeof(double));	
	XCorrelation = (double *)malloc(nElements*sizeof(double));
	fluxSNR = (double *)malloc(nElements*sizeof(double));
	wavelength = (double *)malloc(nElements*sizeof(double));
    if (extended) {
		CreateExtendedvectors(nElements);
	}
}

void operaSpectralElements::CreateExtendedvectors(unsigned nElements) {
    hasExtendedBeamFlux = true;
	tell = (double *)malloc(nElements*sizeof(double));
	rvel = (double *)malloc(nElements*sizeof(double));
	
	rawFlux = new operaFluxVector(nElements);
	normalizedFlux = new operaFluxVector(nElements);
	fcalFlux = new operaFluxVector(nElements);

}

void operaSpectralElements::Createvectors(unsigned nElements, operaSpectralOrder_t format, bool extended) {
	Createvectors(nElements, extended);
	SpectrumType = format;
	switch (format) {
		case RawSpectrum:
		case RawBeamSpectrum:
			setHasRawSpectrum(true);
			break;
		case StandardSpectrum:
		case StandardBeamSpectrum:
			setHasStandardSpectrum(true);
			break;
		case OptimalSpectrum:
		case OptimalBeamSpectrum:
			setHasOptimalSpectrum(true);
			break;
		case OperaOptimalSpectrum:
		case OperaOptimalBeamSpectrum:
			setHasOperaOptimalSpectrum(true);
			break;
		default:
			break;
	}
}

void operaSpectralElements::Resizevector(unsigned nElements, operaSpectralOrder_t format) {
	if (nElements == 0) {
		throw operaException("operaSpectralElements: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (maxnSpectralElements < nElements) {
		Deletevectors();
		Createvectors(nElements, format);	// sets maxnSpectralElements...
	} else {
		nSpectralElements = nElements;
		SpectrumType = format;
	}
}

void operaSpectralElements::Deletevectors(void) {
	if (maxnSpectralElements == 0) {
		return;
	}
	if (photoCenterX)
		free(photoCenterX);
	if (photoCenterY)
		free(photoCenterY);
	if (distd)
		free(distd);	
	if (wavelength)
		free(wavelength);
	if (fluxSNR)
		free(fluxSNR);
	if (XCorrelation)
		free(XCorrelation);    
	if (fluxvector)
		delete fluxvector;
	if (tell)
		free(tell);    
	if (rvel)
		free(rvel);    
	if (rawFlux)
		delete rawFlux;
	if (normalizedFlux)
		delete normalizedFlux;
	if (fcalFlux)
		delete fcalFlux;
	
	fluxvector = NULL;
	photoCenterX = NULL;
	photoCenterY = NULL;
	distd = NULL;
    XCorrelation = NULL;
	wavelength = NULL;
	fluxSNR = NULL;
	tell = NULL;
	rvel = NULL;
	rawFlux = NULL;
	normalizedFlux = NULL;
	fcalFlux = NULL;
	nSpectralElements = 0;
	maxnSpectralElements = 0;
	SpectrumType = None;
    
	setHasRawSpectrum(false);
	setHasStandardSpectrum(false);
	setHasOptimalSpectrum(false);
	setHasOperaOptimalSpectrum(false);
	setHasExtendedBeamFlux(false);
	setHasWavelength(false);   
	setHasDistance(false);
	setHasFluxSNR(false);
    setHasXCorrelation(false);    
}

void operaSpectralElements::setnSpectralElements(unsigned nElem) { 
	if (nElem == 0) {
		throw operaException("operaSpectralElements: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    if (nElem > maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	nSpectralElements = nElem;
}

unsigned operaSpectralElements::getnSpectralElements(void) { 
	return nSpectralElements;
}

void operaSpectralElements::setelementHeight(double Height) { 
	elementHeight = Height;
}

double operaSpectralElements::getelementHeight(void) { 
	return elementHeight;
}

operaSpectralOrder_t operaSpectralElements::getSpectrumType(void) { 
	return SpectrumType;
}

void operaSpectralElements::setSpectrumType(operaSpectralOrder_t format) { 
	SpectrumType = format;
}

/*
 * clone spectralelements
 */
void operaSpectralElements::setSpectralElements(operaSpectralElements &SpectralElements) { 

    if (SpectralElements.nSpectralElements == 0) {
		throw operaException("operaSpectralElements: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    if (SpectralElements.maxnSpectralElements == 0) {
		throw operaException("operaSpectralElements: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (nSpectralElements < SpectralElements.nSpectralElements) {
		throw operaException("operaSpectralLines: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	nSpectralElements = SpectralElements.nSpectralElements;
	maxnSpectralElements = SpectralElements.maxnSpectralElements;
	elementHeight = SpectralElements.elementHeight;
	SpectrumType = SpectralElements.SpectrumType;

	hasRawFlux = SpectralElements.hasRawFlux;
	hasStandardFlux = SpectralElements.hasStandardFlux;
	hasOptimalFlux = SpectralElements.hasOptimalFlux;
	hasOperaOptimalFlux = SpectralElements.hasOperaOptimalFlux;
	hasExtendedBeamFlux = SpectralElements.hasExtendedBeamFlux;
	hasXCorrelation = SpectralElements.hasXCorrelation;
	hasWavelength = SpectralElements.hasWavelength;
	hasDistance = SpectralElements.hasDistance;
	hasFluxSNR = SpectralElements.hasFluxSNR;

	if (fluxvector == NULL && SpectralElements.fluxvector != NULL)
		fluxvector = new operaFluxVector(maxnSpectralElements); 	
	
	if (photoCenterX == NULL && SpectralElements.photoCenterX != NULL)
		photoCenterX = (double *)malloc(maxnSpectralElements*sizeof(double));
	
	if (photoCenterY == NULL && SpectralElements.photoCenterY != NULL)
		photoCenterY = (double *)malloc(maxnSpectralElements*sizeof(double));
	
	if (XCorrelation == NULL && SpectralElements.XCorrelation != NULL)
		XCorrelation = (double *)malloc(maxnSpectralElements*sizeof(double)); 
	
	if (distd == NULL && SpectralElements.distd != NULL)
		distd = (double *)malloc(maxnSpectralElements*sizeof(double));	
	
	if (wavelength == NULL && SpectralElements.wavelength != NULL)
		wavelength = (double *)malloc(maxnSpectralElements*sizeof(double));	
	
	if (fluxSNR == NULL && SpectralElements.fluxSNR != NULL)
		fluxSNR = (double *)malloc(maxnSpectralElements*sizeof(double));	
	
	if (fluxSNR == NULL) {
		throw operaException("operaSpectralElements: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	if (fluxvector != NULL) {
		operaFluxVector *fv = SpectralElements.fluxvector;
		fluxvector->setVectors(fv->getfluxes(), fv->getvariances()); // copies...
	}

	for (unsigned i=0; i<nSpectralElements; i++) {
		if (photoCenterX != NULL && photoCenterY != NULL)
			setphotoCenter(SpectralElements.getphotoCenterX(i), SpectralElements.getphotoCenterY(i), i);
		if (XCorrelation != NULL)
			setXCorrelation(SpectralElements.getXCorrelation(i), i);
		if (distd != NULL)
			setdistd(SpectralElements.getdistd(i), i);
		if (wavelength != NULL)
			setwavelength(SpectralElements.getwavelength(i), i);
		if (fluxSNR != NULL)
			setFluxSNR(SpectralElements.getFluxSNR(i), i);
	}
}

/*
 * flux
 */
operaFluxVector *operaSpectralElements::getFluxVector(void) {
	return fluxvector;
}

// Note that fluxvector is protected against length violations...
double operaSpectralElements::getFlux(unsigned indexElem) {
	return fluxvector->getflux(indexElem);
}

void operaSpectralElements::setFlux(double Flux, unsigned indexElem) {
	fluxvector->setflux(Flux, indexElem);
}
void operaSpectralElements::setFluxVector(operaFluxVector *FluxVector) {
	fluxvector->setVector(*FluxVector);	// copies
}
double operaSpectralElements::getFluxVariance(unsigned indexElem) {
	return fluxvector->getvariance(indexElem);
}

void operaSpectralElements::setFluxVariance(double FluxVariance, unsigned indexElem) {
	fluxvector->setvariance(FluxVariance, indexElem);
}

double operaSpectralElements::getFluxSNR(unsigned indexElem) {
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (!hasFluxSNR) {	// hasFluxSNR is the case where SNR is read from a .sn file
		if (!isnan(fluxvector->geterror(indexElem)) && fluxvector->geterror(indexElem) != 0.0) {
			fluxSNR[indexElem] = fabs(fluxvector->getflux(indexElem)/fluxvector->geterror(indexElem));
		} else {
			fluxSNR[indexElem] = 1.0;
		}
	}
	return fluxSNR[indexElem];
}

void operaSpectralElements::setFluxSNR(double FluxSNR, unsigned indexElem) {
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	fluxSNR[indexElem] = FluxSNR;
}

/*
 * x,y coords
 */
double operaSpectralElements::getphotoCenterX(unsigned indexElem) {							// (flt in counts)
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return photoCenterX[indexElem];
}
double operaSpectralElements::getphotoCenterY(unsigned indexElem) {							// (flt in counts)
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return photoCenterY[indexElem];
}
void operaSpectralElements::setphotoCenter(double x, double y, unsigned indexElem) {						// (flt in counts)
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	photoCenterX[indexElem] = x;
	photoCenterY[indexElem] = y;	
}

/*
 * distance
 */
double operaSpectralElements::getdistd(unsigned indexElem) {							// (flt in counts)
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return distd[indexElem];
}
void operaSpectralElements::setdistd(double Distd, unsigned indexElem) {
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	distd[indexElem] = Distd;
}

/*
 * wavelength
 */
double operaSpectralElements::getwavelength(unsigned indexElem) {
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return wavelength[indexElem];
}
void operaSpectralElements::setwavelength(double Wavelength, unsigned indexElem) {
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	wavelength[indexElem] = Wavelength;
}

void operaSpectralElements::setwavelengthsFromCalibration(operaWavelength *Wavelength) {
    for(unsigned indexElem=0;indexElem<getnSpectralElements();indexElem++) {
        wavelength[indexElem] = Wavelength->evaluateWavelength(distd[indexElem]);
    }
    setHasWavelength(true);
}

/*
 * x-correlation
 */
double operaSpectralElements::getXCorrelation(unsigned indexElem){
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return XCorrelation[indexElem];
}
void operaSpectralElements::setXCorrelation(double Xcorr, unsigned indexElem){
#ifdef RANGE_CHECK
    if (indexElem >= maxnSpectralElements) {
		throw operaException("operaSpectralElements: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	XCorrelation[indexElem] = Xcorr;
}
   
double operaSpectralElements::gettell(unsigned indexElem) { 
	return tell[indexElem];
}
void operaSpectralElements::settell(double value, unsigned indexElem) { 
	tell[indexElem] = value;
}    
void operaSpectralElements::copyTOtell(void) { 
	for (unsigned i=0; i<nSpectralElements; i++) {
		tell[i] = wavelength[i];
	}
}    
void operaSpectralElements::copyFROMtell(void) { 
	for (unsigned i=0; i<nSpectralElements; i++) {
		wavelength[i] = tell[i];
	}
}    
double operaSpectralElements::getrvel(unsigned indexElem) { 
	return rvel[indexElem];
}
void operaSpectralElements::setrvel(double value, unsigned indexElem) { 
	rvel[indexElem] = value;
}    
void operaSpectralElements::copyTOrvel(void) { 
	for (unsigned i=0; i<nSpectralElements; i++) {
		rvel[i] = wavelength[i];
	}
}    
void operaSpectralElements::copyFROMrvel(void) { 
	for (unsigned i=0; i<nSpectralElements; i++) {
		wavelength[i] = rvel[i];
	}
}    
double operaSpectralElements::getnormalizedFlux(unsigned indexElem) { 
	return normalizedFlux->getflux(indexElem);
}
void operaSpectralElements::setnormalizedFlux(double value, unsigned indexElem) { 
	normalizedFlux->setflux(value, indexElem);
}
void operaSpectralElements::setnormalizedFluxVariance(double value, unsigned indexElem) {
	normalizedFlux->setvariance(value, indexElem);
}
void operaSpectralElements::copyTOnormalizedFlux(void) { 
    *normalizedFlux = *fluxvector;
}
void operaSpectralElements::copyFROMnormalizedFlux(void) { 
    *fluxvector = *normalizedFlux;
}
double operaSpectralElements::getfcalFlux(unsigned indexElem) { 
	return fcalFlux->getflux(indexElem);
}
void operaSpectralElements::setfcalFlux(double value, unsigned indexElem) {
	fcalFlux->setflux(value, indexElem);
}
void operaSpectralElements::setfcalFluxVariance(double value, unsigned indexElem) {
	fcalFlux->setvariance(value, indexElem);
}
void operaSpectralElements::copyTOfcalFlux(void) {
    *fcalFlux = *fluxvector;
}
void operaSpectralElements::copyFROMfcalFlux(void) { 
    *fluxvector = *fcalFlux;
}
double operaSpectralElements::getrawFlux(unsigned indexElem) {
	return rawFlux->getflux(indexElem);
}
void operaSpectralElements::setrawFlux(double value, unsigned indexElem) {
	rawFlux->setflux(value, indexElem);
}
void operaSpectralElements::setrawFluxVariance(double value, unsigned indexElem) {
	rawFlux->setvariance(value, indexElem);
}

void operaSpectralElements::copyTOrawFlux(void) {
    *rawFlux = *fluxvector;
}
void operaSpectralElements::copyFROMrawFlux(void) {
    *fluxvector = *rawFlux;
}

/*
 * Other Methods
 */
double operaSpectralElements::getFluxSum(void) {
    double fluxsum = 0;
    
    for(unsigned indexElem=0;indexElem<getnSpectralElements();indexElem++) {
        if(!isnan(fluxvector->getflux(indexElem))) {
            fluxsum += fluxvector->getflux(indexElem);
        }
    }
    
    return fluxsum;
}


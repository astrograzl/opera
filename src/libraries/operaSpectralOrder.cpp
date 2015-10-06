/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralOrder
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

#include <iostream>
#include <fstream>

#include "globaldefines.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaGeometry.h"
#include "libraries/operaWavelength.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaSpectralLines.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaInstrumentProfile.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/PixelSet.h"
#include "libraries/Gaussian.h"

#define DEBUG false
#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

#ifndef SATURATIONLIMIT
#define SATURATIONLIMIT 65535  // this should be retrieved from the config/param file
#endif
#ifndef MAXNUMBEROFLINES
#define MAXNUMBEROFLINES 10000
#endif
#ifndef NSIGMACUT
#define NSIGMACUT 3
#endif
#ifndef NMORESIGMASTOSTARTOPTIMALEXTRACTION
#define NMORESIGMASTOSTARTOPTIMALEXTRACTION 1
#endif
#ifndef DEFAULT_SNR_SMOOTHING
#define DEFAULT_SNR_SMOOTHING 5
#endif

#include "libraries/operaLib.h"     // for itos
#include "libraries/operaStats.h"
#include "libraries/operaFit.h"
#include "libraries/operaCCD.h"
#include "libraries/operaMath.h"
#include "libraries/ladfit.h"
#include "libraries/VLArray.h"

/*!
 * operaSpectralOrder
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the FITS image.
 * \file operaSpectralOrder.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \brief operaSpectralOrder
 * \details A spectral order (SO) consists of a data set containing 
 * \details the information concerned with a full spectral order.
 */

/*
 * Constructors
 */
operaSpectralOrder::operaSpectralOrder() : 
orderNumber(0),
SpectrumType(None),
SNR(0.0),
SpectralElements(NULL),
SkyElements(NULL),
Geometry(NULL),
Wavelength(NULL),
InstrumentProfile(NULL),
SpectralLines(NULL),
Polarimetry(NULL),
SpectralEnergyDistribution(NULL),
numberOfBeams(0),
hasSpectralElements(false), 
hasSkyElements(false),
hasGeometry(false), 
hasWavelength(false), 
hasInstrumentProfile(false), 
hasSpectralLines(false), 
hasExtractionApertures(false),
hasPolarimetry(false),
hasSNR(false),
hasCenterSNROnly(false),
hasSpectralEnergyDistribution(false),
hasWavelengthRange(false)
{
    tiltInDegrees.value = 0.0;
    tiltInDegrees.error = 0.0;
    
    for (unsigned backgroundIndex = 0 ; backgroundIndex < LEFTANDRIGHT ; backgroundIndex++) {
        BackgroundElements[backgroundIndex] = NULL;
        BackgroundApertures[backgroundIndex] = NULL;
    }
    
    for (unsigned beam = 0 ; beam < MAXNUMBEROFBEAMS ; beam++) {
        BeamElements[beam] = NULL;
        BeamProfiles[beam] = NULL;
        ExtractionApertures[beam] = NULL;
    }
}

operaSpectralOrder::operaSpectralOrder(unsigned order) : 
SpectrumType(None),
SNR(0.0),
SpectralElements(NULL),
SkyElements(NULL),
Geometry(NULL),
Wavelength(NULL),
InstrumentProfile(NULL),
SpectralLines(NULL),
Polarimetry(NULL),
SpectralEnergyDistribution(NULL),
numberOfBeams(0),
hasSpectralElements(false), 
hasSkyElements(false),
hasGeometry(false), 
hasWavelength(false), 
hasInstrumentProfile(false), 
hasSpectralLines(false), 
hasExtractionApertures(false),
hasPolarimetry(false),
hasSNR(false),
hasCenterSNROnly(false),
hasSpectralEnergyDistribution(false),
hasWavelengthRange(false)
{ 	
	orderNumber = order;
    tiltInDegrees.value = 0.0;
    tiltInDegrees.error = 0.0;
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
        BackgroundElements[backgroundIndex] = NULL;//new operaSpectralElements();
        BackgroundApertures[backgroundIndex] = NULL;//new operaExtractionAperture();  
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
        BeamElements[beam] = NULL;//new operaSpectralElements();
        ExtractionApertures[beam] = NULL;//new operaExtractionAperture();
        BeamProfiles[beam] = NULL;//new operaInstrumentProfile();
		BeamSED[beam] = NULL;
    }
}

operaSpectralOrder::operaSpectralOrder(unsigned order, unsigned maxdatapoints, unsigned maxValues, unsigned nElements, operaSpectralOrder_t format)  : 
SpectrumType(None),
SNR(0.0),
SpectralElements(NULL),
SkyElements(NULL),
Geometry(NULL),
Wavelength(NULL),
InstrumentProfile(NULL),
SpectralLines(NULL),
Polarimetry(NULL),
SpectralEnergyDistribution(NULL),
numberOfBeams(0),
hasSpectralElements(false), 
hasSkyElements(false), 
hasGeometry(false), 
hasWavelength(false), 
hasInstrumentProfile(false), 
hasSpectralLines(false), 
hasExtractionApertures(false),
hasPolarimetry(false),
hasSNR(false),
hasCenterSNROnly(false),
hasSpectralEnergyDistribution(false),
hasWavelengthRange(false)
{
	if (maxdatapoints == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	// This is OK, it just means you don't need spectralElements...
	if (nElements == 0) {
		//throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (maxValues == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	orderNumber = order;
	SpectrumType = format;
    tiltInDegrees.value = 0.0;
    tiltInDegrees.error = 0.0;
    
	Geometry = new operaGeometry(maxdatapoints, maxValues);
	if (nElements) {
		SpectralElements = new operaSpectralElements(nElements, SpectrumType);
        SkyElements = new operaSpectralElements(nElements, SpectrumType);        
    }
	if (nElements)
		Polarimetry = new operaPolarimetry(nElements);
#if 0   
	if (nElements)
		for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
			BackgroundElements[backgroundIndex] = new operaSpectralElements(nElements, SpectrumType);
			BackgroundApertures[backgroundIndex] = new operaExtractionAperture();  
		}
    
	if (nElements)
		for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
			BeamElements[beam] = new operaSpectralElements(nElements, SpectrumType);
			ExtractionApertures[beam] = new operaExtractionAperture();  
			BeamProfiles[beam] = new operaInstrumentProfile();        
		}    
#else
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
        BackgroundElements[backgroundIndex] = NULL;
        BackgroundApertures[backgroundIndex] = NULL;  
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
        BeamElements[beam] = NULL;
        ExtractionApertures[beam] = NULL;
        BeamProfiles[beam] = NULL;
		BeamSED[beam] = NULL;
    }
#endif
}

/*
 * Destructor
 */
operaSpectralOrder::~operaSpectralOrder() {
	deleteAll();    
}

/*
 * Common Methods
 */

void operaSpectralOrder::deleteAll() {
	
	if (Geometry && hasGeometry) 
		delete Geometry;
	Geometry = NULL;
	if (Wavelength && hasWavelength) 
		delete Wavelength;
	Wavelength = NULL;
	if (SpectralElements && hasSpectralElements) 
		delete SpectralElements;
	SpectralElements = NULL;
	if (SkyElements && hasSkyElements) 
		delete SkyElements;
	SkyElements = NULL;    
	if (InstrumentProfile && hasInstrumentProfile) 
		delete InstrumentProfile;
	InstrumentProfile = NULL;
	if (SpectralLines && hasSpectralLines) 
		delete SpectralLines;
	SpectralLines = NULL;
    if (Polarimetry && hasPolarimetry) 
		delete Polarimetry;
	Polarimetry = NULL;
    if (SpectralEnergyDistribution && hasSpectralEnergyDistribution) 
		delete SpectralEnergyDistribution;
	SpectralEnergyDistribution = NULL;
	
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundElements[backgroundIndex]) 
			delete BackgroundElements[backgroundIndex];
		BackgroundElements[backgroundIndex] = NULL;
		
		if (BackgroundApertures[backgroundIndex])
			delete BackgroundApertures[backgroundIndex];
		BackgroundApertures[backgroundIndex] = NULL;
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (BeamElements[beam]) 
			delete BeamElements[beam];
		BeamElements[beam] = NULL;
		if (ExtractionApertures[beam]) 
			delete ExtractionApertures[beam];
		ExtractionApertures[beam] = NULL;
		if (BeamProfiles[beam]) 
			delete BeamProfiles[beam];
		BeamProfiles[beam] = NULL;
        if (BeamSED[beam]) 
            delete BeamSED[beam];
        BeamSED[beam] = NULL;        
    }    
	hasSpectralElements = false; 
	hasSkyElements = false; 
	hasGeometry = false; 
	hasWavelength = false; 
	hasInstrumentProfile = false; 
	hasSpectralLines = false; 
	hasExtractionApertures = false;
	hasPolarimetry = false;
	hasSNR = false;
	hasSpectralEnergyDistribution = false;
	hasWavelengthRange = false;
}

void operaSpectralOrder::deleteBeamProfiles(void) {
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (BeamProfiles[beam]) 
			delete BeamProfiles[beam];
		BeamProfiles[beam] = NULL;
    }       
}

void operaSpectralOrder::deleteApertures(void) {
    
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundApertures[backgroundIndex]) 
            delete BackgroundApertures[backgroundIndex];
        BackgroundApertures[backgroundIndex] = NULL;
    }    
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (ExtractionApertures[beam]) 
			delete ExtractionApertures[beam];
		ExtractionApertures[beam] = NULL;
    }        
}


void operaSpectralOrder::deleteBeamsAndBackgroundElements(void) {
    for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundElements[backgroundIndex]) 
			delete BackgroundElements[backgroundIndex];
		BackgroundElements[backgroundIndex] = NULL;
    }
    
    for(unsigned beam=0;beam<MAXNUMBEROFBEAMS; beam++) {
		if (BeamElements[beam]) 
			delete BeamElements[beam];
		BeamElements[beam] = NULL;
    }    
}

void operaSpectralOrder::createGeometry(unsigned maxdatapoints, unsigned maxValues) {
	if (Geometry) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	Geometry = new operaGeometry(maxdatapoints, maxValues);
}

void operaSpectralOrder::createWavelength(unsigned maxnumberofcoefficients) {
	if (Wavelength) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	Wavelength = new operaWavelength(maxnumberofcoefficients);    
}

void operaSpectralOrder::createSpectralElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType, bool extended) {
	if (SpectralElements != NULL) {
		// CU Aug 20, 2015 -- Swapped which line is commented out. If this is causing a crash, that's sign that something else is wrong.
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		//SpectralElements->Resizevectors(maxdatapoints, SpectralElements->getSpectrumType());
	} else {
		SpectralElements = new operaSpectralElements(maxdatapoints, SpectrumType, extended);
	}
}

void operaSpectralOrder::createSkyElements(unsigned maxdatapoints, operaSpectralOrder_t SpectrumType) {
	if (SkyElements) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	SkyElements = new operaSpectralElements(maxdatapoints, SpectrumType);
}

void operaSpectralOrder::createBeamsAndBackgrounds(unsigned nElements, unsigned nBeams, operaSpectralOrder_t format, bool extendedbeams) {
	if (nElements == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (nBeams == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	// recycle if we can...
	for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
		if (BackgroundElements[backgroundIndex] != NULL && BackgroundElements[backgroundIndex]->getnSpectralElements() >= nElements) {
			BackgroundElements[backgroundIndex]->setnSpectralElements(nElements);
			BackgroundElements[backgroundIndex]->setSpectrumType(format);
		} else {
			delete BackgroundElements[backgroundIndex];
			BackgroundElements[backgroundIndex] = new operaSpectralElements(nElements, format);
		}
    }
    
    for(unsigned beam=0;beam<nBeams; beam++) {
		if (BeamElements[beam] != NULL && BeamElements[beam]->getnSpectralElements() >= nElements) {
			BeamElements[beam]->setnSpectralElements(nElements);
			BeamElements[beam]->setSpectrumType(format);
		} else {
			delete BeamElements[beam];
			BeamElements[beam] = new operaSpectralElements(nElements, format, extendedbeams);
		}
    }    
}

void operaSpectralOrder::createSpectralEnergyDistributionElements(unsigned nElements) {
	if (nElements == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	// E. Martioli Mar 17 2015 -- this was causing some functions to crash. It should not be a problem if numberOfBeams=0.
	/*if (numberOfBeams == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}*/
    if(SpectralEnergyDistribution)
        delete SpectralEnergyDistribution;
    SpectralEnergyDistribution = new operaSpectralEnergyDistribution(nElements); // DT Jan 2013 added nElements
    SpectralEnergyDistribution->createFluxCalibrationElements(nElements);
    SpectralEnergyDistribution->createThroughputElements(nElements);    
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        if(BeamSED[beam])
            delete BeamSED[beam]; 
        
        BeamSED[beam] = new operaSpectralEnergyDistribution(nElements);// DT Jan 2013 added nElements
        BeamSED[beam]->createFluxCalibrationElements(nElements);
        BeamSED[beam]->createThroughputElements(nElements);         
    }
}


void operaSpectralOrder::deletePolarimetry(void) {
    if (Polarimetry) 
		delete Polarimetry;
	Polarimetry = NULL;
}

void operaSpectralOrder::createPolarimetry(unsigned nElements) {
	if (nElements == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    deletePolarimetry();
	Polarimetry = new operaPolarimetry(nElements);
}

unsigned operaSpectralOrder::getorder(void) const {
	return orderNumber;
}

void operaSpectralOrder::sethasSpectralElements(bool HasSpectralElements) {
	hasSpectralElements = HasSpectralElements;
}

void operaSpectralOrder::sethasSkyElements(bool HasSkyElements) {
	hasSkyElements = HasSkyElements;
}

void operaSpectralOrder::sethasGeometry(bool HasGeometry) {
	hasGeometry = HasGeometry;
}

void operaSpectralOrder::sethasWavelength(bool HasWavelength) {
	hasWavelength = HasWavelength;
}

void operaSpectralOrder::sethasInstrumentProfile(bool HasInstrumentProfile) {
	hasInstrumentProfile = HasInstrumentProfile;
}

void operaSpectralOrder::sethasSpectralLines(bool HasSpectralLines) {
	hasSpectralLines = HasSpectralLines;
}

void operaSpectralOrder::sethasExtractionApertures(bool HasExtractionApertures) {
	hasExtractionApertures = HasExtractionApertures;
}

void operaSpectralOrder::sethasPolarimetry(bool HasPolarimetry) {
	hasPolarimetry = HasPolarimetry;
}

void operaSpectralOrder::sethasSNR(bool HasSNR) {
	hasSNR = HasSNR;
}

void operaSpectralOrder::sethasCenterSNROnly(bool HasCenterSNROnly) {
	hasCenterSNROnly = HasCenterSNROnly;
}

void operaSpectralOrder::sethasSpectralEnergyDistribution(bool HasSpectralEnergyDistribution) {
	hasSpectralEnergyDistribution = HasSpectralEnergyDistribution;
}

bool operaSpectralOrder::gethasSpectralElements(void) const {
	return hasSpectralElements;
}

bool operaSpectralOrder::gethasSkyElements(void) const {
	return hasSkyElements;
}

bool operaSpectralOrder::gethasGeometry(void) const {
	return hasGeometry;
}

bool operaSpectralOrder::gethasWavelength(void) const {
	return hasWavelength;
}

bool operaSpectralOrder::gethasInstrumentProfile(void) const {
	return hasInstrumentProfile;
}

bool operaSpectralOrder::gethasSpectralLines(void) const {
	return hasSpectralLines;
}

bool operaSpectralOrder::gethasExtractionApertures(void) const {
	return hasExtractionApertures;
}

bool operaSpectralOrder::gethasPolarimetry(void) const {
	return hasPolarimetry;
}

bool operaSpectralOrder::gethasSNR(void) const {
	return hasSNR;
}

bool operaSpectralOrder::gethasCenterSNROnly(void) const {
	return hasCenterSNROnly;
}

bool operaSpectralOrder::gethasSpectralEnergyDistribution(void) const {
	return hasSpectralEnergyDistribution;
}

operaSpectralOrder_t operaSpectralOrder::getSpectrumType(void) const { 
	return SpectrumType;
}

void operaSpectralOrder::setSpectrumType(operaSpectralOrder_t format) { 
	SpectrumType = format;
}

/*
 * Common Methods
 */
operaSpectralElements *operaSpectralOrder::getSpectralElements() {
	return SpectralElements;
}

const operaSpectralElements *operaSpectralOrder::getSpectralElements() const {
	return SpectralElements;
}

operaSpectralElements *operaSpectralOrder::getSkyElements() {
	return SkyElements;
}

const operaSpectralElements *operaSpectralOrder::getSkyElements() const {
	return SkyElements;
}

operaGeometry *operaSpectralOrder::getGeometry() {
	return Geometry;
}

const operaGeometry *operaSpectralOrder::getGeometry() const {
	return Geometry;
}

operaWavelength *operaSpectralOrder::getWavelength() {
	return Wavelength;
}

const operaWavelength *operaSpectralOrder::getWavelength() const {
	return Wavelength;
}

operaInstrumentProfile *operaSpectralOrder::getInstrumentProfile() {
	return InstrumentProfile;
}

const operaInstrumentProfile *operaSpectralOrder::getInstrumentProfile() const {
	return InstrumentProfile;
}

operaSpectralLines *operaSpectralOrder::getSpectralLines() {
	return SpectralLines;
}

const operaSpectralLines *operaSpectralOrder::getSpectralLines() const {
	return SpectralLines;
}

operaPolarimetry *operaSpectralOrder::getPolarimetry() {
	return Polarimetry;
}

const operaPolarimetry *operaSpectralOrder::getPolarimetry() const {
	return Polarimetry;
}

operaSpectralEnergyDistribution *operaSpectralOrder::getSpectralEnergyDistribution(void) {
	return SpectralEnergyDistribution;
}

const operaSpectralEnergyDistribution *operaSpectralOrder::getSpectralEnergyDistribution(void) const {
	return SpectralEnergyDistribution;
}

double operaSpectralOrder::getCenterSNR(void) const {
	return SNR;
}

void operaSpectralOrder::setCenterSNR(double Snr) {
	SNR = Snr;
}

unsigned operaSpectralOrder::getnumberOfBeams(void) const {
	return numberOfBeams;
}   

void operaSpectralOrder::setnumberOfBeams(unsigned NumberOfBeams) {
#ifdef RANGE_CHECK
	if (NumberOfBeams >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	numberOfBeams = NumberOfBeams;
}

doubleValue_t operaSpectralOrder::getTiltInDegrees(void) const {
	return tiltInDegrees;
}

void operaSpectralOrder::setTiltInDegrees(doubleValue_t TiltInDegrees) {
	tiltInDegrees = TiltInDegrees;
}

void operaSpectralOrder::setTiltInDegrees(double tilt, double error) {
	tiltInDegrees.value = tilt;
	tiltInDegrees.error = error;    
}

double operaSpectralOrder::getTiltInDegreesValue(void) const {
	return tiltInDegrees.value;    
}

double operaSpectralOrder::getTiltInDegreesError(void) const {
	return tiltInDegrees.error;    
}

operaSpectralElements *operaSpectralOrder::getBeamElements(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamElements[beam];
}

const operaSpectralElements *operaSpectralOrder::getBeamElements(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamElements[beam];
}

operaInstrumentProfile *operaSpectralOrder::getBeamProfiles(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamProfiles[beam]; 
}

const operaInstrumentProfile *operaSpectralOrder::getBeamProfiles(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamProfiles[beam]; 
}

operaSpectralElements *operaSpectralOrder::getBackgroundElements(unsigned LeftOrRight) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundElements[LeftOrRight];
}

const operaSpectralElements *operaSpectralOrder::getBackgroundElements(unsigned LeftOrRight) const {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundElements[LeftOrRight];
}

operaExtractionAperture *operaSpectralOrder::getExtractionApertures(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return ExtractionApertures[beam];
}

const operaExtractionAperture *operaSpectralOrder::getExtractionApertures(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return ExtractionApertures[beam];
}

operaExtractionAperture *operaSpectralOrder::getBackgroundApertures(unsigned LeftOrRight) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundApertures[LeftOrRight];
}

const operaExtractionAperture *operaSpectralOrder::getBackgroundApertures(unsigned LeftOrRight) const {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BackgroundApertures[LeftOrRight];
}

void operaSpectralOrder::setBeamElements(unsigned beam, operaSpectralElements *beamElements) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (BeamElements[beam]) {
		delete BeamElements[beam];
	}
	BeamElements[beam] = beamElements;
}

void operaSpectralOrder::setBeamProfiles(unsigned beam, operaInstrumentProfile *beamProfiles){
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
  	if (BeamProfiles[beam]) {
		delete BeamProfiles[beam];
	}
	BeamProfiles[beam] = beamProfiles;  
}

void operaSpectralOrder::setBackgroundElements(unsigned LeftOrRight, operaSpectralElements *backgroundElements) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (BackgroundElements[LeftOrRight]) {
		delete BackgroundElements[LeftOrRight];
	}
	BackgroundElements[LeftOrRight] = backgroundElements;
}

void operaSpectralOrder::setExtractionApertures(unsigned beam, operaExtractionAperture *extractionApertures) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (ExtractionApertures[beam]) {
		delete ExtractionApertures[beam];
	}
	ExtractionApertures[beam] = extractionApertures;
}

void operaSpectralOrder::setBackgroundApertures(unsigned LeftOrRight, operaExtractionAperture *backgroundApertures) {
#ifdef RANGE_CHECK
	if (LeftOrRight >= LEFTANDRIGHT) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	if (BackgroundApertures[LeftOrRight]) {
		delete BackgroundApertures[LeftOrRight];
	}
	BackgroundApertures[LeftOrRight] = backgroundApertures;
}

double operaSpectralOrder::getminwavelength() const {
	return minwavelength;
}

double operaSpectralOrder::getmaxwavelength() const {
	return maxwavelength;
}

bool operaSpectralOrder::gethasWavelengthRange() const {
	return hasWavelengthRange;
}

void operaSpectralOrder::setminwavelength(double wl) {
	minwavelength = wl;
}

void operaSpectralOrder::setmaxwavelength(double wl) {
	maxwavelength = wl;
}

void operaSpectralOrder::sethasWavelengthRange(bool hasRange) {
	hasWavelengthRange = hasRange;
}

double operaSpectralOrder::getsnrSpectralBinSize() const {
	return snrSpectralBinSize;
}

void operaSpectralOrder::setsnrSpectralBinSize(double spectralbinsize) {
	snrSpectralBinSize = spectralbinsize;
}

void operaSpectralOrder::NormalizeFlat(operaFITSImage &flatMatrix, operaFITSImage &outputMatrix, unsigned nx, unsigned ny, unsigned binsize){
	
	// This function only works for spectralElementHeight=1; IPysize=1; IPyampling=1, so we enforce these below.
	if (InstrumentProfile->getysize()!=1 || InstrumentProfile->getYsampling()!=1 || SpectralElements->getelementHeight()!=1){
		throw operaException("operaSpectralOrder: instrument profile ysize, ysampling and spectralelement height must be 1, not "+itos(InstrumentProfile->getysize())+" and "+dtos(SpectralElements->getelementHeight())+" and "+dtos(SpectralElements->getelementHeight())+" (resp.) as given.",operaErrorInstrumentProfileImproperDimensions, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (SpectralElements->getnSpectralElements() == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned NXPoints = InstrumentProfile->getxsize()*InstrumentProfile->getXsampling();	
	
	unsigned NumberofElementsToBin = binsize;
	
	// DT Jan 213 moved out of loop
	float *SubPixElemSampleMedianCounts = (float *)malloc(NXPoints*sizeof(float));
	if (!SubPixElemSampleMedianCounts) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElementCounts = (float *)malloc(NumberofElementsToBin*sizeof(float));	
	if (!SubPixElementCounts) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixXcoords = (float *)malloc(NXPoints*sizeof(float));
	if (!SubPixXcoords) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *y2 = (float *)malloc(NXPoints*sizeof(float));
	if (!y2) {
		throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	// Below it loops over all spectral elements
	for (unsigned indexElem=0; indexElem < SpectralElements->getnSpectralElements(); indexElem++) {		
		
		// Figure out which is the first and last elements given the number of elements to bin
		unsigned firstElement = indexElem - (unsigned)ceil(NumberofElementsToBin/2);
		unsigned lastElement =  indexElem + (unsigned)ceil(NumberofElementsToBin/2);
		
		// Enforce posive elements
		if (indexElem < ceil(NumberofElementsToBin/2)) {
			firstElement = 0;
		}
		
		// Enforce last element not to be greater than nelements
		if (lastElement > SpectralElements->getnSpectralElements()) {
			lastElement = SpectralElements->getnSpectralElements();
		}
		
		for (unsigned i=0; i<NXPoints; i++) {	
			
			unsigned indexElemInsideBin=0;
			
			for(unsigned iElem=firstElement; iElem < lastElement; iElem++) {
				float xcenter = SpectralElements->getphotoCenterX(iElem);
				float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);				
				unsigned xx = (unsigned)floor(XCenterOfSubPixel);				
				unsigned yy = (unsigned)round(SpectralElements->getphotoCenterY(iElem));
				
				SubPixElementCounts[indexElemInsideBin++] = (float)flatMatrix[yy][xx];
			} 
			SubPixElemSampleMedianCounts[i] = operaArrayMedianQuick(indexElemInsideBin,SubPixElementCounts);		
		}
		free(SubPixElementCounts);
		
		/***** Prep for interpolation **********/
		for (unsigned i=0; i<NXPoints; i++) {
			SubPixXcoords[i] = InstrumentProfile->getIPixXCoordinate(i);
		}
		
		float yp1 = (SubPixElemSampleMedianCounts[1] - SubPixElemSampleMedianCounts[0])/(SubPixXcoords[1] - SubPixXcoords[0]);
		float ypn = (SubPixElemSampleMedianCounts[NXPoints-1] - SubPixElemSampleMedianCounts[NXPoints-2])/(SubPixXcoords[NXPoints-1] - SubPixXcoords[NXPoints-2]);
		
		// Call cubicspline to get second derivatives 
		cubicspline(SubPixXcoords, SubPixElemSampleMedianCounts, NXPoints, yp1, ypn, y2);
		/**** End of prep for interpolation ***/
		
		unsigned xlocalmin = (unsigned)floor(SpectralElements->getphotoCenterX(indexElem) - Geometry->getapertureWidth()/2);
		unsigned xlocalmax = (unsigned)ceil(SpectralElements->getphotoCenterX(indexElem) + Geometry->getapertureWidth()/2);
		
		if (xlocalmin < 0) { xlocalmin = 0;}		
		if (xlocalmax > nx) { xlocalmax = nx;}			
		
		unsigned yy = (unsigned)round(SpectralElements->getphotoCenterY(indexElem));
		
		for(unsigned xx=xlocalmin; xx<xlocalmax; xx++) {
			float ExpectedFlux;
			float xcoordInIPPixunits = (float)xx + 0.5 - SpectralElements->getphotoCenterX(indexElem);
			// Call splineinterpolate for interpolations 
			splineinterpolate(SubPixXcoords, SubPixElemSampleMedianCounts, y2, NXPoints, xcoordInIPPixunits, &ExpectedFlux);			
			outputMatrix[yy][xx] = flatMatrix[yy][xx]/ExpectedFlux;
		}
	}
	// DT Jan 2013 moved outside of loop
	free(y2);
	free(SubPixXcoords);
	free(SubPixElemSampleMedianCounts);
	free(SubPixElementCounts);
}

void operaSpectralOrder::calculateSNR(void) {
    
	unsigned length = SpectralElements->getnSpectralElements();
    
	for (unsigned indexElem=0; indexElem < length; indexElem++){
		SpectralElements->getFluxSNR(indexElem);
	}
	SNR = getCentralSmoothedSNR(length/DEFAULT_SNR_SMOOTHING);
	//SNR = SpectralElements->getFluxSNR(length/2);
	hasSNR = true;
}

float operaSpectralOrder::getCentralSmoothedSNR(int upperlowerbound) const {
    
	float snrs[MAXREFWAVELENGTHSPERORDER];
	unsigned length = SpectralElements->getnSpectralElements();
    
	if (length == 0) {
		return 1.0;
	}
	if (upperlowerbound == 0) {
		return SpectralElements->getFluxSNR(length/2);
	}
	int start = max(length/2-upperlowerbound, 0);
	int stop = min(length/2+upperlowerbound, length-1);
	int range = stop-start+1;
	if ((range > MAXREFWAVELENGTHSPERORDER) || (range == 0)) {
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned i = 0;
	for (int indexElem=start; indexElem <= stop; indexElem++){
		snrs[i++]= SpectralElements->getFluxSNR(indexElem);
	}
	return operaArrayMedianQuick(range, snrs);
}

float operaSpectralOrder::getPeakSmoothedSNR(int upperlowerbound) const {
    
	float snrs[MAXREFWAVELENGTHSPERORDER];
	unsigned length = SpectralElements->getnSpectralElements();
	unsigned maxindex = length/2;
    float maxSNR = 0.0;
	
	for (unsigned indexElem=0; indexElem < length; indexElem++){
		float snr = SpectralElements->getFluxSNR(indexElem);
		if (snr > maxSNR) {
			maxindex = indexElem;
			maxSNR = snr;
		}
	}
	if (length == 0) {
		return 1.0;
	}
	if (upperlowerbound == 0) {
		return SpectralElements->getFluxSNR(length/2);
	}
	int start = max(maxindex-upperlowerbound, 0);
	int stop = min(maxindex+upperlowerbound, length-1);
	int range = stop-start+1;
	if ((range > MAXREFWAVELENGTHSPERORDER) || (range == 0)) {
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned i = 0;
	for (int indexElem=start; indexElem <= stop; indexElem++){
		snrs[i++]= SpectralElements->getFluxSNR(indexElem);
	}
	return operaArrayMedianQuick(range, snrs);
}

void operaSpectralOrder::extractRawSum(operaFITSImage &inputImage, operaFITSProduct &outputSpectrum) {
	
	float distanceInPixelUnits;
	
	for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
		unsigned xlocalmin = (unsigned)floor(SpectralElements->getphotoCenterX(indexElem) - getGeometry()->getapertureWidth()/2);
		unsigned xlocalmax = (unsigned)ceil(SpectralElements->getphotoCenterX(indexElem) + getGeometry()->getapertureWidth()/2);
		unsigned ylocalmin = (unsigned)round(SpectralElements->getphotoCenterY(indexElem) - SpectralElements->getelementHeight()/2);
		unsigned ylocalmax = (unsigned)ceil(SpectralElements->getphotoCenterY(indexElem) + SpectralElements->getelementHeight()/2);
		
		if (xlocalmin < 0) { xlocalmin = 0;}		
		if (xlocalmax > inputImage.getnaxis1()) { xlocalmax = inputImage.getnaxis1();}		
		if (ylocalmin < 0) { ylocalmin = 0;}
		if (ylocalmax > inputImage.getnaxis2()) { ylocalmax = inputImage.getnaxis2();}
		
#ifdef PRINT_DEBUG
		cout << "Elem # " << indexElem << 
		"; PhotoCenter: [" << SpectralElements->getphotoCenterX(indexElem) << ":" << SpectralElements->getphotoCenterY(indexElem) << "]" <<
		"\tXRANGE= [" << xlocalmin <<
		":" << xlocalmax - 1 << "]" <<
		"\tYRANGE= [" << ylocalmin <<
		":" << ylocalmax - 1 << "]" << endl;
#endif				
		
		float sumFlux = 0.0;
		for(unsigned y=ylocalmin;y<ylocalmax;y++) {
			for(unsigned x=xlocalmin;x<xlocalmax;x++) {
				sumFlux += inputImage[y][x];
			}
		}
		
		SpectralElements->setFlux(sumFlux,indexElem);
		distanceInPixelUnits = SpectralElements->getdistd(indexElem);
		
		float table[3] = {(float)getorder(), distanceInPixelUnits, (float)SpectralElements->getFlux(indexElem)};
		outputSpectrum.addRow(table, 3, 0);
	}
	setSpectrumType(RawSpectrum);
	SpectralElements->setSpectrumType(RawSpectrum);
	sethasSpectralElements(true);
}

void operaSpectralOrder::extractRawSum(operaFITSImage &inputImage, ofstream &sout){
	
	float distanceInPixelUnits;
	
	for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
		unsigned xlocalmin = (unsigned)floor(SpectralElements->getphotoCenterX(indexElem) - getGeometry()->getapertureWidth()/2);
		unsigned xlocalmax = (unsigned)ceil(SpectralElements->getphotoCenterX(indexElem) + getGeometry()->getapertureWidth()/2);
		unsigned ylocalmin = (unsigned)round(SpectralElements->getphotoCenterY(indexElem) - SpectralElements->getelementHeight()/2);
		unsigned ylocalmax = (unsigned)ceil(SpectralElements->getphotoCenterY(indexElem) + SpectralElements->getelementHeight()/2);
		
		if (xlocalmin < 0) { xlocalmin = 0;}		
		if (xlocalmax > inputImage.getnaxis1()) { xlocalmax = inputImage.getnaxis1();}		
		if (ylocalmin < 0) { ylocalmin = 0;}
		if (ylocalmax > inputImage.getnaxis2()) { ylocalmax = inputImage.getnaxis2();}
		
#ifdef PRINT_DEBUG
		cout << "Elem # " << indexElem << 
		"; PhotoCenter: [" << SpectralElements->getphotoCenterX(indexElem) << ":" << SpectralElements->getphotoCenterY(indexElem) << "]" <<
		"\tXRANGE= [" << xlocalmin <<
		":" << xlocalmax - 1 << "]" <<
		"\tYRANGE= [" << ylocalmin <<
		":" << ylocalmax - 1 << "]" << endl;
#endif				
		
		float sumFlux = 0.0;
		for(unsigned y=ylocalmin;y<ylocalmax;y++) {
			for(unsigned x=xlocalmin;x<xlocalmax;x++) {
				sumFlux += inputImage[y][x];
			}
		}
		
		SpectralElements->setFlux(sumFlux,indexElem);
		distanceInPixelUnits = SpectralElements->getdistd(indexElem);	
		
#ifdef PRINT_DEBUG
		cout << orderNumber <<
		"\t" << indexElem << 
		"\t" << SpectralElements->getphotoCenterX(indexElem) << 
		"\t" << SpectralElements->getphotoCenterY(indexElem) <<
		"\t" << xlocalmin <<
		"\t" << xlocalmax + 1 <<
		"\t" << ylocalmin <<
		"\t" << ylocalmax + 1 <<				
		"\t" << distanceInPixelUnits <<		
		"\t" << SpectralElements->getFlux(indexElem)<< endl;
#endif				
		
		// DT made to be more like upena...
		sout << distanceInPixelUnits <<
		"\t" << SpectralElements->getFlux(indexElem) <<
		"\t" << getorder() << endl;				
	}
	setSpectrumType(RawSpectrum);
	SpectralElements->setSpectrumType(RawSpectrum);
	sethasSpectralElements(true);
}

void operaSpectralOrder::extractRawSum(operaFITSImage &inputImage, float noise, float gain){
	
	for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
		unsigned xlocalmin = (unsigned)floor(SpectralElements->getphotoCenterX(indexElem) - getGeometry()->getapertureWidth()/2);
		unsigned xlocalmax = (unsigned)ceil(SpectralElements->getphotoCenterX(indexElem) + getGeometry()->getapertureWidth()/2);
		unsigned ylocalmin = (unsigned)round(SpectralElements->getphotoCenterY(indexElem) - SpectralElements->getelementHeight()/2);
		unsigned ylocalmax = (unsigned)ceil(SpectralElements->getphotoCenterY(indexElem) + SpectralElements->getelementHeight()/2);
		
		if (xlocalmin < 0) { xlocalmin = 0;}		
		if (xlocalmax > inputImage.getnaxis1()) { xlocalmax = inputImage.getnaxis1();}		
		if (ylocalmin < 0) { ylocalmin = 0;}
		if (ylocalmax > inputImage.getnaxis2()) { ylocalmax = inputImage.getnaxis2();}
		
#ifdef PRINT_DEBUG
		cout << "Elem # " << indexElem << 
		"; PhotoCenter: [" << SpectralElements->getphotoCenterX(indexElem) << ":" << SpectralElements->getphotoCenterY(indexElem) << "]" <<
		"\tXRANGE= [" << xlocalmin <<
		":" << xlocalmax - 1 << "]" <<
		"\tYRANGE= [" << ylocalmin <<
		":" << ylocalmax - 1 << "]" << endl;
#endif				
		
		float sumFlux = 0.0;
		float sumVar = 0.0;
		for(unsigned y=ylocalmin;y<ylocalmax;y++) {
			for(unsigned x=xlocalmin;x<xlocalmax;x++) {
				float flux = inputImage[y][x];
				sumFlux += flux*gain;
				sumVar += ((noise*noise) + fabs(flux*gain)); // detector noise + photon noise ;
			}
		}
		
		SpectralElements->setFlux(sumFlux,indexElem);
		SpectralElements->setFluxVariance(sumVar,indexElem);
		
#ifdef PRINT_DEBUG
		float distanceInPixelUnits = SpectralElements->getdistd(indexElem);	
		cout << orderNumber <<
		"\t" << indexElem << 
		"\t" << SpectralElements->getphotoCenterX(indexElem) << 
		"\t" << SpectralElements->getphotoCenterY(indexElem) <<
		"\t" << xlocalmin <<
		"\t" << xlocalmax + 1 <<
		"\t" << ylocalmin <<
		"\t" << ylocalmax + 1 <<				
		"\t" << distanceInPixelUnits <<		
		"\t" << SpectralElements->getFlux(indexElem)<< endl;
#endif				
	}
	setSpectrumType(RawSpectrum);
	SpectralElements->setSpectrumType(RawSpectrum);
	sethasSpectralElements(true);
}

void operaSpectralOrder::CalculateWavelengthSolution(void) {
	if (gethasGeometry() && gethasWavelength()) {
		operaWavelength *wavelength = getWavelength();
		operaGeometry *geometry = getGeometry();
		operaSpectralElements *spectralElements = getSpectralElements();
        double dmin = 0.0;
		wavelength->setDmin(dmin);
        double dmax = (double)geometry->CalculateDistance(geometry->getYmin(),geometry->getYmax());
        wavelength->setDmax(dmax);
		for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
			spectralElements->setwavelength(wavelength->evaluateWavelength(spectralElements->getdistd(k)), k);
		}
		spectralElements->setHasWavelength(true);
		sethasWavelength(true);
	}
}

void operaSpectralOrder::measureInstrumentProfileAlongRows(operaFITSImage &masterFlatImage, unsigned binsize, unsigned sampleElementForPlot, ostream *pout){
	
	if (InstrumentProfile->getysize()!=1 || InstrumentProfile->getYsampling()!=1){
		throw operaException("operaSpectralOrder: instrument profile ysize and ysampling must be 1.",operaErrorInstrumentProfileImproperDimensions, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	unsigned NXPoints = InstrumentProfile->getNXPoints();		
	unsigned NYPoints = InstrumentProfile->getNYPoints();
	
	if (NXPoints == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElemSampleMedianCounts = new float[NXPoints];
	float *SubPixElemSampleMedSigCounts = new float[NXPoints];	
	unsigned NumberOfElementSamples;
	unsigned NumberofElementsToBin = binsize;
	
	NumberOfElementSamples = (unsigned)ceil((float)SpectralElements->getnSpectralElements()/(float)NumberofElementsToBin); 
	
	if (NumberOfElementSamples == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElementCounts = new float[NumberOfElementSamples];		
	float *SumSubPixElements = new float[NumberOfElementSamples];
	
	for(unsigned k=0;k<NumberOfElementSamples;k++){
		
		unsigned firstElement = NumberofElementsToBin*(k);
		unsigned lastElement =  NumberofElementsToBin*(k+1);
		if (lastElement > SpectralElements->getnSpectralElements()){
            lastElement = SpectralElements->getnSpectralElements();   
        }
		
		SumSubPixElements[k] = 0;
		
        float distdmidElem = (float)fabs(SpectralElements->getdistd(lastElement-1) + SpectralElements->getdistd(firstElement))/2;		
        
        for (unsigned i=0; i<NXPoints; i++) {	
            
            unsigned indexElemInsideBin=0;
            
            for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
                
                float xcenter =  SpectralElements->getphotoCenterX(indexElem);
                float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);				
                unsigned xx = (unsigned)floor(XCenterOfSubPixel);				
                
                float ycenter =  SpectralElements->getphotoCenterY(indexElem);
                float YCenterOfSubPixel = ycenter + InstrumentProfile->getIPixYCoordinate(0);                
                unsigned yy = (unsigned)floor(YCenterOfSubPixel);
                
                SubPixElementCounts[indexElemInsideBin++] = (float)masterFlatImage[yy][xx]/((float)InstrumentProfile->getXsampling()*(float)InstrumentProfile->getYsampling());
            } 
            
            SubPixElemSampleMedianCounts[i] = operaArrayMedian(indexElemInsideBin,SubPixElementCounts);
            SumSubPixElements[k] += SubPixElemSampleMedianCounts[i];		
            SubPixElemSampleMedSigCounts[i] = operaArrayMedianSigma(indexElemInsideBin, SubPixElementCounts, SubPixElemSampleMedianCounts[i]);
        }
        
        for (unsigned j=0; j<NYPoints; j++) {			
            for (unsigned i=0; i<NXPoints; i++) {
                InstrumentProfile->setdataCubeValues(SubPixElemSampleMedianCounts[i]/SumSubPixElements[k],i,j,k);
                InstrumentProfile->seterrorsCubeValues(SubPixElemSampleMedSigCounts[i]/SumSubPixElements[k],i,j,k);				
            }
        }
        InstrumentProfile->setdistd(distdmidElem,k);
        
#ifdef PRINT_DEBUG        
        unsigned midElement = firstElement + fabs(lastElement - firstElement)/2;        
        unsigned MinimumBinSize = Geometry->CalculateMinimumYBinSize(SpectralElements->getphotoCenterY(midElement));
        cout << k << " " 
        << SpectralElements->getphotoCenterY(midElement) << " " 
        << binsize << " " 
        << MinimumBinSize
        << endl;        
#endif        
        
	}
	
	InstrumentProfile->normalizeCubeData();
	
	InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);
	
	if (pout) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
    
	delete[] SubPixElementCounts;	
	delete[] SubPixElemSampleMedianCounts;	
	delete[] SubPixElemSampleMedSigCounts;	
	delete[] SumSubPixElements;	
}

void operaSpectralOrder::measureInstrumentProfileAlongRowsInto2DWithGaussian(operaFITSImage &masterFlatImage, operaFITSImage &badpix, unsigned binsize, float gaussSig, float tiltInDegrees, bool witherrors, unsigned sampleElementForPlot, ostream *pout, const int minimumLines){
	
	int npar = 3;
	double par[3];
	par[0] = 1.0/((double)gaussSig*sqrt(2.0*M_PI));
	par[1] = 0;
	par[2] = (double)gaussSig;    
	
    float tangentOfTiltInDegrees = tan(tiltInDegrees*M_PI/180.0);
    
	unsigned NXPoints = InstrumentProfile->getNXPoints();		
	unsigned NYPoints = InstrumentProfile->getNYPoints();
	
	if (NXPoints == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElemSampleMedianCounts = new float[NXPoints];
    float *SubPixElemSampleMedSigCounts = new float[NXPoints];
	
	unsigned NumberOfElementSamples;
	unsigned NumberofElementsToBin = binsize;
	
	NumberOfElementSamples = (unsigned)ceil((float)SpectralElements->getnSpectralElements()/(float)NumberofElementsToBin); 
    
	if (NumberofElementsToBin == 0){
		throw operaException("operaSpectralOrder: ",operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	float *SubPixElementCounts = new float[NumberofElementsToBin];
	
    unsigned ndataPoints = 0;
    
	for(unsigned k=0;k<NumberOfElementSamples;k++){
		
		unsigned firstElement = NumberofElementsToBin*(k);
		unsigned lastElement =  NumberofElementsToBin*(k+1);
		if (lastElement > SpectralElements->getnSpectralElements()){
            lastElement = SpectralElements->getnSpectralElements();   
        }
		
        float IPNormalizationFactor = 0;
        float fluxFractionLost = 0; 
        
		for (unsigned i=0; i<NXPoints; i++) {
            
			unsigned indexElemInsideBin=0;
			
			for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
				
                float xcenter =  SpectralElements->getphotoCenterX(indexElem);
				float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);
                unsigned xx = (unsigned)floor(XCenterOfSubPixel);				
                
                unsigned nps = 0;
                float SubPixMeanCounts = 0;
                
                if (xx > 0 && xx < masterFlatImage.getnaxis1()) {
                    float ycenter =  SpectralElements->getphotoCenterY(indexElem);
                    unsigned ylocalmin = (unsigned)round(ycenter -  SpectralElements->getelementHeight()/2);
                    unsigned ylocalmax = (unsigned)ceil(ycenter +  SpectralElements->getelementHeight()/2);
                    
                    for(unsigned yy=ylocalmin; yy<ylocalmax; yy++){
                        if (yy > 0 && yy < masterFlatImage.getnaxis2() &&
                            masterFlatImage[yy][xx] < SATURATIONLIMIT && 
                            badpix[yy][xx] == 1 ) { 
                            
                            SubPixMeanCounts += (float)masterFlatImage[yy][xx]; // in units of ADU					
                            nps++;
                        }
                    }
                }
                if(nps) {
                    SubPixElementCounts[indexElemInsideBin++] = SubPixMeanCounts/(float)(nps*InstrumentProfile->getXsampling());
                }                
			} 
            
            if(indexElemInsideBin > 2) {
                SubPixElemSampleMedianCounts[i] = operaArrayMedian(indexElemInsideBin,SubPixElementCounts);
                IPNormalizationFactor += SubPixElemSampleMedianCounts[i];		                
                if(witherrors) {
                    SubPixElemSampleMedSigCounts[i] = operaArrayMedianSigma(indexElemInsideBin, SubPixElementCounts, SubPixElemSampleMedianCounts[i]);
                }
            } else if (indexElemInsideBin > 0 && indexElemInsideBin < 3) {
                SubPixElemSampleMedianCounts[i] = operaArrayMean(indexElemInsideBin,SubPixElementCounts);
                IPNormalizationFactor += SubPixElemSampleMedianCounts[i];                
                if(witherrors) {
                    SubPixElemSampleMedSigCounts[i] = operaArraySigma(indexElemInsideBin,SubPixElementCounts);
                }
            } else {
                SubPixElemSampleMedianCounts[i] = NAN; 
                fluxFractionLost += 1.0/(float)(InstrumentProfile->getXsampling()*NXPoints);                
                if(witherrors) {
                    SubPixElemSampleMedSigCounts[i] = NAN;
                }
            }
		}
        
        if(IPNormalizationFactor) {
            float distdmidElem = (float)fabs(SpectralElements->getdistd(lastElement-1) + SpectralElements->getdistd(firstElement))/2;		
			
            for (unsigned j=0; j<NYPoints; j++) {	
                double IPPixYCoord = (double)(InstrumentProfile->getIPixYCoordinate(j));				
                for (unsigned i=0; i<NXPoints; i++) {
                    float ipvalue = 0.0;
                    float iperr = 0.0;
                    if(!isnan(SubPixElemSampleMedianCounts[i])) { 
                        ipvalue = (SubPixElemSampleMedianCounts[i]/IPNormalizationFactor)*(1+fluxFractionLost/IPNormalizationFactor);
                        if(witherrors) {
                            float ipvar = SubPixElemSampleMedSigCounts[i]*SubPixElemSampleMedSigCounts[i];
                            iperr = sqrt((ipvar/IPNormalizationFactor)*(1+fluxFractionLost/IPNormalizationFactor)); 
                        }  
                    } else {
                        ipvalue = 1.0/(float)(InstrumentProfile->getXsampling()*NXPoints);
                        if(witherrors) {
                            iperr = NAN;  
                        }
                    }
                    par[1]  = tangentOfTiltInDegrees*(double)(InstrumentProfile->getIPixXCoordinate(i));
                    InstrumentProfile->setdataCubeValues(ipvalue*(float)GaussianFunction(IPPixYCoord,par,npar),i,j,ndataPoints);
                    if(witherrors) {
                        InstrumentProfile->seterrorsCubeValues(iperr,i,j,ndataPoints);	                    
                    }
                }
            }        
            
            InstrumentProfile->setdistd(distdmidElem,ndataPoints);
            ndataPoints++;
        }
	}
    
    InstrumentProfile->setnDataPoints(ndataPoints);
    
    if(ndataPoints > minimumLines) {
        InstrumentProfile->normalizeCubeData();        
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,witherrors);
        
        if (pout != NULL) {
            InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
        }
        sethasInstrumentProfile(true);
    } else {
        sethasInstrumentProfile(false);
    }
    
    delete[] SubPixElementCounts;
	delete[] SubPixElemSampleMedianCounts;
    delete[] SubPixElemSampleMedSigCounts;
}

void operaSpectralOrder::setSpectralElementsByHeight(double Height) {
	
    if (!gethasGeometry()){
        throw operaException("operaSpectralOrder::setSpectralElementsByHeight: ",operaErrorHasNoGeometry, __FILE__, __FUNCTION__, __LINE__);
	}
    if(Height <= 0.0) {
        throw operaException("operaSpectralOrder:",operaErrorDivideByZeroError, __FILE__, __FUNCTION__, __LINE__);
    } 
	unsigned nElements = (unsigned)(ceil(fabs(Geometry->getYmax() - Geometry->getYmin())/Height));
	
    if(nElements == 0) {
        throw operaException("operaSpectralOrder:",operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
    } 
	if (SpectralElements == NULL) {
		SpectralElements = new operaSpectralElements(nElements, getSpectrumType());	
	} else {
		SpectralElements->resize(nElements);
		SpectralElements->setSpectrumType(getSpectrumType());
	}
	
	SpectralElements->setelementHeight(Height);
	
	double ymin = Geometry->getYmin();
	
	double ycenter = ymin + Height/2.0;
	
	for(unsigned indexElem=0; indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
        double d = Geometry->CalculateDistance(ymin, ycenter);
		SpectralElements->setphotoCenter(Geometry->getCenterPolynomial()->Evaluate(ycenter), ycenter, indexElem);
		SpectralElements->setdistd(d, indexElem);
		ycenter += Height;		
	}
	Geometry->CalculateAndSetOrderLength();
	SpectralElements->setHasDistance(true);
	sethasSpectralElements(true);
}

void operaSpectralOrder::setSpectralElementsByStitchingApertures(double effectiveApertureFraction) {
    
    if (!gethasGeometry()){
        throw operaException("operaSpectralOrder::setSpectralElementsByHeight: ",operaErrorHasNoGeometry, __FILE__, __FUNCTION__, __LINE__);
	}
    
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::setSpectralElementsByHeight: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double elementHeight = ExtractionApertures[0]->getLineAperture()->getLineYWidth()*effectiveApertureFraction;
    
    if (SpectralElements == NULL) {
		SpectralElements = new operaSpectralElements(MAXSPECTRALELEMENTSPERORDER, getSpectrumType());
	} else {
		SpectralElements->resize(MAXSPECTRALELEMENTSPERORDER);
		SpectralElements->setSpectrumType(getSpectrumType());
	}
    
    SpectralElements->setelementHeight(elementHeight);
	
    // First remember that the order polynomial is defined as:
    // x(y) = aa + bb*y + cc*y*y
    double aa = (double)Geometry->getCenterPolynomial()->getCoefficient(0);
    double bb = (double)Geometry->getCenterPolynomial()->getCoefficient(1);
    double cc = (double)Geometry->getCenterPolynomial()->getCoefficient(2);
    
    double ymin = (double)Geometry->getYmin();
    double ymax = (double)Geometry->getYmax();
	
    // Pick y coordinate of first element:
    double firstYCenter = ymin + elementHeight/2.0;
    double firstXCenter = (double)Geometry->getCenterPolynomial()->Evaluate((double)firstYCenter);
    // Then calculate the inital aperture line in the image reference system:
    
    // And the aperture central line is defined as: y(x) = dd*x + ff
    double dd = (double)ExtractionApertures[0]->getLineAperture()->getSlope();
    double ff = (double)ExtractionApertures[0]->getLineAperture()->getIntercept();
    
	//    Geometry->getCenterPolynomial()->printEquation(&cout);
    
    double xcenter = firstXCenter;
    double ycenter = firstYCenter;
    
    SpectralElements->setphotoCenter(xcenter, ycenter, 0);
    SpectralElements->setdistd(Geometry->CalculateDistance(ymin, ycenter), 0);
	
    unsigned nElements = 1;
    
    for(unsigned indexElem=1; indexElem<MAXSPECTRALELEMENTSPERORDER; indexElem++) {
#ifdef PRINT_DEBUG
        // Now let's bring the aperture line to the image reference system:
        // y(x) - (firstYCenter + indexElem*elementHeight) = dd*(x - firstXcenter) + ff
        // And let's invert this equation and rearrange the terms:
        // x(y) = ([y - (firstYCenter + indexElem*elementHeight + ff)] / dd) + (firstXcenter)
        // Then we have the new line given by: x(y) = (1/dd)*y + firstXcenter - (firstYCenter + indexElem*elementHeight + ff)/dd
		// And intersect = firstXcenter - (firstYCenter + (double)indexElem*elementHeight + ff)/dd;
#endif
        
        double slope = (1.0/dd);
        double intersect = xcenter - (ycenter + elementHeight - ff)/dd;
		
#ifdef PRINT_DEBUG
		//      Below we print the line equation
		//        cout << "g" << indexElem << "(x)="<< slope << "*x+" << intersect << endl << endl;
        
        // Now we need to find a point {xcenter,ycenter} for the intersection
        // between the order polynomial and the aperture central line
		
        // and let's solve for the y-coordinate of the crossing point:
        // y*slope + intersect = aa + bb*y + cc*y*y
        // rearranging the terms:
        // alpha + beta*y + cc*y*y = 0, where alpha = (aa-intersect), beta = (bb-slope)
#endif
        double alpha = aa - intersect;
        double beta = bb - slope;
		
        // Applying Bhaskara's equation we have:
        // y = (-beta +/- sqrt(beta^2 - 4*alpha*cc)/(2*alpha)
		
        double ycenter_plus, ycenter_minus;
        if((-beta + sqrt(beta*beta - 4*alpha*cc)) != 0 && (beta*beta - 4*alpha*cc) >= 0) {
            ycenter_plus = (2*alpha)/(-beta + sqrt(beta*beta - 4*alpha*cc));
            ycenter_minus = (2*alpha)/(-beta - sqrt(beta*beta - 4*alpha*cc));
        } else {
            //throw operaException("operaSpectralOrder:",operaErrorDivideByZeroError, __FILE__, __FUNCTION__, __LINE__);
			// DT May 20 2014 -- We are likely to run in to problems when indexing through MAXSPECTRALELEMENTSPERORDER
			// so, just break rather than throwing an exception and thus aborting...
			break;
            //throw operaException("operaSpectralOrder:",operaErrorDivideByZeroError, __FILE__, __FUNCTION__, __LINE__);
		}
		//        cout << "ycenter_plus=" << ycenter_plus << " ycenter_minus=" << ycenter_minus << endl;
        
        // pick the branch with closest y to the previous point:
        double dist_yplus = fabs(ycenter_plus - ycenter);
        double dist_yminus = fabs(ycenter_minus - ycenter);
		//        cout << "dist_yplus=" << dist_yplus << " dist_yminus=" << dist_yminus << endl;
		
        if(dist_yminus < dist_yplus) {
            ycenter = ycenter_minus;
        } else {
            ycenter = ycenter_plus;
        }
        
        //ycenter = (2*alpha)/(-beta - sqrt(beta*beta - 4*alpha*cc));
        xcenter = (double)Geometry->getCenterPolynomial()->Evaluate((double)ycenter);
        
        SpectralElements->setphotoCenter(xcenter, ycenter, indexElem);
        SpectralElements->setdistd((double)Geometry->CalculateDistance(ymin, ycenter), indexElem);
        
        nElements++;
        if(ycenter + elementHeight/2 >= ymax) {
            break;
        }
    }
    SpectralElements->resize(nElements);
	
	Geometry->CalculateAndSetOrderLength();
	SpectralElements->setHasDistance(true);
	sethasSpectralElements(true);
}

void operaSpectralOrder::setApertureElements(operaSpectralOrder_t format) {
    if(!hasSpectralElements) {
        throw operaException("operaSpectralOrder:",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);	
    } 
    
	// only reallocate if we need to...
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
		if (BeamElements[beam] == NULL) {
			BeamElements[beam] = new operaSpectralElements(SpectralElements->getnSpectralElements(), format);    			
		} else {
			if (BeamElements[beam]->getnSpectralElements() < SpectralElements->getnSpectralElements()) {
				if (BeamElements[beam]) {
					delete BeamElements[beam];
				}
				BeamElements[beam] = new operaSpectralElements(SpectralElements->getnSpectralElements(), format);    			
			} else {
				BeamElements[beam]->setnSpectralElements(SpectralElements->getnSpectralElements());
			}
		}
    }  
    for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
		if (BackgroundElements[backgroundIndex] == NULL) {
			BackgroundElements[backgroundIndex] = new operaSpectralElements(SpectralElements->getnSpectralElements(), format);    
		} else {
			if (BackgroundElements[backgroundIndex]->getnSpectralElements() < SpectralElements->getnSpectralElements()) {
				if (BackgroundElements[backgroundIndex]) {
					delete BackgroundElements[backgroundIndex];
				}
				BackgroundElements[backgroundIndex] = new operaSpectralElements(SpectralElements->getnSpectralElements(), format); 
			} else {
				BackgroundElements[backgroundIndex]->setnSpectralElements(SpectralElements->getnSpectralElements());
			}
		}
    }    
	
	for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {     
		
        double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
        double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double beamXcenter = elemXcenter + (double)ExtractionApertures[beam]->getLineAperture()->getMidPoint()->getXcoord();
            double beamYcenter = elemYcenter + (double)ExtractionApertures[beam]->getLineAperture()->getMidPoint()->getYcoord();
			
            BeamElements[beam]->setphotoCenter(beamXcenter,beamYcenter,indexElem);    
        }
        for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
            double backgroundXcenter = elemXcenter + (double)BackgroundApertures[backgroundIndex]->getLineAperture()->getMidPoint()->getXcoord();
            double backgroundYcenter = elemYcenter + (double)BackgroundApertures[backgroundIndex]->getLineAperture()->getMidPoint()->getYcoord();
            
            BackgroundElements[backgroundIndex]->setphotoCenter(backgroundXcenter,backgroundYcenter,indexElem);                
        }           
    }    
}

operaExtractionAperture *operaSpectralOrder::calculateMainApertureFromExtractionBeams(bool useIP) {
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::calculateMainApertureFromExtractionBeams: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaPoint point(0,0);
    float mainSlope = ExtractionApertures[0]->getLineAperture()->getSlope();
    float mainWidth = ExtractionApertures[0]->getLineAperture()->getWidth();
    float mainLength = 0;
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        mainLength += ExtractionApertures[beam]->getLineAperture()->getLength();
    }
	
    Line MainExtractionLineAperture(mainSlope, mainWidth, mainLength, &point);
    
    operaExtractionAperture *MainExtractionAperture;
    
    if(useIP) {
        if(!gethasInstrumentProfile()) {
            throw operaException("operaSpectralOrder::calculateMainApertureFromExtractionBeams: ",operaErrorHasNoInstrumentProfile, __FILE__, __FUNCTION__, __LINE__);
        }
		MainExtractionAperture = new operaExtractionAperture(&MainExtractionLineAperture,InstrumentProfile);
    } else {
        unsigned xsampling = ExtractionApertures[0]->getXsampling();
        unsigned ysampling = ExtractionApertures[0]->getYsampling();
        MainExtractionAperture = new operaExtractionAperture(&MainExtractionLineAperture, xsampling, ysampling);
    }
    
    return MainExtractionAperture;
}


void operaSpectralOrder::deleteInstrumentProfile(void) {
	
	if (InstrumentProfile) {
		delete InstrumentProfile;
	}
	InstrumentProfile = NULL;
	sethasInstrumentProfile(false);
}

//
// Caution -- replaces the instrumentProfile pointer!!!!
//
void operaSpectralOrder::setInstrumentProfileVector(unsigned IPxsize, unsigned IPxsampling, unsigned IPysize, unsigned IPysampling, unsigned NDataPoints) {
	
	if (InstrumentProfile) {
		delete InstrumentProfile;
	}
	InstrumentProfile = new operaInstrumentProfile(IPxsize,IPxsampling,IPysize,IPysampling, NDataPoints);
	sethasInstrumentProfile(true);
}

/*
 * interface used by read/write orders
 */
void operaSpectralOrder::setSpectralLines(operaSpectralLines *spectralLines) {
	SpectralLines = spectralLines;
	sethasSpectralLines(true);
}

void operaSpectralOrder::setSpectralLines(operaFITSImage &masterCompImage, operaFITSImage &badpix, operaFITSImage &bias, float noise, float gain, float ReferenceLineWidth,float DetectionThreshold, float LocalMaxFilterWidth, float MinPeakDepth) {
	
    if(!gethasSpectralElements()) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    if(!gethasInstrumentProfile()) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoInstrumentProfile, __FILE__, __FUNCTION__, __LINE__);
    }    
    if(!SpectralElements->getHasXCorrelation()) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoXCorrelation, __FILE__, __FUNCTION__, __LINE__);
    }
	
    if(!SpectralElements->getHasRawSpectrum()) { 
        
        unsigned NXPoints = InstrumentProfile->getNXPoints();		
        unsigned NYPoints = InstrumentProfile->getNYPoints();
        unsigned nMaxDataPoints = SpectralElements->getnSpectralElements();
        
        for(unsigned indexElem=0;indexElem < nMaxDataPoints; indexElem++) {    
            
            float xcenter =  SpectralElements->getphotoCenterX(indexElem);
            float ycenter =  SpectralElements->getphotoCenterY(indexElem);
            float distdElem = SpectralElements->getdistd(indexElem);
            
            float my = 0;
            float myvar = 0; 
            float fluxFractionLost = 0;
            float IPNormalizationFactor = 0;  
            
            for (unsigned j=0; j<NYPoints; j++) {	
                float YCenterOfSubPixel = ycenter + InstrumentProfile->getIPixYCoordinate(j);				
                unsigned yy = (unsigned)floor(YCenterOfSubPixel);             
                for (unsigned i=0; i<NXPoints; i++) {
                    float XCenterOfSubPixel = xcenter + InstrumentProfile->getIPixXCoordinate(i);				
                    unsigned xx = (unsigned)floor(XCenterOfSubPixel); 
                    if (xx > 0 && xx < masterCompImage.getnaxis1() &&
                        yy > 0 && yy < masterCompImage.getnaxis2() &&
                        masterCompImage[yy][xx] < SATURATIONLIMIT && 
                        badpix[yy][xx] == 1 && 
                        (float)masterCompImage[yy][xx] > 0 ) {                    
                        my +=  (float)(masterCompImage[yy][xx]  - bias[yy][xx]) * InstrumentProfile->getipDataFromPolyModel(distdElem,i,j);
                        myvar += (noise/gain)*(noise/gain)/float(NXPoints*NYPoints) + fabs((float)masterCompImage[yy][xx] * InstrumentProfile->getipDataFromPolyModel(distdElem,i,j));
                    } else {
                        fluxFractionLost += InstrumentProfile->getipDataFromPolyModel(distdElem,i,j);
                    }
                    IPNormalizationFactor += InstrumentProfile->getipDataFromPolyModel(distdElem,i,j); //In case IP is not properly normalized. This could occurs after the polynomial fit
                }
            }      
            
            if(my && IPNormalizationFactor) {
                my = (my/IPNormalizationFactor)*(1+fluxFractionLost/IPNormalizationFactor);
                myvar = (myvar/IPNormalizationFactor)*(1+fluxFractionLost/IPNormalizationFactor);
            } else {
                my = NAN;
                myvar = NAN;
            }
            
            SpectralElements->setFlux(my,indexElem);
            SpectralElements->setFluxVariance(myvar,indexElem);
            
#ifdef PRINT_DEBUG	
            //print out extracted cross-correlation spectrum 
            cout << SpectralElements->getphotoCenterY(indexElem) <<
            "\t" << SpectralElements->getFlux(indexElem) <<
            "\t" << sqrt(SpectralElements->getFluxVariance(indexElem)) <<
            "\t" << getorder() << endl;			
#endif	        
        }     
        setSpectrumType(RawSpectrum);   
        SpectralElements->setSpectrumType(RawSpectrum);   
        SpectralElements->setHasRawSpectrum(true);   
    }
    if (SpectralLines) {
		delete SpectralLines;
	}
    SpectralLines = new operaSpectralLines(SpectralElements,(double)ReferenceLineWidth,y_distance_disp);
    
    SpectralLines->detectSpectralFeatures((double)DetectionThreshold,(double)LocalMaxFilterWidth,(double)MinPeakDepth); 
    
    if(SpectralLines->getNFeatures() == 0) { 
        sethasSpectralLines(false);
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    } else {
        sethasSpectralLines(true);	// Note that if the exception was thrown then hasspectrallines is false
    }
}

void operaSpectralOrder::calculateXCorrBetweenIPandImage(operaFITSImage &Image, operaFITSImage &badpix, ostream *pout) {
    if(!gethasInstrumentProfile()) {
        throw operaException("operaSpectralOrder::calculateXCorrBetweenIPandImage: ",operaErrorHasNoInstrumentProfile, __FILE__, __FUNCTION__, __LINE__);
    }
    if(!gethasSpectralElements()) {
        throw operaException("operaSpectralOrder::calculateXCorrBetweenIPandImage: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }    
    
    const unsigned NXPoints = InstrumentProfile->getNXPoints();		
    const unsigned NYPoints = InstrumentProfile->getNYPoints();
    const unsigned nMaxDataPoints = SpectralElements->getnSpectralElements();
    
    for(unsigned indexElem=0;indexElem < nMaxDataPoints; indexElem++) {    
        float xcenter =  SpectralElements->getphotoCenterX(indexElem);
        float ycenter =  SpectralElements->getphotoCenterY(indexElem);
        float distdElem = SpectralElements->getdistd(indexElem);
		
		VLArray<int> xCoords(NXPoints);
        VLArray<int> yCoords(NYPoints);
		for (unsigned j=0; j<NYPoints; j++) yCoords(j) = (int)(ycenter + InstrumentProfile->getIPixYCoordinate(j));
		for (unsigned i=0; i<NXPoints; i++) xCoords(i) = (int)(xcenter + InstrumentProfile->getIPixXCoordinate(i));
        
        unsigned iMin, jMin, iMax, jMax;
        for (iMin=0; iMin < NXPoints; iMin++) if(xCoords(iMin) > 0) break; //find lowest index where xCoords[i] > 0
        for (jMin=0; jMin < NYPoints; jMin++) if(yCoords(jMin) > 0) break; //find lowest index where yCoords[j] > 0
		for (iMax=NXPoints; iMax > 0; iMax--) if(xCoords(iMax-1) < Image.getnaxis1()) break; //find highest index where xCoords[i] < naxis1
		for (jMax=NYPoints; jMax > 0; jMax--) if(yCoords(jMax-1) < Image.getnaxis2()) break; //find highest index where yCoords[j] < naxis2
		
		VLArray<bool> validCoords(NYPoints, NXPoints);
		for (unsigned j=jMin; j<jMax; j++) {	
            for (unsigned i=iMin; i<iMax; i++) {
				validCoords(j, i) = Image[yCoords(j)][xCoords(i)] > 0 && Image[yCoords(j)][xCoords(i)] < SATURATIONLIMIT && badpix[yCoords(j)][xCoords(i)];
			}
		}
        
        VLArray<float> ipvals(NYPoints, NXPoints);
		float avgImg = 0, meanIP = 0;
        unsigned npImg = 0;
		for (unsigned j=jMin; j<jMax; j++) {	
            for (unsigned i=iMin; i<iMax; i++) {
                if (validCoords(j, i)) {                    
                    avgImg += (float)Image[yCoords(j)][xCoords(i)];
                    ipvals(j, i) = InstrumentProfile->getipDataFromPolyModel(distdElem,i,j);
                    meanIP += ipvals(j, i);
                    npImg++;
                }
            }
        }
        avgImg /= (float)npImg;
        meanIP /= (float)npImg;
        
        float Xcorr = 0, imgsqr = 0, ipsqr = 0;
		for (unsigned j=jMin; j<jMax; j++) {	
            for (unsigned i=iMin; i<iMax; i++) {
                if (validCoords(j, i)) {
					const float ip = ipvals(j, i) - meanIP;
					const float imgval = (float)Image[yCoords(j)][xCoords(i)] - avgImg;
                    Xcorr +=  imgval * ip;
                    imgsqr += imgval * imgval;
                    ipsqr += ip * ip;
                }
            }
        }      
        
        if(Xcorr) Xcorr /= sqrt(imgsqr*ipsqr);
        else Xcorr = NAN;
        SpectralElements->setXCorrelation((double)Xcorr,indexElem);
        
        if (pout != NULL) {
            //print out cross-correlation 
            *pout << getorder() << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "            
            << SpectralElements->getXCorrelation(indexElem) << " "
            << endl;			
        }
    }
    SpectralElements->setHasXCorrelation(true);
}

void operaSpectralOrder::measureInstrumentProfile(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, const int minimumLines) {
	
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfile: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
#ifdef PRINT_DEBUG	    
    for(unsigned k=0; k<nSelectedLines; k++) {   
        cout << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
    }   
#endif    
	
    unsigned NXPoints = InstrumentProfile->getNXPoints();		
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    unsigned cubeDataIndex = 0;
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
        
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = (float)Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        
        InstrumentProfile->setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,cubeDataIndex);        
		
        InstrumentProfile->subtractOuterFrame(cubeDataIndex);
		
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex))){
                    SumOfUsefulFluxWithinLine += InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
        
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float ipmatrixvalue = 0;
                if(NormalizationFactor && !isnan(InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex))){
                    ipmatrixvalue = InstrumentProfile->getdataCubeValues(i,j,cubeDataIndex)/NormalizationFactor;
                } else {
                    ipmatrixvalue = InstrumentProfile->getipDataFromPolyModel(distd,i,j);  
                }
                InstrumentProfile->setdataCubeValues(ipmatrixvalue,i,j,cubeDataIndex);
            }
        }  
        
        InstrumentProfile->normalizeCubeData(cubeDataIndex);        
		
        InstrumentProfile->setdistd(distd,cubeDataIndex);            
        
        cubeDataIndex++;  
		
    }
    
    InstrumentProfile->setnDataPoints(cubeDataIndex);  
    
    if(cubeDataIndex < minimumLines) {
#ifdef PRINT_DEBUG
		cerr << "operaSpectralOrder::measureInstrumentProfile: Order " << getorder() << " IP measurements not possible." << endl;
#endif
		sethasInstrumentProfile(false);
    } else if(cubeDataIndex >= minimumLines && cubeDataIndex <= 3*minimumLines) {
        InstrumentProfile->FitMediantoIPDataVector();
    } else if (cubeDataIndex > 3*minimumLines) {
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);        
    }    
	
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
    
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::measureInstrumentProfileWithBinning(operaFITSImage &masterCompImage, operaFITSImage &badpix, double binsize, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, const int minimumLines) {
    
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileWithBinning: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
    if(nSelectedLines == 0) {
		delete[] LinePositionVector;
		delete[] LineSigmaVector;
		delete[] LineAmplitudeVector;
        throw operaException("operaSpectralOrder::measureInstrumentProfileWithBinning: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
#ifdef PRINT_DEBUG	    
    for(unsigned k=0; k<nSelectedLines; k++) {   
        cout << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
    }   
#endif      
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();		
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    unsigned cubeDataIndex = 0;   
    
    unsigned ipxsize = InstrumentProfile->getxsize();	
    unsigned ipysize = InstrumentProfile->getysize();	
    unsigned ipxsampling = InstrumentProfile->getXsampling();  
    unsigned ipysampling = InstrumentProfile->getYsampling();  
    
    operaInstrumentProfile tempIP(ipxsize,ipxsampling,ipysize,ipysampling,nSelectedLines);
    float dlimit = binsize;
    
    unsigned npInBin = 0;
    float averageDistd = 0;
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        averageDistd += distd;
        
        tempIP.setdistd(distd,npInBin);
        
        tempIP.setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,npInBin);
        tempIP.subtractOuterFrame(npInBin);
        
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(tempIP.getdataCubeValues(i,j,npInBin))){
                    SumOfUsefulFluxWithinLine += tempIP.getdataCubeValues(i,j,npInBin);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
        
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float ipmatrixvalue;
                if(NormalizationFactor && !isnan(tempIP.getdataCubeValues(i,j,npInBin))){
                    ipmatrixvalue = tempIP.getdataCubeValues(i,j,npInBin)/NormalizationFactor;
                } else {
                    ipmatrixvalue = InstrumentProfile->getipDataFromPolyModel(distd,i,j);  
                }
                tempIP.setdataCubeValues(ipmatrixvalue,i,j,npInBin);
            }
        }  
        
        tempIP.normalizeCubeData(npInBin);        
        
        npInBin++;
        
        if((npInBin >= 3 && distd > dlimit) || k == nSelectedLines-1)  {
            
			tempIP.setnDataPoints(npInBin);            
            tempIP.FitMediantoIPDataVector();            
            averageDistd /= (float)npInBin;
            
            InstrumentProfile->setdistd(averageDistd,cubeDataIndex);
            
            for (unsigned j=0; j<NYPoints; j++) {
                for (unsigned i=0; i<NXPoints; i++) {
                    float ipmatrixvalue = tempIP.getipDataFromPolyModel(averageDistd,i,j);
                    InstrumentProfile->setdataCubeValues(ipmatrixvalue,i,j,cubeDataIndex);
                }
            }
            InstrumentProfile->normalizeCubeData(cubeDataIndex); 
			
            dlimit += binsize;
            averageDistd = 0;
            cubeDataIndex++;
            npInBin = 0;
        } else if (npInBin < 3 && distd > dlimit) {
            dlimit += binsize;  
        }      
    }
    
    InstrumentProfile->setnDataPoints(cubeDataIndex);        
    
    if(cubeDataIndex < minimumLines) {
#ifdef PRINT_DEBUG	    
		cerr << "operaSpectralOrder::measureInstrumentProfileWithBinning: Order " << getorder() << " IP measurements not possible." << endl;
#endif
        sethasInstrumentProfile(false);
    } else if(cubeDataIndex >= minimumLines) {
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);        
    }    
    
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
	
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::measureInstrumentProfileUsingMedian(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, const int minimumLines) {
    
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileUsingMedian: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
#ifdef PRINT_DEBUG
    if(nSelectedLines == 0) { 
		cerr << "operaSpectralOrder::measureInstrumentProfileUsingMedian: Order " << getorder() << " no spectral lines." << endl;
    }
#endif      
#ifdef PRINT_DEBUG	    
    for(unsigned k=0; k<nSelectedLines; k++) {   
        cout << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
    }   
#endif      
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();		
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        
        InstrumentProfile->setdistd(distd,k);
        InstrumentProfile->setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,k);
        InstrumentProfile->subtractOuterFrame(k);
        
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;  
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    SumOfUsefulFluxWithinLine += InstrumentProfile->getdataCubeValues(i,j,k);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
        
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
		
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float ipmatrixvalue = 0;
                if(NormalizationFactor && !isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    ipmatrixvalue = InstrumentProfile->getdataCubeValues(i,j,k)/NormalizationFactor;
                } else {
                    ipmatrixvalue = InstrumentProfile->getipDataFromPolyModel(distd,i,j);  
                }
                InstrumentProfile->setdataCubeValues(ipmatrixvalue,i,j,k);
            }
        }  
        
        InstrumentProfile->normalizeCubeData(k);        
    }
	
    if(nSelectedLines > minimumLines) {
        InstrumentProfile->FitMediantoIPDataVector();
        sethasInstrumentProfile(true);
    } else {
        sethasInstrumentProfile(false);
    }
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
	
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::measureInstrumentProfileUsingWeightedMean(operaFITSImage &masterCompImage, operaFITSImage &badpix, double MaxContamination, double amplitudeCutOff, unsigned nSigCut, unsigned sampleElementForPlot, ostream *pout, const int minimumLines) {
    
    if(!gethasSpectralLines()) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileUsingWeightedMean: ",operaErrorNoSpectralLineDetected, __FILE__, __FUNCTION__, __LINE__);
    }
    if(SpectralLines->getnLines() == 0) {
        throw operaException("operaSpectralOrder::measureInstrumentProfileUsingWeightedMean: ",operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double *LinePositionVector = new double[SpectralLines->getnLines()];
    double *LineSigmaVector = new double[SpectralLines->getnLines()];
    double *LineAmplitudeVector = new double[SpectralLines->getnLines()];
    
    unsigned nSelectedLines = SpectralLines->selectLines(MaxContamination,nSigCut,amplitudeCutOff,LinePositionVector,LineSigmaVector,LineAmplitudeVector);
    
    double TotalWeight = 0;  
    
    for(unsigned k=0; k<nSelectedLines; k++) {  
        TotalWeight += LineAmplitudeVector[k];
    }  
    
#ifdef PRINT_DEBUG	    
    for(unsigned k=0; k<nSelectedLines; k++) {   
        cout << LinePositionVector[k] << " " << LineAmplitudeVector[k] << " " << LineSigmaVector[k] << endl;
    }   
#endif      
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    unsigned ipxsize = InstrumentProfile->getxsize();
    unsigned ipysize = InstrumentProfile->getysize();
    unsigned ipxsampling = InstrumentProfile->getXsampling();
    unsigned ipysampling = InstrumentProfile->getYsampling();
    
    operaInstrumentProfile tempIP(ipxsize,ipxsampling,ipysize,ipysampling,1);
    
    for (unsigned j=0; j<NYPoints; j++) {	 
        for (unsigned i=0; i<NXPoints; i++) {
            tempIP.setdataCubeValues(0.0,i,j,0);
        }
    }
    
    for(unsigned k=0; k<nSelectedLines; k++) {
        
        float ycenter =  LinePositionVector[k];
		float xcenter = (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        
		float distd = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);        
        
        InstrumentProfile->setdistd(distd,k);
        InstrumentProfile->setdataCubeValues(masterCompImage,badpix,xcenter,ycenter,k);
        InstrumentProfile->subtractOuterFrame(k);
        
        float FractionOfFluxLost = 0;
        float SumOfUsefulFluxWithinLine = 0;  
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                if(!isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    SumOfUsefulFluxWithinLine += InstrumentProfile->getdataCubeValues(i,j,k);                    
                } else {
                    FractionOfFluxLost += InstrumentProfile->getipDataFromPolyModel(distd,i,j);
                }  
            }
        }   
		
        float NormalizationFactor = 0;
        
        if(FractionOfFluxLost<1 && SumOfUsefulFluxWithinLine)
            NormalizationFactor = SumOfUsefulFluxWithinLine/(1.0 - FractionOfFluxLost);
        
        for (unsigned j=0; j<NYPoints; j++) {	 
            for (unsigned i=0; i<NXPoints; i++) {
                float tempipmatrixvalue = tempIP.getdataCubeValues(i,j,0);
                if(NormalizationFactor && !isnan(InstrumentProfile->getdataCubeValues(i,j,k))){
                    tempipmatrixvalue += (InstrumentProfile->getdataCubeValues(i,j,k)/NormalizationFactor)*(float)(LineAmplitudeVector[k]/TotalWeight);
                } else {
                    tempipmatrixvalue += InstrumentProfile->getipDataFromPolyModel(distd,i,j)*(float)(LineAmplitudeVector[k]/TotalWeight);
                }
                tempIP.setdataCubeValues(tempipmatrixvalue,i,j,0);
            }
        }  
    }
    
    InstrumentProfile->setnDataPoints(1);
    
    for (unsigned j=0; j<NYPoints; j++) {
        for (unsigned i=0; i<NXPoints; i++) {
            InstrumentProfile->setdataCubeValues(tempIP.getdataCubeValues(i,j,0),i,j,0);
        }
    }
    
    unsigned zero = 0;
    InstrumentProfile->normalizeCubeData(zero);
    
    InstrumentProfile->FitPolyMatrixtoIPDataVector(1,false); 
    
	if (pout != NULL) {
        InstrumentProfile->printModel(SpectralElements->getdistd(sampleElementForPlot),orderNumber,pout);
	}
    
    delete[] LinePositionVector;
    delete[] LineSigmaVector;
    delete[] LineAmplitudeVector;
}

void operaSpectralOrder::recenterOrderPosition(void) {
	
    for(unsigned i=0; i<Geometry->getNdatapoints(); i++) {    
        
        float ycenter =  Geometry->getCenterY(i);
        float xcenter =  (float)Geometry->getCenterPolynomial()->Evaluate(double(ycenter));
        float d = Geometry->CalculateDistance(Geometry->getYmin(), ycenter);
        
        // Set error as the size of an IP sub-pixel. In fact the error can be calculated,
        // however geometry is not using errors to trace order position, so we leave like this for now.
        
        double xerror = 1.0/(double)InstrumentProfile->getXsampling(); 
		
        double x = (double)(xcenter + InstrumentProfile->getIPphotoCenterX(d));
        
        double y = (double)(ycenter + InstrumentProfile->getIPphotoCenterY(d));
        
        Geometry->resetCenter(x,y,xerror,i);
    }
    
    double chisqr;
    bool witherrors = false;
    unsigned maxOrderofTracePolynomial = 3;	// was 7, too high... DT Apr 25 2013
    
    Geometry->traceOrder(maxOrderofTracePolynomial, chisqr, witherrors);
	
}

/*
 * Print out beam spectra. Useful for plotting.
 */

void operaSpectralOrder::printBeamSpectrum(ostream *pout) {
    if (pout != NULL && gethasSpectralElements()) {
        unsigned NumberofElements = SpectralElements->getnSpectralElements();
        
        for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
            *pout << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getXCorrelation(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *pout << BeamElements[beam]->getphotoCenterX(indexElem) << " "
                << BeamElements[beam]->getphotoCenterY(indexElem) << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " ";
            }
            *pout << endl;
        }
        *pout << endl;
    }
}

void operaSpectralOrder::printBeamSpectrum(string addFirstColumnEntry, ostream *pout) {
    if (pout != NULL && gethasSpectralElements()) {
        unsigned NumberofElements = SpectralElements->getnSpectralElements();
        
        for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
            *pout << addFirstColumnEntry << " "
            << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getXCorrelation(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *pout << BeamElements[beam]->getphotoCenterX(indexElem) << " "
                << BeamElements[beam]->getphotoCenterY(indexElem) << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " ";
            }
            *pout << endl;
        }
        *pout << endl;
    }
}

/*
 * Extract Raw Spectrum - this function uses extraction apertures.
 */
void operaSpectralOrder::extractRawSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, ostream *pout) {
    
	//    operaFITSImage objectCloneImg("clone.fits", objectImage.getnaxis1(), objectImage.getnaxis2(), tfloat, READWRITE, cNone);
	//    objectCloneImg.operaFITSImageCopyHeader(&objectImage);
	//    objectCloneImg = objectImage;
	
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::extractRawSpectrum: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
	
    // Set master spectral elements.
    setSpectralElementsByStitchingApertures(effectiveApertureFraction);
    
    // Set spectral elements for each beam aperture and for the background apertures
    setApertureElements(SpectralElements->getSpectrumType());
    
	unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
    PixelSet *aperturePixels[MAXNUMBEROFBEAMS];
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        aperturePixels[beam] = ExtractionApertures[beam]->getSubpixels();
    }
    
    for (unsigned indexElem=0; indexElem < NumberofElements; indexElem++) {
        // get elem center coordinates
        double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
        double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
        
        double objRawFluxSum = 0;
        double objRawFluxSumVar = 0;
        unsigned NTotalUsefulPoints = 0;
        unsigned NTotalPoints = 0;
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double subpixelArea = (double)(aperturePixels[beam]->getSubpixelArea());
            
            double objBeamRawFluxSum = 0;
            double objBeamFluxSumVar = 0;
            unsigned NPointsInBeam = 0;
            
            // loop over all subpixels within the beam aperture to calculate the sum of subpixel raw values
            for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                // select image col and row of subpixel
                unsigned col = (unsigned)floor(elemXcenter + aperturePixels[beam]->getXcenter(pix));
                unsigned row = (unsigned)floor(elemYcenter + aperturePixels[beam]->getYcenter(pix));
                
                if(col >= 0 && row >= 0 && col < objectImage.getnaxis1() && row < objectImage.getnaxis2()
                   && objectImage[row][col] < SATURATIONLIMIT && badpix[row][col] == 1 && nflatImage[row][col]) {
                    
                    double gain = gainBiasNoise.getGain(col,row);
                    double noise = gainBiasNoise.getNoise(col,row);
                    
                    double objectRawFlux = (double)(gain*(objectImage[row][col] - biasImage[row][col])/nflatImage[row][col]);
                    double objRawFluxVar = 2.0*noise*noise*sqrt(subpixelArea) + fabs(objectRawFlux); // detector noise + photon noise ;
                    
                    //   objectCloneImg[row][col] = objectCloneImg[row][col] - (objectRawFlux/gain)*subpixelArea;
                    
                    //below it calculates the raw (object+background) flux sum in e-/pxl^2
                    objBeamRawFluxSum += objectRawFlux;
                    objBeamFluxSumVar += objRawFluxVar;
                    NPointsInBeam++;
                }
            }
            
            objRawFluxSum += objBeamRawFluxSum*subpixelArea;
            objRawFluxSumVar += objBeamFluxSumVar*subpixelArea;
            NTotalUsefulPoints += NPointsInBeam;
            NTotalPoints += aperturePixels[beam]->getNPixels();
            
            if(NPointsInBeam==0) {
                objBeamRawFluxSum = NAN;
                objBeamFluxSumVar = NAN;
            } else {
                objBeamRawFluxSum *= subpixelArea * (double)aperturePixels[beam]->getNPixels()/(double)NPointsInBeam;
                objBeamFluxSumVar *= subpixelArea * (double)aperturePixels[beam]->getNPixels()/(double)NPointsInBeam;
            }
            
            BeamElements[beam]->setFlux(objBeamRawFluxSum,indexElem);
            BeamElements[beam]->setFluxVariance(objBeamFluxSumVar,indexElem);
        }
        
        if(NTotalUsefulPoints==0) {
            objRawFluxSum = NAN;
            objRawFluxSumVar = NAN;
        } else {
            objRawFluxSum *= (double)NTotalPoints/(double)NTotalUsefulPoints;
            objRawFluxSumVar *= (double)NTotalPoints/(double)NTotalUsefulPoints;
        }
        
        SpectralElements->setFlux(objRawFluxSum,indexElem);
        SpectralElements->setFluxVariance(objRawFluxSumVar,indexElem);
        
	}
	
    printBeamSpectrum(pout);
    
    setSpectrumType(RawBeamSpectrum);
    SpectralElements->setSpectrumType(RawBeamSpectrum);
    SpectralElements->setHasRawSpectrum(true);
    sethasSpectralElements(true);

	//    objectCloneImg.operaFITSImageSave();
	//    objectCloneImg.operaFITSImageClose();
}
/*
 * Extract Standard Spectrum - this function uses extraction apertures.
 */
void operaSpectralOrder::extractStandardSpectrum(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, ostream *pout) {
    
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::extractStandardSpectrum: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    
    // Set master spectral elements.
    setSpectralElementsByStitchingApertures(effectiveApertureFraction);
    
    // Set spectral elements for each beam aperture and for the background apertures
    setApertureElements(SpectralElements->getSpectrumType());
    
	unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
#ifdef PRINT_OUTOFBOUNDS
	if (NumberofElements == 0)
		cerr << "operaSpectralOrder::extractStandardSpectrum Warning: NumberofElements (" << NumberofElements << ") == 0" << endl;
#endif
    
    PixelSet *aperturePixels[MAXNUMBEROFBEAMS];
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        aperturePixels[beam] = ExtractionApertures[beam]->getSubpixels();
    }
    
    PixelSet *backgroundPixels[LEFTANDRIGHT];
    unsigned MAXNPointsInBkg = 0;
    for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
        backgroundPixels[backgroundIndex]  = BackgroundApertures[backgroundIndex]->getSubpixels();
        MAXNPointsInBkg += backgroundPixels[backgroundIndex]->getmaxNPixels();
    }
    
    float *fBackgroundFlux = new float[BackgroundBinsize*MAXNPointsInBkg+1];
    
#ifdef PRINT_OUTOFBOUNDS
	if (BackgroundBinsize == 0)
		cerr << "operaSpectralOrder::extractStandardSpectrum Warning: BackgroundBinsize (" << BackgroundBinsize << ") == 0" << endl;
	if (MAXNPointsInBkg == 0)
		cerr << "operaSpectralOrder::extractStandardSpectrum Warning: MAXNPointsInBkg (" << MAXNPointsInBkg << ") == 0" << endl;
#endif
	unsigned NumberOfElementSamples;
	unsigned NumberofElementsToBin = BackgroundBinsize;
	
	NumberOfElementSamples = (unsigned)ceil((float)SpectralElements->getnSpectralElements()/(float)NumberofElementsToBin);
    
#ifdef PRINT_OUTOFBOUNDS
	if (NumberOfElementSamples == 0)
		cerr << "operaSpectralOrder::extractStandardSpectrum Warning: NumberOfElementSamples (" << NumberOfElementSamples << ") == 0" << endl;
#endif
	float *BackgroundDistd = new float[NumberOfElementSamples];
    float *BackgroundFlux = new float[NumberOfElementSamples];
    
    double BackgroundFluxVar=0;
    unsigned numberOfValidPoints=0;
	
	for(unsigned k=0;k<NumberOfElementSamples;k++){
		
		unsigned firstElement = NumberofElementsToBin*(k);
		unsigned lastElement =  NumberofElementsToBin*(k+1);
		if (lastElement > SpectralElements->getnSpectralElements()){
            lastElement = SpectralElements->getnSpectralElements();
        }
		
        // Loop over background elements
        unsigned NPointsInBkg = 0;
        for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
            // get elem center coordinates
            double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
            double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
            
            for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
                //  double subpixelArea = (double)(backgroundPixels[backgroundIndex]->getSubpixelArea());
                
                for(unsigned pix=0; pix<backgroundPixels[backgroundIndex]->getNPixels(); pix++) {
                    unsigned col = (unsigned)floor(elemXcenter + backgroundPixels[backgroundIndex]->getXcenter(pix));
                    unsigned row = (unsigned)floor(elemYcenter + backgroundPixels[backgroundIndex]->getYcenter(pix));
                    
                    if(col >= 0 && row >= 0 && col < objectImage.getnaxis1() && row < objectImage.getnaxis2()
                       && objectImage[row][col] < SATURATIONLIMIT && badpix[row][col] == 1 && nflatImage[row][col]
                       && NPointsInBkg < BackgroundBinsize*MAXNPointsInBkg+1) {
                        
                        double gain = gainBiasNoise.getGain(col,row);
                        fBackgroundFlux[NPointsInBkg] = (float)(gain*(objectImage[row][col] - biasImage[row][col])/nflatImage[row][col]);
                        NPointsInBkg++;
                        
#ifdef PRINT_OUTOFBOUNDS
                        else {
                            cerr << "operaSpectralOrder::extractStandardSpectrum Warning: NPointsInBkg (" << NPointsInBkg << ") >= BackgroundBinsize*MAXNPointsInBkg(" << (BackgroundBinsize*MAXNPointsInBkg) << endl;;
                        }
#endif
                    }
                }
            }
        }
        
		// DT 7/26/2012 NOTE Change to use MedianQuick, since the values are never used again...
        if (NPointsInBkg > 0 && numberOfValidPoints <= NumberOfElementSamples) {
            float distdmidElem = (float)fabs(SpectralElements->getdistd(lastElement-1) + SpectralElements->getdistd(firstElement))/2;
            BackgroundDistd[numberOfValidPoints] = distdmidElem;
            BackgroundFlux[numberOfValidPoints] = operaArrayMedianQuick(NPointsInBkg,fBackgroundFlux);
            float BackgroundFluxError = operaArrayMedianSigmaQuick(NPointsInBkg,fBackgroundFlux,BackgroundFlux[numberOfValidPoints]);
            BackgroundFluxVar += (double)(BackgroundFluxError*BackgroundFluxError);
            numberOfValidPoints++;
        }
#ifdef PRINT_OUTOFBOUNDS
		else {
            if (NPointsInBkg == 0)
				cerr << "operaSpectralOrder::extractStandardSpectrumWarning: NPointsInBkg (" << NPointsInBkg << ") == 0" << endl;
			else
				cerr << "operaSpectralOrder::extractStandardSpectrumWarning: numberOfValidPoints (" << numberOfValidPoints << ") >= NumberOfElementSamples(" << NumberOfElementSamples <<")" << endl;
		}
#endif
	}
	
    // the model below gives the background in units of e-/pxl^2
    float *BackgroundModelDistd = new float[NumberofElements];
    float *BackgroundModelFlux = new float[NumberofElements];
    
    for (unsigned indexElem=0; indexElem < NumberofElements; indexElem++) {
        BackgroundModelDistd[indexElem] = SpectralElements->getdistd(indexElem);
        BackgroundModelFlux[indexElem] = 0.0;
    }
    
    if(numberOfValidPoints > 0) {
        BackgroundFluxVar /= (double)numberOfValidPoints;
        operaFitSpline(numberOfValidPoints,BackgroundDistd,BackgroundFlux,NumberofElements,BackgroundModelDistd,BackgroundModelFlux);
    } else {
        BackgroundFluxVar = 0;
    }
    
	/* Uncomment below to check background elements
	 for (unsigned indexElem=0; indexElem < NumberofElements; indexElem++) {
	 cout << indexElem << " " << BackgroundModelFlux[indexElem] << " " << BackgroundFluxVar << endl;
	 }
	 exit(EXIT_SUCCESS);
	 */
    
    for (unsigned indexElem=0; indexElem < NumberofElements; indexElem++) {
        // get elem center coordinates
        double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
        double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
        
        for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
            BackgroundElements[backgroundIndex]->setFlux(BackgroundModelFlux[indexElem],indexElem);
            BackgroundElements[backgroundIndex]->setFluxVariance(BackgroundFluxVar,indexElem);
        }
        
        double objStandardFlux = 0;
        double objStandardFluxVar = 0;
        unsigned NTotalUsefulPoints = 0;
        unsigned NTotalPoints = 0;
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double subpixelArea = (double)(aperturePixels[beam]->getSubpixelArea());
            
            double objBeamStandardFlux = 0;
            double objBeamStandardFluxVar = 0;
            unsigned NPointsInBeam = 0;
            
            // loop over all subpixels within the beam aperture to calculate the sum of subpixel raw values
            for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                // select image col and row of subpixel
                unsigned col = (unsigned)floor(elemXcenter + aperturePixels[beam]->getXcenter(pix));
                unsigned row = (unsigned)floor(elemYcenter + aperturePixels[beam]->getYcenter(pix));
                
                if(col >= 0 && row >= 0 && col < objectImage.getnaxis1() && row < objectImage.getnaxis2()
                   && objectImage[row][col] < SATURATIONLIMIT && badpix[row][col] == 1 && nflatImage[row][col]) {
                    
                    double gain = gainBiasNoise.getGain(col,row);
                    double noise = gainBiasNoise.getNoise(col,row);
                    //below is the measured object+background flux in e-/pxl^2
                    double objectRawFlux = (double)(gain*(objectImage[row][col] - biasImage[row][col])/nflatImage[row][col]);
                    double objRawFluxVar = (2.0*noise*noise + (double)BackgroundFluxVar)*sqrt(subpixelArea) + fabs(objectRawFlux); // detector noise + background noise + photon noise ;
                    
                    //below is the raw (object+background) flux sum in e-/pxl^2
                    objBeamStandardFlux += (objectRawFlux - (double)BackgroundModelFlux[indexElem]);
                    objBeamStandardFluxVar += objRawFluxVar;
                    
                    NPointsInBeam++;
                }
            }
            
            objStandardFlux += objBeamStandardFlux*subpixelArea;
            objStandardFluxVar += objBeamStandardFluxVar*subpixelArea;
            NTotalUsefulPoints += NPointsInBeam;
            NTotalPoints += aperturePixels[beam]->getNPixels();
            
            if(NPointsInBeam==0) {
                objBeamStandardFlux = NAN;
                objBeamStandardFluxVar = NAN;
            } else {
                objBeamStandardFlux *= subpixelArea * (double)aperturePixels[beam]->getNPixels()/(double)NPointsInBeam;
                objBeamStandardFluxVar *= subpixelArea * (double)aperturePixels[beam]->getNPixels()/(double)NPointsInBeam;
            }
            
            BeamElements[beam]->setFlux(objBeamStandardFlux,indexElem);
            BeamElements[beam]->setFluxVariance(objBeamStandardFluxVar,indexElem);
        }
        if(NTotalUsefulPoints==0) {
            objStandardFlux = NAN;
            objStandardFluxVar = NAN;
        } else {
            objStandardFlux *= (double)NTotalPoints/(double)NTotalUsefulPoints;
            objStandardFluxVar *= (double)NTotalPoints/(double)NTotalUsefulPoints;
        }
        SpectralElements->setFlux(objStandardFlux,indexElem);
        SpectralElements->setFluxVariance(objStandardFluxVar,indexElem);
	}
    
    printBeamSpectrum(pout);
	
	// DT Sept 4 2012 - free malloc'd storage
	if (NumberofElements > 0) {
		if (BackgroundModelDistd) {
			delete[] BackgroundModelDistd;
			BackgroundModelDistd = NULL;
		}
		if (BackgroundModelFlux) {
			delete[] BackgroundModelFlux;
			BackgroundModelFlux = NULL;
		}
	}
	if (NumberOfElementSamples > 0) {
		if (BackgroundDistd) {
			delete[] BackgroundDistd;
			BackgroundDistd = NULL;
		}
		if (BackgroundFlux) {
			delete[] BackgroundFlux;
			BackgroundFlux = NULL;
		}
	}
	if (BackgroundBinsize*MAXNPointsInBkg > 0) {
		if (fBackgroundFlux) {
			delete[] fBackgroundFlux;
			fBackgroundFlux = NULL;
		}
	}
    
	setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setHasOptimalSpectrum(true);
    sethasSpectralElements(true);
}

/*
 * Extract Standard Spectrum without background subtraction - this function is basically the
 * raw extraction but it also creates the dummy background beams so optimal extraction won't
 * get confused.
 */
void operaSpectralOrder::extractStandardSpectrumNoBackground(operaFITSImage &objectImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, ostream *pout) {
    
    if(!gethasExtractionApertures()) {
        throw operaException("operaSpectralOrder::extractStandardSpectrumNoBackground: ",operaErrorHasNoExtractionAperture, __FILE__, __FUNCTION__, __LINE__);
    }
    // Set master spectral elements.
    setSpectralElementsByStitchingApertures(effectiveApertureFraction);
    // Set spectral elements for each beam aperture and for the background apertures
    setApertureElements(SpectralElements->getSpectrumType());
    
	unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
    PixelSet *aperturePixels[MAXNUMBEROFBEAMS];
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        aperturePixels[beam] = ExtractionApertures[beam]->getSubpixels();
    }
	
    for (unsigned indexElem=0; indexElem < NumberofElements; indexElem++) {
        for(unsigned backgroundIndex = 0; backgroundIndex < LEFTANDRIGHT; backgroundIndex++) {
            BackgroundElements[backgroundIndex]->setFlux(0.0,indexElem);
            BackgroundElements[backgroundIndex]->setFluxVariance(0.0,indexElem);
        }
        
        // get elem center coordinates
        double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
        double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
        
        double objStandardFlux = 0;
        double objStandardFluxVar = 0;
        unsigned NTotalUsefulPoints = 0;
        unsigned NTotalPoints = 0;
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            double subpixelArea = (double)(aperturePixels[beam]->getSubpixelArea());
            
            double objBeamStandardFlux = 0;
            double objBeamStandardFluxVar = 0;
            unsigned NPointsInBeam = 0;
            
            // loop over all subpixels within the beam aperture to calculate the sum of subpixel raw values
            for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                // select image col and row of subpixel
                unsigned col = (unsigned)floor(elemXcenter + aperturePixels[beam]->getXcenter(pix));
                unsigned row = (unsigned)floor(elemYcenter + aperturePixels[beam]->getYcenter(pix));
                
                if(col >= 0 && row >= 0 && col < objectImage.getnaxis1() && row < objectImage.getnaxis2()
                   && objectImage[row][col] < SATURATIONLIMIT && badpix[row][col] == 1 && nflatImage[row][col]) {
                    
                    double gain = gainBiasNoise.getGain(col,row);
                    double noise = gainBiasNoise.getNoise(col,row);
                    //below is the measured object+background flux in e-/pxl^2
                    double objectRawFlux = (double)(gain*(objectImage[row][col] - biasImage[row][col])/nflatImage[row][col]);
                    double objRawFluxVar = 2.0*noise*noise*sqrt(subpixelArea) + fabs(objectRawFlux); // detector noise + photon noise + background noise;
                    
                    //below is the raw (object+background) flux sum in e-/pxl^2
                    objBeamStandardFlux += objectRawFlux;
                    objBeamStandardFluxVar += objRawFluxVar;
                    
                    NPointsInBeam++;
                }
            }
            objStandardFlux += objBeamStandardFlux*subpixelArea;
            objStandardFluxVar += objBeamStandardFluxVar*subpixelArea;
            NTotalUsefulPoints += NPointsInBeam;
            NTotalPoints += aperturePixels[beam]->getNPixels();
            
            if(NPointsInBeam==0) {
                objBeamStandardFlux = NAN;
                objBeamStandardFluxVar = NAN;
            } else {
                objBeamStandardFlux *= subpixelArea * (double)aperturePixels[beam]->getNPixels()/(double)NPointsInBeam;
                objBeamStandardFluxVar *= subpixelArea * (double)aperturePixels[beam]->getNPixels()/(double)NPointsInBeam;
            }
            
            BeamElements[beam]->setFlux(objBeamStandardFlux,indexElem);
            BeamElements[beam]->setFluxVariance(objBeamStandardFluxVar,indexElem);
        }
        if(NTotalUsefulPoints==0) {
            objStandardFlux = NAN;
            objStandardFluxVar = NAN;
        } else {
            objStandardFlux *= (double)NTotalPoints/(double)NTotalUsefulPoints;
            objStandardFluxVar *= (double)NTotalPoints/(double)NTotalUsefulPoints;
        }
        SpectralElements->setFlux(objStandardFlux,indexElem);
        SpectralElements->setFluxVariance(objStandardFluxVar,indexElem);
	}
	
	setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setSpectrumType(StandardBeamSpectrum);
    SpectralElements->setHasOptimalSpectrum(true);
    sethasSpectralElements(true);
}

/*
 * Extract Optimal Spectrum - this function uses extraction apertures.
 */
void operaSpectralOrder::extractOptimalSpectrum(operaFITSImage &objectImage, operaFITSImage &flatImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned BackgroundBinsize, unsigned sigmaClip, unsigned iterations, bool onTargetProfile, bool usePolynomialFit, bool removeBackground, bool verbose, bool calculateXCorrelation, ostream *pout) {
    /*
     * Below it extracts the standard spectrum if it doesn't exist. STEPS #1 thru #4 K. Horne, 1986
     */     
	if(verbose)
		cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating optimal spectrum using apertures ..." << endl;
	
	if (!hasSpectralElements) {
		if(removeBackground == TRUE) {
			if(verbose)
				cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating standard spectrum ..." << endl;
			extractStandardSpectrum(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,BackgroundBinsize,NULL);
		} else if (removeBackground == FALSE) {
			if(verbose)
				cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating standard spectrum without background ..." << endl;
			extractStandardSpectrumNoBackground(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,NULL);
		}
	} else { 
		if(!SpectralElements->getHasStandardSpectrum() && removeBackground == TRUE) {
			if(verbose)
				cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating standard spectrum ..." << endl;
			extractStandardSpectrum(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,BackgroundBinsize,NULL);
		} else if (!SpectralElements->getHasStandardSpectrum() && removeBackground == FALSE) {
			if(verbose)
				cerr << "operaSpectralOrder::extractOptimalSpectrum: calculating standard spectrum without background ..." << endl;
			extractStandardSpectrumNoBackground(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,NULL);
		}
	}
    unsigned NumberofElementsToBin = BackgroundBinsize;
    
    if(calculateXCorrelation) {
        calculateXCorrBetweenIPandImage(objectImage,badpix,NULL);
    }
    
#ifdef PRINT_OUTOFBOUNDS
	if (NumberofElementsToBin == 0)
		cerr << "operaSpectralOrder::extractOptimalSpectrum Warning: NumberofElementsToBin (" << NumberofElementsToBin << ") == 0" << endl;
#endif
    /*
     * Construct spatial profile. STEP #5 K. Horne, 1986
     */   
    if(onTargetProfile) {
        if(verbose)
            cerr << "operaSpectralOrder::extractOptimalSpectrum: measuring spatial profile on object image ..." << endl;        
        measureBeamSpatialProfiles(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,usePolynomialFit);
    } else {
        if(verbose)
            cerr << "operaSpectralOrder::extractOptimalSpectrum: measuring spatial profile on flat-field image ..." << endl;        
        measureBeamSpatialProfiles(flatImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,usePolynomialFit);
    }
    /*
     * The method measureOptimalSpectrum executes the following operations:
     * 1. Revise variance estimates. STEP #6 K. Horne, 1986
     * 2. Mask cosmic ray hits.      STEP #7 K. Horne, 1986
     * 3. Extract optimal spectrum.  STEP #8 K. Horne, 1986     
     */       
    ostream *current_pout = NULL;
    
    if(iterations<=1) {
        current_pout = pout;
    }
    
    if(verbose)
        cerr << "operaSpectralOrder::extractOptimalSpectrum: iter=0 measuring optimal spectrum..." << endl;    
    measureOptimalSpectrum(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction, NMORESIGMASTOSTARTOPTIMALEXTRACTION*sigmaClip, current_pout);
    /*
     * Below it refines profile measurements and iterates optimal extraction. STEP #9 K. Horne, 1986    
     */      
    for(unsigned iter = 1; iter < iterations; iter++){
        if(onTargetProfile) {
            if(verbose)
                cerr << "operaSpectralOrder::extractOptimalSpectrum: iter="<<iter<<" refining spatial profile on object image..." << endl;            
            refineBeamSpatialProfiles(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,NumberofElementsToBin, sigmaClip, usePolynomialFit);
        } else {
            if(verbose)
                cerr << "operaSpectralOrder::extractOptimalSpectrum: iter="<<iter<<" refining spatial profile on flat-field image..." << endl;                        
            refineBeamSpatialProfiles(flatImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,NumberofElementsToBin, sigmaClip, usePolynomialFit);
        }
        if(iter == iterations-1) {
            current_pout = pout;
        }       
        if(verbose)
            cerr << "operaSpectralOrder::extractOptimalSpectrum: iter="<<iter<<" measuring optimal spectrum..." << endl;            
        
        measureOptimalSpectrum(objectImage,nflatImage,biasImage,badpix,gainBiasNoise,effectiveApertureFraction,sigmaClip,current_pout);
    }
	
	//
	// get rid of stuff we don't need anymore
	//
    for (unsigned beam=0; beam < numberOfBeams; beam++) {
		if (BeamProfiles[beam]) {
			delete BeamProfiles[beam];
			BeamProfiles[beam] = NULL;
		}
	}
    if(verbose)
        cerr << "operaSpectralOrder::extractOptimalSpectrum: optimal extraction run successfully! exiting..." << endl;            
    
	setSpectrumType(OptimalBeamSpectrum);
    SpectralElements->setSpectrumType(OptimalBeamSpectrum);
    SpectralElements->setHasOptimalSpectrum(true);
    sethasSpectralElements(true);
}

void operaSpectralOrder::measureBeamSpatialProfiles(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, bool usePolynomialFit) {
    
    // Set master spectral elements either if no elements has been created
    if(!hasSpectralElements) {
        // Set master spectral elements.
        setSpectralElementsByStitchingApertures(effectiveApertureFraction);
        // Set spectral elements for each beam aperture and for the background apertures
        setApertureElements(SpectralElements->getSpectrumType());
    }
    
    unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
#ifdef PRINT_OUTOFBOUNDS
	if (NumberofElements == 0)
		cerr << "operaSpectralOrder::measureSpatialProfileWithinAperture Warning: NumberofElements (" << NumberofElements << ") == 0" << endl;
#endif
	
    // pixel set for each beam aperture
    PixelSet *aperturePixels[MAXNUMBEROFBEAMS];
    
    // The IP dimension is given by the number of beams times the number of subpixels in pixelset   
    unsigned NXPoints = 0;
    unsigned NYPoints = numberOfBeams;
    
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        aperturePixels[beam] = ExtractionApertures[beam]->getSubpixels();
        if(aperturePixels[beam]->getNPixels() > NXPoints) {
            NXPoints = aperturePixels[beam]->getNPixels();
        }
        if(BeamProfiles[beam]) {
            delete BeamProfiles[beam];
        }
        BeamProfiles[beam] = new operaInstrumentProfile(aperturePixels[beam]->getNPixels(),1,1,1, NumberofElements);
    }
    
#ifdef PRINT_OUTOFBOUNDS
    if (NYPoints == 0)
        cerr << "operaSpectralOrder::measureSpatialProfileWithinAperture Warning: sum of aperturePixels[beam]->getNPixels() (" << NYPoints << ") == 0" << endl;
#endif
    
    CMatrix ipmatrix = newCMatrix(NXPoints,NYPoints);
    if (!ipmatrix) {
        throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
    }
    
    if(InstrumentProfile)
        delete InstrumentProfile;
    InstrumentProfile = new operaInstrumentProfile(NXPoints,1,NYPoints,1, NumberofElements);
    
    unsigned ipIndex = 0;
    for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
        float distd = (float)SpectralElements->getdistd(indexElem);
        // get elem center coordinates
        double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
        double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
        
        double BackgroundFlux = 0;
        for(unsigned background=0;background<LEFTANDRIGHT;background++) {
            BackgroundFlux += BackgroundElements[background]->getFlux(indexElem)/(double)LEFTANDRIGHT;
        }
        
        float SumOfUsefulFluxWithinAperture = 0;
        float FractionOfFluxLost = 0;
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            
            double subpixelArea = (double)(aperturePixels[beam]->getSubpixelArea());
            
            for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                
                // select image col and row of subpixel
                unsigned col = (unsigned)floor(elemXcenter + aperturePixels[beam]->getXcenter(pix));
                unsigned row = (unsigned)floor(elemYcenter + aperturePixels[beam]->getYcenter(pix));
                
                if(col >= 0 && row >= 0 && col < inputImage.getnaxis1() && row < inputImage.getnaxis2()
                   && inputImage[row][col] < SATURATIONLIMIT && badpix[row][col] == 1 && nflatImage[row][col]) {
                    
					double gain = gainBiasNoise.getGain(col,row);
					double subPixelRawFlux = (double)(gain*(inputImage[row][col] - biasImage[row][col])/nflatImage[row][col]);
					
					if(subPixelRawFlux > BackgroundFlux && subPixelRawFlux > 0) {
						ipmatrix[beam][pix] = (float)((subPixelRawFlux-BackgroundFlux)*subpixelArea);
						SumOfUsefulFluxWithinAperture += ipmatrix[beam][pix];
					} else {
						ipmatrix[beam][pix] = 0.0;
					}
                    
                } else {
                    ipmatrix[beam][pix] = NAN;
                    FractionOfFluxLost += 1.0/(float)(aperturePixels[beam]->getNPixels());
                }
            }
        }
        
        float NormalizationFactor = 0;
        if(FractionOfFluxLost<1)
            NormalizationFactor = SumOfUsefulFluxWithinAperture/(1.0 - FractionOfFluxLost);
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                if(NormalizationFactor && !isnan(ipmatrix[beam][pix])) {
                    ipmatrix[beam][pix] /= NormalizationFactor;
                } else {
                    ipmatrix[beam][pix] = 1.0/(float)(NXPoints*NYPoints);
                }
                
                BeamProfiles[beam]->setdataCubeValues(ipmatrix[beam][pix],pix,0,ipIndex);
            }
            BeamProfiles[beam]->setdistd(distd,ipIndex);
            BeamProfiles[beam]->normalizeCubeData(ipIndex);
        }
        InstrumentProfile->setdistd(distd,ipIndex);
        InstrumentProfile->setdataCubeValues(ipmatrix,ipIndex);
        InstrumentProfile->normalizeCubeData(ipIndex);
        ipIndex++;
    }
    
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        BeamProfiles[beam]->setnDataPoints(ipIndex);
    }
    InstrumentProfile->setnDataPoints(ipIndex);
    
    if(usePolynomialFit) {
        InstrumentProfile->FitPolyMatrixtoIPDataVector(3,false);
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            BeamProfiles[beam]->FitPolyMatrixtoIPDataVector(3,false);
        }
    } else {
        InstrumentProfile->FitMediantoIPDataVector();
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            BeamProfiles[beam]->FitMediantoIPDataVector();
        }
    }
	
    InstrumentProfile->setdataCubeFromPolyModel();
    InstrumentProfile->normalizeCubeData();
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        BeamProfiles[beam]->setdataCubeFromPolyModel();
        BeamProfiles[beam]->normalizeCubeData();
    }
    
    deleteCMatrix(ipmatrix);
    ipmatrix = NULL;
}

void operaSpectralOrder::measureOptimalSpectrum(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix, GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned sigmaClip, ostream *pout) {
    
    PixelSet *aperturePixels[MAXNUMBEROFBEAMS];
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        aperturePixels[beam] = ExtractionApertures[beam]->getSubpixels();
    }
	
    unsigned NumberofElements = SpectralElements->getnSpectralElements();
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
    CMatrix ipmatrix = newCMatrix(NXPoints,NYPoints);
    if (!ipmatrix) {
        throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
    }
    CMatrix BeamIPmatrix = newCMatrix(NXPoints,1);
    
    for(unsigned indexElem=0;indexElem < NumberofElements; indexElem++) {
        InstrumentProfile->getdataCubeValues(indexElem,ipmatrix);
        
        // get elem center coordinates
        double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
        double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
        
        double BackgroundFlux=0;
        double BackgroundFluxVar=0;
        for(unsigned background=0;background<LEFTANDRIGHT;background++) {
            BackgroundFlux += BackgroundElements[background]->getFlux(indexElem)/(double)LEFTANDRIGHT;
            BackgroundFluxVar += BackgroundElements[background]->getFluxVariance(indexElem)/(double)LEFTANDRIGHT;
        }
        
        double OldFlux = SpectralElements->getFlux(indexElem);
        //double OldVariance = SpectralElements->getFluxVariance(indexElem);
		
        double optimalFluxDenominatorAllBeams = 0;
        double optimalFluxNumeratorAllBeams = 0;
        double SumOfUsefulFluxWithinApertureAllBeams = 0;
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            
            double OldBeamFlux = BeamElements[beam]->getFlux(indexElem);
            //double OldBeamVariance = BeamElements[beam]->getFluxVariance(indexElem);
            
            BeamProfiles[beam]->getdataCubeValues(indexElem,BeamIPmatrix);
            
            double optimalBeamFluxDenominator = 0;
            double optimalBeamFluxNumerator = 0;
            double SumOfUsefulBeamFluxWithinAperture = 0;
            
            double subpixelArea = (double)(aperturePixels[beam]->getSubpixelArea());
			
            // loop over all subpixels within the beam aperture to calculate the sum of subpixel raw values
            for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                // select image col and row of subpixel
                unsigned col = (unsigned)floor(elemXcenter + aperturePixels[beam]->getXcenter(pix));
                unsigned row = (unsigned)floor(elemYcenter + aperturePixels[beam]->getYcenter(pix));
                
                if(col >= 0 && row >= 0 && col < inputImage.getnaxis1() && row < inputImage.getnaxis2()
                   && inputImage[row][col] < SATURATIONLIMIT && badpix[row][col] == 1 && nflatImage[row][col]) {
                    
                    double gain = gainBiasNoise.getGain(col,row);
                    double noise = gainBiasNoise.getNoise(col,row);

                    double subPixelRawFlux = (double)(gain*(inputImage[row][col] - biasImage[row][col])/nflatImage[row][col]);
                    
                    double BeamResidual = fabs((subPixelRawFlux-BackgroundFlux) - (double)(BeamIPmatrix[0][pix])*OldBeamFlux/subpixelArea)/sqrt(fabs(OldBeamFlux)*sqrt(subpixelArea) + noise);
                    double RevisedBeamVariance = 2.0*noise*noise*sqrt(subpixelArea) + fabs((double)BeamIPmatrix[0][pix]*OldBeamFlux + BackgroundFlux*sqrt(subpixelArea));
					
                    if((BeamResidual*BeamResidual/RevisedBeamVariance) < (double)sigmaClip) {
                        optimalBeamFluxDenominator += (double)(BeamIPmatrix[0][pix]*BeamIPmatrix[0][pix])/RevisedBeamVariance;
                        SumOfUsefulBeamFluxWithinAperture += (double)BeamIPmatrix[0][pix];
                        optimalBeamFluxNumerator += (double)(BeamIPmatrix[0][pix])*(subPixelRawFlux-BackgroundFlux)*subpixelArea/RevisedBeamVariance;
                    }
                    
                    double Residual = fabs((subPixelRawFlux-BackgroundFlux) - (double)(ipmatrix[beam][pix])*OldFlux/subpixelArea)/sqrt(fabs(OldFlux)*sqrt(subpixelArea) + noise);
                    double RevisedVariance = 2.0*noise*noise*sqrt(subpixelArea) + fabs((double)ipmatrix[beam][pix]*OldFlux + BackgroundFlux*sqrt(subpixelArea));
					               
                    if((Residual*Residual/RevisedVariance) < ((double)sigmaClip)) {
                        optimalFluxDenominatorAllBeams += (double)(ipmatrix[beam][pix]*ipmatrix[beam][pix])/RevisedVariance;
                        SumOfUsefulFluxWithinApertureAllBeams += (double)ipmatrix[beam][pix];
                        optimalFluxNumeratorAllBeams += (double)(ipmatrix[beam][pix])*(subPixelRawFlux-BackgroundFlux)*subpixelArea/RevisedVariance;
                    }
#ifdef PRINT_DEBUG
                    // Print out all subpixels for a single spectral element in the order
                    if(indexElem == NumberofElements/2) {
                    //if(indexElem > 2475 && indexElem < 2485) {
                        cout<<indexElem << ' ' <<   // 1
                        beam            << ' ' <<   // 2
                        pix             << ' ' <<   // 3
                        col             << ' ' <<   // 4
                        row             << ' ' <<   // 5
                    ipmatrix[beam][pix] << ' ' <<   // 6
                        OldFlux         << ' ' <<   // 7
                        subPixelRawFlux << ' ' <<   // 8
                        Residual        << ' ' <<   // 9
                        RevisedVariance << ' ' <<   // 10
   (Residual*Residual/RevisedVariance)  << ' ' <<   // 11
                      (double)sigmaClip << endl;    // 12
                    }
#endif
                }
            }
            
            if(optimalBeamFluxDenominator) {
                double optimalBeamFlux = optimalBeamFluxNumerator/optimalBeamFluxDenominator;
                double optimalBeamFluxVar = SumOfUsefulBeamFluxWithinAperture/optimalBeamFluxDenominator;
                
                BeamElements[beam]->setFlux(optimalBeamFlux,indexElem);
                BeamElements[beam]->setFluxVariance(optimalBeamFluxVar,indexElem);
            }
        }
        
        if(optimalFluxDenominatorAllBeams) {
            double optimalFlux = optimalFluxNumeratorAllBeams/optimalFluxDenominatorAllBeams;
            double optimalFluxVar = SumOfUsefulFluxWithinApertureAllBeams/optimalFluxDenominatorAllBeams;
            
            SpectralElements->setFlux(optimalFlux,indexElem);
            SpectralElements->setFluxVariance(optimalFluxVar,indexElem);
        }
    }
	
    printBeamSpectrum(pout);
    
    SpectralElements->setHasOptimalSpectrum(true);
    setSpectrumType(OptimalBeamSpectrum);
    SpectralElements->setSpectrumType(OptimalBeamSpectrum);
    
    deleteCMatrix(ipmatrix);
    ipmatrix = NULL;
    deleteCMatrix(BeamIPmatrix);
    BeamIPmatrix = NULL;
}

void operaSpectralOrder::refineBeamSpatialProfiles(operaFITSImage &inputImage, operaFITSImage &nflatImage, operaFITSImage &biasImage, operaFITSImage &badpix,GainBiasNoise &gainBiasNoise, double effectiveApertureFraction, unsigned NumberofElementsToBin, unsigned sigmaClip, bool usePolynomialFit) {
	
    unsigned NumberofElements = SpectralElements->getnSpectralElements();
	
#ifdef PRINT_OUTOFBOUNDS
	if (NumberofElements == 0)
		cerr << "operaSpectralOrder::refineSpatialProfileWithinApertures Warning: NumberofElements (" << NumberofElements << ") == 0" << endl;
#endif
    
 	unsigned NumberOfElementSamples = (unsigned)ceil((float)NumberofElements/(float)NumberofElementsToBin);
    
#ifdef PRINT_OUTOFBOUNDS
	if (NumberOfElementSamples == 0)
		cerr << "operaSpectralOrder::refineSpatialProfileWithinApertures Warning: NumberOfElementSamples (" << NumberOfElementSamples << ") == 0" << endl;
#endif
    
    unsigned NXPoints = InstrumentProfile->getNXPoints();
    unsigned NYPoints = InstrumentProfile->getNYPoints();
    
#ifdef PRINT_OUTOFBOUNDS
    if (NXPoints == 0)
        cerr << "operaSpectralOrder::refineSpatialProfileWithinApertures Warning: NXPoints (" << NXPoints << ") == 0" << endl;
    if (NYPoints == 0)
        cerr << "operaSpectralOrder::refineSpatialProfileWithinApertures Warning: NYPoints (" << NYPoints << ") == 0" << endl;
#endif
    CMatrix ipmatrix = newCMatrix(NXPoints,NYPoints);
    if (!ipmatrix) {
        throw operaException("operaSpectralOrder: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
    }
    CMatrix BeamIPmatrix = newCMatrix(NXPoints,1);
	
    operaInstrumentProfile sampleIP(NXPoints,1,NYPoints,1,NumberofElementsToBin);
    operaInstrumentProfile tempIP(NXPoints,1,NYPoints,1,NumberOfElementSamples);
    
    PixelSet *aperturePixels[MAXNUMBEROFBEAMS];
    operaInstrumentProfile *sampleBeamIP[MAXNUMBEROFBEAMS];
    operaInstrumentProfile *tempBeamIP[MAXNUMBEROFBEAMS];
    
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        aperturePixels[beam] = ExtractionApertures[beam]->getSubpixels();
        sampleBeamIP[beam] = new operaInstrumentProfile(NXPoints,1,1,1,NumberofElementsToBin);
        tempBeamIP[beam] = new operaInstrumentProfile(NXPoints,1,1,1,NumberOfElementSamples);
    }
	
    unsigned ipIndex = 0;
    
    for(unsigned k=0;k<NumberOfElementSamples;k++){
        unsigned firstElement = NumberofElementsToBin*(k);
        unsigned lastElement =  NumberofElementsToBin*(k+1);
        if (lastElement > SpectralElements->getnSpectralElements()) {
            lastElement = SpectralElements->getnSpectralElements();
        }
        
        unsigned zeroBasedIndex=0;
        
        for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
            
            float distd = (float)SpectralElements->getdistd(indexElem);
            InstrumentProfile->getdataCubeValues(indexElem,ipmatrix);
            
            // get elem center coordinates
            double elemXcenter = SpectralElements->getphotoCenterX(indexElem);
            double elemYcenter = SpectralElements->getphotoCenterY(indexElem);
            
            double BackgroundFlux=0;
            double BackgroundFluxVar=0;
            
            for(unsigned background=0;background<LEFTANDRIGHT;background++) {
                BackgroundFlux += BackgroundElements[background]->getFlux(indexElem)/(double)LEFTANDRIGHT;
                BackgroundFluxVar += BackgroundElements[background]->getFluxVariance(indexElem)/(double)LEFTANDRIGHT;
            }
			
            double OldFlux = SpectralElements->getFlux(indexElem);
            //double OldVariance = SpectralElements->getFluxVariance(indexElem);
			
            float FractionOfFluxLost = 0;
            float SumOfUsefulFluxWithinAperture = 0;
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                
                float FractionOfBeamFluxLost = 0;
                float SumOfUsefulBeamFluxWithinAperture = 0;
				
                double OldBeamFlux = BeamElements[beam]->getFlux(indexElem);
                //double OldBeamVariance = BeamElements[beam]->getFluxVariance(indexElem);
                
                BeamProfiles[beam]->getdataCubeValues(indexElem,BeamIPmatrix);
                
                double subpixelArea = (double)(aperturePixels[beam]->getSubpixelArea());
				
                // loop over all subpixels within the beam aperture to calculate the sum of subpixel raw values
                for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                    // select image col and row of subpixel
                    unsigned col = (unsigned)floor(elemXcenter + aperturePixels[beam]->getXcenter(pix));
                    unsigned row = (unsigned)floor(elemYcenter + aperturePixels[beam]->getYcenter(pix));
                    
                    if(col >= 0 && row >= 0 && col < inputImage.getnaxis1() && row < inputImage.getnaxis2()
                       && inputImage[row][col] < SATURATIONLIMIT && badpix[row][col] == 1 && nflatImage[row][col]) {
                        
                        double gain = gainBiasNoise.getGain(col,row);
                        double noise = gainBiasNoise.getNoise(col,row);
                        
                        double subPixelRawFlux = (double)(gain*(inputImage[row][col] - biasImage[row][col])/nflatImage[row][col]);
                        
                        if(!isnan(OldBeamFlux) && OldBeamFlux)  {
                            
                            double BeamResidual = fabs((subPixelRawFlux-BackgroundFlux) - (double)(BeamIPmatrix[0][pix])*OldBeamFlux/subpixelArea)/sqrt(fabs(OldBeamFlux)*sqrt(subpixelArea) + noise);
                            double RevisedBeamVariance = 2.0*noise*noise*sqrt(subpixelArea) + fabs((double)BeamIPmatrix[0][pix]*OldBeamFlux + BackgroundFlux*sqrt(subpixelArea));
                            
                            if((BeamResidual*BeamResidual/RevisedBeamVariance) < ((double)sigmaClip)) {
                                if(subPixelRawFlux > BackgroundFlux) {
									BeamIPmatrix[0][pix] = (float)((subPixelRawFlux - BackgroundFlux)*subpixelArea);
                                    SumOfUsefulBeamFluxWithinAperture += BeamIPmatrix[0][pix];
                                } else {
                                    BeamIPmatrix[0][pix] = 0;
                                }
                            } else {
                                FractionOfBeamFluxLost += BeamIPmatrix[0][pix];
                            }
                        } else {
                            FractionOfBeamFluxLost += BeamIPmatrix[0][pix];
                        }
                        
                        
                        if(!isnan(OldFlux) && OldFlux)  {
                            double Residual = fabs((subPixelRawFlux-BackgroundFlux) - (double)(ipmatrix[beam][pix])*OldFlux/subpixelArea)/sqrt(fabs(OldFlux)*sqrt(subpixelArea) + noise);
                            double RevisedVariance = 2.0*noise*noise*sqrt(subpixelArea) + fabs((double)ipmatrix[beam][pix]*OldFlux + BackgroundFlux*sqrt(subpixelArea));
                                                        
                            if((Residual*Residual/RevisedVariance) < (double)sigmaClip) {
                                if(subPixelRawFlux > BackgroundFlux) {
                                    ipmatrix[beam][pix] = (float)((subPixelRawFlux - BackgroundFlux)*subpixelArea);
                                    SumOfUsefulFluxWithinAperture += ipmatrix[beam][pix];
                                } else {
                                    ipmatrix[beam][pix] = 0;
                                }
                            } else {
                                FractionOfFluxLost += ipmatrix[beam][pix];
                            }
                        } else {
                            FractionOfFluxLost += ipmatrix[beam][pix];
                        }
                    } else {
                        FractionOfFluxLost += ipmatrix[beam][pix];
                        FractionOfBeamFluxLost += BeamIPmatrix[0][pix];
                    }
                }
				
                float BeamNormalizationFactor = 0;
                
                if(FractionOfBeamFluxLost<1)
                    BeamNormalizationFactor = SumOfUsefulBeamFluxWithinAperture/(1.0 - FractionOfBeamFluxLost);
                
                for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                    if(BeamNormalizationFactor && !isnan(BeamIPmatrix[0][pix])){
                        BeamIPmatrix[0][pix] /= BeamNormalizationFactor;
                    } else if (isnan(BeamIPmatrix[0][pix])) {
                        BeamIPmatrix[0][pix] = 1.0/(float)(NXPoints);
                    }
                }
                
                sampleBeamIP[beam]->setdistd(distd,zeroBasedIndex);
                sampleBeamIP[beam]->setdataCubeValues(BeamIPmatrix,zeroBasedIndex);
                sampleBeamIP[beam]->normalizeCubeData(zeroBasedIndex);
            }
			
            float NormalizationFactor = 0;
            
            if(FractionOfFluxLost<1)
                NormalizationFactor = SumOfUsefulFluxWithinAperture/(1.0 - FractionOfFluxLost);
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                for(unsigned pix=0; pix<aperturePixels[beam]->getNPixels(); pix++) {
                    if(NormalizationFactor && !isnan(ipmatrix[beam][pix])){
                        ipmatrix[beam][pix] /= NormalizationFactor;
                    } else if (isnan(ipmatrix[beam][pix])) {
                        ipmatrix[beam][pix] = 1.0/(float)(NXPoints*NYPoints);
                    }
                }
            }
            
            sampleIP.setdistd(distd,zeroBasedIndex);
            sampleIP.setdataCubeValues(ipmatrix,zeroBasedIndex);
            sampleIP.normalizeCubeData(zeroBasedIndex);
            
            zeroBasedIndex++;
        }
        
        float distdmidElem = (float)fabs(SpectralElements->getdistd(lastElement-1) + SpectralElements->getdistd(firstElement))/2;
        
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            sampleBeamIP[beam]->setnDataPoints(zeroBasedIndex);
            sampleBeamIP[beam]->FitMediantoIPDataVector();
			
            tempBeamIP[beam]->setdistd(distdmidElem,ipIndex);
            tempBeamIP[beam]->setdataCubeValues((sampleBeamIP[beam]->getipDataFromPolyModel(distdmidElem,BeamIPmatrix)),ipIndex);
            tempBeamIP[beam]->normalizeCubeData(ipIndex);
        }
        
        sampleIP.setnDataPoints(zeroBasedIndex);
        sampleIP.FitMediantoIPDataVector();
        
        tempIP.setdistd(distdmidElem,ipIndex);
        tempIP.setdataCubeValues((sampleIP.getipDataFromPolyModel(distdmidElem,ipmatrix)),ipIndex);
        tempIP.normalizeCubeData(ipIndex);
        ipIndex++;
    }
    
    tempIP.setnDataPoints(ipIndex);
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        tempBeamIP[beam]->setnDataPoints(ipIndex);
    }
    
    if(usePolynomialFit) {
        tempIP.FitPolyMatrixtoIPDataVector(3,false);
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            tempBeamIP[beam]->FitPolyMatrixtoIPDataVector(3,false);
        }
    } else {
        tempIP.FitMediantoIPDataVector();
        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
            tempBeamIP[beam]->FitMediantoIPDataVector();
        }
    }
    
    /*
     * BAD! memory issue, you may NOT copy a temporary adress into
     * an object which is allocated.
     * FIXED Oct 24 2012 DT
     */
    InstrumentProfile->setipPolyModel(tempIP.getipPolyModel());
    InstrumentProfile->setdataCubeFromPolyModel();
    InstrumentProfile->normalizeCubeData();
    
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        BeamProfiles[beam]->setipPolyModel(tempBeamIP[beam]->getipPolyModel());
        BeamProfiles[beam]->setdataCubeFromPolyModel();
        BeamProfiles[beam]->normalizeCubeData();
    }
    
    deleteCMatrix(ipmatrix);
    ipmatrix = NULL;
	
    deleteCMatrix(BeamIPmatrix);
    BeamIPmatrix = NULL;
    
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        if(sampleBeamIP[beam])
            delete sampleBeamIP[beam];
        if(tempBeamIP[beam])
			delete tempBeamIP[beam];
    }
}

void operaSpectralOrder::createSpectralEnergyDistribution(unsigned binsize) {
    if(!gethasSpectralElements()) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    if(binsize == 0) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    unsigned nDataPoints = (unsigned)ceil((float)(SpectralElements->getFluxVector()->getlength())/(float)binsize) + 2;
    if(nDataPoints == 0 || SpectralElements->getnSpectralElements() == 0) {
        throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    if(SpectralEnergyDistribution)
        delete SpectralEnergyDistribution;
    SpectralEnergyDistribution = new operaSpectralEnergyDistribution(nDataPoints,SpectralElements->getnSpectralElements());
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        unsigned BeamNDataPoints = (unsigned)ceil((float)(BeamElements[beam]->getFluxVector()->getlength())/(float)binsize) + 2;
        if(BeamSED[beam])
            delete BeamSED[beam]; 
		if(BeamNDataPoints == 0 || BeamElements[beam]->getnSpectralElements() == 0) {
			throw operaException("operaSpectralOrder::setSpectralLines: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
		}
        BeamSED[beam] = new operaSpectralEnergyDistribution(BeamNDataPoints,BeamElements[beam]->getnSpectralElements());
    } 
}

void operaSpectralOrder::deleteSpectralEnergyDistribution(void) {
    if(SpectralEnergyDistribution)
        delete SpectralEnergyDistribution;
	SpectralEnergyDistribution = NULL;
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        if(BeamSED[beam])
            delete BeamSED[beam];            
        BeamSED[beam] = NULL;
    } 
}

operaSpectralEnergyDistribution *operaSpectralOrder::getBeamSED(unsigned beam) {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamSED[beam];
}

const operaSpectralEnergyDistribution *operaSpectralOrder::getBeamSED(unsigned beam) const {
#ifdef RANGE_CHECK
	if (beam >= MAXNUMBEROFBEAMS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return BeamSED[beam];
}

void operaSpectralOrder::fitSEDUncalibratedFluxToSample(unsigned uniform_npoints, float *uniform_wl, float *uniform_flux, float **uniform_beamflux) {
	
	// Create arrays to store the flux for our model
	unsigned nElements = SpectralElements->getnSpectralElements();
	float *UncalibratedModelFlux = new float[nElements];
    float *UncalibratedModelBeamFlux[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<numberOfBeams;beam++) {
        UncalibratedModelBeamFlux[beam] = new float[nElements];
    }
    
    // Fit a spline to the uniform sample along the wavelength of our SpectralElements and store the result in UncalibratedModelFlux.
	float *elemWavelength = new float[nElements];
	for(unsigned i=0;i<nElements;i++) {
		elemWavelength[i] = (float)SpectralElements->getwavelength(i);
	}
	operaFitSpline(uniform_npoints,uniform_wl,uniform_flux,nElements,elemWavelength,UncalibratedModelFlux);
	for(unsigned beam=0;beam<numberOfBeams;beam++) {
		operaFitSpline(uniform_npoints,uniform_wl,uniform_beamflux[beam],nElements,elemWavelength,UncalibratedModelBeamFlux[beam]);
	}
	
	// Copy the existing SpectralElements into the UncalibratedFlux of the SEDs.
	operaSpectralElements *uncalibratedFluxElements;
	operaSpectralElements *uncalibratedBeamFluxElements[MAXNUMBEROFBEAMS];
	SpectralEnergyDistribution->setUncalibratedFluxElements(SpectralElements);
	uncalibratedFluxElements = SpectralEnergyDistribution->getUncalibratedFluxElements();
	for(unsigned beam=0;beam<numberOfBeams;beam++) {
		BeamSED[beam]->setUncalibratedFluxElements(BeamElements[beam]);
		uncalibratedBeamFluxElements[beam] = BeamSED[beam]->getUncalibratedFluxElements();
	}
	
	// Set the flux of UncalibratedFlux of the SEDs to the values in UncalibratedModelFlux.
	for(unsigned i=0;i<nElements;i++) {
		uncalibratedFluxElements->setFlux((double)UncalibratedModelFlux[i],i);
		for(unsigned beam=0;beam<numberOfBeams;beam++) {
			uncalibratedBeamFluxElements[beam]->setFlux((double)UncalibratedModelBeamFlux[beam][i],i);
		}
	}
	SpectralEnergyDistribution->setHasUncalibratedFlux(true);
	for(unsigned beam=0;beam<numberOfBeams;beam++) {
		BeamSED[beam]->setHasUncalibratedFlux(true);
	}
	
	// Cleanup
	delete[] elemWavelength;
	delete[] UncalibratedModelFlux;
    for(unsigned beam=0;beam<numberOfBeams;beam++) {
        delete[] UncalibratedModelBeamFlux[beam];
    }
}

void operaSpectralOrder::fitSEDFcalAndThroughputToFlatResp(unsigned npoints, float *frwavelength, float *flatresp) {
	
	// Create arrays to store the flat response flux for our model
    unsigned nElements = SpectralElements->getnSpectralElements();
    float *flatResponseModel = new float[nElements];
    
    // Fit a spline to the flat response along the wavelength of our SpectralElements and store the result in flatResponseModel.
	float *elemWavelength = new float[nElements];
	for (unsigned elemIndex=0; elemIndex<nElements; elemIndex++) {
		elemWavelength[elemIndex] = SpectralElements->getwavelength(elemIndex);
	}
	operaFitSpline(npoints,frwavelength,flatresp,nElements,elemWavelength,flatResponseModel);
	
	// Get the FluxCalibration and InstrumentThroughput of the SEDs
	operaSpectralElements *FluxCalibration = SpectralEnergyDistribution->getFluxCalibrationElements();
	operaSpectralElements *InstrumentThroughput = SpectralEnergyDistribution->getThroughputElements();
	operaSpectralElements *beamFluxcalibration[MAXNUMBEROFBEAMS];
	operaSpectralElements *beamThroughput[MAXNUMBEROFBEAMS];
	for (unsigned beam=0; beam < numberOfBeams; beam++) {
		beamFluxcalibration[beam] = BeamSED[beam]->getFluxCalibrationElements();
		beamThroughput[beam] = BeamSED[beam]->getThroughputElements();
	}
	
	// Set the wl, flux, flux variance, and photoCenter to the values in 1/flatResponseModel.
	for (unsigned elemIndex=0; elemIndex<nElements; elemIndex++) {
		double wl = SpectralElements->getwavelength(elemIndex);
		double fluxcal = 1.0/(double)flatResponseModel[elemIndex];
		double throughput = 1.0/(double)flatResponseModel[elemIndex];
		
		FluxCalibration->setwavelength(wl, elemIndex);
		FluxCalibration->setFlux(fluxcal, elemIndex);
		FluxCalibration->setFluxVariance(fluxcal, elemIndex);
		FluxCalibration->setphotoCenter(0.0, 0.0, elemIndex);
		
		InstrumentThroughput->setwavelength(wl, elemIndex);
		InstrumentThroughput->setFlux(throughput, elemIndex);
		InstrumentThroughput->setFluxVariance(throughput, elemIndex);
		InstrumentThroughput->setphotoCenter(0.0, 0.0, elemIndex);
		
		for (unsigned beam=0; beam < numberOfBeams; beam++) {
			beamFluxcalibration[beam]->setwavelength(wl, elemIndex);
			beamFluxcalibration[beam]->setFlux(fluxcal, elemIndex);
			beamFluxcalibration[beam]->setFluxVariance(fluxcal, elemIndex);
			beamFluxcalibration[beam]->setphotoCenter(0.0, 0.0, elemIndex);
			
			beamThroughput[beam]->setwavelength(wl, elemIndex);
			beamThroughput[beam]->setFlux(throughput, elemIndex);
			beamThroughput[beam]->setFluxVariance(throughput, elemIndex);
			beamThroughput[beam]->setphotoCenter(0.0, 0.0, elemIndex);
		}
	}
	sethasSpectralEnergyDistribution(true);
	FluxCalibration->setHasWavelength(true);
	InstrumentThroughput->setHasWavelength(true);
	SpectralEnergyDistribution->setHasFluxCalibration(true);
	SpectralEnergyDistribution->setHasInstrumentThroughput(true);
	for (unsigned beam=0; beam < numberOfBeams; beam++) {
		beamFluxcalibration[beam]->setHasWavelength(true);
		beamThroughput[beam]->setHasWavelength(true);
		BeamSED[beam]->setHasFluxCalibration(true);
		BeamSED[beam]->setHasInstrumentThroughput(true);
	}
	
	// Cleanup
	delete[] flatResponseModel;
	delete[] elemWavelength;
}

void operaSpectralOrder::applyNormalization(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, ostream *poutspec, ostream *poutcontinuum, bool overwriteUncalFlux, unsigned numberOfprintouts) {
    //unsigned NumberofElements = SpectralElements->getnSpectralElements();	
    
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    createSpectralEnergyDistribution(binsize);      
	
    operaFluxVector *normalizedBeamFlux[MAXNUMBEROFBEAMS];
    operaFluxVector *beamContinuum[MAXNUMBEROFBEAMS];
	operaFluxVector *spectralFlux = SpectralElements->getFluxVector();
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        normalizedBeamFlux[beam]  = new operaFluxVector(BeamElements[beam]->getFluxVector()->getlength());
        beamContinuum[beam] = new operaFluxVector(BeamElements[beam]->getFluxVector()->getlength());
		operaFluxVector *beamFlux = BeamElements[beam]->getFluxVector();
        normalizeSpectrum(*beamFlux, *(normalizedBeamFlux[beam]), *(beamContinuum[beam]), *(BeamSED[beam]), binsize, orderOfPolynomial, usePolynomial);
    }
	
    unsigned numberOfElements = spectralFlux->getlength();
	
	if (numberOfElements == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    operaFluxVector normalizedFlux(numberOfElements);
    operaFluxVector outputContinuum(numberOfElements);
    
    normalizeSpectrum(*spectralFlux, normalizedFlux, outputContinuum, *SpectralEnergyDistribution, binsize, orderOfPolynomial, usePolynomial);
	
    if (gethasWavelength() && SpectralElements->getHasDistance()) {
        SpectralElements->setwavelengthsFromCalibration(getWavelength());
    }
    
    if (poutspec != NULL) {     
        for(unsigned slitview=0; slitview<numberOfprintouts; slitview++) {
            for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
                *poutspec << slitview << " "
                << orderNumber << " "
                << indexElem << " "
                << SpectralElements->getphotoCenterX(indexElem) << " "
                << SpectralElements->getphotoCenterY(indexElem) << " "
                << SpectralElements->getdistd(indexElem) << " "
                << SpectralElements->getwavelength(indexElem) << " "
                << SpectralElements->getFlux(indexElem) << " "
                << outputContinuum.getflux(indexElem) << " "
                << normalizedFlux.getflux(indexElem) << " "
                << SpectralElements->getFluxVariance(indexElem) << " ";
                
                for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                    *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
                    << normalizedBeamFlux[beam]->getflux(indexElem) << " "
                    << beamContinuum[beam]->getflux(indexElem) << " "
                    << BeamElements[beam]->getFluxVariance(indexElem) << " ";
                }
                *poutspec << endl;
            }  
            *poutspec << endl;
        }
    }
	
    if (poutcontinuum != NULL) {
        for(unsigned slitview=0; slitview<numberOfprintouts; slitview++) {
            for(unsigned i=0;i < SpectralEnergyDistribution->getnDataPoints(); i++) {
                *poutcontinuum << slitview << " "
                << orderNumber << " "
                << i << " "
                << SpectralEnergyDistribution->getdistanceData(i) << " "
                << SpectralEnergyDistribution->getwavelengthData(i) << " "
                << SpectralEnergyDistribution->getfluxData(i) << " ";
                
                for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                    *poutcontinuum << beam << " "
                    << BeamSED[beam]->getfluxData(i) << " ";
                }
                *poutcontinuum << endl;     
            }
            *poutcontinuum << endl;
        }
    }
    
    if(overwriteUncalFlux) {
        for(unsigned beam=0; beam < numberOfBeams; beam++) {
            BeamElements[beam]->setFluxVector(normalizedBeamFlux[beam]);
        }    
        SpectralElements->setFluxVector(&normalizedFlux);    
    }
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        if (normalizedBeamFlux[beam])
			delete normalizedBeamFlux[beam];
        normalizedBeamFlux[beam] = NULL;
        if (beamContinuum[beam])
			delete beamContinuum[beam];
        beamContinuum[beam] = NULL;
    }      
    deleteSpectralEnergyDistribution();
}

void operaSpectralOrder::applyNormalizationForEmissionSpectrum(unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial, ostream *poutspec, ostream *poutcontinuum, bool overwriteUncalFlux, unsigned numberOfprintouts) {
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    createSpectralEnergyDistribution(binsize);
    unsigned numberOfElements = SpectralElements->getFluxVector()->getlength();
	
    operaFluxVector UnitVector(numberOfElements);
    UnitVector = 1.0;
    
    operaFluxVector *normalizedBeamFlux[MAXNUMBEROFBEAMS];
    operaFluxVector *beamContinuum[MAXNUMBEROFBEAMS];
	
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        normalizedBeamFlux[beam]  = new operaFluxVector(BeamElements[beam]->getFluxVector()->getlength());
        beamContinuum[beam] = new operaFluxVector(BeamElements[beam]->getFluxVector()->getlength());
        
        double maxbeamflux = -BIG;
        for(unsigned indexElem=0;indexElem < BeamElements[beam]->getFluxVector()->getlength(); indexElem++) {
            if(BeamElements[beam]->getFluxVector()->getflux(indexElem) > maxbeamflux) {
                maxbeamflux = BeamElements[beam]->getFluxVector()->getflux(indexElem);
            }
        }
        
		operaFluxVector *beamFlux = BeamElements[beam]->getFluxVector();
		operaFluxVector *normalizedFlux = normalizedBeamFlux[beam];
		operaFluxVector *continuumFlux = beamContinuum[beam];
		
        *beamFlux = UnitVector - (*beamFlux * (1.0/maxbeamflux));
        
        normalizeSpectrum(*beamFlux, *normalizedFlux, *continuumFlux, *(BeamSED[beam]), binsize, orderOfPolynomial, usePolynomial);
        
        *continuumFlux = UnitVector - *continuumFlux;
        *continuumFlux = UnitVector - *continuumFlux;
        *beamFlux = UnitVector - *beamFlux;
    }
	
    operaFluxVector normalizedFlux(numberOfElements);
    operaFluxVector outputContinuum(numberOfElements);
    
    double maxflux = -BIG;
    for(unsigned indexElem=0;indexElem < numberOfElements; indexElem++) {
        if(SpectralElements->getFluxVector()->getflux(indexElem) > maxflux) {
            maxflux = SpectralElements->getFluxVector()->getflux(indexElem);
        }
    }
    
	operaFluxVector *spectralFlux = SpectralElements->getFluxVector();
    *spectralFlux =  UnitVector - (*spectralFlux * (1.0/maxflux));
    
    normalizeSpectrum(*spectralFlux, normalizedFlux, outputContinuum, *SpectralEnergyDistribution, binsize, orderOfPolynomial, usePolynomial);
	
    normalizedFlux = UnitVector - normalizedFlux;
    outputContinuum = UnitVector - outputContinuum;
    *spectralFlux = UnitVector - *spectralFlux;
    
    if (gethasWavelength() && SpectralElements->getHasDistance()) {
        SpectralElements->setwavelengthsFromCalibration(getWavelength());
    }
    
    if (poutspec != NULL) {
        for(unsigned slitview=0; slitview<numberOfprintouts; slitview++) {
            for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
                *poutspec << slitview << " "
                << orderNumber << " "
                << indexElem << " "
                << SpectralElements->getphotoCenterX(indexElem) << " "
                << SpectralElements->getphotoCenterY(indexElem) << " "
                << SpectralElements->getdistd(indexElem) << " "
                << SpectralElements->getwavelength(indexElem) << " "
                << SpectralElements->getFlux(indexElem) << " "
                << outputContinuum.getflux(indexElem) << " "
                << normalizedFlux.getflux(indexElem) << " "
                << SpectralElements->getFluxVariance(indexElem) << " ";
                
                for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                    *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
                    << normalizedBeamFlux[beam]->getflux(indexElem) << " "
                    << beamContinuum[beam]->getflux(indexElem) << " "
                    << BeamElements[beam]->getFluxVariance(indexElem) << " ";
                }
                *poutspec << endl;
            }
            *poutspec << endl;
        }
    }
	
    if (poutcontinuum != NULL) {
        for(unsigned slitview=0; slitview<numberOfprintouts; slitview++) {
            for(unsigned i=0;i < SpectralEnergyDistribution->getnDataPoints(); i++) {
                *poutcontinuum << slitview << " "
                << orderNumber << " "
                << i << " "
                << SpectralEnergyDistribution->getdistanceData(i) << " "
                << SpectralEnergyDistribution->getwavelengthData(i) << " "
                << SpectralEnergyDistribution->getfluxData(i) << " ";
                
                for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                    *poutcontinuum << beam << " "
                    << BeamSED[beam]->getfluxData(i) << " ";
                }
                *poutcontinuum << endl;
            }
            *poutcontinuum << endl;
        }
    }
    
    if(overwriteUncalFlux) {
        for(unsigned beam=0; beam < numberOfBeams; beam++) {
            BeamElements[beam]->setFluxVector(&(*normalizedBeamFlux[beam]));
        }
        SpectralElements->setFluxVector(&normalizedFlux);
    }
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        if (normalizedBeamFlux[beam])
			delete normalizedBeamFlux[beam];
        normalizedBeamFlux[beam] = NULL;
        if (beamContinuum[beam])
			delete beamContinuum[beam];
        beamContinuum[beam] = NULL;
    }
    deleteSpectralEnergyDistribution();
}


void operaSpectralOrder::applyNormalizationFromExistingContinuum(ostream *poutspec, ostream *poutcontinuum, bool overwriteUncalFlux, bool normalizeBeams, unsigned numberOfprintouts) {
    
    operaFluxVector *normalizedBeamFlux[MAXNUMBEROFBEAMS];
    operaFluxVector *beamContinuum[MAXNUMBEROFBEAMS];
	operaFluxVector *spectralFlux = SpectralElements->getFluxVector();
    
    if(normalizeBeams) {
        for(unsigned beam=0; beam < numberOfBeams; beam++) {
            normalizedBeamFlux[beam] = new operaFluxVector(BeamElements[beam]->getFluxVector()->getlength());
            
            beamContinuum[beam] = BeamSED[beam]->getUncalibratedFluxElements()->getFluxVector();
            operaFluxVector *beamFlux = BeamElements[beam]->getFluxVector();
            
            (*normalizedBeamFlux[beam]) = (*beamFlux)/(*beamContinuum[beam]);
        }
	}
    unsigned numberOfElements = spectralFlux->getlength();
	
	if (numberOfElements == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
	}
    operaFluxVector *outputContinuum = SpectralEnergyDistribution->getUncalibratedFluxElements()->getFluxVector();
    operaFluxVector normalizedFlux = (*spectralFlux)/(*outputContinuum);
    
    if (gethasWavelength() && SpectralElements->getHasDistance()) {
        SpectralElements->setwavelengthsFromCalibration(getWavelength());
    }
    
    if (poutspec != NULL) {
        for(unsigned slitview=0; slitview<numberOfprintouts; slitview++) {
            for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
                *poutspec << slitview << " "
                << orderNumber << " "
                << indexElem << " "
                << SpectralElements->getphotoCenterX(indexElem) << " "
                << SpectralElements->getphotoCenterY(indexElem) << " "
                << SpectralElements->getdistd(indexElem) << " "
                << SpectralElements->getwavelength(indexElem) << " "
                << SpectralElements->getFlux(indexElem) << " "
                << outputContinuum->getflux(indexElem) << " "
                << normalizedFlux.getflux(indexElem) << " "
                << SpectralElements->getFluxVariance(indexElem) << " ";
                
                for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                    *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
                    << normalizedBeamFlux[beam]->getflux(indexElem) << " "
                    << beamContinuum[beam]->getflux(indexElem) << " "
                    << BeamElements[beam]->getFluxVariance(indexElem) << " ";
                }
                *poutspec << endl;
            }
            *poutspec << endl;
        }
    }
	
    if (poutcontinuum != NULL) {
        for(unsigned slitview=0; slitview<numberOfprintouts; slitview++) {
            for(unsigned i=0;i < SpectralEnergyDistribution->getnDataPoints(); i++) {
                *poutcontinuum << slitview << " "
                << orderNumber << " "
                << i << " "
                << SpectralEnergyDistribution->getdistanceData(i) << " "
                << SpectralEnergyDistribution->getwavelengthData(i) << " "
                << SpectralEnergyDistribution->getfluxData(i) << " ";
                
                for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                    *poutcontinuum << beam << " "
                    << BeamSED[beam]->getfluxData(i) << " ";
                }
                *poutcontinuum << endl;
            }
            *poutcontinuum << endl;
        }
    }
    
    if(overwriteUncalFlux) {
        if(normalizeBeams) {
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                BeamElements[beam]->setFluxVector(normalizedBeamFlux[beam]);
            }
        }
        SpectralElements->setFluxVector(&normalizedFlux);
    }
    
    if(normalizeBeams) {
        for(unsigned beam=0; beam < numberOfBeams; beam++) {
            if (normalizedBeamFlux[beam])
                delete normalizedBeamFlux[beam];
            normalizedBeamFlux[beam] = NULL;
        }
    }
}

void operaSpectralOrder::normalizeSpectrum(operaFluxVector &uncalibratedFlux, operaFluxVector &normalizedFlux, operaFluxVector &outputContinuum, operaSpectralEnergyDistribution &spectralEnergyDistribution, unsigned binsize, unsigned orderOfPolynomial, bool usePolynomial) {
    
	
    measureContinuum(uncalibratedFlux,outputContinuum,spectralEnergyDistribution,binsize,NSIGMACUT,orderOfPolynomial,usePolynomial);    
    
    normalizedFlux = uncalibratedFlux/outputContinuum;
    
#ifdef PRINT_DEBUG      
    for(unsigned i=0; i< uncalibratedFlux.getlength(); i++) {
        cout << i << " " << uncalibratedFlux.getflux(i) << " " << outputContinuum.getflux(i) << " " << normalizedFlux.getflux(i) << endl;
    }    
#endif    
}

void operaSpectralOrder::measureContinuum(operaFluxVector &uncalibratedFlux,operaFluxVector &outputContinuum, operaSpectralEnergyDistribution &spectralEnergyDistribution, unsigned binsize, unsigned nsigcut, unsigned orderOfPolynomial, bool usePolynomial) {
	
    unsigned NumberofPoints = uncalibratedFlux.getlength();
    
	unsigned NumberOfSamples = (unsigned)ceil((float)NumberofPoints/(float)binsize); 
    
    float *uncalflux_tmp = new float[3*binsize];
    float *elemindex_tmp = new float[3*binsize];  
    
    float *residuals_tmp = new float[binsize];  
    float *iindex = new float[binsize];     
    int *sindex = new int[binsize];    
	
    double *continuumFluxSample = new double[NumberOfSamples + 2]; 
    double *continuumElemSample = new double[NumberOfSamples + 2];
    double *continuumWavelengthSample = new double[NumberOfSamples + 2];
    
    unsigned actualNumberOfSamples = 0;
    
	for(unsigned k=0;k<NumberOfSamples;k++){
        
        float am,bm,abdevm;
		
		unsigned firstPoint,lastPoint;
        
        if(k==0) {
            firstPoint=0;
            lastPoint=3*binsize;
        } else if (k==NumberOfSamples-1) {
            firstPoint = NumberofPoints - 3*binsize;
            lastPoint = NumberofPoints;
        } else {
            firstPoint = (k-1)*binsize;
            lastPoint = (k+2)*binsize;
        }        
		
        unsigned np=0;
        for(unsigned i=firstPoint;i<lastPoint;i++)
        {
            uncalflux_tmp[np] = uncalibratedFlux.getflux(i);
            elemindex_tmp[np] = (double)i;
            np++;
        }	
        
		// 1st pass: calculate continuum slope      
        ladfit(elemindex_tmp,uncalflux_tmp,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
		
        //--- Clean up spectrum leaving only points within the box that deviate less than abdev from the robust linear fit 
        np=0;
        for(unsigned i=firstPoint;i<lastPoint;i++)
        {
            double fitMedianSlope = (bm*(double)i + am);
            
            if(fabs(uncalibratedFlux.getflux(i) - fitMedianSlope) < abdevm) {
                uncalflux_tmp[np] = uncalibratedFlux.getflux(i);
                elemindex_tmp[np] = (double)i;
                np++;
            }
        }	
        
		// 2nd pass: calculate continuum slope    
        if(np) {
            ladfit(elemindex_tmp,uncalflux_tmp,np,&am,&bm,&abdevm); // robust linear fit: f(x) =  a + b*x
        }
		
        firstPoint = k*binsize;
        lastPoint = (k+1)*binsize;
        if(lastPoint > NumberofPoints) {
            firstPoint = NumberofPoints - binsize;
            lastPoint = NumberofPoints;
        }
		
        // calculate residuals
        np = 0;
        for (unsigned i=firstPoint; i<lastPoint; i++) {
            float fitMedianSlope = (bm*(float)i + am);
            residuals_tmp[np] = uncalibratedFlux.getflux(i) - fitMedianSlope;
            iindex[np] = (float)i;
            np++;
        }
        
        // sort residuals
        operaArrayIndexSort(np,residuals_tmp,sindex);
        
        // select first maximum residual which does not exceed 5x(absolute deviation)         
        float dytop = 0;
        float continuumBinFraction = 0.5;
		// DT May 20 2014 fixed segfault in case that i goes less than zero...
        for(unsigned i=binsize-1; i>(unsigned)((float)binsize*(1.0 - continuumBinFraction)) && i >= 0; i--) {
            if (residuals_tmp[sindex[i]] < nsigcut*abdevm) {
                dytop = residuals_tmp[sindex[i]];
                break;
            }
        }
        
        // if none has been selected then assign 3*abdev       
        if (dytop == 0) { 
            dytop = nsigcut*abdevm; 
        }
        
        // evaluate max at left most edge
        if(k==0) {
            continuumElemSample[actualNumberOfSamples] = firstPoint;
            continuumFluxSample[actualNumberOfSamples] = bm*continuumElemSample[actualNumberOfSamples] + am + dytop; 
			// DT May 20 2014 check for out of bounds...
			// DT May 21 2014 change break to continue...
			if (bm*continuumElemSample[actualNumberOfSamples] + am + dytop > SpectralElements->getnSpectralElements())
				continue;
            continuumWavelengthSample[actualNumberOfSamples] = SpectralElements->getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]);
            
            spectralEnergyDistribution.setdistanceData(SpectralElements->getdistd((unsigned)continuumElemSample[actualNumberOfSamples]),actualNumberOfSamples);
            spectralEnergyDistribution.setwavelengthData(SpectralElements->getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]),actualNumberOfSamples);
            spectralEnergyDistribution.setfluxData(continuumFluxSample[actualNumberOfSamples],actualNumberOfSamples);
            spectralEnergyDistribution.setspectralEnergyData(1.0,actualNumberOfSamples);               
            
            actualNumberOfSamples++;
        }
		
        if (k==NumberOfSamples-1) {
			continuumElemSample[actualNumberOfSamples] = (double)lastPoint-1;
			continuumFluxSample[actualNumberOfSamples] = bm*continuumElemSample[actualNumberOfSamples] + am + dytop;
			// DT May 20 2014 check for out of bounds...
			// DT May 21 2014 change break to continue...
			if (bm*continuumElemSample[actualNumberOfSamples] + am + dytop > SpectralElements->getnSpectralElements())
				continue;
            continuumWavelengthSample[actualNumberOfSamples] = SpectralElements->getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]);
			
            spectralEnergyDistribution.setdistanceData(SpectralElements->getdistd((unsigned)continuumElemSample[actualNumberOfSamples]),actualNumberOfSamples);
            spectralEnergyDistribution.setwavelengthData(SpectralElements->getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]),actualNumberOfSamples);
            spectralEnergyDistribution.setfluxData(continuumFluxSample[actualNumberOfSamples],actualNumberOfSamples);
            spectralEnergyDistribution.setspectralEnergyData(1.0,actualNumberOfSamples);              
			
            actualNumberOfSamples++; 
        } else {
            continuumElemSample[actualNumberOfSamples] = firstPoint + (double)(binsize)/2.0;
            continuumFluxSample[actualNumberOfSamples] = bm*continuumElemSample[actualNumberOfSamples] + am + dytop;
            continuumWavelengthSample[actualNumberOfSamples] = SpectralElements->getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]);
            
            spectralEnergyDistribution.setdistanceData(SpectralElements->getdistd((unsigned)continuumElemSample[actualNumberOfSamples]),actualNumberOfSamples);
            spectralEnergyDistribution.setwavelengthData(SpectralElements->getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]),actualNumberOfSamples);
            spectralEnergyDistribution.setfluxData(continuumFluxSample[actualNumberOfSamples],actualNumberOfSamples);
            spectralEnergyDistribution.setspectralEnergyData(1.0,actualNumberOfSamples);           
            
            actualNumberOfSamples++; 
        }
	}
    
	spectralEnergyDistribution.setnDataPoints(actualNumberOfSamples);  
    
#ifdef PRINT_DEBUG  
    for(unsigned i=0; i<actualNumberOfSamples; i++){
        cout << i << " " << continuumWavelengthSample[i] << " " << continuumElemSample[i] << " " << continuumFluxSample[i] << endl;        
    }
#endif
    
    if(usePolynomial) {
        int npar = (int)orderOfPolynomial;
        double *par = new double[npar];
        double chisqr;
        
        operaLMFitPolynomial(actualNumberOfSamples,continuumElemSample,continuumFluxSample, npar, par, &chisqr);	
        
        Polynomial continuumPolynomial(npar, par);
        
        for(unsigned i=0;i<NumberofPoints;i++) {           
            outputContinuum.setflux(continuumPolynomial.Evaluate((double)i), i);
            outputContinuum.setvariance(uncalibratedFlux.getvariance(i), i);             
#ifdef PRINT_DEBUG           
            cout << i << " " 
            << uncalibratedFlux.getflux(i) << " "
            << outputContinuum.getflux(i) << " " 
            << outputContinuum.getvariance(i) << " " << endl;   
#endif            
        }
 		delete[] par;
    } else {
        double *continuumModelFlux = new double[NumberofPoints];
        double *continuumModelElem = new double[NumberofPoints];        
        for(unsigned i=0;i<NumberofPoints;i++) { 
            continuumModelElem[i] = (double)i;
        }
        
        operaFitSplineDouble(actualNumberOfSamples,continuumElemSample,continuumFluxSample,NumberofPoints,continuumModelElem,continuumModelFlux);
		
        for(unsigned i=0;i<NumberofPoints;i++) {
            outputContinuum.setflux(continuumModelFlux[i], i);
            outputContinuum.setvariance(uncalibratedFlux.getvariance(i), i);             
#ifdef PRINT_DEBUG          
            cout << i << " " 
            << uncalibratedFlux.getflux(i) << " "
            << outputContinuum.getflux(i) << " " 
            << outputContinuum.getvariance(i) << " " << endl;             
#endif            
        }        
        delete[] continuumModelFlux;
        delete[] continuumModelElem;          
    }
    
    delete[] continuumFluxSample; 
    delete[] continuumElemSample;
    delete[] continuumWavelengthSample;
    
    delete[] uncalflux_tmp;
    delete[] elemindex_tmp;  
    delete[] residuals_tmp;
    delete[] iindex;
    delete[] sindex;  
}

void operaSpectralOrder::calculateContinuum(unsigned binsize, unsigned nsigcut, ostream *poutspec, ostream *poutcontinuum) {
    
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    // use info in SpectralElements and binsize to create spectralEnergyDistribution
    
    createSpectralEnergyDistribution(binsize);
    
    SpectralEnergyDistribution->setUncalibratedFluxElements(SpectralElements);
    SpectralEnergyDistribution->calculateUncalibratedElements(binsize,nsigcut);
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            BeamElements[beam]->setdistd(SpectralElements->getdistd(indexElem), indexElem);
        }
        BeamElements[beam]->setwavelengthsFromCalibration(Wavelength);
        
        BeamSED[beam]->setUncalibratedFluxElements(BeamElements[beam]);
        BeamSED[beam]->calculateUncalibratedElements(binsize,nsigcut);
    }
        
    sethasSpectralEnergyDistribution(true);
    
    if (poutspec != NULL) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " "
            << SpectralEnergyDistribution->getFluxCalibrationElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFluxVariance(indexElem) << " " ;
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << beam << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFlux(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getThroughputElements()->getFlux(indexElem) << " "
                << BeamSED[beam]->getThroughputElements()->getFluxVariance(indexElem) << " ";
            }
            *poutspec << endl;
        }
        *poutspec << endl;
    }
	
    if (poutcontinuum != NULL) {
        for(unsigned i=0;i < SpectralEnergyDistribution->getnDataPoints(); i++) {
            *poutcontinuum << orderNumber << " "
			<< i << " "
			<< SpectralEnergyDistribution->getdistanceData(i) << " "
			<< SpectralEnergyDistribution->getwavelengthData(i) << " "
			<< SpectralEnergyDistribution->getfluxData(i) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutcontinuum << beam << " "
				<< BeamSED[beam]->getfluxData(i) << " ";
            }
			
            *poutcontinuum << endl;
        }
        *poutcontinuum << endl;
    }
    
}

void operaSpectralOrder::calculateFluxCalibrationFromExistingContinuum(unsigned nPointsInReference,double *refwl,double *refflux,double referenceFluxForNormalization,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS],double uncalibratedContinuumFluxForNormalization,double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS], ostream *poutspec, ostream *poutcontinuum) {
            
    SpectralEnergyDistribution->calculateCalibratedElements(nPointsInReference,refwl,refflux);

    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            BeamElements[beam]->setdistd(SpectralElements->getdistd(indexElem), indexElem);
        }
        BeamElements[beam]->setwavelengthsFromCalibration(Wavelength);

        BeamSED[beam]->calculateCalibratedElements(nPointsInReference,refwl,refflux);
        BeamSED[beam]->calculateSEDelements(BeamSpectralBinConstant[beam],referenceFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization[beam]);
    }
    
    SpectralEnergyDistribution->calculateSEDelements(spectralBinConstant,referenceFluxForNormalization,uncalibratedContinuumFluxForNormalization);
    
    sethasSpectralEnergyDistribution(true);
    
    if (poutspec != NULL) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " "
            << SpectralEnergyDistribution->getFluxCalibrationElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFluxVariance(indexElem) << " "
			<< referenceFluxForNormalization << " "
			<< uncalibratedContinuumFluxForNormalization << " " ;
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << beam << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFlux(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getThroughputElements()->getFlux(indexElem) << " "
                << BeamSED[beam]->getThroughputElements()->getFluxVariance(indexElem) << " "
                << uncalibratedContinuumBeamFluxForNormalization[beam] << " ";
            }
            *poutspec << endl;
        }
        *poutspec << endl;
    }
	
	
    if (poutcontinuum != NULL) {
        for(unsigned i=0;i < SpectralEnergyDistribution->getnDataPoints(); i++) {
            *poutcontinuum << orderNumber << " "
			<< i << " "
			<< SpectralEnergyDistribution->getdistanceData(i) << " "
			<< SpectralEnergyDistribution->getwavelengthData(i) << " "
			<< SpectralEnergyDistribution->getfluxData(i) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutcontinuum << beam << " "
				<< BeamSED[beam]->getfluxData(i) << " ";
            }
			
            *poutcontinuum << endl;
        }
        *poutcontinuum << endl;
    }
}

void operaSpectralOrder::calculateFluxCalibration(unsigned nPointsInReference,double *refwl,double *refflux, unsigned binsize,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS], ostream *poutspec, ostream *poutcontinuum) {
    
    double referenceFluxForNormalization = 1.0;
    double uncalibratedContinuumFluxForNormalization = 1.0;
    double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS];
    
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    SpectralElements->setwavelengthsFromCalibration(Wavelength);
    
    // use info in SpectralElements and binsize to create spectralEnergyDistribution
    createSpectralEnergyDistribution(binsize);

    unsigned nsigcut = 3;

    SpectralEnergyDistribution->setUncalibratedFluxElements(SpectralElements);
    SpectralEnergyDistribution->calculateCalibratedElements(nPointsInReference,refwl,refflux);
    SpectralEnergyDistribution->calculateUncalibratedElements(binsize,nsigcut);
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        uncalibratedContinuumBeamFluxForNormalization[beam] = 1.0;
        
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            BeamElements[beam]->setdistd(SpectralElements->getdistd(indexElem), indexElem);
        }
        BeamElements[beam]->setwavelengthsFromCalibration(Wavelength);
        
        BeamSED[beam]->setUncalibratedFluxElements(BeamElements[beam]);
        BeamSED[beam]->calculateCalibratedElements(nPointsInReference,refwl,refflux);
        BeamSED[beam]->calculateUncalibratedElements(binsize,nsigcut);
        BeamSED[beam]->calculateSEDelements(BeamSpectralBinConstant[beam],referenceFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization[beam]);
    }
    
    SpectralEnergyDistribution->calculateSEDelements(spectralBinConstant,referenceFluxForNormalization,uncalibratedContinuumFluxForNormalization);
    
    sethasSpectralEnergyDistribution(true);
    
    if (poutspec != NULL) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " "
            << SpectralEnergyDistribution->getFluxCalibrationElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFluxVariance(indexElem) << " "
			<< referenceFluxForNormalization << " "
			<< uncalibratedContinuumFluxForNormalization << " " ;
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << beam << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFlux(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getThroughputElements()->getFlux(indexElem) << " "
                << BeamSED[beam]->getThroughputElements()->getFluxVariance(indexElem) << " "
                << uncalibratedContinuumBeamFluxForNormalization[beam] << " ";
            }
            *poutspec << endl;
        }
        *poutspec << endl;
    }
	
    if (poutcontinuum != NULL) {
        for(unsigned i=0;i < SpectralEnergyDistribution->getnDataPoints(); i++) {
            *poutcontinuum << orderNumber << " "
			<< i << " "
			<< SpectralEnergyDistribution->getdistanceData(i) << " "
			<< SpectralEnergyDistribution->getwavelengthData(i) << " "
			<< SpectralEnergyDistribution->getfluxData(i) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutcontinuum << beam << " "
				<< BeamSED[beam]->getfluxData(i) << " ";
            }
			
            *poutcontinuum << endl;
        }
        *poutcontinuum << endl;
    }
    
}



void operaSpectralOrder::calculateFluxCalibration(unsigned nPointsInReference,double *refwl,double *refflux,double referenceFluxForNormalization, unsigned binsize,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS],double uncalibratedContinuumFluxForNormalization,double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS], ostream *poutspec, ostream *poutcontinuum) {
    
	if (binsize == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    SpectralElements->setwavelengthsFromCalibration(Wavelength);
    
    // use info in SpectralElements and binsize to create spectralEnergyDistribution
    createSpectralEnergyDistribution(binsize);
    
    unsigned nsigcut = 3;    
        
    SpectralEnergyDistribution->setUncalibratedFluxElements(SpectralElements);
    SpectralEnergyDistribution->calculateCalibratedElements(nPointsInReference,refwl,refflux);
    SpectralEnergyDistribution->calculateUncalibratedElements(binsize,nsigcut);

    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            BeamElements[beam]->setdistd(SpectralElements->getdistd(indexElem), indexElem);
        }
        BeamElements[beam]->setwavelengthsFromCalibration(Wavelength);
        
        BeamSED[beam]->setUncalibratedFluxElements(BeamElements[beam]);
        BeamSED[beam]->calculateCalibratedElements(nPointsInReference,refwl,refflux);
        BeamSED[beam]->calculateUncalibratedElements(binsize,nsigcut);
        BeamSED[beam]->calculateSEDelements(BeamSpectralBinConstant[beam],referenceFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization[beam]);
    }     
    
    SpectralEnergyDistribution->calculateSEDelements(spectralBinConstant,referenceFluxForNormalization,uncalibratedContinuumFluxForNormalization);
    
    sethasSpectralEnergyDistribution(true);    
    
    if (poutspec != NULL) {     
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " " 
            << indexElem << " " 
            << SpectralElements->getphotoCenterX(indexElem) << " " 
            << SpectralElements->getphotoCenterY(indexElem) << " " 
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << " "
            << SpectralEnergyDistribution->getFluxCalibrationElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFlux(indexElem) << " "
			<< SpectralEnergyDistribution->getThroughputElements()->getFluxVariance(indexElem) << " "
			<< referenceFluxForNormalization << " "
			<< uncalibratedContinuumFluxForNormalization << " " ;
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << beam << " "
                << BeamElements[beam]->getFlux(indexElem) << " "
                << BeamElements[beam]->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFlux(indexElem) << " "
                << BeamSED[beam]->getFluxCalibrationElements()->getFluxVariance(indexElem) << " "
                << BeamSED[beam]->getThroughputElements()->getFlux(indexElem) << " "            
                << BeamSED[beam]->getThroughputElements()->getFluxVariance(indexElem) << " "
                << uncalibratedContinuumBeamFluxForNormalization[beam] << " ";
            }
            *poutspec << endl;
        }  
        *poutspec << endl;  
    }          
	
    if (poutcontinuum != NULL) {
        for(unsigned i=0;i < SpectralEnergyDistribution->getnDataPoints(); i++) {        
            *poutcontinuum << orderNumber << " " 
			<< i << " "
			<< SpectralEnergyDistribution->getdistanceData(i) << " " 
			<< SpectralEnergyDistribution->getwavelengthData(i) << " " 
			<< SpectralEnergyDistribution->getfluxData(i) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutcontinuum << beam << " " 
				<< BeamSED[beam]->getfluxData(i) << " ";             
            }
			
            *poutcontinuum << endl;     
        }
        *poutcontinuum << endl;  
    }       
    
}

void operaSpectralOrder::normalizeSpectralElementsByConstant(double maxFluxForNormalization, double maxBeamFluxForNormalization[MAXNUMBEROFBEAMS]) {
    if (maxFluxForNormalization == 0) {
		throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    operaFluxVector normalizedFlux(nspecElem);
    normalizedFlux = (*SpectralElements->getFluxVector()) / maxFluxForNormalization;
    SpectralElements->setFluxVector(&normalizedFlux);

    for(unsigned beam = 0; beam < numberOfBeams; beam++) {
        if (maxBeamFluxForNormalization[beam] == 0) {
            throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
        }
        operaFluxVector normalizedBeamFlux(nspecElem);
        normalizedBeamFlux = (*BeamElements[beam]->getFluxVector()) / maxBeamFluxForNormalization[beam];
        BeamElements[beam]->setFluxVector(&normalizedBeamFlux);
    }
}

void operaSpectralOrder::divideSpectralElementsBySEDElements(bool useThroughput, ostream *poutspec, bool StarPlusSky, bool starplusskyInvertSkyFiber) {
    
    if(!gethasSpectralEnergyDistribution()) {
        throw operaException("operaSpectralOrder::applyFluxCalibration: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaSpectralElements *fluxCalibrationElements = SpectralEnergyDistribution->getFluxCalibrationElements();
    operaSpectralElements *thruputElements = SpectralEnergyDistribution->getThroughputElements();
    
    unsigned ncalElem = fluxCalibrationElements->getnSpectralElements();
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    if(ncalElem != nspecElem) {
		throw operaException("operaSpectralOrder:applyFluxCalibration: ncalElem != nspecElem; ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaFluxVector outputCalibratedFlux(nspecElem);
    
    if(useThroughput) {
        outputCalibratedFlux = (*SpectralElements->getFluxVector()) / (*thruputElements->getFluxVector()) ;
    } else {
        outputCalibratedFlux = (*SpectralElements->getFluxVector()) / (*fluxCalibrationElements->getFluxVector()) ;
    }
    
    SpectralElements->setFluxVector(&outputCalibratedFlux);
	
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        operaSpectralElements *fluxCalibrationBeamElements = BeamSED[beam]->getFluxCalibrationElements();
        operaSpectralElements *thruputBeamElements = BeamSED[beam]->getThroughputElements();
        
        operaFluxVector outputCalibratedBeamFlux(nspecElem);
        
        if(useThroughput) {
            outputCalibratedBeamFlux = (*BeamElements[beam]->getFluxVector()) / (*thruputBeamElements->getFluxVector()) ;
        } else {
            outputCalibratedBeamFlux = (*BeamElements[beam]->getFluxVector()) / (*fluxCalibrationBeamElements->getFluxVector()) ;
        }
        
        BeamElements[beam]->setFluxVector(&outputCalibratedBeamFlux);
    }
	
    if(StarPlusSky) {
        calculateStarAndSkyElements(starplusskyInvertSkyFiber,NULL);
    }
    
    if (poutspec != NULL) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
			<< SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
				<< BeamElements[beam]->getFluxVariance(indexElem) << " ";
            }
            *poutspec << endl;
        }
        *poutspec << endl;
    }
}


void operaSpectralOrder::multiplySpectralElementsBySEDElements(bool useThroughput,double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS], ostream *poutspec) {
    
    if(!gethasSpectralEnergyDistribution()) {
        throw operaException("operaSpectralOrder::applyFluxCalibration: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaSpectralElements *fluxCalibrationElements = SpectralEnergyDistribution->getFluxCalibrationElements();
    operaSpectralElements *thruputElements = SpectralEnergyDistribution->getThroughputElements();
    
    unsigned ncalElem = fluxCalibrationElements->getnSpectralElements();
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    if(ncalElem != nspecElem) {
		throw operaException("operaSpectralOrder:applyFluxCalibration: ncalElem != nspecElem; ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaFluxVector outputCalibratedFlux(nspecElem);
    
    if(useThroughput) {
        outputCalibratedFlux = (*SpectralElements->getFluxVector()) * (*thruputElements->getFluxVector()) / spectralBinConstant;
    } else {
        outputCalibratedFlux = (*SpectralElements->getFluxVector()) * (*fluxCalibrationElements->getFluxVector()) / spectralBinConstant;
    }
    SpectralElements->setFluxVector(&outputCalibratedFlux);
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        operaSpectralElements *fluxCalibrationBeamElements = BeamSED[beam]->getFluxCalibrationElements();
        operaSpectralElements *thruputBeamElements = BeamSED[beam]->getThroughputElements();
        
        operaFluxVector outputCalibratedBeamFlux(nspecElem);
        
        if(useThroughput) {
            outputCalibratedBeamFlux = (*BeamElements[beam]->getFluxVector()) * (*thruputBeamElements->getFluxVector()) / BeamSpectralBinConstant[beam];
        } else {
            outputCalibratedBeamFlux = (*BeamElements[beam]->getFluxVector()) * (*fluxCalibrationBeamElements->getFluxVector()) / BeamSpectralBinConstant[beam];
        }
        
        BeamElements[beam]->setFluxVector(&outputCalibratedBeamFlux);
    }
	
    if (poutspec != NULL) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
			<< SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
				<< BeamElements[beam]->getFluxVariance(indexElem) << " ";
            }
            *poutspec << endl;
        }
        *poutspec << endl;
    }
}

/*
 * Same as multiplySpectralElementsBySEDElements but without changing the beams
 */
void operaSpectralOrder::multiplySpectralElementsBySEDElements(bool useThroughput,double spectralBinConstant, ostream *poutspec) {
    
    if(!gethasSpectralEnergyDistribution()) {
        throw operaException("operaSpectralOrder::applyFluxCalibration: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaSpectralElements *fluxCalibrationElements = SpectralEnergyDistribution->getFluxCalibrationElements();
    operaSpectralElements *thruputElements = SpectralEnergyDistribution->getThroughputElements();
    
    unsigned ncalElem = fluxCalibrationElements->getnSpectralElements();
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    if(ncalElem != nspecElem) {
		throw operaException("operaSpectralOrder:applyFluxCalibration: ncalElem != nspecElem; ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaFluxVector outputCalibratedFlux(nspecElem);
    
    if(useThroughput) {
        outputCalibratedFlux = (*SpectralElements->getFluxVector()) * ((*thruputElements->getFluxVector())) / (spectralBinConstant);
    } else {
        outputCalibratedFlux = (*SpectralElements->getFluxVector()) * ((*fluxCalibrationElements->getFluxVector()))  / (spectralBinConstant) ;
    }
    
    SpectralElements->setFluxVector(&outputCalibratedFlux);
	
    if (poutspec != NULL) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
			<< SpectralElements->getFluxVariance(indexElem) << endl;
        }
        *poutspec << endl;
    }
}

void operaSpectralOrder::applyFluxCalibration(double spectralBinConstant,double BeamSpectralBinConstant[MAXNUMBEROFBEAMS],double uncalibratedContinuumFluxForNormalization,double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS], bool absoluteCalibration, ostream *poutspec) {
    
    if(!gethasSpectralEnergyDistribution()) {
        throw operaException("operaSpectralOrder::applyFluxCalibration: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaSpectralElements *fluxCalibrationElements = SpectralEnergyDistribution->getFluxCalibrationElements();
    operaSpectralElements *thruputElements = SpectralEnergyDistribution->getThroughputElements();
    
    unsigned ncalElem = fluxCalibrationElements->getnSpectralElements();
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    if(ncalElem != nspecElem) {
		throw operaException("operaSpectralOrder:applyFluxCalibration: ncalElem != nspecElem; ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaFluxVector outputCalibratedFlux(nspecElem);
    
    if(absoluteCalibration) {
        outputCalibratedFlux = ((*SpectralElements->getFluxVector())*(spectralBinConstant) / (*fluxCalibrationElements->getFluxVector()));
    } else {
        outputCalibratedFlux = (*SpectralElements->getFluxVector()) / ((*thruputElements->getFluxVector())*(uncalibratedContinuumFluxForNormalization)) ;
    }
    
    SpectralElements->setFluxVector(&outputCalibratedFlux);
	
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        operaSpectralElements *fluxCalibrationBeamElements = BeamSED[beam]->getFluxCalibrationElements();
        operaSpectralElements *thruputBeamElements = BeamSED[beam]->getThroughputElements();
        
        operaFluxVector outputCalibratedBeamFlux(nspecElem);
        
        if(absoluteCalibration) {
            outputCalibratedBeamFlux = ((*BeamElements[beam]->getFluxVector())*(BeamSpectralBinConstant[beam]) / (*fluxCalibrationBeamElements->getFluxVector()));
        } else {
            outputCalibratedBeamFlux = (*BeamElements[beam]->getFluxVector()) / ((*thruputBeamElements->getFluxVector())*(uncalibratedContinuumBeamFluxForNormalization[beam])) ;
        }
        
        BeamElements[beam]->setFluxVector(&outputCalibratedBeamFlux);
    }
	
    if (poutspec != NULL) {
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " "
            << indexElem << " "
            << SpectralElements->getphotoCenterX(indexElem) << " "
            << SpectralElements->getphotoCenterY(indexElem) << " "
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
			<< SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
				<< BeamElements[beam]->getFluxVariance(indexElem) << " ";
            }
            *poutspec << endl;
        }
        *poutspec << endl;
    }
}


void operaSpectralOrder::applyFluxCalibration(double exposureTime, ostream *poutspec) {
    
    if(!gethasSpectralEnergyDistribution()) {
        throw operaException("operaSpectralOrder::applyFluxCalibration: ",operaErrorHasNoSpectralElements, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaSpectralElements *fluxCalibrationElements = SpectralEnergyDistribution->getFluxCalibrationElements();
    
    unsigned ncalElem = fluxCalibrationElements->getnSpectralElements();
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    if(ncalElem != nspecElem) {
		throw operaException("operaSpectralOrder:applyFluxCalibration: ncalElem != nspecElem; ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	        
    }
    
    operaFluxVector outputCalibratedFlux(nspecElem);     
    outputCalibratedFlux = ((*SpectralElements->getFluxVector())*(exposureTime) / (*fluxCalibrationElements->getFluxVector()));
    SpectralElements->setFluxVector(&outputCalibratedFlux);
	
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        operaSpectralElements *fluxCalibrationBeamElements = BeamSED[beam]->getFluxCalibrationElements();        
        operaFluxVector outputCalibratedBeamFlux(nspecElem);
        outputCalibratedBeamFlux = ((*BeamElements[beam]->getFluxVector())*(exposureTime) / (*fluxCalibrationBeamElements->getFluxVector()));
        BeamElements[beam]->setFluxVector(&outputCalibratedBeamFlux);
    }    
	
    if (poutspec != NULL) {     
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " " 
            << indexElem << " " 
            << SpectralElements->getphotoCenterX(indexElem) << " " 
            << SpectralElements->getphotoCenterY(indexElem) << " " 
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
			<< SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
				<< BeamElements[beam]->getFluxVariance(indexElem) << " ";             
            }
            *poutspec << endl;
        }  
        *poutspec << endl;  
    }           
}

/*
 *  Apply the LE Moon Flat Response for LE-compatible calibration.
 *  The fluxCalibrationElements contain wavelength and flux factor for all wavelengths,
 *  so the lengths won't match. We need to get the wavelength range and then create
 *  a vector of the matching length, by interpolation.
 */
void operaSpectralOrder::applyFlatResponse(double exposureTime, operaSpectralElements *fluxCalibrationElements, ostream *poutspec) {
    
    
    unsigned ncalElem = fluxCalibrationElements->getnSpectralElements();
    unsigned nspecElem = SpectralElements->getnSpectralElements();
	double wl0 = SpectralElements->getwavelength(0);
	double wlf = SpectralElements->getwavelength(nspecElem-1);
	
	unsigned firstwl = 0;
	unsigned lastwl = 0;
	
	for (unsigned i=0; i<ncalElem; i++) {
		if (fluxCalibrationElements->getwavelength(i) < wl0) {
			firstwl = i;
		}
		if (fluxCalibrationElements->getwavelength(i) < wlf) {
			lastwl = i;
		}
	}
	unsigned length = lastwl - firstwl + 1;
	double *inxs = new double[length];
	double *outxs = new double[nspecElem];
    for (unsigned i=0; i<length; i++) {
		inxs[i] = i;
	}
	double *yin = fluxCalibrationElements->getFluxVector()->getfluxes() + firstwl;
	double *yout = new double[nspecElem];
    operaFitSplineDouble(length, inxs, yin, nspecElem, outxs, yout);
	operaFluxVector calFlux(nspecElem);
    for (unsigned i=0; i<nspecElem; i++) {
		calFlux.setflux(yout[i], i);
	}
    operaFluxVector outputCalibratedFlux(nspecElem);     
    outputCalibratedFlux = ((*SpectralElements->getFluxVector())/calFlux)/exposureTime;
    SpectralElements->setFluxVector(&outputCalibratedFlux);
	
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        operaFluxVector outputCalibratedBeamFlux(nspecElem);     
        outputCalibratedBeamFlux = ((*BeamElements[beam]->getFluxVector())/calFlux)/exposureTime;
        BeamElements[beam]->setFluxVector(&outputCalibratedBeamFlux);
    }
    
	delete yout;
	delete inxs;
	delete outxs;
	
    if (poutspec != NULL) {     
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " " 
            << indexElem << " " 
            << SpectralElements->getphotoCenterX(indexElem) << " " 
            << SpectralElements->getphotoCenterY(indexElem) << " " 
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "
            << SpectralElements->getFlux(indexElem) << " "
			<< SpectralElements->getFluxVariance(indexElem) << " ";
            
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                *poutspec << BeamElements[beam]->getFlux(indexElem) << " "
				<< BeamElements[beam]->getFluxVariance(indexElem) << " ";             
            }
            *poutspec << endl;
        }  
        *poutspec << endl;  
    }           
}

/*
 * Star+Sky Mode, subtract the sky beam
 */
void operaSpectralOrder::calculateStarAndSkyElements(bool starplusskyInvertSkyFiber, ostream *poutspec) {
    
    // Test if number of beams is even.
    if(float(numberOfBeams)/2 - float(round(float(numberOfBeams)/2))) { 
        throw operaException("operaSpectralOrder:calculateStarMinusSky: numberOfBeams must be even;", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    createSkyElements(nspecElem,SpectrumType);    
    
    operaFluxVector skyFlux(nspecElem);
    operaFluxVector starPlusSkyFlux(nspecElem);
    
    // Below it assumes there is only two beams, beam=0 is star+sky and beam=1 is sky
    if(starplusskyInvertSkyFiber) {
        skyFlux = (*BeamElements[0]->getFluxVector());
        starPlusSkyFlux = (*BeamElements[1]->getFluxVector());
    } else {
        starPlusSkyFlux = (*BeamElements[0]->getFluxVector());
        skyFlux = (*BeamElements[1]->getFluxVector());
    }
    
    // Below it assumes that the first half of beams is object+sky 
    // and the second half is sky
    /*
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {   
        if(beam < floor(float(numberOfBeams)/2)) {
            starPlusSkyFlux += (*BeamElements[beam]->getFluxVector());  
        } else if(float(beam) >= float(numberOfBeams)/2) {
            skyFlux += (*BeamElements[beam]->getFluxVector());  
        }
    } */
    
    for(unsigned indexElem=0;indexElem<nspecElem;indexElem++) {
        SkyElements->setphotoCenter(SpectralElements->getphotoCenterX(indexElem),SpectralElements->getphotoCenterY(indexElem),indexElem);        
        SkyElements->setdistd(SpectralElements->getdistd(indexElem),indexElem);
        SkyElements->setwavelength(SpectralElements->getwavelength(indexElem),indexElem);
    }
    
    operaFluxVector starFlux(nspecElem);
        
    // Below is a workaround since operaFluxVector routines can't handle negative values appropriately
    for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
        // E. Martioli Oct 31, 2014. Commented line below cause sky subtraction is causing alias
        //  on main extracted spectra. The beams are saved separately in case someone wants to
        //  subtract sky on a later analysis. The way around it is to use robust statistics to
        //  obtain the sky flux, but this would have to be implemented in the optimal extraction
        //  routines, which we don't wanna touch now. 
        double skysubtractedflux = starPlusSkyFlux.getflux(indexElem) - skyFlux.getflux(indexElem);
        //double skysubtractedflux = starPlusSkyFlux.getflux(indexElem);
        double variance = starPlusSkyFlux.getvariance(indexElem) + skyFlux.getvariance(indexElem);
        starFlux.setflux(skysubtractedflux,indexElem);
        starFlux.setvariance(variance,indexElem);
    }
    
    SpectralElements->setFluxVector(&starFlux);
    
    SkyElements->setFluxVector(&skyFlux);
    
    sethasSkyElements(true);    
    
    if (poutspec != NULL) {     
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " " 
            << indexElem << " " 
            << SpectralElements->getphotoCenterX(indexElem) << " " 
            << SpectralElements->getphotoCenterY(indexElem) << " "             
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "            
            << starPlusSkyFlux.getflux(indexElem) << " "
            << starPlusSkyFlux.getvariance(indexElem) << " "           
            << SkyElements->getFlux(indexElem) << " "
            << SkyElements->getFluxVariance(indexElem) << " "           
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << endl;
        }  
    }      
}

/*
 * Polar Mode, take the sum of the two beams
 */
void operaSpectralOrder::calculatePolarElements(ostream *poutspec) {
    
    // Test if number of beams is even.
    if(float(numberOfBeams)/2 - float(round(float(numberOfBeams)/2))) { 
        throw operaException("operaSpectralOrder:calculateStarMinusSky: numberOfBeams must be even;", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned nspecElem = SpectralElements->getnSpectralElements();
    
    createSkyElements(nspecElem,SpectrumType);    
    
    operaFluxVector Beam1Flux(nspecElem);
    operaFluxVector Beam2Flux(nspecElem);
    
    // This routine assumes that the first half of beams is sky 
    // and the second half is object+sky
    
    for(unsigned beam = 0; beam < numberOfBeams; beam++) {   
        if(beam < floor(float(numberOfBeams)/2)) {
            Beam1Flux += (*BeamElements[beam]->getFluxVector());  
        } else if(float(beam) >= float(numberOfBeams)/2) {
            Beam2Flux += (*BeamElements[beam]->getFluxVector());  
        }
    }        
    
    for(unsigned indexElem=0;indexElem<nspecElem;indexElem++) {
        SkyElements->setphotoCenter(SpectralElements->getphotoCenterX(indexElem),SpectralElements->getphotoCenterY(indexElem),indexElem);        
        SkyElements->setdistd(SpectralElements->getdistd(indexElem),indexElem);
        SkyElements->setwavelength(SpectralElements->getwavelength(indexElem),indexElem);
    }
    
    operaFluxVector starFlux(nspecElem);
    
    starFlux = Beam1Flux + Beam2Flux;
    
    SpectralElements->setFluxVector(&starFlux);
    
    if (poutspec != NULL) {     
        for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
            *poutspec << orderNumber << " " 
            << indexElem << " " 
            << SpectralElements->getphotoCenterX(indexElem) << " " 
            << SpectralElements->getphotoCenterY(indexElem) << " "             
            << SpectralElements->getdistd(indexElem) << " "
            << SpectralElements->getwavelength(indexElem) << " "            
            << SpectralElements->getFlux(indexElem) << " "
            << SpectralElements->getFluxVariance(indexElem) << endl;
        }  
    }      
}

/*
 * Radial Velocity Wavelength Correction
 */
void operaSpectralOrder::applyRvelWavelengthCorrection(double RVcorrectionInKmPerSecond) {
    if(gethasWavelength() && SpectralElements->getHasDistance() && !SpectralElements->getHasWavelength()) {
        SpectralElements->setwavelengthsFromCalibration(getWavelength());
    }
    
    for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
        double currentElemWavelength = SpectralElements->getwavelength(indexElem);
        double wavelengthCorrection = currentElemWavelength*(RVcorrectionInKmPerSecond/SPEED_OF_LIGHT_KMS);
        SpectralElements->setwavelength(currentElemWavelength+wavelengthCorrection,indexElem);
    }
}

void operaSpectralOrder::setExtendedRvelWavelengthCorrection(double RVcorrectionInKmPerSecond) {
    
    if(!gethasWavelength() || !SpectralElements->getHasDistance() || !SpectralElements->getHasExtendedBeamFlux()) {
        throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
        double currentElemWavelength = SpectralElements->getwavelength(indexElem);
        double wavelengthCorrection = currentElemWavelength*(RVcorrectionInKmPerSecond/SPEED_OF_LIGHT_KMS);
        SpectralElements->setrvel(wavelengthCorrection,indexElem);
    }
}

void operaSpectralOrder::applyWavelengthCorrectionFromExtendedRvel(void) {
    
    if(!SpectralElements->getHasExtendedBeamFlux()) {
        throw operaException("operaSpectralOrder: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
        SpectralElements->setwavelength(SpectralElements->getwavelength(indexElem) + SpectralElements->getrvel(indexElem),indexElem);
    }
}

void operaSpectralOrder::CreateExtendedVectors() {
	SpectralElements->CreateExtendedvectors(SpectralElements->getnSpectralElements());
	for(unsigned beam = 0; beam < numberOfBeams; beam++) {
		BeamElements[beam]->CreateExtendedvectors(BeamElements[beam]->getnSpectralElements());
	}
}

void operaSpectralOrder::CopyFluxVectorIntoRawFlux() {
	SpectralElements->copyTOrawFlux();
	for(unsigned beam = 0; beam < numberOfBeams; beam++) BeamElements[beam]->copyTOrawFlux();
}

void operaSpectralOrder::CopyRawFluxIntoFluxVector() {
	SpectralElements->copyFROMrawFlux();
	for(unsigned beam = 0; beam < numberOfBeams; beam++) BeamElements[beam]->copyFROMrawFlux();
}

void operaSpectralOrder::CopyFluxVectorIntoFcalFlux() {
	SpectralElements->copyTOfcalFlux();
	for(unsigned beam = 0; beam < numberOfBeams; beam++) BeamElements[beam]->copyTOfcalFlux();
}

void operaSpectralOrder::CopyFcalFluxIntoFluxVector() {
	SpectralElements->copyFROMfcalFlux();
	for(unsigned beam = 0; beam < numberOfBeams; beam++) BeamElements[beam]->copyFROMfcalFlux();
}

void operaSpectralOrder::CopyFluxVectorIntoNormalizedFlux() {
	SpectralElements->copyTOnormalizedFlux();
	for(unsigned beam = 0; beam < numberOfBeams; beam++) BeamElements[beam]->copyTOnormalizedFlux();
}

void operaSpectralOrder::CopyNormalizedFluxIntoFluxVector() {
	SpectralElements->copyFROMnormalizedFlux();
	for(unsigned beam = 0; beam < numberOfBeams; beam++) BeamElements[beam]->copyFROMnormalizedFlux();
}

void operaSpectralOrder::TrimOrderToWavelengthRange(double wl0, double wlf) {
    unsigned nElements = SpectralElements->getnSpectralElements();
    unsigned newIndexElem = 0;
	for(unsigned indexElem=0; indexElem<nElements; indexElem++) {
		double wl = SpectralElements->getwavelength(indexElem);
		if (wl0 <= wl && wl <= wlf) {
			SpectralElements->setwavelength(wl,newIndexElem);
			double dist = SpectralElements->getdistd(indexElem);
			double flux = SpectralElements->getFlux(indexElem);
			double fluxVariance = SpectralElements->getFluxVariance(indexElem);
			SpectralElements->setdistd(dist,newIndexElem);
			SpectralElements->setFlux(flux,newIndexElem);
			SpectralElements->setFluxVariance(fluxVariance,newIndexElem);
			if(SpectralElements->getHasExtendedBeamFlux()) {
				double tell = SpectralElements->gettell(indexElem);
				double rvel = SpectralElements->getrvel(indexElem);
				double normalizedFlux = SpectralElements->getnormalizedFlux(indexElem);
				double normalizedFluxVariance = SpectralElements->getnormalizedFluxVariance(indexElem);
				double fcalFlux = SpectralElements->getfcalFlux(indexElem);
				double fcalFluxVariance = SpectralElements->getfcalFluxVariance(indexElem);
				double rawFlux = SpectralElements->getrawFlux(indexElem);
				double rawFluxVariance = SpectralElements->getrawFluxVariance(indexElem);
				SpectralElements->settell(tell,newIndexElem);
				SpectralElements->setrvel(rvel,newIndexElem);
				SpectralElements->setnormalizedFlux(normalizedFlux,newIndexElem);
				SpectralElements->setnormalizedFluxVariance(normalizedFluxVariance,newIndexElem);
				SpectralElements->setfcalFlux(fcalFlux,newIndexElem);
				SpectralElements->setfcalFluxVariance(fcalFluxVariance,newIndexElem);
				SpectralElements->setrawFlux(rawFlux,newIndexElem);
				SpectralElements->setrawFluxVariance(rawFluxVariance,newIndexElem);
			}
			for(unsigned i = 0; i < numberOfBeams; i++) {
				double beamFlux = BeamElements[i]->getFlux(indexElem);
				double beamFluxVariance = BeamElements[i]->getFluxVariance(indexElem);
				BeamElements[i]->setFlux(beamFlux,newIndexElem);
				BeamElements[i]->setFluxVariance(beamFluxVariance,newIndexElem);
				if(BeamElements[i]->getHasExtendedBeamFlux()) {
					double normalizedFlux = BeamElements[i]->getnormalizedFlux(indexElem);
					double normalizedFluxVariance = BeamElements[i]->getnormalizedFluxVariance(indexElem);
					double fcalFlux = BeamElements[i]->getfcalFlux(indexElem);
					double fcalFluxVariance = BeamElements[i]->getfcalFluxVariance(indexElem);
					double rawFlux = BeamElements[i]->getrawFlux(indexElem);
					double rawFluxVariance = BeamElements[i]->getrawFluxVariance(indexElem);
					BeamElements[i]->setnormalizedFlux(normalizedFlux,newIndexElem);
					BeamElements[i]->setnormalizedFluxVariance(normalizedFluxVariance,newIndexElem);
					BeamElements[i]->setfcalFlux(fcalFlux,newIndexElem);
					BeamElements[i]->setfcalFluxVariance(fcalFluxVariance,newIndexElem);
					BeamElements[i]->setrawFlux(rawFlux,newIndexElem);
					BeamElements[i]->setrawFluxVariance(rawFluxVariance,newIndexElem);
				}
			}
			if (gethasPolarimetry()) {
				stokes_parameter_t stokesParameter = StokesI;
				if (Polarimetry->getHasStokesV()) stokesParameter = StokesV;
				else if (Polarimetry->getHasStokesQ()) stokesParameter = StokesQ;
				else if (Polarimetry->getHasStokesU()) stokesParameter = StokesU;
				else if (Polarimetry->getHasStokesI()) stokesParameter = StokesI;
				double degree = Polarimetry->getDegreeOfPolarization(stokesParameter)->getflux(indexElem);
				double degreeVariance = Polarimetry->getDegreeOfPolarization(stokesParameter)->getvariance(indexElem);
				double firstNull = Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(indexElem);
				double firstNullVariance = Polarimetry->getFirstNullPolarization(stokesParameter)->getvariance(indexElem);
				double secondNull = Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(indexElem);
				double secondNullVariance = Polarimetry->getSecondNullPolarization(stokesParameter)->getvariance(indexElem);
				Polarimetry->getDegreeOfPolarization(stokesParameter)->setflux(degree, newIndexElem);
				Polarimetry->getDegreeOfPolarization(stokesParameter)->setvariance(degreeVariance, newIndexElem);
				Polarimetry->getFirstNullPolarization(stokesParameter)->setflux(firstNull, newIndexElem);
				Polarimetry->getFirstNullPolarization(stokesParameter)->setvariance(firstNullVariance, newIndexElem);
				Polarimetry->getSecondNullPolarization(stokesParameter)->setflux(secondNull, newIndexElem);
				Polarimetry->getSecondNullPolarization(stokesParameter)->setvariance(secondNullVariance, newIndexElem);
			}
			newIndexElem++;
		} else if (wl > wl0 && wl > wlf) {
			break;
		}
	}
	
	SpectralElements->resize(newIndexElem);
	for(unsigned i = 0; i < numberOfBeams; i++) BeamElements[i]->resize(newIndexElem);
	if (gethasPolarimetry()) Polarimetry->resize(newIndexElem);
}

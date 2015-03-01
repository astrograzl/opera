/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralOrderVector
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

#include <math.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaMEFFITSProduct.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaGeometry.h"		// for calculate order length
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaWavelength.h" // MAXORDEROFWAVELENGTHPOLYNOMIAL and MAXREFWAVELENGTHSPERORDER
#include "libraries/operaSpectralLines.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/GainBiasNoise.h"
#include "libraries/Polynomial.h"
#include "libraries/LaurentPolynomial.h"  
#include "libraries/operaExtractionAperture.h"
#include "libraries/operaGeometricShapes.h"
#include "libraries/operaStokesVector.h"
#include "libraries/operaPolarimetry.h"

#include "libraries/operastringstream.h"
#include "libraries/operaLibCommon.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaFit.h"
#include "libraries/operaStats.h"
#include "libraries/gzstream.h"
#include "libraries/operaFFT.h"    
#include "libraries/ladfit.h" // for ladfit_d

#include "libraries/operaSpectralTools.h"			// void calculateUniformSample, getFluxAtWavelength



/*!
 * operaSpectralOrderVector
 * \author Doug Teeple
 * \brief spectral order vector.
 * \details {This library contains serialization and deserializatin of spectral orders.}
 * \file operaSpectralOrderVector.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectralOrderVector();
 * \brief Base constructor.
 */
operaSpectralOrderVector::operaSpectralOrderVector() :
vector(NULL),
orderSpacingPolynomial(NULL),
length(0),
minorder(0),
maxorder(0),
sequence(0),
instrumentmode(MODE_UNKNOWN),
count(0),
gainBiasNoise(NULL),
BarycentricRadialVelocityCorrection(0.0),
first(NULL),
last(NULL)
{
	unsigned order = 0;
	vector = (operaSpectralOrder **)malloc(sizeof(operaSpectralOrder *) * (MAXORDERS+1));
	if (!vector) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (order=minorder; order<MAXORDERS; order++) {
		vector[order] = new operaSpectralOrder(order);
	}
	vector[order] = NULL;
	length = MAXORDERS;
	orderSpacingPolynomial = new Polynomial();
    for (unsigned dispIndex=0; dispIndex<MAXORDEROFWAVELENGTHPOLYNOMIAL; dispIndex++) {
        dispersionPolynomial[dispIndex] = new LaurentPolynomial();
    }	
	gainBiasNoise = new GainBiasNoise();
}
/* 
 * \class operaSpectralOrderVector(unsigned length, unsigned maxdatapoints, unsigned maxValues, unsigned nElements);
 * \brief Create a NULL-terminated SpectralOrderVector of spectralorders of type "None".
 */
operaSpectralOrderVector::operaSpectralOrderVector(unsigned Length, unsigned maxdatapoints, unsigned maxValues, unsigned nElements) :
vector(NULL),
orderSpacingPolynomial(NULL),
length(0),
minorder(0),
maxorder(0),
sequence(0),
instrumentmode(MODE_UNKNOWN),
count(0),
gainBiasNoise(NULL),
BarycentricRadialVelocityCorrection(0.0),
first(NULL),
last(NULL)
{
	if (Length == 0) {
		throw operaException("operaSpectralOrderVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (Length > MAXORDERS) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	unsigned order = 0;
	length = Length;
	vector = (operaSpectralOrder **)malloc(sizeof(operaSpectralOrder *) * (MAXORDERS+1));
	if (!vector) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (order=minorder; order<MAXORDERS; order++) {
		vector[order] = new operaSpectralOrder(order, maxdatapoints, maxValues, nElements, None);
	}
	vector[order] = NULL;
	orderSpacingPolynomial = new Polynomial();
    for (unsigned dispIndex=0; dispIndex<MAXORDEROFWAVELENGTHPOLYNOMIAL; dispIndex++) {
        dispersionPolynomial[dispIndex] = new LaurentPolynomial();
    }
	gainBiasNoise = new GainBiasNoise();
}
/* 
 * \class operaSpectralOrderVector(string Filename);
 * \brief Base constructor, read a spectral order vector from a filename
 * \brief Filename can contain either a .geom file or a .wcal file as given by the format identifier
 * \brief both of which are stored as vectors of polynomial coefficients.
 * \param Filename - string Filename to read
 */
operaSpectralOrderVector::operaSpectralOrderVector(string Filename) :
vector(NULL),
orderSpacingPolynomial(NULL),
length(0),
minorder(0),
maxorder(0),
sequence(0),
instrumentmode(MODE_UNKNOWN),
count(0),
gainBiasNoise(NULL),
BarycentricRadialVelocityCorrection(0.0),
first(NULL),
last(NULL)
{
	unsigned order = 0;
	vector = (operaSpectralOrder **)malloc(sizeof(operaSpectralOrder *) * (MAXORDERS+1));
	if (!vector) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (order=minorder; order<MAXORDERS; order++) {
		vector[order] = new operaSpectralOrder(order);
	}
	vector[order] = NULL;
	length = MAXORDERS;
	orderSpacingPolynomial = new Polynomial();
    for (unsigned dispIndex=0; dispIndex<MAXORDEROFWAVELENGTHPOLYNOMIAL; dispIndex++) {
        dispersionPolynomial[dispIndex] = new LaurentPolynomial();
    }
	gainBiasNoise = new GainBiasNoise();
	
	ReadSpectralOrders(Filename);
	first = vector;
	last = first + count;
}
/*
 * Destructor
 */
operaSpectralOrderVector::~operaSpectralOrderVector() {
	freeSpectralOrderVector();
	delete orderSpacingPolynomial;
	orderSpacingPolynomial = NULL;
    for (unsigned dispIndex=0; dispIndex<MAXORDEROFWAVELENGTHPOLYNOMIAL; dispIndex++) {
        delete dispersionPolynomial[dispIndex];
        dispersionPolynomial[dispIndex] = NULL;
    }
	
	delete gainBiasNoise;
	gainBiasNoise = NULL;
	first = NULL;
	last = NULL;
}

/*
 * Methods
 */

/*!
 * unsigned getnumberOfDispersionPolynomials(void);
 * \brief returns the number of dispersion polynomials.
 * \return unsigned - numberOfDispersionPolynomials.
 */
unsigned operaSpectralOrderVector::getnumberOfDispersionPolynomials(void) {
    return numberOfDispersionPolynomials;
}

/*!
 * void setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials);
 * \brief sets the number of dispersion polynomials.
 * \return none.
 */
void operaSpectralOrderVector::setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials) {\
	if (NumberOfDispersionPolynomials > MAXORDEROFWAVELENGTHPOLYNOMIAL) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    numberOfDispersionPolynomials = NumberOfDispersionPolynomials;
}

void operaSpectralOrderVector::freeSpectralOrderVector() {
	if (vector) {
		operaSpectralOrder **vec = vector;
		while (operaSpectralOrder *spectralOrder = *vec++) {
			delete spectralOrder;
		}
		free(vector);
		vector = NULL;
		
	}
}
/* 
 * unsigned getGainBiasNoise();
 * \brief returns a pointer to the GainBiasNoise class instance.
 */
GainBiasNoise *operaSpectralOrderVector::getGainBiasNoise() {
	return gainBiasNoise;
}

/* 
 * unsigned getBarycentricRadialVelocityCorrection();
 * \brief returns a double BarycentricRadialVelocityCorrection.
 */
double operaSpectralOrderVector::getBarycentricRadialVelocityCorrection() {
	return BarycentricRadialVelocityCorrection;
}

/* 
 * void getBarycentricRadialVelocityCorrection(double BarycentricRadialVelocityCorrection);
 * \brief sets the double BarycentricRadialVelocityCorrection.
 */
void operaSpectralOrderVector::setBarycentricRadialVelocityCorrection(double theBarycentricRadialVelocityCorrection) {
	BarycentricRadialVelocityCorrection = theBarycentricRadialVelocityCorrection;
}

/* 
 * unsigned getCount();
 * \brief returns the count of spectral orders that have content.
 */
unsigned operaSpectralOrderVector::getCount() {
	return count;
}

/* 
 * unsigned getMinorder();
 * \brief returns the least order number in the vector
 */
unsigned operaSpectralOrderVector::getMinorder() {
	return minorder;
}

/* 
 * unsigned getMaxorder();
 * \brief returns the maximal order number in the vector
 */
unsigned operaSpectralOrderVector::getMaxorder() {
	return maxorder;
}

/* 
 * unsigned setCount();
 * \brief sets the count of spectral orders that have content.
 */
void operaSpectralOrderVector::setCount(unsigned Count) {
#ifdef RANGE_CHECK
	if (Count > length) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	count = Count;
}

/* 
 * unsigned setMinorder();
 * \brief sets the least order number in the vector
 */
void operaSpectralOrderVector::setMinorder(unsigned Minorder) {
#ifdef RANGE_CHECK
	if (Minorder > length) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	minorder = Minorder;
}

/* 
 * void getMaxorder(unsigned Maxorder);
 * \brief sets the maximal order number in the vector
 */
void operaSpectralOrderVector::setMaxorder(unsigned Maxorder) {
#ifdef RANGE_CHECK
	if (Maxorder > length) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	maxorder = Maxorder;
}

/*
 * \sa void setObject(string object)
 * \brief sets the object name
 * \return none.
 */
void operaSpectralOrderVector::setObject(string Object) {
	object = Object;
}

/*
 * \sa string getObject(void);
 * \brief get the object name
 * \return none.
 */
string operaSpectralOrderVector::getObject(void) {
	return object;
}

/* 
 * \sa setSequence(unsigned sequence)
 * \brief sets the sequence number
 * \return none.
 */
void operaSpectralOrderVector::setSequence(unsigned Sequence) {
	sequence = Sequence;
}

/* 
 * \sa unsigned getSequence(void);
 * \brief get the sequence number
 * \return none.
 */
unsigned operaSpectralOrderVector::getSequence(void) {
	return sequence;
}
/* 
 * Polynomial *getOrderSpacingPolynomial(void);
 * \brief gets the order spacing polynomial
 */
Polynomial *operaSpectralOrderVector::getOrderSpacingPolynomial(void) {
	return orderSpacingPolynomial;
}

/* 
 * setOrderSpacingPolynomial(PolynomialCoeffs_t *pc);
 * \brief sets the order spacing polynomial
 */
void operaSpectralOrderVector::setOrderSpacingPolynomial(PolynomialCoeffs_t *pc) {
	delete orderSpacingPolynomial;
	orderSpacingPolynomial = new Polynomial(pc);
}


/*
 * \sa method Polynomial *getDispersionPolynomial(unsigned index);
 * \brief gets the dispersion polynomial
 */
LaurentPolynomial *operaSpectralOrderVector::getDispersionPolynomial(unsigned index) {
    if (index >= numberOfDispersionPolynomials) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	return dispersionPolynomial[index];
}

/*
 * setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc);
 * \brief sets the dispersion polynomial
 */
void operaSpectralOrderVector::setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc) {
    if (index >= MAXORDEROFWAVELENGTHPOLYNOMIAL) {
		throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	delete dispersionPolynomial[index];
	dispersionPolynomial[index] = new LaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial,pc);
}

/* 
 * void ReadSpectralOrders(string Filename);
 * \brief augment an existing vector with information from a file
 */
void operaSpectralOrderVector::ReadSpectralOrders(string Filename) {
	operaSpectralOrder_t format = None;
	operaistream fin(Filename.c_str());
	if (fin.is_open()) {
		string dataline;
		if (fin.good()) {
			getline(fin, dataline);
			if (!dataline.compare("#!geom")) {
				format = Geom;
			}
			if (!dataline.compare("#!wave")) {
				format = Wave;
			}
			if (!dataline.compare("#!prof")) {
				format = Prof;
			}
			if (!dataline.compare("#!aper")) {
				format = Aperture;
			}
			if (!dataline.compare("#!line")) {
				format = Lines;
			}
			if (!dataline.compare("#!fcal")) {
				format = Fcal;
			}            
			if (!dataline.compare("#!rvel")) {
				format = RVel;
			}            
			if (!dataline.compare("#!rawspectrum")) {
				format = RawSpectrum;
			}
			if (!dataline.compare("#!standardspectrum")) {
				format = StandardSpectrum;
			}
			if (!dataline.compare("#!optimalspectrum")) {
				format = OptimalSpectrum;
			}
			if (!dataline.compare("#!operaoptimalspectrum")) {
				format = OperaOptimalSpectrum;
			}
			if (!dataline.compare("#!rawbeamspectrum")) {
				format = RawBeamSpectrum;
			}
			if (!dataline.compare("#!standardbeamspectrum")) {
				format = StandardBeamSpectrum;
			}
			if (!dataline.compare("#!optimalbeamspectrum")) {
				format = OptimalBeamSpectrum;
			}
			if (!dataline.compare("#!operaoptimalbeamspectrum")) {
				format = OperaOptimalBeamSpectrum;
			}
			if (!dataline.compare("#!calibratedrawspectrum")) {
				format = CalibratedRawSpectrum;
			}
			if (!dataline.compare("#!calibratedstandardspectrum")) {
				format = CalibratedStandardSpectrum;
			}
			if (!dataline.compare("#!calibratedoptimalspectrum")) {
				format = CalibratedOptimalSpectrum;
			}
			if (!dataline.compare("#!calibratedoperaoptimalspectrum")) {
				format = CalibratedOperaOptimalSpectrum;
			}
			if (!dataline.compare("#!calibratedrawbeamspectrum")) {
				format = CalibratedRawBeamSpectrum;
			}
			if (!dataline.compare("#!calibratedstandardbeamspectrum")) {
				format = CalibratedStandardBeamSpectrum;
			}
			if (!dataline.compare("#!calibratedoptimalbeamspectrum")) {
				format = CalibratedOptimalBeamSpectrum;
			}
			if (!dataline.compare("#!calibratedoperaoptimalbeamspectrum")) {
				format = CalibratedOperaOptimalBeamSpectrum;
			}
			if (!dataline.compare("#!calibratedextendedbeamspectrum")) {
				format = CalibratedExtendedBeamSpectrum;
			}
			if (!dataline.compare("#!extendedpolarimetry")) {
				format = ExtendedPolarimetry;
			}
			if (!dataline.compare("#!referencespectrum")) {
				format = ReferenceSpectrum;
			}
			if (!dataline.compare("#!ordp")) {
				format = Orderspacing;
			}
			if (!dataline.compare("#!disp")) {
				format = Disp;
			}
			if (!dataline.compare("#!polar")) {
				format = Polarimetry;
			}
			if (!dataline.compare("#!gain")) {
				format = GainNoise;
			}
			if (!dataline.compare("#!SNR")) {
				format = SNR;
			}
			if (startsWith(dataline.c_str(), "#!csv")) {
				format = CSV;
			}
			if (startsWith(dataline.c_str(), "***")) {
				format = LibreEspritSpectrum;
			}
		}
	}
	fin.close();
	unsigned count2 = 0, minorder2 = 0, maxorder2 = 0;
	switch (format) {
		case GainNoise:
			readGainNoise(Filename);
			break;
		case Polarimetry:
			readOrdersFromPolar(Filename, count2, minorder2, maxorder2);
			break;
		case Aperture:
			readOrdersFromAperture(Filename, count2, minorder2, maxorder2);
			break;
		case SNR:
			readOrdersFromSNR(Filename, count2, minorder2, maxorder2);
			break;
		case Geom:
			readOrdersFromGeometry(Filename, count2, minorder2, maxorder2);
			break;
		case Wave:
			readOrdersFromWavelength(Filename, count2, minorder2, maxorder2);
			break;
		case Prof:
			readOrdersFromProfile(Filename, count2, minorder2, maxorder2);
			break;
		case Lines:
			readOrdersFromLines(Filename, count2, minorder2, maxorder2);
			break;
		case Fcal:
			readOrdersFromFluxCalibrationBeamSpectrum(Filename, format, count2, minorder2, maxorder2);
			break;            
		case RVel:
			readRadialVelocityCorrection(Filename);
			break;
		case RawSpectrum:
		case StandardSpectrum:
		case OptimalSpectrum:
		case OperaOptimalSpectrum:
			readOrdersFromSpectrum(Filename, format, count2, minorder2, maxorder2);
			break;
		case RawBeamSpectrum:
		case StandardBeamSpectrum:
		case OptimalBeamSpectrum:
		case OperaOptimalBeamSpectrum:
			readOrdersFromBeamSpectrum(Filename, format, count2, minorder2, maxorder2);
			break;
		case CalibratedRawSpectrum:
		case CalibratedStandardSpectrum:
		case CalibratedOptimalSpectrum:
		case CalibratedOperaOptimalSpectrum:
			readOrdersFromCalibratedSpectrum(Filename, format, count2, minorder2, maxorder2);
			break;
		case CalibratedRawBeamSpectrum:
		case CalibratedStandardBeamSpectrum:
		case CalibratedOptimalBeamSpectrum:
		case CalibratedOperaOptimalBeamSpectrum:
			readOrdersFromCalibratedBeamSpectrum(Filename, format, count2, minorder2, maxorder2);
			break;
		case CalibratedExtendedBeamSpectrum:
			readOrdersFromCalibratedExtendedBeamSpectrum(Filename, format, count2, minorder2, maxorder2);
			break;
		case ExtendedPolarimetry:
			readOrdersFromExtendedPolarimetry(Filename, count2, minorder2, maxorder2);
			break;
		case ReferenceSpectrum:
			readOrdersFromReferenceSpectrum(Filename, count2, minorder2, maxorder2);
			break;
		case LibreEspritSpectrum:
			readOrdersFromLibreEspritSpectrum(Filename, count2, minorder2, maxorder2);
			break;
		case Orderspacing:
			readOrderSpacingPolynomial(Filename);
			break;
		case Disp:
			readDispersionPolynomial(Filename);
			break;
		case CSV:
			readOrdersFromCSV(Filename);
			break;
		default:
			throw operaException("operaSpectralOrderVector: unkown content type in "+Filename+' ', operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);
			break;
	}
	if (format != GainNoise) {
#ifdef RANGE_CHECK
		if (maxorder2 > length) {
			throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (maxorder2 > maxorder) {
			maxorder = maxorder2;
		}
#ifdef RANGE_CHECK
		if (count2 > length) {
			throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
#ifdef RANGE_CHECK
		if (minorder2 > length) {
			throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
#endif
		if (minorder2 < minorder || minorder == 0) {
			minorder = minorder2;
		}
		count = maxorder - minorder + 1;
		first = vector + minorder;
		last = vector + maxorder;		
	}
}

/* 
 * operaSpectralOrder* operaSpectralOrderVector::GetSpectralOrder(unsigned order);
 * \brief Gets an operaSpectralOrder* to a given order, else NULL
 */
operaSpectralOrder* operaSpectralOrderVector::GetSpectralOrder(unsigned order) {
	operaSpectralOrder **vec = vector;
	operaSpectralOrder *spectralOrder = NULL;
	while ((spectralOrder = *vec++)) {
		if (spectralOrder->getorder() == order) {
			return spectralOrder;
		}
	}
	return NULL;
}

/* 
 * void operaSpectralOrderVector::ReadSpectralOrders(string Filename, operaSpectralOrder_t format)
 * \brief Reads from an m.fits product to create spectralorders. Also added CSV support.
 */
bool operaSpectralOrderVector::ReadSpectralOrders(string Filename, operaSpectralOrder_t format) {
	bool Gotit = false;
	if (!Filename.empty()) {
		if (format == CSV) {
			readOrdersFromCSV(Filename);
			return true;
		}
		if (Filename.find("m.fits") == string::npos) {
			throw operaException("operaSpectralOrderVector: "+Filename+": ", operaErrorInvalidInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		operaMEFFITSProduct in(Filename, READONLY);
		unsigned extensions = (unsigned)in.getNExtensions();
		switch (format) {
			case RawBeamSpectrum:
			case StandardBeamSpectrum:
			case OptimalBeamSpectrum:
			case OperaOptimalBeamSpectrum:
			case CalibratedRawBeamSpectrum:
			case CalibratedStandardBeamSpectrum:
			case CalibratedOptimalBeamSpectrum:
			case CalibratedOperaOptimalBeamSpectrum: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "BEAMFLUX") {
						//	<orders><order number><nElements><nBeams><elementindex><SpectralElements photoCenterX><SpectralElements photoCenterY><SpectralElements dist><SpectralElements wl><SpectralElements flux><SpectralElements flux variance><XCorrelation><nBeams>[repeat <beam><BeamElements[beam] photoCenterX><BeamElements[beam] photoCenterY><BeamElements[beam] flux><BeamElements[beam] flux variance>]
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							unsigned nElements = (unsigned)Product[row][column++];
							unsigned beams = (unsigned)Product[row][column++];
							column++; // skip index 
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							spectralOrder->createSpectralElements(nElements, format);
							operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
							spectralOrder->createBeamsAndBackgrounds(nElements, beams, format);
							spectralOrder->setnumberOfBeams(beams);
							spectralOrder->setSpectrumType(format);
							spectralOrder->sethasSpectralElements(true);
							spectralElements->setHasDistance(true);
							if (format == CalibratedRawBeamSpectrum ||
								format == CalibratedStandardBeamSpectrum ||
								format == CalibratedOptimalBeamSpectrum ||
								format == CalibratedOperaOptimalBeamSpectrum) {
								spectralElements->setHasWavelength(true);
							}
							spectralElements->setHasXCorrelation(true);
							for (unsigned element=0; element<nElements; element++) {
								double photoCenterX = Product[row][column++];
								double photoCenterY = Product[row][column++];
								double d = Product[row][column++];
								double wl = Product[row][column++];
								double flux = Product[row][column++];
								double variance = Product[row][column++];
								double xcorrelation = Product[row][column++];
								spectralElements->setdistd(d, element);
								spectralElements->setwavelength(wl, element);
								spectralElements->setFlux(flux, element);
								spectralElements->setFluxVariance(variance, element);
								spectralElements->setXCorrelation(xcorrelation, element);
								spectralElements->setphotoCenter(photoCenterX, photoCenterY, element);
								column++;	// skip nbeams
								// beams
								for (unsigned beam=0; beam < beams; beam++) {
									column++;	// skip beam
									double beamphotoCenterX = Product[row][column++];
									double beamphotoCenterY = Product[row][column++];
									double beamflux = Product[row][column++];
									double beamvariance = Product[row][column++];
									operaSpectralElements *beamElements = spectralOrder->getBeamElements(beam);
									beamElements->setFlux(beamflux, element);
									beamElements->setFluxVariance(beamvariance, element);
									beamElements->setphotoCenter(beamphotoCenterX, beamphotoCenterY, element);
								}
								column = 5;
								row++;
							}
							column = 2;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case GainNoise: {	// also bias
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "GAIN") {
						//<amps><amp><gain><noise><gainerror><bias><datasec x1><datasec x2><datasec y1><datasec y2><amp><gain><noise><gainerror><bias><datasec x1><datasec x2><datasec y1><datasec y2>
						unsigned row = 0;
						unsigned column = 0;
						unsigned amps = (unsigned)Product[row][column++];
						getGainBiasNoise()->setAmps(amps);
						for (unsigned amp=0; amp<amps; amp++) {
							column++; // skip amp
							float gain = Product[row][column++];
							float noise = Product[row][column++];
							float gainerror = Product[row][column++];
							float bias = Product[row][column++];
							DATASEC_t datasec;
							datasec.x1 = (unsigned)Product[row][column++];
							datasec.x2 = (unsigned)Product[row][column++];
							datasec.y1 = (unsigned)Product[row][column++];
							datasec.y2 = (unsigned)Product[row][column++];
							getGainBiasNoise()->setGain(amp, gain);
							getGainBiasNoise()->setGainError(amp, gainerror);
							getGainBiasNoise()->setNoise(amp, noise);
							getGainBiasNoise()->setBias(amp, bias);
							getGainBiasNoise()->setDatasec(amp, datasec);
							column = 1;
							row++;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Geom: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "GEOMETRY") {
						// <orders><order number><number of coefficients><ndatapoints>[<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL <chisqr><YBinning><miny><maxy>
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							unsigned npar = (unsigned)Product[row][column++];
							unsigned ndatapoints = (unsigned)Product[row][column++];
							float chisqr = 0.0;
							float miny = 0.0;
							float maxy = 0.0;
							unsigned ybinning = 0;
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							operaGeometry *geometry = spectralOrder->getGeometry();
							if (geometry == NULL) {
								spectralOrder->createGeometry(ndatapoints, ndatapoints);
							}
							geometry = spectralOrder->getGeometry();
							Polynomial *p = geometry->getCenterPolynomial();
							p->setOrderOfPolynomial(npar);
							float coeff = 0.0;
							float coefferr = 0.0;
							for (unsigned i=0; i<MAXPOLYNOMIAL; i++) {
								if (i < npar) {
									coeff = Product[row][column++];
									coefferr = Product[row][column++];
									p->setCoefficient(i, coeff);
									p->setCoefficientError(i, coefferr);									
								} else {
									column+=2;
								}
							}
							chisqr = Product[row][column++];
							ybinning = (unsigned)Product[row][column++];
							miny = Product[row][column++];
							maxy = Product[row][column++];
							p->setChisqr(chisqr);
							geometry->setYmin(miny);
							geometry->setYmax(maxy);
							geometry->setNumberofPointsToBinInYDirection(ybinning);
							spectralOrder->sethasGeometry(true);
							column = 2;
							row++;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Wave: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "WAVELENGTH") {
						// <orders><order number><number of coefficients><polynomial coefficient><polynomial coefficient error>...
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							unsigned npar = (unsigned)Product[row][column++];
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							if (!spectralOrder->getWavelength()) {
								spectralOrder->createWavelength(MAXORDEROFWAVELENGTHPOLYNOMIAL);
							}
							operaWavelength *wavelength = spectralOrder->getWavelength();         
							Polynomial *p = wavelength->getWavelengthPolynomial();
							p->setOrderOfPolynomial(npar);
							for (unsigned i=0; i<npar; i++) {
								float coeff = Product[row][column++];
								float coefferr = Product[row][column++];
								p->setCoefficient(i, coeff);
								p->setCoefficientError(i, coefferr);
							}
							spectralOrder->sethasWavelength(true);
							column = 2;
							row++;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Prof: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "PROFILE") {
						// <orders><order><number of columns i><number of rows j><xsize><xsampling><ysize><ysampling><i><j><number of coefficients><ndatapoints><polynomial coefficients><chisqr>
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							unsigned ii = (unsigned)Product[row][column++];
							unsigned jj = (unsigned)Product[row][column++];
							unsigned xsize = (unsigned)Product[row][column++];
							unsigned xsampling = (unsigned)Product[row][column++];
							unsigned ysize = (unsigned)Product[row][column++];
							unsigned ysampling = (unsigned)Product[row][column++];
							unsigned i = (unsigned)Product[row][column++]; i=i;
							unsigned j = (unsigned)Product[row][column++]; j=j;
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							spectralOrder->setInstrumentProfileVector(xsize, xsampling, ysize, ysampling, 1);
							operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
							spectralOrder->sethasInstrumentProfile(true);
							for (unsigned j=0; j<jj; j++) {
								for (unsigned i=0; i<ii; i++) {
									unsigned npar  = (unsigned)Product[row][column++];
									column++; // ndatapoints not used
									PolynomialCoeffs_t *pp = (PolynomialCoeffs_t *)malloc(sizeof(PolynomialCoeffs_t));
									if (!pp) {
										throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
									}
									pp->orderofPolynomial = npar;
									for	(unsigned coeff=0; coeff<MAXPOLYNOMIAL; coeff++) {
										if (coeff < npar) {
											pp->p[coeff] = Product[row][column++];
										} else {
											column++;
										}
									}
									instrumentProfile->setipPolyModelCoefficients(pp, i, j);
									float chisqr = Product[row][column++];
									instrumentProfile->setchisqrMatrixValue(chisqr, i, j);
									row++;
									column = 10;
								}
							}
							column = 2;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Disp: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "DISPERSION") {
						// <NumberOfDispersionPolynomials> <PolynomialIndex> <min coefficient> <max coefficient> [<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL
						unsigned row = 0;
						unsigned column = 0;
						unsigned NumberOfDispersionPolynomials = (unsigned)Product[row][column++];
						for (unsigned PolynomialIndex = 0; PolynomialIndex < NumberOfDispersionPolynomials; PolynomialIndex++) {
							unsigned PolyIndex = (unsigned)Product[row][column++];
							int MinorderOfLaurentPolynomial = (int)Product[row][column++];
							int MaxorderOfLaurentPolynomial = (int)Product[row][column++];
							setnumberOfDispersionPolynomials(NumberOfDispersionPolynomials);
							LaurentPolynomial *p = getDispersionPolynomial(PolyIndex);
							p->setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);
							p->setChisqr(0.0);
							for (unsigned i=0; i<p->getNumberOfCoefficients(); i++) {
								p->setCoefficient(i, Product[row][column++]);
								p->setCoefficientError(i, Product[row][column++]);
							}
							column = 1;
							row++;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case SNR: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "SNR") {
						// central SNR
						// <orders><Columns><order number><Center SNR>
						// All SNR
						// <orders><Columns><order number><nElements><Center SNR><wavelength><SNR>
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned columns = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
                        bool centralsnr = (columns <= 5);
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							spectralOrder->sethasSNR(true);
                            if (centralsnr) {
								float centersnr = Product[row][column++];
								spectralOrder->sethasCenterSNROnly(true);
								spectralOrder->setCenterSNR(centersnr);
								row++;
                            } else {
								unsigned nElements = (unsigned)Product[row][column++];
								spectralOrder->createSpectralElements(nElements, None);
                                operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                                spectralOrder->sethasSpectralElements(true);
                                spectralElements->setHasWavelength(true);
                                spectralElements->setHasFluxSNR(true);
                                for (unsigned element=0; element<nElements; element++) {
                                    float centersnr = Product[row][column++];
                                    float wl = Product[row][column++];
                                    float snr = Product[row][column++];
									spectralOrder->setCenterSNR(centersnr);
									spectralElements->setwavelength(wl, element);
                                    spectralElements->setFluxSNR(snr, element);
                                    row++;
                                    column = 4;
                                }
                            }
							column = 3;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Polarimetry: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "POLAR") {
						// <orders><order number><StokesParameter_t><method_t><length><distance><wavelength><Stokes(Q,U,V) flux><Stokes(Q,U,V) variance><StokesI flux><StokesI variance><degree of polarization flux><degree of polarization variance> <first null polarization><first null polarization variance><second null polarization><second null polarization variance>
						
						unsigned row = 0;
						unsigned column = 0;
						unsigned lastorder = 0;
						unsigned index = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							operaPolarimetry *Polarimetry = NULL;
							operaSpectralElements *spectralElements = NULL;
							unsigned StokesParameter = (unsigned)Product[row][column++];
							unsigned method = (unsigned)Product[row][column++];
							unsigned length = (unsigned)Product[row][column++];
							double distance = Product[row][column++];
							double wavelength = Product[row][column++];
							double crosscorrelation = Product[row][column++];
							double QUVFlux = Product[row][column++];
							double QUVVariance = Product[row][column++];
							double IFlux = Product[row][column++];
							double IVariance = Product[row][column++];
							double DegPolarFlux = Product[row][column++];
							double DegPolarVariance = Product[row][column++];
							double FirstNullPolarization = Product[row][column++];
							double FirstNullPolarizationVariance = Product[row][column++];
							double SecondNullPolarization = Product[row][column++];
							double SecondNullPolarizationVariance = Product[row][column++];
							
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							if (order != lastorder) {
								spectralOrder->createPolarimetry(length);
								Polarimetry = spectralOrder->getPolarimetry();
								spectralOrder->sethasPolarimetry(true);
								spectralOrder->createSpectralElements(length, None);
								spectralElements = spectralOrder->getSpectralElements();
								spectralOrder->sethasSpectralElements(true);
								spectralOrder->sethasWavelength(true);
								spectralElements->setHasWavelength(true);
								Polarimetry->setmethod((method_t)method);
								index = 0;
							}
							Polarimetry->setwavelength(wavelength, index);
							spectralElements->setwavelength(wavelength, index);
							spectralElements->setdistd(distance, index);
							spectralElements->setFlux(IFlux, index);
							spectralElements->setFluxVariance(IVariance, index);
							spectralElements->setXCorrelation(crosscorrelation, index);
							if ((stokes_parameter_t)StokesParameter == StokesQ) {
								spectralElements->setwavelength(wavelength, index);
								Polarimetry->setStokesParameter(StokesQ, QUVFlux, QUVVariance, index);
								Polarimetry->setStokesParameter(StokesI, IFlux, IVariance, index);
								Polarimetry->setDegreeOfPolarization(StokesQ, DegPolarFlux, DegPolarVariance, index);
								Polarimetry->setFirstNullPolarization(StokesQ, FirstNullPolarization, FirstNullPolarizationVariance, index);
								Polarimetry->setSecondNullPolarization(StokesQ, SecondNullPolarization, SecondNullPolarizationVariance, index);
								Polarimetry->setHasWavelength(true);
								Polarimetry->setHasStokesI(true);
								Polarimetry->setHasStokesQ(true);
							}
							if ((stokes_parameter_t)StokesParameter == StokesU) {
								spectralElements->setwavelength(wavelength, index);
								Polarimetry->setStokesParameter(StokesU, QUVFlux, QUVVariance, index);
								Polarimetry->setStokesParameter(StokesI, IFlux, IVariance, index);
								Polarimetry->setDegreeOfPolarization(StokesU, DegPolarFlux, DegPolarVariance, index);
								Polarimetry->setFirstNullPolarization(StokesU, FirstNullPolarization, FirstNullPolarizationVariance, index);
								Polarimetry->setSecondNullPolarization(StokesU, SecondNullPolarization, SecondNullPolarizationVariance, index);
								Polarimetry->setHasWavelength(true);
								Polarimetry->setHasStokesI(true);
								Polarimetry->setHasStokesU(true);
							}
							if ((stokes_parameter_t)StokesParameter == StokesV) {
								spectralElements->setwavelength(wavelength, index);
								Polarimetry->setStokesParameter(StokesV, QUVFlux, QUVVariance, index);
								Polarimetry->setStokesParameter(StokesI, IFlux, IVariance, index);
								Polarimetry->setDegreeOfPolarization(StokesV, DegPolarFlux, DegPolarVariance, index);
								Polarimetry->setFirstNullPolarization(StokesV, FirstNullPolarization, FirstNullPolarizationVariance, index);
								Polarimetry->setSecondNullPolarization(StokesV, SecondNullPolarization, SecondNullPolarizationVariance, index);
								Polarimetry->setHasWavelength(true);
								Polarimetry->setHasStokesI(true);
								Polarimetry->setHasStokesV(true);
							}
							column = 2;
							row++;
							index++;
							lastorder = order;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case ExtendedPolarimetry: 
				break;
			case Orderspacing: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "ORDERPOLY") {
						// <number of coefficients> [<polynomial coefficient><polynomial coefficienterror>]*MAXPOLYNOMIAL
						unsigned row = 0;
						unsigned column = 0;
						unsigned npar = (unsigned)Product[row][column++];
						Polynomial *p = getOrderSpacingPolynomial();
						p->setOrderOfPolynomial(npar);
						p->setChisqr(0.0);
						for (unsigned i=0; i<npar; i++) {
							p->setCoefficient(i, Product[row][column++]);
							p->setCoefficientError(i, Product[row][column++]);
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Aperture: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "APERTURE") {
						// <orders><order number><number of beams><measured tilt><tilt error><leftBackgroundIndex><xsampling><ysampling><lb height><lb width><lb slope><lb xcenter><lb ycenter><lb fluxFraction><rightBackgroundIndex><xsampling><ysampling><rb height><rb width><rb slope><rb xcenter><rb ycenter><rb fluxFraction><beam><xsampling><ysampling><beam height><beam width><beam slope><beam xcenter><beam ycenter><beam fluxFraction>
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						unsigned beams = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							float width, length, slope, midpointx, midpointy, fluxfraction;
							unsigned xsampling, ysampling;
							float tiltInDegreesValue = Product[row][column++];
							float tiltInDegreesError = Product[row][column++];
							
							operaSpectralOrder *spectralOrder = NULL;
							spectralOrder = GetSpectralOrder(order);
							spectralOrder->setnumberOfBeams(beams);
							spectralOrder->setTiltInDegrees(tiltInDegreesValue, tiltInDegreesError);
							spectralOrder->sethasExtractionApertures(true);
							
							for (unsigned background=0; background < 2; background++) {
								column++; // get rid of background index
								xsampling = (unsigned)Product[row][column++];
								ysampling = (unsigned)Product[row][column++];
								width = Product[row][column++];
								length = Product[row][column++];
								slope = Product[row][column++];
								midpointx  = Product[row][column++];
								midpointy = Product[row][column++];
								fluxfraction = Product[row][column++];
								operaPoint point(midpointx, midpointy);
								Line backgroundLineAperture(slope, width, length, &point);                    
								operaExtractionAperture *backgroundAperture = new operaExtractionAperture(&backgroundLineAperture, xsampling, ysampling);
								spectralOrder->setBackgroundApertures(background, backgroundAperture);
								backgroundAperture->setFluxFraction(fluxfraction); 
							}
							
							// beams...
							
							for (unsigned beam=0; beam < beams; beam++) {
								column++; // get rid of beam
								xsampling = (unsigned)Product[row][column++];
								ysampling = (unsigned)Product[row][column++];
								width = Product[row][column++];
								length = Product[row][column++];
								slope = Product[row][column++];
								midpointx  = Product[row][column++];
								midpointy = Product[row][column++];
								fluxfraction = Product[row][column++];
								operaPoint point(midpointx, midpointy);
								Line extractionLineAperture(slope, width, length, &point);
								operaExtractionAperture *extractionAperture = new operaExtractionAperture(&extractionLineAperture, xsampling, ysampling);
								spectralOrder->setExtractionApertures(beam, extractionAperture);
								extractionAperture->setFluxFraction(fluxfraction);
							}
							column = 3;
							row++;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Fcal: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					unsigned sequence = 0;
					if (extname == "FLUXCALIBRATION" || extname == "FLUXCALIBRATION1" || extname == "FLUXCALIBRATION2" || extname == "FLUXCALIBRATION3" || extname == "FLUXCALIBRATION4") {
						if (extname == "FLUXCALIBRATION1") {
							sequence = 1;
						} else if (extname == "FLUXCALIBRATION2") {
							sequence = 2;
						} else if (extname == "FLUXCALIBRATION3") {
							sequence = 3;
						} else if (extname == "FLUXCALIBRATION4") {
							sequence = 4;
						}
						setSequence(sequence);
						// <orders><order number><nElements><nBeams><elementindex><wavelength><SpectralElements flux conversion><flux conversion variance><SpectralElements throughput><throughput variance><nElements><nBeams><beam><BeamSED[beam] flux conversion><flux conversion variance><BeamSED[beam] throughput><throughput variance>
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							unsigned nElements = (unsigned)Product[row][column++];
							unsigned beams = (unsigned)Product[row][column++];
                            double wavelengthForNormalization = (unsigned)Product[row][column++];
							
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							spectralOrder->setnumberOfBeams(beams);
							spectralOrder->createSpectralEnergyDistributionElements(nElements);
                            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                            spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                            operaSpectralElements *FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements();
							operaSpectralElements *InstrumentThroughput = spectralEnergyDistribution->getThroughputElements(); 
							spectralOrder->sethasSpectralEnergyDistribution(true);
							FluxCalibration->setHasWavelength(true);
							InstrumentThroughput->setHasWavelength(true);                        
							spectralEnergyDistribution->setHasFluxCalibration(true);    
							spectralEnergyDistribution->setHasInstrumentThroughput(true);                                                 
							for (unsigned indexElem=0; indexElem < nElements; indexElem++) {
								column = 5;
								double wl = Product[row][column++];
								double fluxcal = Product[row][column++];
								double fcalvariance = Product[row][column++];
								double throughput = Product[row][column++];
								double throughputvariance = Product[row][column++];

								FluxCalibration->setwavelength(wl, indexElem);
								FluxCalibration->setFlux(fluxcal, indexElem);
								FluxCalibration->setFluxVariance(fcalvariance, indexElem);
								FluxCalibration->setphotoCenter(0.0, 0.0, indexElem);
								
								InstrumentThroughput->setwavelength(wl, indexElem);
								InstrumentThroughput->setFlux(throughput, indexElem);
								InstrumentThroughput->setFluxVariance(throughputvariance, indexElem);
								InstrumentThroughput->setphotoCenter(0.0, 0.0, indexElem);
								
								// beams
								for (unsigned beam=0; beam < beams; beam++) {
									unsigned b = (unsigned)Product[row][column++]; b = b;
									double beamfluxcal = Product[row][column++];
									double beamfcalvariance = Product[row][column++];
									double beamthroughput = Product[row][column++];
									double beamthrouputvariance = Product[row][column++];
									
									operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
									
									operaSpectralElements *beamFluxcalibration = BeamSED->getFluxCalibrationElements();                    
									beamFluxcalibration->setwavelength(wl, indexElem);
									beamFluxcalibration->setFlux(beamfluxcal, indexElem);
									beamFluxcalibration->setFluxVariance(beamfcalvariance, indexElem);
									beamFluxcalibration->setphotoCenter(0.0, 0.0, indexElem);
									
									operaSpectralElements *beamThroughput = BeamSED->getThroughputElements();                    
									beamThroughput->setwavelength(wl, indexElem);
									beamThroughput->setFlux(beamthroughput, indexElem);
									beamThroughput->setFluxVariance(beamthrouputvariance, indexElem);
									beamThroughput->setphotoCenter(0.0, 0.0, indexElem);
								}
								row++;
							}
							column = 2;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case RVel: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "RADIALVELOCITY") {
						// <BarycentricRadialVelocityCorrection>
						unsigned row = 0;
						unsigned column = 0;
						float rvel = (unsigned)Product[row][column++];
                        setBarycentricRadialVelocityCorrection(rvel);
						Gotit = true;
						break;
					}
				}
			}
				break;
			case PRVel: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "PRADIALVELOCITY") {
						// <BarycentricRadialVelocityCorrection>
						unsigned row = 0;
						unsigned column = 0;
						float rvel = (unsigned)Product[row][column++];
                        setBarycentricRadialVelocityCorrection(rvel);
						Gotit = true;
						break;
					}
				}
			}
				break;
			case Tell: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "TELLCORR") {
						// <orders><order number><number of coefficients><polynomial coefficient><polynomial coefficient error>...
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							unsigned npar = (unsigned)Product[row][column++];
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							if (!spectralOrder->getWavelength()) {
								spectralOrder->createWavelength(MAXORDEROFWAVELENGTHPOLYNOMIAL);
							}
							operaWavelength *wavelength = spectralOrder->getWavelength();         
							Polynomial *p = wavelength->getWavelengthPolynomial();
							p->setOrderOfPolynomial(npar);
							for (unsigned coeff=0; coeff<npar; coeff++) {
								float coefficient = Product[row][column++];
								float coefferr = Product[row][column++];
								p->setCoefficient(coeff, coefficient);
								p->setCoefficientError(coeff, coefferr);
							}
							spectralOrder->sethasWavelength(true);
							column = 2;
							row++;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			case PTell: {
				for (unsigned extension=1; extension<=extensions; extension++) {
					in.readExtension(extension);
					operaFITSProduct Product(in);
					string extname = in.operaFITSGetHeaderValue("EXTNAME", extension);
					if (extname == "PTELLCORR") {
						// <orders><order number><number of coefficients><polynomial coefficient><polynomial coefficient error>...
						unsigned row = 0;
						unsigned column = 0;
						unsigned orders = (unsigned)Product[row][column++];
						unsigned firstorder = (unsigned)Product[row][column++];
						for (unsigned order=firstorder; order<orders+firstorder; order++) {
							unsigned npar = (unsigned)Product[row][column++];
							operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
							if (!spectralOrder->getWavelength()) {
								spectralOrder->createWavelength(MAXORDEROFWAVELENGTHPOLYNOMIAL);
							}
							operaWavelength *wavelength = spectralOrder->getWavelength();         
							Polynomial *p = wavelength->getWavelengthPolynomial();
							p->setOrderOfPolynomial(npar);
							for (unsigned coeff=0; coeff<npar; coeff++) {
								float coefficient = Product[row][column++];
								float coefferr = Product[row][column++];
								p->setCoefficient(coeff, coefficient);
								p->setCoefficientError(coeff, coefferr);
							}
							spectralOrder->sethasWavelength(true);
							column = 2;
							row++;
						}
						Gotit = true;
						break;
					}
				}
			}
				break;
			default:
				break;
		}
	}
	return Gotit;
}
/* 
* void WriteSpectralOrders(string Filename, operaSpectralOrder_t Format), unsigned order=0, unsigned min=0;
* \brief Writes a SpectralOrder to a File
* \brief in the right place such that the vector remains ordered.
* \brief the optional order argument permits incremental addition to the output file, where zero means write all.
*/
void operaSpectralOrderVector::WriteSpectralOrders(string Filename, operaSpectralOrder_t Format, unsigned order, unsigned min) {
	operaostream fout;
	if (!Filename.empty()) {
		fout.open(Filename.c_str());
		switch (Format) {
			case GainNoise: {
				fout << "#!gain\n";
				fout << "######################################################################\n";
				fout << "# Gain Noise format is:\n";
				fout << "# <number of amps> <newline>\n";
				fout << "# <amp> <gain> <noise> <gainerror> <bias> <datasec x1> <datasec x2> <datasec y1> <datasec y2> <newline>\n";
				fout << "# ...\n";
				fout << "# Note that <amp> is zero-based.\n";
				fout << "# Note that gain is in units e/ADU, noise in e and bias in ADU.\n";
				fout << "#\n";
				fout << "######################################################################\n";
				fout << getGainBiasNoise()->getAmps() << endl;
				DATASEC_t datasec;
				for (unsigned i=0; i<getGainBiasNoise()->getAmps(); i++) {
					getGainBiasNoise()->getDatasec(i, datasec);
					fout << i << ' ' << getGainBiasNoise()->getGain(i) << ' '  << getGainBiasNoise()->getNoise(i) << ' ' << getGainBiasNoise()->getGainError(i)  << ' ' << getGainBiasNoise()->getBias(i) << ' ' << datasec.x1  << ' ' << datasec.x2  << ' ' << datasec.y1  << ' ' << datasec.y2 << endl;
				}
			}
				break;
			case Aperture: {
				fout << "#!aper\n";
				fout << "######################################################################\n";
				fout << "# Extraction Aperture format is:\n";
				fout << "#\n";
				fout << "# <number of orders><newline>\n";
				fout << "# <order number> <number of beams> <measured tilt> <tilt error>\n";
				fout << "# <leftBackgroundIndex> <xsampling> <ysampling> <lb height> <lb width> <lb slope> <lb xcenter> <lb ycenter> <lb fluxFraction>\n";
				fout << "# <rightBackgroundIndex> <xsampling> <ysampling>  <rb height> <rb width> <rb slope> <rb xcenter> <rb ycenter> <rb fluxFraction>\n";
				fout << "# <beam> <xsampling> <ysampling>  <beam height> <beam width> <beam slope> <beam xcenter> <beam ycenter>  <beam fluxFraction> <newline>\n";
				fout << "#\n";
				fout << "# Note that leftBackgroundIndex = 0\n";                
				fout << "# Note that rightBackgroundIndex = 1\n";                
				fout << "# Note that <beam> is zero-based.\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasExtractionApertures()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1) << ' ' << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
                    if (spectralOrder->gethasExtractionApertures()) {
                        fout << spectralOrder->getorder() << ' ' 
                        << spectralOrder->getnumberOfBeams() << ' '
                        << spectralOrder->getTiltInDegreesValue() << ' '
                        << spectralOrder->getTiltInDegreesError() << ' '; 
                        
                        for(unsigned backgroundIndex=0;backgroundIndex<LEFTANDRIGHT;backgroundIndex++) {
                            operaExtractionAperture *backgroundAperture = spectralOrder->getBackgroundApertures(backgroundIndex);
                            Line *backgroundLineAperture = backgroundAperture->getLineAperture();
                            
                            fout << backgroundIndex << ' '
                            << backgroundAperture->getXsampling() << ' '
                            << backgroundAperture->getYsampling() << ' '                            
                            << backgroundLineAperture->getWidth() << ' '
                            << backgroundLineAperture->getLength() << ' '
                            << backgroundLineAperture->getSlope() << ' '
                            << backgroundLineAperture->getMidPoint()->getXcoord() << ' '
                            << backgroundLineAperture->getMidPoint()->getYcoord() << ' '
                            << backgroundAperture->getFluxFraction() << ' ';                     
                        }
                        
                        for(unsigned beam=0;beam<spectralOrder->getnumberOfBeams(); beam++) {
                            operaExtractionAperture *beamAperture = spectralOrder->getExtractionApertures(beam);
                            Line *beamLineAperture = beamAperture->getLineAperture();
                            
                            fout << beam << ' '
                            << beamAperture->getXsampling() << ' '
                            << beamAperture->getYsampling() << ' '                            
                            << beamLineAperture->getWidth() << ' '
                            << beamLineAperture->getLength() << ' '
                            << beamLineAperture->getSlope() << ' '
                            << beamLineAperture->getMidPoint()->getXcoord() << ' '
                            << beamLineAperture->getMidPoint()->getYcoord() << ' '
                            << beamAperture->getFluxFraction() << ' ';                     
                        }     
						fout << endl;
                    }
				}
			}
				break;
			case Fcal: {
				fout << "#!fcal\n";
				fout << "######################################################################\n";
				fout << "# Flux Calibration Beam Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <wavelengthForNormalization> <elementindex> <wavelength> <SpectralElements flux conversion> <flux conversion variance> <SpectralElements throughput> <throughput variance> \n";
				fout << "# <beam> <BeamSED[beam] flux conversion> <flux conversion variance> <BeamSED[beam] throughput> <throughput variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralEnergyDistribution()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1) << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralEnergyDistribution()) {
						
						operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        
                        if (spectralEnergyDistribution->getHasFluxCalibration() && spectralEnergyDistribution->getHasInstrumentThroughput()){
                            
                            operaSpectralElements *FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements();                        
                            operaSpectralElements *InstrumentThroughput = spectralEnergyDistribution->getThroughputElements();
                            
                            for (unsigned indexElem=0; indexElem < FluxCalibration->getnSpectralElements(); indexElem++) {
                                
                                fout << spectralOrder->getorder() << ' ' << FluxCalibration->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
                                fout << fixed << setprecision(4) << spectralEnergyDistribution->getwavelengthForNormalization() << ' ';
                                fout << indexElem << ' ' << fixed << setprecision(4) << FluxCalibration->getwavelength(indexElem) << ' ';
                                fout << scientific << setprecision(6) << FluxCalibration->getFlux(indexElem) << ' '
                                << FluxCalibration->getFluxVariance(indexElem) << ' '
                                << InstrumentThroughput->getFlux(indexElem) << ' '                                
                                << InstrumentThroughput->getFluxVariance(indexElem) << ' ';
                                
                                for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                                    fout << beam << ' ' 
                                    << spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFlux(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFluxVariance(indexElem) << ' '             
                                    << spectralOrder->getBeamSED(beam)->getThroughputElements()->getFlux(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getThroughputElements()->getFluxVariance(indexElem) << ' ';                                                 
                                }
								fout << endl;
                            }  
                        }
					}						
				}
			}
				break;                
			case PRVel:
			case RVel: {
				fout << "#!rvel\n";
				fout << "######################################################################\n";
				fout << "# Radial Velocity Correction (km/s) format is:\n";
				fout << "# <radialvelocity> <newline>\n";
				fout << "#\n";
				fout << "######################################################################\n";
				fout << BarycentricRadialVelocityCorrection << endl;
			}
				break;
			case Polarimetry: {
				fout << "#!polar\n";
				fout << "######################################################################\n";
				fout << "# Polarimetry format is:\n";
				fout << "# <number of orders> <cols> <method> <newline>\n";
				fout << "# <order number> <StokesParameter_t> <length> <distance> <wavelength> <crosscorrelation> <Stokes(Q,U,V) flux> <Stokes(Q,U,V) variance> <StokesI flux> <StokesI variance> <degree of polarization flux> <degree of polarization variance> <first null polarization> <first null polarization variance> <second null polarization> <second null polarization variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				const unsigned columns = 14;
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasPolarimetry()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				bool firstline = true;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasPolarimetry() && spectralOrder->gethasSpectralElements()) {
						operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
						operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
						if (SpectralElements->getHasWavelength()) {
							operaStokesVector *StokesVector = Polarimetry->getStokesVector();
							operaStokesVector *DegreeOfPolarization = Polarimetry->getDegreeOfPolarization();
							operaStokesVector *FirstNullPolarization = Polarimetry->getFirstNullPolarization();
							operaStokesVector *SecondNullPolarization = Polarimetry->getSecondNullPolarization();
							unsigned length = Polarimetry->getLength();
							if (firstline) {
								fout << (maxorder - minorder + 1) << ' ' << columns << ' ' << Polarimetry->getmethod() << endl;
								firstline = false;
							}
							if (Polarimetry->getHasStokesQ()) {
								for (unsigned index = 0 ; index < StokesVector->getLength() ; index++) {
									fout << spectralOrder->getorder() << ' ' << StokesQ  << ' ' << length << ' ' 
									<< SpectralElements->getdistd(index) << ' ' ;
									fout << fixed << setprecision(4) << SpectralElements->getwavelength(index) << ' ';
									fout << scientific << SpectralElements->getXCorrelation(index) << ' '
									<< StokesVector->getStokesParameterFlux(StokesQ, index) << ' '
									<< StokesVector->getStokesParameterVariance(StokesQ, index) << ' '
									<< StokesVector->getStokesParameterFlux(StokesI, index) << ' '
									<< StokesVector->getStokesParameterVariance(StokesI, index) << ' '
									<< DegreeOfPolarization->getStokesParameterFlux(StokesQ, index) << ' '
									<< DegreeOfPolarization->getStokesParameterVariance(StokesQ, index) << ' ';
									if (Polarimetry->getHasFirstNullPolarization()) {
										fout  
										<< FirstNullPolarization->getStokesParameterFlux(StokesQ, index) << ' '
										<< FirstNullPolarization->getStokesParameterVariance(StokesQ, index) << ' ';
									} else {
										fout << " 0.0 0.0 ";
									}
									if (Polarimetry->getHasSecondNullPolarization()) {
										fout  
										<< SecondNullPolarization->getStokesParameterFlux(StokesQ, index) << ' '
										<< SecondNullPolarization->getStokesParameterVariance(StokesQ, index) << ' ';
										
									} else {
										fout << " 0.0 0.0 ";
									}
									fout << endl;
								}
							}
							if (Polarimetry->getHasStokesU()) {
								for (unsigned index = 0 ; index < StokesVector->getLength() ; index++) {
									fout << spectralOrder->getorder() << ' ' << StokesU << ' ' << length << ' '
									<< SpectralElements->getdistd(index) << ' ' ;
									fout << fixed << setprecision(4) << SpectralElements->getwavelength(index) << ' ';
									fout << scientific << SpectralElements->getXCorrelation(index) << ' '
									<< StokesVector->getStokesParameterFlux(StokesU, index) << ' '
									<< StokesVector->getStokesParameterVariance(StokesU, index) << ' '
									<< StokesVector->getStokesParameterFlux(StokesI, index) << ' '
									<< StokesVector->getStokesParameterVariance(StokesI, index) << ' '
									<< DegreeOfPolarization->getStokesParameterFlux(StokesU, index) << ' '
									<< DegreeOfPolarization->getStokesParameterVariance(StokesU, index) << ' ';
									if (Polarimetry->getHasFirstNullPolarization()) {
										fout 
										<< FirstNullPolarization->getStokesParameterFlux(StokesU, index) << ' '
										<< FirstNullPolarization->getStokesParameterVariance(StokesU, index) << ' ';
									} else {
										fout << " 0.0 0.0 ";
									}
									if (Polarimetry->getHasSecondNullPolarization()) {
										fout  
										<< SecondNullPolarization->getStokesParameterFlux(StokesU, index) << ' '
										<< SecondNullPolarization->getStokesParameterVariance(StokesU, index) << ' ';
										
									} else {
										fout << " 0.0 0.0 ";
									}
									fout << endl;
								}
							}
							if (Polarimetry->getHasStokesV()) {
								for (unsigned index = 0 ; index < StokesVector->getLength() ; index++) {
									fout << spectralOrder->getorder() << ' ' << StokesV << ' ' << length << ' ' 
									<< SpectralElements->getdistd(index) << ' ' ;
									fout << fixed << setprecision(4) << SpectralElements->getwavelength(index) << ' ';
									fout << scientific << SpectralElements->getXCorrelation(index) << ' '
									<< StokesVector->getStokesParameterFlux(StokesV, index) << ' '
									<< StokesVector->getStokesParameterVariance(StokesV, index) << ' '
									<< StokesVector->getStokesParameterFlux(StokesI, index) << ' '
									<< StokesVector->getStokesParameterVariance(StokesI, index) << ' '
									<< DegreeOfPolarization->getStokesParameterFlux(StokesV, index) << ' '
									<< DegreeOfPolarization->getStokesParameterVariance(StokesV, index) << ' ';
									if (Polarimetry->getHasFirstNullPolarization()) {
										fout 
										<< FirstNullPolarization->getStokesParameterFlux(StokesV, index) << ' '
										<< FirstNullPolarization->getStokesParameterVariance(StokesV, index) << ' ';
									} else {
										fout << " 0.0 0.0 ";
									}
									if (Polarimetry->getHasSecondNullPolarization()) {
										fout 
										<< SecondNullPolarization->getStokesParameterFlux(StokesV, index) << ' '
										<< SecondNullPolarization->getStokesParameterVariance(StokesV, index) << ' ';
									} else {
										fout << " 0.0 0.0 ";
									}
									fout << endl;
								}
							}
						}
					}
				}
			}
				break;
                
			case ExtendedPolarimetry: {
				fout << "#!extendedpolarimetry\n";
				fout << "######################################################################\n";
				fout << "# Extended Polarimetry format is:\n";
				fout << "# <number of orders> <StokesParameter_t> <method> <newline>\n";
				fout << "# <order number> <nElements> <elementindex> <wavelength> <wavelength telluric corrected> <barycentric wavelength correction> <crosscorrelation>\n";
				fout << "# <StokesI flux> <StokesI variance> <normalized StokesI flux> <calibrated StokesI flux>\n";
				fout << "# <degree of polarization> <degree of polarization variance> <first null polarization> <second null polarization> <newline>\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasPolarimetry()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				bool firstline = true;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasPolarimetry() && spectralOrder->gethasSpectralElements()) {
						operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
						operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
						if (SpectralElements->getHasWavelength()) {
							operaStokesVector *StokesVector = Polarimetry->getStokesVector();
							operaStokesVector *DegreeOfPolarization = Polarimetry->getDegreeOfPolarization();
							operaStokesVector *FirstNullPolarization = Polarimetry->getFirstNullPolarization();
							operaStokesVector *SecondNullPolarization = Polarimetry->getSecondNullPolarization();
							unsigned length = Polarimetry->getLength();

							if (Polarimetry->getHasStokesQ()) {
                                if (firstline) {
                                    fout << (maxorder - minorder + 1) << ' ' << StokesQ  << ' ' << Polarimetry->getmethod() << endl;
                                    firstline = false;
                                }
								for (unsigned index = 0 ; index < StokesVector->getLength() ; index++) {
									fout << spectralOrder->getorder() << ' ' << length << ' ' << index << ' ' ;
                                    fout << fixed << setprecision(8) << SpectralElements->getwavelength(index) << ' ';
                                    fout << SpectralElements->gettell(index) << ' ' << SpectralElements->getrvel(index) << ' ';
                                    fout << scientific << SpectralElements->getXCorrelation(index) << ' ';
                                    fout << SpectralElements->getFlux(index) << ' ' << SpectralElements->getFluxVariance(index) << ' ';
                                    fout << SpectralElements->getnormalizedFlux(index) << ' ' << SpectralElements->getfcalFlux(index) << ' ';
                                    fout << DegreeOfPolarization->getStokesParameterFlux(StokesQ, index) << ' ';
									fout << DegreeOfPolarization->getStokesParameterVariance(StokesQ, index) << ' ';
                                    
									if (Polarimetry->getHasFirstNullPolarization()) {
										fout  << FirstNullPolarization->getStokesParameterFlux(StokesQ, index) << ' ';
									} else {
										fout << " 0.0";
									}
									if (Polarimetry->getHasSecondNullPolarization()) {
										fout << SecondNullPolarization->getStokesParameterFlux(StokesQ, index) << ' ';
										
									} else {
										fout << " 0.0";
									}
									fout << endl;
								}
							}
							if (Polarimetry->getHasStokesU()) {
                                if (firstline) {
                                    fout << (maxorder - minorder + 1) << ' ' << StokesU  << ' ' << Polarimetry->getmethod() << endl;
                                    firstline = false;
                                }
								for (unsigned index = 0 ; index < StokesVector->getLength() ; index++) {
									fout << spectralOrder->getorder() << ' ' << length << ' ' << index << ' ' ;
                                    fout << fixed << setprecision(8) << SpectralElements->getwavelength(index) << ' ';
                                    fout << SpectralElements->gettell(index) << ' ' << SpectralElements->getrvel(index) << ' ';
                                    fout << scientific << SpectralElements->getXCorrelation(index) << ' ';
                                    fout << SpectralElements->getFlux(index) << ' ' << SpectralElements->getFluxVariance(index) << ' ';
                                    fout << SpectralElements->getnormalizedFlux(index) << ' ' << SpectralElements->getfcalFlux(index) << ' ';
                                    fout << DegreeOfPolarization->getStokesParameterFlux(StokesU, index) << ' ';
									fout << DegreeOfPolarization->getStokesParameterVariance(StokesU, index) << ' ';
                                    
									if (Polarimetry->getHasFirstNullPolarization()) {
										fout  << FirstNullPolarization->getStokesParameterFlux(StokesU, index) << ' ';
									} else {
										fout << " 0.0";
									}
									if (Polarimetry->getHasSecondNullPolarization()) {
										fout << SecondNullPolarization->getStokesParameterFlux(StokesU, index) << ' ';
										
									} else {
										fout << " 0.0";
									}
									fout << endl;
								}
							}
							if (Polarimetry->getHasStokesV()) {
                                if (firstline) {
                                    fout << (maxorder - minorder + 1) << ' ' << StokesV  << ' ' << Polarimetry->getmethod() << endl;
                                    firstline = false;
                                }
								for (unsigned index = 0 ; index < StokesVector->getLength() ; index++) {
									fout << spectralOrder->getorder() << ' ' << length << ' ' << index << ' ' ;
                                    fout << fixed << setprecision(8) << SpectralElements->getwavelength(index) << ' ';
                                    fout << SpectralElements->gettell(index) << ' ' << SpectralElements->getrvel(index) << ' ';
                                    fout << scientific << SpectralElements->getXCorrelation(index) << ' ';
                                    fout << SpectralElements->getFlux(index) << ' ' << SpectralElements->getFluxVariance(index) << ' ';
                                    fout << SpectralElements->getnormalizedFlux(index) << ' ' << SpectralElements->getfcalFlux(index) << ' ';
                                    fout << DegreeOfPolarization->getStokesParameterFlux(StokesV, index) << ' ';
									fout << DegreeOfPolarization->getStokesParameterVariance(StokesV, index) << ' ';
                                    
									if (Polarimetry->getHasFirstNullPolarization()) {
										fout  << FirstNullPolarization->getStokesParameterFlux(StokesV, index) << ' ';
									} else {
										fout << " 0.0";
									}
									if (Polarimetry->getHasSecondNullPolarization()) {
										fout << SecondNullPolarization->getStokesParameterFlux(StokesV, index) << ' ';
										
									} else {
										fout << " 0.0";
									}
									fout << endl;
								}
							}
						}
					}
				}
			}
				break;
			case SNR: {
				fout << "#!SNR\n";
				fout << "######################################################################\n";
				fout << "# SNR format is:\n";
				fout << "# <number of orders> <cols> <newline>\n";
				fout << "# <order number> <nElements> <Center SNR> <wavelength> <SNR> <newline>\n";
				fout << "# or\n";
				fout << "# <order number> <wavelength> <Center SNR> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder = *v;
				bool center = false;
				unsigned rows = 0;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasCenterSNROnly()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
						center = true;
						rows++;
					} else if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
						if (spectralOrder->gethasSpectralElements()) {
							if (spectralOrder->getorder() < minorder) {
								minorder = spectralOrder->getorder();
							}
							if (spectralOrder->getorder() > maxorder) {
								maxorder = spectralOrder->getorder();
							}
							rows += spectralOrder->getSpectralElements()->getnSpectralElements();
						}
					}
				}
				if (center) {
					fout << (maxorder - minorder + 1) << " 2" << endl;
				} else {
					fout << rows << " 5" << endl;
				}
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						if (center) {
							fout << spectralOrder->getorder() << ' ' << fixed << setprecision(4) << spectralElements->getwavelength(spectralElements->getnSpectralElements()/2) << ' ';
							fout << scientific << spectralOrder->getCenterSNR() << endl;
						} else {
							for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
								fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getCenterSNR() << ' ';
								fout << fixed << setprecision(4) << spectralElements->getwavelength(k) << ' ';
								fout << scientific << spectralElements->getFluxSNR(k) << endl;
							}
						}
					}
				}
			}
				break;
			case LibreEspritSNR: {
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				bool center = false;
				unsigned rows = 0;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasCenterSNROnly()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
						if (!center) {
							rows = 0;
						}
						center = true;
						rows++;
					} else if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
						if (spectralOrder->gethasSpectralElements()) {
							if (spectralOrder->getorder() < minorder) {
								minorder = spectralOrder->getorder();
							}
							if (spectralOrder->getorder() > maxorder) {
								maxorder = spectralOrder->getorder();
							}
							rows += spectralOrder->getSpectralElements()->getnSpectralElements();
						}
					}
				}
				fout << "***SNR of '" << object << "'" << endl;
				if (center) {
					fout << (maxorder - minorder + 1) << " 1" << endl;
				} else {
					fout << rows << " 1" << endl;
				}
				for (order=maxorder; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
					unsigned count = 0;
					if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength() && spectralOrder->getSpectralElements()->getnSpectralElements() > 0) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						if (center) {
							if (!isnan(spectralOrder->getCenterSNR())) {
								fout << fixed << setprecision(4) << spectralElements->getwavelength(spectralElements->getnSpectralElements()/2) << ' ';
								fout << scientific << spectralOrder->getCenterSNR() << endl;
							}
						} else {
							for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
								if (!isnan(spectralElements->getFluxSNR(k))) {
									fout << fixed << setprecision(4) << spectralElements->getwavelength(k) << ' ';
									fout << scientific << spectralElements->getFluxSNR(k) << ' ' << endl;
									count++;
								}
							}
							if (count > 0) {
								if (NEWLINES_BETWEEN_ORDERS) fout << endl;	// for plotting, split the orders
							}
						}
					}
				}
			}
				break;
			case RawSpectrum: {
				fout << "#!rawspectrum\n";
				fout << "######################################################################\n";
				fout << "# Raw Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <distance> <flux> <variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ' << spectralElements->getdistd(k) << ' ' << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case CalibratedRawSpectrum: {
				fout << "#!calibratedrawspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Raw Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <wavelength> <flux> <variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(k) <<  ' ';
							fout << scientific << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case StandardSpectrum: {
				fout << "#!standardspectrum\n";
				fout << "######################################################################\n";
				fout << "# Standard Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <distance> <flux> <variance> <integrated distance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ' << spectralElements->getdistd(k) << ' ' << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case CalibratedStandardSpectrum: {
				fout << "#!calibratedstandardspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Standard Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <wavelength> <flux> <variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(k) <<  ' ';
							fout << scientific << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case OptimalSpectrum: {
				fout << "#!optimalspectrum\n";
				fout << "######################################################################\n";
				fout << "# Optimal Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <distance> <flux> <variance> <integrated distance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ' << spectralElements->getdistd(k) << ' ' << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case CalibratedOptimalSpectrum: {
				fout << "#!calibratedoptimalspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Optimal Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <wavelength> <flux> <variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(k) <<  ' ';
							fout << scientific << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case OperaOptimalSpectrum: {
				fout << "#!operaoptimalspectrum\n";
				fout << "######################################################################\n";
				fout << "# Opera Optimal Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <distance> <flux> <variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(k) <<  ' ';
							fout << scientific << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case CalibratedOperaOptimalSpectrum: {
				fout << "#!calibratedoperaoptimalspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Opera Optimal Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <wavelength> <flux> <variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned k=0; k<spectralElements->getnSpectralElements(); k++) {
							fout << spectralOrder->getorder() << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(k) <<  ' ';
							fout << scientific << spectralElements->getFlux(k) << ' ' << spectralElements->getFluxVariance(k) << endl;
						}
					}						
				}
			}
				break;
			case RawBeamSpectrum: {
				fout << "#!rawbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Raw Beam Spectrum (as output from operaExtraction) format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <SpectralElements photoCenterX> <SpectralElements photoCenterY> <SpectralElements dist> <SpectralElements flux> <SpectralElements flux variance> <XCorrelation>\n";
				fout << "# <beam> <BeamElements[beam] photoCenterX> <BeamElements[beam] photoCenterY> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralelements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralelements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' '
							<< spectralelements->getphotoCenterX(indexElem) << ' ' 
							<< spectralelements->getphotoCenterY(indexElem) << ' ' 
							<< spectralelements->getdistd(indexElem) << ' '
							<< spectralelements->getFlux(indexElem) << ' '
							<< spectralelements->getFluxVariance(indexElem) << ' '
							<< spectralelements->getXCorrelation(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterX(indexElem) << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterY(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}						
				}
			}
				break;
			case CalibratedRawBeamSpectrum: {
				fout << "#!calibratedrawbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Raw Beam Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <wavelength> <SpectralElements flux> <SpectralElements flux variance> \n";
				fout << "# <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(indexElem) << ' ';
							fout << scientific << spectralElements->getFlux(indexElem) << ' '
							<< spectralElements->getFluxVariance(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}						
				}
			}
				break;
			case StandardBeamSpectrum: {
				fout << "#!standardbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Standard Beam Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <SpectralElements photoCenterX> <SpectralElements photoCenterY> <SpectralElements dist> <SpectralElements flux> <SpectralElements flux variance> <XCorrelation>\n";
				fout << "# <beam> <BeamElements[beam] photoCenterX> <BeamElements[beam] photoCenterY> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements() ) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' '
							<< spectralElements->getphotoCenterX(indexElem) << ' ' 
							<< spectralElements->getphotoCenterY(indexElem) << ' ' 
							<< spectralElements->getdistd(indexElem) << ' '
							<< spectralElements->getFlux(indexElem) << ' '
							<< spectralElements->getFluxVariance(indexElem) << ' '
							<< spectralElements->getXCorrelation(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterX(indexElem) << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterY(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}
				}
			}
				break;
			case CalibratedStandardBeamSpectrum: {
				fout << "#!calibratedstandardbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Standard Beam Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <wavelength> <SpectralElements flux> <SpectralElements flux variance> \n";
				fout << "# <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(indexElem) << ' ';
							fout << scientific << spectralElements->getFlux(indexElem) << ' '
							<< spectralElements->getFluxVariance(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}
				}
			}
				break;
			case OptimalBeamSpectrum: {
				fout << "#!optimalbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Optimal Beam Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <SpectralElements photoCenterX> <SpectralElements photoCenterY> <SpectralElements dist> <SpectralElements flux> <SpectralElements flux variance> <XCorrelation>\n";
				fout << "# <beam> <BeamElements[beam] photoCenterX> <BeamElements[beam] photoCenterY> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' '
							<< spectralElements->getphotoCenterX(indexElem) << ' ' 
							<< spectralElements->getphotoCenterY(indexElem) << ' ' 
							<< spectralElements->getdistd(indexElem) << ' '
							<< spectralElements->getFlux(indexElem) << ' '
							<< spectralElements->getFluxVariance(indexElem) << ' '
							<< spectralElements->getXCorrelation(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterX(indexElem) << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterY(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}
				}
			}
				break;
			case CalibratedOptimalBeamSpectrum: {
				fout << "#!calibratedoptimalbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Optimal Beam Spectrum format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <wavelength> <SpectralElements flux> <SpectralElements flux variance><cross correlation> \n";
				fout << "# <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' '
							<< fixed << setprecision(4) << spectralElements->getwavelength(indexElem) << ' '
							<< scientific << spectralElements->getFlux(indexElem) << ' '
							<< spectralElements->getFluxVariance(indexElem) << ' '
							<< spectralElements->getXCorrelation(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}
				}
			}
				break;
			case OperaOptimalBeamSpectrum: {
				fout << "#!operaoptimalbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Opera Optimal Beam Spectrum (as output from operaExtraction) format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <SpectralElements photoCenterX> <SpectralElements photoCenterY> <SpectralElements dist> <SpectralElements flux> <SpectralElements flux variance> <XCorrelation>\n";
				fout << "# <beam> <BeamElements[beam] photoCenterX> <BeamElements[beam] photoCenterY> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' '
							<< spectralElements->getphotoCenterX(indexElem) << ' ' 
							<< spectralElements->getphotoCenterY(indexElem) << ' ' 
							<< spectralElements->getdistd(indexElem) << ' '
							<< spectralElements->getFlux(indexElem) << ' '
							<< spectralElements->getFluxVariance(indexElem) << ' '
							<< spectralElements->getXCorrelation(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterX(indexElem) << ' ' 
								<< spectralOrder->getBeamElements(beam)->getphotoCenterY(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}
				}
			}
				break;
			case CalibratedOperaOptimalBeamSpectrum: {
				fout << "#!calibratedoperaoptimalbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Opera Optimal Beam Spectrum (as output from operaExtraction) format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <wavelength> <SpectralElements flux> <SpectralElements flux variance> \n";
				fout << "# <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' ';
							fout << fixed << setprecision(4) << spectralElements->getwavelength(indexElem) << ' ';
							fout << scientific << spectralElements->getFlux(indexElem) << ' ' 
							<< spectralElements->getFluxVariance(indexElem) << ' ';
							
							for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
								fout << beam << ' ' 
								<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
								<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
							}
							fout << endl;
						}  
					}
				}
			}
				break;
			case CalibratedExtendedBeamSpectrum: {
				fout << "#!calibratedextendedbeamspectrum\n";
				fout << "######################################################################\n";
				fout << "# Calibrated Extended Beam Spectrum (as output from operaExtraction) format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nElements> <nBeams> <elementindex> <wavelength> <wavelength telluric corrected> <barycentric wavelength correction> <crosscorrelation>\n";
				fout << "# <SpectralElements flux> <SpectralElements flux variance> <SpectralElements normalizedFlux> <SpectralElements fcalFlux>\n";
				fout << "# <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance> ... <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
						for (unsigned indexElem=0;indexElem < spectralElements->getnSpectralElements(); indexElem++) {
							
							fout << spectralOrder->getorder() << ' ' << spectralElements->getnSpectralElements() << ' ' << spectralOrder->getnumberOfBeams() << ' ';
							fout << indexElem << ' ';
							fout << fixed << setprecision(8) << spectralElements->getwavelength(indexElem) << ' ';
							fout << spectralElements->gettell(indexElem) << ' ' << spectralElements->getrvel(indexElem) << ' ';
                            fout << scientific << spectralElements->getXCorrelation(indexElem) << ' '
                            << spectralElements->getFlux(indexElem) << ' '
                            << spectralElements->getFluxVariance(indexElem) << ' ';
							fout << spectralElements->getnormalizedFlux(indexElem) << ' ' << spectralElements->getfcalFlux(indexElem) << ' ';
							
							if (spectralOrder->getnumberOfBeams() > 1) {
								for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
									fout << beam << ' ' 
									<< spectralOrder->getBeamElements(beam)->getFlux(indexElem) << ' '
									<< spectralOrder->getBeamElements(beam)->getFluxVariance(indexElem) << ' ';             
								}
							}
							fout << endl;
						}  
					}
				}
			}
				break;
			case LibreEspritpolarimetry:{
                unsigned rows = 0;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
                    if (spectralOrder->gethasPolarimetry() && spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
						if (SpectralElements->getHasWavelength()) {
							double  DegreeOfPolarization, DegreeOfPolarizationVariance;
							double Intensity, NullSpectrum1, NullSpectrum2;
							operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
							unsigned length = spectralOrder->getPolarimetry()->getLength();
							stokes_parameter_t stokesParameter = StokesI;
							for (unsigned index = 0 ; index < length ; index++) {
								Intensity = SpectralElements->getFlux(index);
								//Variance = SpectralElements->getFluxVariance(index);
								//Polarization = Polarimetry->getStokesParameter(stokesParameter)->getflux(index);
								//PolarizationVariance = Polarimetry->getStokesParameter(stokesParameter)->getvariance(index);
								DegreeOfPolarization = Polarimetry->getDegreeOfPolarization(stokesParameter)->getflux(index);
								DegreeOfPolarizationVariance = Polarimetry->getDegreeOfPolarization(stokesParameter)->getvariance(index);
								NullSpectrum1 = Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(index);
								//NullSpectrum1Variance = Polarimetry->getFirstNullPolarization(stokesParameter)->getvariance(index);
								NullSpectrum2 = Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(index);
								//NullSpectrum2Variance = Polarimetry->getSecondNullPolarization(stokesParameter)->getvariance(index);
								if (!isnan(DegreeOfPolarization) && !isnan(NullSpectrum1) && !isnan(NullSpectrum2) && !isnan(DegreeOfPolarizationVariance)) {
									rows++;
								}
							}
						}
                    }
				}
				fout << "***Reduced spectrum of '" << object << "'" << endl;
				fout << rows << " 5" << endl;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
						if (spectralOrder->gethasPolarimetry() && spectralOrder->gethasSpectralElements() && SpectralElements->getHasWavelength()) {
							double  DegreeOfPolarization, DegreeOfPolarizationVariance;
							double Intensity, NullSpectrum1, NullSpectrum2;
							operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
							unsigned length = spectralOrder->getPolarimetry()->getLength();
							stokes_parameter_t stokesParameter = StokesI;
							if (Polarimetry->getHasStokesV()) {
								stokesParameter = StokesV;
							} else if (Polarimetry->getHasStokesQ()) {
								stokesParameter = StokesQ;
							} else if (Polarimetry->getHasStokesU()) {
								stokesParameter = StokesU;
							} else if (Polarimetry->getHasStokesI()) {
								stokesParameter = StokesI;
							}
							unsigned count = 0;
							for (unsigned index = 0 ; index < length ; index++) {
								Intensity = SpectralElements->getFlux(index);
								//Variance = SpectralElements->getFluxVariance(index);
								//Polarization = Polarimetry->getStokesParameter(stokesParameter)->getflux(index);
								//PolarizationVariance = Polarimetry->getStokesParameter(stokesParameter)->getvariance(index);
								DegreeOfPolarization = Polarimetry->getDegreeOfPolarization(stokesParameter)->getflux(index);
								DegreeOfPolarizationVariance = Polarimetry->getDegreeOfPolarization(stokesParameter)->getvariance(index);
								NullSpectrum1 = Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(index);
								//NullSpectrum1Variance = Polarimetry->getFirstNullPolarization(stokesParameter)->getvariance(index);
								NullSpectrum2 = Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(index);
								//NullSpectrum2Variance = Polarimetry->getSecondNullPolarization(stokesParameter)->getvariance(index);
								if (!isnan(DegreeOfPolarization) && !isnan(NullSpectrum1) && !isnan(NullSpectrum2) && !isnan(DegreeOfPolarizationVariance)) {
									fout << fixed << setprecision(4) << SpectralElements->getwavelength(index) << ' ';
									fout << scientific << Intensity << ' '
									<< DegreeOfPolarization << ' '
									<< NullSpectrum1 << ' '
									<< NullSpectrum2 << ' '
									<< sqrt(DegreeOfPolarizationVariance) << endl;
									count++;
								}
							}
							if (count > 0) {
								if (NEWLINES_BETWEEN_ORDERS) fout << endl; // split the orders for plotting
							}
						}
					}
				}
			}
				break;
			case LibreEspritpolSpectrum:{
                unsigned rows = 0;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
                    if (spectralOrder->gethasSpectralElements()) {
                        operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						if (spectralelements->getHasWavelength()) {
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								if (!isnan(spectralelements->getFlux(i)) && !isnan(spectralelements->getFluxVariance(i))) {
									rows++;
								}
							}
						}
                    }
				}
				fout << "***Reduced spectrum of '" << object << "'" << endl;
				fout << rows << " 2" << endl;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						if (spectralelements->getHasWavelength() && spectralelements->getnSpectralElements() > 0) {
							unsigned count = 0;
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								if (!isnan(spectralelements->getFlux(i)) && !isnan(spectralelements->getFluxVariance(i))) {
									fout << fixed << setprecision(4) << spectralelements->getwavelength(i) << ' ';
									fout << scientific << spectralelements->getFlux(i) << ' '
									<< sqrt(spectralelements->getFluxVariance(i)) << endl;
									count++;
								}
							}
							if (count > 0) {
								if (NEWLINES_BETWEEN_ORDERS) fout << endl; // split the orders for plotting
							}
						}
					}						
				}
			}
				break;
			case LibreEspritsp1Spectrum:{
                unsigned rows = 0;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
                    if (spectralOrder->gethasSpectralElements()) {
                        operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						if (spectralelements->getHasWavelength()) {
							operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
							operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								if (!isnan(spectralelements->getFlux(i)) && !isnan(beamElements1->getFlux(i)) && !isnan(beamElements0->getFluxVariance(i) && !isnan(beamElements1->getFluxVariance(i)))) {
									rows++;
								}
							}
						}
                    }
				}
				fout << "***Reduced spectrum of '" << object << "'" << endl;
				fout << rows << " 6" << endl;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						if (spectralelements->getHasWavelength() && spectralelements->getnSpectralElements() > 0) {
							operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
							operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
							unsigned count = 0;
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								if (!isnan(spectralelements->getFlux(i)) && !isnan(beamElements1->getFlux(i)) && !isnan(beamElements0->getFluxVariance(i) && !isnan(beamElements1->getFluxVariance(i)))) {
									fout << fixed << setprecision(4) << spectralelements->getwavelength(i) << ' ';
									fout << scientific << spectralelements->getFlux(i) << ' '
									<< beamElements0->getFlux(i) << ' '
									<< beamElements1->getFlux(i) << ' '
									<< sqrt(beamElements0->getFluxVariance(i)+beamElements1->getFluxVariance(i)) << ' '
									<< sqrt(beamElements0->getFluxVariance(i)) << ' '
									<< sqrt(beamElements1->getFluxVariance(i)) << endl;
									count++;
								}
							}
							if (count > 0) {
								if (NEWLINES_BETWEEN_ORDERS) fout << endl; // split the orders for plotting
							}
						}
					}						
				}
			}
				break;
			case LibreEspritsp2Spectrum:{
                unsigned rows = 0;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
                    if (spectralOrder->gethasSpectralElements()) {
                        operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						if (spectralelements->getHasWavelength()) {
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								if (!isnan(spectralelements->getFlux(i)) && !isnan(spectralelements->getFluxVariance(i))) {
									rows++;
								}
							}
						}
                    }
				}
				fout << "***Reduced spectrum of '" << object << "'" << endl;
				fout << rows << " 2" << endl;
				for (order=MAXORDERS-1; order>=minorder; order--) {
					operaSpectralOrder *spectralOrder = vector[order];
					if (spectralOrder->gethasSpectralElements()) {
						operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
						if (spectralelements->getHasWavelength() && spectralelements->getnSpectralElements() > 0) {
							unsigned count = 0;
							for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
								if (!isnan(spectralelements->getFlux(i)) && !isnan(spectralelements->getFluxVariance(i))) {
									fout << fixed << setprecision(4) << spectralelements->getwavelength(i) << ' ';
									fout << scientific << spectralelements->getFlux(i) << ' '
									<< sqrt(spectralelements->getFluxVariance(i)) << endl;
									count++;
								}
							}
							if (count > 0) {
								if (NEWLINES_BETWEEN_ORDERS) fout << endl; // split the orders for plotting
							}
						}
					}						
				}
			}
				break;
			case Orderspacing: {
				fout << "#!ordp\n";
				fout << "######################################################################\n";
				fout << "# Order Spacing Polynomial, the format is:\n";
				fout << "# <number of coefficients> <polynomial coefficient> <polynomial coefficienterror> ... <newline>\n";
				fout << "#\n";
				fout << "######################################################################\n";
				Polynomial *polynomial = getOrderSpacingPolynomial();
				unsigned npar = polynomial->getOrderOfPolynomial();
				fout << npar << ' ';
				for (unsigned coeff=0; coeff<npar; coeff++) {
					fout << polynomial->getCoefficient(coeff) << ' ' << polynomial->getCoefficientError(coeff) << ' ';
				}
				fout << endl;
			}
				break;
			case Disp: {
				fout << "#!disp\n";
				fout << "######################################################################\n";
				fout << "# Dispersion Polynomial, the format is:\n";
				fout << "# <numberOfDispersionPolynomials> <newline>\n";                
				fout << "# <PolynomialIndex> <MinorderOfPolynomial> <MaxorderOfPolynomial> <polynomial coefficient> <polynomial coefficienterror> ... <newline>\n";
				fout << "#\n";
				fout << "######################################################################\n";
				fout << getnumberOfDispersionPolynomials() << endl;
                for (unsigned dispIndex=0; dispIndex<getnumberOfDispersionPolynomials(); dispIndex++) {
                    LaurentPolynomial *polynomial = getDispersionPolynomial(dispIndex);
                    int minordOfLaurent = polynomial->getMinorderOfLaurentPolynomial();
                    int maxordOfLaurent = polynomial->getMaxorderOfLaurentPolynomial();
                    unsigned npar = polynomial->getNumberOfCoefficients();
                    
                    fout << dispIndex << ' ' << minordOfLaurent << ' ' << maxordOfLaurent << ' ';
                    for (unsigned coeff=0; coeff<npar; coeff++) {
                        fout << polynomial->getCoefficient(coeff) << ' ' << polynomial->getCoefficientError(coeff) << ' ';
                    }
                    fout << endl;
                }
			}
				break;
			case Geom: {
				fout << "#!geom\n";
				fout << "######################################################################\n";
				fout << "# Geometry Calibration, the format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <number of coefficients> <ndatapoints> <polynomial coefficient> <polynomial coefficienterror>... <chisqr> <YBinning> <miny><maxy> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasGeometry()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasGeometry()) {
						Polynomial *polynomial = spectralOrder->getGeometry()->getCenterPolynomial();
						unsigned npar = polynomial->getOrderOfPolynomial();
						fout << spectralOrder->getorder() << ' ' << npar << ' ' << spectralOrder->getGeometry()->getNdatapoints() << ' ';
						for (unsigned coeff=0; coeff<npar; coeff++) {
							fout << polynomial->getCoefficient(coeff) << ' ' << polynomial->getCoefficientError(coeff) << ' ';
						}
						float miny = spectralOrder->getGeometry()->getYmin();
						float maxy = spectralOrder->getGeometry()->getYmax();
						fout << polynomial->getChisqr() << ' ' << spectralOrder->getGeometry()->getNumberofPointsToBinInYDirection() << ' ' << miny << ' ' << maxy << endl;
						
					}						
				}
			}
				break;
			case Wave: {
				fout << "#!wave\n";
				fout << "######################################################################\n";
				fout << "# Wavelength Calibration, the format is:\n";
				fout << "# <count of orders> <newline>\n";
				fout << "# <order number> <number of coefficients> <polynomial coefficients> <polynomial coefficient error> ... <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasWavelength()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) { 
					if (spectralOrder->gethasWavelength()) {
						Polynomial *polynomial = spectralOrder->getWavelength()->getWavelengthPolynomial() ;
						unsigned npar = polynomial->getOrderOfPolynomial();
						fout << spectralOrder->getorder() << ' ' << npar << ' ';
						for	(unsigned coeff=0; coeff<npar; coeff++) {
							fout << scientific 
							<< polynomial->getCoefficient(coeff) << ' '
                            << polynomial->getCoefficientError(coeff) << ' ';
						}
						fout << endl;						
					}
				}
			}
				break;
			case Prof: {
				fout << "#!prof\n";
				fout << "######################################################################\n";
				fout << "# Instrument Profile Calibration the format is:\n";
				fout << "# <number of orders> <number of columns i> <number of rows j> <xsize> <xsampling> <ysize> <ysampling> <newline>\n";
				fout << "# <order number> <i> <j> <number of coefficients> <ndatapoints> <polynomial coefficients> <chisqr> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasInstrumentProfile()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				v = vector;
				unsigned index = 0;
				while ((spectralOrder = *v++)) { 
					if (spectralOrder->gethasInstrumentProfile()) {
						operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
						if (index == 0) {
							fout << (maxorder - minorder + 1) << ' ' << instrumentProfile->getNXPoints() << ' ' << instrumentProfile->getNYPoints() << ' ' << instrumentProfile->getxsize() << ' ' << instrumentProfile->getXsampling() << ' ' << instrumentProfile->getysize() << ' ' << instrumentProfile->getYsampling() << endl;
						}
						for (unsigned j=0; j<instrumentProfile->getNYPoints(); j++) {
							for (unsigned i=0; i<instrumentProfile->getNXPoints(); i++) {
								PolynomialCoeffs_t *pp = instrumentProfile->getipPolyModelCoefficients(i,j);
								unsigned npar = pp->orderofPolynomial;
								fout << spectralOrder->getorder() << ' ' << i << ' ' << j << ' ' << npar << ' ' << '0'/*instrumentProfile->getnDataPoints()*/ << ' ';
								for	(unsigned coeff=0; coeff<npar; coeff++) {
									fout << pp->p[coeff] << ' ';
								}
								fout << instrumentProfile->getchisqrMatrixValue(i, j) << endl;
							}
						}	
						index++;
					}
				}
			}
				break;
			case Lines: {
				fout << "#!line\n";
				fout << "######################################################################\n";
				fout << "# Spectral Lines format is:\n";
				fout << "# <number of orders> <newline>\n";
				fout << "# <order number> <nfeatures> <feature number> <nlines in feature> <line number in feature> <center> <error> <sigma> <error> <amplitude> <error> <chisqr> <newline>\n";
				fout << "# ...\n";
				fout << "#\n";
				fout << "######################################################################\n";
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralLines()) {
						if (spectralOrder->getorder() < minorder || minorder == 0) {
							minorder = spectralOrder->getorder();
						}
						if (spectralOrder->getorder() > maxorder || maxorder == 0) {
							maxorder = spectralOrder->getorder();
						}
					}
				}
				fout << (maxorder - minorder + 1)  << endl;
				v = vector;
				while ((spectralOrder = *v++)) {
					if (spectralOrder->gethasSpectralLines()) {
						operaSpectralLines *spectralLines = spectralOrder->getSpectralLines();
						unsigned nFeatures = spectralLines->getNFeatures();
						for(unsigned featurenumber=0;featurenumber<nFeatures;featurenumber++) {
							operaSpectralFeature *spectralFeature = spectralLines->getSpectralFeature(featurenumber);
							double *center = spectralFeature->getGaussianFit()->getCenterVector();
							double *sigma = spectralFeature->getGaussianFit()->getSigmaVector();
							double *amplitude = spectralFeature->getGaussianFit()->getAmplitudeVector();        
							double *centerError = spectralFeature->getGaussianFit()->getCenterErrorVector();
							double *sigmaError = spectralFeature->getGaussianFit()->getSigmaErrorVector();
							double *amplitudeError = spectralFeature->getGaussianFit()->getAmplitudeErrorVector(); 
							
							for(unsigned line=0; line<spectralFeature->getnLines(); line++) {
								fout << spectralOrder->getorder() << 
								' ' << nFeatures <<             
								' ' << featurenumber <<
								' ' << spectralFeature->getnLines() <<
								' ' << line <<
								' ' << *center++ << 
								' ' << *centerError++ <<             
								' ' << *sigma++ <<
								' ' << *sigmaError++ <<            
								' ' << *amplitude++ << 
								' ' << *amplitudeError++ << 
								' ' << spectralFeature->getGaussianFit()->getGaussianChisqr() << endl;
							}        
						}    
					}						
				}
			}
				break;
			case CSV: {
				operaSpectralOrder **v = vector;
				operaSpectralOrder *spectralOrder;
				while ((spectralOrder = *v++)) {
					switch (instrumentmode) {
						case MODE_POLAR:
							if (spectralOrder->gethasPolarimetry()) {
								fout << "#!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance,DegreeOfPolarization,StokesParameter,NullSpectrum1,NullSpectrum2,PolarizationError" << endl;
							} else {
								fout << "#!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance" << endl;
							}
							for (unsigned order=minorder; order<=maxorder; order++) {
								if (spectralOrder->gethasSpectralElements()) {
									operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
									for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
										
										if (spectralOrder->gethasPolarimetry()) {
											double PolarizationVariance, DegreeOfPolarization, DegreeOfPolarizationVariance, NullSpectrum1Variance, NullSpectrum2Variance;
											double Intensity, Variance, Polarization, NullSpectrum1, NullSpectrum2;
											operaPolarimetry *Polarimetry = spectralOrder->getPolarimetry();
											stokes_parameter_t stokesParameter = StokesI;
											if (Polarimetry->getHasStokesV()) {
												stokesParameter = StokesV;
											} else if (Polarimetry->getHasStokesQ()) {
												stokesParameter = StokesQ;
											} else if (Polarimetry->getHasStokesU()) {
												stokesParameter = StokesU;
											} else if (Polarimetry->getHasStokesI()) {
												stokesParameter = StokesI;
											}
											Intensity = Polarimetry->getStokesParameter(StokesI)->getflux(i);
											Variance = Polarimetry->getStokesParameter(StokesI)->getvariance(i);
											Polarization = Polarimetry->getStokesParameter(stokesParameter)->getflux(i);
											PolarizationVariance = Polarimetry->getStokesParameter(stokesParameter)->getvariance(i);
											DegreeOfPolarization = Polarimetry->getDegreeOfPolarization(stokesParameter)->getflux(i);
											DegreeOfPolarizationVariance = Polarimetry->getDegreeOfPolarization(stokesParameter)->getvariance(i);
											NullSpectrum1 = Polarimetry->getFirstNullPolarization(stokesParameter)->getflux(i);
											NullSpectrum1Variance = Polarimetry->getFirstNullPolarization(stokesParameter)->getvariance(i);
											NullSpectrum2 = Polarimetry->getSecondNullPolarization(stokesParameter)->getflux(i);
											NullSpectrum2Variance = Polarimetry->getSecondNullPolarization(stokesParameter)->getvariance(i);
											fout << ',' 
                                            << ",,polar" << ','
                                            << spectralelements->getwavelength(i) << ','
                                            << spectralelements->getFlux(i) << ','
                                            << spectralelements->getFlux(i) << ',' // fixme
                                            << sqrt(spectralelements->getFluxVariance(i)) << ','
                                            << spectralelements->getFluxSNR(i) << ','
											<< Polarization << ','
											<< stokesParameter << ','
											<< NullSpectrum1 << ','
											<< NullSpectrum2 << ','
											<< sqrt(PolarizationVariance)
                                            << endl;
										} else {	// gethasPolarimetry()
                                            fout << ",,pol" << ','
                                            << spectralelements->getwavelength(i) << ','
                                            << spectralelements->getFlux(i) << ','
                                            << spectralelements->getFlux(i) << ',' // fixme
                                            << sqrt(spectralelements->getFluxVariance(i)) << ','
                                            << spectralelements->getFluxSNR(i)
                                            << endl;
                                        }
									}
								}
							}
							break;
						case MODE_STAR_ONLY:
							fout << "#!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR" << endl;
							for (unsigned order=minorder; order<=maxorder; order++) {
								if (spectralOrder->gethasSpectralElements()) {
									operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
									for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
										fout << ",,sp2" << ','
										<< spectralelements->getwavelength(i) << ','
										<< spectralelements->getFlux(i) << ','
										<< spectralelements->getFlux(i) << ',' // fixme
										<< sqrt(spectralelements->getFluxVariance(i)) << ','
										<< spectralelements->getFluxSNR(i)
										<< endl;
									}
								}
							}
							break;
						case MODE_STAR_PLUS_SKY:
							fout << "#!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance" << endl;
							for (unsigned order=minorder; order<=maxorder; order++) {
								operaSpectralElements *beamElements0 = spectralOrder->getBeamElements(0);
								operaSpectralElements *beamElements1 = spectralOrder->getBeamElements(1);
								if (spectralOrder->gethasSpectralElements()) {
									operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
									for (unsigned i=0; i<spectralelements->getnSpectralElements(); i++) {
										fout << ",,sp1" << ','
										<< spectralelements->getwavelength(i) << ','
										<< (beamElements0->getFlux(i) + beamElements1->getFlux(i)) << ','
										<< (beamElements0->getFlux(i) + beamElements1->getFlux(i)) << ','// fixme
										<< sqrt(beamElements0->getFluxVariance(i)+beamElements1->getFluxVariance(i)) << ','
										<< spectralelements->getFluxSNR(i) << ','
										<< beamElements0->getFlux(i) << ','
										<< beamElements1->getFlux(i) << ','
										<< sqrt(beamElements0->getFluxVariance(i)) << ','
										<< sqrt(beamElements1->getFluxVariance(i)) << endl;
									}
								}
							}
							break;
						default:
							break;
					}
				}
			}
				break;
			default:
				break;
		}
		fout.close();
	}
}

/* 
 * void readOrdersFromPolar(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .p.s file.
 */
void operaSpectralOrderVector::readOrdersFromPolar(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	// <number of orders> <cols>\n";
	// <order number> <StokesParameter_t> <length> <distance> <wavelength> <crosscorrelation> <Stokes(Q,U,V) flux> <Stokes(Q,U,V) variance> <StokesI flux> <StokesI variance> <degree of polarization flux> <degree of polarization variance> <first null polarization> <first null polarization variance> <second null polarization> <second null polarization variance> <newline>\n";
	operaistream fpolar(filename.c_str());
	if (fpolar.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned method = 0;
		unsigned columns = 0;
		unsigned length = 0;
		unsigned index = 0;
		unsigned StokesParameter = 0;
		
		// Note use of class "Double" here, to handle nans and infs
		Double wavelength = 0.0;
		Double distance = 0.0;
		Double crosscorrelation = 0.0;
		Double QUVFlux = 0.0;
		Double QUVVariance = 0.0;
		Double IFlux = 0.0;
		Double IVariance = 0.0;
		Double DegPolarFlux = 0.0;
		Double DegPolarVariance = 0.0;
		Double FirstNullPolarization = 0.0;
		Double FirstNullPolarizationVariance = 0.0;
		Double SecondNullPolarization = 0.0;
		Double SecondNullPolarizationVariance = 0.0;
		
		operaPolarimetry *Polarimetry = NULL;
		operaSpectralElements *spectralElements = NULL;
		
		while (fpolar.good()) {
			getline(fpolar, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					sscanf(dataline.c_str(), "%u %u %u", &count, &columns, &method);
					line++;
				} else {
					//sscanf(dataline.c_str(), "%u %u %u %u %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &order, &StokesParameter, &length, &distance, &wavelength, &crosscorrelation, &QUVFlux, &QUVVariance, &IFlux, &IVariance, &DegPolarFlux, &DegPolarVariance, &FirstNullPolarization, &FirstNullPolarizationVariance, &SecondNullPolarization, &SecondNullPolarizationVariance);                    
					stringstream ss (stringstream::in | stringstream::out);
					ss << dataline.c_str();
					ss >> order;
					ss >> StokesParameter;
					ss >> length;
					ss >> distance;
					ss >> wavelength;
					ss >> crosscorrelation;
					ss >> QUVFlux;
					ss >> QUVVariance;
					ss >> IFlux;
					ss >> IVariance;
					ss >> DegPolarFlux;
					ss >> DegPolarVariance;
					ss >> FirstNullPolarization;
					ss >> FirstNullPolarizationVariance;
					ss >> SecondNullPolarization;
					ss >> SecondNullPolarizationVariance;
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order != lastorder) {
						operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
                        spectralOrder->createPolarimetry(length);
                        Polarimetry = spectralOrder->getPolarimetry();
						spectralOrder->sethasPolarimetry(true);
						Polarimetry->setmethod((method_t)method);
                        lastorder = order;
						index = 0;
						spectralOrder->createSpectralElements(length, None);
                        spectralElements = spectralOrder->getSpectralElements();
						spectralOrder->sethasSpectralElements(true);
						Polarimetry->setHasWavelength(true);
						spectralElements->setHasWavelength(true);
						spectralElements->setHasDistance(true);
						spectralElements->setHasXCorrelation(true);
					}
					// populate spectral elements with the combined Stokes Intensity flux
					//spectralElements->setnSpectralElements(index+1);
					Polarimetry->setwavelength(wavelength.d, index);
					spectralElements->setwavelength(wavelength.d, index);
					spectralElements->setdistd(distance.d, index);
					spectralElements->setFlux(IFlux.d, index);
					spectralElements->setFluxVariance(IVariance.d, index);
					spectralElements->setXCorrelation(crosscorrelation.d, index);
					if ((stokes_parameter_t)StokesParameter == StokesQ) {
						Polarimetry->setStokesParameter(StokesQ, QUVFlux.d, QUVVariance.d, index);
						Polarimetry->setStokesParameter(StokesI, IFlux.d, IVariance.d, index);
						Polarimetry->setDegreeOfPolarization(StokesQ, DegPolarFlux.d, DegPolarVariance.d, index);
						Polarimetry->setFirstNullPolarization(StokesQ, FirstNullPolarization.d, FirstNullPolarizationVariance.d, index);
						Polarimetry->setSecondNullPolarization(StokesQ, SecondNullPolarization.d, SecondNullPolarizationVariance.d, index);
						Polarimetry->setHasStokesI(true);
						Polarimetry->setHasStokesQ(true);
					}
					if ((stokes_parameter_t)StokesParameter == StokesU) {
						Polarimetry->setStokesParameter(StokesU, QUVFlux.d, QUVVariance.d, index);
						Polarimetry->setStokesParameter(StokesI, IFlux.d, IVariance.d, index);
						Polarimetry->setDegreeOfPolarization(StokesU, DegPolarFlux.d, DegPolarVariance.d, index);
						Polarimetry->setFirstNullPolarization(StokesU, FirstNullPolarization.d, FirstNullPolarizationVariance.d, index);
						Polarimetry->setSecondNullPolarization(StokesU, SecondNullPolarization.d, SecondNullPolarizationVariance.d, index);
						Polarimetry->setHasStokesI(true);
						Polarimetry->setHasStokesU(true);
					}
					if ((stokes_parameter_t)StokesParameter == StokesV) {
						Polarimetry->setStokesParameter(StokesV, QUVFlux.d, QUVVariance.d, index);
						Polarimetry->setStokesParameter(StokesI, IFlux.d, IVariance.d, index);
						Polarimetry->setDegreeOfPolarization(StokesV, DegPolarFlux.d, DegPolarVariance.d, index);
						Polarimetry->setFirstNullPolarization(StokesV, FirstNullPolarization.d, FirstNullPolarizationVariance.d, index);
						Polarimetry->setSecondNullPolarization(StokesV, SecondNullPolarization.d, SecondNullPolarizationVariance.d, index);
						Polarimetry->setHasStokesI(true);
						Polarimetry->setHasStokesV(true);
					}
					index++;
					line++;
				}
			}
		}
		fpolar.close();
	}
}
/* 
 * void readOrdersFromAperture(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .aper file.
 */
void operaSpectralOrderVector::readOrdersFromAperture(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
    operaistream faperture(filename.c_str());
	if (faperture.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned beams = 1;
		unsigned b = 0;
		Double tiltInDegreesValue = 0.0;
		Double tiltInDegreesError = 0.0;
		unsigned line = 0;
		operaSpectralOrder *spectralOrder = NULL;
		
		while (faperture.good()) {
			getline(faperture, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					ss >> count;
					line++;
				} else {
					Double width, length, slope, midpointx, midpointy, fluxfraction;
                    unsigned xsampling, ysampling;
					
					ss >> order;
					ss >> beams;
					ss >> tiltInDegreesValue;
					ss >> tiltInDegreesError;
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					spectralOrder = GetSpectralOrder(order);
					spectralOrder->setnumberOfBeams(beams);
					spectralOrder->setTiltInDegrees(tiltInDegreesValue.d, tiltInDegreesError.d);
					spectralOrder->sethasExtractionApertures(true);
					
					// backgrounds...
					
					for (unsigned background=0; background < 2; background++) {
						ss >> b;
						ss >> xsampling;
						ss >> ysampling;
						ss >> width;
						ss >> length;
						ss >> slope;
						ss >> midpointx;
						ss >> midpointy;
						ss >> fluxfraction;
						operaPoint point(midpointx.d, midpointy.d);
						Line backgroundLineAperture(slope.d, width.d, length.d, &point);                    
						operaExtractionAperture *backgroundAperture = new operaExtractionAperture(&backgroundLineAperture, xsampling, ysampling);
						spectralOrder->setBackgroundApertures(background, backgroundAperture);
						backgroundAperture->setFluxFraction(fluxfraction.d); 
					}
					
					// beams...
					
					for (unsigned beam=0; beam < beams; beam++) {
						ss >> b;
						ss >> xsampling;
						ss >> ysampling;
						ss >> width;
						ss >> length;
						ss >> slope;
						ss >> midpointx;
						ss >> midpointy;
						ss >> fluxfraction;
						operaPoint point(midpointx.d, midpointy.d);
						Line extractionLineAperture(slope.d, width.d, length.d, &point);
						operaExtractionAperture *extractionAperture = new operaExtractionAperture(&extractionLineAperture, xsampling, ysampling);
						spectralOrder->setExtractionApertures(beam, extractionAperture);
						extractionAperture->setFluxFraction(fluxfraction.d);
					}
					line++;
				}
			}
		}
		faperture.close();
	}
}

/* 
 * void readGainNoise(string filename)
 * \brief Reads gain, gain error and noise.
 */
void operaSpectralOrderVector::readGainNoise(string filename) {
	operaistream fgain(filename.c_str());
	if (fgain.is_open()) {
		string dataline;
		Double gain = 0.0;
		Double gainerror = 0.0;
		Double noise = 0.0;
		Double bias = 0.0;
		unsigned line = 0;
		unsigned count = 0;
		unsigned amp = 0;
		unsigned x1 = 0, y1 = 0, x2 = 0, y2 = 0;
		DATASEC_t datasec;
		
		while (fgain.good()) {
			getline(fgain, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					sscanf(dataline.c_str(), "%u", &count);
					getGainBiasNoise()->setAmps(count);
					line++;
				} else {
					stringstream ss (stringstream::in | stringstream::out);
					ss << dataline.c_str();
					//sscanf(dataline.c_str(), "%u %lf %lf %lf %lf %u %u %u %u", &amp, &gain, &noise, &gainerror, &bias, &x1, &x2, &y1, &y2);
					ss >> amp;
					ss >> gain;
					ss >> noise;
					ss >> gainerror;
					ss >> bias;
					ss >> x1;
					ss >> x2;
					ss >> y1;
					ss >> y2;
					datasec.x1 = x1;
					datasec.x2 = x2;
					datasec.y1 = y1;
					datasec.y2 = y2;
					getGainBiasNoise()->setGain(amp, gain.d);
					getGainBiasNoise()->setGainError(amp, gainerror.d);
					getGainBiasNoise()->setNoise(amp, noise.d);
					getGainBiasNoise()->setBias(amp, bias.d);
					getGainBiasNoise()->setDatasec(amp, datasec);
					line++;
				}
			}
		}
		fgain.close();
	}
}

/* 
 * void readOrdersFromSNR(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .geom file.
 */
void operaSpectralOrderVector::readOrdersFromSNR(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the SNR table
	//
	operaistream fsnr(filename.c_str());
	if (fsnr.is_open()) {
		Float snr = 0.0, wl = 0.0, centersnr = 0.0;
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned cols = 0;
		unsigned rows = 0;
		unsigned count = 0;
		operaSpectralElements *spectralElements = NULL;
		// <number of orders> <cols>\n";
		// <order number> <rows> <Center SNR> <wavelength> <SNR>\n";
		// or
		// <order number> <wavelength> <Center SNR>\n";
		
		while (fsnr.good()) {
			getline(fsnr, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					sscanf(dataline.c_str(), "%u %u", &rows, &cols);
					line++;
				} else {
					stringstream ss (stringstream::in | stringstream::out);
					ss << dataline.c_str();
					if (cols < 4) {
						//sscanf(dataline.c_str(), "%u %f %f", &order, &wl, &centersnr);
						ss >> order >> wl >> centersnr;
					} else {
						//sscanf(dataline.c_str(), "%u %u %f %f %f", &order, &rows, &centersnr, &wl, &snr);
						ss >> order >> rows >> centersnr >> wl >> snr;
					}
					if (order < minorder || minorder == 0)
						minorder = order;
					if (order > maxorder || maxorder == 0)
						maxorder = order;
					if (lastorder != order) {
						count = 0;
					}
					operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
					spectralOrder->setCenterSNR(centersnr.f);
					spectralOrder->sethasSNR(true);
					if (cols < 4) {
						spectralOrder->sethasCenterSNROnly(true);
					} else {
						if (spectralOrder->gethasSpectralElements() == false) {
							spectralOrder->createSpectralElements(rows, None);
							spectralOrder->sethasSpectralElements(true);
							count = 0;
						}
						spectralElements = spectralOrder->getSpectralElements();
						if (count >= spectralElements->getnSpectralElements()) 
							continue;
						spectralElements->setwavelength(wl.f, count);
						spectralElements->setFluxSNR(snr.f, count);
						spectralElements->setHasWavelength(true);
						spectralElements->setHasFluxSNR(true);
					}
					count++;
					line++;
					lastorder = order;
				}
			}
		}
		count = line;
		fsnr.close();
	}
}

/* 
 * void readOrdersFromGeometry(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .geom file.
 */
void operaSpectralOrderVector::readOrdersFromGeometry(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the geometry polynomial coefficients as written by operaGeometry
	//
	operaistream fgeom(filename.c_str());
	if (fgeom.is_open()) {
		Float chisqr = 0.0;
		Float miny = 0.0;
		Float maxy = 0.0;
		unsigned ybinning = 0;
		string dataline;
		unsigned npar;
		int ndatapoints = 0;
		unsigned order = 0;
		unsigned line = 0;
		
		while (fgeom.good()) {
			getline(fgeom, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					ss >> count;
					line++;
				} else {
					ss >> order;
					ss >> npar;
					ss >> ndatapoints;
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
					operaGeometry *geometry = spectralOrder->getGeometry();
					if (geometry == NULL) {
						spectralOrder->createGeometry(ndatapoints, ndatapoints);
					}
					geometry = spectralOrder->getGeometry();
					Polynomial *p = geometry->getCenterPolynomial();
					p->setOrderOfPolynomial(npar);
					float coeff = 0.0;
					float coefferr = 0.0;
					for (unsigned i=0; i<npar; i++) {
						ss >> coeff;
						ss >> coefferr;
						p->setCoefficient(i, coeff);
						p->setCoefficientError(i, coefferr);
					}
					ss >> chisqr;
					ss >> ybinning;
					ss >> miny;
					ss >> maxy;
					p->setChisqr(chisqr.f);
					geometry->setYmin(miny.f);
					geometry->setYmax(maxy.f);
					geometry->setNumberofPointsToBinInYDirection(ybinning);
					geometry->CalculateAndSetOrderLength();
					spectralOrder->sethasGeometry(true);
					line++;
				}
			}
		}
		fgeom.close();
	}
}

/* 
 * void readOrdersFromSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .s file.
 * \brief The raw spectrum has distances rather than wavelength.
 */
void operaSpectralOrderVector::readOrdersFromSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from .s file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned ordercount = 0;
		Double d, flux, variance;
		operaSpectralOrder *spectralOrder = NULL;		
		operaSpectralElements *spectralElements = NULL;
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					ss >> count;
					line++;
				} else {
					ss >> order;
					ss >> d;
					ss >> flux;
					ss >> variance;
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					if (lastorder != order) {
						spectralOrder = GetSpectralOrder(order);
						
						if (!spectralOrder->getSpectralElements()) {
							spectralOrder->createSpectralElements(MAXSPECTRALELEMENTSPERORDER, format);
						} else {
							spectralOrder->getSpectralElements()->Resizevector(MAXSPECTRALELEMENTSPERORDER, format);
						}
                        spectralElements = spectralOrder->getSpectralElements();
                        
						spectralOrder->sethasSpectralElements(true);
						spectralElements->setHasDistance(true);
						ordercount = 0;
					}
					spectralElements->setnSpectralElements(ordercount+1);
					spectralElements->setdistd(d.d, ordercount);
					spectralElements->setFlux(flux.d, ordercount);
					spectralElements->setFluxVariance(variance.d, ordercount);
					ordercount++;
					lastorder = order;
					line++;
				}
			}
		}
		fspectrum.close();
	}
}

/* 
 * void readOrdersFromBeamSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .s file.
 * \brief The raw spectrum has distances rather than wavelength.
 */
void operaSpectralOrderVector::readOrdersFromBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from .e file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned count = 0;
		unsigned beams = 1;
		unsigned b = 0;
		Double d, photoCenterX, photoCenterY, flux, variance, xcorrelation;
		Double beamphotoCenterX, beamphotoCenterY, beamflux, beamvariance;
		unsigned elementindex = 0, nElements = 0;
		operaSpectralOrder *spectralOrder = NULL;		
		operaSpectralElements *spectralElements = NULL;
		operaSpectralElements *beamElements = NULL;
		
		// <order number> <nElements> <nBeams> <elementindex> <SpectralElements photoCenterX> <SpectralElements photoCenterX> <SpectralElements dist> <SpectralElements flux> <SpectralElements flux variance> <XCorrelation>
		// <beam> <BeamElements[beam] photoCenterX> <BeamElements[beam] photoCenterY> <BeamElements[beam] flux> <BeamElements[beam] flux variance>
		// ...
		
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					ss >> count;
					line++;
				} else {
					
					ss >> order;
					ss >> nElements;
					ss >> beams;
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					if (lastorder != order) {
						spectralOrder = GetSpectralOrder(order);
						
						spectralOrder->createSpectralElements(nElements, format);
						spectralElements = spectralOrder->getSpectralElements();
                        
                        spectralOrder->createBeamsAndBackgrounds(nElements, beams, format);
						spectralOrder->setnumberOfBeams(beams);
						spectralOrder->setSpectrumType(format);
						spectralOrder->sethasSpectralElements(true);
						spectralElements->setHasDistance(true);
						spectralElements->setHasXCorrelation(true);
					}
					ss >> elementindex;
					ss >> photoCenterX;
					ss >> photoCenterY;
					ss >> d;
					ss >> flux;
					ss >> variance;
					ss >> xcorrelation;
#ifdef RANGE_CHECK
					if (elementindex > nElements) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					spectralElements->setdistd(d.d, elementindex);
					spectralElements->setFlux(flux.d, elementindex);
					spectralElements->setFluxVariance(variance.d, elementindex);
					spectralElements->setXCorrelation(xcorrelation.d, elementindex);
					spectralElements->setphotoCenter(photoCenterX.d, photoCenterY.d, elementindex);
					
					// beams
					for (unsigned beam=0; beam < beams; beam++) {
						ss >> b;
						ss >> beamphotoCenterX;
						ss >> beamphotoCenterY;
						ss >> beamflux;
						ss >> beamvariance;
						beamElements = spectralOrder->getBeamElements(beam);
						beamElements->setFlux(beamflux.d, elementindex);
						beamElements->setFluxVariance(beamvariance.d, elementindex);
						beamElements->setphotoCenter(beamphotoCenterX.d, beamphotoCenterY.d, elementindex);
					}
					lastorder = order;
					line++;
				}
			}
		}
		fspectrum.close();
	}
}

/* 
 * void readOrdersFromCalibratedSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .s file.
 */
void operaSpectralOrderVector::readOrdersFromCalibratedSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from .s file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned ordercount = 0;
		Double w, flux, fluxvariance;
		operaSpectralOrder *spectralOrder = NULL;		
		operaSpectralElements *spectralElements = NULL;
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					sscanf(dataline.c_str(), "%u", &count);
					line++;
				} else {
					//sscanf(dataline.c_str(), "%u %lf %lf %lf", &order, &w, &flux, &fluxvariance);
					ss >> order >> w >> flux >> fluxvariance;
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					if (lastorder != order) {
						spectralOrder = GetSpectralOrder(order);
						spectralOrder->setSpectrumType(format);
						spectralOrder->sethasWavelength(true);
                        
						spectralOrder->createSpectralElements(MAXSPECTRALELEMENTSPERORDER, format);
						spectralElements = spectralOrder->getSpectralElements();
                        
						spectralOrder->sethasSpectralElements(true);
						spectralElements->setHasWavelength(true);
						ordercount = 0;
					}
					spectralElements->setnSpectralElements(ordercount+1);
					spectralElements->setwavelength(w.d, ordercount);
					spectralElements->setFlux(flux.d, ordercount);
					spectralElements->setFluxVariance(fluxvariance.d, ordercount);
					ordercount++;
					lastorder = order;
					line++;
				}
			}
		}
		fspectrum.close();
	}
}

/* 
 * void readOrdersFromCalibratedBeamSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .e file.
 */
void operaSpectralOrderVector::readOrdersFromCalibratedBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from .s file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned count = 0;
		unsigned beams = 1;
		unsigned b = 0;
		Double wl, flux, variance;
		Double beamflux, beamvariance;
		Double xcorrelation;
		unsigned elementindex = 0, nElements = 0;
		operaSpectralOrder *spectralOrder = NULL;		
		operaSpectralElements *spectralElements = NULL;
		operaSpectralElements *beamElements = NULL;
		
		// <order number> <nElements> <nBeams> <elementindex> <wavelength> <SpectralElements flux> <SpectralElements flux variance> <cross correlation>
		// <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance>
		// ...
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					ss >> count;
					line++;
				} else {
					ss >> order;
					ss >> nElements;
					ss >> beams;
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					if (lastorder != order) {
						spectralOrder = GetSpectralOrder(order);
						
						spectralOrder->createSpectralElements(nElements, format);
						spectralElements = spectralOrder->getSpectralElements();
                        
                        spectralOrder->createBeamsAndBackgrounds(nElements, beams, format);
						spectralOrder->setnumberOfBeams(beams);
						spectralOrder->setSpectrumType(format);
						spectralOrder->sethasSpectralElements(true);
						spectralElements->setHasWavelength(true);
						spectralElements->setHasXCorrelation(true);
					}
					ss >> elementindex;
					ss >> wl;
					ss >> flux;
					ss >> variance;
					ss >> xcorrelation;
#ifdef RANGE_CHECK
					if (elementindex > nElements) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					spectralElements->setwavelength(wl.d, elementindex);
					spectralElements->setFlux(flux.d, elementindex);
					spectralElements->setFluxVariance(variance.d, elementindex);
					spectralElements->setXCorrelation(xcorrelation.d, elementindex);
					
					// beams
					for (unsigned beam=0; beam < beams; beam++) {
						ss >> b;
						ss >> beamflux;
						ss >> beamvariance;
						beamElements = spectralOrder->getBeamElements(beam);
						beamElements->setFlux(beamflux.d, elementindex);
						beamElements->setFluxVariance(beamvariance.d, elementindex);
					}
					lastorder = order;
					line++;
				}
			}
		}
		fspectrum.close();
	}
}
/* 
 * void readOrdersFromCalibratedExtendedBeamSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the i.e file.
 */
void operaSpectralOrderVector::readOrdersFromCalibratedExtendedBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from i.e file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned count = 0;
		unsigned beams = 1;
		unsigned b = 0;
		Double wl, flux, variance;
		Double beamflux, beamvariance;
		Double xcorrelation;
		Double tell, rvel;
		Double normalizedFlux, fcalFlux;
		unsigned elementindex = 0, nElements = 0;
		operaSpectralOrder *spectralOrder = NULL;		
		operaSpectralElements *spectralElements = NULL;
		operaSpectralElements *beamElements = NULL;
		
		// <number of orders> <newline>\n";
		// <order number> <nElements> <nBeams> <elementindex> <wavelength> <wavelength telluric corrected> <barycentric wavelength correction> <crosscorrelation>
		// <SpectralElements flux> <SpectralElements flux variance> <SpectralElements normalizedFlux> <SpectralElements fcalFlux>
		// <beam> <BeamElements[beam] flux> <BeamElements[beam] flux variance> <newline>
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					ss >> count;
					line++;
				} else {
					ss >> order;
					ss >> nElements;
					ss >> beams;
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					if (lastorder != order) {
						spectralOrder = GetSpectralOrder(order);
						
						spectralOrder->createSpectralElements(nElements, format, true);
						spectralElements = spectralOrder->getSpectralElements();
                        
                        spectralOrder->createBeamsAndBackgrounds(nElements, beams, format);
						spectralOrder->setnumberOfBeams(beams);
						spectralOrder->setSpectrumType(format);
						spectralOrder->sethasSpectralElements(true);
						spectralElements->setHasWavelength(true);
						spectralElements->setHasXCorrelation(true);
						spectralElements->setHasExtendedBeamFlux(true);
					}
					ss >> elementindex;
					ss >> wl >> tell >> rvel;
					ss >> xcorrelation;
					ss >> flux;
					ss >> variance;
					ss >> normalizedFlux >> fcalFlux;
#ifdef RANGE_CHECK
					if (elementindex > nElements) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					spectralElements->setwavelength(wl.d, elementindex);
					spectralElements->settell(tell.d, elementindex);
					spectralElements->setrvel(rvel.d, elementindex);
                    spectralElements->setXCorrelation(xcorrelation.d, elementindex);

					spectralElements->setFlux(flux.d, elementindex);
                    spectralElements->setFluxVariance(variance.d, elementindex);
                    
                    spectralElements->setrawFlux(flux.d,elementindex);
                    spectralElements->setrawFluxVariance(variance.d, elementindex);

					spectralElements->setnormalizedFlux(normalizedFlux.d, elementindex);
					spectralElements->setnormalizedFluxVariance(variance.d/(flux.d*flux.d), elementindex);

					spectralElements->setfcalFlux(fcalFlux.d, elementindex);
					spectralElements->setfcalFluxVariance(variance.d*(fcalFlux.d/flux.d)*(fcalFlux.d/flux.d), elementindex);
					
					// beams
					for (unsigned beam=0; beam < beams; beam++) {
						if (beams > 1) {
							ss >> b;
							ss >> beamflux;
							ss >> beamvariance;
						}
						beamElements = spectralOrder->getBeamElements(beam);
						beamElements->setFlux(flux.d, elementindex);
						beamElements->setFluxVariance(variance.d, elementindex);
					}
					lastorder = order;
					line++;
				}
			}
		}
		fspectrum.close();
	}
}

/* 
 * void readOrdersFromExtendedPolarimetry(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the i.e file.
 */
void operaSpectralOrderVector::readOrdersFromExtendedPolarimetry(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {

    // <number of orders> <StokesParameter_t> <method>
	// <order number> <nElements> <elementindex> <wavelength> <wavelength telluric corrected> <barycentric wavelength correction> <crosscorrelation>
    // <StokesI flux> <StokesI variance> <normalized StokesI flux> <calibrated StokesI flux>
    // <degree of polarization> <degree of polarization variance> <first null polarization> <second null polarization>

	operaistream fpolar(filename.c_str());
	if (fpolar.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned method = 0;
		unsigned length = 0;
		unsigned index = 0;
		unsigned StokesParameter = 0;
		
		// Note use of class "Double" here, to handle nans and infs
		Double wavelength = 0.0;
		Double crosscorrelation = 0.0;
		Double IFlux = 0.0;
		Double IVariance = 0.0;
		Double tell, rvel;
		Double normalizedFlux, fcalFlux;
		Double DegPolarFlux = 0.0;
		Double DegPolarVariance = 0.0;
		Double FirstNullPolarization = 0.0;
		Double SecondNullPolarization = 0.0;
		
		operaPolarimetry *Polarimetry = NULL;
		operaSpectralElements *spectralElements = NULL;
		
		while (fpolar.good()) {
			getline(fpolar, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					sscanf(dataline.c_str(), "%u %u %u", &count, &StokesParameter, &method);
					line++;
				} else {
					stringstream ss (stringstream::in | stringstream::out);
					ss << dataline.c_str();
					ss >> order;
					ss >> length;
					ss >> index;
					ss >> wavelength;
                    ss >> tell;
                    ss >> rvel;
					ss >> crosscorrelation;
					ss >> IFlux;
					ss >> IVariance;
					ss >> normalizedFlux;
					ss >> fcalFlux;
					ss >> DegPolarFlux;
					ss >> DegPolarVariance;
					ss >> FirstNullPolarization;
					ss >> SecondNullPolarization;
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order != lastorder) {
						operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
                        spectralOrder->createPolarimetry(length);
                        Polarimetry = spectralOrder->getPolarimetry();
						spectralOrder->sethasPolarimetry(true);
						Polarimetry->setmethod((method_t)method);
                        lastorder = order;
						index = 0;
						spectralOrder->createSpectralElements(length, None, true);
                        spectralElements = spectralOrder->getSpectralElements();
						spectralOrder->sethasSpectralElements(true);
						Polarimetry->setHasWavelength(true);
						spectralElements->setHasWavelength(true);
						spectralElements->setHasDistance(true);
						spectralElements->setHasXCorrelation(true);
					}
					// populate spectral elements with the combined Stokes Intensity flux
					//spectralElements->setnSpectralElements(index+1);
					Polarimetry->setwavelength(wavelength.d, index);
					spectralElements->setwavelength(wavelength.d, index);
					spectralElements->settell(tell.d, index);
					spectralElements->setrvel(rvel.d, index);

					spectralElements->setFlux(IFlux.d, index);
					spectralElements->setFluxVariance(IVariance.d, index);
                    spectralElements->setrawFlux(IFlux.d,index);
                    spectralElements->setrawFluxVariance(IVariance.d, index);

					spectralElements->setnormalizedFlux(normalizedFlux.d, index);
                    spectralElements->setnormalizedFluxVariance(IVariance.d/(IFlux.d*IFlux.d), index);
        
					spectralElements->setfcalFlux(fcalFlux.d, index);
					spectralElements->setfcalFluxVariance(IVariance.d*(fcalFlux.d/IFlux.d)*(fcalFlux.d/IFlux.d), index);
                    
					spectralElements->setXCorrelation(crosscorrelation.d, index);
					if ((stokes_parameter_t)StokesParameter == StokesQ) {
						Polarimetry->setStokesParameter(StokesQ, 0.0, 0.0, index);
						Polarimetry->setStokesParameter(StokesI, IFlux.d, IVariance.d, index);
						Polarimetry->setDegreeOfPolarization(StokesQ, DegPolarFlux.d, DegPolarVariance.d, index);
						Polarimetry->setFirstNullPolarization(StokesQ, FirstNullPolarization.d, 0.0, index);
						Polarimetry->setSecondNullPolarization(StokesQ, SecondNullPolarization.d, 0.0, index);
						Polarimetry->setHasStokesI(true);
						Polarimetry->setHasStokesQ(true);
					}
					if ((stokes_parameter_t)StokesParameter == StokesU) {
						Polarimetry->setStokesParameter(StokesU, 0.0, 0.0, index);
						Polarimetry->setStokesParameter(StokesI, IFlux.d, IVariance.d, index);
						Polarimetry->setDegreeOfPolarization(StokesU, DegPolarFlux.d, DegPolarVariance.d, index);
						Polarimetry->setFirstNullPolarization(StokesU, FirstNullPolarization.d, 0.0, index);
						Polarimetry->setSecondNullPolarization(StokesU, SecondNullPolarization.d, 0.0, index);
						Polarimetry->setHasStokesI(true);
						Polarimetry->setHasStokesU(true);
					}
					if ((stokes_parameter_t)StokesParameter == StokesV) {
						Polarimetry->setStokesParameter(StokesV, 0.0, 0.0, index);
						Polarimetry->setStokesParameter(StokesI, IFlux.d, IVariance.d, index);
						Polarimetry->setDegreeOfPolarization(StokesV, DegPolarFlux.d, DegPolarVariance.d, index);
						Polarimetry->setFirstNullPolarization(StokesV, FirstNullPolarization.d,0.0, index);
						Polarimetry->setSecondNullPolarization(StokesV, SecondNullPolarization.d, 0.0, index);
						Polarimetry->setHasStokesI(true);
						Polarimetry->setHasStokesV(true);
					}
					line++;
				}
			}
		}
		fpolar.close();
	}
}

/* 
 * void readOrdersFromFluxCalibrationBeamSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .s file.
 */
void operaSpectralOrderVector::readOrdersFromFluxCalibrationBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from .fcal file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = 0;
		unsigned lastorder = 0;
		unsigned line = 0;
		unsigned count = 0;
		unsigned beams = 1;
        double wavelengthForNormalization = 548;
		unsigned b = 0;
		Double wl, fluxcal, fcalvariance, throughput, throughputvariance;
		Double beamfluxcal, beamfcalvariance, beamthroughput, beamthrouputvariance;
		unsigned elementindex = 0, nElements = 0;
		operaSpectralOrder *spectralOrder = NULL;		
        
		operaSpectralEnergyDistribution *spectralEnergyDistribution = NULL;
		operaSpectralElements *FluxCalibration = NULL; 
        operaSpectralElements *InstrumentThroughput = NULL; 
        
		operaSpectralEnergyDistribution *BeamSED = NULL;        
		operaSpectralElements *beamFluxcalibration = NULL;
		operaSpectralElements *beamThroughput = NULL;
		
		// <order number> <nElements> <nBeams> <wavelengthForNormalization> <elementindex> <wavelength> <flux calibration> <flux calibration variance> <instrument throughput> <throughput variance>
		// <beam> <BeamElements[beam] flux calibration> <BeamElements[beam] flux cal variance> <BeamElements[beam] throughput> <BeamElements[beam] throughput variance>
		// ...
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {
					ss >> count;
					line++;
				} else {
					
					ss >> order;
					ss >> nElements;
					ss >> beams;
                    ss >> wavelengthForNormalization;
                    
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					if (lastorder != order) {
						spectralOrder = GetSpectralOrder(order);
                        spectralOrder->setnumberOfBeams(beams);
                        spectralOrder->createSpectralEnergyDistributionElements(nElements);
                        spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                        FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements(); 
                        InstrumentThroughput = spectralEnergyDistribution->getThroughputElements(); 
						spectralOrder->sethasSpectralEnergyDistribution(true);
						FluxCalibration->setHasWavelength(true);
                        InstrumentThroughput->setHasWavelength(true);                        
                        spectralEnergyDistribution->setHasFluxCalibration(true);    
                        spectralEnergyDistribution->setHasInstrumentThroughput(true);                                                 
					}
					ss >> elementindex;
					ss >> wl;
					ss >> fluxcal;
					ss >> fcalvariance;
					ss >> throughput;
					ss >> throughputvariance;
#ifdef RANGE_CHECK
					if (elementindex > nElements) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					FluxCalibration->setwavelength(wl.d, elementindex);
					FluxCalibration->setFlux(fluxcal.d, elementindex);
					FluxCalibration->setFluxVariance(fcalvariance.d, elementindex);
					FluxCalibration->setphotoCenter(0.0, 0.0, elementindex);
					
					InstrumentThroughput->setwavelength(wl.d, elementindex);
					InstrumentThroughput->setFlux(throughput.d, elementindex);
					InstrumentThroughput->setFluxVariance(throughputvariance.d, elementindex);
					InstrumentThroughput->setphotoCenter(0.0, 0.0, elementindex);
					
					// beams
					for (unsigned beam=0; beam < beams; beam++) {
						ss >> b;
						ss >> beamfluxcal;
						ss >> beamfcalvariance;
						ss >> beamthroughput;
						ss >> beamthrouputvariance;
						
						BeamSED = spectralOrder->getBeamSED(beam);
						
						beamFluxcalibration = BeamSED->getFluxCalibrationElements();                    
						beamFluxcalibration->setwavelength(wl.d, elementindex);
						beamFluxcalibration->setFlux(beamfluxcal.d, elementindex);
						beamFluxcalibration->setFluxVariance(beamfcalvariance.d, elementindex);
						beamFluxcalibration->setphotoCenter(0.0, 0.0, elementindex);
						
						beamThroughput = BeamSED->getThroughputElements();                    
						beamThroughput->setwavelength(wl.d, elementindex);
						beamThroughput->setFlux(beamthroughput.d, elementindex);
						beamThroughput->setFluxVariance(beamthrouputvariance.d, elementindex);
						beamThroughput->setphotoCenter(0.0, 0.0, elementindex);
					}
					lastorder = order;
					line++;
				}
			}
		}
		fspectrum.close();
	}
}



/* 
 * void readOrdersFromReferenceSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief reads a reference (solar or moon) spectrum.
 */
void operaSpectralOrderVector::readOrdersFromReferenceSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from reference .s file
	//
	// <Wavelength (nm)> <ETR (W*m-2*nm-1)> <Global Tilt (W*m-2*nm-1)> <Direct + Circumsolar (W*m-2*nm-1)>
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned line = 0;
		double wavelength, etr, globaltilt, circumsolar;
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf %lf %lf", &wavelength, &etr, &globaltilt, &circumsolar);
					line++;
				}
			}
		}
		fspectrum.close();
	}
}

/* 
 * void readOrdersFromLibreEspritSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder)
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .s file.
 */
void operaSpectralOrderVector::readOrdersFromLibreEspritSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the spectrum from .s file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = minorder;
		float lastwavelength = BIG;
		unsigned line = 0;
		unsigned ordercount = 0;
		double w, flux, fluxvariance;
		operaSpectralOrder *spectralOrder = NULL;		
		operaSpectralElements *spectralElements = NULL;
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '*') {
					// skip comments
				} else if (line == 0) {
					sscanf(dataline.c_str(), "%u", &count);
					line++;
				} else {
					sscanf(dataline.c_str(), "%lf %lf %lf", &w, &flux, &fluxvariance);
					if (lastwavelength > w) {
						spectralOrder = GetSpectralOrder(order);
						
						spectralOrder->createSpectralElements(MAXSPECTRALELEMENTSPERORDER, None);
						spectralElements = spectralOrder->getSpectralElements();
                        
						spectralElements->setHasWavelength(true);
						spectralOrder->sethasSpectralElements(true);
						order++;
						ordercount = 0;
					}
					spectralElements->setnSpectralElements(ordercount+1);
					spectralElements->setwavelength(w, ordercount);
					spectralElements->setFlux(flux, ordercount);
					spectralElements->setFluxVariance(fluxvariance, ordercount);
					lastwavelength = w;
					ordercount++;
					count = ordercount;
					line++;
				}
			}
		}
		fspectrum.close();
	}
}

/* 
 * \brief Creates a null terminated vector of pointers to orders as initialized from the .s file.
 */
void operaSpectralOrderVector::readOrdersFromLibreEspritPolarimetry(string filename, stokes_parameter_t StokesParameter, unsigned &count, unsigned &maxorder) {
	//
	// read in the spectrum from .s file
	//
	operaistream fspectrum(filename.c_str());
	if (fspectrum.is_open()) {
		string dataline;
		unsigned order = maxorder;
		float LastWavelength = BIG;
		unsigned line = 0;
		unsigned TotalNumberOfLines;
		unsigned ordercount = 0;
		double Wavelength, Intensity, Polarization, NullSpectrum1, NullSpectrum2, Variance;
		operaSpectralOrder *spectralOrder = NULL;		
		operaPolarimetry *Polarimetry = NULL;
		
		while (fspectrum.good()) {
			getline(fspectrum, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '*') {
					// skip comments
				} else if (line == 0) {
					sscanf(dataline.c_str(), "%u", &TotalNumberOfLines);
					line++;
				} else {
					sscanf(dataline.c_str(), "%lf %lf %lf %lf %lf %lf", &Wavelength, &Intensity, &Polarization, &NullSpectrum1, &NullSpectrum2, &Variance);
					if (LastWavelength > Wavelength) {
                        if (Polarimetry)
                            Polarimetry->setLength(ordercount);
						spectralOrder = GetSpectralOrder(order);
                        spectralOrder->createPolarimetry(MAXSPECTRALELEMENTSPERORDER);
                        Polarimetry = spectralOrder->getPolarimetry();
						spectralOrder->sethasPolarimetry(true);
                        
						order--;
						ordercount = 0;
					}
                    Polarimetry->setStokesParameter(StokesI, Intensity, Variance, ordercount);
                    Polarimetry->setStokesParameter(StokesParameter, -Polarization, -Variance, ordercount);
                    Polarimetry->setDegreeOfPolarization(StokesParameter, (-Polarization/Intensity), (-Variance/Intensity), ordercount);
                    Polarimetry->setFirstNullPolarization(StokesParameter, (-NullSpectrum1/Intensity), (-Variance/Intensity), ordercount);
                    Polarimetry->setSecondNullPolarization(StokesParameter, (-NullSpectrum2/Intensity), (-Variance/Intensity), ordercount);
                    Polarimetry->setHasStokesI(true);
                    if (StokesParameter == StokesQ) {
                        Polarimetry->setHasStokesQ(true);
                        Polarimetry->setHasDegreeOfStokesQ(true);
                    }
                    if (StokesParameter == StokesU) {
                        Polarimetry->setHasStokesU(true);
                        Polarimetry->setHasDegreeOfStokesU(true);
                    }
                    if (StokesParameter == StokesV) {
                        Polarimetry->setHasStokesV(true);
                        Polarimetry->setHasDegreeOfStokesV(true);
                    }
                    Polarimetry->setHasFirstNullPolarization(true);
                    Polarimetry->setHasSecondNullPolarization(true);
                    
					LastWavelength = Wavelength;
					ordercount++;
					line++;
				}
			}
		}
        if (Polarimetry)
            Polarimetry->setLength(ordercount);
        
        count = maxorder - order + 1;
        
		fspectrum.close();
        
        double PolarizationVariance, DegreeOfPolarization, DegreeOfPolarizationVariance, NullSpectrum1Variance, NullSpectrum2Variance;
        
        for (unsigned temporder = maxorder; temporder > order ; temporder--) {
            spectralOrder = GetSpectralOrder(temporder);
            unsigned length = spectralOrder->getPolarimetry()->getLength();
            
            operaPolarimetry TemporaryPolarimetry(length);
            
            for (unsigned index = 0 ; index < length ; index++) {
                Intensity = Polarimetry->getStokesParameter(StokesI)->getflux(index);
                Variance = Polarimetry->getStokesParameter(StokesI)->getvariance(index);
                Polarization = Polarimetry->getStokesParameter(StokesParameter)->getflux(index);
                PolarizationVariance = Polarimetry->getStokesParameter(StokesParameter)->getvariance(index);
                DegreeOfPolarization = Polarimetry->getDegreeOfPolarization(StokesParameter)->getflux(index);
                DegreeOfPolarizationVariance = Polarimetry->getDegreeOfPolarization(StokesParameter)->getvariance(index);
                NullSpectrum1 = Polarimetry->getFirstNullPolarization(StokesParameter)->getflux(index);
                NullSpectrum1Variance = Polarimetry->getFirstNullPolarization(StokesParameter)->getvariance(index);
                NullSpectrum2 = Polarimetry->getSecondNullPolarization(StokesParameter)->getflux(index);
                NullSpectrum2Variance = Polarimetry->getSecondNullPolarization(StokesParameter)->getvariance(index);
                
                TemporaryPolarimetry.setStokesParameter(StokesI, Intensity, Variance, index);
                TemporaryPolarimetry.setStokesParameter(StokesParameter, Polarization, PolarizationVariance, index);
                TemporaryPolarimetry.setDegreeOfPolarization(StokesParameter, DegreeOfPolarization, DegreeOfPolarizationVariance, index);
                TemporaryPolarimetry.setFirstNullPolarization(StokesParameter, NullSpectrum1, NullSpectrum1Variance, index);
                TemporaryPolarimetry.setSecondNullPolarization(StokesParameter, NullSpectrum2, NullSpectrum2Variance, index);
            }
            
            unsigned reversedIndex = length;
            for (unsigned index = 0 ; index < length ; index++) {
                Intensity = TemporaryPolarimetry.getStokesParameter(StokesI)->getflux(reversedIndex);
                Variance = TemporaryPolarimetry.getStokesParameter(StokesI)->getvariance(reversedIndex);
                Polarization = TemporaryPolarimetry.getStokesParameter(StokesParameter)->getflux(reversedIndex);
                PolarizationVariance = TemporaryPolarimetry.getStokesParameter(StokesParameter)->getvariance(reversedIndex);
                DegreeOfPolarization = TemporaryPolarimetry.getDegreeOfPolarization(StokesParameter)->getflux(reversedIndex);
                DegreeOfPolarizationVariance = TemporaryPolarimetry.getDegreeOfPolarization(StokesParameter)->getvariance(reversedIndex);
                NullSpectrum1 = TemporaryPolarimetry.getFirstNullPolarization(StokesParameter)->getflux(reversedIndex);
                NullSpectrum1Variance = TemporaryPolarimetry.getFirstNullPolarization(StokesParameter)->getvariance(reversedIndex);
                NullSpectrum2 = TemporaryPolarimetry.getSecondNullPolarization(StokesParameter)->getflux(reversedIndex);
                NullSpectrum2Variance = TemporaryPolarimetry.getSecondNullPolarization(StokesParameter)->getvariance(reversedIndex);
                
                Polarimetry->setStokesParameter(StokesI, Intensity, Variance, index);
                Polarimetry->setStokesParameter(StokesParameter, Polarization, PolarizationVariance, index);
                Polarimetry->setDegreeOfPolarization(StokesParameter, DegreeOfPolarization, DegreeOfPolarizationVariance, index);
                Polarimetry->setFirstNullPolarization(StokesParameter, NullSpectrum1, NullSpectrum1Variance, index);
                Polarimetry->setSecondNullPolarization(StokesParameter, NullSpectrum2, NullSpectrum2Variance, index);
                
                reversedIndex--;
            }
        }
	}
}

void operaSpectralOrderVector::readOrdersFromWavelength(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	//
	// read in the wavelength polynomial coefficients as written by operaWavelengthCalibration
	//
	operaistream fwcalref(filename.c_str());		
	if (fwcalref.is_open()) {
		unsigned line = 0;
		unsigned order, npar;
		string dataline;
		while (fwcalref.good()) {
			getline(fwcalref, dataline);
			if (strlen(dataline.c_str())) {
				stringstream ss (stringstream::in | stringstream::out);
				ss << dataline.c_str();
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else if (line == 0) {				// number of entries
					ss >> count;
					line++;
				} else {
					ss >> order;
					ss >> npar;									
#ifdef RANGE_CHECK
					if (order > length) {
						throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
					}
#endif
					if (order < minorder || minorder == 0) {
						minorder = order;
					}
					if (order > maxorder || maxorder == 0) {
						maxorder = order;
					}
					operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
					if (!spectralOrder->getWavelength()) {
						spectralOrder->createWavelength(MAXORDEROFWAVELENGTHPOLYNOMIAL);
					}
					operaWavelength *wavelength = spectralOrder->getWavelength();         
					Polynomial *p = wavelength->getWavelengthPolynomial();
					p->setOrderOfPolynomial(npar);
					Float coeff = 0.0;
					Float coefferr = 0.0;
					for (unsigned i=0; i<npar; i++) {
						ss >> coeff;
						ss >> coefferr;
						p->setCoefficient(i, coeff.f);
						p->setCoefficientError(i, coefferr.f);
					}
					
					spectralOrder->sethasWavelength(true);
					line++;
				}
			}
		}	
		fwcalref.close();
	}
}

/*
 * Read Orders From Instrument Profile
 */
void operaSpectralOrderVector::readOrdersFromProfile(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	
	string dataline;
	unsigned line = 0;
	unsigned order = 0, npar = 0, nx = 0, ny = 0, xsize = 0, xsampling = 0, ysize = 0, ysampling = 0;
	unsigned col = 0, row = 0, ndatapoints = 0;
	unsigned lastorder = 0;
	Float chisqr;
	Float ProfilePolynomial[MAXPOLYNOMIAL];
	operaSpectralOrder *spectralOrder = NULL;
	operaInstrumentProfile *instrumentProfile = NULL;
	if (!filename.empty()) {
		operaistream flist(filename.c_str());		
		if (flist.is_open()) {
			while (flist.good()) {
				getline (flist,dataline);
				if (strlen(dataline.c_str())) {
					stringstream ss (stringstream::in | stringstream::out);
					ss << dataline.c_str();
					if (dataline.c_str()[0] == '#') {
						// skip comments
					} else if (line == 0) {				// number of entries
						sscanf(dataline.c_str(), "%u %u %u %u %u %u %u", &count, &nx, &ny, &xsize, &xsampling, &ysize, &ysampling);
						line++;
					} else {
						sscanf(dataline.c_str(), "%u %u %u %u", &order, &col, &row, &npar);
						switch (npar) {
							case 1:
								//sscanf(dataline.c_str(), "%u %u %u %u %u %g %g", &order, &col, &row, &npar, &ndatapoints, &ProfilePolynomial[0], &chisqr);
								ss >> order >> col >> row >> npar >> ndatapoints >> ProfilePolynomial[0] >> chisqr;
								break;
							case 2:
								//sscanf(dataline.c_str(), "%u %u %u %u %u %g %g %g", &order, &col, &row, &npar, &ndatapoints, &ProfilePolynomial[0], &ProfilePolynomial[1], &chisqr);
								ss >> order >> col >> row >> npar >> ndatapoints >> ProfilePolynomial[0] >> ProfilePolynomial[1] >> chisqr;
								break;
							case 3:
								//sscanf(dataline.c_str(), "%u %u %u %u %u %g %g %g %g", &order, &col, &row, &npar, &ndatapoints, &ProfilePolynomial[0], &ProfilePolynomial[1], &ProfilePolynomial[2], &chisqr);
								ss >> order >> col >> row >> npar >> ndatapoints >> ProfilePolynomial[0] >> ProfilePolynomial[1] >> ProfilePolynomial[2] >> chisqr;
								break;
							case 4:
								//sscanf(dataline.c_str(), "%u %u %u %u %u %g %g %g %g %g", &order, &col, &row, &npar, &ndatapoints, &ProfilePolynomial[0], &ProfilePolynomial[1], &ProfilePolynomial[2], &ProfilePolynomial[3], &chisqr);
								ss >> order >> col >> row >> npar >> ndatapoints >> ProfilePolynomial[0] >> ProfilePolynomial[1] >> ProfilePolynomial[2] >> ProfilePolynomial[3] >> chisqr;
								break;
							case 5:
								//sscanf(dataline.c_str(), "%u %u %u %u %u %g %g %g %g %g %g", &order, &col, &row, &npar, &ndatapoints, &ProfilePolynomial[0], &ProfilePolynomial[1], &ProfilePolynomial[2], &ProfilePolynomial[3], &ProfilePolynomial[4], &chisqr);
								ss >> order >> col >> row >> npar >> ndatapoints >> ProfilePolynomial[0] >> ProfilePolynomial[1] >> ProfilePolynomial[2] >> ProfilePolynomial[3] >> ProfilePolynomial[4] >> chisqr;
								break;
							case 6:
								//sscanf(dataline.c_str(), "%u %u %u %u %u %g %g %g %g %g %g %g", &order, &col, &row, &npar, &ndatapoints, &ProfilePolynomial[0], &ProfilePolynomial[1], &ProfilePolynomial[2], &ProfilePolynomial[3], &ProfilePolynomial[4], &ProfilePolynomial[5], &chisqr);
								ss >> order >> col >> row >> npar >> ndatapoints >> ProfilePolynomial[0] >> ProfilePolynomial[1] >> ProfilePolynomial[2] >> ProfilePolynomial[3] >> ProfilePolynomial[4] >> ProfilePolynomial[5] >> chisqr;
								break;
							default:
								throw operaException("operaSpectralOrder: profile polynomial order must be between 1 and 6, got "+itos(npar), operaErrorInvalidParameter, __FILE__, __FUNCTION__, __LINE__);	
								break;
						}
						if (order != lastorder) {
							spectralOrder = GetSpectralOrder(order);
							spectralOrder->setInstrumentProfileVector(xsize, xsampling, ysize, ysampling, 1);
							instrumentProfile = spectralOrder->getInstrumentProfile();
							spectralOrder->sethasInstrumentProfile(true);
						}
						PolynomialCoeffs_t *pp = (PolynomialCoeffs_t *)malloc(sizeof(PolynomialCoeffs_t));
						if (!pp) {
							throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
						}
						pp->orderofPolynomial = npar;
						for (unsigned coeff=0; coeff<npar; coeff++) {
							pp->p[coeff] = ProfilePolynomial[coeff].f;
						}
						instrumentProfile->setipPolyModelCoefficients(pp, col, row);
						instrumentProfile->setchisqrMatrixValue(chisqr.f, col, row);
#ifdef RANGE_CHECK
						if (order > length) {
							throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
						}
#endif
						if (order < minorder || minorder == 0) {
							minorder = order;
						}
						if (order > maxorder || maxorder == 0) {
							maxorder = order;
						}
						lastorder = order;
						line++;
					}
				}
			}
		}	
		flist.close();
	}
}

/*
 * NOT COMPLETE
 */
void operaSpectralOrderVector::readOrdersFromLines(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder) {
	
	string dataline;
	unsigned line = 0;
	unsigned order = 0;
	unsigned nFeatures = 0;
	unsigned featurenumber = 0;
	unsigned nlinesinfeature = 0;
	unsigned lineinfeature = 0;
	Double center = 0.0;
	Double centererror = 0.0;
	Double sigma = 0.0;
	Double sigmaerror = 0.0;
	Double amplitude = 0.0;
	Double amplitudeerror = 0.0;
	Double chisqr = 0.0;
	double *centervector = NULL;
	double *sigmavector = NULL;
	double *amplitudevector = NULL;        
	double *centererrorvector = NULL;
	double *sigmaerrorvector = NULL;
	double *amplitudeerrorvector = NULL; 
	unsigned lastorder = 0;
	int lastfeature = -1;
	operaSpectralLines *spectralLines = NULL;
	operaSpectralOrder *spectralOrder = NULL;
	operaSpectralFeature *spectralFeature = NULL;
	if (!filename.empty()) {
		operaistream flist(filename.c_str());		
		if (flist.is_open()) {
			while (flist.good()) {
				getline (flist,dataline);
				if (strlen(dataline.c_str())) {
					if (dataline.c_str()[0] == '#') {
						// skip comments
					} else if (line == 0) {				// number of entries
						sscanf(dataline.c_str(), "%u", &count);
						line++;
					} else {
						stringstream ss (stringstream::in | stringstream::out);
						ss << dataline.c_str();
						//sscanf(dataline.c_str(), "%u %u %u %u %u %lf %lf %lf %lf %lf %lf %lf", &order, &nFeatures, &featurenumber, &nlinesinfeature, &lineinfeature, &center, &centererror, &sigma, &sigmaerror, &amplitude, &amplitudeerror, &chisqr);
						ss >> order >> nFeatures >> featurenumber >> nlinesinfeature >> lineinfeature >> center >> centererror >> sigma >> sigmaerror >> amplitude >> amplitudeerror >> chisqr;
						// <number of orders>;
						// <order number> <nfeatures> <feature number> <nlines in feature> <line number in feature> <center> <error> <sigma> <error> <amplitude> <error> <chisqr>
						if (order != lastorder) {
							spectralOrder = GetSpectralOrder(order);
							operaSpectralLines *SpectralLines = spectralOrder->getSpectralLines();
							if (SpectralLines)
								delete SpectralLines;
							spectralLines = new operaSpectralLines(nFeatures, MAXLINESPERFEATURE);
							spectralOrder->setSpectralLines(spectralLines);
							spectralOrder->sethasSpectralLines(true);
							lastfeature = -1;
						}
						if (lastfeature != (int)featurenumber) {
							spectralFeature = spectralLines->getSpectralFeature(featurenumber);
							spectralFeature->setnLines(nlinesinfeature);
							centervector = spectralFeature->getGaussianFit()->getCenterVector();
							sigmavector = spectralFeature->getGaussianFit()->getSigmaVector();
							amplitudevector = spectralFeature->getGaussianFit()->getAmplitudeVector();        
							centererrorvector = spectralFeature->getGaussianFit()->getCenterErrorVector();
							sigmaerrorvector = spectralFeature->getGaussianFit()->getSigmaErrorVector();
							amplitudeerrorvector = spectralFeature->getGaussianFit()->getAmplitudeErrorVector(); 
							spectralFeature->getGaussianFit()->setGaussianChisqr(chisqr.d);
						}
						*centervector++ = center.d; 
						*centererrorvector++ = centererror.d; 
						*sigmavector++ = sigma.d; 
						*sigmaerrorvector++ = sigmaerror.d; 
						*amplitudevector++ = amplitude.d; 
						*amplitudeerrorvector++ = amplitudeerror.d; 
#ifdef RANGE_CHECK
						if (order > length) {
							throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
						}
#endif
						if (order < minorder || minorder == 0) {
							minorder = order;
						}
						if (order > maxorder || maxorder == 0) {
							maxorder = order;
						}
						lastorder = order;
						lastfeature = featurenumber;
						line++;
					}
				}
			}
		}	
		flist.close();
	}
}

void operaSpectralOrderVector::readOrdersFromCSV(string filename) {
	string dataline;
	unsigned line = 0;
	unsigned order = 0;
	unsigned lastorder = 0;
    char Object[48];
    char Mode[24];
	float Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance,DegreeOfPolarization,StokesParameter,NullSpectrum1,NullSpectrum2,PolarizationError;

	//operaSpectralOrder *spectralOrder = NULL;
	if (!filename.empty()) {
		operaistream flist(filename.c_str());
		if (flist.is_open()) {
			while (flist.good()) {
				getline (flist,dataline);
				if (strlen(dataline.c_str())) {
					if (dataline.c_str()[0] == '#') {
						sscanf(dataline.c_str(), ",%s", Object);
					} else {
						sscanf(dataline.c_str(), ",,%s", Mode);
                        // pol
                        // #!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance,DegreeOfPolarization,StokesParameter,NullSpectrum1,NullSpectrum2,PolarizationError
                        // #!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance
						if (startsWith(Mode, "polar")) {
                            sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR, &LeftBeamFlux, &RightBeamFlux, &LeftBeamFluxVariance, &RightBeamFluxvariance, &DegreeOfPolarization, &StokesParameter, &NullSpectrum1, &NullSpectrum2, &PolarizationError);
                        } else if (startsWith(Mode, "pol")) {
                            sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR, &LeftBeamFlux, &RightBeamFlux, &LeftBeamFluxVariance, &RightBeamFluxvariance);
                        }
                        // star plus sky
                        // #!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR,LeftBeamFlux,RightBeamFlux,LeftBeamFluxVariance,RightBeamFluxvariance
						if (startsWith(Mode, "sp1")) {
                            sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR, &LeftBeamFlux, &RightBeamFlux, &LeftBeamFluxVariance, &RightBeamFluxvariance);
                        }
                        // star only
                        // #!csv,"+object+",Mode,Order,Wavelength,Flux,NormalizedFlux,FluxVariance,SNR
						if (startsWith(Mode, "sp2")) {
                            sscanf(dataline.c_str(), ",,,%u,%f,%f,%f,%f,%f", &order, &Wavelength, &Flux, &NormalizedFlux, &FluxVariance, &SNR);
                        }
                        if (order != lastorder) {
							//spectralOrder = GetSpectralOrder(order);
						}
#ifdef RANGE_CHECK
						if (order > length) {
							throw operaException("operaSpectralOrderVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
						}
#endif
						if (order < minorder || minorder == 0) {
							minorder = order;
						}
						if (order > maxorder || maxorder == 0) {
							maxorder = order;
						}
						lastorder = order;
						line++;
					}
				}
			}
		}
		flist.close();
	}
}

void operaSpectralOrderVector::readOrderSpacingPolynomial(string filename) {
	string dataline;
	if (!filename.empty()) {
		operaistream flist(filename.c_str());		
		if (flist.is_open()) {
			unsigned line = 0;
			unsigned npar = 0;
			double orderSpacingPolynomial[MAXPOLYNOMIAL];
			double orderSpacingPolynomialErrors[MAXPOLYNOMIAL];
			double chisqr = 0.0;
			while (flist.good()) {
				getline (flist,dataline);
				if (strlen(dataline.c_str())) {
					if (dataline.c_str()[0] == '#') {
						// skip comments
					} else {
                        sscanf(dataline.c_str(), "%u", &npar);
                        switch (npar) {
                            case 2:
                                sscanf(dataline.c_str(), "%u %lg %lg %lg %lg", &npar, &orderSpacingPolynomial[0], &orderSpacingPolynomialErrors[0], &orderSpacingPolynomial[1], &orderSpacingPolynomialErrors[1]);
                                break;
                            case 3:
                                sscanf(dataline.c_str(), "%u %lg %lg %lg %lg %lg %lg", &npar, &orderSpacingPolynomial[0], &orderSpacingPolynomialErrors[0], &orderSpacingPolynomial[1], &orderSpacingPolynomialErrors[1], &orderSpacingPolynomial[2], &orderSpacingPolynomialErrors[2]);
                                break;
                            case 4:
                                sscanf(dataline.c_str(), "%u %lg %lg %lg %lg %lg %lg %lg %lg", &npar, &orderSpacingPolynomial[0], &orderSpacingPolynomialErrors[0], &orderSpacingPolynomial[1], &orderSpacingPolynomialErrors[1], &orderSpacingPolynomial[2], &orderSpacingPolynomialErrors[2], &orderSpacingPolynomial[3], &orderSpacingPolynomialErrors[3]);
                                break;
                            case 5:
                                sscanf(dataline.c_str(), "%u %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &npar, &orderSpacingPolynomial[0], &orderSpacingPolynomialErrors[0], &orderSpacingPolynomial[1], &orderSpacingPolynomialErrors[1], &orderSpacingPolynomial[2], &orderSpacingPolynomialErrors[2], &orderSpacingPolynomial[3], &orderSpacingPolynomialErrors[3], &orderSpacingPolynomial[4], &orderSpacingPolynomialErrors[4]);
                                break;
                            case 6:
                                sscanf(dataline.c_str(), "%u %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &npar, &orderSpacingPolynomial[0], &orderSpacingPolynomialErrors[0], &orderSpacingPolynomial[1], &orderSpacingPolynomialErrors[1], &orderSpacingPolynomial[2], &orderSpacingPolynomialErrors[2], &orderSpacingPolynomial[3], &orderSpacingPolynomialErrors[3], &orderSpacingPolynomial[4], &orderSpacingPolynomialErrors[4], &orderSpacingPolynomial[5], &orderSpacingPolynomialErrors[5]);
                                break;
                            default:
                                throw operaException("operaSpectralOrder: polynomial order must be between 2 and 6", operaErrorInvalidParameter, __FILE__, __FUNCTION__, __LINE__);
                                break;
                        }
                        Polynomial *p = getOrderSpacingPolynomial();
                        p->setOrderOfPolynomial(npar);
                        p->setChisqr(chisqr);
                        for (unsigned i=0; i<npar; i++) {
                            p->setCoefficient(i, orderSpacingPolynomial[i]);
                            p->setCoefficientError(i, orderSpacingPolynomialErrors[i]);
                        }
						line++;
					}
				}									
			}	
			flist.close();
		}
	}		
}

void operaSpectralOrderVector::readDispersionPolynomial(string filename) {
	string dataline;
	if (!filename.empty()) {
		operaistream flist(filename.c_str());
		if (flist.is_open()) {
			unsigned line = 0;
			int MinorderOfLaurentPolynomial = 0;
            int MaxorderOfLaurentPolynomial = 0;
			double dispersionPolynomialCoefficients[MAXPOLYNOMIAL];
			double dispersionPolynomialErrors[MAXPOLYNOMIAL];
			double chisqr = 0.0;
            unsigned dispIndex,NumberOfDispersionPolynomials;
            
			while (flist.good()) {
				getline (flist,dataline);
				if (strlen(dataline.c_str())) {
					if (dataline.c_str()[0] == '#') {
						// skip comments
					} else if (line == 0) {	  // number of dispersion polynomials
						sscanf(dataline.c_str(), "%u", &NumberOfDispersionPolynomials);
                        setnumberOfDispersionPolynomials(NumberOfDispersionPolynomials);
						line++;
					} else {
                        sscanf(dataline.c_str(), "%u %d %d", &dispIndex, &MinorderOfLaurentPolynomial, &MaxorderOfLaurentPolynomial);
                        unsigned npar = (unsigned)(MaxorderOfLaurentPolynomial - MinorderOfLaurentPolynomial) + 1;
                        switch (npar) {
                            case 1:
                                sscanf(dataline.c_str(), "%u %d %d %lg %lg", &dispIndex, &MinorderOfLaurentPolynomial, &MaxorderOfLaurentPolynomial, &dispersionPolynomialCoefficients[0], &dispersionPolynomialErrors[0]);
                                break;
                            case 2:
                                sscanf(dataline.c_str(), "%u %d %d %lg %lg %lg %lg", &dispIndex, &MinorderOfLaurentPolynomial, &MaxorderOfLaurentPolynomial, &dispersionPolynomialCoefficients[0], &dispersionPolynomialErrors[0], &dispersionPolynomialCoefficients[1], &dispersionPolynomialErrors[1]);
                                break;
                            case 3:
                                sscanf(dataline.c_str(), "%u %d %d %lg %lg %lg %lg %lg %lg", &dispIndex, &MinorderOfLaurentPolynomial, &MaxorderOfLaurentPolynomial, &dispersionPolynomialCoefficients[0], &dispersionPolynomialErrors[0], &dispersionPolynomialCoefficients[1], &dispersionPolynomialErrors[1], &dispersionPolynomialCoefficients[2], &dispersionPolynomialErrors[2]);
                                break;
                            case 4:
                                sscanf(dataline.c_str(), "%u %d %d %lg %lg %lg %lg %lg %lg %lg %lg", &dispIndex, &MinorderOfLaurentPolynomial, &MaxorderOfLaurentPolynomial, &dispersionPolynomialCoefficients[0], &dispersionPolynomialErrors[0], &dispersionPolynomialCoefficients[1], &dispersionPolynomialErrors[1], &dispersionPolynomialCoefficients[2], &dispersionPolynomialErrors[2], &dispersionPolynomialCoefficients[3], &dispersionPolynomialErrors[3]);
                                break;
                            case 5:
                                sscanf(dataline.c_str(), "%u %d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &dispIndex, &MinorderOfLaurentPolynomial, &MaxorderOfLaurentPolynomial, &dispersionPolynomialCoefficients[0], &dispersionPolynomialErrors[0], &dispersionPolynomialCoefficients[1], &dispersionPolynomialErrors[1], &dispersionPolynomialCoefficients[2], &dispersionPolynomialErrors[2], &dispersionPolynomialCoefficients[3], &dispersionPolynomialErrors[3], &dispersionPolynomialCoefficients[4], &dispersionPolynomialErrors[4]);
                                break;
                            case 6:
                                sscanf(dataline.c_str(), "%u %d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &dispIndex, &MinorderOfLaurentPolynomial, &MaxorderOfLaurentPolynomial, &dispersionPolynomialCoefficients[0], &dispersionPolynomialErrors[0], &dispersionPolynomialCoefficients[1], &dispersionPolynomialErrors[1], &dispersionPolynomialCoefficients[2], &dispersionPolynomialErrors[2], &dispersionPolynomialCoefficients[3], &dispersionPolynomialErrors[3], &dispersionPolynomialCoefficients[4], &dispersionPolynomialErrors[4], &dispersionPolynomialCoefficients[5], &dispersionPolynomialErrors[5]);
                                break;
                            default:
                                throw operaException("operaSpectralOrder: polynomial order must be between 2 and 6", operaErrorInvalidParameter, __FILE__, __FUNCTION__, __LINE__);
                                break;
                        }
                        LaurentPolynomial *p = getDispersionPolynomial(dispIndex);
                        p->setMinMaxOrderOfLaurentPolynomial(MinorderOfLaurentPolynomial,MaxorderOfLaurentPolynomial);
                        p->setChisqr(chisqr);
                        for (unsigned i=0; i<p->getNumberOfCoefficients(); i++) {
                            p->setCoefficient(i, dispersionPolynomialCoefficients[i]);
                            p->setCoefficientError(i, dispersionPolynomialErrors[i]);
                        }
						line++;
					}
				}									
			}	
			
		}
	}		
}

void operaSpectralOrderVector::readRadialVelocityCorrection(string filename) {
	string dataline;
	if (!filename.empty()) {
		operaistream flist(filename.c_str());
		if (flist.is_open()) {
			double rvel;
            
			while (flist.good()) {
				getline (flist,dataline);
				if (strlen(dataline.c_str())) {
					if (dataline.c_str()[0] == '#') {
						// skip comments
					} else {
						sscanf(dataline.c_str(), "%lf", &rvel);
                        setBarycentricRadialVelocityCorrection(rvel);
					}
				}									
			}	
			
		}
	}		
}

void operaSpectralOrderVector::fitOrderSpacingPolynomial(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, unsigned sampleCenterPosition, unsigned referenceOrderNumber, float referenceOrderSeparation, int detectionMethod, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2, unsigned cleanbinsize, float nsigcut, ostream *pout) {
    
    unsigned nx = x2 - x1;
    unsigned ny = y2 - y1;
    
    float *fx = (float *) malloc (nx * sizeof(float));
    float *fy = (float *) malloc (ny * sizeof(float));
    float *fytmp = (float *) malloc (ny * sizeof(float));
    float *fyerr = (float *) malloc (ny * sizeof(float));
	if (!fyerr) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
	}
    float xmean[MAXORDERS],xmeanerr[MAXORDERS],ymean[MAXORDERS];
    
    float TemporaryOrderPosition[MAXORDERS], TemporaryOrderSeparation[MAXORDERS],TemporaryOrderSeparationError[MAXORDERS];
    
    unsigned nords=0;
    float threshold = DETECTTHRESHOLD;
    float sigma = slit/4;
    
    float *fys = (float *) malloc ((nsamples+1)* sizeof(float));
	if (!fys) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);
	}
    
    // Sample nsamples (default 3) equally spaced regions of the detector at e.g. 1/4, 1/2, 3/4 in the vertical direction
    // what we are doing here is taking three(NSAMPLE) samples to try to determine:
    // what is the polynomial that descibes the function of spacing between orders in the x direction
    unsigned y0,yf;
    if(sampleCenterPosition - nsamples/2 < y1) {
        y0 = y1;
        yf = y1 + nsamples;
    } else {
        y0 = sampleCenterPosition - nsamples/2;
        yf = y0 + nsamples;
    }
    if(yf > y2) {
        yf=y2;
    }
    
    unsigned np = 0;
    
    for (unsigned x=x1; x<nx; x++) {
        fx[np] = (float)x + 0.5;
        
        unsigned ns=0;
        for (unsigned y=y0; y<yf; y++) {
            if(badpixImage[y][x] == 1) {
                fys[ns++] = masterFlatImage[y][x];
            }
        }
        
        if (ns == 0) {
            if(FFTfilter){
                fytmp[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
                fyerr[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
            } else {
                fy[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
                fyerr[np] = 0.0;	// SHOULD BE NAN, but library doesn't handle it...
            }
        } else {
            if(FFTfilter){
                fytmp[np] = operaArrayMedian(ns,fys);
                fyerr[np] = operaArrayMedianSigma(ns, fys, fytmp[np]);
            } else {
                fy[np] = operaArrayMedian(ns,fys);
                fyerr[np] = operaArrayMedianSigma(ns, fys, fy[np]);
            }
        }
        
#ifdef PRINT_DEBUG
        cerr << np << ' ' << fx[np] << ' ' << fy[np] << ' ' << endl;
#endif
        np++;
    }
    
    if(FFTfilter){
        operaFFTLowPass(np,fytmp,fy, 0.1);
    }
    
    if(detectionMethod == 1) {
#ifdef WITH_ERRORS
        nords = operaCCDDetectPeaksWithErrorsUsingGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean,xmeanerr);
#else
        nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean);
#endif
    } else if (detectionMethod == 2) {
#ifdef WITH_ERRORS
        nords = operaCCDDetectPeaksWithErrorsUsingGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean,xmeanerr);
#else
        nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean);
#endif
    } else if (detectionMethod == 3) {
#ifdef WITH_ERRORS
        nords = operaCCDDetectPeaksWithErrorsUsingTopHat(np,fx,fy,(unsigned)slit,noise,gain,threshold,xmean,ymean,xmeanerr);
#else
        nords = operaCCDDetectPeaksWithTopHat(np,fx,fy,(unsigned)slit,noise,gain,threshold,xmean,ymean);
#endif
    }
    
    unsigned npspc = 0;
    
    for(unsigned i=1;i<nords;i++) {
        TemporaryOrderSeparation[npspc] = fabs(xmean[i] - xmean[i-1]);
        TemporaryOrderSeparationError[npspc] = sqrt(xmeanerr[i]*xmeanerr[i] + xmeanerr[i-1]*xmeanerr[i-1]);
        TemporaryOrderPosition[npspc] = xmean[i-1] + (xmean[i] + xmean[i-1])/2;
        
        if(i > 1) {
            float ratio = TemporaryOrderSeparation[npspc]/TemporaryOrderSeparation[npspc - 1];
            int roundedRatio = (int)round(ratio);
            
            // If ratio>1 it means the current separation is at least twice larger than the size of
            // previous separation. Therefore the current has skipped one (or more) order(s).
            if(roundedRatio > 1
               && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc]
               && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc - 1]) {
                TemporaryOrderSeparation[npspc] /= float(roundedRatio);
                
                // Else if ratio<1 it means the current separation is at least twice smaller than the size
                // of previous separations. In this case we need to fix all previous separations.
            } else if (roundedRatio < 1
                       && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc]
                       && TemporaryOrderSeparationError[npspc] < 0.5*TemporaryOrderSeparation[npspc - 1]) {
                unsigned j = npspc - 1;
                
                // Fix previous separation
                TemporaryOrderSeparation[j] /= 1.0/float(roundedRatio);
                
                //Check for separations before previous one
                for (unsigned ii = 0; ii < i-2; ii++) {
                    ratio = TemporaryOrderSeparation[j]/TemporaryOrderSeparation[j-1];
                    roundedRatio = (int)round(ratio);
                    if(roundedRatio > 1
                       && TemporaryOrderSeparationError[j] < 0.5*TemporaryOrderSeparation[j]
                       && TemporaryOrderSeparationError[j] < 0.5*TemporaryOrderSeparation[j - 1]) {
                        // Fix separation value.
                        TemporaryOrderSeparation[j] /= float(roundedRatio);
                    }
                    j--;
                }
            }
        }
        npspc++;
    }
    
    int SortIndex[MAXORDERS];
    operaArrayIndexSort((int)npspc,TemporaryOrderPosition,SortIndex);
    
    double OrderPosition[MAXORDERS];
    double OrderSeparation[MAXORDERS],OrderSeparationError[MAXORDERS];
    
    // identify reference order, which will be the order for which OrderSeparation ~ referenceOrderSeparation
    float minResidualSeparation = BIG;
    unsigned minj = 0;
    
    for(unsigned j=0; j<npspc; j++) {
        if(fabs(TemporaryOrderSeparation[SortIndex[j]] - referenceOrderSeparation) < minResidualSeparation) {
            minResidualSeparation = fabs(TemporaryOrderSeparation[SortIndex[j]] - referenceOrderSeparation);
            minj = j;
        }
        
        //OrderPosition[j] = (double)TemporaryOrderPosition[SortIndex[j]];
        OrderSeparation[j] = (double)TemporaryOrderSeparation[SortIndex[j]];
        OrderSeparationError[j] = (double)TemporaryOrderSeparationError[SortIndex[j]];
    }
    
    // Below we are replacing the variable detector position to order number, which
    // is based on the reference order provided by user -> referenceOrderNumber
    
    for(unsigned j=0; j<npspc; j++) {
        OrderPosition[j] = double(referenceOrderNumber + j - minj);
    }
	
    // here is the polynomial fit of spacing between orders across the array in the x direction
    // npar is different in this case than order tracing
    unsigned npars = 3;
    double par[3] = {1,1,1};
    double ecoeffs[3] = {0,0,0};
    double chisqr = 0.0;
    
#ifdef WITH_ERRORS
    int errorcode = operaMPFitPolynomial(npspc, OrderPosition, OrderSeparation, OrderSeparationError, npars, par, ecoeffs, &chisqr);
    if (errorcode <= 0) {
        throw operaException("operaSpectralOrderVector: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);
    }
    
#else
    operaLMFitPolynomial(npspc, OrderPosition, OrderSeparation, npars, par, &chisqr);
#endif
    
    
    /*
     *  Below it applies a sigma clip cleaning
     */
    unsigned binsize = cleanbinsize;
    float nsig = nsigcut;
    
    if (binsize==0 || npspc==0) {
        throw operaException("operaSpectralOrderVector: binsize=0 or nDataPoints=0", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);
    }
    unsigned nDataPoints = npspc;
    
    double *cleanxdataVector = new double[nDataPoints];
    double *cleanydataVector = new double[nDataPoints];
    double *cleanyerrorVector = new double[nDataPoints];

    unsigned numberOfCleanPoints = 0;
    
    float *xtmp = new float[nDataPoints];
    float *ytmp = new float[nDataPoints];
    
    for(unsigned i=0; i<nDataPoints; i++) {
        
        int firstPoint = (int)i - (int)binsize;
        int lastPoint = (int)i + (int)binsize + 1;
        
        if(firstPoint < 0) {
            firstPoint = 0;
            lastPoint = 2*(int)binsize + 1;
        }
        if(lastPoint > (int)nDataPoints) {
            lastPoint = (int)nDataPoints;
            firstPoint = (int)nDataPoints - 2*(int)binsize - 1;
            if(firstPoint < 0) {
                firstPoint = 0;
            }
        }
        
        unsigned np = 0;
        for(unsigned ii=(unsigned)firstPoint; ii<(unsigned)lastPoint; ii++) {
            xtmp[np] = (float)OrderPosition[ii];
            ytmp[np] = (float)OrderSeparation[ii];
            np++;
        }
        
        float am,bm,abdevm;
        
        //--- Robust linear fit
        ladfit(xtmp,ytmp,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
        
        //--- Clean up
        float fitMedianSlope = (bm*(float)OrderPosition[i] + am);
        
        if(fabs((float)OrderSeparation[i] - fitMedianSlope) < nsig*abdevm) {
            cleanxdataVector[numberOfCleanPoints] = OrderPosition[i];
            cleanydataVector[numberOfCleanPoints] = OrderSeparation[i];
#ifdef WITH_ERRORS

            cleanyerrorVector[numberOfCleanPoints] = OrderSeparationError[i];  
#else
            cleanyerrorVector[numberOfCleanPoints] = 0.0;
#endif
            numberOfCleanPoints++;
        }
    }
    
    for(unsigned i=0; i<numberOfCleanPoints; i++) {
        OrderPosition[i] = cleanxdataVector[i];
        OrderSeparation[i] = cleanydataVector[i];
        OrderSeparationError[i] = cleanyerrorVector[i];
    }
    nDataPoints = numberOfCleanPoints;
    
    delete[] cleanxdataVector;
    delete[] cleanydataVector;
    delete[] cleanyerrorVector;
    
    npspc = nDataPoints;
    /*
     *  End of sigma clip cleaning
     */
    
    
    /*
     * Perform polynomial fit on clean data
     */
#ifdef WITH_ERRORS
    errorcode = operaMPFitPolynomial(npspc, OrderPosition, OrderSeparation, OrderSeparationError, npars, par, ecoeffs, &chisqr);
    if (errorcode <= 0) {
        throw operaException("operaGeometryCalibration: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);
    }
#else
    operaLMFitPolynomial(npspc, OrderPosition, OrderSeparation, npars, par, &chisqr);
#endif
    
    if (pout != NULL) {
        for(unsigned i=0;i<npspc;i++) {
            *pout <<  OrderPosition[i] << ' ' << OrderSeparation[i] << ' ' << OrderSeparationError[i] << ' ' << PolynomialFunction(OrderPosition[i],par,npars) << endl;
        }
    }
    
    // set in the order spacing polynomial
    PolynomialCoeffs_t orderSpacing;
    orderSpacing.orderofPolynomial = npars;
    orderSpacing.polychisqr = chisqr;
    for (unsigned i=0; i<npars; i++) {
        orderSpacing.p[i] = par[i];
        orderSpacing.e[i] = ecoeffs[i];
    }
    setOrderSpacingPolynomial(&orderSpacing);
    free(fx);
    free(fy);
    free(fytmp);
    free(fyerr);
}

void operaSpectralOrderVector::measureIPAlongRowsFromSamples(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2,float *ipfunc, float *ipx, float *iperr) {
	
    unsigned nx = x2 - x1;
    unsigned ny = y2 - y1;    
    
    float *fx = (float *) malloc (nx * sizeof(float));
    float *fy = (float *) malloc (ny * sizeof(float));
    float *fytmp = (float *) malloc (ny * sizeof(float));		
	if (!fytmp) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    float xmean[MAXORDERS],ymean[MAXORDERS];	
    
    unsigned nords=0;
    float threshold = DETECTTHRESHOLD;
    float sigma = slit/4;
    unsigned uslit = (unsigned)slit;
    
    unsigned numberofpointsinydirectiontomerge = NPINSAMPLE;
    float *fys = (float *) malloc (NPINSAMPLE * sizeof(float)); 
	if (!fys) {
		throw operaException("operaSpectralOrderVector: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
    float *ipsample = new float[uslit];
    float *ipxsample = new float[uslit];
    float *ipvar = new float[uslit];
	
    memset(ipfunc, 0, sizeof(float)*uslit);
    memset(ipx, 0, sizeof(float)*uslit);
    memset(ipvar, 0, sizeof(float)*uslit);
    
    for(unsigned k=0; k<nsamples; k++) {
        unsigned np=0;
        for (unsigned x=x1; x<nx; x++) {			
            unsigned ns=0;
            for (unsigned y=(ny-y1)*(k+1)/(nsamples + 1) - numberofpointsinydirectiontomerge/2; y < (ny-y1)*(k+1)/(nsamples + 1) + numberofpointsinydirectiontomerge/2; y++) {
                if(badpixImage[y][x] == 1)
                    fys[ns++] = masterFlatImage[y][x];
            }
            
            fx[np] = (float)x + 0.5;
            
			if (ns == 0) {
				if(FFTfilter){
					fytmp[np] = 0.0;
				} else {
					fy[np] = 0.0;	
				}				
			} else {
				if(FFTfilter){
					fytmp[np] = operaArrayMedian(ns,fys);
				} else {
					fy[np] = operaArrayMedian(ns,fys);	
				}				
			}
			
#ifdef PRINT_DEBUG
            cerr << k << ' ' << fx[np] << ' ' << fy[np] << ' ' << endl;
#endif
            
            np++;
        }
        
        if(FFTfilter){
            operaFFTLowPass(np,fytmp,fy, 0.1);
        }
        
        nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,noise,gain,threshold,xmean,ymean);
        operaCCDFitIP(np,fx,fy,nords,xmean,ymean,ipsample,ipxsample,iperr,uslit);
        
        for(unsigned i=0;i<uslit;i++) {
            ipfunc[i] += ipsample[i];
            ipx[i] += ipxsample[i];
            ipvar[i] += (iperr[i]*iperr[i]);
        } 
    }
    float IPNormalizationFactor = 0;
    for(unsigned i=0;i<uslit;i++) {
        IPNormalizationFactor += ipfunc[i];
        ipx[i] /= (double)nsamples;
    }     
    for(unsigned i=0;i<uslit;i++) {
        ipfunc[i] /= IPNormalizationFactor;
        iperr[i] = sqrt(ipvar[i]/IPNormalizationFactor);
    }      
    delete[] ipsample;
    delete[] ipxsample;
    delete[] ipvar;
    free(fx);
    free(fy);
    free(fytmp);
}



unsigned operaSpectralOrderVector::getElemIndexAndOrdersByWavelength(int *orderForWavelength, unsigned *elemIndexForWavelength, double wavelength) {
    
    /* E.Martioli -- Jul 1 2013
     * It probably needs a range check for the size of *orderForWavelength and *elemIndexForWavelength, 
     * which must be greater than (maxorder - minorder + 1)
     * But I don't know what is the best way to do that.
     */
    unsigned nOrdersSelected = 0;
    
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
            
            double wl0 = SpectralElements->getwavelength(0);
            double wlf = SpectralElements->getwavelength(SpectralElements->getnSpectralElements()-1);

#ifdef PRINT_DEBUG
            cout << wl0 << ' ' << wlf << ' ' << wavelength << ' ' << order << endl;
#endif

            if((wl0 <= wlf && wavelength >= wl0 && wavelength < wlf) ||
               (wl0 >  wlf && wavelength > wlf && wavelength <= wl0) ) { // for ascending wavelengths
                
                orderForWavelength[nOrdersSelected] = order;
                
                for (unsigned elemIndex=0; elemIndex<SpectralElements->getnSpectralElements()-1; elemIndex++) {
                    double elem_wl0 = SpectralElements->getwavelength(elemIndex);
                    double elem_wlf = SpectralElements->getwavelength(elemIndex+1);
                    
                    if((elem_wl0 <= elem_wlf && wavelength >= elem_wl0 && wavelength < elem_wlf) ||
                       (elem_wl0 >  elem_wlf && wavelength > elem_wlf && wavelength <= elem_wl0) ) {
                        if(fabs(wavelength - elem_wl0) <= fabs(wavelength - elem_wlf)) {
                            elemIndexForWavelength[nOrdersSelected] = elemIndex;
                            break;
                        } else if (fabs(wavelength - elem_wl0) > fabs(wavelength - elem_wlf)) {
                            elemIndexForWavelength[nOrdersSelected] = elemIndex+1;
                            break;
                        }
                    }
                }
                
                nOrdersSelected++;
            }
        }
    }
    return nOrdersSelected;
}


unsigned operaSpectralOrderVector::getOrdersByWavelengthRange(int *orderForWavelengthRange, double Range_wl0, double Range_wlf) {
    
    unsigned nOrdersSelected = 0;
    
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
            
            double wl0 = SpectralElements->getwavelength(0);
            double wlf = SpectralElements->getwavelength(SpectralElements->getnSpectralElements()-1);
                        
#ifdef PRINT_DEBUG
            cout << wl0 << ' ' << wlf << ' ' << Range_wl0 << ' ' << Range_wlf << ' ' << order << endl;
#endif
            
            if((wl0 <= wlf && Range_wl0 >= wl0 && Range_wl0 < wlf) ||
               (wl0 >  wlf && Range_wl0 > wlf && Range_wl0 <= wl0) ||
               (wl0 <= wlf && Range_wlf >= wl0 && Range_wlf < wlf) ||
               (wl0 >  wlf && Range_wlf > wlf && Range_wlf <= wl0)) { // for ascending wavelengths
                
                orderForWavelengthRange[nOrdersSelected] = order;
                nOrdersSelected++;
            }
        }
    }
    return nOrdersSelected;
}


void operaSpectralOrderVector::measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut) {
    
    unsigned maxNDataPoints = 0;
    unsigned numberOfBeams = 0;
    
    int usefulMinorder = (int)minorder;
    int usefulMaxorder = (int)maxorder;
    
    bool hasUsefulMinorder = false;
    
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            
            if(!hasUsefulMinorder) {
                usefulMinorder = order;
                hasUsefulMinorder = true;
            }
            usefulMaxorder = order;
            
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            // Below it populates wavelengths in SpectralElement from wavelength calibration
            SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
            
            numberOfBeams = spectralOrder->getnumberOfBeams();

            spectralOrder->calculateContinuum(binsize,nsigcut,NULL, NULL);
            
            if(spectralOrder->gethasSpectralEnergyDistribution()) {
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                
                if(maxNDataPoints < spectralEnergyDistribution->getnDataPoints()) {
                    maxNDataPoints = spectralEnergyDistribution->getnDataPoints();
                }
            }
        }
    }
    
    operaFluxVector *FluxCalibrationFitVector[MAXORDERS];
    operaFluxVector *BeamFluxCalibrationFitVector[MAXNUMBEROFBEAMS][MAXORDERS];
    
    float *fluxData = new float[MAXORDERS];
    float *orderData = new float[MAXORDERS];
        
    float *fluxBeamData[MAXNUMBEROFBEAMS];
    float *orderBeamData[MAXNUMBEROFBEAMS];
    float orderBeamFlux[MAXNUMBEROFBEAMS];
    unsigned npBeam[MAXNUMBEROFBEAMS];
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        npBeam[beam] = 0;
        fluxBeamData[beam] = new float[MAXORDERS];
        orderBeamData[beam] = new float[MAXORDERS];
    }
    
    float fluxvariance = 0;
    float beamfluxvariance = 0;

    unsigned order_count = 0;
    
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        int loword = order - orderBin;
        int hiord = order + orderBin;
        if(loword < usefulMinorder) {
            loword = usefulMinorder;
        }
        if(hiord > usefulMaxorder) {
            hiord = usefulMaxorder;
        }
        
        int nord = (hiord - loword + 1);
        
        if(!(nord%2)) {
            if(hiord == usefulMaxorder) {
                loword--;
            } else {
                hiord++;
            }
        }
        
        nord = (hiord - loword + 1);
        
        if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution() && GetSpectralOrder(order)->gethasWavelength()) {
            
            unsigned nDataPointsOfCurrentOrder = GetSpectralOrder(order)->getSpectralEnergyDistribution()->getnDataPoints();
            
            FluxCalibrationFitVector[order_count]  = new operaFluxVector(nDataPointsOfCurrentOrder);
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                BeamFluxCalibrationFitVector[beam][order_count] = new operaFluxVector(nDataPointsOfCurrentOrder);
            }
            
            for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
                if(dataIndex>=nDataPointsOfCurrentOrder) {
                    break;
                }
                
                unsigned np = 0;
                float orderFlux = NAN;
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    npBeam[beam] = 0;
                    orderBeamFlux[beam] = NAN;
                }
                
                for(int o=loword; o<=hiord; o++) {
                    operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                    
                    if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()  && spectralOrder->gethasWavelength()) {
                        
                        operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        
                        if(dataIndex < spectralEnergyDistribution->getnDataPoints() && dataIndex < nDataPointsOfCurrentOrder) {
                            float flux = (float)spectralEnergyDistribution->getfluxData(dataIndex);
                            fluxvariance = flux;
                            
                            if(!isnan(flux)) {
                                orderData[np] = (float)o;
                                fluxData[np] = flux;
                                if(o==order) {
                                    orderFlux = flux;
                                }
                                np++;
                            }
                            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                                float beamflux = beamSED->getfluxData(dataIndex);
                                beamfluxvariance = beamflux;
                                
                                if(!isnan(beamflux)) {
                                    orderBeamData[beam][npBeam[beam]] = (float)o;
                                    fluxBeamData[beam][npBeam[beam]] = beamflux;
                                    if(o==order) {
                                        orderBeamFlux[beam] = beamflux;
                                    }
                                    npBeam[beam]++;
                                }
                            }
                        }
                        
                    }
                }
                
                if(np && !(np%2)) {
                    np--;
                }
                if(np>=3) {
                    float am,bm,abdevm;

                    ladfit(orderData,fluxData,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
                    
                    float dytop = 0;
                    
                    for(unsigned i=0;i<np;i++){
                        float fitMedianSlope = (bm*orderData[i] + am);
                        
                        if(fabs(fluxData[i] - fitMedianSlope) < abdevm &&
                           fabs(fluxData[i] - fitMedianSlope) > dytop) {
                            dytop = fabs(fluxData[i] - fitMedianSlope);
                        }
                    }
                    if(dytop == 0) {
                        dytop = abdevm;
                    }
                    
                    FluxCalibrationFitVector[order_count]->setflux(double(bm*(float)order + am + dytop),dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance(double(0.674433*abdevm)*double(0.674433*abdevm),dataIndex);
                } else {
                    FluxCalibrationFitVector[order_count]->setflux((double)orderFlux,dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance((double)fluxvariance,dataIndex);
                }
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam] && !(npBeam[beam]%2)) {
                        npBeam[beam]--;
                    }
                }
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam]>=3) {
                        float amBeam,bmBeam,abdevmBeam;
                        ladfit(orderBeamData[beam],fluxBeamData[beam],npBeam[beam],&amBeam,&bmBeam,&abdevmBeam); /* robust linear fit: f(x) =  a + b*x */
                        float dytop = 0;
                        
                        for(unsigned i=0;i<npBeam[beam];i++){
                            float fitMedianSlope = (bmBeam*orderBeamData[beam][i] + amBeam);
                            
                            if(fabs(fluxBeamData[beam][i] - fitMedianSlope) < abdevmBeam &&
                               (fluxBeamData[beam][i] - fitMedianSlope) > dytop) {
                                dytop = fluxBeamData[beam][i] - fitMedianSlope;
                            }
                        }
                        
                        if(dytop == 0) {
                            dytop = abdevmBeam;
                        }
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(bmBeam*(float)order + amBeam + dytop),dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(0.674433*abdevmBeam)*double(0.674433*abdevmBeam),dataIndex);
                    } else {
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(orderBeamFlux[beam]),dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(beamfluxvariance),dataIndex);
                    }
                } // for(unsigned beam=0; beam < numberOfBeams; beam++) {
            } // for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
            order_count++;
        } // if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution()) {
    } // for(int order=(int)minorder; order<=(int)maxorder; order++) {

    order_count = 0;
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution() && spectralOrder->gethasWavelength()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            unsigned nDataPointsOfCurrentOrder = spectralEnergyDistribution->getnDataPoints();
            
            for(unsigned dataIndex = 0; dataIndex < nDataPointsOfCurrentOrder; dataIndex++) {
                spectralEnergyDistribution->setfluxData(FluxCalibrationFitVector[order_count]->getflux(dataIndex),dataIndex);
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    beamSED->setfluxData(BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex),dataIndex);
                }
#ifdef PRINT_DEBUG
                cout << order << ' ' << dataIndex << ' '
                << spectralEnergyDistribution->getdistanceData(dataIndex) << ' '
                << spectralEnergyDistribution->getfluxData(dataIndex) << ' '
                << FluxCalibrationFitVector[order_count]->getflux(dataIndex) << ' ';
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    cout << beam << ' ' << beamSED->getfluxData(dataIndex) << ' ' << BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex) << ' ';
                }
                cout << endl;
#endif
            }
#ifdef PRINT_DEBUG
            cout << endl;
#endif
            order_count++;
        }
    }

    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution() && spectralOrder->gethasWavelength()) {
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            spectralEnergyDistribution->populateUncalibratedElementsFromContinuumData();
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                beamSED->populateUncalibratedElementsFromContinuumData();
            }
        }
    } //for(int order=(int)minorder; order<=(int)maxorder; order++) {
}


void operaSpectralOrderVector::measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut, unsigned nOrdersPicked, int *orderForWavelength) {
    
    int usefulMinorder = (int)minorder;
    int usefulMaxorder = (int)maxorder;
    
    bool hasUsefulMinorder = false;

    unsigned maxNDataPoints = 0;
    unsigned numberOfBeams = 0;
    
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            
            if(!hasUsefulMinorder) {
                usefulMinorder = order;
                hasUsefulMinorder = true;
            }
            usefulMaxorder = order;
            
            operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
            // Below it populates wavelengths in SpectralElement from wavelength calibration
            SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
            
            numberOfBeams = spectralOrder->getnumberOfBeams();
            
            spectralOrder->calculateContinuum(binsize,nsigcut,NULL, NULL);
            
            if(spectralOrder->gethasSpectralEnergyDistribution()) {
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                
                if(maxNDataPoints < spectralEnergyDistribution->getnDataPoints()) {
                    maxNDataPoints = spectralEnergyDistribution->getnDataPoints();
                }
            }
        }
    }
    
    int actualMinorder = usefulMinorder;
    int actualMaxorder = usefulMaxorder;
    
    if(orderForWavelength[0] < usefulMinorder) {
        actualMinorder = orderForWavelength[0];
    }
    
    if(orderForWavelength[nOrdersPicked-1] > usefulMaxorder) {
        actualMaxorder = orderForWavelength[nOrdersPicked-1];
    }
    
    operaFluxVector *FluxCalibrationFitVector[MAXORDERS];
    operaFluxVector *BeamFluxCalibrationFitVector[MAXNUMBEROFBEAMS][MAXORDERS];
    
    float *fluxData = new float[MAXORDERS];
    float *orderData = new float[MAXORDERS];
    
    float *fluxBeamData[MAXNUMBEROFBEAMS];
    float *orderBeamData[MAXNUMBEROFBEAMS];
    float orderBeamFlux[MAXNUMBEROFBEAMS];
    unsigned npBeam[MAXNUMBEROFBEAMS];
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        npBeam[beam] = 0;
        fluxBeamData[beam] = new float[MAXORDERS];
        orderBeamData[beam] = new float[MAXORDERS];
    }
    
    float fluxvariance = 0;
    float beamfluxvariance = 0;
    
    unsigned order_count = 0;
    
    for(int order=actualMinorder; order<=actualMaxorder; order++) {
        
        int loword = order - orderBin;
        int hiord = order + orderBin;
        if(loword < actualMinorder) {
            loword = actualMinorder;
        }
        if(hiord > actualMaxorder) {
            hiord = actualMaxorder;
        }
        
        int nord = (hiord - loword + 1);
        
        if(!(nord%2)) {
            if(hiord == actualMaxorder) {
                loword--;
            } else {
                hiord++;
            }
        }
        
        nord = (hiord - loword + 1);
                
        if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution()) {
            
            unsigned nDataPointsOfCurrentOrder = GetSpectralOrder(order)->getSpectralEnergyDistribution()->getnDataPoints();
            
            FluxCalibrationFitVector[order_count]  = new operaFluxVector(nDataPointsOfCurrentOrder);
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                BeamFluxCalibrationFitVector[beam][order_count] = new operaFluxVector(nDataPointsOfCurrentOrder);
            }
            
            for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
                if(dataIndex>=nDataPointsOfCurrentOrder) {
                    break;
                }
                
                unsigned np = 0;
                float orderFlux = NAN;
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    npBeam[beam] = 0;
                    orderBeamFlux[beam] = NAN;
                }
                
                for(int o=loword; o<=hiord; o++) {
                    operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                    
                    if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
                        
                        operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        
                        if(dataIndex < spectralEnergyDistribution->getnDataPoints() && dataIndex < nDataPointsOfCurrentOrder) {
                            float flux = spectralEnergyDistribution->getfluxData(dataIndex);
                            fluxvariance = flux;
                            
                            if(!isnan(flux)) {
                                orderData[np] = (float)o;
                                fluxData[np] = flux;
                                if(o==order) {
                                    orderFlux = flux;
                                }
                                np++;
                            }
                            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                                float beamflux = beamSED->getfluxData(dataIndex);
                                beamfluxvariance = beamflux;
                                
                                if(!isnan(beamflux)) {
                                    orderBeamData[beam][npBeam[beam]] = (float)o;
                                    fluxBeamData[beam][npBeam[beam]] = beamflux;
                                    if(o==order) {
                                        orderBeamFlux[beam] = beamflux;
                                    }
                                    npBeam[beam]++;
                                }
                            }
                        }
                        
                    }
                }
                
                if(np && !(np%2)) {
                    np--;
                }
                
                if(np>=3) {
                    float am,bm,abdevm;
                    ladfit(orderData,fluxData,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
                    
                    float dytop = 0;
                    
                    for(unsigned i=0;i<np;i++){
                        float fitMedianSlope = (bm*orderData[i] + am);
                        
                        if(fabs(fluxData[i] - fitMedianSlope) < abdevm &&
                           fabs(fluxData[i] - fitMedianSlope) > dytop) {
                            dytop = fabs(fluxData[i] - fitMedianSlope);
                        }
                    }
                    
                    if(dytop == 0) {
                        dytop = abdevm;
                    }
                    
                    FluxCalibrationFitVector[order_count]->setflux(double(bm*(float)order + am + dytop),dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance(double(0.674433*abdevm)*double(0.674433*abdevm),dataIndex);
                } else {
                    FluxCalibrationFitVector[order_count]->setflux((double)orderFlux,dataIndex);
                    FluxCalibrationFitVector[order_count]->setvariance((double)fluxvariance,dataIndex);
                }
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam] && !(npBeam[beam]%2)) {
                        npBeam[beam]--;
                    }
                }
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam]>=3) {
                        float amBeam,bmBeam,abdevmBeam;
                        ladfit(orderBeamData[beam],fluxBeamData[beam],npBeam[beam],&amBeam,&bmBeam,&abdevmBeam); /* robust linear fit: f(x) =  a + b*x */
                        
                        float dytop = 0;
                        
                        for(unsigned i=0;i<npBeam[beam];i++){
                            double fitMedianSlope = (bmBeam*orderBeamData[beam][i] + amBeam);
                            
                            if(fabs(fluxBeamData[beam][i] - fitMedianSlope) < abdevmBeam &&
                               (fluxBeamData[beam][i] - fitMedianSlope) > dytop) {
                                dytop = fluxBeamData[beam][i] - fitMedianSlope;
                            }
                        }
                        
                        if(dytop == 0) {
                            dytop = abdevmBeam;
                        }
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(bmBeam*(float)order + amBeam + dytop),dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(0.674433*abdevmBeam)*double(0.674433*abdevmBeam),dataIndex);
                    } else {
                        BeamFluxCalibrationFitVector[beam][order_count]->setflux((double)orderBeamFlux[beam],dataIndex);
                        BeamFluxCalibrationFitVector[beam][order_count]->setvariance((double)beamfluxvariance,dataIndex);
                    }
                }
            } // for(unsigned dataIndex=0;dataIndex < maxNDataPoints; dataIndex++) {
			order_count++;
        } // if(GetSpectralOrder(order)->gethasSpectralElements() && GetSpectralOrder(order)->gethasSpectralEnergyDistribution()) {
    } // for(int order=(int)minorder; order<=(int)maxorder; order++) {
    
    order_count = 0;
                       
    for(int order=(int)actualMinorder; order<=(int)actualMaxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if(spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            unsigned nDataPointsOfCurrentOrder = spectralEnergyDistribution->getnDataPoints();
            
            for(unsigned dataIndex = 0; dataIndex < nDataPointsOfCurrentOrder; dataIndex++) {
                spectralEnergyDistribution->setfluxData(FluxCalibrationFitVector[order_count]->getflux(dataIndex),dataIndex);
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    beamSED->setfluxData(BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex),dataIndex);
                }
#ifdef PRINT_DEBUG
                cout << order << ' ' << dataIndex << ' '
                << spectralEnergyDistribution->getdistanceData(dataIndex) << ' '
                << spectralEnergyDistribution->getfluxData(dataIndex) << ' '
                << FluxCalibrationFitVector[order_count]->getflux(dataIndex) << ' ';
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    cout << beam << ' ' << beamSED->getfluxData(dataIndex) << ' ' << BeamFluxCalibrationFitVector[beam][order_count]->getflux(dataIndex) << ' ';
                }
                cout << endl;
#endif
            }
#ifdef PRINT_DEBUG
            cout << endl;
#endif
            order_count++;
        }
    }
    
    for(int order=(int)actualMinorder; order<=(int)actualMaxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            spectralEnergyDistribution->populateUncalibratedElementsFromContinuumData();
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                beamSED->populateUncalibratedElementsFromContinuumData();
            }
        }
    } //for(int order=(int)minorder; order<=(int)maxorder; order++) {

}

void operaSpectralOrderVector::FitFluxCalibrationAcrossOrders(int lowOrderToClip, int highOrderToClip, int orderBin, bool throughput) {
    
    int usefulMinorder = (int)minorder;
    int usefulMaxorder = (int)maxorder;
    
    bool hasUsefulMinorder = false;
    
    unsigned maxNElements = 0;
    unsigned numberOfBeams = 0;
    
    operaSpectralElements* fluxCalibrationElements = NULL;
    
    // Step 1. Figure out maximum number of elements and number of beams
    for(int order=(int)minorder; order<=(int)maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralEnergyDistribution() &&
            spectralOrder->gethasSpectralElements() &&
            spectralOrder->gethasWavelength()) {
            
            if(!hasUsefulMinorder) {
                usefulMinorder = order;
                hasUsefulMinorder = true;
            }
            
            usefulMaxorder = order;
            
            numberOfBeams = spectralOrder->getnumberOfBeams();
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            
            if(throughput == true) {
                fluxCalibrationElements = spectralEnergyDistribution->getThroughputElements();
            } else {
                fluxCalibrationElements = spectralEnergyDistribution->getFluxCalibrationElements();
            }
            if(maxNElements < fluxCalibrationElements->getnSpectralElements()) {
                maxNElements = fluxCalibrationElements->getnSpectralElements();
            }
        }
    }

    operaFluxVector *FluxCalibrationFitVector[MAXORDERS];
    float *fluxData = new float[MAXORDERS];
    float *orderData = new float[MAXORDERS];
    float *fluxBeamData[MAXNUMBEROFBEAMS];
    
    operaFluxVector *BeamFluxCalibrationFitVector[MAXNUMBEROFBEAMS][MAXORDERS];
    float *orderBeamData[MAXNUMBEROFBEAMS];
    float orderBeamFlux[MAXNUMBEROFBEAMS];
    unsigned npBeam[MAXNUMBEROFBEAMS];
    
    for(unsigned beam=0; beam < numberOfBeams; beam++) {
        npBeam[beam] = 0;
        fluxBeamData[beam] = new float[MAXORDERS];
        orderBeamData[beam] = new float[MAXORDERS];
    }
    
    float fluxvariance = 0;
    float beamfluxvariance = 0;
    
    unsigned order_count = 0;
    
    // Step 2. Calculate robust fit between neighbor orders
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        
        int loword = order - orderBin;
        int hiord = order + orderBin;
        if(loword < usefulMinorder) {
            loword = usefulMinorder;
        }
        if(hiord > usefulMaxorder) {
            hiord = usefulMaxorder;
        }
        
        int nord = (hiord - loword + 1);
        
        if(!(nord%2)) {
            if(hiord == usefulMaxorder) {
                loword--;
            } else {
                hiord++;
            }
        }
        
        nord = (hiord - loword + 1);
        
        unsigned nElementsOfCurrentOrder = GetSpectralOrder(order)->getSpectralElements()->getnSpectralElements();
        
        FluxCalibrationFitVector[order_count]  = new operaFluxVector(nElementsOfCurrentOrder);
        
        for(unsigned beam=0; beam < numberOfBeams; beam++) {
            BeamFluxCalibrationFitVector[beam][order_count] = new operaFluxVector(nElementsOfCurrentOrder);
        }

        for(unsigned elemIndex=0;elemIndex < maxNElements; elemIndex++) {
            
            if(elemIndex>=nElementsOfCurrentOrder){break;}
            
            unsigned np = 0;
            float orderFlux = NAN;
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                npBeam[beam] = 0;
                orderBeamFlux[beam] = NAN;
            }
            
            if (order > lowOrderToClip && order < highOrderToClip) {
                for(int o=loword; o<=hiord; o++) {
                    operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                    
                    if(spectralOrder->gethasSpectralElements() &&
                       spectralOrder->gethasSpectralEnergyDistribution() &&
                       spectralOrder->gethasWavelength()) {
                        
                        operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                        
                        if(throughput == true) {
                            fluxCalibrationElements =  spectralEnergyDistribution->getThroughputElements();
                        } else {
                            fluxCalibrationElements = spectralEnergyDistribution->getFluxCalibrationElements();
                        }
                        
                        if(elemIndex < fluxCalibrationElements->getnSpectralElements() && elemIndex < nElementsOfCurrentOrder) {
                            float flux = (float)fluxCalibrationElements->getFlux(elemIndex);
                            fluxvariance = (float)fluxCalibrationElements->getFluxVariance(elemIndex);
                            
                            if ((o > lowOrderToClip && o < highOrderToClip) || o==order) {
                                
                                if(!isnan(flux)) {
                                    orderData[np] = (float)o;
                                    fluxData[np] = flux;
                                    if(o==order) {
                                        orderFlux = flux;
                                    }
                                    np++;
                                }
                                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                                    operaSpectralElements* BeamFluxCalibrationElements;
                                    
                                    if(throughput == true) {
                                        BeamFluxCalibrationElements = beamSED->getThroughputElements();
                                    } else {
                                        BeamFluxCalibrationElements = beamSED->getFluxCalibrationElements();
                                    }
                                    
                                    float beamflux = (float)BeamFluxCalibrationElements->getFlux(elemIndex);
                                    beamfluxvariance = (float)BeamFluxCalibrationElements->getFluxVariance(elemIndex);
                                    
                                    if(!isnan(beamflux)) {
                                        orderBeamData[beam][npBeam[beam]] = (float)o;
                                        fluxBeamData[beam][npBeam[beam]] = beamflux;
                                        if(o==order) {
                                            orderBeamFlux[beam] = beamflux;
                                        }
                                        npBeam[beam]++;
                                    }
                                }
                            }
                        }
                        
                    }
                }
                if(np && !(np%2)) {
                    np--;
                }
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    if(npBeam[beam] && !(npBeam[beam]%2)) {
                        npBeam[beam]--;
                    }
                }
            } else {
                int o = order;
                operaSpectralOrder *spectralOrder = GetSpectralOrder(o);
                
                if(spectralOrder->gethasSpectralElements() &&
                   spectralOrder->gethasSpectralEnergyDistribution() &&
                   spectralOrder->gethasWavelength()) {
                    
                    operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                    
                    if(throughput == true) {
                        fluxCalibrationElements =  spectralEnergyDistribution->getThroughputElements();
                    } else {
                        fluxCalibrationElements = spectralEnergyDistribution->getFluxCalibrationElements();
                    }
                    
                    if(elemIndex < fluxCalibrationElements->getnSpectralElements() && elemIndex < nElementsOfCurrentOrder) {
                        float flux = (float)fluxCalibrationElements->getFlux(elemIndex);
                        fluxvariance = (float)fluxCalibrationElements->getFluxVariance(elemIndex);
                        
                        if(!isnan(flux)) {
                            orderData[np] = (float)o;
                            fluxData[np] = flux;
                            orderFlux = flux;
                            np++;
                        }
                        for(unsigned beam=0; beam < numberOfBeams; beam++) {
                            operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                            operaSpectralElements* BeamFluxCalibrationElements;
                            if(throughput == true) {
                                BeamFluxCalibrationElements =  beamSED->getThroughputElements();
                            } else {
                                BeamFluxCalibrationElements = beamSED->getFluxCalibrationElements();
                            }
                            float beamflux = (float)BeamFluxCalibrationElements->getFlux(elemIndex);
                            beamfluxvariance = (float)BeamFluxCalibrationElements->getFluxVariance(elemIndex);
                            if(!isnan(beamflux)) {
                                orderBeamData[beam][npBeam[beam]] = (float)o;
                                fluxBeamData[beam][npBeam[beam]] = beamflux;
                                orderBeamFlux[beam] = beamflux;
                                npBeam[beam]++;
                            }
                        }
                    }
                }
            }
            
            if(np>=3) {
                float am,bm,abdevm;
                ladfit(orderData,fluxData,np,&am,&bm,&abdevm); // robust linear fit: f(x) =  a + b*x 
                
                FluxCalibrationFitVector[order_count]->setflux(double(bm*(double)order + am),elemIndex);
                FluxCalibrationFitVector[order_count]->setvariance(double(0.674433*abdevm)*double(0.674433*abdevm),elemIndex);
            } else {
                FluxCalibrationFitVector[order_count]->setflux((double)orderFlux,elemIndex);
                FluxCalibrationFitVector[order_count]->setvariance((double)fluxvariance,elemIndex);
            }
            
            for(unsigned beam=0; beam < numberOfBeams; beam++) {
                if(npBeam[beam]>=3) {
                    float amBeam,bmBeam,abdevmBeam;
                    ladfit(orderBeamData[beam],fluxBeamData[beam],npBeam[beam],&amBeam,&bmBeam,&abdevmBeam); // robust linear fit: f(x) =  a + b*x 
                    BeamFluxCalibrationFitVector[beam][order_count]->setflux(double(bmBeam*(double)order + amBeam),elemIndex);
                    BeamFluxCalibrationFitVector[beam][order_count]->setvariance(double(0.674433*abdevmBeam)*double(0.674433*abdevmBeam),elemIndex);
                } else {
                    BeamFluxCalibrationFitVector[beam][order_count]->setflux((double)orderBeamFlux[beam],elemIndex);
                    BeamFluxCalibrationFitVector[beam][order_count]->setvariance((double)beamfluxvariance,elemIndex);
                }
            }
        } //  for(unsigned elemIndex=0;elemIndex < maxNElements; elemIndex++) {
        order_count++;
    } // for(int order=(int)minorder; order<=(int)maxorder; order++) {
    
    order_count = 0;
    // Step 3. Feed back fit quantities to order SED
    for(int order=usefulMinorder; order<=usefulMaxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if(spectralOrder->gethasSpectralElements() &&
           spectralOrder->gethasSpectralEnergyDistribution() && 
           spectralOrder->gethasWavelength()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            
            if(throughput == true) {
                fluxCalibrationElements = spectralEnergyDistribution->getThroughputElements();
            } else {
                fluxCalibrationElements = spectralEnergyDistribution->getFluxCalibrationElements();
            }
            
            for(unsigned elemIndex=0;elemIndex < fluxCalibrationElements->getnSpectralElements(); elemIndex++) {
                
                fluxCalibrationElements->setFlux(FluxCalibrationFitVector[order_count]->getflux(elemIndex),elemIndex);
                fluxCalibrationElements->setFluxVariance(FluxCalibrationFitVector[order_count]->getvariance(elemIndex),elemIndex);
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    operaSpectralElements* BeamFluxCalibrationElements;
                    if(throughput == true) {
                        BeamFluxCalibrationElements = beamSED->getThroughputElements();
                    } else {
                        BeamFluxCalibrationElements = beamSED->getFluxCalibrationElements();
                    }
                    BeamFluxCalibrationElements->setFlux(BeamFluxCalibrationFitVector[beam][order_count]->getflux(elemIndex),elemIndex);
                    BeamFluxCalibrationElements->setFluxVariance(BeamFluxCalibrationFitVector[beam][order_count]->getvariance(elemIndex),elemIndex);
                }
            }
            order_count++;
        }
    }
}

void operaSpectralOrderVector::getContinuumFluxesForNormalization(double *uncalibratedContinuumFluxForNormalization, double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS],unsigned binsize, int orderBin, unsigned nsigcut) {
    
    /*
     * Note: This function returns the continuum flux measured for a given wavelength obtained from 
     *       spectral energy distribution class.  If the "wavelengthfornormalization" given in the SED  
     *       is not in the range covevered by the orders, then it returns NaNs.
     *
     */
    
    unsigned NumberofBeams = 0;
    double wavelengthForNormalization=0;
    
    for (int order=(int)minorder; order<=(int)maxorder; order++) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasSpectralEnergyDistribution() &&
            spectralOrder->gethasSpectralElements() &&
            spectralOrder->gethasWavelength()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            wavelengthForNormalization = spectralEnergyDistribution->getwavelengthForNormalization();
            NumberofBeams = spectralOrder->getnumberOfBeams();
            break;
        }
    }
    
    if(wavelengthForNormalization == 0 || NumberofBeams == 0) {
        //throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		return;	// ths case of a constant DT Aug 9 2013
	}
    
    int *orderWithReferenceFluxForNormalization = new int[MAXORDERS];
    unsigned *elemIndexWithReferenceFluxForNormalization = new unsigned[MAXORDERS];
    
    unsigned nOrdersPicked = getElemIndexAndOrdersByWavelength(orderWithReferenceFluxForNormalization,elemIndexWithReferenceFluxForNormalization,wavelengthForNormalization);
    
    measureContinuumAcrossOrders(binsize,orderBin,nsigcut,nOrdersPicked,orderWithReferenceFluxForNormalization);
    
    // DT May 20 2014 -- notused --int orderpicked = 0;
    unsigned elemIndexPicked = 0;
    
    *uncalibratedContinuumFluxForNormalization = -BIG;
    bool hasfluxforNormalization = false;
    
    for(unsigned i=0;i<nOrdersPicked;i++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(orderWithReferenceFluxForNormalization[i]);
        
        if (spectralOrder->gethasSpectralEnergyDistribution()) {
            
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            double tmp_uncalibratedContinuumFluxForNormalization = spectralEnergyDistribution->getUncalibratedFluxElements()->getFlux(elemIndexWithReferenceFluxForNormalization[i]);
            
            if(*uncalibratedContinuumFluxForNormalization < tmp_uncalibratedContinuumFluxForNormalization) {
                *uncalibratedContinuumFluxForNormalization = tmp_uncalibratedContinuumFluxForNormalization;
                //orderpicked = orderWithReferenceFluxForNormalization[i];
                elemIndexPicked = elemIndexWithReferenceFluxForNormalization[i];
                for(unsigned beam=0; beam < NumberofBeams; beam++) {
                    operaSpectralEnergyDistribution *beamSED = spectralOrder->getBeamSED(beam);
                    uncalibratedContinuumBeamFluxForNormalization[beam] = beamSED->getUncalibratedFluxElements()->getFlux(elemIndexPicked);
                }
                hasfluxforNormalization = true;
            }
        }
    }
    
    if(!hasfluxforNormalization) {
        *uncalibratedContinuumFluxForNormalization = NAN;
        for(unsigned beam=0; beam < NumberofBeams; beam++) {
            uncalibratedContinuumBeamFluxForNormalization[beam] = NAN;
        }
    }
    
    delete[] orderWithReferenceFluxForNormalization;
    delete[] elemIndexWithReferenceFluxForNormalization;
}

unsigned operaSpectralOrderVector::getLEElementCount(string LEfluxCalibration) {
	ifstream astream;
	string dataline;
	unsigned length = 0;
    
	astream.open(LEfluxCalibration.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
                    length++;  
                }
            }
		} // while (astream.good())
		astream.close();
	}	// if (astream.open()
	return length;
}

void operaSpectralOrderVector::readLEFluxCalibration(string LEfluxCalibration, operaSpectralElements *fluxCalibrationElements) {
	ifstream astream;
	string dataline;
	unsigned length = 0;
	double wl, fcal;
    
	astream.open(LEfluxCalibration.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &wl, &fcal);
					fluxCalibrationElements->setwavelength(wl, length);
					fluxCalibrationElements->setFlux(fcal, length);
					fluxCalibrationElements->setHasWavelength(true);
                    length++;  
                }
            }
		} // while (astream.good())
		astream.close();
	}	// if (astream.open()
}

void operaSpectralOrderVector::calculateCleanUniformSampleOfContinuum(int Minorder, int Maxorder, unsigned binsize, double delta_wl, string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, float *uniform_wl, float *uniform_flux,float *uniform_Beamflux[MAXNUMBEROFBEAMS], bool useBeams) {
    bool debug = false;
    /*
     * Notes: This function returns a clean uniform sample of the continuum flux for
     *        both the main elements and the beam elements of all orders.
     *        The spectral regions for the continuum is obtained from an input mask
     *        inputWavelengthMaskForUncalContinuum.
     *        The wavelength vector of the sample is output to uniform_wl, 
     *        The flux vectors are output to uniform_flux and uniform_Beamflux[beam]
     *        The size of final sample is given by numberOfPointsInUniformSample
     *        Binsize defines how many points should be binned to get rid of noise/features
     *        delta_wl defines which is the minimum wavelength range to merge data
     *        This function uses a number of routines in the operaSpectralTools library
     */
    
    //---------------------------------
    // Loop over orders to set maximum number of elements and number of beams
    // --> maxNElements & NumberofBeams
    unsigned NumberofBeams = getNumberofBeams(Minorder, Maxorder);
    unsigned maxNElements = getMaxNumberOfElementsInOrder(Minorder, Maxorder);

    if(!useBeams){
        NumberofBeams = 0;
    }
    
   /* It's commented cause I'm not sure this would work. 
      It would be a way to create memory for input NULL vectors
    
    if(uniform_wl==NULL) {
        uniform_wl = new float[numberOfPointsInUniformSample];
    }
    if(uniform_flux==NULL){
        uniform_flux = new float[numberOfPointsInUniformSample];
    }
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        if(uniform_Beamflux[beam]==NULL) {
            uniform_Beamflux[beam] = new float[numberOfPointsInUniformSample];
        }
    }
    */
    //---------------------------------
    // Calculate a clean sample of the continuum from the uncalibrated spectrum
    unsigned maxNumberoOfTotalPoints = (unsigned)(ceil((float)maxNElements/(float)binsize) + 1)*MAXORDERS;

    //---------------------------------
    // Collect continuum sample using input mask
    
    float *uncal_wl = new float[maxNumberoOfTotalPoints];
    float *uncal_flux = new float[maxNumberoOfTotalPoints];
    float *uncal_Beamflux[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        uncal_Beamflux[beam] = new float[maxNumberoOfTotalPoints];
    }
    
    double *wl0_vector = new double[MAXNUMBEROFWLRANGES];
    double *wlf_vector = new double[MAXNUMBEROFWLRANGES];
    
    unsigned nRangesInWLMask = readContinuumWavelengthMask(inputWavelengthMaskForUncalContinuum,wl0_vector,wlf_vector);

    float *wl_tmp = new float[2*binsize];
    float *flux_tmp = new float[2*binsize];
    float *beamflux_tmp[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        beamflux_tmp[beam] = new float[2*binsize];
    }
    
    unsigned nTotalPoints = 0;
    double minwl = BIG;
    double maxwl = -BIG;
    
    for (int order=Minorder; order<=Maxorder; order++) {

        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            
            if(debug) {
                for(unsigned i=0; i<spectralElements->getnSpectralElements(); i++) {
                    cout << order << " "
                    << spectralElements->getwavelength(i) << " "
                    << spectralElements->getFlux(i) << endl;
                }
            }
            unsigned NumberOfElementSamples = (unsigned)ceil((float)spectralElements->getnSpectralElements()/(float)binsize);
            
            for(unsigned k=0;k<NumberOfElementSamples;k++){
                unsigned firstElement = binsize*(k);
                unsigned lastElement =  binsize*(k+1);
                
                if (lastElement > spectralElements->getnSpectralElements()){
                    lastElement = spectralElements->getnSpectralElements();
                    if(binsize > lastElement) {
                        firstElement = 0;
                    } else {
                        firstElement = lastElement - binsize;
                    }
                }
                
                unsigned np=0;
                for(unsigned elemIndex=firstElement;elemIndex < lastElement; elemIndex++) {
                    
                    for(unsigned rangeElem=0; rangeElem<nRangesInWLMask; rangeElem++) {
                        if(spectralElements->getwavelength(elemIndex) >= wl0_vector[rangeElem] &&
                           spectralElements->getwavelength(elemIndex) <= wlf_vector[rangeElem] ) {
                            
                            if(spectralElements->getFlux(elemIndex) && !isnan((float)spectralElements->getFlux(elemIndex))) {
                                
                                wl_tmp[np] = spectralElements->getwavelength(elemIndex);
                                if(minwl > spectralElements->getwavelength(elemIndex)) {
                                    minwl = spectralElements->getwavelength(elemIndex);
                                }
                                if(maxwl < spectralElements->getwavelength(elemIndex)) {
                                    maxwl = spectralElements->getwavelength(elemIndex);
                                }
                                
                                flux_tmp[np] = spectralElements->getFlux(elemIndex);
                                
                                //cout << wl_tmp[np] << " " << flux_tmp[np] << endl;
                                
                                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                                    beamflux_tmp[beam][np] = spectralOrder->getBeamElements(beam)->getFlux(elemIndex);
                                }
                                np++;
                                break;
                            }
                        }
                    }
                }
                
                if(np > MINNUMBEROFPOINTSINSIDEBIN) {
                    uncal_wl[nTotalPoints] = operaArrayMean(np,wl_tmp);
                    uncal_flux[nTotalPoints] = operaArrayMedian(np,flux_tmp);
                    for(unsigned beam=0;beam<NumberofBeams;beam++) {
                        uncal_Beamflux[beam][nTotalPoints] = operaArrayMedian(np,beamflux_tmp[beam]);
                    }
                    
                    nTotalPoints++;
                }
            }
        }
    }

    //---------------------------------
    // Sort sample
    int *sindex = new int[nTotalPoints];
    operaArrayIndexSort((int)nTotalPoints,uncal_wl,sindex);
   // cout << "minwl=" << minwl << " maxwl=" << maxwl << endl;
   // cout << uncal_wl[sindex[0]] << " " << uncal_wl[sindex[nTotalPoints-1]] << endl;
    
   /* for(unsigned index=0; index<nTotalPoints; index++) {
        cout << uncal_wl[sindex[index]] << " " << uncal_flux[sindex[index]] << endl;
    }*/
    
    //---------------------------------
    // Merge data from all orders
    
    float *uncal_final_wl = new float[maxNumberoOfTotalPoints+2];
    float *uncal_final_flux = new float[maxNumberoOfTotalPoints+2];
    float *uncal_final_Beamflux[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        uncal_final_Beamflux[beam] = new float[maxNumberoOfTotalPoints+2];
    }

    unsigned npFinal = 0;
    unsigned np = 0;
    for(unsigned index=0; index<nTotalPoints; index++) {
        wl_tmp[np] = uncal_wl[sindex[index]];
        flux_tmp[np] = uncal_flux[sindex[index]];
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            beamflux_tmp[beam][np] = uncal_Beamflux[beam][sindex[index]];
        }
        np++;
        
        if(index) {
            if(fabs(uncal_wl[sindex[index]] - uncal_wl[sindex[index-1]]) > delta_wl) {
                np--;
                
                if(np==1) {
                    uncal_final_wl[npFinal] = uncal_wl[sindex[index-1]];
                    uncal_final_flux[npFinal] = uncal_flux[sindex[index-1]];
                    for(unsigned beam=0;beam<NumberofBeams;beam++) {
                        uncal_final_Beamflux[beam][npFinal] = uncal_Beamflux[beam][sindex[index-1]];
                    }
                } else {
                    uncal_final_wl[npFinal] = operaArrayMean(np,wl_tmp);
                    uncal_final_flux[npFinal] = operaArrayMean(np,flux_tmp);
                    for(unsigned beam=0;beam<NumberofBeams;beam++) {
                        uncal_final_Beamflux[beam][npFinal] = operaArrayMean(np,beamflux_tmp[beam]);
                    }
                }

                npFinal++;
                
                np = 0;
                wl_tmp[np] = uncal_wl[sindex[index]];
                flux_tmp[np] = uncal_flux[sindex[index]];
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    beamflux_tmp[beam][np] = uncal_Beamflux[beam][sindex[index]];
                }
                np++;
            }
            
            if(index == nTotalPoints-1) {
                if(np==1) {
                    uncal_final_wl[npFinal] = uncal_wl[sindex[index-1]];
                    uncal_final_flux[npFinal] = uncal_flux[sindex[index-1]];
                    for(unsigned beam=0;beam<NumberofBeams;beam++) {
                        uncal_final_Beamflux[beam][npFinal] = uncal_Beamflux[beam][sindex[index-1]];
                    }
                } else {
                    uncal_final_wl[npFinal] = operaArrayMean(np,wl_tmp);
                    uncal_final_flux[npFinal] = operaArrayMean(np,flux_tmp);
                    for(unsigned beam=0;beam<NumberofBeams;beam++) {
                        uncal_final_Beamflux[beam][npFinal] = operaArrayMean(np,beamflux_tmp[beam]);
                    }
                }
                npFinal++;
            }
            
        }
        
        if(npFinal==1) {
            uncal_final_wl[npFinal] = uncal_final_wl[npFinal-1];
            uncal_final_flux[npFinal] = uncal_final_flux[npFinal-1];
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                uncal_final_Beamflux[beam][npFinal] = uncal_final_Beamflux[beam][npFinal-1];
            }
            uncal_final_wl[npFinal-1] = minwl;
            npFinal++;
        }

        //cout << uncal_final_wl[npFinal-1] << " " << uncal_final_flux[npFinal-1] << endl;
    }

    uncal_final_wl[npFinal] = maxwl;
    if(npFinal) {
        uncal_final_flux[npFinal] = uncal_final_flux[npFinal-1];
        
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            uncal_final_Beamflux[beam][npFinal] = uncal_final_Beamflux[beam][npFinal-1];
        }
    }

    npFinal++;

 /*   for (unsigned i=0; i<npFinal; i++) {
        cout << uncal_final_wl[i] << " " << uncal_final_flux[i] << endl;
    }
   */     
    //---------------------------------
    // Calculate an uniform sample.
    // This is a necessary step in order to make further spline interpolations to work.

    calculateUniformSample(npFinal,uncal_final_wl,uncal_final_flux,numberOfPointsInUniformSample,uniform_wl,uniform_flux);
    
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        calculateUniformSample(npFinal,uncal_final_wl,uncal_final_Beamflux[beam],numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam]);
    }
    
    if(debug) {
        // original sample
        for (unsigned i=0; i<npFinal; i++) {
            cout << i << " " << uncal_final_wl[i] << " " << uncal_final_flux[i] << " ";
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                cout << uncal_final_Beamflux[beam][i] << " ";
            }
            cout << endl;
        }
        // uniform sample
        for (unsigned i=0; i<numberOfPointsInUniformSample; i++) {
            cout << uniform_wl[i] << " " << uniform_flux[i] << " ";
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                cout << uniform_Beamflux[beam][i] << " ";
            }
            cout << endl;
        }
    }

    delete[] uncal_wl;
    delete[] uncal_flux;
    delete[] wl0_vector;
    delete[] wlf_vector;
    delete[] wl_tmp;
    delete[] flux_tmp;
    delete[] uncal_final_wl;
    delete[] uncal_final_flux;
    
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        delete[] uncal_Beamflux[beam];
        delete[] beamflux_tmp[beam];
        delete[] uncal_final_Beamflux[beam];
    }
    
    delete[] sindex;
}

// Load telluric corrected wavelength calibration into extended spectra
void operaSpectralOrderVector::readTelluricWavelengthINTOExtendendSpectra(string telluriccorrection, int Minorder, int Maxorder) {
    if (telluriccorrection.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    ReadSpectralOrders(telluriccorrection);
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            operaWavelength *wavelength = spectralOrder->getWavelength();
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralElements->setwavelengthsFromCalibration(wavelength);
            spectralElements->copyTOtell();	// Save the tell
            spectralElements->setHasWavelength(true);
        }
    }
}

void operaSpectralOrderVector::readRVCorrectionINTOExtendendSpectra(string Radialvelocitycorrection, string WavelengthCalibration, int Minorder, int Maxorder) {
    if (Radialvelocitycorrection.empty() || WavelengthCalibration.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    readRadialVelocityCorrection(Radialvelocitycorrection);

    ReadSpectralOrders(WavelengthCalibration);

    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            operaWavelength *wavelength = spectralOrder->getWavelength();
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralElements->setwavelengthsFromCalibration(wavelength);
            spectralOrder->setExtendedBarycentricWavelengthCorrection(getBarycentricRadialVelocityCorrection());
            spectralElements->setHasWavelength(true);
        }
    }
}

void operaSpectralOrderVector::correctFlatField(string inputFlatFluxCalibration, int Minorder, int Maxorder, bool StarPlusSky) {
    if (inputFlatFluxCalibration.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }

    unsigned initialNumberofBeams = getNumberofBeams(Minorder, Maxorder);        
    ReadSpectralOrders(inputFlatFluxCalibration);
    unsigned flatCalNumberofBeams = getNumberofBeams(Minorder, Maxorder);
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        //Below it ignores all Beams in case when the number of beams are different
        if(initialNumberofBeams != flatCalNumberofBeams) {
            spectralOrder->setnumberOfBeams(0);
        }
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            spectralOrder->divideSpectralElementsBySEDElements(true, NULL,StarPlusSky,false);
        } else {
            spectralOrder->sethasSpectralElements(false);
        }
    }
}

void operaSpectralOrderVector::correctFlatField(string inputFlatFluxCalibration, int Minorder, int Maxorder, bool StarPlusSky, bool starplusskyInvertSkyFiber) {
    if (inputFlatFluxCalibration.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned initialNumberofBeams = getNumberofBeams(Minorder, Maxorder);    
    ReadSpectralOrders(inputFlatFluxCalibration);
    unsigned flatCalNumberofBeams = getNumberofBeams(Minorder, Maxorder);
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        //Below it ignores all Beams in case when the number of beams are different
        if(initialNumberofBeams != flatCalNumberofBeams) {
            spectralOrder->setnumberOfBeams(0);
        }
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            spectralOrder->divideSpectralElementsBySEDElements(true, NULL,StarPlusSky,starplusskyInvertSkyFiber);
        } else {
            spectralOrder->sethasSpectralElements(false);
        }
    }
}

void operaSpectralOrderVector::saveExtendedRawFlux(int Minorder, int Maxorder) {
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            spectralOrder->getSpectralElements()->copyTOnormalizedFlux();
            spectralOrder->getSpectralElements()->copyTOfcalFlux();
        }
    }
}

void operaSpectralOrderVector::normalizeFluxINTOExtendendSpectra(string inputWavelengthMaskForUncalContinuum, unsigned numberOfPointsInUniformSample, unsigned normalizationBinsize, double delta_wl, int Minorder, int Maxorder, bool normalizeBeams) {
    if (inputWavelengthMaskForUncalContinuum.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned NumberofBeams = getNumberofBeams(Minorder, Maxorder);
    unsigned maxNElements = getMaxNumberOfElementsInOrder(Minorder, Maxorder);
    if(!normalizeBeams) {
        NumberofBeams = 0;
    }
    float *uniform_wl = new float[numberOfPointsInUniformSample];
    float *uniform_flux = new float[numberOfPointsInUniformSample];
    float *uniform_Beamflux[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        uniform_Beamflux[beam] = new float[numberOfPointsInUniformSample];
    }

    calculateCleanUniformSampleOfContinuum(Minorder,Maxorder,normalizationBinsize,delta_wl,inputWavelengthMaskForUncalContinuum,numberOfPointsInUniformSample,uniform_wl,uniform_flux,uniform_Beamflux, normalizeBeams);
    
    float *UncalibratedModelFlux = new float[maxNElements];
    float *UncalibratedModelBeamFlux[MAXNUMBEROFBEAMS];
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            UncalibratedModelBeamFlux[beam] = new float[maxNElements];
        }
    float *elemWavelength = new float[maxNElements];
    operaSpectralEnergyDistribution *BeamSED[MAXNUMBEROFBEAMS];
    operaSpectralElements *uncalibratedBeamFluxElements[MAXNUMBEROFBEAMS];
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() &&
            spectralOrder->gethasSpectralElements() &&
            spectralOrder->gethasSpectralEnergyDistribution()) {
            
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            
            unsigned nElements = spectralElements->getnSpectralElements();
            
            spectralEnergyDistribution->setUncalibratedFluxElements(spectralElements);
            
            operaSpectralElements *uncalibratedFluxElements = spectralEnergyDistribution->getUncalibratedFluxElements();
            
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    BeamSED[beam] = spectralOrder->getBeamSED(beam);
                    BeamSED[beam]->setUncalibratedFluxElements(spectralOrder->getBeamElements(beam));
                    uncalibratedBeamFluxElements[beam] = BeamSED[beam]->getUncalibratedFluxElements();
                }
            for(unsigned i=0;i<nElements;i++) {
                elemWavelength[i] =  (float)spectralElements->getwavelength(i);
            }

            operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_flux,nElements,elemWavelength,UncalibratedModelFlux);
            
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam],nElements,elemWavelength,UncalibratedModelBeamFlux[beam]);
                }
            for(unsigned i=0;i<nElements;i++) {
                uncalibratedFluxElements->setFlux((double)UncalibratedModelFlux[i],i);
                    for(unsigned beam=0;beam<NumberofBeams;beam++) {
                        // uncalibratedBeamFluxElements[beam]->setFlux((double)UncalibratedModelBeamFlux[beam][i],i);
                        uncalibratedBeamFluxElements[beam]->setFlux(1.0,i);
                    }
            }
            spectralEnergyDistribution->setHasUncalibratedFlux(true);
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    BeamSED[beam]->setHasUncalibratedFlux(true);
                }
            spectralOrder->applyNormalizationFromExistingContinuum(NULL,NULL,TRUE,normalizeBeams,2);

            spectralElements->copyTOnormalizedFlux();
            spectralElements->copyFROMrawFlux();
        }
    }
}

void operaSpectralOrderVector::normalizeAndCalibrateFluxINTOExtendendSpectra(string inputWavelengthMaskForUncalContinuum,string fluxCalibration, double exposureTime, bool AbsoluteCalibration, unsigned numberOfPointsInUniformSample, unsigned normalizationBinsize, double delta_wl, int Minorder, int Maxorder, bool normalizeBeams, bool StarPlusSky) {
    if (inputWavelengthMaskForUncalContinuum.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    if (fluxCalibration.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
   
    ReadSpectralOrders(fluxCalibration);
    
    unsigned NumberofBeams = getNumberofBeams(Minorder, Maxorder);
    
    if(StarPlusSky) {
        for (int order=Minorder; order<=Maxorder; order++) {
            operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
            if (spectralOrder->gethasWavelength() &&
                spectralOrder->gethasSpectralElements() &&
                spectralOrder->gethasSpectralEnergyDistribution()) {
                
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                unsigned nspecElem = spectralEnergyDistribution->getFluxCalibrationElements()->getnSpectralElements();
                operaFluxVector fluxCalibrationFluxVector(nspecElem);
                fluxCalibrationFluxVector = (*spectralEnergyDistribution->getFluxCalibrationElements()->getFluxVector());
                operaFluxVector thruputFluxVector(nspecElem);
                thruputFluxVector = (*spectralEnergyDistribution->getFluxCalibrationElements()->getFluxVector());
                
                for(unsigned beam=0; beam < NumberofBeams; beam++) {
                    if(float(beam) >= float(NumberofBeams)/2) {// Sky Fiber
                        operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
                        operaSpectralElements *fluxCalibrationBeamElements = BeamSED->getFluxCalibrationElements();
                        operaSpectralElements *thruputBeamElements = BeamSED->getThroughputElements();
                        
                        fluxCalibrationBeamElements->setFluxVector(&fluxCalibrationFluxVector);
                        thruputBeamElements->setFluxVector(&thruputFluxVector);
                    }
                }
            }
        }
    }
    
    if(StarPlusSky) {
        normalizeBeams = false;
    }
    
    if(!normalizeBeams) {
        NumberofBeams = 0;
    }
    
    unsigned maxNElements = getMaxNumberOfElementsInOrder(Minorder, Maxorder);

    float *uniform_wl = new float[numberOfPointsInUniformSample];
    float *uniform_flux = new float[numberOfPointsInUniformSample];
    float *uniform_Beamflux[MAXNUMBEROFBEAMS];

    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        uniform_Beamflux[beam] = new float[numberOfPointsInUniformSample];
    }

    calculateCleanUniformSampleOfContinuum(Minorder,Maxorder,normalizationBinsize,delta_wl,inputWavelengthMaskForUncalContinuum,numberOfPointsInUniformSample,uniform_wl,uniform_flux,uniform_Beamflux, normalizeBeams);
    /*
    for(unsigned i=0;i<numberOfPointsInUniformSample;i++) {
        cout << uniform_wl[i] << " " << uniform_flux[i] << endl;
    }
    exit(1);
    */
    float *UncalibratedModelFlux = new float[maxNElements];
    float *UncalibratedModelBeamFlux[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        UncalibratedModelBeamFlux[beam] = new float[maxNElements];
    }

    float *elemWavelength = new float[maxNElements];
    operaSpectralEnergyDistribution *BeamSED[MAXNUMBEROFBEAMS];
    operaSpectralElements *uncalibratedBeamFluxElements[MAXNUMBEROFBEAMS];

    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() &&
            spectralOrder->gethasSpectralElements() &&
            spectralOrder->gethasSpectralEnergyDistribution()) {
            
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            
            unsigned nElements = spectralElements->getnSpectralElements();
            
            spectralEnergyDistribution->setUncalibratedFluxElements(spectralElements);
            
            operaSpectralElements *uncalibratedFluxElements = spectralEnergyDistribution->getUncalibratedFluxElements();
            
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                BeamSED[beam] = spectralOrder->getBeamSED(beam);
                BeamSED[beam]->setUncalibratedFluxElements(spectralOrder->getBeamElements(beam));
                uncalibratedBeamFluxElements[beam] = BeamSED[beam]->getUncalibratedFluxElements();
            }
            for(unsigned i=0;i<nElements;i++) {
                elemWavelength[i] =  (float)spectralElements->getwavelength(i);
                
            }
            
            operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_flux,nElements,elemWavelength,UncalibratedModelFlux);
            
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam],nElements,elemWavelength,UncalibratedModelBeamFlux[beam]);
            }

            for(unsigned i=0;i<nElements;i++) {
                uncalibratedFluxElements->setFlux((double)UncalibratedModelFlux[i],i);
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    // uncalibratedBeamFluxElements[beam]->setFlux((double)UncalibratedModelBeamFlux[beam][i],i);
                    uncalibratedBeamFluxElements[beam]->setFlux(1.0,i);
                }
            }
            spectralEnergyDistribution->setHasUncalibratedFlux(true);
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                BeamSED[beam]->setHasUncalibratedFlux(true);
            }
            spectralOrder->applyNormalizationFromExistingContinuum(NULL,NULL,TRUE,normalizeBeams,2);
            
            spectralElements->copyTOnormalizedFlux();
            spectralElements->copyFROMfcalFlux();
        }
    }
    
    double wavelengthForNormalization=0;
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            if(spectralEnergyDistribution->getwavelengthForNormalization()>0) {
                wavelengthForNormalization = spectralEnergyDistribution->getwavelengthForNormalization();
                break;
            }
        }
    }
    
    double spectralBinConstant = 1.0;
    double BeamSpectralBinConstant[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        BeamSpectralBinConstant[beam] = 1.0;
    }

    //-- Calculate Sky over Star fiber area ratio to compensate for different apertures.
    double SkyOverStarFiberAreaRatio = (2.2*2.2)/(1.6*1.6);
    
    if(AbsoluteCalibration) {
        spectralBinConstant = exposureTime;
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            if(StarPlusSky && float(beam) >= float(NumberofBeams)/2) {
                // Sky Fiber
                BeamSpectralBinConstant[beam] = exposureTime/SkyOverStarFiberAreaRatio;
            } else {
                BeamSpectralBinConstant[beam] = exposureTime; 
            }
        }
    } else {
        spectralBinConstant = (double)getFluxAtWavelength(numberOfPointsInUniformSample,uniform_wl,uniform_flux,wavelengthForNormalization);
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            BeamSpectralBinConstant[beam] = (double)getFluxAtWavelength(numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam],wavelengthForNormalization);
            if(StarPlusSky && float(beam) >= float(NumberofBeams)/2) {
                // Sky Fiber
                BeamSpectralBinConstant[beam] /= SkyOverStarFiberAreaRatio;
            }
        }
    }

    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            
            if(AbsoluteCalibration) { // use fluxcalibation
                spectralOrder->multiplySpectralElementsBySEDElements(false, spectralBinConstant, NULL);
            } else { // use throughput
                spectralOrder->multiplySpectralElementsBySEDElements(true, spectralBinConstant, NULL);
            }
            spectralElements->copyTOfcalFlux();
            spectralElements->copyFROMrawFlux();
        }
    }
    
    delete[] uniform_wl;
    delete[] uniform_flux;
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        delete[] uniform_Beamflux[beam];
    }
    delete[] UncalibratedModelFlux;
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        delete[] UncalibratedModelBeamFlux[beam];
    }
    delete[] elemWavelength;
}



/*
 *  The function below normalizes the spectrum to the continuum AND 
 *  applies a flat response flux calibration. Both spectra are saved separately into
 *  the extended spectra.
 */

void operaSpectralOrderVector::normalizeAndApplyFlatResponseINTOExtendendSpectra(string inputWavelengthMaskForUncalContinuum, string flatResponse, unsigned numberOfPointsInUniformSample, unsigned normalizationBinsize, double delta_wl, int Minorder, int Maxorder, bool normalizeBeams, bool StarPlusSky) {

    if (inputWavelengthMaskForUncalContinuum.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    if (flatResponse.empty()) {
        throw operaException("operaSpectralOrderVector: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
    }
    
    // Function below reads flat response spectrum, interpolate values and
    // feed them into spectralEnergyDistribution for all orders
    readLibreEspritFlatResponseIntoSED(flatResponse,Minorder,Maxorder);
    
    unsigned NumberofBeams = getNumberofBeams(Minorder, Maxorder);
    
    if(StarPlusSky) {
        for (int order=Minorder; order<=Maxorder; order++) {
            operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
            if (spectralOrder->gethasWavelength() &&
                spectralOrder->gethasSpectralElements() &&
                spectralOrder->gethasSpectralEnergyDistribution()) {
                
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                unsigned nspecElem = spectralEnergyDistribution->getFluxCalibrationElements()->getnSpectralElements();
                operaFluxVector fluxCalibrationFluxVector(nspecElem);
                fluxCalibrationFluxVector = (*spectralEnergyDistribution->getFluxCalibrationElements()->getFluxVector());
                operaFluxVector thruputFluxVector(nspecElem);
                thruputFluxVector = (*spectralEnergyDistribution->getFluxCalibrationElements()->getFluxVector());
                
                for(unsigned beam=0; beam < NumberofBeams; beam++) {
                    if(float(beam) >= float(NumberofBeams)/2) {// Sky Fiber
                        operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
                        operaSpectralElements *fluxCalibrationBeamElements = BeamSED->getFluxCalibrationElements();
                        operaSpectralElements *thruputBeamElements = BeamSED->getThroughputElements();
                        
                        fluxCalibrationBeamElements->setFluxVector(&fluxCalibrationFluxVector);
                        thruputBeamElements->setFluxVector(&thruputFluxVector);
                    }
                }
            }
        }
    }
    
    if(StarPlusSky) {
        normalizeBeams = false;
    }
    
    if(!normalizeBeams) {
        NumberofBeams = 0;
    }
    
    unsigned maxNElements = getMaxNumberOfElementsInOrder(Minorder, Maxorder);
    
    float *uniform_wl = new float[numberOfPointsInUniformSample];
    float *uniform_flux = new float[numberOfPointsInUniformSample];
    float *uniform_Beamflux[MAXNUMBEROFBEAMS];
    
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        uniform_Beamflux[beam] = new float[numberOfPointsInUniformSample];
    }
    
    calculateCleanUniformSampleOfContinuum(Minorder,Maxorder,normalizationBinsize,delta_wl,inputWavelengthMaskForUncalContinuum,numberOfPointsInUniformSample,uniform_wl,uniform_flux,uniform_Beamflux, normalizeBeams);
    
    /*
     for(unsigned i=0;i<numberOfPointsInUniformSample;i++) {
     cout << uniform_wl[i] << " " << uniform_flux[i] << endl;
     }
     exit(1);
     */
    float *UncalibratedModelFlux = new float[maxNElements];
    float *UncalibratedModelBeamFlux[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        UncalibratedModelBeamFlux[beam] = new float[maxNElements];
    }
    
    float *elemWavelength = new float[maxNElements];
    operaSpectralEnergyDistribution *BeamSED[MAXNUMBEROFBEAMS];
    operaSpectralElements *uncalibratedBeamFluxElements[MAXNUMBEROFBEAMS];
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasWavelength() &&
            spectralOrder->gethasSpectralElements() &&
            spectralOrder->gethasSpectralEnergyDistribution()) {
            
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            
            unsigned nElements = spectralElements->getnSpectralElements();
            
            spectralEnergyDistribution->setUncalibratedFluxElements(spectralElements);
            
            operaSpectralElements *uncalibratedFluxElements = spectralEnergyDistribution->getUncalibratedFluxElements();
            
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                BeamSED[beam] = spectralOrder->getBeamSED(beam);
                BeamSED[beam]->setUncalibratedFluxElements(spectralOrder->getBeamElements(beam));
                uncalibratedBeamFluxElements[beam] = BeamSED[beam]->getUncalibratedFluxElements();
            }
            
            for(unsigned i=0;i<nElements;i++) {
                elemWavelength[i] =  (float)spectralElements->getwavelength(i);
                
            }
            
            operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_flux,nElements,elemWavelength,UncalibratedModelFlux);
            
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam],nElements,elemWavelength,UncalibratedModelBeamFlux[beam]);
            }
            
            for(unsigned i=0;i<nElements;i++) {
                uncalibratedFluxElements->setFlux((double)UncalibratedModelFlux[i],i);
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    // uncalibratedBeamFluxElements[beam]->setFlux((double)UncalibratedModelBeamFlux[beam][i],i);
                    uncalibratedBeamFluxElements[beam]->setFlux(1.0,i);
                }
            }
            spectralEnergyDistribution->setHasUncalibratedFlux(true);
            for(unsigned beam=0;beam<NumberofBeams;beam++) {
                BeamSED[beam]->setHasUncalibratedFlux(true);
            }
            spectralOrder->applyNormalizationFromExistingContinuum(NULL,NULL,TRUE,normalizeBeams,2);
            
            spectralElements->copyTOnormalizedFlux();
            spectralElements->copyFROMfcalFlux();
        }
    }
    
    double wavelengthForNormalization=0;
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
            if(spectralEnergyDistribution->getwavelengthForNormalization()>0) {
                wavelengthForNormalization = spectralEnergyDistribution->getwavelengthForNormalization();
                break;
            }
        }
    }
    
    double BeamSpectralBinConstant[MAXNUMBEROFBEAMS];
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        BeamSpectralBinConstant[beam] = 1.0;
    }
    double spectralBinConstant =1.0;
    //-- Calculate Sky over Star fiber area ratio to compensate for different apertures.
    /*double SkyOverStarFiberAreaRatio = (2.2*2.2)/(1.6*1.6);
    
    double spectralBinConstant = (double)getFluxAtWavelength(numberOfPointsInUniformSample,uniform_wl,uniform_flux,wavelengthForNormalization);
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        BeamSpectralBinConstant[beam] = (double)getFluxAtWavelength(numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam],wavelengthForNormalization);
        if(StarPlusSky && float(beam) >= float(NumberofBeams)/2) {
            // Sky Fiber
            BeamSpectralBinConstant[beam] /= SkyOverStarFiberAreaRatio;
        }
    }
    */
    // The flux calibration is effectively calculated below, but it expects that SED elements are in there.
    bool AbsoluteCalibration = false;
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasSpectralEnergyDistribution()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            
            if(AbsoluteCalibration) { // use fluxcalibation
                spectralOrder->multiplySpectralElementsBySEDElements(false, spectralBinConstant, NULL);
            } else { // use throughput
                spectralOrder->multiplySpectralElementsBySEDElements(true, spectralBinConstant, NULL);
            }
            spectralElements->copyTOfcalFlux();
            spectralElements->copyFROMrawFlux();
        }
    }
    
    delete[] uniform_wl;
    delete[] uniform_flux;
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        delete[] uniform_Beamflux[beam];
    }
    delete[] UncalibratedModelFlux;
    for(unsigned beam=0;beam<NumberofBeams;beam++) {
        delete[] UncalibratedModelBeamFlux[beam];
    }
    delete[] elemWavelength;
}

unsigned operaSpectralOrderVector::getMaxNumberOfElementsInOrder(int Minorder, int Maxorder) {
    unsigned maxNElements = 0;
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            if(maxNElements < spectralElements->getnSpectralElements()) {
                maxNElements = spectralElements->getnSpectralElements();
            }
        }
    }
    return maxNElements;
}

unsigned operaSpectralOrderVector::getNumberofBeams(int Minorder, int Maxorder) {
    unsigned NumberofBeams = 0;
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            if(NumberofBeams==0 && spectralOrder->getnumberOfBeams()) {
                NumberofBeams = spectralOrder->getnumberOfBeams();
                break;
            }
        }
    }
    return NumberofBeams;
}

unsigned operaSpectralOrderVector::getSpectrumWithinTelluricMask(string inputWavelengthMaskForTelluric, int Minorder, int Maxorder, bool normalized, unsigned normalizationBinsize, double *wavelength, double *spectrum, double *variance) {
    
    double *wl0_vector = new double[MAXNUMBEROFWLRANGES];
    double *wlf_vector = new double[MAXNUMBEROFWLRANGES];
    
    unsigned nRangesInWLMask = readContinuumWavelengthMask(inputWavelengthMaskForTelluric,wl0_vector,wlf_vector);
    
    float *wl = new float[MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
    float *flux = new float[MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
    float *var = new float[MAXORDERS*MAXSPECTRALELEMENTSPERORDER];
    unsigned nelem = 0;
    
    for (int order=Maxorder; order>=Minorder; order--) {
        
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        
        if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
            
            spectralOrder->getSpectralElements()->CreateExtendedvectors(spectralOrder->getSpectralElements()->getnSpectralElements());
            spectralOrder->getSpectralElements()->copyTOrawFlux();

            operaWavelength *wavelength =  spectralOrder->getWavelength();
            
            if(normalized && normalizationBinsize>0) {
                spectralOrder->applyNormalization(normalizationBinsize,0,FALSE,NULL,NULL,TRUE,0);
            }
            
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            spectralElements->setwavelengthsFromCalibration(wavelength);
            
            for (unsigned elemIndex=0; elemIndex<spectralElements->getnSpectralElements(); elemIndex++) {
                for(unsigned rangeElem=0; rangeElem<nRangesInWLMask; rangeElem++) {
                    if(spectralElements->getwavelength(elemIndex) >= wl0_vector[rangeElem] &&
                       spectralElements->getwavelength(elemIndex) <= wlf_vector[rangeElem] ) {
                        
                        if(!isnan((float)spectralElements->getFlux(elemIndex))) {
                            flux[nelem] = (float)spectralElements->getFlux(elemIndex);
                            var[nelem] = (float)spectralElements->getFluxVariance(elemIndex);
                            wl[nelem] = (float)spectralElements->getwavelength(elemIndex);
                            nelem++;
                        }
                        
                        break;
                    }
                }
            }
        }
    }
    
    int *sindex = new int[nelem];
    
    operaArrayIndexSort((int)nelem,wl,sindex);
    
    for(unsigned index=0; index<nelem; index++) {
        wavelength[index] = (double)wl[sindex[index]];
        spectrum[index] = (double)flux[sindex[index]];
        variance[index] = (double)var[sindex[index]];
    }
    
    delete[] wl;
    delete[] flux;
    delete[] var;
    delete[] wl0_vector;
    delete[] wlf_vector;

    return nelem;
}

void operaSpectralOrderVector::calculateRawFluxQuantities(int Minorder, int Maxorder, double *integratedFlux, double *meanFlux, double *maxSNR, double *meanSNR) {
    
    unsigned totalNelem = 0;
    *integratedFlux = 0;
    *maxSNR = -BIG;
    double sumSNR = 0;
    
    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);
        if (spectralOrder->gethasSpectralElements()) {
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
            
            for (unsigned elemIndex=0; elemIndex<spectralElements->getnSpectralElements(); elemIndex++) {
                if(!isnan(spectralElements->getFlux(elemIndex))) {
                    
                    *integratedFlux += spectralElements->getFlux(elemIndex);
                    
                    sumSNR += spectralElements->getFluxSNR(elemIndex);
                    
                    if(spectralElements->getFluxSNR(elemIndex) > *maxSNR) {
                        *maxSNR = spectralElements->getFluxSNR(elemIndex);
                    }
                    
                    totalNelem++;
                }
            }
        }
    }
    if(totalNelem) {
        *meanFlux = *integratedFlux/(double)totalNelem;
        *meanSNR = sumSNR/(double)totalNelem;
    } else {
        *meanFlux = 0;
        *meanSNR = 0;
    }
}

/*
 * void readLibreEspritFlatResponseIntoSED(string filename)
 * \brief Read Libre-Esprit flat response from the .s file.
 */
void operaSpectralOrderVector::readLibreEspritFlatResponseIntoSED(string filename,int Minorder, int Maxorder) {
    
    float *flatresp_tmp = new float[MAXNUMBEROFPOINTSINFLATRESPONSE];
    float *frwavelength_tmp = new float[MAXNUMBEROFPOINTSINFLATRESPONSE];
    //
    // first read in the flat response spectrum from .s file
    //
    unsigned np = 0;
    operaistream fspectrum(filename.c_str());
    if (fspectrum.is_open()) {
        string dataline;
        unsigned line = 0;
        float w, fresp;
        
        while (fspectrum.good()) {
            getline(fspectrum, dataline);
            if (strlen(dataline.c_str())) {
                if (dataline.c_str()[0] == '*') {
                    // skip comments
                } else if (line == 0) {
                    sscanf(dataline.c_str(), "%u", &count);
                    line++;
                } else {
                    sscanf(dataline.c_str(), "%f %f", &w, &fresp);
                    flatresp_tmp[np] = fresp;
                    frwavelength_tmp[np] = w;
                    np++;
                    
                    line++;
                }
            }
        }
        fspectrum.close();
    }
    
    float *flatresp = new float[MAXNUMBEROFPOINTSINFLATRESPONSE];
    float *frwavelength = new float[MAXNUMBEROFPOINTSINFLATRESPONSE];

    calculateUniformSample(np,frwavelength_tmp,flatresp_tmp,np,frwavelength,flatresp);
    
    double wavelengthForNormalization = 0;
    float maxflatresp = -BIG;
    // get wavelength at maximum flat response, which should be around 1 since this quantity is normalized.
    for (unsigned i=0; i<np; i++) {
        if (flatresp[i] > maxflatresp) {
            maxflatresp = flatresp[i];
            wavelengthForNormalization = frwavelength[i];
        }
    }
    
    // create memory space for flat response interpolations within orders
    unsigned maxNElements = getMaxNumberOfElementsInOrder(Minorder, Maxorder);
    float *flatResponseModel = new float[maxNElements];
    float *elemWavelength = new float[maxNElements];

    operaSpectralEnergyDistribution *BeamSED = NULL;
    operaSpectralElements *beamFluxcalibration = NULL;
    operaSpectralElements *beamThroughput = NULL;

    for (int order=Minorder; order<=Maxorder; order++) {
        operaSpectralOrder *spectralOrder = GetSpectralOrder(order);

        if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {

            operaWavelength *wavelength =  spectralOrder->getWavelength();
            operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();

            spectralElements->setwavelengthsFromCalibration(wavelength);

            unsigned nElements = spectralElements->getnSpectralElements();
            
            spectralOrder->createSpectralEnergyDistributionElements(nElements);
           
            operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();

            spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
            
            operaSpectralElements *FluxCalibration = spectralEnergyDistribution->getFluxCalibrationElements();
            operaSpectralElements *InstrumentThroughput = spectralEnergyDistribution->getThroughputElements();
            
            for (unsigned beam=0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                BeamSED = spectralOrder->getBeamSED(beam);
                beamFluxcalibration = BeamSED->getFluxCalibrationElements();
                beamThroughput = BeamSED->getThroughputElements();
            }
            
            for (unsigned elemIndex=0; elemIndex<nElements; elemIndex++) {
                elemWavelength[elemIndex] =  spectralElements->getwavelength(elemIndex);
            }

            operaFitSpline(np,frwavelength,flatresp,nElements,elemWavelength,flatResponseModel);
            
            for (unsigned elemIndex=0; elemIndex<nElements; elemIndex++) {
                
                double wl = spectralElements->getwavelength(elemIndex);
                double fluxcal = 1.0/(double)flatResponseModel[elemIndex];
                double throughput = 1.0/(double)flatResponseModel[elemIndex];
                double beamfluxcal = 1.0/(double)flatResponseModel[elemIndex];
                double beamthroughput = 1.0/(double)flatResponseModel[elemIndex];
                
                FluxCalibration->setwavelength(wl, elemIndex);
                FluxCalibration->setFlux(fluxcal, elemIndex);
                FluxCalibration->setFluxVariance(fluxcal, elemIndex);
                FluxCalibration->setphotoCenter(0.0, 0.0, elemIndex);
                
                InstrumentThroughput->setwavelength(wl, elemIndex);
                InstrumentThroughput->setFlux(throughput, elemIndex);
                InstrumentThroughput->setFluxVariance(throughput, elemIndex);
                InstrumentThroughput->setphotoCenter(0.0, 0.0, elemIndex);
                
                for (unsigned beam=0; beam < spectralOrder->getnumberOfBeams(); beam++) {
        
                    beamFluxcalibration->setwavelength(wl, elemIndex);
                    beamFluxcalibration->setFlux(beamfluxcal, elemIndex);
                    beamFluxcalibration->setFluxVariance(beamfluxcal, elemIndex);
                    beamFluxcalibration->setphotoCenter(0.0, 0.0, elemIndex);
                    
                    beamThroughput->setwavelength(wl, elemIndex);
                    beamThroughput->setFlux(beamthroughput, elemIndex);
                    beamThroughput->setFluxVariance(beamthroughput, elemIndex);
                    beamThroughput->setphotoCenter(0.0, 0.0, elemIndex);
                }
            }
            
            spectralOrder->sethasSpectralEnergyDistribution(true);
            FluxCalibration->setHasWavelength(true);
            InstrumentThroughput->setHasWavelength(true);
            spectralEnergyDistribution->setHasFluxCalibration(true);
            spectralEnergyDistribution->setHasInstrumentThroughput(true);
        }
    }
    delete[] frwavelength;
    delete[] flatresp;
    delete[] frwavelength_tmp;
    delete[] flatresp_tmp;
}

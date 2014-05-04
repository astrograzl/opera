#ifndef OPERASPECTRALELEMENTS_H
#define OPERASPECTRALELEMENTS_H

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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include "libraries/operaLibCommon.h"			// for CMatrix
#include "libraries/Polynomial.h"	
#include "libraries/operaFluxVector.h"			// for operaFluxVector
#include "libraries/operaSpectralOrder.h"		// for operaSpectralOrder_t
#include "libraries/operaPolarimetry.h"			// for operaInstrumentProfile
#include "libraries/operaInstrumentProfile.h"	// for operaInstrumentProfile
#include "libraries/operaFluxVector.h"          // for operaFluxVectors
#include "libraries/operaWavelength.h"          // for operaWavelength

#define MAXSPECTRALELEMENTSPERORDER 8000

/*! 
 * \sa class operaSpectralElements
 * \brief operaSpectralElements
 * \details A spectral element (SE) consists of a spectrally resolved 
 * \details subdivision of a spectral order taken in the dispersion direction i.e. along the trace polynomial
 * \details more or less in the y direction for espadons, but is actually curved. 
 * \details The minimal width of the spectral element is defined by the resolution 
 * \details of the spectrograph (RE). The actual width could be defined as any size 
 * \details smaller than the order size and as long as RE < RS, where RE is the 
 * \details resolution of the element.
 * \return none
 * \file operaSpectralElements.h
 * \ingroup libraries
 */

class operaSpectralOrder;

class operaSpectralElements {
	
private:
	unsigned nSpectralElements;			// this is a count, the vector is not retained
	unsigned maxnSpectralElements;		// this is a count, the vector is not retained
	double elementHeight;				// height of spectral element in pixel units measured along the dispersion direction (CenterPolynomial)
	operaSpectralOrder_t SpectrumType;	// what kind of spectrum is this

	operaFluxVector *fluxvector;		// (dbl in counts)
	double *XCorrelation;               // cross-correlation vector. Used for detection of spectral lines.
	double *photoCenterX, *photoCenterY;// x,y coordinates of the photocenter in the image reference frame (pixel units) 
	double *distd;						// distance in pixel units measured along the dispersion, starting at the first element d[0] = 0 
	double *wavelength;					// wavelength in nm 
	double *fluxSNR;					//SNR at each spectralelement 

	// Extended spectrum includes various variations of wl and flu data...
	double *tell;						// telluric wl information
	double *rvel;						// barycentric wl information
	
	operaFluxVector *rawFlux;			// raw flux
	operaFluxVector *normalizedFlux;	// normalized flux
	operaFluxVector *fcalFlux;			// flux calibrated flux
	
	bool hasRawFlux;
	bool hasStandardFlux;
	bool hasOptimalFlux;
	bool hasOperaOptimalFlux;
	bool hasExtendedBeamFlux;
    
	bool hasXCorrelation;    
	bool hasWavelength;    
	bool hasDistance;    
	bool hasFluxSNR;    
	
public:
	
	/*
	 * Constructor
	 */
	operaSpectralElements();
	operaSpectralElements(unsigned nElements);
	operaSpectralElements(unsigned nElements, operaSpectralOrder_t format, bool extended = false);
	/*
	 * Destructor
	 */
	~operaSpectralElements();
	
	/*
	 * getters/setters
	 */
	
	void Createvectors(unsigned nElements, bool extended = false);
	
	void CreateExtendedvectors(unsigned nElements); 
	
	void Createvectors(unsigned nElements, operaSpectralOrder_t format, bool extended = false);
	
	void Deletevectors(void);
	
	void Resizevector(unsigned nElements, operaSpectralOrder_t format);
	
	void setnSpectralElements(unsigned nElem);
	
	unsigned getnSpectralElements(void);
	
	void setelementHeight(double Height);
	
	double getelementHeight(void);
	
	bool getHasRawSpectrum() { return hasRawFlux; };
	
	void setHasRawSpectrum(bool HasRawFlux) { hasRawFlux = HasRawFlux; };
	
	bool getHasStandardSpectrum() { return hasStandardFlux; };
	
	void setHasStandardSpectrum(bool HasStandardFlux) { hasStandardFlux = HasStandardFlux; };
	
	bool getHasOptimalSpectrum() { return hasOptimalFlux; };
	
	void setHasOptimalSpectrum(bool HasOptimalFlux) { hasOptimalFlux = HasOptimalFlux; };
	
	bool getHasOperaOptimalSpectrum() { return hasOperaOptimalFlux; };
	
	void setHasOperaOptimalSpectrum(bool HasOperaOptimalFlux) { hasOperaOptimalFlux = HasOperaOptimalFlux; };
	
	bool getHasXCorrelation() { return hasXCorrelation; };
	
	void setHasXCorrelation(bool HasXCorrelation) { hasXCorrelation = HasXCorrelation; };
    
	bool getHasWavelength() { return hasWavelength; };
	
	void setHasWavelength(bool HasWavelength) { hasWavelength = HasWavelength; };
    
	bool getHasDistance() { return hasDistance; };
	
	void setHasDistance(bool HasDistance) { hasDistance = HasDistance; };
    
	bool getHasFluxSNR() { return hasFluxSNR; };
	
	void setHasFluxSNR(bool HasFluxSNR) { hasFluxSNR = HasFluxSNR; };
    
	bool getHasExtendedBeamFlux() { return hasExtendedBeamFlux; };
	
	void setHasExtendedBeamFlux(bool HasExtendedBeamFlux) { hasExtendedBeamFlux = HasExtendedBeamFlux; };
    
	/*
	 * Flux
	 */
	
	operaFluxVector *getFluxVector(void);
	
	double getFlux(unsigned indexElem);
	
	void setFluxVector(operaFluxVector *Flux);
	
	void setFlux(double Flux, unsigned indexElem);
	
	double getFluxVariance(unsigned indexElem);
		
	void setFluxVariance(double FluxVariance, unsigned indexElem);
	
	double getFluxSNR(unsigned indexElem);
	
	void setFluxSNR(double FluxSNR, unsigned indexElem);
	
	/*
	 * others
	 */
	double getphotoCenterX(unsigned indexElem);
	double getphotoCenterY(unsigned indexElem);
	void setphotoCenter(double x, double y, unsigned indexElem);
	double getdistd(unsigned indexElem);
	void setdistd(double Distd, unsigned indexElem);
	double getwavelength(unsigned indexElem);
	void setwavelength(double Wavelength, unsigned indexElem);
    void setwavelengthsFromCalibration(operaWavelength *Wavelength);
	double getXCorrelation(unsigned indexElem);  
	void setXCorrelation(double Xcorr, unsigned indexElem);    

	double gettell(unsigned indexElem);  
	void settell(double value, unsigned indexElem);    
	void copyTOtell(void);
	void copyFROMtell(void);
	double getrvel(unsigned indexElem);  
	void setrvel(double value, unsigned indexElem);    
	void copyTOrvel(void);
	void copyFROMrvel(void);
	double getnormalizedFlux(unsigned indexElem);  
	void setnormalizedFlux(double value, unsigned indexElem);    
	void copyTOnormalizedFlux(void);
	void copyFROMnormalizedFlux(void);
	double getfcalFlux(unsigned indexElem);  
	void setfcalFlux(double value, unsigned indexElem);    
	void copyTOfcalFlux(void);
	void copyFROMfcalFlux(void);
	double getrawFlux(unsigned indexElem);  
	void setrawFlux(double value, unsigned indexElem);    
	void copyTOrawFlux(void);
	void copyFROMrawFlux(void);
	
	operaSpectralOrder_t getSpectrumType(void);	
	
	void setSpectrumType(operaSpectralOrder_t format);
	
	/*
	 * clone spectralelements
	 */
	void setSpectralElements(operaSpectralElements &SpectralElements);

	/*
	 * Other Methods
	 */
    double getFluxSum(void);
    
};
#endif

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                     ****
 ********************************************************************
 Library name: operaWavelength
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

#ifndef OPERAWAVELENGTH_H
#define OPERAWAVELENGTH_H

#include "libraries/operaLibCommon.h"			// for doubleValue_t
#include "libraries/Polynomial.h"	
#include "libraries/operaWavelength.h"	

#define MAXORDEROFWAVELENGTHPOLYNOMIAL 10

/*! 
 * \sa class operaWavelength
 * \brief Encapsulation of Wavelength information.
 * \return none
 * \file operaWavelength.h
 * \ingroup libraries
 */

class operaSpectralOrder;

class operaWavelength {
	
private:
	double dmin;
	double dmax;		
    
    unsigned maxNDataPoints;    
    unsigned nDataPoints;
    
    double *distanceData;
    double *wavelengthData;
    double *wavelengthErrors; 
    unsigned *matchAtlasindex;
    unsigned *matchComparisonindex;

    unsigned maxNAtlasLines;
    unsigned nAtlasLines;

    double *atlasLinesflux;
    double *atlasLineswl;
    double *atlasLineswlError;    

    unsigned maxNComparisonLines;   
    unsigned nComparisonLines;   

    double *comparisonLinesflux;
    double *comparisonLinespix;
    double *comparisonLinespixError;    
    double *comparisonLineswl;    

    doubleValue_t spectralResolution;
    double radialVelocityPrecision;    
    
    double xcorrelation;
    
    Polynomial *wavelengthPolynomial; // lambda(d)
	
public:
		
	/*
	 * Constructors
	 */
	operaWavelength();
	operaWavelength(unsigned Coeffs);
	~operaWavelength();
	
	/*
	 * Methods
	 */
	double getDmin(void) const;
	
	double getDmax(void) const;
	
	void setDmin(double Dmin);
	
	void setDmax(double Dmax);
	
	double getDistance(unsigned index) const;
	
	double getWavelength(unsigned index) const;
	
	double getWavelengthError(unsigned index) const;

	unsigned getMatchAtlasIndex(unsigned index) const;
    
	unsigned getMatchComparisonIndex(unsigned index) const;
    
	double getcentralWavelength(void) const;
	
	double getinitialWavelength(void) const;
	
	double getfinalWavelength(void) const;
	
	Polynomial *getWavelengthPolynomial(void);
	
	const Polynomial *getWavelengthPolynomial(void) const;
	
	void CalculateWavelengthSolution(unsigned maxcoeffs, bool witherrors);
    
    void RefineWavelengthSolution(unsigned ncoeffs, bool witherrors);    
	    
    /*
	 * Clean wavelength and distance data for fitting
	 */
    
    void setnDataPoints(unsigned NDataPoints);
    
    unsigned getnDataPoints(void) const;
	
	void createDataVectors(unsigned NDataPoints);
    
    void createDataVectors(unsigned NDataPoints, double *WavelengthData, double *WavelengthErrors, double *DistanceData);

	void deleteDataVectors(void);    
	/*
	 * Atlas wavelength data
	 */
	
    void setnAtlasLines(unsigned NAtlasLines);
    
    unsigned getnAtlasLines(void) const;   

    double getatlasLinesflux(unsigned index) const;
    
    double getatlasLineswl(unsigned index) const;
    
    double getatlasLineswlError(unsigned index) const;    
    
    void setatlasLinesflux(double AtlasLinesflux, unsigned index);
    
    void setatlasLineswl(double AtlasLineswl, unsigned index);     
    
    void setatlasLineswl(double AtlasLineswl, double AtlasLineswlError, unsigned index);     

    void createAtlasDataVectors(unsigned NAtlasLines);
    
    void createAtlasDataVectors(unsigned NAtlasLines, double *AtlasLineswl, double *AtlasLineswlError, double *AtlasLinesflux);
    
    void deleteAtlasDataVectors(void);
        
	/*
	 *  Comparison pixel data
	 */    
    
    void setnComparisonLines(unsigned NComparisonLines);

    unsigned getnComparisonLines(void) const;

    double getcomparisonLinesflux(unsigned index) const;
    
    double getcomparisonLinespix(unsigned index) const;
    
    double getcomparisonLinespixError(unsigned index) const;
    
    double getcomparisonLineswl(unsigned index) const;
    
    void setcomparisonLinesflux(double ComparisonLinesflux, unsigned index);
    
    void setcomparisonLinespix(double ComparisonLinespix, double ComparisonLinespixError, unsigned index); 
    
    void setcomparisonLineswl(double ComparisonLineswl, unsigned index);  
    
    void recalculateComparisonLineswlVector(void);   
    
    void createComparisonDataVectors(unsigned NComparisonLines);
    
    void createComparisonDataVectors(unsigned NComparisonLines, double *ComparisonLinespix, double *ComparisonLinespixError, double *ComparisonLinesflux);    
    
    void deleteComparisonDataVectors(void); 
    
    
	 /* 
      * Other Methods     
      */
    void setSpectralResolution(doubleValue_t Resolution);
    
    doubleValue_t getSpectralResolution(void) const;
	
    void calculateSpectralResolution(doubleValue_t ResolutionElementInPixels);
    
    double evaluateWavelength(double distanceValue) const;
    
    double convertPixelToWavelength(double DeltaDistanceInPixels) const;
    
    void setRadialVelocityPrecision(double radialvelocityprecision);    
    
    double getRadialVelocityPrecision(void) const;
    
    void calculateRadialVelocityPrecision(void);
    
    double calculateWavelengthRMSPrecision(void);
    
    double calculateWavelengthMedianPrecision(void);   
    
    void refineWavelengthSolutionOfSecondOrderByXCorrelation(unsigned nPointsPerParameter, double parameterRangetoSearch);
    
    void refineWavelengthSolutionByFindingMaxMatching(unsigned NpointsPerPar, double ParRangeSizeInPerCent, double acceptableMismatch);
    
    unsigned createAtlasSimulatedSpectrumWithConstantFlux(double *outputwl, double *outputSpectrum, unsigned nstepspersigma);
    
    unsigned createComparisonSimulatedSpectrumWithConstantFlux(double *outputwl, double *outputSpectrum, unsigned nstepspersigma);
    
    unsigned createAtlasSimulatedSpectrum(double *outputwl, double *outputSpectrum, unsigned nstepspersigma);
        
    unsigned createComparisonSimulatedSpectrum(double *outputwl, double *outputSpectrum, unsigned nstepspersigma);

	double getxcorrelation(void) const;
	
	void setxcorrelation(double Xcorrelation); 
    
    void calculateXCorrelation(void);
    
    void matchAtlaswithComparisonLines(double acceptableMismatch);    
            
    double getPerCentageOfComparisonMatch(void) const;
    
    double getPerCentageOfAtlasMatch(void) const;
    
    void filterDataPointsBySigmaClip(double nsig);
    
    void filterDataPointsByErrorClip(double nsig);
    
    void applyRadialVelocityCorrection(double rvshift_InKPS);

};
#endif

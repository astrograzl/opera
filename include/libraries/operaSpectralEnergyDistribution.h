#ifndef OPERASPECTRALENERGYDISTRIBUTION_H
#define OPERASPECTRALENERGYDISTRIBUTION_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralEnergyDistribution
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

#include "libraries/Polynomial.h"	

#ifndef MAXNUMBEROFMASKEDREGIONS
#define MAXNUMBEROFMASKEDREGIONS 300
#endif

/*! 
 * \sa class operaSpectralEnergyDistribution
 * \brief Encapsulation of Spectral Energy Distribution information.
 * \ingroup libraries
 */

class operaSpectralOrder;

class operaSpectralEnergyDistribution {
	
private:    
    
    unsigned nDataPoints;
    
    unsigned maxnDataPoints;
    
    double *distanceData;
    
    double *wavelengthData;
   	
    double wavelengthForNormalization;
    
    operaFluxVector *fluxData;                      // flux measurements in instrumental units (ADU/pixel^2/spectralElement)
   	
    operaFluxVector *spectralEnergyData;            // energy flux density obtained from standard measurements
                                                    // flux ( ergs/(cm*cm*s*A) * 10**16 )
      
    operaSpectralElements *uncalibratedFluxElements;         // operaSpectralElements for uncalibrated flux
    
    operaSpectralElements *calibratedFluxElements;         // operaSpectralElements for calibrated flux

    
    operaSpectralElements *fluxCalibration;         // operaSpectralElements for flux calibration
    
    operaSpectralElements *instrumentThroughput;    // operaSpectralElements for instrument throughput

    unsigned numberOfMaskedRegions;                     // number of wavelength regions to mask
    double MaskWavelength1[MAXNUMBEROFMASKEDREGIONS];   // vector of initial wavelengths for masked regions
    double MaskWavelength2[MAXNUMBEROFMASKEDREGIONS];   // vector of final wavelengths for masked regions
    
	bool hasFluxData;
	bool hasSpectralEnergyData;   
    
	bool hasUncalibratedFlux;
	bool hasCalibratedFlux;
	bool hasFluxCalibration;
   	bool hasInstrumentThroughput;      

public:
		
	/*
	 * Constructors
	 */
    
	operaSpectralEnergyDistribution();
    
	operaSpectralEnergyDistribution(unsigned NDataPoints);   

    operaSpectralEnergyDistribution(unsigned NDataPoints, unsigned nElements);  
    
    operaSpectralEnergyDistribution(operaSpectralElements &inputFluxCalibrationElements);   

	/*
	 * Destructor
	 */   
    
	~operaSpectralEnergyDistribution();
	
	/*
	 * Methods for managing data
	 */    
    
    void setnDataPoints(unsigned NDataPoints);
    
    unsigned getnDataPoints(void) const;
    
    void setwavelengthForNormalization(double WavelengthForNormalization);
    
    double getwavelengthForNormalization(void) const;
        
	void createDataVectors(unsigned NDataPoints);

	void deleteDataVectors(void);
    
    void setdistanceData(double Distance, unsigned index);
    
    void setwavelengthData(double Wavelength, unsigned index);
    
    void setfluxData(double Flux, unsigned index); 
    
    void setspectralEnergyData(double SpectralEnergy, unsigned index);    
   
    double getdistanceData(unsigned index) const;
    
    double getwavelengthData(unsigned index) const;
    
    double getfluxData(unsigned index) const;
    
    double getspectralEnergyData(unsigned index) const;
    
	bool getHasFluxData(void) const { return hasFluxData; };
    
	void setHasFluxData(bool HasFluxData) { hasFluxData = HasFluxData; };
    
	bool getHasSpectralEnergyData(void) const { return hasSpectralEnergyData; };
    
	void setHasSpectralEnergyData(bool HasSpectralEnergyData) { hasSpectralEnergyData = HasSpectralEnergyData; };
    
	/*
	 * Methods for managing flux calibration and throughput elements
	 */
	  
	void createFluxCalibrationElements(unsigned nElements);     
    
	void deleteFluxCalibrationElements(void);    
	
	void createThroughputElements(unsigned nElements);     
    
	void deleteThroughputElements(void); 
    
	void createUncalibratedFluxElements(unsigned nElements);     
    
	void deleteUncalibratedFluxElements(void); 
    
	void createCalibratedFluxElements(unsigned nElements);     
    
	void deleteCalibratedFluxElements(void);     
    
    void setUncalibratedFluxElements(operaSpectralElements *UncalibratedFluxElements);

    void setCalibratedFluxElements(operaSpectralElements *CalibratedFluxElements);    
    
    operaSpectralElements* getUncalibratedFluxElements(void){ return uncalibratedFluxElements; }
    
    const operaSpectralElements* getUncalibratedFluxElements(void) const { return uncalibratedFluxElements; }
    
    operaSpectralElements* getCalibratedFluxElements(void){ return calibratedFluxElements; }
    
    const operaSpectralElements* getCalibratedFluxElements(void) const { return calibratedFluxElements; }
    
    void setFluxCalibrationElements(operaSpectralElements *FluxCalibrationElements);
    
    void setThroughputElements(operaSpectralElements *ThroughputElements);
    
    operaSpectralElements* getFluxCalibrationElements(void){ return fluxCalibration; }
    
    const operaSpectralElements* getFluxCalibrationElements(void) const { return fluxCalibration; }
    
    operaSpectralElements* getThroughputElements(void){ return instrumentThroughput; }
    
    const operaSpectralElements* getThroughputElements(void) const { return instrumentThroughput; }
        
	bool getHasUncalibratedFlux(void) const { return hasUncalibratedFlux; };

	void setHasUncalibratedFlux(bool HasUncalibratedFlux) { hasUncalibratedFlux = HasUncalibratedFlux; };

	bool getHasCalibratedFlux(void) const { return hasCalibratedFlux; };

	void setHasCalibratedFlux(bool HasCalibratedFlux) { hasCalibratedFlux = HasCalibratedFlux; };

	bool getHasFluxCalibration(void) const { return hasFluxCalibration; };

	void setHasFluxCalibration(bool HasFluxCalibration) { hasFluxCalibration = HasFluxCalibration; };
    
	bool getHasInstrumentThroughput(void) const { return hasInstrumentThroughput; };

	void setHasInstrumentThroughput(bool HasInstrumentThroughput) { hasInstrumentThroughput = HasInstrumentThroughput; };
    
    /*
     * wavelength mask regions
     */
    
    void setnumberOfMaskedRegions(unsigned nMaskedRegions);
    
    unsigned getnumberOfMaskedRegions(void) const;
    
    void setMaskRegion(unsigned index, double Wavelength1, double Wavelength2);
    
    double getMaskRegionWavelength1(unsigned index) const;
    
    double getMaskRegionWavelength2(unsigned index) const;
 
    void readMaskedRegionsFromFile(string MaskedRegionsFileName);

    bool isItMasked(double Wavelength);
    
    /*
     * Other Methods     
     */
    void measureUncalibratedContinuum(operaSpectralElements &uncalibratedFluxElements, unsigned binsize, unsigned nsigcut);
    
    void calculateUncalibratedElements(unsigned binsize, unsigned nsigcut);

    void populateUncalibratedElementsFromContinuumData(void);

    void calculateCalibratedElements(unsigned nPointsInReference, double *refwl, double *refflux);
    
    void calculateSEDelements(double spectralBinConstant,double referenceFluxForNormalization,double uncalibratedContinuumFluxForNormalization);
    
};
#endif

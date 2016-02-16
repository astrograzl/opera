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
#include <stdio.h>
#include <fstream>

#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralElements.h"            // for operaSpectralElements
#include "libraries/operaSpectralEnergyDistribution.h"  // for operaSpectralEnergyDistribution
#include "libraries/ladfit.h"                           // for ladfit
#include "libraries/operaFit.h"                         // for operaFitSpline

/*!
 * operaSpectralEnergyDistribution
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the spectral energy distribution object.
 * \file operaSpectralEnergyDistribution.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaSpectralEnergyDistribution
 * \brief Encapsulation of Wavelength information.
 * \return none
 */

/*
 * Constructors
 */

operaSpectralEnergyDistribution::operaSpectralEnergyDistribution() :
nDataPoints(0),
maxnDataPoints(0),
distanceData(NULL),
wavelengthData(NULL),
wavelengthForNormalization(0),
fluxData(NULL),
spectralEnergyData(NULL),
uncalibratedFluxElements(NULL),
calibratedFluxElements(NULL),
fluxCalibration(NULL),
instrumentThroughput(NULL),
numberOfMaskedRegions(0),
hasFluxData(false),   
hasSpectralEnergyData(false),
hasUncalibratedFlux(false),
hasCalibratedFlux(false),
hasFluxCalibration(false),
hasInstrumentThroughput(false) 
{
	throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
}

operaSpectralEnergyDistribution::operaSpectralEnergyDistribution(unsigned NDataPoints) :
nDataPoints(0),
maxnDataPoints(0),
distanceData(NULL),
wavelengthData(NULL),
wavelengthForNormalization(0),
fluxData(NULL),
spectralEnergyData(NULL),
uncalibratedFluxElements(NULL),
calibratedFluxElements(NULL),
fluxCalibration(NULL),
instrumentThroughput(NULL),
numberOfMaskedRegions(0),
hasFluxData(false),   
hasSpectralEnergyData(false),
hasUncalibratedFlux(false),
hasCalibratedFlux(false),
hasFluxCalibration(false),
hasInstrumentThroughput(false) 
{
    if (NDataPoints == 0) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	createDataVectors(NDataPoints);
    for(unsigned index=0;index<MAXNUMBEROFMASKEDREGIONS;index++) {
        MaskWavelength1[index] = 0;
        MaskWavelength2[index] = 0;
    }
}

operaSpectralEnergyDistribution::operaSpectralEnergyDistribution(operaSpectralElements &inputFluxCalibrationElements) :
nDataPoints(0),
maxnDataPoints(0),
distanceData(NULL),
wavelengthData(NULL),
wavelengthForNormalization(0),
fluxData(NULL),
spectralEnergyData(NULL),
uncalibratedFluxElements(NULL),
calibratedFluxElements(NULL),
fluxCalibration(NULL),
instrumentThroughput(NULL),
numberOfMaskedRegions(0),
hasFluxData(false),   
hasSpectralEnergyData(false),
hasUncalibratedFlux(false),
hasCalibratedFlux(false),
hasFluxCalibration(false),
hasInstrumentThroughput(false) 
{
    unsigned nElements = inputFluxCalibrationElements.getnSpectralElements();
    if (nElements == 0) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	createFluxCalibrationElements(nElements);
    for(unsigned index=0;index<MAXNUMBEROFMASKEDREGIONS;index++) {
        MaskWavelength1[index] = 0;
        MaskWavelength2[index] = 0;
    }

}

operaSpectralEnergyDistribution::operaSpectralEnergyDistribution(unsigned NDataPoints, unsigned nElements) :
nDataPoints(0),
maxnDataPoints(0),
distanceData(NULL),
wavelengthData(NULL),
wavelengthForNormalization(0),
fluxData(NULL),
spectralEnergyData(NULL),
uncalibratedFluxElements(NULL),
calibratedFluxElements(NULL),
fluxCalibration(NULL),
instrumentThroughput(NULL),
numberOfMaskedRegions(0),
hasFluxData(false),   
hasSpectralEnergyData(false),
hasUncalibratedFlux(false),
hasCalibratedFlux(false),
hasFluxCalibration(false),
hasInstrumentThroughput(false) 
{
    if (NDataPoints == 0 || nElements == 0) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    createDataVectors(NDataPoints); 
    createUncalibratedFluxElements(nElements);       
	createCalibratedFluxElements(nElements);         
	createFluxCalibrationElements(nElements);    
	createThroughputElements(nElements);
    for(unsigned index=0;index<MAXNUMBEROFMASKEDREGIONS;index++) {
        MaskWavelength1[index] = 0;
        MaskWavelength2[index] = 0;
    }
}

/*
 * Destructor
 */  

operaSpectralEnergyDistribution::~operaSpectralEnergyDistribution() {
    deleteDataVectors();
    
    deleteUncalibratedFluxElements();         
	deleteCalibratedFluxElements();     
	deleteFluxCalibrationElements();
    deleteThroughputElements();    
}

/*
 * Methods for managing data
 */  

// DT modified to only reduce the number of datapoints to avoid memory overwrites
// Jan 2013
void operaSpectralEnergyDistribution::setnDataPoints(unsigned NDataPoints) {
    if (NDataPoints > maxnDataPoints) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    nDataPoints = NDataPoints;
}

unsigned operaSpectralEnergyDistribution::getnDataPoints(void) const {
    return nDataPoints;
}

void operaSpectralEnergyDistribution::setwavelengthForNormalization(double WavelengthForNormalization) {
    wavelengthForNormalization = WavelengthForNormalization;
}

double operaSpectralEnergyDistribution::getwavelengthForNormalization(void) const {
    return wavelengthForNormalization;
}

void operaSpectralEnergyDistribution::createDataVectors(unsigned NDataPoints) {
	// DT Jan 2013 only reallocated storage if you need to, otherwise reuse
	if (NDataPoints > maxnDataPoints) {
		deleteDataVectors();
		distanceData = new double[NDataPoints];  
		wavelengthData = new double[NDataPoints];
		fluxData = new operaFluxVector(NDataPoints); 
		spectralEnergyData = new operaFluxVector(NDataPoints);         
	}
	nDataPoints = maxnDataPoints = NDataPoints;
}

void operaSpectralEnergyDistribution::deleteDataVectors(void) {
    if(distanceData) {
        free(distanceData);
        distanceData = NULL;
    }
    if(wavelengthData) {
        free(wavelengthData);  
        wavelengthData = NULL;
    }
    if(fluxData){
        delete fluxData;    
        fluxData = NULL;
    }
    if(spectralEnergyData) {
        delete spectralEnergyData;      
        spectralEnergyData = NULL;
    }
    nDataPoints = maxnDataPoints = 0;
}

void operaSpectralEnergyDistribution::setdistanceData(double Distance, unsigned index) {
    if (index > nDataPoints) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}    
    distanceData[index] = Distance;
}

void operaSpectralEnergyDistribution::setwavelengthData(double Wavelength, unsigned index) {
    if (index > nDataPoints) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}    
    wavelengthData[index] = Wavelength;   
}

void operaSpectralEnergyDistribution::setfluxData(double Flux, unsigned index) {
    if (index > fluxData->getlength()) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}    
    fluxData->setflux(Flux,index);    
}

void operaSpectralEnergyDistribution::setspectralEnergyData(double SpectralEnergy, unsigned index) {
    if (index > spectralEnergyData->getlength()) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}    
    spectralEnergyData->setflux(SpectralEnergy,index);
}

double operaSpectralEnergyDistribution::getdistanceData(unsigned index) const {
    if (index > nDataPoints) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}       
    return distanceData[index];
}

double operaSpectralEnergyDistribution::getwavelengthData(unsigned index) const {
    if (index > nDataPoints) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}    
    return wavelengthData[index];    
}

double operaSpectralEnergyDistribution::getfluxData(unsigned index) const {
    if (index > fluxData->getlength()) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}    
    return fluxData->getflux(index);    
}

double operaSpectralEnergyDistribution::getspectralEnergyData(unsigned index) const {
    if (index > spectralEnergyData->getlength()) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}    
    return spectralEnergyData->getflux(index);    
}

/*
 * Methods for managing flux calibration and throughput elements
 */

void operaSpectralEnergyDistribution::createFluxCalibrationElements(unsigned nElements) {
    deleteFluxCalibrationElements();
    fluxCalibration = new operaSpectralElements(nElements);
}

void operaSpectralEnergyDistribution::createThroughputElements(unsigned nElements) {
    deleteThroughputElements();
    instrumentThroughput = new operaSpectralElements(nElements);    
}

void operaSpectralEnergyDistribution::deleteFluxCalibrationElements(void) {
    if(fluxCalibration) {
        delete fluxCalibration;
        fluxCalibration = NULL;                
    }      
}

void operaSpectralEnergyDistribution::deleteThroughputElements(void) {
    if(instrumentThroughput) {
        delete instrumentThroughput;
        instrumentThroughput = NULL;        
    }
}

void operaSpectralEnergyDistribution::createUncalibratedFluxElements(unsigned nElements) {
    deleteUncalibratedFluxElements();
    uncalibratedFluxElements = new operaSpectralElements(nElements);     
}

void operaSpectralEnergyDistribution::deleteUncalibratedFluxElements(void) {
    if(uncalibratedFluxElements) {
        delete uncalibratedFluxElements;
        uncalibratedFluxElements = NULL;
    }    
    
}

void operaSpectralEnergyDistribution::createCalibratedFluxElements(unsigned nElements) {
    deleteCalibratedFluxElements();
    calibratedFluxElements = new operaSpectralElements(nElements);     
}

void operaSpectralEnergyDistribution::deleteCalibratedFluxElements(void) {
    if(calibratedFluxElements) {
        delete calibratedFluxElements;
        calibratedFluxElements = NULL;                
    }    
}

void operaSpectralEnergyDistribution::setUncalibratedFluxElements(operaSpectralElements *UncalibratedFluxElements) {
    createUncalibratedFluxElements(UncalibratedFluxElements->getnSpectralElements());
	*uncalibratedFluxElements = *UncalibratedFluxElements;
}

void operaSpectralEnergyDistribution::setCalibratedFluxElements(operaSpectralElements *CalibratedFluxElements) {
    createCalibratedFluxElements(CalibratedFluxElements->getnSpectralElements());
	*calibratedFluxElements = *CalibratedFluxElements;
}

void operaSpectralEnergyDistribution::setFluxCalibrationElements(operaSpectralElements *FluxCalibrationElements) {
    createFluxCalibrationElements(FluxCalibrationElements->getnSpectralElements());
	*fluxCalibration = *FluxCalibrationElements;
}

void operaSpectralEnergyDistribution::setThroughputElements(operaSpectralElements *ThroughputElements) {
    createThroughputElements(ThroughputElements->getnSpectralElements());
	*instrumentThroughput = *ThroughputElements;
}

/*
 * wavelength mask regions
 */

void operaSpectralEnergyDistribution::setnumberOfMaskedRegions(unsigned nMaskedRegions) {
    if (numberOfMaskedRegions > MAXNUMBEROFMASKEDREGIONS) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    numberOfMaskedRegions = nMaskedRegions;
}

unsigned operaSpectralEnergyDistribution::getnumberOfMaskedRegions(void) const {
    return numberOfMaskedRegions;
}

void operaSpectralEnergyDistribution::setMaskRegion(unsigned index, double Wavelength1, double Wavelength2) {
    if (index > MAXNUMBEROFMASKEDREGIONS) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    MaskWavelength1[index] = Wavelength1;
    MaskWavelength2[index] = Wavelength2;
}


double operaSpectralEnergyDistribution::getMaskRegionWavelength1(unsigned index) const {
    if (index > MAXNUMBEROFMASKEDREGIONS) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    return MaskWavelength1[index];
}

double operaSpectralEnergyDistribution::getMaskRegionWavelength2(unsigned index) const {
    if (index > MAXNUMBEROFMASKEDREGIONS) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
    return MaskWavelength2[index];
}


/*
 * Read wavelength mask regions from file
 */
void operaSpectralEnergyDistribution::readMaskedRegionsFromFile(string MaskedRegionsFileName) {
	ifstream astream;
	string dataline;
    
	double tmpwl1 = 0.0;
	double tmpwl2 = 0.0;
	unsigned np = 0;
	
	astream.open(MaskedRegionsFileName.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl1, &tmpwl2);
                    setMaskRegion(np,tmpwl1,tmpwl2);
                    np++;
                }	// skip comments
            }
		} // while (astream.good())
        astream.close();
	}	// if (astream.open())
    
	setnumberOfMaskedRegions(np);
}

bool operaSpectralEnergyDistribution::isItMasked(double Wavelength) {
    for(unsigned index=0;index<getnumberOfMaskedRegions();index++) {
        if(Wavelength > MaskWavelength1[index] && Wavelength < MaskWavelength2[index]) {
            return true;
        }
    }
    return false;
}


/*
 * Other Methods
 */

void operaSpectralEnergyDistribution::measureUncalibratedContinuum(operaSpectralElements &uncalibratedFluxElements, unsigned binsize, unsigned nsigcut) {

    unsigned NumberofPoints = uncalibratedFluxElements.getnSpectralElements();
    
	unsigned NumberOfSamples = (unsigned)ceil((float)NumberofPoints/(float)binsize); 
    
    float *uncalflux_tmp = new float[3*binsize];
    float *elemindex_tmp = new float[3*binsize];  
    
    float *residuals_tmp = new float[binsize];  
    float *iindex = new float[binsize];     
    int *sindex = new int[binsize];    
	
    float *continuumFluxSample = new float[NumberOfSamples + 2]; 
    float *continuumElemSample = new float[NumberOfSamples + 2];   

    createDataVectors(NumberOfSamples + 2);
    
    unsigned actualNumberOfSamples = 0;
    
	for(unsigned k=0;k<NumberOfSamples;k++){
        
        float am,bm,abdevm;
		
		unsigned firstPoint,lastPoint;
        
        if(k==0) {
            firstPoint=0;
            lastPoint=3*binsize;
        } else if (k>NumberOfSamples-3) {
            firstPoint = NumberofPoints - 3*binsize;
            lastPoint = NumberofPoints;
        } else {
            firstPoint = (k-1)*binsize;
            lastPoint = (k+2)*binsize;
        }        
		
        unsigned np=0;
        for(unsigned i=firstPoint;i<lastPoint;i++)
        {
            uncalflux_tmp[np] = (float)uncalibratedFluxElements.getFlux(i);
            elemindex_tmp[np] = (float)i;
            np++;
        }	
        
		// 1st pass: calculate continuum slope      
        ladfit(elemindex_tmp,uncalflux_tmp,np,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
		
        //--- Clean up spectrum leaving only points within the box that deviate less than abdev from the robust linear fit 
        np=0;
        for(unsigned i=firstPoint;i<lastPoint;i++)
        {
            float fitMedianSlope = (bm*(float)i + am);
            
            if(fabs((float)uncalibratedFluxElements.getFlux(i) - fitMedianSlope) < abdevm) {
                uncalflux_tmp[np] = (float)uncalibratedFluxElements.getFlux(i);
                elemindex_tmp[np] = (float)i;
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
        for(unsigned i=firstPoint;i<lastPoint;i++) {
            float fitMedianSlope = (bm*(float)i + am);
            residuals_tmp[np] = (float)uncalibratedFluxElements.getFlux(i) - fitMedianSlope;
            iindex[np] = (float)i;
            np++;
        }
        
        // sort residuals
        operaArrayIndexSort(np,residuals_tmp,sindex);
        
        // select first maximum residual which does not exceed 5x(absolute deviation)         
        float dytop = 0;
        float continuumBinFraction = 0.5;
        for(unsigned i=binsize-1;i>(unsigned)((float)binsize*(1.0 - continuumBinFraction));i--) {
            if(residuals_tmp[sindex[i]] < float(nsigcut)*abdevm) {
                dytop = residuals_tmp[sindex[i]];
                break;
            }
        }
        
        // if none has been selected then assign 3*abdev       
        if(dytop == 0) { 
            dytop = float(nsigcut)*abdevm; 
        }
        
        // evaluate max at left most edge
        if(k==0) {            
            continuumElemSample[actualNumberOfSamples] = (float)firstPoint;
            continuumFluxSample[actualNumberOfSamples] = bm*continuumElemSample[actualNumberOfSamples] + am + dytop;            
            
            distanceData[actualNumberOfSamples] = (double)uncalibratedFluxElements.getdistd((unsigned)continuumElemSample[actualNumberOfSamples]);
            wavelengthData[actualNumberOfSamples] = (double)uncalibratedFluxElements.getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]);
            fluxData->setflux(continuumFluxSample[actualNumberOfSamples],actualNumberOfSamples);
            spectralEnergyData->setflux(1.0,actualNumberOfSamples);              
            setwavelengthForNormalization(wavelengthData[actualNumberOfSamples]);	// DT Sep 6 2013 is this right?
            actualNumberOfSamples++;
        }
		
        if (k==NumberOfSamples-1) {
			continuumElemSample[actualNumberOfSamples] = (float)lastPoint-1;
			continuumFluxSample[actualNumberOfSamples] = bm*continuumElemSample[actualNumberOfSamples] + am + dytop;
            
            distanceData[actualNumberOfSamples] = (double)uncalibratedFluxElements.getdistd((unsigned)continuumElemSample[actualNumberOfSamples]);
            wavelengthData[actualNumberOfSamples] = (double)uncalibratedFluxElements.getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]);
            fluxData->setflux(continuumFluxSample[actualNumberOfSamples],actualNumberOfSamples);            
            spectralEnergyData->setflux(1.0,actualNumberOfSamples);            
			setwavelengthForNormalization(wavelengthData[actualNumberOfSamples]);	// DT Sep 6 2013 is this right?
            actualNumberOfSamples++; 
        } else {
            continuumElemSample[actualNumberOfSamples] = (float)firstPoint + (float)(binsize)/2.0;
            continuumFluxSample[actualNumberOfSamples] = bm*continuumElemSample[actualNumberOfSamples] + am + dytop;
            
            distanceData[actualNumberOfSamples] = (double)uncalibratedFluxElements.getdistd((unsigned)continuumElemSample[actualNumberOfSamples]);
            wavelengthData[actualNumberOfSamples] = (double)uncalibratedFluxElements.getwavelength((unsigned)continuumElemSample[actualNumberOfSamples]);
            fluxData->setflux(continuumFluxSample[actualNumberOfSamples],actualNumberOfSamples);
            spectralEnergyData->setflux(1.0,actualNumberOfSamples);              
            setwavelengthForNormalization(wavelengthData[actualNumberOfSamples]);	// DT Sep 6 2013 is this right?
            actualNumberOfSamples++; 
        }
	}
    
    if(actualNumberOfSamples <= nDataPoints){
        nDataPoints = actualNumberOfSamples;
        if(actualNumberOfSamples) {
            setHasFluxData(true);
        }
    } else {
        createDataVectors(actualNumberOfSamples);
        nDataPoints = actualNumberOfSamples;
        if(actualNumberOfSamples) {
            setHasFluxData(true);
        }        
    }
    
    delete[] continuumFluxSample; 
    delete[] continuumElemSample;
    
    delete[] uncalflux_tmp;
    delete[] elemindex_tmp;  
    delete[] residuals_tmp;
    delete[] iindex;
    delete[] sindex;      
}

void operaSpectralEnergyDistribution::calculateCalibratedElements(unsigned nPointsInReference, double *refwl, double *refflux) {

    if(!uncalibratedFluxElements || !calibratedFluxElements) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	        
    }
    
    float *wavelengthData_f = new float[nPointsInReference];
    float *spectralEnergyData_f = new float[nPointsInReference];
    
    for(unsigned i=0;i<nPointsInReference;i++) { 
        wavelengthData_f[i] = (float)refwl[i];
        spectralEnergyData_f[i] =  (float)refflux[i];        
    }    
    
    unsigned nElements = uncalibratedFluxElements->getnSpectralElements();
    
    float *SEDModelFlux = new float[nElements];    
    float *elemWavelength = new float[nElements];  
    
    for(unsigned i=0;i<nElements;i++) { 
        elemWavelength[i] =  (float)uncalibratedFluxElements->getwavelength(i);
    }    
    
    operaFitSpline(nPointsInReference,wavelengthData_f,spectralEnergyData_f,nElements,elemWavelength,SEDModelFlux);
    
    for(unsigned i=0;i<nElements;i++) {
        calibratedFluxElements->setwavelength((double)elemWavelength[i],i);
        calibratedFluxElements->setFlux((double)SEDModelFlux[i],i);
    }
    
    setHasCalibratedFlux(true);
    
    delete[] wavelengthData_f;
    delete[] spectralEnergyData_f;
    
    delete[] SEDModelFlux;
    delete[] elemWavelength;
}

void operaSpectralEnergyDistribution::calculateUncalibratedElements(unsigned binsize, unsigned nsigcut) {
    if(!uncalibratedFluxElements) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	        
    }
    
    measureUncalibratedContinuum(*uncalibratedFluxElements,binsize,nsigcut);
    
    float *wavelengthData_f = new float[nDataPoints];
    float *fluxData_f = new float[nDataPoints];
        
    for(unsigned i=0;i<nDataPoints;i++) { 
        wavelengthData_f[i] = (float)wavelengthData[i];
        fluxData_f[i] =  (float)fluxData->getflux(i);
    }
    
    unsigned nElements = uncalibratedFluxElements->getnSpectralElements();
    
    float *continuumModelFlux = new float[nElements];
    float *elemWavelength = new float[nElements];  
        
    for(unsigned i=0;i<nElements;i++) { 
        elemWavelength[i] =  uncalibratedFluxElements->getwavelength(i);
    }
        
    operaFitSpline(nDataPoints,wavelengthData_f,fluxData_f,nElements,elemWavelength,continuumModelFlux);
    
    for(unsigned i=0;i<nElements;i++) { 
        uncalibratedFluxElements->setFlux((double)continuumModelFlux[i],i);
//        uncalibratedFluxElements->setFluxVariance((double)continuumModelvariance[i],i);
#ifdef PRINT_DEBUG          
        cout << i << " " 
        << uncalibratedFluxElements->getwavelength(i) << " "
        << uncalibratedFluxElements->getFlux(i) << " " 
        << calibratedFluxElements->getFlux(i) << " "
        << endl;             
#endif            
    }        
    
    setHasUncalibratedFlux(true);
    
    delete[] continuumModelFlux;
    delete[] elemWavelength;  
    
    delete[] wavelengthData_f;
    delete[] fluxData_f;  
}

void operaSpectralEnergyDistribution::populateUncalibratedElementsFromContinuumData(void) {
    if(!uncalibratedFluxElements) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    if(!getHasFluxData()) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
        
    float *wavelengthData_f = new float[nDataPoints];
    float *fluxData_f = new float[nDataPoints];
    
    for(unsigned i=0;i<nDataPoints;i++) {
        wavelengthData_f[i] = (float)wavelengthData[i];
        fluxData_f[i] =  (float)fluxData->getflux(i);
    }
    
    unsigned nElements = uncalibratedFluxElements->getnSpectralElements();
    
    float *continuumModelFlux = new float[nElements];
    float *elemWavelength = new float[nElements];
    
    for(unsigned i=0;i<nElements;i++) {
        elemWavelength[i] =  uncalibratedFluxElements->getwavelength(i);
    }
    
    operaFitSpline(nDataPoints,wavelengthData_f,fluxData_f,nElements,elemWavelength,continuumModelFlux);
    
    for(unsigned i=0;i<nElements;i++) {
        uncalibratedFluxElements->setFlux((double)continuumModelFlux[i],i);
        //        uncalibratedFluxElements->setFluxVariance((double)continuumModelvariance[i],i);
#ifdef PRINT_DEBUG
        cout << order << " " <<  i << " "
        << uncalibratedFluxElements->getwavelength(i) << " "
        << uncalibratedFluxElements->getFlux(i) << " "
        << calibratedFluxElements->getFlux(i) << " "
        << endl;
#endif
    }

    setHasUncalibratedFlux(true);
    
    delete[] continuumModelFlux;
    delete[] elemWavelength;
    
    delete[] wavelengthData_f;
    delete[] fluxData_f;
}

void operaSpectralEnergyDistribution::calculateSEDelements(double spectralBinConstant,double referenceFluxForNormalization,double uncalibratedContinuumFluxForNormalization) {
    
    if(!hasUncalibratedFlux || !hasCalibratedFlux) {
		throw operaException("operaSpectralEnergyDistribution: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned nElements = uncalibratedFluxElements->getnSpectralElements();
        
    for(unsigned i=0;i<nElements;i++) {
        fluxCalibration->setwavelength(uncalibratedFluxElements->getwavelength(i),i);
        instrumentThroughput->setwavelength(uncalibratedFluxElements->getwavelength(i),i);
        
        double uncalibratedFlux = uncalibratedFluxElements->getFlux(i);
        double calibratedFlux = calibratedFluxElements->getFlux(i);
        
        if(calibratedFlux && spectralBinConstant && referenceFluxForNormalization && uncalibratedContinuumFluxForNormalization) {
            fluxCalibration->setFlux(uncalibratedFlux/(calibratedFlux*spectralBinConstant),i);
            instrumentThroughput->setFlux((uncalibratedFlux/uncalibratedContinuumFluxForNormalization)/(calibratedFlux/referenceFluxForNormalization),i);
        } else {
            fluxCalibration->setFlux(NAN,i);
            instrumentThroughput->setFlux(NAN,i);
        }
        
#ifdef PRINT_DEBUG
        cout << i << " "
        << uncalibratedFluxElements->getwavelength(i) << " "
        << uncalibratedFluxElements->getFlux(i) << " "
        << calibratedFluxElements->getFlux(i) << " "
        << fluxCalibration->getFlux(i) << " "
        << instrumentThroughput->getFlux(i) << " "
        << endl;
#endif
    }
    setHasFluxCalibration(true);
    setHasInstrumentThroughput(true);
}

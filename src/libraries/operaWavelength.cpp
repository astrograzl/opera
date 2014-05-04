/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaWavelength.h"
#include "libraries/Gaussian.h"
#include "libraries/operaFit.h"	 // for operaMPFitPolynomial and operaLMFitPolynomial
#include "libraries/operaMath.h"	 // for DiffPolynomialFunction

#define NPOINTPERSIGMA 4

/*!
 * operaWavelength
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the wavelength object.
 * \file operaWavelength.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaWavelength
 * \brief Encapsulation of Wavelength information.
 * \return none
 */
/*
 * Constructors
 */
operaWavelength::operaWavelength() :
dmin(0.0),
dmax(0.0),
maxNDataPoints(0),
nDataPoints(0),
distanceData(NULL),
wavelengthData(NULL),
wavelengthErrors(NULL),
matchAtlasindex(NULL),
matchComparisonindex(NULL),
maxNAtlasLines(0),
nAtlasLines(0),
atlasLinesflux(NULL),
atlasLineswl(NULL),
atlasLineswlError(NULL),
maxNComparisonLines(0), 
nComparisonLines(0), 
comparisonLinesflux(NULL),
comparisonLinespix(NULL),
comparisonLinespixError(NULL),
comparisonLineswl(NULL),
radialVelocityPrecision(0.0),
xcorrelation(0.0),
wavelengthPolynomial(NULL)	   
{
	wavelengthPolynomial = new Polynomial();
}

operaWavelength::operaWavelength(unsigned Coeffs) :
dmin(0.0),
dmax(0.0),
maxNDataPoints(0),
nDataPoints(0),
distanceData(NULL),
wavelengthData(NULL),
wavelengthErrors(NULL),
matchAtlasindex(NULL),
matchComparisonindex(NULL),
maxNAtlasLines(0),
nAtlasLines(0),
atlasLinesflux(NULL),
atlasLineswl(NULL),
atlasLineswlError(NULL),
maxNComparisonLines(0), 
nComparisonLines(0), 
comparisonLinesflux(NULL),
comparisonLinespix(NULL),
comparisonLinespixError(NULL),
comparisonLineswl(NULL),
radialVelocityPrecision(0.0),
xcorrelation(0.0),
wavelengthPolynomial(NULL)	   
{
	wavelengthPolynomial = new Polynomial(Coeffs);
}

operaWavelength::~operaWavelength() {
	deleteDataVectors();
	deleteAtlasDataVectors();
	deleteComparisonDataVectors();
	//if (wavelengthPolynomial != NULL)
	//	delete wavelengthPolynomial;    
	wavelengthPolynomial = NULL;
}

/*
 * Methods
 */

double operaWavelength::getDmin(void) {
	return dmin;
}

double operaWavelength::getDmax(void) {
	return dmax;
}

void operaWavelength::setDmin(double Dmin) {
	dmin = Dmin;
}

void operaWavelength::setDmax(double Dmax) {
	dmax = Dmax;
}

double operaWavelength::getDistance(unsigned index) {
	if (index >= maxNDataPoints) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return distanceData[index];
}

double operaWavelength::getWavelength(unsigned index) {
	if (index >= maxNDataPoints) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return wavelengthData[index];
}


double operaWavelength::getWavelengthError(unsigned index) {
	if (index >= maxNDataPoints) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return wavelengthErrors[index];
}

unsigned operaWavelength::getMatchAtlasIndex(unsigned index) {
	if (index >= maxNDataPoints) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	return matchAtlasindex[index];
}

unsigned operaWavelength::getMatchComparisonIndex(unsigned index) {
	if (index >= maxNDataPoints) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	return matchComparisonindex[index];
}

double operaWavelength::getcentralWavelength(void) {
	return wavelengthPolynomial->Evaluate(dmin + (dmax - dmin)/2);
}

double operaWavelength::getinitialWavelength(void) {
	return wavelengthPolynomial->Evaluate(dmin);
}

double operaWavelength::getfinalWavelength(void) {
	return wavelengthPolynomial->Evaluate(dmax);
}

Polynomial *operaWavelength::getWavelengthPolynomial(void) {
	return wavelengthPolynomial;
}

void operaWavelength::setnDataPoints(unsigned NDataPoints) {
	if (NDataPoints == 0) {
		// DT Apr 2013, causes wcals to throw an exception when ndatapoints is zero...
		//throw operaException("operaWavelength: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (NDataPoints > maxNDataPoints) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    nDataPoints = NDataPoints;
}

unsigned operaWavelength::getnDataPoints(void) {
    return nDataPoints;
}

void operaWavelength::createDataVectors(unsigned NDataPoints) {
	if (NDataPoints == 0) {
		throw operaException("operaWavelength: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	distanceData = (double *)malloc(sizeof(double)*NDataPoints);  
	wavelengthData = (double *)malloc(sizeof(double)*NDataPoints);      
	wavelengthErrors = (double *)malloc(sizeof(double)*NDataPoints);        
	matchAtlasindex = (unsigned *)malloc(sizeof(unsigned)*NDataPoints); 
	matchComparisonindex = (unsigned *)malloc(sizeof(unsigned)*NDataPoints); 
    maxNDataPoints = NDataPoints;
}

void operaWavelength::deleteDataVectors(void) {
    if(distanceData) {
        free(distanceData);
        distanceData = NULL;
    }
    if(wavelengthData) {
        free(wavelengthData);  
        wavelengthData = NULL;
    }
    if(wavelengthErrors) {
        free(wavelengthErrors); 
        wavelengthErrors = NULL;  
    }
    if(matchAtlasindex) {
        free(matchAtlasindex);  
        matchAtlasindex = NULL; 
    }
    if(matchComparisonindex) {
        free(matchComparisonindex); 
        matchComparisonindex = NULL;
    }
    maxNDataPoints = 0;
    nDataPoints = 0;
}

void operaWavelength::CalculateWavelengthSolution(unsigned maxcoeffs, bool witherrors) {
    
	doublePolynomialCoeffs_t coeffs;
	int nparbestfit = maxcoeffs;
	double bestchisqr = BIG;
    double chisqr = BIG;
    
	for (unsigned currentfit=1; currentfit<=maxcoeffs; currentfit++) {
        if (wavelengthPolynomial) {
            int npar = wavelengthPolynomial->getOrderOfPolynomial();
            double *currentpar = (double *)wavelengthPolynomial->getVector();
            double *currenterrs = (double *)wavelengthPolynomial->getErrorVector();
            
            for(unsigned i=0;i<currentfit;i++) {
                if(i < (unsigned)npar) {
                    coeffs.p[i] = currentpar[i];
                    coeffs.e[i] = currenterrs[i];
                } 
            }       
        } else {
            for	(unsigned i=0; i<maxcoeffs; i++) {
                coeffs.p[i] = 1.0;
                coeffs.e[i] = 0.0;
            }
        }
		
		if (witherrors) {
			int errorcode = operaMPFitPolynomial(nDataPoints, distanceData, wavelengthData, wavelengthErrors, currentfit, coeffs.p, coeffs.e, &chisqr);
			if (errorcode <= 0) {
				throw operaException("operaWavelength: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);	
			}	
		} else {
			operaLMFitPolynomial(nDataPoints, distanceData, wavelengthData, currentfit, coeffs.p, &chisqr);
		}
		if (chisqr < bestchisqr) {
			bestchisqr = chisqr;
			nparbestfit = currentfit;
		}
	}
	
    if (wavelengthPolynomial) {
        int npar = wavelengthPolynomial->getOrderOfPolynomial();
        double *currentpar = (double *)wavelengthPolynomial->getVector();
        double *currenterrs = (double *)wavelengthPolynomial->getErrorVector();
        
        for(unsigned i=0;i<(unsigned)nparbestfit;i++) {
            if(i < (unsigned)npar) {
				coeffs.p[i] = currentpar[i];
				coeffs.e[i] = currenterrs[i];
            }
        }       
    } else {
        for	(unsigned i=0; i<(unsigned)nparbestfit; i++) {
			coeffs.p[i] = 1.0;
			coeffs.e[i] = 0.0;
        }
    }    
	if (witherrors) {
		int errorcode = operaMPFitPolynomial(nDataPoints, distanceData, wavelengthData, wavelengthErrors, nparbestfit, coeffs.p, coeffs.e, &chisqr);
		if (errorcode <= 0) {
			throw operaException("operaWavelength: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);	
		}
	} else {
		operaLMFitPolynomial(nDataPoints, distanceData, wavelengthData, nparbestfit, coeffs.p, &chisqr);
	}		
	coeffs.orderofPolynomial = nparbestfit;
	coeffs.polychisqr = chisqr;
	PolynomialCoeffs_t fcoeffs;
	PolynomialCoeffsToFloat(&fcoeffs, &coeffs);
	// wants floats...
	wavelengthPolynomial->setPolynomialCoeffs(&fcoeffs);    
}

void operaWavelength::RefineWavelengthSolution(unsigned ncoeffs, bool witherrors) {
    if (!wavelengthPolynomial) {
        throw operaException("operaWavelength: no polynomial found: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);	        
    }
	
    int npar = wavelengthPolynomial->getOrderOfPolynomial();
    
	doublePolynomialCoeffs_t coeffs;
    
	double *currentpar = (double *)wavelengthPolynomial->getVector();
    double *currenterrs = (double *)wavelengthPolynomial->getErrorVector();
	
    for(unsigned i=0;i<(unsigned)npar;i++) {
		coeffs.p[i] = currentpar[i];
		coeffs.e[i] = currenterrs[i];
    } 
    
    if(ncoeffs > (unsigned)npar) {
        for(unsigned i=(unsigned)npar;i<ncoeffs;i++) {
			coeffs.p[i] = 1.0;
			coeffs.e[i] = 0.0;
        } 
    }
    
	double currentchisqr = wavelengthPolynomial->getChisqr();
    double newchisqr = 0.0;
    
    if (witherrors) {
        int errorcode = operaMPFitPolynomial(nDataPoints, distanceData, wavelengthData, wavelengthErrors, ncoeffs, coeffs.p, coeffs.e, &newchisqr);
        if (errorcode <= 0) {
            throw operaException("operaWavelength: ", operaErrorGeometryBadFit, __FILE__, __FUNCTION__, __LINE__);	
        }	
    } else {
        operaLMFitPolynomial(nDataPoints, distanceData, wavelengthData, ncoeffs, coeffs.p, &newchisqr);
    }
    if (newchisqr < currentchisqr) {
		coeffs.orderofPolynomial = ncoeffs;
		coeffs.polychisqr = newchisqr;
		PolynomialCoeffs_t fcoeffs;
		PolynomialCoeffsToFloat(&fcoeffs, &coeffs);
		// wants floats...
		wavelengthPolynomial->setPolynomialCoeffs(&fcoeffs);    
    }
}


void operaWavelength::setSpectralResolution(doubleValue_t Resolution) {
    spectralResolution = Resolution;
}

doubleValue_t operaWavelength::getSpectralResolution(void) {
    return spectralResolution;
}


/*
 * Function to calculate spectral resolution
 */
void operaWavelength::calculateSpectralResolution(doubleValue_t ResolutionElementInPixels) {
    
    double deltawl = convertPixelToWavelength(ResolutionElementInPixels.value);
    double errorwl = convertPixelToWavelength(ResolutionElementInPixels.error);
    
    spectralResolution.value = getcentralWavelength()/deltawl;   
    spectralResolution.error = errorwl * getcentralWavelength()/(deltawl*deltawl);     
}

/*
 * Function to evaluate the wavelength (in nm) correspondent to a given distance value (in pixels)
 */
double operaWavelength::evaluateWavelength(double distanceValue) {
    return wavelengthPolynomial->Evaluate(distanceValue);
}

/*
 * Function to evaluate the wavelength (in nm) correspondent to a given distance value (in pixels)
 */
double operaWavelength::convertPixelToWavelength(double DeltaDistanceInPixels) {
	int npar = wavelengthPolynomial->getOrderOfPolynomial();
	double *par = (double *)wavelengthPolynomial->getVector();
    double deltawl = DeltaDistanceInPixels * (double)DiffPolynomialFunction(getcentralWavelength(),par,npar);
    return deltawl;
}


void operaWavelength::setRadialVelocityPrecision(double radialvelocityprecision) {
    radialVelocityPrecision = radialvelocityprecision;
}

double operaWavelength::getRadialVelocityPrecision(void) {
    return radialVelocityPrecision;
}

double operaWavelength::getxcorrelation(void) {
    return xcorrelation;
}

void operaWavelength::setxcorrelation(double Xcorrelation) {
    xcorrelation = Xcorrelation;
}

/*
 * Function to calculate radial velocity precision in m/s
 */
void operaWavelength::calculateRadialVelocityPrecision(void) {
    
    double speedoflight = 299792458; // in m/s    
    double radialvelocity_residuals = 0;
    
    for(unsigned i=0;i<nDataPoints;i++) {
		radialvelocity_residuals += (speedoflight*(wavelengthData[i] - evaluateWavelength(distanceData[i]))/evaluateWavelength(distanceData[i]))*(speedoflight*(wavelengthData[i] - evaluateWavelength(distanceData[i]))/evaluateWavelength(distanceData[i]));
    }
	
    radialVelocityPrecision = sqrt(radialvelocity_residuals/(double)nDataPoints);
}

/*
 * Function to calculate radial velocity precision in m/s
 */
double operaWavelength::calculateWavelengthRMSPrecision(void) {
    
    double rms_of_wlresiduals = 0;
    
    for(unsigned i=0;i<nDataPoints;i++) {
        rms_of_wlresiduals += (wavelengthData[i] - evaluateWavelength(distanceData[i]))*(wavelengthData[i] - evaluateWavelength(distanceData[i]));
    }
    
    rms_of_wlresiduals = sqrt(rms_of_wlresiduals/(double)nDataPoints);
    
    return rms_of_wlresiduals;
}

/*
 * Function to calculate radial velocity precision in m/s
 */
double operaWavelength::calculateWavelengthMedianPrecision(void) {
    
    float *wlresiduals = new float[nDataPoints];
    
    for(unsigned i=0;i<nDataPoints;i++) {
        wlresiduals[i] = (float)fabs(wavelengthData[i] - evaluateWavelength(distanceData[i]));
    }
    
    double MedianPrecision = (double)operaArrayMedianSigmaQuick(nDataPoints,wlresiduals,operaArrayMedianQuick(nDataPoints,wlresiduals));
    
    delete[] wlresiduals;
    
    return MedianPrecision;
}

/*
 * Atlas wavelength data
 */

void operaWavelength::setnAtlasLines(unsigned NAtlasLines) {
    maxNAtlasLines = NAtlasLines;
    nAtlasLines = NAtlasLines;
}

unsigned operaWavelength::getnAtlasLines(void) {
    return nAtlasLines;
}    
double operaWavelength::getatlasLinesflux(unsigned index) {
	if (index > nAtlasLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return atlasLinesflux[index];
} 

double operaWavelength::getatlasLineswl(unsigned index) {
	if (index > nAtlasLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return atlasLineswl[index];
}

double operaWavelength::getatlasLineswlError(unsigned index) {
	if (index > nAtlasLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return atlasLineswlError[index];
}

void operaWavelength::setatlasLinesflux(double AtlasLinesflux, unsigned index) {
	if (index > nAtlasLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    atlasLinesflux[index] = AtlasLinesflux;
}

void operaWavelength::setatlasLineswl(double AtlasLineswl, unsigned index) {
	if (index > nAtlasLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    atlasLineswl[index] = AtlasLineswl;   
}

void operaWavelength::setatlasLineswl(double AtlasLineswl, double AtlasLineswlError, unsigned index) {
	if (index > nAtlasLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    atlasLineswl[index] = AtlasLineswl;
    atlasLineswlError[index] = AtlasLineswlError;    
}

void operaWavelength::createAtlasDataVectors(unsigned NAtlasLines, double *AtlasLineswl, double *AtlasLineswlError, double *AtlasLinesflux) {
    nAtlasLines = NAtlasLines; 
    maxNAtlasLines = NAtlasLines; 
    if(nAtlasLines) {
        atlasLinesflux = (double *)malloc(sizeof(double)*nAtlasLines);  
        atlasLineswl = (double *)malloc(sizeof(double)*nAtlasLines);  
        atlasLineswlError = (double *)malloc(sizeof(double)*nAtlasLines);         
        for(unsigned i=0; i< nAtlasLines; i++ ) {
            atlasLinesflux[i] = AtlasLinesflux[i];
            atlasLineswl[i] = AtlasLineswl[i];
            atlasLineswlError[i] = AtlasLineswlError[i];            
        }
    }
}

void operaWavelength::createAtlasDataVectors(unsigned NAtlasLines) {
    nAtlasLines = NAtlasLines; 
    maxNAtlasLines = NAtlasLines; 
    if(nAtlasLines) {
        atlasLinesflux = (double *)malloc(sizeof(double)*nAtlasLines);  
        atlasLineswl = (double *)malloc(sizeof(double)*nAtlasLines);  
        atlasLineswlError = (double *)malloc(sizeof(double)*nAtlasLines);        
    }
}

void operaWavelength::deleteAtlasDataVectors(void) {
    if(atlasLinesflux)
        free(atlasLinesflux);
	atlasLinesflux = NULL;
    if(atlasLineswl)
        free(atlasLineswl);
	atlasLineswl = NULL;
    if(atlasLineswlError)
        free(atlasLineswlError); 
	atlasLineswlError = NULL;
    nAtlasLines = 0;
    maxNAtlasLines = 0;
}

/*
 *  Comparison pixel data
 */    

void operaWavelength::setnComparisonLines(unsigned NComparisonLines) {
    nComparisonLines = NComparisonLines;
    maxNComparisonLines = NComparisonLines;
}    

unsigned operaWavelength::getnComparisonLines(void) {
    return nComparisonLines;
}  

double operaWavelength::getcomparisonLinesflux(unsigned index) {
	if (index > nComparisonLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return comparisonLinesflux[index];
} 

double operaWavelength::getcomparisonLinespix(unsigned index) {
	if (index > nComparisonLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return comparisonLinespix[index];
} 

double operaWavelength::getcomparisonLinespixError(unsigned index) {
 	if (index > nComparisonLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return comparisonLinespixError[index];
} 

double operaWavelength::getcomparisonLineswl(unsigned index) {
	if (index > nComparisonLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    return comparisonLineswl[index];
} 

void operaWavelength::setcomparisonLinesflux(double ComparisonLinesflux, unsigned index) {
	if (index > nComparisonLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    comparisonLinesflux[index] = ComparisonLinesflux;
}

void operaWavelength::setcomparisonLinespix(double ComparisonLinespix, double ComparisonLinespixError, unsigned index) {
	if (index > nComparisonLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    comparisonLinespix[index] = ComparisonLinespix;
    comparisonLinespixError[index] = ComparisonLinespixError;    
}

void operaWavelength::setcomparisonLineswl(double ComparisonLineswl, unsigned index) {
	if (index > nComparisonLines) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    comparisonLineswl[index] = ComparisonLineswl;
}

void operaWavelength::createComparisonDataVectors(unsigned NComparisonLines) {
    nComparisonLines = NComparisonLines; 
    maxNComparisonLines = NComparisonLines; 
    if(nComparisonLines) {
        comparisonLinesflux = (double *)malloc(sizeof(double)*nComparisonLines);  
        comparisonLinespix = (double *)malloc(sizeof(double)*nComparisonLines); 
        comparisonLinespixError = (double *)malloc(sizeof(double)*nComparisonLines);        
        comparisonLineswl = (double *)malloc(sizeof(double)*nComparisonLines);         
    }
}

void operaWavelength::createComparisonDataVectors(unsigned NComparisonLines, double *ComparisonLinespix, double *ComparisonLinespixError, double *ComparisonLinesflux){
    nComparisonLines = NComparisonLines; 
    maxNComparisonLines = NComparisonLines; 
    if(nComparisonLines) {
        comparisonLinesflux = (double *)malloc(sizeof(double)*nComparisonLines);  
        comparisonLinespix = (double *)malloc(sizeof(double)*nComparisonLines); 
        comparisonLinespixError = (double *)malloc(sizeof(double)*nComparisonLines);         
        comparisonLineswl = (double *)malloc(sizeof(double)*nComparisonLines);   
        for(unsigned i=0; i< nComparisonLines; i++ ) {
            comparisonLinesflux[i] = ComparisonLinesflux[i];
            comparisonLinespix[i] = ComparisonLinespix[i];
            comparisonLinespixError[i] = ComparisonLinespixError[i];            
            comparisonLineswl[i] = evaluateWavelength(comparisonLinespix[i]);            
        }        
    }
}

void operaWavelength::recalculateComparisonLineswlVector(void) {
    for(unsigned i=0; i< nComparisonLines; i++ ) {
        comparisonLineswl[i] = evaluateWavelength(comparisonLinespix[i]);            
    }        
}

void operaWavelength::deleteComparisonDataVectors(void) {
    if(comparisonLinesflux)
        free(comparisonLinesflux);
	comparisonLinesflux = NULL;
    if(comparisonLinespix)
        free(comparisonLinespix);  
	comparisonLinespix = NULL;
    if(comparisonLinespixError)
        free(comparisonLinespixError);    
	comparisonLinespixError = NULL;
    if(comparisonLineswl)
        free(comparisonLineswl);      
	comparisonLineswl = NULL;
    nComparisonLines = 0;
    maxNComparisonLines = 0;
}

void operaWavelength::refineWavelengthSolutionByXCorrelation(unsigned nPointsPerParameter, double parameterRangetoSearch, unsigned maxpolyorder) {
	
    double simulwl[MAXPOINTSINSIMULATEDSPECTRUM]; 
    double atlasSimulSpectrum[MAXPOINTSINSIMULATEDSPECTRUM]; 
    double comparisonSimulSpectrum[MAXPOINTSINSIMULATEDSPECTRUM]; 
    unsigned NSimulatedPoints;
	
    calculateXCorrelation();
	
    double Maxcorrelation = getxcorrelation();     
    
    int npar = wavelengthPolynomial->getOrderOfPolynomial();
    double *par = (double *)wavelengthPolynomial->getVector();   
    if(maxpolyorder > (unsigned)npar) {
        maxpolyorder = (unsigned)npar;
    }
    /*
     * Below it attempts to find the coefficients that gives the highest correlation 
     * between the raw and atlas simulated spectra.
     */    
    for(unsigned k=0;k<maxpolyorder;k++) {    
        /*
         * Note that nPointsPerParameter and parameterRangetoSearch are input parameters
         * that determine the step and range for which the coefficients will be searched
         */                           
        double dpar = fabs(par[k] * parameterRangetoSearch/100.0)/double(nPointsPerParameter);
        double par0 = par[k] -  fabs(par[k] * parameterRangetoSearch/100.0)/2.0;
        
        double maxpar0 = par[k];
        par[k] = par0;
        
        for(unsigned l=0; l<nPointsPerParameter; l++) {
            recalculateComparisonLineswlVector();
			
            NSimulatedPoints = createAtlasSimulatedSpectrum(simulwl, atlasSimulSpectrum, NPOINTPERSIGMA);
			
            NSimulatedPoints = createComparisonSimulatedSpectrum(simulwl, comparisonSimulSpectrum, NPOINTPERSIGMA);  
			
            double crosscorrelation = operaCrossCorrelation(NSimulatedPoints,atlasSimulSpectrum,comparisonSimulSpectrum);                       
			
            if(crosscorrelation > Maxcorrelation) {
                Maxcorrelation = crosscorrelation;
                maxpar0 = par[k];
            }
            par[k] += dpar;
        }
        par[k] = maxpar0;
    }
    recalculateComparisonLineswlVector();
    
    setxcorrelation(Maxcorrelation);
}

void operaWavelength::refineWavelengthSolutionByFindingMaxMatching(unsigned NpointsPerPar, double ParRangeSizeInPerCent, double acceptableMismatch) {
    double maxpercentage = 0;
    double *par = (double *)(getWavelengthPolynomial()->getVector());
    double parcentre = par[0];
    double dpar = fabs(par[0] * ParRangeSizeInPerCent/100.0)/double(NpointsPerPar);
    double par0 = par[0] -  fabs(par[0] * ParRangeSizeInPerCent/100.0)/2.0;
    double maxpar0 = par[0];
    //double dlambdamax = 0;
    par[0] = par0;
    
    for(unsigned i = 0; i < NpointsPerPar; i++) {
        matchAtlaswithComparisonLines(acceptableMismatch);
        //double dlambda = par[0] - parcentre;
        double MatchPercentage = (getPerCentageOfComparisonMatch() + getPerCentageOfAtlasMatch())/2;
        if(MatchPercentage > maxpercentage) {
            maxpercentage = MatchPercentage;
            //dlambdamax = dlambda;
            maxpar0 = par[0];
        }
        
        //    cout << order << " " <<  getWavelengthPolynomial()->getVector()[0] << " " << MatchPercentage << " " << dlambda << endl;
        par[0] += dpar;
    }
    par[0] = maxpar0;
}


void operaWavelength::calculateXCorrelation(void) {
    
    double simulwl[MAXPOINTSINSIMULATEDSPECTRUM];
    double atlasSimulSpectrum[MAXPOINTSINSIMULATEDSPECTRUM];
    double comparisonSimulSpectrum [MAXPOINTSINSIMULATEDSPECTRUM];
    
    unsigned NSimulatedPoints;
    NSimulatedPoints = createAtlasSimulatedSpectrum(simulwl, atlasSimulSpectrum, NPOINTPERSIGMA);
    NSimulatedPoints = createComparisonSimulatedSpectrum(simulwl, comparisonSimulSpectrum, NPOINTPERSIGMA);
    double crosscorrelation = operaCrossCorrelation(NSimulatedPoints,atlasSimulSpectrum,comparisonSimulSpectrum);                       
    setxcorrelation(crosscorrelation);
}
/*
 * Function to create an Atlas simulated spectrum
 */
unsigned operaWavelength::createAtlasSimulatedSpectrum(double *outputwl, double *outputSpectrum, unsigned nstepspersigma) {
    unsigned nLines = nAtlasLines;
    double *lineCenters = atlasLineswl;
    double *lineAmplitudes = atlasLinesflux;
    
    double wlc = getcentralWavelength();
    double minwl = getinitialWavelength();
    double maxwl = getfinalWavelength();    
    
    double lineSigma = wlc/(spectralResolution.value);
    double wlstep = lineSigma/(double)nstepspersigma;
    
    unsigned npoints = (unsigned)ceil(fabs(maxwl - minwl)/wlstep);
    
    double *sigmaVector = new double[nLines];
    
    for(unsigned i=0;i<nLines;i++) {
        sigmaVector[i] = lineSigma;
    }
    
    Gaussian spectrumModel(nLines,lineAmplitudes,sigmaVector,lineCenters);
	
    double wl = minwl;
    
    for(unsigned i=0;i<npoints;i++) {
        outputSpectrum[i] = spectrumModel.EvaluateGaussian(wl);
        outputwl[i] = wl;
        wl += wlstep;
    }
    
    delete[] sigmaVector;
    
    return npoints;
}

/*
 * Function to create a comparison simulated spectrum
 */
unsigned operaWavelength::createComparisonSimulatedSpectrum(double *outputwl, double *outputSpectrum, unsigned nstepspersigma) {
    unsigned nLines = nComparisonLines;
    
	if (nLines >= MAXPOINTSINSIMULATEDSPECTRUM) {
		throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    double *lineAmplitudes = comparisonLinesflux;
    double wlc = getcentralWavelength();
    double minwl = getinitialWavelength();
    double maxwl = getfinalWavelength();
    double lineSigma = wlc/(spectralResolution.value);
    double wlstep = lineSigma/(double)nstepspersigma;
    unsigned npoints = (unsigned)ceil(fabs(maxwl - minwl)/wlstep);
    double *sigmaVector = new double[nLines];
    for(unsigned i=0;i<nLines;i++) {
        sigmaVector[i] = lineSigma;
    }
    recalculateComparisonLineswlVector();
    Gaussian spectrumModel(nLines,lineAmplitudes,sigmaVector,comparisonLineswl);
    double wl = minwl;
    for(unsigned i=0;i<npoints;i++) {
        outputSpectrum[i] = spectrumModel.EvaluateGaussian(wl);
        outputwl[i] = wl;
        wl += wlstep;
    }
    delete[] sigmaVector;
    return npoints;
}


void operaWavelength::matchAtlaswithComparisonLines(double acceptableMismatch) {
    
    double wlc = getcentralWavelength();
    double lineSigma = wlc/(spectralResolution.value);    
    
    double acceptMismatchInwlUnits = acceptableMismatch*lineSigma;
    
    recalculateComparisonLineswlVector();
    
    /*
     * Below it identifies and select the set of lines that match both comparison and atlas.
     * The criteria for matching is that the difference between centers must be < acceptableMismatch x sigma 
     */          
    unsigned nmatch = 0;
    unsigned nextfirstline = 0;
	
    for (unsigned i=0; i<getnComparisonLines(); i++) {
        
        unsigned bestAtlasMatchIndex = 0;
        double mindifference = acceptMismatchInwlUnits;
        
        for(unsigned l=nextfirstline;l<getnAtlasLines();l++) {
            
            double difference = fabs(comparisonLineswl[i] - atlasLineswl[l]);
#ifdef PRINT_DEBUG	    
            cout << "mindiff=" << mindifference << " diff=" << difference << " comp[" << i << "]=" << comparisonLineswl[i] << " atlas[" << l << "]=" << atlasLineswl[l] << endl; 
#endif               
            
            if(comparisonLineswl[i] > atlasLineswl[l]  && difference < mindifference) {
                mindifference = difference;
                bestAtlasMatchIndex = l;
            } else if (comparisonLineswl[i] <= atlasLineswl[l] && difference < mindifference) {
                distanceData[nmatch] = comparisonLinespix[i];
                wavelengthData[nmatch] = atlasLineswl[l];
                wavelengthErrors[nmatch] = atlasLineswlError[l];
                matchAtlasindex[nmatch] = l;
                matchComparisonindex[nmatch] = i;
                nextfirstline = l+1; 
                nmatch++;
                break;
            } else if (comparisonLineswl[i] <= atlasLineswl[l]  && difference > mindifference) {
                if(bestAtlasMatchIndex) {
                    distanceData[nmatch] = comparisonLinespix[i];
                    wavelengthData[nmatch] = atlasLineswl[bestAtlasMatchIndex];
                    wavelengthErrors[nmatch] = atlasLineswlError[bestAtlasMatchIndex];
                    matchAtlasindex[nmatch] = bestAtlasMatchIndex;
                    matchComparisonindex[nmatch] = i;
                    nextfirstline = bestAtlasMatchIndex+1;                    
                    nmatch++;                     
                }
                break;
            }
            if(l==nAtlasLines-1) {
                if(bestAtlasMatchIndex) {
                    distanceData[nmatch] = comparisonLinespix[i];
                    wavelengthData[nmatch] = atlasLineswl[bestAtlasMatchIndex];
                    wavelengthErrors[nmatch] = atlasLineswlError[bestAtlasMatchIndex];
                    matchAtlasindex[nmatch] = bestAtlasMatchIndex;
                    matchComparisonindex[nmatch] = i;
                    nextfirstline = bestAtlasMatchIndex+1;                    
                    nmatch++;                     
                }  
            }
			if (nmatch >= getnAtlasLines()) {
				throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (nmatch >= getnComparisonLines()) {
				throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
			}
			if (nmatch >= maxNDataPoints) {
				throw operaException("operaWavelength: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
			}
        }
    }
    setnDataPoints(nmatch);
}

double operaWavelength::getPerCentageOfComparisonMatch(void) {
    return 100*(double)getnDataPoints()/(double)nComparisonLines;
}

double operaWavelength::getPerCentageOfAtlasMatch(void) {
    return 100*(double)getnDataPoints()/(double)nAtlasLines;
}

void operaWavelength::filterDataPointsBySigmaClip(double nsig) {
    double sigmaclip = (nsig)*calculateWavelengthRMSPrecision();
    unsigned np = 0;
    for(unsigned i=0;i<getnDataPoints();i++) {
        if(fabs(wavelengthData[i] - evaluateWavelength(distanceData[i])) < sigmaclip) {
            distanceData[np] = distanceData[i];
            wavelengthData[np] = wavelengthData[i];
            wavelengthErrors[np] = wavelengthErrors[i];
            matchAtlasindex[np] = matchAtlasindex[i];
            matchComparisonindex[np] = matchComparisonindex[i];
            np++; 
        }
    }
    setnDataPoints(np);
}

void operaWavelength::filterDataPointsByErrorClip(double nsig) {
    double sigmaclip = (nsig)*calculateWavelengthRMSPrecision();
    unsigned np = 0;
    for(unsigned i=0;i<getnDataPoints();i++) {
        if(wavelengthErrors[i] < sigmaclip) {
            distanceData[np] = distanceData[i];
            wavelengthData[np] = wavelengthData[i];
            wavelengthErrors[np] = wavelengthErrors[i];
            matchAtlasindex[np] = matchAtlasindex[i];
            matchComparisonindex[np] = matchComparisonindex[i];
            np++; 
        }
    }
    setnDataPoints(np);
}

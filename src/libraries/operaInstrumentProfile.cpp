/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaInstrumentProfile
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
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaInstrumentProfile.h"
#include "libraries/operaException.h"
#include "libraries/operaExtractionAperture.h"
#include "libraries/PixelSet.h"

#include "libraries/operaFit.h"
#include "libraries/operaMath.h"
#include "libraries/operaStats.h"

#ifndef SATURATIONLIMIT
#define SATURATIONLIMIT 65535  // this should be retrieved from the config/param file
#endif

/*! 
 * operaInstrumentProfile
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the Instrument Profile.
 * \file operaInstrumentProfile.cpp
 * \ingroup libraries
 */

using namespace std;

/*! 
 * \brief operaInstrumentProfile
 * \details The instrument profile (IP) consists of a data set (pixelized image) 
 * \details representing the distribution of the fraction of flux for a uniformly illuminated 
 * \details monochromatic image of the entrance slit projected on the spectrograph focal plane. 
 * \details The fiducial set of coordinates chosen is such that the ordinate is oriented along 
 * \details with the dispersion direction and the abscissa with the spatial direction.
 */

/*
 * Constructor
 */

operaInstrumentProfile::operaInstrumentProfile(void) :
NXPoints(1), // number of points in x direction
NYPoints(1), // number of points in y directon
NTotalPoints(1), // total number of points
Xsampling(1),		// (sampling factor for the PSF in x direction)
Ysampling(1),		// (sampling factor for the  PSF in y direction)
xsize(1),			// (size along the x direction in pix units) 
ysize(1),			// (size along the y direction in pix units)
geometricCenterX(0.5),
geometricCenterY(0.5),
nDataPoints(0),
maxnDataPoints(0),
dataCube(NULL),
errorsCube(NULL),
distd(NULL)
{
	ipPolyModel = newPolynomialMatrix(NXPoints, NYPoints);
	chisqrMatrix = newCMatrix(NXPoints, NYPoints);
}

operaInstrumentProfile::operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling) :
NXPoints(1),		// number of points in x direction
NYPoints(1),		// number of points in y directon
NTotalPoints(1),	// total number of points
Xsampling(1),		// (sampling factor for the PSF in x direction)
Ysampling(1),		// (sampling factor for the  PSF in y direction)
xsize(1),			// (size along the x direction in pix units) 
ysize(1),			// (size along the y direction in pix units)
geometricCenterX(0.5),
geometricCenterY(0.5),
nDataPoints(0),
maxnDataPoints(0),
dataCube(NULL),
errorsCube(NULL),
distd(NULL),
ipPolyModel(NULL),
chisqrMatrix(NULL)
{
	if (NXPoints == 0 || NYPoints == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	setsize(ipxsize,ipysize);
	setsampling(ipxsampling,ipysampling);
	setNXPoints(ipxsize*ipxsampling);	
	setNYPoints(ipysize*ipysampling);
	setNTotalPoints(ipxsize*ipxsampling*ipysize*ipysampling);
	setGeometricCenter();	
	
	maxnDataPoints = 1;
	setnDataPoints(1);	
	setCubes(1);	
	
	ipPolyModel = newPolynomialMatrix(NXPoints, NYPoints);
	chisqrMatrix = newCMatrix(NXPoints, NYPoints);	
}

operaInstrumentProfile::operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling,unsigned NDataPoints) :
NXPoints(1),		// number of points in x direction
NYPoints(1),		// number of points in y directon
NTotalPoints(1),	// total number of points
Xsampling(1),		// (sampling factor for the PSF in x direction)
Ysampling(1),		// (sampling factor for the  PSF in y direction)
xsize(1),			// (size along the x direction in pix units) 
ysize(1),			// (size along the y direction in pix units)
geometricCenterX(0.5),
geometricCenterY(0.5),
nDataPoints(0),
maxnDataPoints(0),
dataCube(NULL),
errorsCube(NULL),
distd(NULL),
ipPolyModel(NULL),
chisqrMatrix(NULL)
{
	if (NXPoints == 0 || NYPoints == 0 || NDataPoints == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	setsize(ipxsize,ipysize);
	setsampling(ipxsampling,ipysampling);
	setNXPoints(ipxsize*ipxsampling);	
	setNYPoints(ipysize*ipysampling);
	setNTotalPoints(ipxsize*ipxsampling*ipysize*ipysampling);
	setGeometricCenter();	
	
	maxnDataPoints = NDataPoints;
	setnDataPoints(NDataPoints);	
	setCubes(nDataPoints);	
	
	ipPolyModel = newPolynomialMatrix(NXPoints, NYPoints);
	chisqrMatrix = newCMatrix(NXPoints, NYPoints);	
}

/*
 * Destructor
 */
operaInstrumentProfile::~operaInstrumentProfile(void) {
	
	if (dataCube) deleteCCube(dataCube);
	if (errorsCube) deleteCCube(errorsCube);
	
	dataCube = NULL;
	errorsCube = NULL;
	
	if (distd) free(distd);
	distd = NULL;
	
	if (ipPolyModel) deletePolynomialMatrix(ipPolyModel);
	if (chisqrMatrix) deleteCMatrix(chisqrMatrix);
	
	ipPolyModel = NULL;
	chisqrMatrix = NULL;	

	maxnDataPoints = nDataPoints = 0;
}

/*
 * Setter/Getters
 */

unsigned operaInstrumentProfile::getNXPoints(void){
	return NXPoints;
}

void operaInstrumentProfile::setNXPoints(unsigned npx){
	NXPoints = npx;
}

unsigned operaInstrumentProfile::getNYPoints(void){
	return NYPoints;
}

void operaInstrumentProfile::setNYPoints(unsigned npy){
	NYPoints = npy;
}

unsigned operaInstrumentProfile::getNTotalPoints(void){
	return NTotalPoints;
}

void operaInstrumentProfile::setNTotalPoints(unsigned np){
	NTotalPoints = np;
}

// (fraction of pixel to sample PSF in dispersion direction)
unsigned operaInstrumentProfile::getYsampling(void) {
	return Ysampling;
}

// (fraction of pixel to sample PSF in spatial direction)
unsigned operaInstrumentProfile::getXsampling(void) {
	return Xsampling;
}

void operaInstrumentProfile::setYsampling(unsigned ysamp) {
	Ysampling = ysamp;
}

void operaInstrumentProfile::setXsampling(unsigned xsamp) {
	Xsampling = xsamp;
}

// (fraction of pixel to sample PSF in spatial direction)
void operaInstrumentProfile::setsampling(unsigned xsamp, unsigned ysamp) {
	Xsampling = xsamp;
	Ysampling = ysamp;
}

// (size along the spatial direction in pix units) 
unsigned operaInstrumentProfile::getxsize(void) {
	return xsize;
}

// (size along the dispersion direction in pix units)
unsigned operaInstrumentProfile::getysize(void) {
	return ysize;
}

// (size along the dispersion direction in pix units)
void operaInstrumentProfile::setxsize(unsigned xs) {
	xsize = xs;
}

// (size along the dispersion direction in pix units)
void operaInstrumentProfile::setysize(unsigned ys) {
	ysize = ys;
}

// (size along the dispersion direction in pix units)
void operaInstrumentProfile::setsize(unsigned xs, unsigned ys) {
	xsize = xs;
	ysize = ys;
}

void operaInstrumentProfile::subtractOuterFrame(unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    unsigned mp = 2*NYPoints+2*NXPoints;
	if (mp == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	double *xx = (double*) malloc(mp * sizeof(double));
	double *yy = (double*) malloc(mp * sizeof(double));	
	double *z = (double*) malloc(mp * sizeof(double));	
    float *fz = (float*) malloc(mp * sizeof(float));
	if (!fz) {
		throw operaException("operaInstrumentProfile: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
	
    unsigned k=0;
    unsigned numberofNANs = 0;
    
    for (unsigned j=0; j<NYPoints; j++) {	
        for (unsigned i=0; i<NXPoints; i++) {
            if((j==0 || j==NYPoints-1) && (i==0 || i==NXPoints-1)) {
                if(!isnan(getdataCubeValues(i,j,index))) {
                    yy[k] = (double)j;
                    xx[k] = (double)i;
                    fz[k] = getdataCubeValues(i,j,index);
                    z[k] = (double)getdataCubeValues(i,j,index);
                    k++;
                } else {
                    numberofNANs++;
                }
            }
        }
    }		
    
    if(k > numberofNANs) {
        unsigned npar = 4;
        double pars[4] = {1,1,1,1};
        double chisqr;
        
        operaLMFit2DPolynomial(k,xx,yy,z,npar,pars,&chisqr);
        
        for (unsigned j=0; j<NYPoints; j++) {		
            for (unsigned i=0; i<NXPoints; i++) {	
                float value = getdataCubeValues(i,j,index) - (float)Polynomial2DFunction((double)i,(double)j,pars,npar);
                setdataCubeValues(value,i,j,index);            
            }
        }	
    } else if (k > 3) {
        float minValue = operaArrayMinValue(k,fz);
        for (unsigned j=0; j<NYPoints; j++) {		
            for (unsigned i=0; i<NXPoints; i++) {	
                float value = getdataCubeValues(i,j,index) - minValue;
                setdataCubeValues(value,i,j,index);            
            }
        }	
    }
    free(xx);
    free(yy);
    free(z);
    free(fz);	// DT one more...
}


void operaInstrumentProfile::normalizeCubeData(unsigned index) {
#ifdef RANGE_CHECK
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    float sum=0;
    for (unsigned j=0; j<NYPoints; j++) {		
        for (unsigned i=0; i<NXPoints; i++) {	
            sum += getdataCubeValues(i,j,index);
        }
    }		
    for (unsigned j=0; j<NYPoints; j++) {		
        for (unsigned i=0; i<NXPoints; i++) {	
            if(sum) {
                float value = getdataCubeValues(i,j,index)/sum;
                setdataCubeValues(value,i,j,index);
                float errorvalue = geterrorsCubeValues(i,j,index)/sum;	
                seterrorsCubeValues(errorvalue,i,j,index);
            }
        }
    }			
}

void operaInstrumentProfile::normalizeCubeData(void) {
	for(unsigned index = 0; index < getnDataPoints(); index++) {
		float sum=0;
		for (unsigned j=0; j<NYPoints; j++) {		
			for (unsigned i=0; i<NXPoints; i++) {	
				sum += getdataCubeValues(i,j,index);
			}
		}		
		for (unsigned j=0; j<NYPoints; j++) {		
			for (unsigned i=0; i<NXPoints; i++) {	
				if(sum) {
					float value = getdataCubeValues(i,j,index)/sum;
					setdataCubeValues(value,i,j,index);
					float errorvalue = geterrorsCubeValues(i,j,index)/sum;	
					seterrorsCubeValues(errorvalue,i,j,index);
				}
			}
		}			
	}
}


void operaInstrumentProfile::normalizeCubeData(PixelSet *aperturePixels, unsigned index) { 
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (NXPoints == 0 || NYPoints == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (aperturePixels->getNPixels() == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    CMatrix ipmatrix = newCMatrix(NXPoints, NYPoints);    
    CMatrix errormatrix = newCMatrix(NXPoints, NYPoints);   
    
    float sum=0;
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);
            ipmatrix[j][i] = getdataCubeValues(i,j,index);
            sum += ipmatrix[j][i]; 
            setdataCubeValues(0.0,i,j,index);
            errormatrix[j][i] = geterrorsCubeValues(i,j,index);
            seterrorsCubeValues(NAN,i,j,index);
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
    }	
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);
            if(sum) {
                float value = ipmatrix[j][i]/sum;
                setdataCubeValues(value,i,j,index);
                float errorvalue = errormatrix[j][i]/sum;	
                seterrorsCubeValues(errorvalue,i,j,index);
            }
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                        
    }
    deleteCMatrix(ipmatrix);
    deleteCMatrix(errormatrix); 
}

void operaInstrumentProfile::normalizeCubeData(PixelSet *aperturePixels) {
	if (getnDataPoints() == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (NXPoints == 0 || NYPoints == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (aperturePixels->getNPixels() == 0) {
		throw operaException("operaInstrumentProfile: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    CMatrix ipmatrix = newCMatrix(NXPoints, NYPoints);    
    CMatrix errormatrix = newCMatrix(NXPoints, NYPoints); 
    
	for(unsigned index = 0; index < getnDataPoints(); index++) {
        float sum=0;
        for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
            if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
               aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
                
                unsigned i = (unsigned)aperturePixels->getiIndex(pix);
                unsigned j = (unsigned)aperturePixels->getjIndex(pix);
                ipmatrix[j][i] = getdataCubeValues(i,j,index);
                sum += ipmatrix[j][i];
                setdataCubeValues(0.0,i,j,index);
                errormatrix[j][i] = geterrorsCubeValues(i,j,index);
                seterrorsCubeValues(NAN,i,j,index);
            } else {
#ifdef PRINT_OUTOFBOUNDS
                cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
                cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
            }                        
        }	
        for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
            if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
               aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
                
                unsigned i = (unsigned)aperturePixels->getiIndex(pix);
                unsigned j = (unsigned)aperturePixels->getjIndex(pix);
                if(sum) {
                    float value = ipmatrix[j][i]/sum;
                    setdataCubeValues(value,i,j,index);
                    float errorvalue = errormatrix[j][i]/sum;	
                    seterrorsCubeValues(errorvalue,i,j,index);
                }
            } else {
#ifdef PRINT_OUTOFBOUNDS
                cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
                cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
            }                        
        }		
	}
    deleteCMatrix(ipmatrix);
    deleteCMatrix(errormatrix);     
}



// set x,y coordinates of geometric center (requires x,y size and x,y sampling)
void operaInstrumentProfile::setGeometricCenter(void){
	geometricCenterX = (float)xsize/2.0;
	geometricCenterY = (float)ysize/2.0;	
}

float operaInstrumentProfile::getGeometricCenterX(void){
	return geometricCenterX;
}

float operaInstrumentProfile::getGeometricCenterY(void){	
	return geometricCenterY;
}

float operaInstrumentProfile::getIPixXCoordinate(unsigned i){
    float halfXgridstep = 1.0/float(2*getXsampling());
    return (float(1+2*i)*halfXgridstep - geometricCenterX);
	//return (((float)i+0.5)/(float)Xsampling - geometricCenterX); 
}

float operaInstrumentProfile::getIPixYCoordinate(unsigned j){
    float halfYgridstep = 1.0/float(2*getYsampling());
    return (float(1+2*j)*halfYgridstep - geometricCenterY);    
	//return (((float)j+0.5)/(float)Ysampling - geometricCenterY);	
}

unsigned operaInstrumentProfile::getIPixiIndex(float xcoord){
	int i = (int)round((xcoord + geometricCenterX)*(float)Xsampling - 0.5);
	
	if(i < 0 || i >= (int)NXPoints) {
		throw operaException("operaInstrumentProfile: out of bound x-coordinate.",operaErrorInstrumentProfileImproperCoordinateRequest, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	return (unsigned)i;
}

unsigned operaInstrumentProfile::getIPixjIndex(float ycoord){
	int j = (int)round((ycoord + geometricCenterY)*(float)Ysampling - 0.5);
	
	if(j < 0 || j >= (int)NYPoints) {
		throw operaException("operaInstrumentProfile: out of bound y-coordinate.",operaErrorInstrumentProfileImproperCoordinateRequest, __FILE__, __FUNCTION__, __LINE__);	
	}
	
	return (unsigned)j;
}

unsigned operaInstrumentProfile::getnDataPoints(void) {
	return nDataPoints;
}

void operaInstrumentProfile::setnDataPoints(unsigned NDataPoints) {
	if (NDataPoints > maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	nDataPoints = NDataPoints;
}

CMatrix* operaInstrumentProfile::getdataCube(void) {
	return dataCube;
}

CMatrix operaInstrumentProfile::getdataCube(unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return dataCube[index];	
}

CMatrix operaInstrumentProfile::getdataCubeValues(unsigned index, CMatrix dataSlice) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			 dataSlice[j][i] = dataCube[index][j][i];
		}
	}	
    return dataSlice;    
}

CMatrix operaInstrumentProfile::getdataCubeValues(PixelSet *aperturePixels, unsigned index, CMatrix dataSlice) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);
            dataSlice[j][i] = dataCube[index][j][i];
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
	}	
    return dataSlice;
}

void operaInstrumentProfile::setdataCubeValues(operaFITSImage &image, operaFITSImage &badpix, float xcenter, float ycenter, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for (unsigned j=0; j<NYPoints; j++) {	
        float YCenterOfSubPixel = ycenter + getIPixYCoordinate(j);				
        unsigned yy = (unsigned)floor(YCenterOfSubPixel);             
        for (unsigned i=0; i<NXPoints; i++) {
            float XCenterOfSubPixel = xcenter + getIPixXCoordinate(i);				
            unsigned xx = (unsigned)floor(XCenterOfSubPixel); 
            if (xx > 0 && xx < image.getnaxis1() && yy > 0 && yy < image.getnaxis2()){
                if (image[yy][xx] < SATURATIONLIMIT && badpix[yy][xx] == 1 && (float)image[yy][xx] > 0) {                    
                    dataCube[index][j][i] = (float)image[yy][xx];                       
                } else {
                    dataCube[index][j][i] = NAN;
                }
            }
        }
    }    
}

void operaInstrumentProfile::setdataCubeValues(operaFITSImage &image, operaFITSImage &badpix, PixelSet *aperturePixels, float xcenter, float ycenter, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);
            
            float YCenterOfSubPixel = ycenter + getIPixYCoordinate(j);				
            unsigned yy = (unsigned)floor(YCenterOfSubPixel);             
            
            float XCenterOfSubPixel = xcenter + getIPixXCoordinate(i);				
            unsigned xx = (unsigned)floor(XCenterOfSubPixel); 
            
            if (xx > 0 && xx < image.getnaxis1() && yy > 0 && yy < image.getnaxis2()) {
                if (image[yy][xx] < SATURATIONLIMIT && badpix[yy][xx] == 1 && (float)image[yy][xx] > 0) {                    
                    dataCube[index][j][i] = (float)image[yy][xx];                       
                } else {
                    dataCube[index][j][i] = NAN;
                }
            }
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
    }    
}

void operaInstrumentProfile::setdataCubeValues(CMatrix DataMatrix, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			dataCube[index][j][i] = (float)(DataMatrix[j][i]);
		}
	}	
}

void operaInstrumentProfile::setdataCubeValues(PixelSet *aperturePixels, CMatrix DataMatrix, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);
            dataCube[index][j][i] = (float)(DataMatrix[j][i]);
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
	}	
}

void operaInstrumentProfile::setdataCubeValues(float DataValue, unsigned i, unsigned j, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	dataCube[index][j][i] = DataValue;
}

float operaInstrumentProfile::getdataCubeValues(unsigned i, unsigned j, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return dataCube[index][j][i];
}

float operaInstrumentProfile::getdataGivenCoords(float xcoord, float ycoord, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return dataCube[index][getIPixjIndex(ycoord)][getIPixiIndex(xcoord)];
}

CMatrix* operaInstrumentProfile::geterrorsCube(void) {
	return errorsCube;
}

CMatrix operaInstrumentProfile::geterrorsCube(unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return errorsCube[index];
}

void operaInstrumentProfile::seterrorsCube(CMatrix ErrorsMatrix, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	errorsCube[index] = ErrorsMatrix;	
}

void operaInstrumentProfile::seterrorsCubeValues(CMatrix ErrorsMatrix, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			errorsCube[index][j][i] = (float)(ErrorsMatrix[j][i]);
		}
	}	
}

void operaInstrumentProfile::seterrorsCubeValues(float ErrorValue, unsigned i, unsigned j, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	errorsCube[index][j][i] = ErrorValue;
}

float operaInstrumentProfile::geterrorsCubeValues(unsigned i, unsigned j, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return errorsCube[index][j][i];
}

float* operaInstrumentProfile::getdistdVector(void){
	return distd;
}

float operaInstrumentProfile::getdistd(unsigned index){
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return distd[index];
}

void operaInstrumentProfile::setdistdVector(float *Distd){
	if (distd) {
		free(distd);
	}
    distd = (float *)malloc(sizeof(float)*nDataPoints);
	if (!distd) {
		throw operaException("operaInstrumentProfile: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned i=0;i<nDataPoints;i++) {
        distd[i] = Distd[i];
    }
}

void operaInstrumentProfile::setdistd(float DistdValue, unsigned index){
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	distd[index] = DistdValue;
}

void operaInstrumentProfile::setCubes(unsigned NDataPoints){
	
	setnDataPoints(NDataPoints);
    
	if (dataCube) {
		deleteCCube(dataCube);
	}
	dataCube = newCCube(NXPoints, NYPoints, NDataPoints);
	if (errorsCube) {
		deleteCCube(errorsCube);
	}
	errorsCube = newCCube(NXPoints, NYPoints, NDataPoints);
	if (distd) {
		free(distd);
	}
	distd = (float *)malloc(sizeof(float)*NDataPoints);
	if (!distd) {
		throw operaException("operaInstrumentProfile: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
}

void operaInstrumentProfile::printData(unsigned index, ostream *pout) {	
	if (pout) {
		*pout << "#i\tj\tk\tx\ty\td\tz\tzerr\n" << endl;
		for (unsigned j=0; j<NYPoints; j++) {
			for (unsigned i=0; i<NXPoints; i++) {
				*pout << i << "\t" << j << "\t" << index << "\t";
				*pout << getIPixXCoordinate(i) << 
				"\t" << getIPixYCoordinate(j) << 
				"\t" << (float)(distd[index]) << 			
				"\t" << getdataCubeValues(i,j,index) <<
				"\t" << geterrorsCubeValues(i,j,index) << endl;
			}
			*pout << endl;
		}	
		
	}
}

void operaInstrumentProfile::printData(unsigned index, int ordernumber, ostream *pout) {	

	if (pout) {
		*pout << "#o\ti\tj\tk\tx\ty\td\tz\tzerr\n" << endl;
		for (unsigned j=0; j<NYPoints; j++) {
			for (unsigned i=0; i<NXPoints; i++) {
				*pout << ordernumber << "\t"
                      << i << "\t" 
                      << j << "\t" 
                      << index << "\t";
				*pout << getIPixXCoordinate(i) << 
				"\t" << getIPixYCoordinate(j) << 
				"\t" << (float)(distd[index]) << 			
				"\t" << getdataCubeValues(i,j,index) <<
				"\t" << geterrorsCubeValues(i,j,index) << endl;
			}
			*pout << endl;
		}	
		
	}
}

void operaInstrumentProfile::printModel(float DistanceInPixels, unsigned ordernumber, ostream *pout) {
    
	if (pout) {
		*pout << "#o\ti\tj\tx\ty\td\tz\tzerr\n" << endl;
		for (unsigned j=0; j<NYPoints; j++) {
			for (unsigned i=0; i<NXPoints; i++) {
				*pout << ordernumber << "\t"
                << i << "\t"
                << j << "\t";
				*pout << getIPixXCoordinate(i) <<
				"\t" << getIPixYCoordinate(j) <<
				"\t" << DistanceInPixels <<
                "\t" << getipDataFromPolyModel(DistanceInPixels,i,j) <<
				"\t" << sqrt(getchisqrMatrixValue(i,j)) << endl;
			}
			*pout << endl;
		}
	}
}

/*
 * DT/EM Fixed memory leak by copying. Oct 24 2012
 */
void operaInstrumentProfile::setipPolyModel(PolynomialMatrix IPPolyModel) {
	unsigned cols = getPolynomialMatrixCols(IPPolyModel);
	unsigned rows = getPolynomialMatrixRows(IPPolyModel);
	
	for (unsigned j=0; j<rows; j++) {
		for (unsigned i=0; i<cols; i++) {
			PolynomialCoeffs_t *pp = ipPolyModel[j][i];
			PolynomialCoeffs_t *PP = IPPolyModel[j][i];
			pp->orderofPolynomial = PP->orderofPolynomial;
			unsigned coeffs = PP->orderofPolynomial;
			for (unsigned coeff=0; coeff<coeffs; coeff++) {
				pp->p[coeff] = PP->p[coeff];
				pp->e[coeff] = PP->e[coeff];
			}
			pp->polychisqr = PP->polychisqr;
		}
	}
}

PolynomialMatrix operaInstrumentProfile::getipPolyModel(void) {
	return ipPolyModel;
}

PolynomialCoeffs_t *operaInstrumentProfile::getipPolyModelCoefficients(unsigned i,unsigned j) {
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return ipPolyModel[j][i];
}

// DT -- BAD INTERFACE -- copies a pointer!!!! Mar 2013
void operaInstrumentProfile::setipPolyModelCoefficients(PolynomialCoeffs_t *PolyModelCoeffs,unsigned i,unsigned j) {
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	ipPolyModel[j][i] = PolyModelCoeffs;
}

CMatrix operaInstrumentProfile::getipDataFromPolyModel(float d, CMatrix ipmatrix) {
	
	Polynomial polyn;  
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
			polyn.setPolynomialCoeffs(pp);
			ipmatrix[j][i] = (float)polyn.Evaluate((double)d);			
		}
	}		
	return ipmatrix;
}

CMatrix operaInstrumentProfile::getipDataFromPolyModel(PixelSet *aperturePixels, float d, CMatrix ipmatrix) {
	
	Polynomial polyn;  
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);    
            PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
            polyn.setPolynomialCoeffs(pp);
            ipmatrix[j][i] = (float)polyn.Evaluate((double)d);	
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
	}		
	return ipmatrix;
}

/*
CMatrix operaInstrumentProfile::getipDataFromPolyModel(float d) {
	
	CMatrix ipmatrix = newCMatrix(NXPoints, NYPoints);
	
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
			Polynomial polyn(pp);  
			ipmatrix[j][i] = polyn.Evaluate((double)d);			
		}
	}		
	
	return ipmatrix;
}
*/

float operaInstrumentProfile::getipDataFromPolyModel(float d, unsigned i, unsigned j) {
	
	PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
	
	//Commented out the original way of doing this
	/*Polynomial polyn(pp);  
	float modelValue = (float)polyn.Evaluate((double)d);	
	return modelValue;*/
	
	//Switched to inline polynomial evaluation to optimize for operaExtraction
	int polyOrder = pp->orderofPolynomial;
	float* polyCoeffs = pp->p;
	
	double eval = polyCoeffs[0];
	double x_to_i = d;
	for (unsigned i=1; i<polyOrder; i++) {
		eval += polyCoeffs[i]*x_to_i;
		x_to_i *= d;
	}
	return (float)eval;
}

float operaInstrumentProfile::getipDataFromPolyModel(unsigned i, unsigned j, unsigned index) {
	
	PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
	Polynomial polyn(pp);  
	float modelValue = (float)polyn.Evaluate((double)getdistd(index));	
	
	return modelValue;
}

CMatrix operaInstrumentProfile::getipDataFromPolyModel(unsigned index, CMatrix ipmatrix) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
			Polynomial polyn(pp);  
			ipmatrix[j][i] = (float)polyn.Evaluate((double)getdistd(index));			
		}
	}		
	
	return ipmatrix;
}


CMatrix operaInstrumentProfile::getipDataFromPolyModel(PixelSet *aperturePixels, unsigned index, CMatrix ipmatrix) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix); 
            
            PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
            Polynomial polyn(pp);  
            ipmatrix[j][i] = (float)polyn.Evaluate((double)getdistd(index));
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        } 
	}		
	
	return ipmatrix;
}

void operaInstrumentProfile::setipDataFromPolyModel(unsigned index){
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
			Polynomial polyn(pp);  
			setdataCubeValues(polyn.Evaluate((double)getdistd(index)),i,j,index);
		}
	}		
}

void operaInstrumentProfile::setipDataFromPolyModel(float d, unsigned index){
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
			Polynomial polyn(pp);  
			setdataCubeValues(polyn.Evaluate((double)d),i,j,index);
		}
	}		
}

void operaInstrumentProfile::setdataCubeFromPolyModel(void) {
	for(unsigned index = 0; index < getnDataPoints(); index++) {
		for (unsigned j=0; j<NYPoints; j++) {		
			for (unsigned i=0; i<NXPoints; i++) {			
				setdataCubeValues(getipDataFromPolyModel(distd[index],i,j),i,j,index);
			}
		}			
	}	
}

void operaInstrumentProfile::setdataCubeFromPolyModel(unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for (unsigned j=0; j<NYPoints; j++) {		
        for (unsigned i=0; i<NXPoints; i++) {	
            setdataCubeValues(getipDataFromPolyModel(distd[index],i,j),i,j,index);
        }
    }			
}

void operaInstrumentProfile::setipDataFromPolyModel(PixelSet *aperturePixels, unsigned index){
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix); 
            PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
            Polynomial polyn(pp);  
            setdataCubeValues(polyn.Evaluate((double)getdistd(index)),i,j,index);
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
	}		
}

void operaInstrumentProfile::setipDataFromPolyModel(PixelSet *aperturePixels, float d, unsigned index){
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix); 
            PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
            Polynomial polyn(pp);  
            setdataCubeValues(polyn.Evaluate((double)d),i,j,index);
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        } 
	}		
}

void operaInstrumentProfile::setdataCubeFromPolyModel(PixelSet *aperturePixels) {
	for(unsigned index = 0; index < getnDataPoints(); index++) {
        for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
            if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
               aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
                
                unsigned i = (unsigned)aperturePixels->getiIndex(pix);
                unsigned j = (unsigned)aperturePixels->getjIndex(pix);		
                setdataCubeValues(getipDataFromPolyModel(distd[index],i,j),i,j,index);
            } else {
#ifdef PRINT_OUTOFBOUNDS
                cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
                cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
            }                        
		}			
	}	
}

void operaInstrumentProfile::setdataCubeFromPolyModel(PixelSet *aperturePixels, unsigned index) {
	if (index >= maxnDataPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix); 	
            setdataCubeValues(getipDataFromPolyModel(distd[index],i,j),i,j,index);
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        } 
    }			
}

void operaInstrumentProfile::FitPolyMatrixtoIPDataVector(unsigned coeffs, bool witherrors) {
    
    if(ipPolyModel)
        deletePolynomialMatrix(ipPolyModel);
    if(chisqrMatrix)
        deleteCMatrix(chisqrMatrix);
    
	ipPolyModel = NULL;
	chisqrMatrix = NULL;
	
	ipPolyModel = newPolynomialMatrix(NXPoints, NYPoints);
	chisqrMatrix = newCMatrix(NXPoints, NYPoints);	
	
	double *par = (double *)malloc(coeffs*sizeof(double));
	double *parError = (double *)malloc(coeffs*sizeof(double)); 
	double *ytmp = (double *)malloc(nDataPoints*sizeof(double));
	double *ytmpError = (double *)malloc(nDataPoints*sizeof(double));        
	double *xtmp = (double *)malloc(nDataPoints*sizeof(double));
	if (!xtmp) {
		throw operaException("operaInstrumentProfile: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {
			
			int nparbestfit = coeffs;
			double bestchisqr = BIG;
			double chisqr;
			
            unsigned npts = 0;
			for(unsigned index=0;index<nDataPoints;index++){	
                if(!isnan(getdataCubeValues(i,j,index))) {
                    ytmp[npts] = (double)getdataCubeValues(i,j,index);				
                    xtmp[npts] = (double)getdistd(index);
                    ytmpError[npts] = (double)geterrorsCubeValues(i,j,index);  
                    npts++;
                }
			}
			
            if(npts < coeffs) {
                coeffs = npts;
            } 
            
            if(coeffs) {
                for (unsigned currentfit=1; currentfit<=coeffs; currentfit++) {
                    for	(unsigned k=0; k<currentfit; k++) {
                        par[k] = 1.0;
                    }
                    if (witherrors) {
                        operaMPFitPolynomial(npts, xtmp, ytmp, ytmpError, currentfit, par, parError, &chisqr);
                        if (fabs(chisqr-1.0) < bestchisqr) {
                            bestchisqr = chisqr;
                            nparbestfit = currentfit;
                        }
                        
                    } else {
                        operaLMFitPolynomial(npts, xtmp, ytmp, currentfit, par, &chisqr);
                        if (chisqr < bestchisqr) {                
                            bestchisqr = chisqr;
                            nparbestfit = currentfit;
                        }
                    }
                }
                for	(int k=0; k<nparbestfit; k++) {
                    par[k] = 1.0;
                }
                
                if (witherrors) {
                    operaMPFitPolynomial(npts, xtmp, ytmp, ytmpError, nparbestfit, par, parError, &chisqr);                        
                } else {
                    operaLMFitPolynomial(npts, xtmp, ytmp, nparbestfit, par, &chisqr);
                }
                
                PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
                
                pp->orderofPolynomial = nparbestfit;
                for	(int k=0; k<nparbestfit; k++) {
                    pp->p[k] = par[k];
                }
                
                setchisqrMatrixValue((float)chisqr,i,j);
                
#ifdef PRINT_DEBUG
                cout << "i="<< i << " j=" << j << " np=" << nparbestfit;
                for	(unsigned k=0; k<nparbestfit; k++) {
                    cout << " p["<<k<<"]=" << par[k];
                }
                cout << " chi2=" << chisqr << endl;	
                
                if(i==(0*NXPoints/8) || i==(1*NXPoints/8) || i==(2*NXPoints/8) || i==(3*NXPoints/8) || i==(4*NXPoints/8) || 
                   i==(5*NXPoints/8) || i==(6*NXPoints/8) || i==(7*NXPoints/8) || i==(8*NXPoints/8)){
                    for(unsigned jj=1;jj<npts;jj++){
                        cout << i << " " << xtmp[jj] << " " << ytmp[jj] << " " << ytmpError[jj] << " " << getipDataFromPolyModel((float)xtmp[jj],i,j) << " "  << fabs(xtmp[jj] - xtmp[jj-1]) << endl;
                    }
                }  
#endif      
            }
        }
	}		
	free(par);
	free(parError);    
	free(ytmp);
    free(ytmpError);    
	free(xtmp);			
}


void operaInstrumentProfile::FitPolyMatrixtoIPDataVector(PixelSet *aperturePixels, unsigned coeffs, bool witherrors) {
    
    if(ipPolyModel)
        deletePolynomialMatrix(ipPolyModel);
    if(chisqrMatrix)
        deleteCMatrix(chisqrMatrix);
    
	ipPolyModel = NULL;
	chisqrMatrix = NULL;
	
	ipPolyModel = newPolynomialMatrix(NXPoints, NYPoints);
	chisqrMatrix = newCMatrix(NXPoints, NYPoints);	
	
	double *par = (double *)malloc(coeffs*sizeof(double));
	double *parError = (double *)malloc(coeffs*sizeof(double)); 
	double *ytmp = (double *)malloc(nDataPoints*sizeof(double));
	double *ytmpError = (double *)malloc(nDataPoints*sizeof(double));        
	double *xtmp = (double *)malloc(nDataPoints*sizeof(double));
	if (!xtmp) {
		throw operaException("operaInstrumentProfile: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);
            
            int nparbestfit = coeffs;
            double bestchisqr = BIG;
            double chisqr;
            
            unsigned npts = 0;
            for(unsigned index=0;index<nDataPoints;index++){	
                if(!isnan(getdataCubeValues(i,j,index))) {
                    ytmp[npts] = (double)getdataCubeValues(i,j,index);				
                    xtmp[npts] = (double)getdistd(index);
                    ytmpError[npts] = (double)geterrorsCubeValues(i,j,index);  
                    npts++;
                }
            }
            
            if(npts < coeffs) {
                coeffs = npts;
            } 
            
            if(coeffs) {
                for (int currentfit=1; currentfit<=(int)coeffs; currentfit++) {
                    for	(int k=0; k<(int)currentfit; k++) {
                        par[k] = 1.0;
                    }
                    if (witherrors) {
                        operaMPFitPolynomial(npts, xtmp, ytmp, ytmpError, currentfit, par, parError, &chisqr);
                        if (fabs(chisqr-1.0) < bestchisqr) {
                            bestchisqr = chisqr;
                            nparbestfit = currentfit;
                        }
                        
                    } else {
                        operaLMFitPolynomial(npts, xtmp, ytmp, currentfit, par, &chisqr);
                        if (chisqr < bestchisqr) {                
                            bestchisqr = chisqr;
                            nparbestfit = currentfit;
                        }
                    }
                }
                for	(int k=0; k<nparbestfit; k++) {
                    par[k] = 1.0;
                }
                
                if (witherrors) {
                    operaMPFitPolynomial(npts, xtmp, ytmp, ytmpError, nparbestfit, par, parError, &chisqr);                        
                } else {
                    operaLMFitPolynomial(npts, xtmp, ytmp, nparbestfit, par, &chisqr);
                }
                
                PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
                
                pp->orderofPolynomial = nparbestfit;
                for	(int k=0; k<nparbestfit; k++) {
                    pp->p[k] = par[k];
                }
                
                setchisqrMatrixValue((float)chisqr,i,j);
                
#ifdef PRINT_DEBUG
                cout << "i="<< i << " j=" << j << " np=" << nparbestfit;
                for	(unsigned k=0; k<nparbestfit; k++) {
                    cout << " p["<<k<<"]=" << par[k];
                }
                cout << " chi2=" << chisqr << endl;	
                
                if(i==(0*NXPoints/8) || i==(1*NXPoints/8) || i==(2*NXPoints/8) || i==(3*NXPoints/8) || i==(4*NXPoints/8) || 
                   i==(5*NXPoints/8) || i==(6*NXPoints/8) || i==(7*NXPoints/8) || i==(8*NXPoints/8)){
                    for(unsigned jj=1;jj<npts;jj++){
                        cout << i << " " << xtmp[jj] << " " << ytmp[jj] << " " << ytmpError[jj] << " " << getipDataFromPolyModel((float)xtmp[jj],i,j) << " "  << fabs(xtmp[jj] - xtmp[jj-1]) << endl;
                    }
                }  
#endif      
            }
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
	}		
	free(par);
	free(parError);    
	free(ytmp);
    free(ytmpError);    
	free(xtmp);			
}

void operaInstrumentProfile::FitMediantoIPDataVector(void) {
    unsigned coeffs = 1;
    double *par = (double *)malloc(coeffs*sizeof(double));
    
    if(ipPolyModel)
        deletePolynomialMatrix(ipPolyModel);
    if(chisqrMatrix)
        deleteCMatrix(chisqrMatrix);
    
	ipPolyModel = NULL;
	chisqrMatrix = NULL;
	
	ipPolyModel = newPolynomialMatrix(NXPoints, NYPoints);
	chisqrMatrix = newCMatrix(NXPoints, NYPoints);	
	
	float *ytmp = (float *)malloc(nDataPoints*sizeof(float));
	if (!ytmp) {
		throw operaException("operaInstrumentProfile: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    
	for (unsigned j=0; j<NYPoints; j++) {
		for (unsigned i=0; i<NXPoints; i++) {

			for(unsigned index=0;index<nDataPoints;index++) {	
				ytmp[index] = getdataCubeValues(i,j,index);
			}

            par[0] = (double)operaArrayMedian(nDataPoints,ytmp);
            
            PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
            
            pp->orderofPolynomial = coeffs;
            for	(unsigned k=0; k<coeffs; k++) {
                pp->p[k] = par[k];
            }
            
            float chisqr = operaArrayChisqr(nDataPoints, ytmp, (float)par[0], nDataPoints - 1);
            setchisqrMatrixValue(chisqr,i,j);
        }
	}		
	free(ytmp);
    free(par);
}

void operaInstrumentProfile::FitMediantoIPDataVector(PixelSet *aperturePixels) {
    unsigned coeffs = 1;
    double *par = (double *)malloc(coeffs*sizeof(double));

    if(ipPolyModel)
        deletePolynomialMatrix(ipPolyModel);
    if(chisqrMatrix)
        deleteCMatrix(chisqrMatrix);
    
	ipPolyModel = NULL;
	chisqrMatrix = NULL;
	
	ipPolyModel = newPolynomialMatrix(NXPoints, NYPoints);
	chisqrMatrix = newCMatrix(NXPoints, NYPoints);	
	
	float *ytmp = (float *)malloc(nDataPoints*sizeof(float));
	if (!ytmp) {
		throw operaException("operaInstrumentProfile: ", operaErrorNoMemory, __FILE__, __FUNCTION__, __LINE__);	
	}
    
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {    
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);
            
            for(unsigned index=0;index<nDataPoints;index++) {	
                ytmp[index] = getdataCubeValues(i,j,index);				
            }
            
            par[0] = (double)operaArrayMedian(nDataPoints,ytmp);
            
            PolynomialCoeffs_t *pp = getipPolyModelCoefficients(i,j);
            
            pp->orderofPolynomial = coeffs;
            for	(unsigned k=0; k<coeffs; k++) {
                pp->p[k] = par[k];
            }
            
            float chisqr = operaArrayChisqr(nDataPoints, ytmp, (float)par[0], nDataPoints - 1);
            setchisqrMatrixValue(chisqr,i,j);
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }
	}		
	free(ytmp);
    free(par);
}

CMatrix operaInstrumentProfile::getchisqrMatrix(void){
	return chisqrMatrix;
}

void operaInstrumentProfile::setchisqrMatrix(CMatrix ChisqrMatrix){
	chisqrMatrix = ChisqrMatrix;
}

float operaInstrumentProfile::getchisqrMatrixValue(unsigned i, unsigned j){
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return chisqrMatrix[j][i];
}

void operaInstrumentProfile::setchisqrMatrixValue(float Value, unsigned i, unsigned j){
	if (j >= NYPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (i >= NXPoints) {
		throw operaException("operaInstrumentProfile: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	chisqrMatrix[j][i] = Value;
}

void operaInstrumentProfile::deleteDataCubes(void) {
	
    if (dataCube) {
		deleteCCube(dataCube);
	}
	dataCube = NULL;    
	if (errorsCube) {
		deleteCCube(errorsCube);
	}
	errorsCube = NULL;    
	if (distd) {
		free(distd);
	}
    distd = NULL;
	
    nDataPoints = 0;    
}

float operaInstrumentProfile::getIPphotoCenterX(float d) {
    
    float sumflux = 0;
    float xc = 0;
    
    for (unsigned j=0; j<NYPoints; j++) {		
        for (unsigned i=0; i<NXPoints; i++) {				
            sumflux += getipDataFromPolyModel(d,i,j);
            xc += getipDataFromPolyModel(d,i,j)*getIPixXCoordinate(i);
        }
    }	    
    
    return xc/sumflux;
}

float operaInstrumentProfile::getIPphotoCenterY(float d) {
    
    float sumflux = 0;
    float yc = 0;
    
    for (unsigned i=0; i<NXPoints; i++) {    
        for (unsigned j=0; j<NYPoints; j++) {
            sumflux += getipDataFromPolyModel(d,i,j);
            yc += getipDataFromPolyModel(d,i,j)*getIPixYCoordinate(j);
        }
    }	    
    
    return yc/sumflux;    
}

float operaInstrumentProfile::getIPphotoCenterX(PixelSet *aperturePixels, float d) {
    
    float sumflux = 0;
    float xc = 0;
    
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);				
            sumflux += getipDataFromPolyModel(d,i,j);
            xc += getipDataFromPolyModel(d,i,j)*getIPixXCoordinate(i);    
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }
    }	    
    
    return xc/sumflux;
}

float operaInstrumentProfile::getIPphotoCenterY(PixelSet *aperturePixels, float d) {
    
    float sumflux = 0;
    float yc = 0;
    
    for(unsigned pix=0; pix<aperturePixels->getNPixels(); pix++) {        
        if(aperturePixels->getiIndex(pix) >=0 && (unsigned)aperturePixels->getiIndex(pix) < NXPoints &&  
           aperturePixels->getjIndex(pix) >=0 && (unsigned)aperturePixels->getjIndex(pix) < NYPoints) { 
            
            unsigned i = (unsigned)aperturePixels->getiIndex(pix);
            unsigned j = (unsigned)aperturePixels->getjIndex(pix);				
            sumflux += getipDataFromPolyModel(d,i,j);
            yc += getipDataFromPolyModel(d,i,j)*getIPixYCoordinate(j);   
        } else {
#ifdef PRINT_OUTOFBOUNDS
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getiIndex(pix) (" << aperturePixels->getiIndex(pix) << ") >= 0 || < " << NXPoints << endl;
            cerr << "operaInstrumentProfile:: Warning: aperturePixels->getjIndex(pix) (" << aperturePixels->getjIndex(pix) << ") >= 0 || < " << NYPoints << endl;
#endif            
        }                    
    }	    
    
    return yc/sumflux;    
}

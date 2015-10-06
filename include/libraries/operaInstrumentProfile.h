#ifndef OPERAINSTRUMENTPROFILE_H
#define OPERAINSTRUMENTPROFILE_H

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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include <ostream>

#include "libraries/operaLibCommon.h"			// for CMatrix

#include "libraries/operaFITSImage.h"
#include "libraries/Polynomial.h"
#include "libraries/PixelSet.h"

class operaSpectralOrder;

/*! 
 * \brief class encapsulating the operaInstrumentProfile.
 * \details The instrument profile (IP) consists of a data set (pixelized image) 
 * \details representing the distribution of the fraction of flux for a uniformly illuminated 
 * \details monochromatic image of the entrance slit projected on the spectrograph focal plane. 
 * \details The fiducial set of coordinates chosen is such that the ordinate is oriented along 
 * \details with the dispersion direction and the abscissa with the spatial direction.
 * \return none
 * \file operaInstrumentProfile.h
 * \ingroup libraries
 */
class operaInstrumentProfile {
	
private:
	unsigned NXPoints;
	unsigned NYPoints;
	unsigned NTotalPoints;
	unsigned Xsampling;		// (sampling factor for the PSF in x direction)
	unsigned Ysampling;		// (sampling factor for the  PSF in y direction)
	unsigned xsize;			// (size along the x direction in pix units) 
	unsigned ysize;			// (size along the y direction in pix units)
	float geometricCenterX; // geometric x center in pixel units 
	float geometricCenterY; // geometric y center in pixel units 
	
	unsigned nDataPoints;
	unsigned maxnDataPoints;
	CCube dataCube;
	CCube errorsCube;
	
	float *distd;
	
	PolynomialMatrix ipPolyModel;
    CMatrix chisqrMatrix;	
	
public:
	
	/*
	 * Constructor
	 */
	operaInstrumentProfile(void);
	
	operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling);
	
	operaInstrumentProfile(unsigned ipxsize,unsigned ipxsampling,unsigned ipysize,unsigned ipysampling,unsigned NDataPoints);
	
	/*
	 * Destructor
	 */
	~operaInstrumentProfile(void);
	
	/*
	 * Methods
	 */
	
	void setGlobalInstrumentProfileDimensions(unsigned IPxsize, unsigned IPxsampling, unsigned IPysize, unsigned IPysampling);
	
	void setGlobalInstrumentProfileDimensions(unsigned nDataPoints, unsigned IPxsize, unsigned IPxsampling, unsigned IPysize, unsigned IPysampling);
	
	unsigned getNXPoints(void) const;
	
	void setNXPoints(unsigned npx);
	
	unsigned getNYPoints(void) const;
	
	void setNYPoints(unsigned npy);
	
	unsigned getNTotalPoints(void) const;
	
	void setNTotalPoints(unsigned np);
	
	unsigned getYsampling(void) const;		// (fraction of pixel to sample PSF in dispersion direction)
	
	unsigned getXsampling(void) const;		// (fraction of pixel to sample PSF in spatial direction)
	
	void setYsampling(unsigned ysamp);		// (fraction of pixel to sample PSF in dispersion direction)
	
	void setXsampling(unsigned xsamp);		// (fraction of pixel to sample PSF in spatial direction)
	
	void setsampling(unsigned xsamp, unsigned ysamp);		// (fraction of pixel to sample PSF in spatial direction)
	
	unsigned getxsize(void) const;			// (size along the spatial direction in pix units) 
	
	unsigned getysize(void) const;			// (size along the dispersion direction in pix units)
	
	void setxsize(unsigned xs);   // (size along the spatial direction in pix units) 
	
	void setysize(unsigned ys);			// (size along the dispersion direction in pix units)
	
	void setsize(unsigned xs, unsigned ys);			// (size along the dispersion direction in pix units)
	
	void normalizeCubeData(void); 			// normalize data to Sum = 1
    
	void normalizeCubeData(unsigned index);   
    
    void normalizeCubeData(PixelSet *aperturePixels, unsigned index);
    
    void normalizeCubeData(PixelSet *aperturePixels);
    
	void subtractOuterFrame(unsigned index);    
    
	void setGeometricCenter(void); // set x,y coordinates of the geometric center (requires x,y size and x,y sampling)
	
	float getGeometricCenterX(void) const;
	
	float getGeometricCenterY(void) const;
	
	float getIPixXCoordinate(unsigned i) const;
	
	float getIPixYCoordinate(unsigned j) const;
	
	unsigned getIPixiIndex(float xcoord) const;
	
	unsigned getIPixjIndex(float ycoord) const;
	
	void printData(unsigned index, ostream *pout);
	
	void printData(unsigned index, int ordernumber, ostream *pout);
	
    void printModel(float DistanceInPixels, unsigned ordernumber, ostream *pout);
    
	void setnDataPoints(unsigned NDataPoints);
	
	unsigned getnDataPoints(void) const;
	
	CCube getdataCube(void);
	
	CMatrix getdataCube(unsigned index);
	
    void setdataCubeValues(operaFITSImage &image, operaFITSImage &badpix, float xcenter, float ycenter, unsigned index);
    
    void setdataCubeValues(operaFITSImage &image, operaFITSImage &badpix, PixelSet *aperturePixels, float xcenter, float ycenter, unsigned index);
    
	void setdataCubeValues(CMatrix DataMatrix, unsigned index);
    
    void setdataCubeValues(PixelSet *aperturePixels, CMatrix DataMatrix, unsigned index);
	
	void setdataCubeValues(float DataValue, unsigned i, unsigned j, unsigned index);
	
	float getdataCubeValues(unsigned i, unsigned j, unsigned index) const;
	
	float getdataGivenCoords(float xcoord, float ycoord, unsigned index) const;

	CCube geterrorsCube(void);
	
	CMatrix geterrorsCube(unsigned index);
	
	void seterrorsCube(CMatrix ErrorsMatrix, unsigned index);	
	
	void seterrorsCubeValues(CMatrix ErrorsMatrix, unsigned index);
	
	void seterrorsCubeValues(float ErrorValue, unsigned i, unsigned j, unsigned index);
	
	float geterrorsCubeValues(unsigned i, unsigned j, unsigned index) const;
	
	float *getdistdVector(void);
	
	float getdistd(unsigned index) const;
	
	void setdistdVector(float *Distd);
	
	void setdistd(float DistdValue, unsigned index);
	
	void setCubes(unsigned NDataPoints);	
	
	void setipPolyModel(PolynomialMatrix IPPolyModel);
	
	PolynomialMatrix getipPolyModel(void);
	
	PolynomialCoeffs_t* getipPolyModelCoefficients(unsigned i,unsigned j);
	
	const PolynomialCoeffs_t* getipPolyModelCoefficients(unsigned i,unsigned j) const;

	void setipPolyModelCoefficients(PolynomialCoeffs_t PolyModelCoeffs,unsigned i,unsigned j);

	CMatrix getipDataFromPolyModel(float d, CMatrix output);

    CMatrix getipDataFromPolyModel(PixelSet *aperturePixels, float d, CMatrix output);

	CMatrix getipDataFromPolyModel(unsigned index, CMatrix output);
    
    CMatrix getipDataFromPolyModel(PixelSet *aperturePixels, unsigned index, CMatrix output);

	float getipDataFromPolyModel(float d, unsigned i, unsigned j) const;
	
	float getipDataFromPolyModel(unsigned i, unsigned j, unsigned index) const;

	void setipDataFromPolyModel(unsigned index);	
    
    void setipDataFromPolyModel(PixelSet *aperturePixels, unsigned index);  

	void setipDataFromPolyModel(float d, unsigned index);
    
    void setipDataFromPolyModel(PixelSet *aperturePixels, float d, unsigned index);
    
    void setdataCubeFromPolyModel(PixelSet *aperturePixels);

    void setdataCubeFromPolyModel(PixelSet *aperturePixels, unsigned index);
    
    void FitPolyMatrixtoIPDataVector(unsigned coeffs, bool witherrors);		

    void FitPolyMatrixtoIPDataVector(PixelSet *aperturePixels, unsigned coeffs, bool witherrors);    
    
    void FitMediantoIPDataVector(void);
    
    void FitMediantoIPDataVector(PixelSet *aperturePixels);
    
	CMatrix getchisqrMatrix(void);
	
	void setchisqrMatrix(CMatrix ChisqrMatrix);
	
	float getchisqrMatrixValue(unsigned i, unsigned j) const;
	
	void setchisqrMatrixValue(float Value, unsigned i, unsigned j);

	void setdataCubeFromPolyModel(void);
    
    void setdataCubeFromPolyModel(unsigned index);
    
	void deleteDataCubes(void);  
    
    float getIPphotoCenterX(float d) const;
    
    float getIPphotoCenterY(float d) const;
    
    float getIPphotoCenterX(PixelSet *aperturePixels, float d) const;
    
    float getIPphotoCenterY(PixelSet *aperturePixels, float d) const;
    
    CMatrix getdataCubeValues(PixelSet *aperturePixels, unsigned index, CMatrix dataSlice);
    
    CMatrix getdataCubeValues(unsigned index, CMatrix dataSlice);    
	
};

#endif

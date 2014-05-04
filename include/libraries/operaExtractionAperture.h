#ifndef OPERAEXTRACTIONAPERTURE_H
#define OPERAEXTRACTIONAPERTURE_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaExtractionAperture
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

#ifndef CIRCLE
#define CIRCLE 999
#endif
#ifndef RECTANGLE
#define RECTANGLE 998
#endif
#ifndef POLYGON
#define POLYGON 997
#endif
#ifndef LINE
#define LINE 996
#endif
#ifndef PIXELSET
#define PIXELSET 995
#endif

#include <ostream>
#include <math.h>

#include "libraries/operaLibCommon.h"			// for MAXNPOLYGONSIDES
#include "libraries/operaStats.h"				// for operaArrayIndexSort

#include "libraries/operaGeometricShapes.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaInstrumentProfile.h"
#include "libraries/PixelSet.h"

enum apertureShape {circle=CIRCLE, rectangle=RECTANGLE, polygon=POLYGON, line=LINE, pixelset=PIXELSET};

class operaSpectralOrder;

/*! 
 * \sa class operaExtractionAperture
 * \brief The Extraction Aperture.
 * \details The Extraction Aperture defines the set of pixels or subpixels where  
 * \details the flux information is extracted from. The aperture shape can be defined 
 * \details as one of the available options.
 * \return none
 * \file operaExtractionAperture.h
 * \ingroup libraries
 */
class operaExtractionAperture {
	
private:    
    apertureShape shape;	// (CIRCLE,RECTANGLE,POLYGON,LINE,PIXELSET)
	
    Circle circleAperture;         // Circle operaGeometricShape
    Rectangle rectangleAperture;   // Rectangle operaGeometricShape 
    Polygon polygonAperture;       // Polygon operaGeometricShape
    Line lineAperture;             // Line operaGeometricShape
    
    Rectangle boundingBox;  // Box size that englobes aperture shape
    
	unsigned xsampling;		// sampling factor for the extraction aperture in x direction
	unsigned ysampling;		// sampling factor for the extraction aperture in y direction
    
    PixelSet *subpixels;    // set of subpixels with coordinates, values, and area of each subpixel
    
    float fluxFraction;     // fraction of flux (with respect to the entire IP) contained in aperture  
    
public:
	
	/*
	 * Constructor
	 */
	
    operaExtractionAperture(void);
    
    operaExtractionAperture(Circle *CircleAperture, unsigned XSampling, unsigned YSampling);
    
    operaExtractionAperture(Circle *CircleAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image);

    operaExtractionAperture(Circle *CircleAperture, operaInstrumentProfile *instrumentProfile, float distd);

    operaExtractionAperture(Circle *CircleAperture, operaInstrumentProfile *instrumentProfile);
    
    operaExtractionAperture(Rectangle *RectangleAperture, unsigned XSampling, unsigned YSampling); 
    
    operaExtractionAperture(Rectangle *RectangleAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image);
    
    operaExtractionAperture(Rectangle *RectangleAperture, operaInstrumentProfile *instrumentProfile, float distd);

    operaExtractionAperture(Rectangle *RectangleAperture, operaInstrumentProfile *instrumentProfile);

    operaExtractionAperture(Polygon *PolygonAperture, unsigned XSampling, unsigned YSampling);
    
    operaExtractionAperture(Polygon *PolygonAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image);
    
    operaExtractionAperture(Polygon *PolygonAperture, operaInstrumentProfile *instrumentProfile, float distd);
    
    operaExtractionAperture(Polygon *PolygonAperture, operaInstrumentProfile *instrumentProfile);    
    
    operaExtractionAperture(Line *LineAperture, unsigned XSampling, unsigned YSampling);
    
    operaExtractionAperture(Line *LineAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image);
    
    operaExtractionAperture(Line *LineAperture, operaInstrumentProfile *instrumentProfile, float distd);    
        
    operaExtractionAperture(Line *LineAperture, operaInstrumentProfile *instrumentProfile);    
#if 0    
    operaExtractionAperture(PixelSet *Subpixels, unsigned XSampling, unsigned YSampling);
#endif    
	
    /*
	 * Destructor
	 */
	
    ~operaExtractionAperture(void);
	
	/*
	 * Setters/Getters
	 */      
    
	void operaSetExtractionAperture(Line *LineAperture, unsigned XSampling, unsigned YSampling);
	
    void setShape(apertureShape Shape);  
    
    apertureShape getShape(void);      
    
	PixelSet* getSubpixels(void);     
    
    void setCircleAperture(Circle *CircleAperture);  
    
    Circle *getCircleAperture(void);     
    
    void setRectangleAperture(Rectangle *RectangleAperture);  
    
    Rectangle* getRectangleAperture(void); 
    
    void setPolygonAperture(Polygon *PolygonAperture);  
    
    Polygon* getPolygonAperture(void); 
    
    void setLineAperture(Line *LineAperture);  
    
    Line* getLineAperture(void);     
    
    void setBoundingBox(Rectangle *BoundingBox);  
    
    Rectangle *getBoundingBox(void);
    
    void setSampling(unsigned Xsampling, unsigned Ysampling);  
    
    unsigned getXsampling(void);
    
    unsigned getYsampling(void);

    void setSubpixels(void);    
    
    void setSubpixelsWithRedundancy(void);    
    
    void setSubpixels(PixelSet *Subpixels); 
    
    void setSubpixels(operaFITSImage &Image); 
    
    void setSubpixels(int naxis1, int naxis2);
    
    void setSubpixelsWithRedundancy(int naxis1, int naxis2);
    
    void setSubpixels(operaInstrumentProfile *instrumentProfile, float d); 
    
    void setSubpixels(operaInstrumentProfile *instrumentProfile);     

    float getFluxFraction(void);
    
    void setFluxFraction(float FluxFraction);
    
	/*
	 * Methods
	 */    

    void shiftAperture(float xshift, float yshift); 
    
    void shiftAperture(float xshift, float yshift, operaFITSImage &Image);
    
    void recenterAperture(operaPoint &NewCenter);
    
    void recenterAperture(operaPoint &NewCenter, operaFITSImage &Image); 
    
    void recenterAperture(operaPoint &NewCenter,int naxis1, int naxis2);    
    
    void recenterApertureWithRedundancy(operaPoint &NewCenter);     
    
    void recenterApertureWithRedundancy(operaPoint &NewCenter,int naxis1, int naxis2);    
    
	unsigned calculateMAXnSubPixels(void);    
    
    unsigned calculateMAXnSubPixels(int naxis1, int naxis2);
    
    unsigned calculateMAXnSubPixels(operaFITSImage &Image);
    
    unsigned calculateMAXnSubPixels(operaInstrumentProfile *instrumentProfile);    
	
};

#endif

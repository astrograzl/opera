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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaExtractionAperture.h"
#include "libraries/PixelSet.h"
#include "libraries/operaStats.h"

/*!
 * operaExtractionAperture
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the Extraction Aperture.
 * \file operaExtractionAperture.cpp
 * \ingroup libraries
 */

using namespace std;

/* 
 * \class operaExtractionAperture
 * \brief The Extraction Aperture defines the set of pixels or subpixels where  
 * \brief the flux information is extracted from
 * \return none
 */

/*
 * operaExtractionAperture Constructor
 */

operaExtractionAperture::operaExtractionAperture(void) :
shape(pixelset),
xsampling(1),
ysampling(1)
{
    
}

operaExtractionAperture::operaExtractionAperture(Circle *CircleAperture, unsigned XSampling, unsigned YSampling) :
shape(circle),
xsampling(1),
ysampling(1)
{
    setCircleAperture(CircleAperture);         
    xsampling = XSampling;
    ysampling = YSampling;
	
    float radius = circleAperture.getRadius();
    operaPoint *center = circleAperture.getCenter();
	
    boundingBox.setRectangle(2*radius, 2*radius, 0.0, *center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(); 
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels();
}

operaExtractionAperture::operaExtractionAperture(Circle *CircleAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image) :
shape(circle),
xsampling(1),
ysampling(1)
{
    setCircleAperture(CircleAperture);         
    xsampling = XSampling;
    ysampling = YSampling;
    
    float radius = circleAperture.getRadius();
    operaPoint *center = circleAperture.getCenter();
    
    boundingBox.setRectangle(2*radius, 2*radius, 0.0, *center);    
	
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels();
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(Image);
}

operaExtractionAperture::operaExtractionAperture(Circle *CircleAperture, operaInstrumentProfile *instrumentProfile, float distd) :
shape(circle),
xsampling(1),
ysampling(1)
{
    setCircleAperture(CircleAperture);         
    xsampling = instrumentProfile->getXsampling();
    ysampling = instrumentProfile->getYsampling();
    
    float radius = circleAperture.getRadius();
    operaPoint *center = circleAperture.getCenter();
    
    boundingBox.setRectangle(2*radius, 2*radius, 0.0, *center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(instrumentProfile,distd);
}

operaExtractionAperture::operaExtractionAperture(Circle *CircleAperture, operaInstrumentProfile *instrumentProfile) :
shape(circle),
xsampling(1),
ysampling(1)
{
    setCircleAperture(CircleAperture);         
    xsampling = instrumentProfile->getXsampling();
    ysampling = instrumentProfile->getYsampling();
    
    float radius = circleAperture.getRadius();
    operaPoint *center = circleAperture.getCenter();
    
    boundingBox.setRectangle(2*radius, 2*radius, 0.0, *center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(instrumentProfile);
}

operaExtractionAperture::operaExtractionAperture(Rectangle *RectangleAperture, unsigned XSampling, unsigned YSampling) :
shape(rectangle),
xsampling(1),
ysampling(1)
{
    setRectangleAperture(RectangleAperture); 
    xsampling = XSampling;
    ysampling = YSampling;      
    
    operaPoint *center = RectangleAperture->getCenter();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<FOURSIDES; i++) {
        if(RectangleAperture->getCorner(i)->getXcoord() < minx) {
            minx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getXcoord() > maxx) {
            maxx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() < miny) {
            miny = RectangleAperture->getCorner(i)->getYcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() > maxy) {
            maxy = RectangleAperture->getCorner(i)->getYcoord();
        }        
	}
	
    boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,*center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(); 
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels();    
}

operaExtractionAperture::operaExtractionAperture(Rectangle *RectangleAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image) :
shape(rectangle),
xsampling(1),
ysampling(1)
{
    setRectangleAperture(RectangleAperture); 
    xsampling = XSampling;
    ysampling = YSampling;      
    
    operaPoint *center = RectangleAperture->getCenter();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<FOURSIDES; i++) {
        if(RectangleAperture->getCorner(i)->getXcoord() < minx) {
            minx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getXcoord() > maxx) {
            maxx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() < miny) {
            miny = RectangleAperture->getCorner(i)->getYcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() > maxy) {
            maxy = RectangleAperture->getCorner(i)->getYcoord();
        }        
    }
    
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,*center);    
	
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels();
	
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(Image);  
}

operaExtractionAperture::operaExtractionAperture(Rectangle *RectangleAperture, operaInstrumentProfile *instrumentProfile, float distd) :
shape(rectangle),
xsampling(1),
ysampling(1)
{
    setRectangleAperture(RectangleAperture); 
    xsampling = instrumentProfile->getXsampling();
    ysampling = instrumentProfile->getYsampling();      
    
    operaPoint *center = RectangleAperture->getCenter();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<FOURSIDES; i++) {
        if(RectangleAperture->getCorner(i)->getXcoord() < minx) {
            minx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getXcoord() > maxx) {
            maxx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() < miny) {
            miny = RectangleAperture->getCorner(i)->getYcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() > maxy) {
            maxy = RectangleAperture->getCorner(i)->getYcoord();
        }        
	}
	
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,*center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(instrumentProfile,distd);  
}

operaExtractionAperture::operaExtractionAperture(Rectangle *RectangleAperture, operaInstrumentProfile *instrumentProfile) :
shape(rectangle),
xsampling(1),
ysampling(1)
{
    setRectangleAperture(RectangleAperture); 
    xsampling = instrumentProfile->getXsampling();
    ysampling = instrumentProfile->getYsampling();      
    
    operaPoint *center = RectangleAperture->getCenter();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<FOURSIDES; i++) {
        if(RectangleAperture->getCorner(i)->getXcoord() < minx) {
            minx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getXcoord() > maxx) {
            maxx = RectangleAperture->getCorner(i)->getXcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() < miny) {
            miny = RectangleAperture->getCorner(i)->getYcoord();
        }
        if(RectangleAperture->getCorner(i)->getYcoord() > maxy) {
            maxy = RectangleAperture->getCorner(i)->getYcoord();
        }        
	}
	
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,*center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(instrumentProfile);  
}


operaExtractionAperture::operaExtractionAperture(Polygon *PolygonAperture, unsigned XSampling, unsigned YSampling) :
shape(polygon),
xsampling(1),
ysampling(1)
{
    setPolygonAperture(PolygonAperture);
    xsampling = XSampling;
    ysampling = YSampling;      
    
    unsigned nSides = PolygonAperture->getNSides();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<nSides; i++) {
        if(PolygonAperture->getVertex(i)->getXcoord() < minx) {
            minx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getXcoord() > maxx) {
            maxx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() < miny) {
            miny = PolygonAperture->getVertex(i)->getYcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() > maxy) {
            maxy = PolygonAperture->getVertex(i)->getYcoord();
        }        
    }
    
    operaPoint center(minx + fabs(maxx-minx)/2 , miny + fabs(maxy-miny)/2);
    
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(); 
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels();    
}

operaExtractionAperture::operaExtractionAperture(Polygon *PolygonAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image) :
shape(polygon),
xsampling(1),
ysampling(1)
{
    setPolygonAperture(PolygonAperture);
    xsampling = XSampling;
    ysampling = YSampling;      
    
    unsigned nSides = PolygonAperture->getNSides();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<nSides; i++) {
        if(PolygonAperture->getVertex(i)->getXcoord() < minx) {
            minx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getXcoord() > maxx) {
            maxx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() < miny) {
            miny = PolygonAperture->getVertex(i)->getYcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() > maxy) {
            maxy = PolygonAperture->getVertex(i)->getYcoord();
        }        
    }
    
    operaPoint center(minx + fabs(maxx-minx)/2 , miny + fabs(maxy-miny)/2);
    
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels();
	
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(Image);      
}

operaExtractionAperture::operaExtractionAperture(Polygon *PolygonAperture, operaInstrumentProfile *instrumentProfile, float distd) :
shape(polygon),
xsampling(1),
ysampling(1)
{
    setPolygonAperture(PolygonAperture);
	xsampling = instrumentProfile->getXsampling();
	ysampling = instrumentProfile->getYsampling();    
	
    unsigned nSides = PolygonAperture->getNSides();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<nSides; i++) {
        if(PolygonAperture->getVertex(i)->getXcoord() < minx) {
            minx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getXcoord() > maxx) {
            maxx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() < miny) {
            miny = PolygonAperture->getVertex(i)->getYcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() > maxy) {
            maxy = PolygonAperture->getVertex(i)->getYcoord();
        }        
    }
    
    operaPoint center(minx + fabs(maxx-minx)/2 , miny + fabs(maxy-miny)/2);
    
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);
    
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    
    setSubpixels(instrumentProfile, distd);      
}

operaExtractionAperture::operaExtractionAperture(Polygon *PolygonAperture, operaInstrumentProfile *instrumentProfile) :
shape(polygon),
xsampling(1),
ysampling(1)
{
    setShape(polygon);    
    setPolygonAperture(PolygonAperture);
	xsampling = instrumentProfile->getXsampling();
	ysampling = instrumentProfile->getYsampling();    
    
    unsigned nSides = PolygonAperture->getNSides();
    
    float minx = BIG;
    float maxx = -BIG;
    float miny = BIG;
    float maxy = -BIG;
    
    for (unsigned i=0; i<nSides; i++) {
        if(PolygonAperture->getVertex(i)->getXcoord() < minx) {
            minx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getXcoord() > maxx) {
            maxx = PolygonAperture->getVertex(i)->getXcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() < miny) {
            miny = PolygonAperture->getVertex(i)->getYcoord();
        }
        if(PolygonAperture->getVertex(i)->getYcoord() > maxy) {
            maxy = PolygonAperture->getVertex(i)->getYcoord();
        }        
    }
    
    operaPoint center(minx + fabs(maxx-minx)/2 , miny + fabs(maxy-miny)/2);
    
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,center);    
    
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    
    subpixels = new PixelSet(SubPixelArea);
    
    setSubpixels(instrumentProfile);      
}

void operaExtractionAperture::operaSetExtractionAperture(Line *LineAperture, unsigned XSampling, unsigned YSampling)
{
    setShape(line);    
    setLineAperture(LineAperture);
    xsampling = XSampling;
    ysampling = YSampling;    
	
    Line topline(*LineAperture, line_top);
    Line bottomline(*LineAperture, line_bottom);
    Line leftline(*LineAperture, line_left);
    Line rightline(*LineAperture, line_right);
    
    float maxx = lineAperture.getMidPoint()->getXcoord() + lineAperture.getLength()/2;
    float minx = lineAperture.getMidPoint()->getXcoord() - lineAperture.getLength()/2;
    float maxy = lineAperture.getMidPoint()->getYcoord() + lineAperture.getWidth()/2;
    float miny = lineAperture.getMidPoint()->getYcoord() - lineAperture.getWidth()/2;    
	
    operaPoint Corner1(rightline, bottomline);
    operaPoint Corner2(topline, rightline);    
    operaPoint Corner3(leftline, topline);
    operaPoint Corner4(bottomline, leftline);     
    
    if(lineAperture.getSlope() > 0) {
        maxx = Corner1.getXcoord();
        maxy = Corner2.getYcoord();
        minx = Corner3.getXcoord();
        miny = Corner4.getYcoord();
    } else if (lineAperture.getSlope() < 0){
        maxx = Corner2.getXcoord();
        maxy = Corner3.getYcoord();
        minx = Corner4.getXcoord();
        miny = Corner1.getYcoord();        
    }
	boundingBox.setRectangle(fabs(maxx-minx),fabs(maxy-miny),0.0,*(lineAperture.getMidPoint()));    
}

operaExtractionAperture::operaExtractionAperture(Line *LineAperture, unsigned XSampling, unsigned YSampling) :
shape(line),
xsampling(1),
ysampling(1)
{
	operaSetExtractionAperture(LineAperture, XSampling, YSampling);    
    unsigned MAXnPixels = calculateMAXnSubPixels(); 
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    setSubpixels();      
}

operaExtractionAperture::operaExtractionAperture(Line *LineAperture, unsigned XSampling, unsigned YSampling, operaFITSImage &Image) 
{
	operaSetExtractionAperture(LineAperture, XSampling, YSampling);    
    unsigned MAXnPixels = calculateMAXnSubPixels(Image);
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
	setSubpixels(Image);  
}

operaExtractionAperture::operaExtractionAperture(Line *LineAperture, operaInstrumentProfile *instrumentProfile, float distd)
{
	operaSetExtractionAperture(LineAperture, instrumentProfile->getXsampling(), instrumentProfile->getYsampling());    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    setSubpixels(instrumentProfile, distd);  
}

operaExtractionAperture::operaExtractionAperture(Line *LineAperture, operaInstrumentProfile *instrumentProfile)
{
	operaSetExtractionAperture(LineAperture, instrumentProfile->getXsampling(), instrumentProfile->getYsampling());    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);
    float SubPixelArea = 1.0/(float)(getXsampling()*getYsampling());
    subpixels = new PixelSet(MAXnPixels,SubPixelArea);
    setSubpixels(instrumentProfile);  
}

/*
 * operaExtractionAperture Destructor
 */
operaExtractionAperture::~operaExtractionAperture(void) {
	if (subpixels)
		delete subpixels;
	subpixels = NULL;
}

/*
 * operaExtractionAperture Setters/Getters
 */

void operaExtractionAperture::setShape(apertureShape Shape) {
    shape = Shape;
}

apertureShape operaExtractionAperture::getShape(void) const {
    return shape;
}

void operaExtractionAperture::setSubpixels(void) {
    
    Rectangle *box = getBoundingBox();
	
    int i0 = (int)floor(box->getCenter()->getXcoord() - ceil(box->getWidth()/2));
    int nx = (int)ceil(box->getCenter()->getXcoord() + ceil(box->getWidth()/2));
    int j0 = (int)floor(box->getCenter()->getYcoord() - ceil(box->getHeight()/2));
    int ny = (int)ceil(box->getCenter()->getYcoord() + ceil(box->getHeight()/2));
    
    unsigned MAXnPixels = calculateMAXnSubPixels();
    
    subpixels->resize(MAXnPixels);    
    
    unsigned nPixels = 0;
	
    switch (shape) {
        case circle:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                             
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getCircleAperture()->pointInCircle(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case rectangle:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                            
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getRectangleAperture()->pointInRectangle(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case polygon:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                             
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getPolygonAperture()->pointInPolygon(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case line:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());  
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getLineAperture()->pointInLineWidth(TestPoint,*getLineAperture())) {  
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }            
            break;                    
        default:
            break;
    }
	subpixels->setNPixels(nPixels);
}

void operaExtractionAperture::setSubpixels(PixelSet *Subpixels) {
    subpixels->setSubPixels(*Subpixels);
}

void operaExtractionAperture::setSubpixelsWithRedundancy(void) {
    
    Rectangle *box = getBoundingBox();
	
    int i0 = (int)floor(box->getCenter()->getXcoord() - ceil(box->getWidth()/2));
    int nx = (int)ceil(box->getCenter()->getXcoord() + ceil(box->getWidth()/2));
    int j0 = (int)floor(box->getCenter()->getYcoord() - ceil(box->getHeight()/2));
    int ny = (int)ceil(box->getCenter()->getYcoord() + ceil(box->getHeight()/2));
    
    unsigned MAXnPixels = calculateMAXnSubPixels();
    
    subpixels->resize(MAXnPixels);    
    
    unsigned nPixels = 0;
	
    switch (shape) {
        case circle:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getCircleAperture()->pointInCircle(TestPoint)) {                               
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }
            break;
        case rectangle:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getRectangleAperture()->pointInRectangle(TestPoint)) {
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }            
            break;
        case polygon:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getPolygonAperture()->pointInPolygon(TestPoint)) {
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }               
            break;
        case line:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getLineAperture()->pointInLineWidth(TestPoint,*getLineAperture())) {                                 
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }                         
            break;                    
        default:
            break;
    } 
	subpixels->setNPixels(nPixels);
}

void operaExtractionAperture::setSubpixels(operaFITSImage &Image) {
    
    Rectangle *box = getBoundingBox();
    
    int i0=0,j0=0,nx=(int)Image.getnaxis1(),ny=(int)Image.getnaxis2();
    
    int ii0 = (int)floor(box->getCenter()->getXcoord() - ceil(box->getWidth()/2));
    if(ii0 < 0) {
        i0 = 0;
    } else if (ii0 >= 0 && ii0 <= (int)Image.getnaxis1()) {
        i0 = ii0;
    } else if (ii0 > (int)Image.getnaxis1()) {
        i0 = (int)Image.getnaxis1();
    }
        
    int inx = (int)ceil(box->getCenter()->getXcoord() + ceil(box->getWidth()/2));
    if(inx < 0) {
        nx = 0;
    } else if (inx >= 0 && inx <= (int)Image.getnaxis1()){
        nx = inx;  
    } else if (inx > (int)Image.getnaxis1()) {
        nx = (int)Image.getnaxis1();
    }
    
    int ij0 = (int)floor(box->getCenter()->getYcoord() - ceil(box->getHeight()/2));
    if(ij0 < 0) {
        j0 = 0; 
    } else if(ij0 >= 0 && ij0 <= (int)Image.getnaxis2()) {
        j0 = ij0;
    } else if (ij0 > (int)Image.getnaxis2()) {
        j0 = (int)Image.getnaxis2();
    }
    
    int iny = (int)ceil(box->getCenter()->getYcoord() + ceil(box->getHeight()/2));
    if(iny < 0) {
        ny = 0; 
    } else if (iny >= 0 && iny <= (int)Image.getnaxis2()){
        ny = iny;  
    } else if(iny > (int)Image.getnaxis2()) {
        ny = (int)Image.getnaxis2(); 
    }
    
    unsigned MAXnPixels = calculateMAXnSubPixels(Image);
    
	subpixels->resize(MAXnPixels);
    
    unsigned nPixels = 0;
    
    switch (shape) {
        case circle:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getCircleAperture()->pointInCircle(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                subpixels->setPixelValue((float)Image[j][i], nPixels);
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case rectangle:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getRectangleAperture()->pointInRectangle(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                subpixels->setPixelValue((float)Image[j][i], nPixels);
								nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case polygon:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                          
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getPolygonAperture()->pointInPolygon(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                subpixels->setPixelValue((float)Image[j][i], nPixels);
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case line:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                            
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getLineAperture()->pointInLineWidth(TestPoint,*getLineAperture())) {                                
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                subpixels->setPixelValue((float)Image[j][i], nPixels);
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }            
            break;                    
        default:
            break;
    } 
	subpixels->setNPixels(nPixels);
}

void operaExtractionAperture::setSubpixels(int naxis1, int naxis2)  {
    
    Rectangle *box = getBoundingBox();
    
    int i0=0,j0=0,nx=naxis1,ny=naxis2;
    int ii0 = (int)floor(box->getCenter()->getXcoord() - ceil(box->getWidth()/2));
    
    if(ii0 < 0) {
        i0 = 0;
    } else if (ii0 >= 0 && ii0 <= naxis1) {
        i0 = ii0;
    } else if (ii0 > naxis1) {
        i0 = naxis1;
    }
    
    int inx = (int)ceil(box->getCenter()->getXcoord() + ceil(box->getWidth()/2));
    if(inx < 0) {
        nx = 0;
    } else if (inx >= 0 && inx <= naxis1){
        nx = inx;  
    } else if (inx > naxis1) {
        nx = naxis1;
    }
    
    int ij0 = (int)floor(box->getCenter()->getYcoord() - ceil(box->getHeight()/2));
    if(ij0 < 0) {
        j0 = 0; 
    } else if(ij0 >= 0 && ij0 <= naxis2) {
        j0 = ij0;
    } else if (ij0 > naxis2) {
        j0 = naxis2;
    }
    
    int iny = (int)ceil(box->getCenter()->getYcoord() + ceil(box->getHeight()/2));
    if(iny < 0) {
        ny = 0; 
    } else if (iny >= 0 && iny <= naxis2){
        ny = iny;  
    } else if(iny > naxis2) {
        ny = naxis2; 
    }
    
    unsigned MAXnPixels = calculateMAXnSubPixels(naxis1,naxis2);
    
	subpixels->resize(MAXnPixels);
    
    unsigned nPixels = 0;
    
    switch (shape) {
        case circle:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getCircleAperture()->pointInCircle(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case rectangle:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getRectangleAperture()->pointInRectangle(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
								nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case polygon:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                          
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getPolygonAperture()->pointInPolygon(TestPoint)) {
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }
            break;
        case line:
            for(int i=i0;i<nx;i++){
                for(unsigned ii=0;ii<getXsampling();ii++) {
                    float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                    for(int j=j0;j<ny;j++){
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                            
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getLineAperture()->pointInLineWidth(TestPoint,*getLineAperture())) {                                
                                subpixels->setXcenter(xsubpix, nPixels);
                                subpixels->setYcenter(ysubpix, nPixels);
                                subpixels->setiIndex(i, nPixels);
                                subpixels->setjIndex(j, nPixels);
                                subpixels->setredundancy(1,nPixels);                                
                                nPixels++;
                            }
                            
                        }
                    }
                }
            }            
            break;                    
        default:
            break;
    } 
	subpixels->setNPixels(nPixels);
}

void operaExtractionAperture::setSubpixelsWithRedundancy(int naxis1, int naxis2) {
    
    Rectangle *box = getBoundingBox();
    
    int i0=0,j0=0,nx=naxis1,ny=naxis2;
    int ii0 = (int)floor(box->getCenter()->getXcoord() - ceil(box->getWidth()/2));
    
    if(ii0 < 0) {
        i0 = 0;
    } else if (ii0 >= 0 && ii0 <= naxis1) {
        i0 = ii0;
    } else if (ii0 > naxis1) {
        i0 = naxis1;
    }
    
    int inx = (int)ceil(box->getCenter()->getXcoord() + ceil(box->getWidth()/2));
    if(inx < 0) {
        nx = 0;
    } else if (inx >= 0 && inx <= naxis1){
        nx = inx;  
    } else if (inx > naxis1) {
        nx = naxis1;
    }
    
    int ij0 = (int)floor(box->getCenter()->getYcoord() - ceil(box->getHeight()/2));
    if(ij0 < 0) {
        j0 = 0; 
    } else if(ij0 >= 0 && ij0 <= naxis2) {
        j0 = ij0;
    } else if (ij0 > naxis2) {
        j0 = naxis2;
    }
    
    int iny = (int)ceil(box->getCenter()->getYcoord() + ceil(box->getHeight()/2));
    if(iny < 0) {
        ny = 0; 
    } else if (iny >= 0 && iny <= naxis2){
        ny = iny;  
    } else if(iny > naxis2) {
        ny = naxis2; 
    }

    unsigned MAXnPixels = calculateMAXnSubPixels(naxis1,naxis2);
    
	subpixels->resize(MAXnPixels);
    
    unsigned nPixels = 0;
    
    switch (shape) {
        case circle:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getCircleAperture()->pointInCircle(TestPoint)) {                               
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }
            break;
        case rectangle:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getRectangleAperture()->pointInRectangle(TestPoint)) {
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }            
            break;
        case polygon:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getPolygonAperture()->pointInPolygon(TestPoint)) {
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }               
            break;
        case line:
            for(int i=i0;i<nx;i++){
                for(int j=j0;j<ny;j++){
                    unsigned redundancy = 0;
                    for(unsigned ii=0;ii<getXsampling();ii++) {
                        float xsubpix = (float)i + (float)(1+2*ii)/float(2*getXsampling());                
                        for(unsigned jj=0;jj<getYsampling();jj++) {
                            float ysubpix = (float)j + (float)(1+2*jj)/float(2*getYsampling());                           
                            operaPoint TestPoint(xsubpix,ysubpix);
                            if(getLineAperture()->pointInLineWidth(TestPoint,*getLineAperture())) {                                
                                redundancy++;
                            }                    
                        }
                    }
                    if(redundancy) {
                        subpixels->setXcenter((float)i+0.5, nPixels);
                        subpixels->setYcenter((float)j+0.5, nPixels);
                        subpixels->setiIndex(i, nPixels);
                        subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(redundancy,nPixels);
                        nPixels++;    
                    }
                }
            }                         
            break;                    
        default:
            break;
    } 
	subpixels->setNPixels(nPixels);
}


void operaExtractionAperture::setSubpixels(operaInstrumentProfile *instrumentProfile, float d) {
    
    Rectangle *box = getBoundingBox();
    
    int i0 = 0;
    if(box->getCenter()->getXcoord() - ceil(box->getWidth()/2) > instrumentProfile->getIPixXCoordinate(0)) {
        if((int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() - ceil(box->getWidth()/2)) > 0) {
            i0=(int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() - ceil(box->getWidth()/2)) - 1;
        }
    }
	
    int nx = (int)instrumentProfile->getNXPoints();
    if(box->getCenter()->getXcoord() + ceil(box->getWidth()/2) < instrumentProfile->getIPixXCoordinate(instrumentProfile->getNXPoints()-1)) {
        if((int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() + ceil(box->getWidth()/2)) + 1 < nx) {
            nx = (int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() + ceil(box->getWidth()/2)) + 1;
        }
    }
    
    int j0 = 0;
    if(box->getCenter()->getYcoord() - ceil(box->getHeight()/2) > instrumentProfile->getIPixYCoordinate(0)) {
        if((int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() - ceil(box->getHeight()/2)) > 0) {
            j0=(int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() - ceil(box->getHeight()/2)) - 1;
        }
    }
    
    int ny = (int)instrumentProfile->getNYPoints();
    if(box->getCenter()->getYcoord() + ceil(box->getHeight()/2) < instrumentProfile->getIPixYCoordinate(instrumentProfile->getNYPoints()-1)) {
        if((int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() + ceil(box->getHeight()/2)) + 1 < ny) {
            ny = (int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() + ceil(box->getHeight()/2)) + 1;
        }
    }    
    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);    
    
	subpixels->resize(MAXnPixels);
    
    unsigned nPixels = 0;
    
    switch (shape) {
        case circle:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                            
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getCircleAperture()->pointInCircle(TestPoint)) {                                
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
						subpixels->setPixelValue(instrumentProfile->getipDataFromPolyModel(d,(unsigned)i,(unsigned)j), nPixels);
                        nPixels++;
                    }
                }
            }   
            break;
        case rectangle:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                            
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getRectangleAperture()->pointInRectangle(TestPoint)) {                                
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
						subpixels->setPixelValue(instrumentProfile->getipDataFromPolyModel(d,(unsigned)i,(unsigned)j), nPixels);
                        nPixels++;
                    }
                }
            }   
            break;
        case polygon:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                         
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getPolygonAperture()->pointInPolygon(TestPoint)) {                                
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
						subpixels->setPixelValue(instrumentProfile->getipDataFromPolyModel(d,(unsigned)i,(unsigned)j), nPixels);
                        nPixels++;
                    }
                }
            }   
            break;
        case line:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                         
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getLineAperture()->pointInLineWidth(TestPoint,*getLineAperture())) {                                 
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
						subpixels->setPixelValue(instrumentProfile->getipDataFromPolyModel(d,(unsigned)i,(unsigned)j), nPixels);
                        nPixels++;
                    }
                }
            }            
            break;                    
        default:
            break;
    } 
	subpixels->setNPixels(nPixels);
}

void operaExtractionAperture::setSubpixels(operaInstrumentProfile *instrumentProfile) {
    
    Rectangle *box = getBoundingBox();
    
    int i0 = 0;
    if(box->getCenter()->getXcoord() - ceil(box->getWidth()/2) > instrumentProfile->getIPixXCoordinate(0)) {
        if((int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() - ceil(box->getWidth()/2)) > 0) {
            i0=(int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() - ceil(box->getWidth()/2)) - 1;
        }
    }
	
    int nx = (int)instrumentProfile->getNXPoints();
    if(box->getCenter()->getXcoord() + ceil(box->getWidth()/2) < instrumentProfile->getIPixXCoordinate(instrumentProfile->getNXPoints()-1)) {
        if((int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() + ceil(box->getWidth()/2)) + 1 < nx) {
            nx = (int)instrumentProfile->getIPixiIndex(box->getCenter()->getXcoord() + ceil(box->getWidth()/2)) + 1;
        }
    }
    
    int j0 = 0;
    if(box->getCenter()->getYcoord() - ceil(box->getHeight()/2) > instrumentProfile->getIPixYCoordinate(0)) {
        if((int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() - ceil(box->getHeight()/2)) > 0) {
            j0=(int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() - ceil(box->getHeight()/2)) - 1;
        }
    }
    
    int ny = (int)instrumentProfile->getNYPoints();
    if(box->getCenter()->getYcoord() + ceil(box->getHeight()/2) < instrumentProfile->getIPixYCoordinate(instrumentProfile->getNYPoints()-1)) {
        if((int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() + ceil(box->getHeight()/2)) + 1 < ny) {
            ny = (int)instrumentProfile->getIPixjIndex(box->getCenter()->getYcoord() + ceil(box->getHeight()/2)) + 1;
        }
    }    
    
    unsigned MAXnPixels = calculateMAXnSubPixels(instrumentProfile);    
    
	subpixels->resize(MAXnPixels);
    
    unsigned nPixels = 0;
    
    switch (shape) {
        case circle:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                                                   
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getCircleAperture()->pointInCircle(TestPoint)) {                                
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
                        nPixels++;
                    }
                }
            }   
            break;
        case rectangle:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                          
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getRectangleAperture()->pointInRectangle(TestPoint)) {                                
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
                        nPixels++;
                    }
                }
            }   
            break;
        case polygon:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                         
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getPolygonAperture()->pointInPolygon(TestPoint)) {                                
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
                        nPixels++;
                    }
                }
            }   
            break;
        case line:
            for(int i=i0;i<nx;i++){
                float xsubpix = instrumentProfile->getIPixXCoordinate((unsigned)i);                
                for(int j=j0;j<ny;j++){
                    float ysubpix = instrumentProfile->getIPixYCoordinate((unsigned)j);                           
                    operaPoint TestPoint(xsubpix,ysubpix);
                    if(getLineAperture()->pointInLineWidth(TestPoint,*getLineAperture())) {                                
						subpixels->setXcenter(xsubpix, nPixels);
						subpixels->setYcenter(ysubpix, nPixels);
						subpixels->setiIndex(i, nPixels);
						subpixels->setjIndex(j, nPixels);
                        subpixels->setredundancy(1,nPixels);                        
                        nPixels++;
                    }
                }
            }            
            break;                    
        default:
            break;
    } 
	subpixels->setNPixels(nPixels);
}

PixelSet* operaExtractionAperture::getSubpixels(void) {
    return subpixels;
}

const PixelSet* operaExtractionAperture::getSubpixels(void) const {
    return subpixels;
}

void operaExtractionAperture::setCircleAperture(Circle *CircleAperture) {
    circleAperture.setCircle(*CircleAperture);
} 

Circle* operaExtractionAperture::getCircleAperture(void) {
    return &circleAperture;
}     

const Circle* operaExtractionAperture::getCircleAperture(void) const {
    return &circleAperture;
}     

void operaExtractionAperture::setRectangleAperture(Rectangle *RectangleAperture) {
    rectangleAperture.setRectangle(*RectangleAperture);
}  

Rectangle* operaExtractionAperture::getRectangleAperture(void) {
    return &rectangleAperture;
} 

const Rectangle* operaExtractionAperture::getRectangleAperture(void) const {
    return &rectangleAperture;
} 

void operaExtractionAperture::setPolygonAperture(Polygon *PolygonAperture) {
    polygonAperture.setPolygon(*PolygonAperture);
}  

Polygon* operaExtractionAperture::getPolygonAperture(void) {
    return &polygonAperture;
} 

const Polygon* operaExtractionAperture::getPolygonAperture(void) const {
    return &polygonAperture;
} 

void operaExtractionAperture::setLineAperture(Line *LineAperture) {
    lineAperture.setLine(*LineAperture);
}  

Line* operaExtractionAperture::getLineAperture(void) {
    return &lineAperture;
}

const Line* operaExtractionAperture::getLineAperture(void) const {
    return &lineAperture;
}

void operaExtractionAperture::setBoundingBox(Rectangle *BoundingBox) {
    boundingBox.setRectangle(*BoundingBox);
}  

Rectangle *operaExtractionAperture::getBoundingBox(void) {
    return &boundingBox;
}

const Rectangle *operaExtractionAperture::getBoundingBox(void) const {
    return &boundingBox;
}

void operaExtractionAperture::setSampling(unsigned Xsampling, unsigned Ysampling) {
    xsampling = Xsampling;
    ysampling = Ysampling;    
}

unsigned operaExtractionAperture::getXsampling(void) const {
    return xsampling;
}

unsigned operaExtractionAperture::getYsampling(void) const {
    return ysampling;
}

float operaExtractionAperture::getFluxFraction(void) const {
    return fluxFraction;
}

void operaExtractionAperture::setFluxFraction(float FluxFraction) {
    fluxFraction = FluxFraction;
}

/*
 * operaExtractionAperture Methods
 */

void operaExtractionAperture::shiftAperture(float xshift, float yshift) {
    
    float boxWidth = boundingBox.getWidth();
    float boxHeight = boundingBox.getHeight(); 
    float boxAngle = 0.0;    
    
    if(shape == circle) {
        operaPoint newCenter(circleAperture.getCenter()->getXcoord() + xshift, circleAperture.getCenter()->getYcoord() + yshift);
        Circle newcircle(circleAperture.getRadius(), newCenter);                       
        circleAperture.setCircle(newcircle);
        
        boundingBox.setRectangle(boxWidth,boxHeight,boxAngle,newCenter); 
        
    } else if (shape == rectangle) {
        operaPoint newCenter(rectangleAperture.getCenter()->getXcoord() + xshift,rectangleAperture.getCenter()->getYcoord() + yshift);            
        Rectangle newrectangle(rectangleAperture.getWidth(),rectangleAperture.getHeight(),rectangleAperture.getAngle(),newCenter);
		rectangleAperture.setRectangle(newrectangle);
        
        boundingBox.setRectangle(boxWidth,boxHeight,boxAngle,newCenter);  		
        
    } else if (shape == polygon) {
        unsigned nsides = polygonAperture.getNSides();
        operaPoint polygoncorners[MAXNPOLYGONSIDES];
        for(unsigned i=0;i<nsides;i++) {
            polygoncorners[i].setPoint(polygonAperture.getVertex(i)->getXcoord() + xshift,polygonAperture.getVertex(i)->getYcoord() + yshift);
        }
        Polygon newpolygon(nsides, polygoncorners);
        
		polygonAperture.setPolygon(newpolygon);
        
        operaPoint newCenter(boundingBox.getCenter()->getXcoord() + xshift,boundingBox.getCenter()->getYcoord() + yshift);            
        
        boundingBox.setRectangle(boxWidth,boxHeight,boxAngle,newCenter);   
		
    } else if (shape == line) {
        operaPoint newMidPoint(lineAperture.getMidPoint()->getXcoord() + xshift,lineAperture.getMidPoint()->getYcoord() + yshift);            
        Line newline(lineAperture.getSlope(),lineAperture.getWidth(),lineAperture.getLength(),newMidPoint);
		lineAperture.setLine(newline);
        
        boundingBox.setRectangle(boxWidth,boxHeight,boxAngle,newMidPoint);         
    }
}

void operaExtractionAperture::shiftAperture(float xshift, float yshift, operaFITSImage &Image) {
    shiftAperture(xshift,yshift);
    setSubpixels(Image);    
}

void operaExtractionAperture::recenterAperture(operaPoint &NewCenter) {
    
    float boxWidth = boundingBox.getWidth();
    float boxHeight = boundingBox.getHeight(); 
    float boxAngle = 0.0;  
    
    boundingBox.setRectangle(boxWidth,boxHeight,boxAngle,NewCenter);
    
    if(shape == circle) {
        Circle newcircle(circleAperture.getRadius(),NewCenter);                       
        circleAperture.setCircle(newcircle);
    } else if (shape == rectangle) {
        Rectangle newrectangle(rectangleAperture.getWidth(),rectangleAperture.getHeight(),rectangleAperture.getAngle(),NewCenter);
		rectangleAperture.setRectangle(newrectangle);
        
    } else if (shape == polygon) {
        unsigned nsides = polygonAperture.getNSides();
        operaPoint polygoncorners[MAXNPOLYGONSIDES];
        
        float xshift = NewCenter.getXcoord() - polygonAperture.getVertex(0)->getXcoord();
        float yshift = NewCenter.getYcoord() - polygonAperture.getVertex(0)->getYcoord();
        
        polygoncorners[0].setPoint(NewCenter.getXcoord(), NewCenter.getYcoord());
        
        for(unsigned i=1;i<nsides;i++) {
            polygoncorners[i].setPoint(polygonAperture.getVertex(i)->getXcoord() + xshift,polygonAperture.getVertex(i)->getYcoord() + yshift);
        }
        Polygon newpolygon(nsides, polygoncorners);        
		polygonAperture.setPolygon(newpolygon);
    } else if (shape == line) {
        Line newline(lineAperture.getSlope(),lineAperture.getWidth(),lineAperture.getLength(),NewCenter);
		lineAperture.setLine(newline);
    }
}

void operaExtractionAperture::recenterAperture(operaPoint &NewCenter, operaFITSImage &Image) {
    recenterAperture(NewCenter);
    setSubpixels(Image);    
}

void operaExtractionAperture::recenterAperture(operaPoint &NewCenter,int naxis1, int naxis2) {
    recenterAperture(NewCenter);
    setSubpixels(naxis1,naxis2);
}

void operaExtractionAperture::recenterApertureWithRedundancy(operaPoint &NewCenter,int naxis1, int naxis2) {
    recenterAperture(NewCenter);
    setSubpixelsWithRedundancy(naxis1,naxis2);    
}

void operaExtractionAperture::recenterApertureWithRedundancy(operaPoint &NewCenter) {
    recenterAperture(NewCenter);
    setSubpixelsWithRedundancy();    
}

unsigned operaExtractionAperture::calculateMAXnSubPixels(void) {
    int i0 = (int)floor(boundingBox.getCenter()->getXcoord() - boundingBox.getWidth()/2);
    int nx = (int)ceil(boundingBox.getCenter()->getXcoord() + boundingBox.getWidth()/2);  
    int j0 = (int)floor(boundingBox.getCenter()->getYcoord() - boundingBox.getHeight()/2);
    int ny = (int)ceil(boundingBox.getCenter()->getYcoord() + boundingBox.getHeight()/2);
    
    return (unsigned)(fabs((float)(nx-i0))*fabs((float)(ny-j0)))*getXsampling()*getYsampling();
}

unsigned operaExtractionAperture::calculateMAXnSubPixels(int naxis1, int naxis2) {
    
    Rectangle *box = getBoundingBox();
    
    int i0=0,j0=0,nx=naxis1,ny=naxis2;
    
    int ii0 = (int)floor(box->getCenter()->getXcoord() - ceil(box->getWidth())/2);
    if(ii0 < 0) {
        i0 = 0;
    } else if (ii0 >= 0 && ii0 <= naxis1) {
        i0 = ii0;
    } else if (ii0 > naxis1) {
        i0 = naxis1;
    }
    
    int inx = (int)ceil(box->getCenter()->getXcoord() + ceil(box->getWidth())/2);
    if(inx < 0) {
        nx = 0;
    } else if (inx >= 0 && inx <= naxis1){
        nx = inx;  
    } else if (inx > naxis1) {
        nx = naxis1;
    }
    
    int ij0 = (int)floor(box->getCenter()->getYcoord() - ceil(box->getHeight())/2);
    if(ij0 < 0) {
        j0 = 0; 
    } else if(ij0 >= 0 && ij0 <= naxis2) {
        j0 = ij0;
    } else if (ij0 > naxis2) {
        j0 = naxis2;
    }
    
    int iny = (int)ceil(box->getCenter()->getYcoord() + ceil(box->getHeight())/2);
    if(iny < 0) {
        ny = 0; 
    } else if (iny >= 0 && iny <= naxis2){
        ny = iny;  
    } else if(iny > naxis2) {
        ny = naxis2; 
    }
    
    return (unsigned)(fabs((float)(nx-i0))*fabs((float)(ny-j0)))*getXsampling()*getYsampling();    
}

unsigned operaExtractionAperture::calculateMAXnSubPixels(operaFITSImage &Image) {
    return calculateMAXnSubPixels((int)Image.getnaxis1(), (int)Image.getnaxis2());
}

unsigned operaExtractionAperture::calculateMAXnSubPixels(operaInstrumentProfile *instrumentProfile) {
    return instrumentProfile->getNXPoints()*instrumentProfile->getNYPoints();
}

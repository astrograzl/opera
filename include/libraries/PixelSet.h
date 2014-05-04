#ifndef PIXELSET_H
#define PIXELSET_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: PixelSet
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

/*! 
 * \sa class PixelSet
 * \brief Encapsulation of the pixel set.
 * \details Defines the set of pixels or subpixels where  
 * \details the flux information is extracted from. 
 * \file PixelSet.h
 * \ingroup libraries
 */

class PixelSet {
private:
    unsigned nPixels;       // number of subpixels
    unsigned maxnPixels;    // maximum number of subpixels allocated
    float *xcenter;         // subpixel x-center positions
    float *ycenter;         // subpixel y-center positions
    int *iIndex;			// original image pixel i-index (col) from which the subpixel was obtained
    int *jIndex;			// original image pixel j-index (row) from which the subpixel was obtained
    
    unsigned *redundancy;   // redundancy  
    
    float *pixelValue;      // pixel value
    float subpixelArea;     // area of subpixel in pixel^2, value depends on pixelation
public:
	
	/*
	 * Constructor
	 */
    
    PixelSet(void);
    
    PixelSet(float SubPixelArea);  
    
    PixelSet(unsigned NPixels, float SubPixelArea);
    
    /*
	 * Destructor
	 */
    
    ~PixelSet(void);  
	
	/*
	 * Setters/Getters
	 */    
	
    unsigned getNPixels(void);
	unsigned getmaxNPixels(void);
    float getXcenter(unsigned index);
    float getYcenter(unsigned index); 
    int getiIndex(unsigned index); 
    int getjIndex(unsigned index); 
    unsigned getredundancy(unsigned index);     
    float getPixelValue(unsigned index);   
    void setNPixels(unsigned NPixels);
	void setXcenter(float Xcenter, unsigned index);
	void setYcenter(float Ycenter, unsigned index);
	void setiIndex(int iindex, unsigned k);
	void setjIndex(int jindex, unsigned k);
	void setredundancy(unsigned Redundancy, unsigned k);
	void setPixelValue(float PixelValue, unsigned index);
	
    float getSubpixelArea(void);
    void setSubpixelArea(float SubpixelArea);  
    
	void resize(unsigned newsize);
	void setSubPixels(PixelSet &pixels);
	/*
	 * Methods
	 */      
    void createVectors(unsigned NPixels);
    void deleteVectors(void);
    
    float getMinXcoord(void);
    float getMaxXcoord(void);
    float getMinYcoord(void);
    float getMaxYcoord(void);

};

#endif

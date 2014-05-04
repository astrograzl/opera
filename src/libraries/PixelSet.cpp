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

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaLib.h"		// for itos
#include "libraries/operaExtractionAperture.h"
#include "libraries/operaException.h"
#include "libraries/operaStats.h"
#include "libraries/PixelSet.h"

/*! 
 * PixelSet
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates a set of pixels or subpixels
 * \file PixelSet.cpp
 * \ingroup libraries
 */

using namespace std;


/* 
 * PixelSet
 * \brief It defines a set of pixels or subpixels where  
 * \brief the flux information is extracted from
 * \return none
 */

/*
 * PixelSet Constructor
 */

PixelSet::PixelSet(void) :
nPixels(0),
maxnPixels(0),
xcenter(NULL),
ycenter(NULL),
iIndex(NULL),
jIndex(NULL),
redundancy(NULL),
pixelValue(NULL),
subpixelArea(1)
{
    
}

PixelSet::PixelSet(float SubPixelArea) :
nPixels(0),
maxnPixels(0),
xcenter(NULL),
ycenter(NULL),
iIndex(NULL),
jIndex(NULL),
redundancy(NULL),
pixelValue(NULL),
subpixelArea(1)
{
    subpixelArea = SubPixelArea;
}

PixelSet::PixelSet(unsigned NPixels, float SubPixelArea) :
nPixels(0),
maxnPixels(0),
xcenter(NULL),
ycenter(NULL),
iIndex(NULL),
jIndex(NULL),
redundancy(NULL),
pixelValue(NULL),
subpixelArea(1)
{   
	if (NPixels == 0) {
		throw operaException("PixelSet: npixels == 0 ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    subpixelArea = SubPixelArea;
    createVectors(NPixels);
    nPixels = NPixels;         
    maxnPixels = NPixels;    
}

/*
 * PixelSet Destructor
 */
PixelSet::~PixelSet(void) {
    deleteVectors();
    nPixels = 0;
	maxnPixels = 0;
} 

/*
 * PixelSet Setters/Getters
 */

unsigned PixelSet::getNPixels(void) {
    return nPixels;
}
unsigned PixelSet::getmaxNPixels(void) {
    return maxnPixels;
}
float PixelSet::getXcenter(unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return xcenter[index];
}
float PixelSet::getYcenter(unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return ycenter[index];
}
int PixelSet::getiIndex(unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return iIndex[index];    
}
int PixelSet::getjIndex(unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return jIndex[index];    
}
unsigned PixelSet::getredundancy(unsigned index){
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return redundancy[index];    
}     
float PixelSet::getPixelValue(unsigned index) {
#ifdef RANGE_CHECK
    if (index >= nPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(nPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    return pixelValue[index];
}
float PixelSet::getSubpixelArea(void) {
    return subpixelArea;
}
void PixelSet::setNPixels(unsigned NPixels) {
#ifdef RANGE_CHECK
    if (NPixels > maxnPixels) {
		throw operaException("PixelSet: NPixels="+itos(NPixels)+" maxnPixels="+itos(maxnPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    nPixels = NPixels;
}

void PixelSet::setXcenter(float Xcenter, unsigned index) {
#ifdef RANGE_CHECK
    if (index >= maxnPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(maxnPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    xcenter[index] = Xcenter;
}

void PixelSet::setYcenter(float Ycenter, unsigned index) {
#ifdef RANGE_CHECK
    if (index >= maxnPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(maxnPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    ycenter[index] = Ycenter;
}

void PixelSet::setiIndex(int iindex, unsigned k) {
#ifdef RANGE_CHECK
    if (k >= maxnPixels) {
		throw operaException("PixelSet: k="+itos(k)+" nPixels="+itos(maxnPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    iIndex[k] = iindex;
}

void PixelSet::setjIndex(int jindex, unsigned k) {
#ifdef RANGE_CHECK
    if (k >= maxnPixels) {
		throw operaException("PixelSet: k="+itos(k)+" nPixels="+itos(maxnPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    jIndex[k] = jindex;
}

void PixelSet::setredundancy(unsigned Redundancy, unsigned k) {
#ifdef RANGE_CHECK
    if (k >= maxnPixels) {
		throw operaException("PixelSet: k="+itos(k)+" nPixels="+itos(maxnPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
    redundancy[k] = Redundancy;
}

void PixelSet::setPixelValue(float PixelValue, unsigned index) {
#ifdef RANGE_CHECK
    if (index >= maxnPixels) {
		throw operaException("PixelSet: index="+itos(index)+" nPixels="+itos(maxnPixels)+"  ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	pixelValue[index] = PixelValue;
}

void PixelSet::setSubpixelArea(float SubpixelArea) {
    subpixelArea = SubpixelArea;
}


/*
 * PixelSet Methods
 */

void PixelSet::createVectors(unsigned NPixels) {
	if (NPixels == 0) {
		throw operaException("PixelSet: NPixels == 0 ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	xcenter = new float[NPixels];
    ycenter = new float[NPixels];
    iIndex = new int[NPixels];
    jIndex = new int[NPixels];
    redundancy = new unsigned[NPixels];
    pixelValue = new float[NPixels];
    nPixels = NPixels;         
    maxnPixels = NPixels;         
}

void PixelSet::deleteVectors(void) {
	if (maxnPixels) {
		if(xcenter) {
			delete[] xcenter;
		}
		xcenter = NULL;
		if(ycenter) {
			delete[] ycenter;
		}    
		ycenter = NULL;
		if(pixelValue) {
			delete[] pixelValue;
		}
		pixelValue = NULL;
		if(iIndex) {
			delete[] iIndex;
		}
		iIndex = NULL;
		if(jIndex) {
			delete[] jIndex;
		}       
		jIndex = NULL;
		if(redundancy) {
			delete[] redundancy;
		}       
		redundancy = NULL;        

		nPixels = 0;
		maxnPixels = 0;
	}
}

void PixelSet::setSubPixels(PixelSet &pixels) {
	
	subpixelArea = pixels.subpixelArea;
	if (maxnPixels >= pixels.nPixels) { // we have enough room, just copy
		nPixels = pixels.nPixels;
		for (unsigned i=0; i<nPixels; i++) {
			setXcenter(pixels.getXcenter(i), i);
			setYcenter(pixels.getYcenter(i), i);
			setiIndex(pixels.getiIndex(i), i);
			setjIndex(pixels.getjIndex(i), i);
			setredundancy(pixels.getredundancy(i), i);                        
			setPixelValue(pixels.getPixelValue(i), i);
		}
	} else {
		deleteVectors();
		createVectors(pixels.nPixels);
		for (unsigned i=0; i<nPixels; i++) {
			setXcenter(pixels.getXcenter(i), i);
			setYcenter(pixels.getYcenter(i), i);
			setiIndex(pixels.getiIndex(i), i);
			setjIndex(pixels.getjIndex(i), i);
			setredundancy(pixels.getredundancy(i), i);              
			setPixelValue(pixels.getPixelValue(i), i);
		}
	}
}

void PixelSet::resize(unsigned newsize) {
	if (maxnPixels == 0) { // make more room
		createVectors(newsize);
	} else if (maxnPixels < newsize) { // make more room
		deleteVectors();
		createVectors(newsize);
	} else {
		nPixels = newsize;
	}
}

float PixelSet::getMinXcoord(void) {    
    return operaArrayMinValue(nPixels,xcenter);
}

float PixelSet::getMaxXcoord(void) {
    return operaArrayMaxValue(nPixels,xcenter);               
}

float PixelSet::getMinYcoord(void) {
    return operaArrayMinValue(nPixels,ycenter);               
}

float PixelSet::getMaxYcoord(void) {
    return operaArrayMaxValue(nPixels,ycenter);               
}

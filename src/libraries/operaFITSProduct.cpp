/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaFITSProduct
 Version: 1.0
 Description: class encapsulates a FITS image.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 Contact: teeple@cfht.hawaii.edu
 
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
 * operaFITSProduct
 * \author Doug Teeple
 * \brief This class encapsulates the OPERA products.
 * \file operaFITSProduct.cpp
 * \ingroup libraries
 */

#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaMEFFITSProduct.h"

using namespace std;

/*!
 * \verbatim
 *	   upena-compatible headers for column names in extension 1
 *	 
 *     (1) Spectroscopy Star only mode
 *         First column = wavelength in nanometres
 *         Second column = intensity
 *         Third column = error bar
 *     
 *     (2) Polarimetry
 *         1st col = wavelength in nanometres
 *         2d  col = intensity
 *         3rd col = polarisation (Q or U or V or W)
 *         4th col = Check Spectra #1
 *         5th col = Check Spectra #2
 *         6th col = error bar
 *     
 *     (3) Spectroscopy Star + Sky
 *         1st col = wavelength
 *         2d  col = star spectra (sky subtracted)
 *         3rd col = star + sky spectra
 *         4th col = sky spectra
 *         5, 6, 7 = error bars for each column 2, 3, 4
 *
 * The headers spectoscopy star-only look like this:
 *
 * REDUCTIO= 'Intensity'          / Type of reduction                              
 * NORMAL  =                    2 / Normalized and Un-normalized Data              
 * COMMENT  File contains automatic wavelength correction and uncorrected data.    
 * COL1    = 'Wavelength'         / Normalized                                     
 * COL2    = 'Intensity'          / Normalized                                     
 * COL3    = 'ErrorBar'           / Normalized                                     
 * COL4    = 'Wavelength'         / UnNormalized                                   
 * COL5    = 'Intensity'          / UnNormalized                                   
 * COL6    = 'ErrorBar'           / UnNormalized                                   
 * COL7    = 'Wavelength'         / Normalized, no autowave correction             
 * COL8    = 'Intensity'          / Normalized, no autowave correction             
 * COL9    = 'ErrorBar'           / Normalized, no autowave correction             
 * COL10   = 'Wavelength'         / UnNormalized, no autowave correction           
 * COL11   = 'Intensity'          / UnNormalized, no autowave correction           
 * COL12   = 'ErrorBar'           / UnNormalized, no autowave correction  
 *
 * or, for polar:
 *
 * REDUCTIO= 'Polar   '           / Type of reduction                              
 * NORMAL  =                    2 / Normalized and Un-normalized Data              
 * COMMENT  File contains automatic wavelength correction and uncorrected data.    
 * COL1    = 'Wavelength'         / Normalized                                     
 * COL2    = 'Intensity'          / Normalized                                     
 * COL3    = 'Stokes  '           / Normalized                                     
 * COL4    = 'CheckN1 '           / Normalized                                     
 * COL5    = 'CheckN2 '           / Normalized                                     
 * COL6    = 'ErrorBar'           / Normalized                                     
 * COL7    = 'Wavelength'         / UnNormalized                                   
 * COL8    = 'Intensity'          / UnNormalized                                   
 * COL9    = 'Stokes  '           / UnNormalized                                   
 * COL10   = 'CheckN1 '           / UnNormalized                                   
 * COL11   = 'CheckN2 '           / UnNormalized                                   
 * COL12   = 'ErrorBar'           / UnNormalized                                   
 * COL13   = 'Distance'           / Normalized, no autowave correction             
 * COL14   = 'Intensity'          / Normalized, no autowave correction             
 * COL15   = 'Stokes  '           / Normalized, no autowave correction             
 * COL16   = 'CheckN1 '           / Normalized, no autowave correction             
 * COL17   = 'CheckN2 '           / Normalized, no autowave correction             
 * COL18   = 'ErrorBar'           / Normalized, no autowave correction             
 * COL19   = 'Distance'           / UnNormalized, no autowave correction           
 * COL20   = 'Intensity'          / UnNormalized, no autowave correction           
 * COL21   = 'Stokes  '           / UnNormalized, no autowave correction           
 * COL22   = 'CheckN1 '           / UnNormalized, no autowave correction           
 * COL23   = 'CheckN2 '           / UnNormalized, no autowave correction           
 * COL24   = 'ErrorBar'           / UnNormalized, no autowave correction           
 * COMMENT For Stokes Q, V, and W, keep the Stokes parameter sign as is            
 * COMMENT For Stokes U, invert the sign of the Stokes parameter 
 *
 * or for spectoscopy, star+sky"
 * COL1    = 'Wavelength'         / Normalized                                     
 * COL2    = 'Star    '           / Normalized                                     
 * COL3    = 'Star+sky'           / Normalized                                     
 * COL4    = 'Sky     '           / Normalized                                     
 * COL5    = 'ErrorBar1'          / Normalized                                     
 * COL6    = 'ErrorBar2'          / Normalized                                     
 * COL7    = 'ErrorBar3'          / Normalized                                     
 * COL8    = 'Wavelength'         / UnNormalized                                   
 * COL9    = 'Star    '           / UnNormalized                                   
 * COL10   = 'Star+sky'           / UnNormalized                                   
 * COL11   = 'Sky     '           / UnNormalized                                   
 * COL12   = 'ErrorBar1'          / UnNormalized                                   
 * COL13   = 'ErrorBar2'          / UnNormalized                                   
 * COL14   = 'ErrorBar3'          / UnNormalized                                   
 * COL15   = 'Distance'           / Normalized, no autowave correction             
 * COL16   = 'Star    '           / Normalized, no autowave correction             
 * COL17   = 'Star+sky'           / Normalized, no autowave correction             
 * COL18   = 'Sky     '           / Normalized, no autowave correction             
 * COL19   = 'ErrorBar1'          / Normalized, no autowave correction             
 * COL20   = 'ErrorBar2'          / Normalized, no autowave correction             
 * COL21   = 'ErrorBar3'          / Normalized, no autowave correction             
 * COL22   = 'Distance'           / UnNormalized, no autowave correction           
 * COL23   = 'Star    '           / UnNormalized, no autowave correction           
 * COL24   = 'Star+sky'           / UnNormalized, no autowave correction           
 * COL25   = 'Sky     '           / UnNormalized, no autowave correction           
 * COL26   = 'ErrorBar1'          / UnNormalized, no autowave correction           
 * COL27   = 'ErrorBar2'          / UnNormalized, no autowave correction           
 * COL28   = 'ErrorBar3'          / UnNormalized, no autowave correction           
 * 
 * 
 * \endverbatim
 */

const string Wavelength("Wavelength");
const string Intensity("Intensity");
const string Polar("Polar");
const string Star("Star");
const string Starplussky("Star+sky");
const string Sky("Sky");
const string ErrorBar("ErrorBar");
const string ErrorBar1("ErrorBar1");
const string ErrorBar2("ErrorBar2");
const string ErrorBar3("ErrorBar3");
const string Stokes("Stokes");
const string CheckN1("CheckN1");
const string CheckN2("CheckN2");

const string Normalized("Normalized");
const string UnNormalized("UnNormalized");
const string NormalizedNoAuto("Normalized, no autowave correction");
const string UnNormalizedNoAuto("UnNormalized, no autowave correction");
const string COLS[MAXPRODUCTCOLS] = {
	"COL1",
	"COL2",
	"COL3",
	"COL4",
	"COL5",
	"COL6",
	"COL7",
	"COL8",
	"COL9",
	"COL10",
	"COL11",
	"COL12",
	"COL13",
	"COL14",
	"COL15",
	"COL16",
	"COL17",
	"COL18",
	"COL19",
	"COL20",
	"COL21",
	"COL22",
	"COL23",
	"COL24",
	"COL25",
	"COL26",
	"COL27",
	"COL28"
};

/* 
 * \class operaFITSProduct()
 * \brief Basic operaFITSProduct constructor.
 * \extends operaFITSImage
 * \return none
 */
operaFITSProduct::operaFITSProduct(void) : operaFITSImage(),
instrumentmode(MODE_UNKNOWN), columns(0), datapoints(0), next_row(0)
{
}

/* 
 * \class operaFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector)
 * \brief operaFITSProduct constructconsturcted from a spectral order vector.
 * \extends operaFITSImage
 * \oaram Filename
 * \oaram spectralOrerVector
 * \return none
 */
operaFITSProduct::operaFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector) :
instrumentmode(MODE_UNKNOWN), columns(0), datapoints(0), next_row(0)
{
}
/* 
 * \class operaFITSProduct(operaMEFFITSProduct &mefproduct)
 * \brief operaFITSProduct constructed from a spectral order vector.
 * \extends operaFITSImage
 * \oaram Filename
 * \oaram spectralOrerVector
 * \return none
 */
operaFITSProduct::operaFITSProduct(operaMEFFITSProduct &mefproduct)
{
	filename = mefproduct.filename;
	hdu = mefproduct.hdu;	// this is wrong?
	naxis = mefproduct.naxis;
	naxis1 = mefproduct.naxis1;
	naxis2 = mefproduct.naxis2;
	naxis3 = 1;
	naxes[0] = naxis1;
	naxes[1] = naxis2;
	naxes[2] = naxis3;
	extensions = mefproduct.extensions;
	npixels = naxis1 * naxis2;
	npixels_per_slice = npixels_per_slice;
	npixels_per_extension = npixels_per_slice;
	current_extension = 1;
	current_slice = 1;
    mode = mefproduct.mode;
    imageType = mefproduct.imageType;
	
	fptr = mefproduct.fptr;
	bzero = mefproduct.bzero;
	bscale = mefproduct.bscale;
	datatype = mefproduct.datatype;
	bitpix = mefproduct.bitpix;
	isLazy = mefproduct.isLazy;
    isClone = true;				// i.e.do not close fptr
	viewOnly = true;
	AllExtensions = mefproduct.AllExtensions;
	if (AllExtensions) {
		current_extension = 1;
	}
	AllSlices = mefproduct.AllSlices;
	if (AllSlices) {
		current_slice = 1;
	}
	
	// The memcpy may not be needed
	if (mefproduct.pixptr != NULL) {
		pixptr = mefproduct.pixptr;
	} else {
		size_t size = getsize();
        pixptr = malloc(MAX(size, toSize(float_img, npixels))); 
		if (mefproduct.getpixels()) {
			memcpy(pixptr, (void *)((float *)mefproduct.pixptr+(mefproduct.current_extension-1)*naxis1*naxis2*naxis3+(mefproduct.current_slice-1)*naxis1*naxis2), size);
		}
		if (mefproduct.pixptr == NULL) {
			mefproduct.pixptr = pixptr;
		}
	}
}
/* 
 * \class operaFITSProduct()
 * \brief Basic unnames operaFITSProduct constructor with a size.
 * \extends operaFITSImage
 * \oaram Columns
 * \oaram Rows
 * \return none
 */
operaFITSProduct::operaFITSProduct(unsigned Columns, unsigned Rows) : 
operaFITSImage(Columns, Rows, tfloat),
instrumentmode(MODE_UNKNOWN), columns(0), datapoints(0), next_row(0)
{
	columns = Columns;
	datapoints = Rows;
}

/* 
 * \class operaFITSProduct(string Filename, int mode=READWRITE|READONLY)
 * \brief Constructor for readng a FITS file and creating the corresponding object.
 * \param Filename containing a product (i.fits or p.fits)
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaFITSProduct::operaFITSProduct(string Filename, int mode, unsigned Compression) : 
operaFITSImage(Filename, tfloat, mode, Compression, false),
instrumentmode(MODE_UNKNOWN), columns(0), datapoints(0), next_row(0)
{
	if (getnaxis1() > getnaxis2()) {
		rotate90();
	}
	columns = getnaxis1();
	datapoints = getnaxis2();

	try {
		for (unsigned i=0; i<columns; i++) {
			colnames[i] = operaFITSGetHeaderValue(COLS[i]);
		}
	} catch (operaException e) {
		// do nothing, sometimes the column names are not available
	}
	switch (columns) {
		case MODE_STAR_ONLY_COLS:
			instrumentmode = MODE_STAR_ONLY;
			break;
		case MODE_POLAR_COLS:
			instrumentmode = MODE_POLAR;
			break;
		case MODE_STAR_PLUS_SKY_COLS:
			instrumentmode = MODE_STAR_PLUS_SKY;
			break;
		default:
			break;
	}
}

/* 
 * \class operaFITSProduct(string Filename, string baseOnFilename,  int mode, unsigned Columns, unsigned Rows, unsigned Compression)
 * \brief Constructor for readng a FITS object file and creating the corresponding product.
 * \param Filename the product file to create
 * \param baseOnFilename the object file from which to get the headers
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaFITSProduct::operaFITSProduct(string Filename, string baseOnFilename, instrumentmode_t Instrumentmode, unsigned Columns, unsigned Rows, unsigned Compression) : 
operaFITSImage(Filename, Columns, Rows, tfloat, Compression, false),
instrumentmode(MODE_UNKNOWN), columns(0), datapoints(0), next_row(0)
{
	operaFITSImage basedon(baseOnFilename, tfloat, READONLY, cNone, true);
	operaFITSImageCopyHeader(&basedon);
	
	columns = Columns;
	datapoints = Rows;
	instrumentmode = Instrumentmode;
	updateColumnNames();
	setHeaderColumnNames();
}
/* 
 * \class operaFITSProduct(string Filename, unsigned Columns, unsigned Rows, int mode, unsigned Compression)
 * \brief Constructor for readng a FITS object file and creating the corresponding product.
 * \param Filename the product file to create
 * \param baseOnFilename the object file from which to get the headers
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaFITSProduct::operaFITSProduct(string Filename, instrumentmode_t Instrumentmode, unsigned Columns, unsigned Rows, unsigned Compression) : 
operaFITSImage(Filename, Columns, Rows, tfloat, Compression, false),
instrumentmode(MODE_UNKNOWN), columns(0), datapoints(0), next_row(0)
{	
	columns = Columns;
	datapoints = Rows;
	instrumentmode = Instrumentmode;
	updateColumnNames();
	setHeaderColumnNames();
}
/* 
 * \class operaFITSProduct(string Filename, instrumentmode_t, Instrumentmode, int mode=READWRITE|READONLY);
 * \brief Constructor for creating a FITS product file with the correct column width.
 * \param Filename
 * \param mode
 * \param Instrumentmode
 * \throws operaException operaErrorHeaderProblem
 * \return none
 */
operaFITSProduct::operaFITSProduct(string Filename, instrumentmode_t Instrumentmode, unsigned Compression) : 
operaFITSImage(Filename, tfloat, READWRITE, Compression, false),
instrumentmode(MODE_UNKNOWN), columns(0), datapoints(0), next_row(0)
{
	instrumentmode = Instrumentmode;
	updateColumnNames();
	setHeaderColumnNames();
	resize(columns, 0);
}
/* 
 * \class operaFITSProduct()
 * \brief destructor
 * \return void
 */
operaFITSProduct::~operaFITSProduct(void) {
}

/*
 * Helper functions
 */
/*
 * Set the default column names
 */
void operaFITSProduct::updateColumnNames(void) {
	unsigned col = 0;
	switch (instrumentmode) {
		case MODE_STAR_ONLY:
			// Normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = ErrorBar;
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = ErrorBar;
			// un Normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = ErrorBar;
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = ErrorBar;
			columns = col;
			break;
		case MODE_STAR_PLUS_SKY:
			// Normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Star;
			colnames[col++] = Starplussky;
			colnames[col++] = Sky;
			colnames[col++] = ErrorBar2;
			colnames[col++] = ErrorBar3;
			
			// un Normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Star;
			colnames[col++] = Starplussky;
			colnames[col++] = Sky;
			colnames[col++] = ErrorBar1;
			colnames[col++] = ErrorBar2;
			colnames[col++] = ErrorBar3;
			
			// Not Wavelength corrected, normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Star;
			colnames[col++] = Starplussky;
			colnames[col++] = Sky;
			colnames[col++] = ErrorBar1;
			colnames[col++] = ErrorBar2;
			colnames[col++] = ErrorBar3;
			
			// Not Wavelength corrected, un normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Star;
			colnames[col++] = Starplussky;
			colnames[col++] = Sky;
			colnames[col++] = ErrorBar1;
			colnames[col++] = ErrorBar2;
			colnames[col++] = ErrorBar3;
			columns = col;
			break;
		case MODE_POLAR:
			// Normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = Stokes;
			colnames[col++] = CheckN1;
			colnames[col++] = CheckN2;
			colnames[col++] = ErrorBar;
			
			// un Normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = Wavelength;
			colnames[col++] = CheckN1;
			colnames[col++] = CheckN2;
			colnames[col++] = ErrorBar;
			
			// Not Wavelength corrected, normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = Stokes;
			colnames[col++] = CheckN1;
			colnames[col++] = CheckN2;
			colnames[col++] = ErrorBar;
			
			// Not Wavelength corrected, un normalized
			colnames[col++] = Wavelength;
			colnames[col++] = Intensity;
			colnames[col++] = Stokes;
			colnames[col++] = CheckN1;
			colnames[col++] = CheckN2;
			colnames[col++] = ErrorBar;
			columns = col;
			break;
		case MODE_UNKNOWN:
			break;
		default:
			break;
	}
}

/*
 * Set the column names into the fits header
 */
void operaFITSProduct::setHeaderColumnNames(void) {
	unsigned col = 0;
	switch (instrumentmode) {
		case MODE_STAR_ONLY:
			operaFITSSetHeaderValue("REDUCTIO", Intensity, "Type of reduction");                               
			operaFITSSetHeaderValue("NORMAL", "2", "Normalized and Un-normalized Data");  
			operaFITSAddComment("File contains automatic wavelength correction and uncorrected data.");
			operaFITSSetHeaderValue(COLS[col++], Wavelength, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], Intensity, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], Wavelength, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], Intensity, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, UnNormalized);                               
			
			operaFITSSetHeaderValue(COLS[col++], Wavelength, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Intensity, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Wavelength, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Intensity, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, UnNormalizedNoAuto);                               
			columns = col;
			break;
		case MODE_STAR_PLUS_SKY:
			operaFITSSetHeaderValue("REDUCTIO", Intensity, "Type of reduction");                               
			operaFITSSetHeaderValue("NORMAL", "2", "Normalized and Un-normalized Data");                               
			operaFITSAddComment("File contains automatic wavelength correction and uncorrected data.");
			operaFITSSetHeaderValue(COLS[col++], Wavelength, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], Star, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], Starplussky, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], Sky, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar2, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar2, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar3, Normalized);                               
			
			operaFITSSetHeaderValue(COLS[col++], Wavelength, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], Star, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], Starplussky, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], Sky, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar1, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar2, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar3, UnNormalized);                               
			
			operaFITSSetHeaderValue(COLS[col++], Wavelength, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Star, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Starplussky, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Sky, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar1, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar2, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar3, NormalizedNoAuto);                               
			
			operaFITSSetHeaderValue(COLS[col++], Wavelength, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Star, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Starplussky, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Sky, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar2, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar2, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar3, UnNormalizedNoAuto);                               
			columns = col;
			break;
		case MODE_POLAR:
			operaFITSSetHeaderValue("REDUCTIO", Polar, "Type of reduction");                               
			operaFITSSetHeaderValue("NORMAL", "2", "Normalized and Un-normalized Data");                               
			operaFITSAddComment("File contains automatic wavelength correction and uncorrected data.");
			operaFITSAddComment("For Stokes Q, V, and W, keep the Stokes parameter sign as is");
			operaFITSAddComment("For Stokes U, invert the sign of the Stokes parameter");
			operaFITSSetHeaderValue(COLS[col++], Wavelength, Normalized); 
			operaFITSSetHeaderValue(COLS[col++], Intensity, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], Stokes, Normalized);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN1, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN2, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, UnNormalized);                               
			
			operaFITSSetHeaderValue(COLS[col++], Wavelength, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], Intensity, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], Stokes, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN1, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN2, UnNormalized);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, UnNormalized);                               
			
			operaFITSSetHeaderValue(COLS[col++], Wavelength, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Intensity, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Stokes, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN1, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN2, NormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, NormalizedNoAuto);                               
			
			operaFITSSetHeaderValue(COLS[col++], Wavelength, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Intensity, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], Stokes, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN1, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], CheckN2, UnNormalizedNoAuto);                               
			operaFITSSetHeaderValue(COLS[col++], ErrorBar, UnNormalizedNoAuto); 
			columns = col;
			break;
		case MODE_UNKNOWN:
			break;
		default:
			break;
	}
}

/*
 * Interface to spectra in text files compatible with Libre-Esprit
 */

void operaFITSProduct::addRow(float *rowdata, unsigned Columns, unsigned startColumn) {
	if (next_row < datapoints) {
		for (unsigned c=startColumn; c<Columns && c <columns; c++) {
			*this[next_row][c] = rowdata[c];
		}
		next_row++;		
	}
}

/* 
 * \class operaFITSProduct()
 * \brief getNextRow - get the current row
 * \return void
 */
unsigned operaFITSProduct::getNextRow(void) {
	return next_row;		
}

/* 
 * \class operaFITSProduct()
 * \brief setNextRow - set the current row
 * \param Row - unsigned
 * \return void
 */
void operaFITSProduct::setNextRow(unsigned Row) {
	if (Row < datapoints)
		next_row = Row;		
}

/* 
 * \class operaFITSProduct()
 * \brief readText - read a .s spectrum
 * \param Filename - string
 * \return void
 */
void operaFITSProduct::readText(string Filename) {
	instrumentmode = MODE_UNKNOWN;
	string title;
	//bool normalized = false;
	//bool wavelengthcorrected = false;
	
	if (Filename.find(".s") == string::npos)
		throw operaException("operaFITSProduct: spectrum must end in .s ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);

	if (Filename.find("in.s") != string::npos) {
		instrumentmode = MODE_STAR_ONLY;	// actually we don't know...
		//normalized = true;
	}
	if (Filename.find("pn.s") != string::npos) {
		instrumentmode = MODE_POLAR;
		//normalized = true;
	}
	if (Filename.find("iu.s") != string::npos) {
		instrumentmode = MODE_STAR_ONLY;	// actually we don't know...
	}
	if (Filename.find("pu.s") != string::npos) {
		instrumentmode = MODE_POLAR;
	}
	if (Filename.find("iuw.s") != string::npos) {
		instrumentmode = MODE_STAR_ONLY;	// actually we don't know...
		//wavelengthcorrected = true;
	}
	if (Filename.find("puw.s") != string::npos) {
		instrumentmode = MODE_POLAR;
		//wavelengthcorrected = true;
	}
	if (Filename.find("inw.s") != string::npos) {
		instrumentmode = MODE_STAR_ONLY;	// actually we don't know...
		//normalized = true;
		//wavelengthcorrected = true;
	}
	if (Filename.find("pnw.s") != string::npos) {
		instrumentmode = MODE_POLAR;
		//normalized = true;
		//wavelengthcorrected = true;
	}
	/*
	 *     (1) Spectroscopy Star only mode
	 *         First column = wavelength in nanometres
	 *         Second column = intensity
	 *         Third column = error bar )order in the cae of opera raw spactrum)
	 *     
	 *     (2) Polarimetry
	 *         1st col = wavelength in nanometres
	 *         2d  col = intensity
	 *         3rd col = polarisation (Q or U or V or W)
	 *         4th col = Check Spectra #1
	 *         5th col = Check Spectra #2
	 *         6th col = error bar
	 *     
	 *     (3) Spectroscopy Star + Sky
	 *         1st col = wavelength
	 *         2d  col = star spectra (sky subtracted)
	 *         3rd col = star + sky spectra
	 *         4th col = sky spectra
	 *         5, 6, 7 = error bars for each column 2, 3, 4
	 */
	
	ifstream fspectrum(Filename.c_str());	
	
	if (fspectrum.is_open()) {
		int line = 0;
		string dataline;
		float data[10];
		while (fspectrum.good()) {
			getline (fspectrum, dataline);
			if (strlen(dataline.c_str())){
				if (line == 0) {
					title = dataline; 			 
				} else if (line == 1) {
					sscanf(dataline.c_str(), "%u %u", &datapoints, &columns);
					columns++;
					naxis1 = columns;
					naxis2 = datapoints;
					npixels = naxis1*naxis2;
					naxes[0] = naxis1;
					naxes[1] = naxis2;
					compression = 0;
					hdu = 1;
					
					datatype = tfloat;
					bitpix = tobitpix(datatype);
					long size = toSize(bitpix, npixels);	
					pixptr = malloc(size); 
					memset(pixptr, 0, size);
					//resize(columns, datapoints);
				} else {
					if (columns == 3) {	// upena Spectroscopy Star only mode and opera raw spectrum
						sscanf(dataline.c_str(), "%f %f %f", &data[0],  &data[1],  &data[2]);
						for (unsigned i=0; i<columns; i++) {
							*this[i][next_row] = data[i];
						}
					} else if (columns == 6) {	// upena Polarimetry
						sscanf(dataline.c_str(), "%f %f %f %f %f %f", &data[0],  &data[1],  &data[2],  &data[3],  &data[4], &data[5]);
						for (unsigned i=0; i<columns; i++) {
							*this[i][next_row] = data[i];
						}
					} else if (columns == 7) {	// Spectroscopy Star + Sky
						sscanf(dataline.c_str(), "%f %f %f %f %f %f %f", &data[0],  &data[1],  &data[2],  &data[3],  &data[4], &data[5], &data[6]);
						for (unsigned i=0; i<columns; i++) {
							*this[i][next_row] = data[i];
						}
					}
					next_row++;
				}
			}
			line++;
		}
		fspectrum.close();
	}
}

/* 
 * \class operaFITSProduct()
 * \brief writeText - write a .s spectrum
 * \param Filename - string
 * \return void
 */
void operaFITSProduct::writeText(string Filename) {
	ofstream sout(Filename.c_str());
	// write the header for Libre-Esprit compatibility
	sout << "*** Reduced Spectrum of " << (strrchr(Filename.c_str(), '/')?strrchr(Filename.c_str()+1, '/'):Filename) << endl;		
	sout << datapoints << " " << columns << endl;
	for (unsigned line=0; line<datapoints; line++) {
		for (unsigned col=0; col<columns; col++) {
			sout << setprecision(4) << fixed << *this[col][line];
		}
		sout << endl;
	}
	sout.close();
}

/*
 * getters / setters
 */

/* 
 * string getcolname(unsigned col) 
 * \brief returns the string column name of col.
 * \param col
 * \return string
 */
string operaFITSProduct::getcolname(unsigned col) {
	return colnames[col];
}

/* 
 * unsigned getcols() 
 * \brief returns number of columns.
 * \return string
 */
unsigned operaFITSProduct::getcols() {
	return columns;
}

/* 
 * unsigned getrows() 
 * \brief returns number of rows.
 * \return string
 */
unsigned operaFITSProduct::getrows() {
	return datapoints;
}

/* 
 * instrumentmode_t getmode() 
 * \brief returns intrument mode.
 * \return string
 */
instrumentmode_t operaFITSProduct::getmode() {
	return instrumentmode;
}

/* 
 * instrumentmode_t getmode() 
 * \brief returns intrument mode.
 * \return string
 */
void operaFITSProduct::setmode(instrumentmode_t Mode) {
	instrumentmode = Mode;
}


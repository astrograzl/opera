/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaPixelSensitivityMap
 Version: 1.0
 Description: Create a normalized map of pixel-by-pixel sensitivity variations  
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
 Contact: opera@cfht.hawaii.edu
 
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

#include <stdio.h>
#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaPixelSensitivityMap.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaEspadonsImage.h"			// for imtype_t
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaInstrumentProfile.h"		// for operaInstrumentProfile

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"					// for MAXORDERS
#include "libraries/operaFit.h"	
#include "libraries/operaFFT.h"	

/*! \file operaPixelSensitivityMap.cpp */

using namespace std;

/*! 
 * operaPixelSensitivityMap
 * \author Eder Martioli
 * \brief Create a normalized map of pixel-by-pixel sensitivity variations.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;
	string outputfilename;
	string masterbias; 
	string masterflat; 
	string badpixelmask;	
    eCompression compression = cNone;

	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"outputfilename",              1, NULL, 'o'},
		{"masterbias",                  1, NULL, 'b'},
		{"masterflat",                  1, NULL, 'f'},
		{"badpixelmask",                1, NULL, 'm'},
        {"compressiontype",             1, NULL, 'C'},
		{"plot",		optional_argument, NULL, 'p'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",              no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "o:b:f:m:C:v::d::t::p::h", longopts, NULL))  != -1)
    {
		switch(opt) 
		{
			case 'o':
				outputfilename = optarg;
				break;   
			case 'b':		// output
				masterbias = optarg;
				break;
			case 'f':		// masterflat
				masterflat = optarg;
				break;
			case 'm':		// badpixelmask
				badpixelmask = optarg;
				break;
            case 'C':
                compression = (eCompression)atoi(optarg);
                break;
			case 'v':
				verbose = 1;
				break;
			case 'p':
				plot = 1;
				break;
			case 'd':
				debug = 1;
				break;
			case 't':
				trace = 1;
				break;         
			case 'h':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
			case '?':
				printUsageSyntax(argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}	
	
	/*Start the module here*/
	
	try {
		// we need a masterflat...
		if (masterflat.empty()) {
			throw operaException("operaPixelSensitivityMap: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}			
        // we need an output...
		if (outputfilename.empty()) {
			throw operaException("operaPixelSensitivityMap: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		
		if (verbose) {
			cout << "operaPixelSensitivityMap: outputfilename = " << outputfilename << endl; 
			cout << "operaPixelSensitivityMap: masterflat = " << masterflat << endl; 	
			cout << "operaPixelSensitivityMap: masterbias = " << masterbias << endl; 		
			cout << "operaPixelSensitivityMap: badpixelmask = " << badpixelmask << endl;
            cout << "operaPixelSensitivityMap: compression = " << compression << endl;
		}
        
		operaFITSImage *flat = new operaFITSImage(masterflat, tfloat, READONLY);
		
        operaFITSImage *bias = NULL;
        operaFITSImage *badpix = NULL;
        
		if (!masterbias.empty()){
			bias = new operaFITSImage(masterbias, tfloat, READONLY);
		} else {
            bias = new operaFITSImage(flat->getnaxis1(),flat->getnaxis2(),tfloat);
            *bias = 0.0;
        }
        
		if (!badpixelmask.empty()){
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat->getnaxis1(),flat->getnaxis2(),tfloat);
            *badpix = 1.0;
        }

		*flat -= *bias;			// remove bias from masterflat
		
        long npixels = flat->getnpixels();
        
        float *flatData = (float *)flat->getpixels();
        float *badpixData = (float *)badpix->getpixels();
        
        float maxvalue = -BIG;
        
        for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
            if(maxvalue < flatData[pixIndex] && badpixData[pixIndex]) {
                maxvalue = flatData[pixIndex];
            }
        }
        
        operaFITSImage outputImage(outputfilename, flat->getnaxis1(), flat->getnaxis2(), tfloat, compression);
		outputImage.operaFITSImageCopyHeader(flat);
        
		for (unsigned y=0; y<flat->getnaxis2(); y++) {
			for (unsigned x=0; x<flat->getnaxis1(); x++) {
				outputImage[y][x] = ( (*flat)[y][x] + (*bias)[y][x] ) / maxvalue;
			}
		}
        
        outputImage.operaFITSImageSave();
		outputImage.operaFITSImageClose();

        flat->operaFITSImageClose();
		        
        if(bias)
            delete bias;
        if(badpix)
            delete badpix;
        if(flat)
            delete flat;
	}
	catch (operaException e) {
		cerr << "operaPixelSensitivityMap: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPixelSensitivityMap: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	" --inputGeometryFile=<GEOM_FILE>"
	" --outputNormFlat=<FITS_IMAGE>"
	" --inputMasterBias=<FITS_IMAGE>"	
	" --inputMasterFlat=<FITS_IMAGE>"	
	" --badpixelmask=<FITS_IMAGE>"
	" --aperture=<FLT_VALUE>"	
	" --binsize=<UNS_VALUE>  \n\n"
	
	" Example: "+string(modulename)+" -g OLAPAa_sp2_Normal.geom -A 30 -B 20 -f masterflat_OLAPAa_sp2_Normal.fits -b masterbias_OLAPAa_sp2_Normal.fits -o normflat_OLAPAa_sp2_Normal.fits -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	" -g, --inputGeometryFile=<GEOM_FILE>, Input geometry file\n"
	" -o, --outputNormFlat=<FITS_IMAGE>, Ouput pixel-by-pixel sensitivity map FITS image\n"
	" -b, --inputMasterBias=<FITS_IMAGE>, Input Master Bias FITS image\n"	
	" -f, --inputMasterFlat=<FITS_IMAGE>, Input Master Flat-Field FITS image\n"		
	" -m, --badpixelmask=<FITS_IMAGE>, FITS image with badpixel mask\n"
	" -A, --aperture=<FLT_VALUE>, Aperture for extraction in X-direction\n"	
	" -B, --binsize=<UNS_VALUE>, Bin size in Y-direction \n\n";
}


/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterFlat
 Version: 1.0
 Description: THis module creates a Master Flat from a series of flats.
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

#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaMasterFlat.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLib.h"     // for itos

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"

/*! \file operaMasterFlat.cpp */

using namespace std;

/*!
 * operaMasterFlat
 * \author Doug Teeple
 * \brief Creates a master Flat FITS image from a list of input Flat FITS file names.
 * \arg argc
 * \arg argv
 * \note [--images=...]* --output=...[ --pick=\<posint\>0\> ]
 * \note Pick one or median stack images.
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \ingroup core
 * \return EXIT_STATUS
 */
int main(int argc, char *argv[])
{
	int opt;
	string badpixelmask;
	unsigned pick = 0;
	string images[MAXIMAGES];
	unsigned short *flats[MAXIMAGES];
	string listofimages;
	string output;
	string version = "OPERA-1.0";
	string date = "";
	unsigned imageIndex = 0;
	eCompression compression = cNone;
	
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"images",			1, NULL, 'i'},	// series of input flats
		{"list",			1, NULL, 'l'},	// list of input flats		
		{"output",			1, NULL, 'o'},	// a single master flat output fits file
		{"badpixelmask",	1, NULL, 'm'},	// bad pixel mask fits file to ignore pixels
		{"pick",			1, NULL, 'k'},	// don't average, just pick this image index
		{"compressiontype", 1, NULL, 'C'},			
		{"version",			1, NULL, 'V'},
		{"date",			1, NULL, 'a'},
		
		{"plot",			optional_argument, NULL, 'p'},       
		{"verbose",			optional_argument, NULL, 'v'},
		{"debug",			optional_argument, NULL, 'd'},
		{"trace",			optional_argument, NULL, 't'},
		{"help",			no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:l:o:m:k:C:v::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'i':		// images
					images[imageIndex++] = optarg;
					break;
				case 'l':		// list of images
					listofimages = optarg;
					break;					
				case 'o':		// output
					output = optarg;
					break;
				case 'm':		// badpixelmask
					badpixelmask = optarg;
					break;            
				case 'k':		// just pick one - one based not zero based
					pick = atoi(optarg);
					break;    
				case 'C':
					compression = (eCompression)atoi(optarg);	
					break;    
				case 'V':
					version = optarg;	
					break;    
				case 'a':
					date = optarg;	
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
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
				default:
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
			}	// switch
		}	// while
		
		/* 
		 * Read image path from input file and append to list of input images	
		 */ 
		if (!listofimages.empty()) {
			ifstream flist(listofimages.c_str());		
			if (flist.is_open())
			{
				while (flist.good()) {
					getline (flist,images[imageIndex++]);
					if (images[imageIndex-1].size() == 0 || images[imageIndex-1][0] == '#')
						imageIndex--;					
				}	
				flist.close();
			}
		}
		if (debug) {
			for(unsigned i=0;i<imageIndex;i++)
				cout << images[i] << " " << i << endl;
		}
		/*
		 * end of reading list of images
		 */
		if (imageIndex == 0) {
			throw operaException("operaMasterFlat: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (output.empty()) {
			throw operaException("operaMasterFlat: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// if there are not enough images to median combine, just pick one
		if ((imageIndex < 2 && pick == 0) || (imageIndex < pick)) {
			pick = 1;	
			if (verbose)
				cout << "operaMasterFlat: too few images (" << imageIndex << "), picking " << pick << endl;
		}
		if (verbose) {
			cout << "operaMasterFlat: output= " << output << endl;
			cout << "operaMasterFlat: imageIndex= " << imageIndex << endl;
			cout << "operaMasterFlat: pick= " << pick << endl;
			cout << "operaMasterFlat: OPERA version= " << version << endl;
			cout << "operaMasterFlat: Reduction date= " << date << endl;
		}
		/*
		 * choose one input flat -- one-based
		 */
		if (pick) {
			if (pick > imageIndex) {
				throw (operaErrorPickOutofRange);
			}
			if (verbose)
				cout << "operaMasterFlat: picking " << images[pick-1] << endl;
			operaFITSImage *flatIn = new operaFITSImage(images[pick-1], READONLY);
			operaFITSImage *masterFlat = new operaFITSImage(output, flatIn->getnaxis1(), flatIn->getnaxis2(), tushort, compression);
			masterFlat->operaFITSImageCopyHeader(flatIn);
			*masterFlat = *flatIn;	// copy the pixels
			masterFlat->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
			masterFlat->operaFITSAddComment(version);
			masterFlat->operaFITSAddComment("Picking a single flat "+images[pick-1]);
			masterFlat->operaFITSSetHeaderValue("FILENAME", output, "Filename");
			masterFlat->operaFITSDeleteHeaderKey("EXPNUM");
			masterFlat->operaFITSDeleteHeaderKey("OBSID");
			masterFlat->operaFITSImageSave();
			masterFlat->operaFITSImageClose();
			delete masterFlat;
			delete flatIn;
			return EXIT_SUCCESS;
		} else {
			if (verbose)
				cout << "operaMasterFlat: median of " << imageIndex << " images." << endl;
		}
		/*
		 * else do a median
		 */
		long npixels = 0;
		unsigned i;
		operaFITSImage *masterFlat = NULL;
		unsigned short *masterData = NULL;
		for (i=0; i<imageIndex; i++) {
			operaFITSImage *flatIn = new operaFITSImage(images[i], READONLY);
			if (i == 0) {
				masterFlat = new operaFITSImage(output, flatIn->getnaxis1(), flatIn->getnaxis2(), tushort, compression);
				masterFlat->operaFITSImageCopyHeader(flatIn);
				masterData = (unsigned short *)masterFlat->getpixels();
				npixels = masterFlat->getnpixels();
			}
			flats[i] = flatIn->operaFITSImageClonePixelsUSHORT();
			flatIn->operaFITSImageClose();
			delete flatIn;
		}
		flats[i] = NULL;
		masterData = operaArrayMedianCombineUSHORT(imageIndex, npixels, masterData, flats);
		masterFlat->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
		masterFlat->operaFITSAddComment(version);
		masterFlat->operaFITSAddComment("A median of "+itos(imageIndex)+" images.");
		for (i=0; i<imageIndex; i++) {
			masterFlat->operaFITSAddComment("Using flat image "+images[i]);
		}
		masterFlat->operaFITSSetHeaderValue("FILENAME", output, "Filename");
		masterFlat->operaFITSDeleteHeaderKey("EXPNUM");
		masterFlat->operaFITSDeleteHeaderKey("OBSID");
		masterFlat->operaFITSImageSave();
		masterFlat->operaFITSImageClose();
		delete masterFlat;
	}
	catch (operaException e) {
		cerr << "operaMasterFlat: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMasterFlat: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n"
	" Usage: operaMasterFlat [--images=<flat filename>]+ --output=<master flat file name> [badpixlemask=<bad pixel mask file name>] [--pick=<n>] -[dvth]\n";
}	

/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterBias
 Version: 1.0
 Description: This module creates a master bias.
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

/*
 * Algorithm
 *
 * Pick one or median stack images.
 *
 */

#include <getopt.h>
#include <iostream>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaMasterBias.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLib.h"     // for itos

#define MAXDIRNAMESIZE 1000

#include "libraries/operaStats.h"
#include "libraries/operaImage.h"

/*! \file operaMasterBias.cpp */

using namespace std;

/*! 
 * operaMasterBias
 * \author Doug Teeple
 * \brief Creates a master bias FITS image from a list of input bias FITS file names.
 * \arg argc
 * \arg argv
 * \note operaMasterBias [--images=...]* --output=...[ --pick=\<posint\>0\> ]
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
	unsigned short *biases[MAXIMAGES];
	string listofimages;
	string output;
	string version = "OPERA-1.0";
	string date = "";
	unsigned imageIndex = 0;
	eCompression compression = cNone;
	unsigned rotate = 0;
	
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"images",			1, NULL, 'i'},	// series of input flats
		{"list",			1, NULL, 'l'},	// list of input flats		
		{"output",			1, NULL, 'o'},	// a single master flat output fits file
		{"pick",			1, NULL, 'k'},	// a single master flat output fits file
		{"badpixelmask",	1, NULL, 'm'},	// bad pixel mask fits file to ignore pixels
		{"rotate",			1, NULL, 'r'},	// rotate output by 90 degrees
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
		while ((opt = getopt_long(argc, argv, "i:l:o:k:m:r:C:V:a:v::d::t::p::h", longopts, NULL))  != -1) {
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
				case 'r':		// rotate by 90 degrees
					rotate = atoi(optarg);
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
					if(images[imageIndex-1].size() == 0 || images[imageIndex-1][0] == '#')
						imageIndex--;					
				}	
				flist.close();
			}
		}
		if(verbose) {
			for(unsigned i=0;i<imageIndex;i++)
				cout << images[i] << " " << i << endl;
		}
		/*
		 * end of reading list of images
		 */
		if (imageIndex == 0) {
			throw operaException("operaMasterBias: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (output.empty()) {
			throw operaException("operaMasterBias: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (verbose) {
			cout << "operaMasterBias: output= " << output << endl;
			cout << "operaMasterBias: imageIndex= " << imageIndex << endl;
			cout << "operaMasterBias: pick= " << pick << endl;
			cout << "operaMasterBias: OPERA version= " << version << endl;
			cout << "operaMasterBias: Reduction date= " << date << endl;
		}
		// if there are not enough images to median combine, just pick one
		if ((imageIndex < 3 && pick == 0) || (imageIndex < pick)) {
			pick = 1;	
			if (verbose)
				cout << "operaMasterBias: too few images (" << imageIndex << "), picking " << pick << endl;
		}
		/*
		 * choose one input bias -- one-based
		 */
		if (pick) {
			if (pick > imageIndex) {
				throw (operaErrorPickOutofRange);
			}
			if (verbose)
				cout << "operaMasterBias: picking " << images[pick-1] << endl;
			operaFITSImage *biasIn = new operaFITSImage(images[pick-1], READONLY);
			operaFITSImage *masterBias = new operaFITSImage(output, biasIn->getnaxis1(), biasIn->getnaxis2(), tushort, compression);
			masterBias->operaFITSImageCopyHeader(biasIn);
			*masterBias = *biasIn;	// copy the pixels
			masterBias->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
			masterBias->operaFITSAddComment(version);
			masterBias->operaFITSAddComment("Picking a single bias "+images[pick-1]);
			masterBias->operaFITSSetHeaderValue("FILENAME", output, "Filename");
			masterBias->operaFITSDeleteHeaderKey("EXPNUM");
			masterBias->operaFITSDeleteHeaderKey("OBSID");
			if (rotate) {
				masterBias->rotate90();
			}
			masterBias->operaFITSImageSave();
			masterBias->operaFITSImageClose();
			delete biasIn;
			delete masterBias;
			return EXIT_SUCCESS;
		} else {
			if (verbose)
				cout << "operaMasterBias: median of " << imageIndex << " images." << endl;
		}
		/*
		 * else median combine stack
		 */
		operaFITSImage *masterBias = NULL;
		unsigned short *masterData = NULL;
		long npixels = 0;
		unsigned i;
		for (i=0; i<imageIndex; i++) {
			operaFITSImage *biasIn = new operaFITSImage(images[i], READONLY);
			if (i == 0) {
				masterBias = new operaFITSImage(output, biasIn->getnaxis1(), biasIn->getnaxis2(), tushort, compression);
				masterBias->operaFITSImageCopyHeader(biasIn);
				masterData = (unsigned short *)masterBias->getpixels();
				npixels = masterBias->getnpixels();
			}
			biases[i] = biasIn->operaFITSImageClonePixelsUSHORT();
			biasIn->operaFITSImageClose();
			delete biasIn;
		}
		biases[i] = NULL;
		masterData = operaArrayMedianCombineUSHORT(imageIndex, npixels, masterData, biases);
		masterBias->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
		masterBias->operaFITSAddComment(version);
		masterBias->operaFITSAddComment("A median of "+itos(imageIndex)+" images.");
		for (i=0; i<imageIndex; i++) {
			masterBias->operaFITSAddComment("Using bias image "+images[i]);
		}
		masterBias->operaFITSSetHeaderValue("FILENAME", output, "Filename");
		masterBias->operaFITSDeleteHeaderKey("EXPNUM");
		masterBias->operaFITSDeleteHeaderKey("OBSID");
		if (rotate) {
			masterBias->rotate90();
		}
		masterBias->operaFITSImageSave();
		masterBias->operaFITSImageClose();

		delete masterBias;
	}
	catch (operaException e) {
		cerr << "operaMasterBias: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMasterBias: " << operaStrError(errno) << endl;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n"
	" Usage: operaMasterBias [--images=<flat filename>]+ --output=<master flat file name> [badpixlemask=<bad pixel mask file name>] [--pick=<n>] -[dvth]\n";

}	

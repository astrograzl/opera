/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterFabPerot
 Version: 1.0
 Description: Create a Master Fabry-Perot image.
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
#include "core-espadons/operaMasterFabPerot.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLib.h"     // for itos

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"

/*! \file operaMasterFabPerot.cpp */

using namespace std;

/*!
 * operaMasterFabPerot
 * \author Doug Teeple
 * \brief Creates a master Align FITS image from a list of input Align FITS file names.
 * \arg argc
 * \arg argv
 * \note [--images=...]* --output=... [ --pick=\<posint\>0\> ]
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
	unsigned short *aligns[MAXIMAGES];
	string listofimages;
	string output;
	string version = "OPERA-1.0";
	string date = "";
	unsigned imageIndex = 0;
	eCompression compression = cNone;
	
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"images",1, NULL, 'i'},
		{"output",1, NULL, 'o'},
		{"badpixelmask",1, NULL, 'm'},
		{"pick",1, NULL, 'k'},
		{"compressiontype", 1, NULL, 'C'},			
		{"version",			1, NULL, 'V'},
		{"date",			1, NULL, 'a'},
		
		{"plot",		optional_argument, NULL, 'p'},       
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "i:o:k:m:C:V:a:v::d::t::p::p::h", longopts, NULL))  != -1) {
			switch(opt) 
			{
				case 'i':		// images
					images[imageIndex++] = optarg;
					break;
				case 'o':		// output
					output = optarg;
					break;
				case 'k':		// just pick one
					pick = atoi(optarg);
					break;    
				case 'm':		// badpixelmask
					badpixelmask = optarg;
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
				case '?':
					printUsageSyntax();
					exit(EXIT_SUCCESS);
					break;
			}
		}	
		
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
			throw operaException("operaMasterFabPerot: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (output.empty()) {
			throw operaException("operaMasterFabPerot: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// if there are not enough images to median combine, just pick one
		if ((imageIndex < 3 && pick == 0) || (imageIndex < pick)) {
			pick = 1;	
			if (verbose)
				cout << "operaMasterFabPerot: too few images (" << imageIndex << "), picking " << pick << endl;
		}		
		if (verbose) {
			cout << "operaMasterFabPerot: output= " << output << endl;
			cout << "operaMasterFabPerot: imageIndex= " << imageIndex << endl;
			cout << "operaMasterFabPerot: pick= " << pick << endl;
			cout << "operaMasterFabPerot: OPERA version= " << version << endl;
			cout << "operaMasterFabPerot: Reduction date= " << date << endl;
		}
		/*
		 * choose one input Fabry-Perot -- one-based
		 */
		if (pick) {
			if (pick > imageIndex) {
				throw (operaErrorPickOutofRange);
			}
			if (verbose)
				cout << "operaMasterFabPerot: picking " << pick << endl;
			operaFITSImage *alignIn = new operaFITSImage(images[pick-1], READONLY);
			operaFITSImage *masterAlign = new operaFITSImage(output, alignIn->getnaxis1(), alignIn->getnaxis2(), tushort, compression);
			masterAlign->operaFITSImageCopyHeader(alignIn);
			*masterAlign = *alignIn;	// copy the pixels
			masterAlign->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
			masterAlign->operaFITSAddComment(version);
			masterAlign->operaFITSAddComment("Picking a single align "+images[pick-1]);
			masterAlign->operaFITSSetHeaderValue("FILENAME", output, "Filename");
			masterAlign->operaFITSDeleteHeaderKey("EXPNUM");
			masterAlign->operaFITSDeleteHeaderKey("OBSID");
			masterAlign->operaFITSImageSave();
			masterAlign->operaFITSImageClose();
			delete masterAlign;
			delete alignIn;
			return EXIT_SUCCESS;
		} else {
			if (verbose)
				cout << "operaMasterFabPerot: median of " << imageIndex << " images." << endl;
		}
		/*
		 * else median combine stack
		 */
		operaFITSImage *masterAlign = NULL;
		unsigned short *masterData = NULL;
		long npixels = 0;
		unsigned i;
		for (i=0; i<imageIndex; i++) {
			operaFITSImage *alignIn = new operaFITSImage(images[i], READONLY);
			if (i == 0) {
				masterAlign = new operaFITSImage(output, alignIn->getnaxis1(), alignIn->getnaxis2(), tushort, compression);
				masterAlign->operaFITSImageCopyHeader(alignIn);
				masterData = (unsigned short *)masterAlign->getpixels();
				npixels = masterAlign->getnpixels();
			}
			aligns[i] = alignIn->operaFITSImageClonePixelsUSHORT();
			alignIn->operaFITSImageClose();
			delete alignIn;
		}
		aligns[i] = NULL;
		masterData = operaArrayMedianCombineUSHORT(imageIndex, npixels, masterData, aligns);
		masterAlign->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
		masterAlign->operaFITSAddComment(version);
		masterAlign->operaFITSAddComment("A median of "+itos(imageIndex)+" images.");
		for (i=0; i<imageIndex; i++) {
			masterAlign->operaFITSAddComment("Using align image "+images[i]);
		}
		masterAlign->operaFITSSetHeaderValue("FILENAME", output, "Filename");
		masterAlign->operaFITSDeleteHeaderKey("EXPNUM");
		masterAlign->operaFITSDeleteHeaderKey("OBSID");
		masterAlign->operaFITSImageSave();
		masterAlign->operaFITSImageClose();
		
		delete masterAlign;
	}
	catch (operaException e) {
		cerr << "operaMasterFabPerot: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMasterFabPerot: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	
	cout <<
	"\n"
	" Usage: operaMasterFabPerot [--images=<flat filename>]+ --output=<master flat file name> [badpixlemask=<bad pixel mask file name>] [--pick=<n>] -[dvth]\n";

}	

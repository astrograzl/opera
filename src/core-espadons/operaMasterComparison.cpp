/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterComparison
 Version: 1.0
 Description: Create a master comparison.
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jan/2011
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

#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaMasterComparison.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaLib.h"     // for itos

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"

#ifndef SATURATIONLIMIT
#define SATURATIONLIMIT (unsigned short)65535  // this should be retrieved from the config/param file
#endif

/*! \file operaMasterComparison.cpp */

using namespace std;

/*! 
 * operaMasterComparison
 * \author Doug Teeple
 * \brief Creates a master Comparison FITS image from a list of input Comparison FITS file names.
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
	string listofimages;
	string output;
	string version = "OPERA-1.0";
	string date = "";
	unsigned imageIndex = 0;
	eCompression compression = cNone;
	
	string expTimeFITSKeyword = "EXPTIME";
    string masterbias;
    unsigned combineMethod = 1;
    /*
     * combineMethod 1: median combine.
     *
     * combineMethod 2: stack and save total exposure time per pixel. Note that
     * a given pixel may be saturated in some images and not saturated in others
     * so the only way to bring all pixel counts to the same units is
     * by knowing the total exposure time accounted for each pixel. The other trick
     * here is that the masterData should be initiallized as biasData to avoid negative
     * assignments to the unsigned variable.
     *
     * combineMethod 3: a simple pixel-by-pixel sum ignoring saturated pixels
     */

    unsigned short saturationLimit = (unsigned short)(SATURATIONLIMIT);
    float outputExposureTime = 60;
    
    bool biasConstant = true;
    
    bool truncateOuputFluxToSaturation = true;
    
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"images",				1, NULL, 'i'},
		{"output",				1, NULL, 'o'},
		{"badpixelmask",		1, NULL, 'm'},
		{"masterbias",			1, NULL, 'b'},
		{"pick",				1, NULL, 'k'},
		{"compressiontype",		1, NULL, 'C'},
		{"expTimeFITSKeyword",	1, NULL, 'T'},
		{"combineMethod",		1, NULL, 'M'},
		{"saturationLimit",		1, NULL, 'S'},
		{"outputExposureTime",	1, NULL, 'E'},
		{"biasConstant",		1, NULL, 'B'},
		{"truncateOuputFluxToSaturation", 1, NULL, 's'}, 
		{"version",				1, NULL, 'V'},
		{"date",				1, NULL, 'a'},

		{"plot",				optional_argument, NULL, 'p'},       
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	try {
		while((opt = getopt_long(argc, argv, "i:l:o:k:m:b:C:T:M:S:E:B:V:a:s:v::d::t::p::h", 
								 longopts, NULL))  != -1)
		{
			switch(opt) 
			{
				case 'i':		// images
					images[imageIndex++] = optarg;
					break;
				case 'l':		// list of images
					listofimages = optarg;
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
				case 'b':		// master bias
					masterbias = optarg;
					break;
				case 'C':
					compression = (eCompression)atoi(optarg);
					break;
				case 'T':
					expTimeFITSKeyword = optarg;
					break;
                case 'M':		// method to combine images 1. Median, 2. Stack
                    combineMethod = atoi(optarg);
                    break;
                case 'S':
                    saturationLimit = (unsigned short)atoi(optarg);
                    break;
                case 'E':
                    outputExposureTime = atoi(optarg);
                    break;
                case 'B':
                    biasConstant = (atoi(optarg)?true:false);
                    break;
                case 's':
                    truncateOuputFluxToSaturation = (atoi(optarg)?true:false);
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
					printUsageSyntax(argv[0]);
					exit(EXIT_SUCCESS);
					break;
				case '?':
					printUsageSyntax(argv[0]);
					exit(EXIT_SUCCESS);
					break;
			}
		}	
		
        
		if (verbose) {
            cout << "operaMasterComparison: listofimages = " << listofimages << endl;
            cout << "operaMasterComparison: output = " << output << endl;
            cout << "operaMasterComparison: compression = " << compression << endl;
            cout << "operaMasterComparison: pick = " << pick << endl;
            cout << "operaMasterComparison: badpixelmask = " << badpixelmask << endl;
            cout << "operaMasterComparison: masterbias = " << masterbias << endl;
 			cout << "operaMasterComparison: combineMethod = " << combineMethod << endl;
 			cout << "operaMasterComparison: saturationLimit = " << saturationLimit << endl;
 			cout << "operaMasterComparison: outputExposureTime = " << outputExposureTime << endl;
 			cout << "operaMasterComparison: biasConstant = " << biasConstant << endl;
 			cout << "operaMasterComparison: truncateOuputFluxToSaturation = " << truncateOuputFluxToSaturation << endl;
			cout << "operaMasterComparison: expTimeFITSKeyword = " << expTimeFITSKeyword << endl;
			cout << "operaMasterComparison: OPERA version= " << version << endl;
			cout << "operaMasterComparison: Reduction date= " << date << endl;
        }
        
		// we need a masterbias for methods 2 ... 
		if (combineMethod == 2 || combineMethod==3) {
            if (masterbias.empty()) {
                throw operaException("operaMasterComparison: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
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
        
		if (verbose) {
            for(unsigned i=0;i<imageIndex;i++) {
                cout << "operaMasterComparison: image[" << i << "] = " << images[i] << endl;
            }
		}
		/*
		 * end of reading list of images
		 */
		if (imageIndex == 0) {
			throw operaException("operaMasterComparison: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (output.empty()) {
			throw operaException("operaMasterComparison: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// if there are not enough images to median combine, just pick one
		if ((imageIndex < 3 && pick == 0) || (imageIndex < pick)) {
			pick = 1;	
			if (verbose)
				cout << "operaMasterComparison: too few images (" << imageIndex << "), picking " << pick << endl;
		}		
		/*
		 * choose one input image -- one-based
		 */
		if (pick) {
			if (pick > imageIndex) {
				throw (operaErrorPickOutofRange);
			}
			if (verbose)
				cout << "operaMasterComparison: picking " << images[pick-1] << endl;
			operaFITSImage *compIn = new operaFITSImage(images[pick-1], READONLY);
			operaFITSImage *masterComparison = new operaFITSImage(output, compIn->getnaxis1(), compIn->getnaxis2(), tushort, compression);
			masterComparison->operaFITSImageCopyHeader(compIn);
			*masterComparison = *compIn;	// copy the pixels
			masterComparison->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
			masterComparison->operaFITSAddComment(version);
			masterComparison->operaFITSAddComment("Picking a single comparison "+images[pick-1]);
			masterComparison->operaFITSSetHeaderValue("FILENAME", output, "Filename");
			masterComparison->operaFITSDeleteHeaderKey("EXPNUM");
			masterComparison->operaFITSDeleteHeaderKey("OBSID");
			masterComparison->operaFITSImageSave();
			masterComparison->operaFITSImageClose();
			delete masterComparison;
			delete compIn;
			return EXIT_SUCCESS;
		} else {
			if (verbose) {
				switch (combineMethod) {
					case 1:
						cout << "operaMasterComparison: median of " << imageIndex << " images." << endl;
						break;
					case 2:
						cout << "operaMasterComparison: etime-based combination of " << imageIndex << " images." << endl;
						break;
					case 3:
						cout << "operaMasterComparison: sum ignoring saturated pixels of " << imageIndex << " images." << endl;
						break;
					default:
						break;
				}
			}
		}

        string exposuretimes[MAXIMAGES];
        
        operaFITSImage *masterComparison = NULL;
        operaFITSImage *badpix = NULL;
        operaFITSImage *bias = NULL;
        
        unsigned short *masterData = NULL;
        unsigned short *comparisons[MAXIMAGES];
        
        unsigned short *biasData = NULL;        
        unsigned short *badpixData = NULL;
        float *exptimePerPixel = NULL;
        
        unsigned short biasConstantValue = 0;
        
        long npixels = 0;
        unsigned i;
        for (i=0; i<imageIndex; i++) {
            operaFITSImage *comparisonIn = new operaFITSImage(images[i], READONLY);

            if (i == 0) {
                masterComparison = new operaFITSImage(output, comparisonIn->getnaxis1(), comparisonIn->getnaxis2(), tushort, compression);
                masterComparison->operaFITSImageCopyHeader(comparisonIn);
                masterData = (unsigned short *)masterComparison->getpixels();
                npixels = masterComparison->getnpixels();
                if(verbose)
					cout << "operaMasterComparison: base image =" << images[i]  << " NAXIS1=" << comparisonIn->getnaxis1() << " NAXIS2=" << comparisonIn->getnaxis2()<< " npixels=" << npixels << endl;

                if (!masterbias.empty()) {
                    bias = new operaFITSImage(masterbias, tushort, READONLY);
                } else {
                    bias = new operaFITSImage(comparisonIn->getnaxis1(),comparisonIn->getnaxis2(),tushort);
                    //biasData = new unsigned short[npixels];
                }
                biasData = (unsigned short *)bias->getpixels();
                if(biasConstant) {
                    biasConstantValue = operaArrayMedianQuickUSHORT(bias->getnpixels(),bias->operaFITSImageClonePixelsUSHORT());
                }
                if(debug)
                    cerr << "operaMasterComparison: biasConstantValue=" << biasConstantValue << endl;
                
                if (!badpixelmask.empty()){
                    badpix = new operaFITSImage(badpixelmask, tushort, READONLY);
                } else {
                    badpix = new operaFITSImage(comparisonIn->getnaxis1(),comparisonIn->getnaxis2(),tushort);
                }
                badpixData = (unsigned short *)badpix->getpixels();
                
                if(combineMethod==2 || combineMethod==3) {
                    exptimePerPixel = new float[npixels];
                    for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
                        if(badpixelmask.empty()){
                            badpixData[pixIndex] = 1;
                        }
                        exptimePerPixel[pixIndex] = 0.0;
                        if(biasConstant) {
                            masterData[pixIndex] =  biasConstantValue;
                        } else {
                            masterData[pixIndex] = biasData[pixIndex]; // use biasData if biasConstant is zero
                        }
                    }
                }
            }
            comparisons[i] = comparisonIn->operaFITSImageClonePixelsUSHORT();
            
            if(combineMethod==2 || combineMethod==3) {
                /*
                 * Method 2: stack and save total exposure time per pixel. Note that 
                 * a given pixel may be saturated in some images and not saturated in others
                 * so the only way to bring all pixel counts to the same units is
                 * by knowing the total exposure time accounted for each pixel. The other trick
                 * here is that the masterData should be initiallized as biasData to avoid negative
                 * assignments to the unsigned variable.
                 */
                
                /*
                 * Method 3: This method is a simple pixel-by-pixel stack ignoring saturated pixels
                 */
                
                exposuretimes[i] = comparisonIn->operaFITSGetHeaderValue(expTimeFITSKeyword.c_str());
                float exptime = atof(exposuretimes[i].c_str());
                if(verbose) {
                    cout << "operaMasterComparison: image= " << images[i] << " -> " << expTimeFITSKeyword << "= " << exptime << endl;
                }
                for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
                    if(debug) {
                        cerr << "operaMasterComparison: badpixData[pixIndex]=" << badpixData[pixIndex] << " comparisons[i][pixIndex]=" << comparisons[i][pixIndex] << " saturationLimit=" << saturationLimit << " npixels=" << npixels << " exptimePerPixel[" <<pixIndex<<"]=" << exptimePerPixel[pixIndex] << " exptime=" << exptime << " masterData[pixIndex]=" << masterData[pixIndex] << endl;
                    }
                    
                    if(comparisons[i][pixIndex] < saturationLimit && badpixData[pixIndex] == 1) {
                        masterData[pixIndex] += int(comparisons[i][pixIndex]) - int(biasData[pixIndex]); // WARNING: if *master is not initialized as *bias or biasConstant this could be negative
                        exptimePerPixel[pixIndex] += exptime;
                    }
                }
            }
            comparisonIn->operaFITSImageClose();
            delete comparisonIn;
        }
        comparisons[i] = NULL;

		masterComparison->operaFITSAddComment("Created by the OPERA Open Source Pipeline "+date);
		masterComparison->operaFITSAddComment(version);
        if(combineMethod==1) {
            /*
             * Method 1: median combine stack
             */
            masterData = operaArrayMedianCombineUSHORT(imageIndex, npixels, masterData, comparisons);

            masterComparison->operaFITSAddComment("A median of "+itos(imageIndex)+" images.");
            for (i=0; i<imageIndex; i++) {
                masterComparison->operaFITSAddComment("Using comparison image "+images[i]);
            }
        } else if(combineMethod==2) {

            for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
                if(debug)
                    cerr << "operaMasterComparison: exptimePerPixel[pixIndex]=" << exptimePerPixel[pixIndex] << " masterData[pixIndex]=" << masterData[pixIndex] << endl;
                
                if(exptimePerPixel[pixIndex] > 0) {
                    if(biasConstant) {
                        masterData[pixIndex] = (unsigned)(((float(masterData[pixIndex]) - float(biasConstantValue))
                                                           *outputExposureTime/exptimePerPixel[pixIndex]) + (float)biasConstantValue);
                    } else {
                        masterData[pixIndex] = (unsigned)(((float(masterData[pixIndex]) - float(biasData[pixIndex]))
                                                           *outputExposureTime/exptimePerPixel[pixIndex]) + (float)biasData[pixIndex]);
                    }

                    if(masterData[pixIndex] > saturationLimit && truncateOuputFluxToSaturation) {
                        masterData[pixIndex] = saturationLimit;
                    }
                } else {
                    masterData[pixIndex] = saturationLimit;
                }
            }
            masterComparison->operaFITSAddComment("A stack of "+itos(imageIndex)+" images.");
            masterComparison->operaFITSAddComment("Combined flux is equivalent to an exptime of "+ftos(outputExposureTime)+" s");
            for (i=0; i<imageIndex; i++) {
                masterComparison->operaFITSAddComment("Using comparison image "+images[i]+" exptime="+exposuretimes[i]+" s");
            }
        }  else if(combineMethod==3) {
            masterComparison->operaFITSAddComment("A stack of "+itos(imageIndex)+" images.");
            masterComparison->operaFITSAddComment("Combined flux by the sum of all images");
            
            for(unsigned pixIndex=0;pixIndex<(unsigned)npixels;pixIndex++) {
                if(exptimePerPixel[pixIndex] > 0) {
                    if(masterData[pixIndex] > saturationLimit && truncateOuputFluxToSaturation) {
                        masterData[pixIndex] = saturationLimit;
                    }
                } else {
                    masterData[pixIndex] = saturationLimit;
                }
            }
            for (i=0; i<imageIndex; i++) {
                masterComparison->operaFITSAddComment("Using comparison image "+images[i]+" exptime="+exposuretimes[i]+" s");
            }
        } else {
            throw operaException("operaMasterComparison: combineMethod not supported. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
        
        masterComparison->operaFITSSetHeaderValue("FILENAME", output, "Filename");
        masterComparison->operaFITSDeleteHeaderKey("EXPNUM");
        masterComparison->operaFITSDeleteHeaderKey("OBSID");
        masterComparison->operaFITSImageSave();
        masterComparison->operaFITSImageClose();
        
        delete masterComparison;
        if(badpix)
            delete badpix;
        if(bias)
            delete bias;
	}
	catch (operaException e) {
		cerr << "operaMasterComparison: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMasterComparison: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
    " --images=<FITS_IMAGES>"
    " --output=<FITS_IMAGE>"
	" --badpixelmask=<FITS_IMAGE>"
	" --masterbias=<FITS_IMAGE>"
    " --pick=<BOOL>"
    " --compressiontype=<COMPRESSION_TYPE>"
    " --expTimeFITSKeyword=<STRING>"
    " --combineMethod=<UNS_VALUE>"
    " --saturationLimit=<UNSIGNED SHORT>"
    " --outputExposureTime=<FLOAT>"
    " --biasConstant=<BOOL>"
    " --truncateOuputFluxToSaturation=<BOOL>\n\n"
	" Example: "+string(modulename)+" --output=mastercomparison_Stack_Normal.fits.fz --images=/data/espadons/51Peg-12BQ10-Nov27/1598185c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598186c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598187c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598188c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598189c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598190c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598191c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598192c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598193c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598194c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598195c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598196c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598197c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598198c.fits --images=/data/espadons/51Peg-12BQ10-Nov27/1598249c.fits  --compressiontype=11 --pick=0 --expTimeFITSKeyword=\"EXPTIME\" --combineMethod=2 --masterbias=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Nov27/masterbias_OLAPAa_pol_Normal.fits.fz --saturationLimit=65500 --outputExposureTime=60 --truncateOuputFluxToSaturation=1 -v --biasConstant=1\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
    "  -i, --images=<FITS_IMAGE>, Input individual comparison FITS files\n"
    "  -o, --output=<FITS_IMAGE>, Output master comparison\n"
	"  -m, --badpixelmask=<FITS_IMAGE>, FITS image with badpixel mask\n"
	"  -b, --masterbias=<FITS_IMAGE>, FITS image with masterbias\n"
    "  -k, --pick=<BOOL>, Choose this option to pick only one image\n"
    "  -C, --compressiontype=<COMPRESSION_TYPE>, Ouput compression type\n"
    "  -T, --expTimeFITSKeyword=<STRING>, String to identify exposure time header keyword\n"
    "  -M, --combineMethod=<UNS_VALUE>, Method for combining images\n"
    "                              Available options are = 1, 2, and  3, where: \n"
    "                              1. Median (default)\n"
    "                              2. Mean weighted by exposure time \n"
    "                              3. Sum \n"
    "  -S, --saturationLimit=<UNSIGN>, Define saturation or linearity limit\n"
    "  -E, --outputExposureTime=<FLOAT>, Exposure time to reset output image (only used in combineMethod=2) \n"    
    "  -B, --biasConstant=<BOOL>, Use median bias constant to add to output image. Useful to avoid negative numbers and when bias has low gradient. \n"
    "  -s, --truncateOuputFluxToSaturation=<BOOL>, Set flux value to saturation if final output is above saturation\n\n";
}

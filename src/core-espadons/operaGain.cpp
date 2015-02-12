/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaGain
 Version: 1.0
 Description: Calculate gain and noise.
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
#include <iostream>
#include <iomanip>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaGain.h"
#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaConfigurationAccess.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"
#include "libraries/operaFit.h"	

/*! \file operaGain.cpp */

using namespace std;

/*!
 * operaGain
 * \author Doug Teeple & Eder Martioli
 * \brief Output gain and noise based on at least 2 flat images.
 * \arg argc
 * \arg argv
 * \arg [--biasimgs=...]*
 * \arg [--flatimgs=...]* 
 * \arg --listofbiasimgs=...
 * \arg --listofflatimgs=... 
 * \arg --badpixelmask=...
 * \arg --output=...
 * \arg --defaultgain=...
 * \arg --defaultnoise=...
 * \note --subwindow="x1 nx y1 ny"
 * \note --defaultgain=...
 * \note --defaultnoise=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 */

/**
 * \defgroup core Core Modules
 */

int main(int argc, char *argv[])
{
	int opt;
    
	string biasimgs[MAXIMAGES];			
	string flatimgs[MAXIMAGES];	
	string listofbiasimgs;
	string listofflatimgs;
	string output; 
	string badpixelmask;
    
	int debug=0, verbose=0, trace=0, plot=0;
	
	float noise[MAXNAMPS];
	float gain[MAXNAMPS];
	float gainError[MAXNAMPS];
	float bias[MAXNAMPS];
	
	for(unsigned i=0;i<MAXNAMPS;i++) {
		noise[i] = 0.0;
		gain[i] = 0.0;
		gainError[i] = 0.0;	
		bias[i] = 0.0;	
	}
	
	float defaultgain = 1.0;
	float defaultnoise = 0.0;
	
    DATASEC_t datasec = {1,2068,1,4608};
    DATASEC_t dseca   = {21,1044,1,4608};
    DATASEC_t dsecb   = {1045,2068,1,4608};
    
	struct subwindow {
		int x0, nx;
		int y0, ny;
	} subwindow = {0,0,0,0};
	
	unsigned biasimgIndex = 0;
	unsigned flatimgIndex = 0;
	
	unsigned gainMinPixPerBin = 100;
	unsigned gainMaxNBins = 1000;
	unsigned gainLowestCount = 1000;
	unsigned gainHighestCount = 25000;
    
    unsigned maximages = 12;    // otherwise we run out of memory o 32 bit systems...
	
	unsigned namps = 1;         // 1 for EEV1 and OLAPAa, 2 for OLAPAab
	
	struct option longopts[] = {
		{"biasimgs",            1, NULL, 'b'},
		{"flatimgs",            1, NULL, 'f'},		
		{"listofbiasimgs",      1, NULL, 'L'},				
		{"listofflatimgs",      1, NULL, 'F'},		
		{"output",              1, NULL, 'o'},
		{"badpixelmask",        1, NULL, 'm'},		
		{"maximages",           1, NULL, 's'},
		{"numberofamplifiers",  1, NULL, 'a'},
		{"defaultgain",         1, NULL, 'G'},
		{"defaultnoise",        1, NULL, 'N'},
		{"gainMinPixPerBin",    1, NULL, 'n'},
		{"gainMaxNBins",        1, NULL, 'x'},
		{"gainLowestCount",     1, NULL, 'l'},
		{"gainHighestCount",    1, NULL, 'c'},
		{"subwindow",           1, NULL, 'w'},
		{"DATASEC",             1, NULL, 'D'},
		{"DSECA",               1, NULL, 'A'},
		{"DSECB",               1, NULL, 'B'},
        
		{"plot",                optional_argument, NULL, 'p'},
		{"verbose",             optional_argument, NULL, 'v'},
		{"debug",               optional_argument, NULL, 'd'},
		{"trace",               optional_argument, NULL, 't'},
		{"help",                no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "b:f:L:F:o:m:a:s:G:N:n:x:l:c:w:D:A:B:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'b':		// bias images
				if (biasimgIndex < maximages)
                    biasimgs[biasimgIndex++] = optarg;
				break;
			case 'f':		// flat images
                if (flatimgIndex < maximages)
                    flatimgs[flatimgIndex++] = optarg;
				break;		
			case 'L':		// list of bias images
				listofbiasimgs = optarg;
				break;	
			case 'F':		// list of flat images
				listofflatimgs = optarg;
				break;		
			case 'o':		// output
				output = optarg;
				break;
			case 'm':		// badpixelmask
                badpixelmask = optarg;
				break;
			case 's':		// max images
                maximages = atoi(optarg);
				break;
			case 'a':		// number of amplifiers = 1 or 2 (for Espadons)
				namps = atoi(optarg);
				break;	
			case 'G':		// default gain
				defaultgain = atof(optarg);
				for(unsigned i=0;i<MAXNAMPS;i++) {
					gain[i] = defaultgain;	
				}				
				break;      
			case 'N':		// default noise
				defaultnoise = atof(optarg);
				for(unsigned i=0;i<MAXNAMPS;i++) {
					noise[i] = defaultnoise;	
				}		
				break;               
			case 'n':		// gainMinPixPerBin
				gainMinPixPerBin = atoi(optarg);
				break;            
			case 'x':		// gainMaxNBins
				gainMaxNBins = atoi(optarg);
				break;            
			case 'l':		// gainLowestCount
				gainLowestCount = atoi(optarg);
				break;            
			case 'c':		// gainHighestCount
				gainHighestCount = atoi(optarg);
				break;            
			case 'w':		// subwindow - subwindow within which data will used for the gain/noise calculation
				if (strlen(optarg))
					sscanf(optarg, "%d %d %d %d", &subwindow.x0, &subwindow.nx, &subwindow.y0, &subwindow.ny);
				break;
			case 'D':
				if (strlen(optarg))
					sscanf(optarg, "%u %u %u %u", &datasec.x1, &datasec.x2, &datasec.y1, &datasec.y2);
				break;
			case 'A':
				if (strlen(optarg))
					sscanf(optarg, "%u %u %u %u", &dseca.x1, &dseca.x2, &dseca.y1, &dseca.y2);
				break;
			case 'B':
				if (strlen(optarg))
					sscanf(optarg, "%u %u %u %u", &dsecb.x1, &dsecb.x2, &dsecb.y1, &dsecb.y2);
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
	
	/* 
	 * Read image path from input file and append to list of input images	
	 */ 
	// Read list of bias images
	if (!listofbiasimgs.empty()) {
		ifstream flist(listofbiasimgs.c_str());		
		if (flist.is_open())
		{
			while (flist.good()) {
				getline (flist,biasimgs[biasimgIndex++]);
				if (biasimgs[biasimgIndex-1].size() == 0 || biasimgs[biasimgIndex-1][0] == '#')
					biasimgIndex--;					
                if (biasimgIndex >= maximages)
                    break;
			}
			flist.close();
		}
	}
	// Read list of flat images	
	if (!listofflatimgs.empty()) {
		ifstream flist(listofflatimgs.c_str());		
		if (flist.is_open())
		{
			while (flist.good()) {
				getline (flist,flatimgs[flatimgIndex++]);
				if (flatimgs[flatimgIndex-1].size() == 0 || flatimgs[flatimgIndex-1][0] == '#')
					flatimgIndex--;					
                if (flatimgIndex >= maximages)
                    break;
			}
			flist.close();
		}
	}
	/*
	 * end of reading list of images
	 */
	
	// we need at least 2 biases...
	if (biasimgIndex < 1) {
		throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
	}
	// we need at least 2 flats to calc the gain...
	if (flatimgIndex < 2) {
		throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
	}
	// we need an output...
	if (output.empty()) {
		throw operaException("operaGain: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
	}
	// only accept number of amplifiers = 1 or 2...
	if (namps != 1 && namps != 2) {
		throw operaException("operaGain: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
	}	
	
	if (debug) {
		for(unsigned i=0; i<biasimgIndex; i++)
			cout << "operaGain: input bias: " << i << ' ' << biasimgs[i] << endl;
		for(unsigned i=0; i<flatimgIndex;i ++)
			cout << "operaGain: input flat: " << i << ' ' << flatimgs[i] << endl;		
	}
	
	if (verbose) {
		cout << "operaGain: subwindow = " << subwindow.x0 << " " << subwindow.nx << " "<< subwindow.y0  << " "<< subwindow.ny << endl; 
		cout << "operaGain: gainMinPixPerBin = " << gainMinPixPerBin << " gainMinPixPerBin = " << gainMinPixPerBin << " gainLowestCount = " << gainLowestCount << " gainHighestCount = " << gainHighestCount << endl; 
	}
	
	try {
		/*
		 * open badpixelmask and load data into a vector
		 */
		
		long npixels = subwindow.nx * subwindow.ny;
		
		int x1=0,x2=0,y1=0,y2=0;
		int amp_x0=0,amp_y0=0,amp_xf=0,amp_yf=0;
		
		operaSpectralOrderVector *orders = new operaSpectralOrderVector();
		
        float *flatdata[MAXIMAGES];
        float *biasdata[MAXIMAGES];
        /*
         * open bias images and load data into vectors
         */
        operaFITSImage *badpix = NULL;
        operaFITSImage *biasIn = NULL;
        operaFITSImage *flatIn = NULL;
        
        float *badpixdata = NULL;
        
		for (unsigned amp=1; amp<=namps; amp++) {	// loop over all possible amplifiers
			for (unsigned i=0; i<biasimgIndex; i++) {	
                biasIn = new operaFITSImage(biasimgs[i], tfloat, READONLY);					
				if (i == 0) {
                    if (!badpixelmask.empty()){
                        badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
                    } else {
                        badpix = new operaFITSImage(biasIn->getnaxis1(),biasIn->getnaxis2(),tfloat);
                        *badpix = 1.0;
                    }
                    badpixdata = (float *)badpix->operaFITSImageClonePixels(amp_x0, amp_y0, amp_xf, amp_yf);
					// delete badpix;
                    
                    if (namps == 2) {
                        if (amp == 1) {
                            x1 = dseca.x1;
                            x2 = dseca.x2;
                            y1 = dseca.y1;
                            y2 = dseca.y2;
                            if (verbose) {
                                cout << "operaGain: amp " << amp << " dseca = " << x1 << " : " << x2 << " , "<< y1  << " : "<< y2 << endl;
                            }					
						} else {
                            x1 = dsecb.x1;
                            x2 = dsecb.x2;
                            y1 = dsecb.y1;
                            y2 = dsecb.y2;
							if (verbose) {
								cout << "operaGain: amp " << amp << " dsecb = " << x1 << " : " << x2 << " , "<< y1  << " : "<< y2 << endl;
							}
                        }
                        amp_x0 = x1 + subwindow.x0;
                        amp_y0 = y1 + subwindow.y0;
                        amp_xf = amp_x0 + subwindow.nx;
                        amp_yf = amp_y0 + subwindow.ny;
						datasec.x1 = x1;
						datasec.x2 = x2;
						datasec.y1 = y1;
						datasec.y2 = y2;
                        npixels = (amp_yf-amp_y0)*(amp_xf-amp_x0);
					} else if (namps == 1) {
                        x1 = datasec.x1;
                        x2 = datasec.x2;
                        y1 = datasec.y1;
                        y2 = datasec.y2;
                        amp_x0 = x1 + subwindow.x0;
                        amp_y0 = y1 + subwindow.y0;
                        amp_xf = amp_x0 + subwindow.nx;
                        amp_yf = amp_y0 + subwindow.ny;
                        npixels = (amp_yf-amp_y0)*(amp_xf-amp_x0);
                        if (verbose) {
                            cout << "operaGain: amp " << amp << " datasec = " << x1 << " : " << x2 << " , "<< y1  << " : "<< y2 << endl;
                        }
                    }
					
					if(amp_x0 < x1 || amp_y0 < y1 || amp_xf > x2 || amp_xf > y2)
					{
						throw "operaGain: error: subwindow exceeds area of datasec\n";
					}
					
					if (verbose) {
						cout << "operaGain: amp = " << amp << ": xw1 = " << amp_x0 << " xw2 = " << amp_xf << " yw1 = "<< amp_y0  << " yw2 = "<< amp_yf   << " npixels = "<< npixels << endl; 
					}					
				}
				
				biasdata[i] = (float *)biasIn->operaFITSImageClonePixels(amp_x0, amp_y0, amp_xf, amp_yf);
			}
			
			/*
			 * open flat images and load data into vectors
			 */	
			
			for (unsigned i=0; i<flatimgIndex; i++) {	
				flatIn = new operaFITSImage(flatimgs[i], tfloat, READONLY);			
				flatdata[i] = (float *)flatIn->operaFITSImageClonePixels(amp_x0, amp_y0, amp_xf, amp_yf);
			}
			
			operaCCDGainNoise(npixels, biasimgIndex, biasdata, flatimgIndex, flatdata, badpixdata, gainLowestCount, gainHighestCount, gainMaxNBins, gainMinPixPerBin, &gain[amp], &gainError[amp], &bias[amp], &noise[amp]);
			
			if (isnan(gain[amp]))
				throw operaException("operaGain: gain: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(gain[amp]))
				throw operaException("operaGain: gain: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(noise[amp]))
				throw operaException("operaGain: noise: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(noise[amp]))
				throw operaException("operaGain: noise: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(gainError[amp]))
				throw operaException("operaGain: gainError: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(gainError[amp]))
				throw operaException("operaGain: gainError: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			
			// Note that ordersstores amps as zero-based....
			orders->getGainBiasNoise()->setGain(amp-1, gain[amp]);
			orders->getGainBiasNoise()->setGainError(amp-1, gainError[amp]);
			orders->getGainBiasNoise()->setNoise(amp-1, noise[amp]);
			orders->getGainBiasNoise()->setBias(amp-1, bias[amp]);
			orders->getGainBiasNoise()->setDatasec(amp-1, datasec);
			if (verbose) {
				cout << "operaGain: amp " << amp << " gain " << gain[amp] << " gainError " << gainError[amp] << " Noise " << noise[amp]  << " Bias " << bias[amp] << endl; 
			}
		}
		orders->getGainBiasNoise()->setAmps(namps);
		orders->WriteSpectralOrders(output, GainNoise);
	}
	catch (operaException e) {
		cout << "operaGain: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cout << "operaGain: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


/* Print out the proper program usage syntax */
void printUsageSyntax(char * modulename) {
	
	cout <<
	" Usage: "+string(modulename)+"  [-vdth] \n"
	"  -b, --biasimgs=<FITS_FILE>, bias fits file name  \n"
	"  -f, --flatimgs=<FITS_FILE>, flat fits file name  \n"
	"  -B, --listofbiasimgs=<LIST_FILE>, list of bias fits file paths  \n"
	"  -F, --listofflatimgs=<LIST_FILE>, list of flat fits file paths  \n"
	"  -o, --output=<FILENAME>, output file  \n"
	"  -m, --badpixelmask=<FITS_FILE>, badpixel mask fits file name \n"
	"  -a, --numberofamplifiers=<INT_VALUE>, amplifier mode 1 or 2 amps \n"
	"  -G, --defaultgain=<FLT_VALUE>, default gain value  \n"
	"  -N, --defaultnoise=<FLT_VALUE>, default noise value  \n"
	"  -n, --gainMinPixPerBin=<INT_VALUE>\n"
	"  -x, --gainMaxNBins=<INT_VALUE>\n"
	"  -l, --gainLowestCount=<UNS_VALUE>\n"
	"  -c, --gainHighestCount=<UNS_VALUE>\n"
	"  -s, --subwindow=<X0 NX Y0 NY>, subwindow to estimate the gain/noise \n"
    "  -D, --useDATASECKey=<BOOL>, Whether or not to use DATASEC keyword\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"\n\n"
	" Example: "+string(modulename)+"  -B bias.list -F flat.list -m badpixmask.fits -a 1 -n 1000 -x 100 -l 1000 -c 30000 --subwindow=\"100 800 100 3000\" \n\n"
    ;
}


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

#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaCCD.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

/*! \file operaGain.cpp */

using namespace std;

/*!
 * operaGain
 * \author Doug Teeple & Eder Martioli
 * \brief Output gain and noise based on at least 2 flat images.
 * \arg argc
 * \arg argv
 * \arg --biasimgs=...
 * \arg --flatimgs=...
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
	const unsigned MAXIMAGES=1000;
	const unsigned MAXNAMPS=2;
	
	operaArgumentHandler args;
    
	string listofbiasimgs;
	string listofflatimgs;
	string biaslistfile;
	string flatlistfile;
	string output; 
	string badpixelmask;
	unsigned namps = 1; // 1 for EEV1 and OLAPAa, 2 for OLAPAab	
	double defaultgain = 1.0;
	double defaultnoise = 0.0;
	unsigned gainMinPixPerBin = 100;
	unsigned gainMaxNBins = 1000;
	unsigned gainLowestCount = 1000;
	unsigned gainHighestCount = 25000;
	string subwindow_str;
	string datasec_str;
	string dseca_str;
	string dsecb_str;
	//unsigned maximages = 12;    // otherwise we run out of memory on 32 bit systems...
	
	args.AddOptionalArgument("biasimgs", listofbiasimgs, "", "List of bias fits files, seperated by spaces");
	args.AddOptionalArgument("flatimgs", listofflatimgs, "", "List of flat fits files, seperated by spaces");
	args.AddOptionalArgument("biaslistfile", biaslistfile, "", "File containing list of bias fits files");
	args.AddOptionalArgument("flatlistfile", flatlistfile, "", "File containing list of flat fits files");
	args.AddRequiredArgument("output", output, "Output file");
	args.AddRequiredArgument("badpixelmask", badpixelmask, "Badpixel mask fits file");	
	args.AddRequiredArgument("numberofamplifiers", namps, "Number of amplifiers (1 or 2)");
	args.AddRequiredArgument("defaultgain", defaultgain, "Default gain value");
	args.AddRequiredArgument("defaultnoise", defaultnoise, "Default noise value");
	args.AddRequiredArgument("gainMinPixPerBin", gainMinPixPerBin, "Minimum number of pixels allowed per bin");
	args.AddRequiredArgument("gainMaxNBins", gainMaxNBins, "Maximum number of bins to use");
	args.AddRequiredArgument("gainLowestCount", gainLowestCount, "Lowest pixel count value to be considered");
	args.AddRequiredArgument("gainHighestCount", gainHighestCount, "Highest pixel count value to be considered");
	args.AddOptionalArgument("subwindow", subwindow_str, "", "Subwindow where data will used for the gain/noise calculation \"x1 x2 y1 y2\"");
	args.AddOptionalArgument("DATASEC", datasec_str, "", "Valid data region on the amplifier (1 amp mode) \"x1 x2 y1 y2\"");
	args.AddOptionalArgument("DSECA", dseca_str, "", "Valid data region on amplifier A (2 amp mode) \"x1 x2 y1 y2\"");
	args.AddOptionalArgument("DSECB", dsecb_str, "", "Valid data region on amplifier B (2 amp mode) \"x1 x2 y1 y2\"");
	
	try {
		args.Parse(argc, argv);

		struct subwindow {
			int x0, nx;
			int y0, ny;
		} subwindow = {0,0,0,0};
		DATASEC_t datasec = {1,2068,1,4608};
		DATASEC_t dseca   = {21,1044,1,4608};
		DATASEC_t dsecb   = {1045,2068,1,4608};
		if (!subwindow_str.empty()) sscanf(subwindow_str.c_str(), "%d %d %d %d", &subwindow.x0, &subwindow.nx, &subwindow.y0, &subwindow.ny);
		if (!datasec_str.empty()) sscanf(datasec_str.c_str(), "%u %u %u %u", &datasec.x1, &datasec.x2, &datasec.y1, &datasec.y2);
		if (!dseca_str.empty()) sscanf(dseca_str.c_str(), "%u %u %u %u", &dseca.x1, &dseca.x2, &dseca.y1, &dseca.y2);
		if (!dsecb_str.empty()) sscanf(dsecb_str.c_str(), "%u %u %u %u", &dsecb.x1, &dsecb.x2, &dsecb.y1, &dsecb.y2);

		string biasimgs[MAXIMAGES];
		string flatimgs[MAXIMAGES];
		unsigned biasimgIndex = 0;
		unsigned flatimgIndex = 0;
		SplitStringIntoArray(listofbiasimgs, biasimgs, biasimgIndex, MAXIMAGES); // Split list of bias images into array
		SplitStringIntoArray(listofflatimgs, flatimgs, flatimgIndex, MAXIMAGES); // Split list of flat images into array
		ReadStringsFromFileIntoArray(biaslistfile, biasimgs, biasimgIndex, MAXIMAGES); // Read list of images from file
		ReadStringsFromFileIntoArray(flatlistfile, flatimgs, flatimgIndex, MAXIMAGES); // Read list of images from file
		
		float gain[MAXNAMPS];
		float noise[MAXNAMPS];
		float gainError[MAXNAMPS];
		float bias[MAXNAMPS];
		for(unsigned i=0;i<MAXNAMPS;i++) {
			gain[i] = defaultgain;	
			noise[i] = defaultnoise;
			gainError[i] = 0.0;	
			bias[i] = 0.0;	
		}
		
		// we need at least 2 biases...
		if (biasimgIndex < 2) {
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
		
		if (args.debug) {
			for(unsigned i=0; i<biasimgIndex; i++)
				cout << "operaGain: input bias: " << i << ' ' << biasimgs[i] << endl;
			for(unsigned i=0; i<flatimgIndex;i ++)
				cout << "operaGain: input flat: " << i << ' ' << flatimgs[i] << endl;		
		}
		
		if (args.verbose) {
			cout << "operaGain: subwindow = " << subwindow.x0 << " " << subwindow.nx << " "<< subwindow.y0  << " "<< subwindow.ny << endl; 
			cout << "operaGain: gainMinPixPerBin = " << gainMinPixPerBin << " gainMinPixPerBin = " << gainMinPixPerBin << " gainLowestCount = " << gainLowestCount << " gainHighestCount = " << gainHighestCount << endl; 
		}
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
                            if (args.verbose) {
                                cout << "operaGain: amp " << amp << " dseca = " << x1 << " : " << x2 << " , "<< y1  << " : "<< y2 << endl;
                            }					
						} else {
                            x1 = dsecb.x1;
                            x2 = dsecb.x2;
                            y1 = dsecb.y1;
                            y2 = dsecb.y2;
							if (args.verbose) {
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
                        if (args.verbose) {
                            cout << "operaGain: amp " << amp << " datasec = " << x1 << " : " << x2 << " , "<< y1  << " : "<< y2 << endl;
                        }
                    }
					
					if(amp_x0 < x1 || amp_y0 < y1 || amp_xf > x2 || amp_xf > y2) {
						throw "operaGain: error: subwindow exceeds area of datasec\n";
					}
					
					if (args.verbose) {
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
			
			if (isnan(gain[amp])) throw operaException("operaGain: gain: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(gain[amp])) throw operaException("operaGain: gain: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(noise[amp])) throw operaException("operaGain: noise: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(noise[amp])) throw operaException("operaGain: noise: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(gainError[amp])) throw operaException("operaGain: gainError: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(gainError[amp])) throw operaException("operaGain: gainError: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			
			// Note that ordersstores amps as zero-based....
			orders->getGainBiasNoise()->setGain(amp-1, gain[amp]);
			orders->getGainBiasNoise()->setGainError(amp-1, gainError[amp]);
			orders->getGainBiasNoise()->setNoise(amp-1, noise[amp]);
			orders->getGainBiasNoise()->setBias(amp-1, bias[amp]);
			orders->getGainBiasNoise()->setDatasec(amp-1, datasec);
			if (args.verbose) {
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

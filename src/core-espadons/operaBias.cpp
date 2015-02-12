/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaBias
 Version: 1.0
 Description: This module calculates the chip bias.
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

#include <string.h>
#include <getopt.h>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaBias.h"
#include "core-espadons/operaGain.h"							// for MAXIMAGES
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaFITSImage.h"

#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaLibCommon.h"

/*! \file operaBias.cpp */

using namespace std;

/*!
 * operaBias
 * \author Doug Teeple / Eder Martioli
 * \brief Calculate bias levels of each amp and update the GainNoiseBias stats file.
 * \arg argc
 * \arg argv
 * \note --bias=...
 * \note 
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOutput
 * \ingroup tools
 * \return EXIT_STATUS
 */

int main(int argc, char *argv[])
{
	int opt;
	string gainNoiseBias;
	string output;
	string biasimgs[MAXIMAGES];			
	unsigned namps = 0;
	unsigned index = 0;
	unsigned nampsexpected = 1;
	const unsigned overscan = 20;
	const unsigned AmpBStartCol = 1045;
	const unsigned AmpBEndCol = 2068;
	const unsigned CCDX = 2080;
	const unsigned CCDY = 4608;
	bool dooverscan = false;			// whether to use the median of the whole array or just the ast column
										// of ampA and the first column of ampB (reduce troughing)
	
	int debug=0, verbose=0, trace=0, plot=0;
	
	struct option longopts[] = {
		{"bias",                1, NULL, 'b'},	// bias		
		{"gain",                1, NULL, 'g'},	// gain file		
		{"output",              1, NULL, 'o'},	// output .bias file		
		{"overscan",            0, NULL, 's'},	// do overscan only	
		{"numberofamplifiers",	1, NULL, 'a'},	// how many amps do we expect		
		
		{"plot",                optional_argument, NULL, 'p'},       
		{"verbose",             optional_argument, NULL, 'v'},
		{"debug",               optional_argument, NULL, 'd'},
		{"trace",               optional_argument, NULL, 't'},
		{"help",                no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	try  {
		while ((opt = getopt_long(argc, argv, "g:b:a:o:sv::d::t::p::h", longopts, NULL))  != -1) {
			switch (opt) {
				case 'b':		// bias name
					biasimgs[index++] = optarg;
					break;						
				case 's':
					dooverscan = true;
					break;						
				case 'g':		// gain / noise / bias file
					gainNoiseBias = optarg;
					break;						
				case 'a':		// setting expectations
					nampsexpected = atoi(optarg);
					break;						
				case 'o':
					output = optarg;
					break;						
					
				case 'v':
					verbose = 1;
					break;
				case 'p':
					plot = 1;
					break;
				case 'd':
					debug = true;
					break;
				case 't':
					trace = true; 
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
		
		
		if (output.empty()) {
			throw operaException("operaBias: please specify an output ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (biasimgs[0].empty()) {
			throw operaException("operaBias: please specify a bias ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		operaSpectralOrderVector *spectralOrders = new operaSpectralOrderVector();
		if (!gainNoiseBias.empty()) {
			spectralOrders->readGainNoise(gainNoiseBias);
			namps = spectralOrders->getGainBiasNoise()->getAmps();
		}
		if (verbose) {
			cout << "operaBias: gain " << gainNoiseBias << endl;
			cout << "operaBias: namps " << namps << endl;
			cout << "operaBias: nampsexpected " << nampsexpected << endl;
			cout << "operaBias: do overscan only " << dooverscan << endl;
		}
		
		for (unsigned i=0; i<index; i++) {
			if (verbose) {
				cout << "operaBias: bias " << biasimgs[i] << endl;
			}
			// get the two median bias levels
			operaFITSImage bias(biasimgs[i], READONLY);
			unsigned short medianBiasB = 0;
			unsigned short medianBiasA = 0;
			
			if (dooverscan) {
				// ampA
				unsigned short *pixels = new unsigned short[CCDY * overscan];
				unsigned short *p = pixels;
				for (unsigned y=0; y<CCDY; y++) {
					for (unsigned x=0; x<overscan; x++) {
						*p++ = bias.getpixelUSHORT(x, y);
					}
				}
				medianBiasA = operaArrayMedianQuickUSHORT(CCDY, pixels);	// scrambles the pixels
																			// ampB
				p = pixels;
				for (unsigned y=0; y<CCDY; y++) {
					for (unsigned x=AmpBEndCol; x<CCDX; x++) {
						*p++ = bias.getpixelUSHORT(x, y);
					}
				}
				medianBiasB = operaArrayMedianQuickUSHORT(CCDY, pixels);	// scrambles the pixels
				free(pixels);
			} else {
				// ampA
				unsigned short *pixels = new unsigned short[(AmpBStartCol-overscan) * CCDY];
				unsigned short *p = pixels;
				for (unsigned y=0; y<CCDY; y++) {
					for (unsigned x=overscan; x<AmpBStartCol; x++) {
						*p++ = bias.getpixelUSHORT(x, y);
					}
				}
				medianBiasA = operaArrayMedianQuickUSHORT((AmpBStartCol-overscan) * CCDY, pixels);	// scrambles the pixels
																									// ampB
				p = pixels;
				for (unsigned y=0; y<CCDY; y++) {
					for (unsigned x=AmpBStartCol; x<AmpBEndCol; x++) {
						*p++ = bias.getpixelUSHORT(x, y);
					}
				}
				medianBiasB = operaArrayMedianQuickUSHORT((AmpBEndCol-AmpBStartCol) * CCDY, pixels);	// scrambles the pixels
				free(pixels);
			}
			
			if (isnan(medianBiasA))
				throw operaException("operaBias: bias 0 : ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(medianBiasA))
				throw operaException("operaBias: bias 0: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (isnan(medianBiasB))
				throw operaException("operaBias: bias 1: ", operaErrorIsNaN, __FILE__, __FUNCTION__, __LINE__);	
			if (isinf(medianBiasB))
				throw operaException("operaBias: bias 1: ", operaErrorIsInf, __FILE__, __FUNCTION__, __LINE__);	
			if (verbose) {
				cout << "operaBias: Changing Bias ampA from " << spectralOrders->getGainBiasNoise()->getBias(0) << " to " << medianBiasA  << endl;					
				if (nampsexpected > 1) {
					cout << "operaBias: Changing Bias ampB from " << spectralOrders->getGainBiasNoise()->getBias(1) << " to " << medianBiasB  << endl;					
				}
			}
			// set the gain / noise / bias
			spectralOrders->getGainBiasNoise()->setBias(0, medianBiasA);
			spectralOrders->getGainBiasNoise()->setBias(1, medianBiasB);
			bias.operaFITSImageClose();
		}
		if (namps < 2)
			spectralOrders->getGainBiasNoise()->setAmps(nampsexpected);
		spectralOrders->WriteSpectralOrders(output, GainNoise);
		delete spectralOrders;
	}
	catch (operaException e) {
		cerr << "operaBias: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaBias: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax() {
	cout <<
	"\n"
	" Usage: operaBias --dooverscan=1|0 --gainNoiseBias=...gain.gz --nampsexpected=1|2 --output=... [--bias=image]* -[dvpth]\n";
}	


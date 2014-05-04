/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaSNR
 Version: 1.0
 Description: Calculate SNR stats.\.
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

#include <stdio.h>
#include <getopt.h>
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t

#include "core-espadons/operaSNR.h"

/*! \file operaSNR.cpp */

using namespace std;

/*! 
 * operaSNR
 * \author Doug Teeple
 * \brief Signal to Noise calculation.
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
	
	string inputfilename; 
	string outputfilename; 
	string wcalfilename; 
	string object;
	
	operaSpectralOrder_t spectralOrderType = SNR;
	
	bool debug=false, verbose=false, trace=false, plot=false;
	
	bool centralsnr = false;
	
	struct option longopts[] = {
		{"input",						1,	NULL, 'i'},
		{"output",						1,	NULL, 'o'},
		{"wavelengthCalibration",		1,	NULL, 'w'},
		{"object",						1,	NULL, 'O'},	// needed for Libre-Esprit output
		{"spectrumtype",				1,  NULL, 'T'},	
		{"centralsnr",					1,  NULL, 'C'},	
		
		{"plot",		optional_argument, NULL, 'p'},       
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:w:T:O:C:v::d::p::t::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputfilename = optarg;
				break;    
			case 'o':
				outputfilename = optarg;
				break; 
			case 'w':
				wcalfilename = optarg;
				break; 
			case 'T':		// spectrum type
				spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'O':
				object = optarg;
				break;   
			case 'C':
				centralsnr = atoi(optarg) == 1;
				break;   
				
			case 'p':
				plot = true;
				break;
			case 'v':
				verbose = true;
				break;
			case 'd':
				debug = true;
				break;
			case 't':
				trace = true;
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
	try {
		if (inputfilename.empty()) {
			throw operaException("operaSNR: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputfilename.empty()) {
			throw operaException("operaSNR: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (wcalfilename.empty()) {
			throw operaException("operaSNR: wcal: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cout << "operaSNR: input = " << inputfilename << endl; 
			cout << "operaSNR: object = " << object << endl; 
			cout << "operaSNR: output = " << outputfilename << endl;
			cout << "operaSNR: spectrum type = " << spectralOrderType << endl;							
			cout << "operaSNR: centralsnr = " << centralsnr << endl;							
			cout << "operaSNR: wavelength calibration file = " << wcalfilename << endl;
		}
		
		operaSpectralOrderVector spectralOrderVector(inputfilename);
        spectralOrderVector.ReadSpectralOrders(wcalfilename);
		spectralOrderVector.setObject(object);
		
		unsigned minorder = spectralOrderVector.getMinorder();
		unsigned maxorder = spectralOrderVector.getMaxorder();
		for (unsigned order=minorder; order <= maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrderVector.GetSpectralOrder(order);
			if (spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
				operaSpectralElements *spectralElements = spectralOrder->getSpectralElements(); 
				if (spectralElements->getnSpectralElements() > 0) {
					operaWavelength *wavelength = spectralOrder->getWavelength();  
					Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
					unsigned elements = spectralElements->getnSpectralElements();
					while (elements--) {
						spectralElements->setwavelength(wavelengthPolynomial->Evaluate(spectralElements->getdistd(elements)), elements);
					}
					spectralOrder->sethasWavelength(true);   
					spectralElements->setHasWavelength(true);   
					spectralElements->setHasDistance(false);
					spectralOrder->sethasCenterSNROnly(centralsnr);
					spectralOrder->calculateSNR();	// has side effect of retaining center SNR					
				}
			}
		}
		spectralOrderVector.WriteSpectralOrders(outputfilename, spectralOrderType);
	}
	catch (operaException e) {
		cerr << "operaSNR: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaSNR: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth] --outputfilename=... --inputfilename=... \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -p, --plot, Turn on plotting \n"
	"  -o, --outputfilename=$(spectradir)$*i.sn\n"
	"  -i, --inputfilename=$(spectradir)$*iu.s\n\n";
}


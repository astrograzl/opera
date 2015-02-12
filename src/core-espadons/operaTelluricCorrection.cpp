/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ****
 ********************************************************************
 Module name: operaTelluricCorrection
 Version: 1.0
 Description: This is a template for helping developers 
 to start up with an OPERA module.
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
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>
#include <fitsio.h>
#include <getopt.h>
#include <errno.h>
#include <string>
#include <iostream>

#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "core-espadons/operaTelluricWavelengthCorrection.h"

/*! \brief telluric correction. */
/*! \file operaTelluricCorrection.cpp */
/*! \package operaTelluricCorrection */

using namespace std;
#if 0
static int debug=0, verbose=0, trace=0;
#endif
/*! 
 * operaTelluricCorrection
 * \author Doug Teeple
 * \brief telluric correction.
 * \arg argc
 * \arg argv
 * \note --output=...
 * \note --input=...
 * \note --wave=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */


int main(int argc, char *argv[])
{
#if 0
	
???? DT May 8 2014 this doesnt seem right at all????
	int opt;
	string wavelengthcalibrationfilename;
	string inputfilename;
	string outputfilename;
	
	operaSpectrum_t spectrumType = RawFluxSum;
	operaSpectralOrder_t spectralOrderType = RawSpectrum;
	
	struct option longopts[] = {
		{"input",1, NULL, 'i'},
		{"output",1, NULL, 'o'},
		{"wave",1, NULL, 'w'},
		{"spectrumtype",1, NULL, 'T'},	
		
		{"verbose",0, NULL, 'v'},
		{"debug",0, NULL, 'd'},
		{"trace",0, NULL, 't'},
		{"help",0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:w:T:vdth",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputfilename = optarg;	
				break;    
			case 'w':
				wavelengthcalibrationfilename = optarg;	
				break;    
			case 'o':		// output
				outputfilename = optarg;
				break;
			case 'T':		// spectrum type
				spectrumType = (operaSpectrum_t)atoi(optarg);
				switch (spectrumType) {
					case RawFluxSum:
						spectralOrderType = RawSpectrum;
						break;
					case StandardFlux:
						spectralOrderType = StandardSpectrum;
						break;
					case OptimalFlux:
						spectralOrderType = OptimalSpectrum;
						break;
					case OperaOptimalFlux:
						spectralOrderType = OperaOptimalSpectrum;
						break;
					default:
						break;
				}
				break;
				
			case 'v':
				verbose = 1;
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
	
	try {
		// we need an input...
		if (inputfilename.empty()) {
			throw operaException("operaNormalize: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		operaSpectralOrderVector spectralOrdervector(inputfilename);
		spectralOrdervector.ReadSpectralOrders(wavelengthcalibrationfilename);
		unsigned minorder = spectralOrdervector.getMinorder();
		unsigned maxorder = spectralOrdervector.getMaxorder();
		for (unsigned order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrdervector.GetSpectralOrder(order);
			if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
				operaSpectralElements *spectralelements = spectralOrder->getSpectralElements();
				operaWavelength *wavelength = spectralOrder->getWavelength();
				Polynomial *p = wavelength->getWavelengthPolynomial();
				double *wl = spectralelements->getwavelengthVector();
				double *dv = spectralelements->getdistdVector();
				for (unsigned k=0; k<spectralelements->getnSpectralElements(); k++) {
					*wl++ = p->Evaluate(*dv++);
				}
			}						
		}
		spectralOrdervector.WriteSpectralOrders(outputfilename, CalibratedSpectrum);
	}
	catch (operaException e) {
		cerr << "operaTelluricCorrection: " << e.getFormattedMessage() << endl;
	}
	catch (...) {
		cerr << "operaTelluricCorrection: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
#endif	
	return EXIT_SUCCESS;
} 
#if 0
/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cerr <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth] --param=<PARAMETER_VALUE_1> --param=<PARAMETER_VALUE_2> ... --output=<PRODUCT_FILE_NAME> --input=<INPUT_FILE_1> --input=<INPUT_FILE_2> ... \n\n"
	" Example: "+string(modulename)+"  -v -p 10 -p 2 --output=o.fits -i 001.fits -i 002.fits -i bad_pix.dat \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -p, --param=<PARAMETER_VALUE>, Input parameters  \n"
	"  -o, --output=<PRODUCT_FILE_NAME>, Output product file  \n"
	"  -i, --input=<INPUT_FILE_NAME>, Input files  \n\n";
}
#endif

/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaGenerateLEFormats
 Version: 1.0
 Description:  Module to generate LE compatible formats.
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
#include "core-espadons/operaGenerateLEFormats.h"

#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaWavelength.h"				// for wavelength polynomial
#include "libraries/operaLib.h"						// systemf
#include "libraries/operaLibCommon.h"               // for anglecoord_t and timecoord_t
#include "libraries/operaHelio.h"					// for sexigesimal conversion
#include "libraries/operaSpectralTools.h"			// 
#include "libraries/operaFit.h"						// for operaLMFitPolynomial

/*! \file operaGenerateLEFormats.cpp */

using namespace std;

/*! 
 * operaGenerateLEFormats
 * \author Eder Martioli
 * \brief  Module to generate LE compatible formats
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


#define NOTPROVIDED -999

int main(int argc, char *argv[])
{
	int opt;
    
	string inputOperaSpectrum; 
	string outputLEfilename;
	string object = "Nowhere";

    operaSpectralOrder_t LibreEspritSpectrumType = LibreEspritsp2Spectrum;
    /*  Available LibreEspritSpectrumType options for LE formats are:
        LibreEspritpolarimetry
        LibreEspritpolSpectrum
        LibreEspritsp1Spectrum
        LibreEspritsp2Spectrum
     */
    
    operaFluxType_t fluxType = RawFluxInElectronsPerElement;
    /*  Available operaFluxType_t options are:
		1 = RawFluxInElectronsPerElement
		2 = NormalizedFluxToContinuum
		3 = CalibratedFluxNormalizedToRefWavelength
    */
    
    operaWavelengthType_t wavelengthType = ThArCalibratedInNM;
    /*  Available operaWavelengthType_t options are:
		1 = ThArCalibratedInNM
		2 = TelluricCorrectedWavelengthInNM
		3 = RVCorrectedWavelengthInNM
        4 = RVAndTelluricCorrectedWavelengthInNM
     */
    
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    int maxorder = 62;    
    bool minorderprovided = false;
    bool maxorderprovided = false; 
    
	bool debug=false, verbose=false, trace=false, plot=false;
	
	struct option longopts[] = {
		{"inputOperaSpectrum",          1, NULL, 'i'},	// .spc
		{"outputLEfilename",            1, NULL, 'l'},  // .s
		{"object",                      1, NULL, 'o'},  // object name
		{"LibreEspritSpectrumType",		1, NULL, 'S'},	// spectrum type
		{"fluxType",                    1, NULL, 'F'},	// flux type
		{"wavelengthType",				1, NULL, 'W'},	// wavelength type

		{"ordernumber",					1, NULL, 'O'},	// just do a particular order
		{"minorder",					1, NULL, 'M'},	// only consider this order range
		{"maxorder",					1, NULL, 'X'}, 	// only consider this order range

		{"plot",						0, NULL, 'p'},
		{"verbose",						0, NULL, 'v'},
		{"debug",						0, NULL, 'd'},
		{"trace",						0, NULL, 't'},
		{"help",						0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:l:o:S:F:W:O:M:X:p::v::d::t::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputOperaSpectrum = optarg;
				break;   
			case 'l':
				outputLEfilename = optarg;
                break; 		
			case 'o':
				object = optarg;
                break;
			case 'S':		// spectrum type
				LibreEspritSpectrumType = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'F':		// flux type
				fluxType = (operaFluxType_t)atoi(optarg);
				break;
			case 'W':		// wavelength type
				wavelengthType = (operaWavelengthType_t)atoi(optarg);
				break;
			case 'O':
				ordernumber = atoi(optarg);
				break;				
			case 'M':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
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
		// we need an input .e spectrum...
		if (inputOperaSpectrum.empty()) {
			throw operaException("operaGenerateLEFormats: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputLEfilename.empty()) {
			throw operaException("operaGenerateLEFormats: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
        
		if (verbose) {
			cout << "operaGenerateLEFormats: inputOperaSpectrum = " << inputOperaSpectrum << endl; 
			cout << "operaGenerateLEFormats: outputLEfilename = " << outputLEfilename << endl;
			cout << "operaGenerateLEFormats: object = " << object << endl;
			cout << "operaGenerateLEFormats: LibreEspritSpectrumType = " << LibreEspritSpectrumType << endl;
			cout << "operaGenerateLEFormats: fluxType = " << fluxType << endl;
			cout << "operaGenerateLEFormats: wavelengthType = " << wavelengthType << endl;
		}
        
		/*
		 * Down to business, read in all the source and calibration data.
		 */        
		operaSpectralOrderVector spectralOrders(inputOperaSpectrum);
		if(!minorderprovided) {
            minorder = spectralOrders.getMinorder();
        }
        
        if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}        
		
		if(!maxorderprovided) {
            maxorderprovided = spectralOrders.getMaxorder();
        }
        
        if(ordernumber != NOTPROVIDED) {
			minorder = ordernumber;
			maxorder = ordernumber;
		}
        
        if (verbose)
			cout << "operaGenerateLEFormats: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			if (spectralOrder->gethasSpectralElements() ) {
                operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                if (spectralElements->getHasExtendedBeamFlux()){
                    
                    switch (fluxType) {
                        case RawFluxInElectronsPerElement: {
                            //spectralElements->copyFROMrawFlux();
                        }
                            break;
                        case NormalizedFluxToContinuum: {
                            spectralElements->copyFROMnormalizedFlux();
                        }
                            break;
                        case CalibratedFluxNormalizedToRefWavelength: {
                            spectralElements->copyFROMfcalFlux();
                        }
                            break;
                        default:
                            break;
                    }
                    
                    switch (wavelengthType) {
                        case ThArCalibratedInNM:
                            break;
                        case TelluricCorrectedWavelengthInNM: {
                            spectralElements->copyFROMtell();
                        }
                            break;
                        case RVCorrectedWavelengthInNM: {
                            spectralOrder->applyBarycentricWavelengthCorrectionFromExtendedRvel();
                        }
                        case RVAndTelluricCorrectedWavelengthInNM: {
                            spectralElements->copyFROMtell();
                            spectralOrder->applyBarycentricWavelengthCorrectionFromExtendedRvel();
                        }
                            break;
                        default:
                            break;
                    }
                }
                
			}
		}
        //---------------------------------
        
 		// output wavelength/flux calibrated spectrum...
		spectralOrders.setObject(object);
		spectralOrders.WriteSpectralOrders(outputLEfilename, LibreEspritSpectrumType);
	}
	catch (operaException e) {
		cerr << "operaGenerateLEFormats: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaGenerateLEFormats: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	"  -i, --inputOperaSpectrum=<EXTENDED OPERA .SPC>"
	"  -l, --outputLEfilename=<FILE_NAME>"
	"  -o, --object=<UNS_VALUE>"
	"  -S, --LibreEspritSpectrumType=<FITS_IMAGE>, .e\n"
	"  -F, --fluxType=<UNS_VALUE>, Method for extraction\n"
	"  -W, --wavelengthType=<FILE_NAME>, wavelength calibration polynomials\n\n";
}

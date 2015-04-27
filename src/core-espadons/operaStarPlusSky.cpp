/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaStarPlusSky
 Version: 1.0
 Description: Extract raw spectrum in Star Plus Sky mode
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
#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaStarPlusSky.h"

#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaWavelength.h"				// for wavelength polynomial
#include "libraries/operaLib.h"						// systemf
#include "libraries/operaLibCommon.h"               // for anglecoord_t and timecoord_t
#include "libraries/operaHelio.h"					// for sexigesimal conversion

#define NOTPROVIDED -999

/*! \file operaStarPlusSky.cpp */

using namespace std;

/*! 
 * operaStarPlusSky
 * \author Eder Martioli
 * \brief Module to extract raw spectra in Star Plus Sky mode.
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
    
	string input; 
	string object;
	string outputSpectraFile;
	string wavelengthCalibration;
	string inputFlatFluxCalibration;
	operaSpectralOrder_t spectralOrderType = CalibratedExtendedBeamSpectrum;

    string inputWavelengthMaskForUncalContinuum;
    unsigned numberOfPointsInUniformSample = 150;
    
    /*
     * Parameters for normalization
     */
    unsigned normalizationBinsize = 100;
    
    /*
     * Barycentric radial velocity correction
     */    
	string radialvelocitycorrection;
    /*
     * Telluric correction
     */    
	string telluriccorrection;
	
    /*
     * Parameters for flux calibration
     */    
    string fluxCalibration;
    string flatResponse;
    
	float exposureTime = 0.0;
    
    bool AbsoluteCalibration = false;
    
    double delta_wl = 1.0; // wavelength (in nm) range for stiching non-overlapping orders

 	bool starplusskyInvertSkyFiber = false;
    
    int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    int maxorder = 62;    
    bool minorderprovided = false;
    bool maxorderprovided = false; 
    
    bool interactive = false;
    
	bool debug=false, verbose=false, trace=false, plot=false;
    
    string plotfilename;	
	string spectrumDataFilename, continuumDataFilename;	
	string scriptfilename;	

	struct option longopts[] = {
		{"inputUncalibratedSpectrum",	1, NULL, 'i'},	// .e
        {"object",						1, NULL, 'o'},	// needed for Libre-Esprit output
		{"outputCalibratedSpectrum",	1, NULL, 's'},  // i.e
		{"spectrumtype",				1, NULL, 'y'},	// spectrum type
		{"wavelengthCalibration",		1, NULL, 'w'},	// wavelength calibration file (.wcal or .wcar)
 		{"radialvelocitycorrection",	1, NULL, 'V'},  // Barycentric wavelength correction file (.rvel)
 		{"telluriccorrection",			1, NULL, 'T'},  // Telluric wavelength correction file (.tell)
 		{"inputFlatFluxCalibration",                1, NULL, 'm'},  // flat field spectrum ff_
        {"inputWavelengthMaskForUncalContinuum",    1, NULL, 'u'},
		{"numberOfPointsInUniformSample",           1, NULL, 'l'},

 		{"normalizationBinsize",		1, NULL, 'b'},	// binsize for normalization
        
        {"fluxCalibration",				1, NULL, 'C'},	// flux calibration; file (.fcal)
        {"flatResponse",				1, NULL, 'R'},	// apply flat response calibration; file (LE .s)
        
        {"etime",						1, NULL, 'E'},	// needed for flux calibration        
		{"AbsoluteCalibration",         1, NULL, 'A'},  // absolute or relative flux calibration

        {"starplusskyInvertSkyFiber",   1, NULL, 'K'},

		{"ordernumber",					1, NULL, 'O'},	// just do a particular order
		{"minorder",					1, NULL, 'M'},	// only consider this order range
		{"maxorder",					1, NULL, 'X'}, 	// only consider this order range
        
        {"plotfilename",				1, NULL, 'P'},
		{"spectrumDataFilename",		1, NULL, 'F'},
		{"continuumDataFilename",		1, NULL, 'c'},
        
		{"scriptfilename",				1, NULL, 'S'},  
		{"interactive",					0, NULL, 'I'},
		
		{"plot",						0, NULL, 'p'},
		{"verbose",						0, NULL, 'v'},
		{"debug",						0, NULL, 'd'},
		{"trace",						0, NULL, 't'},
		{"help",						0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:s:y:w:V:T:m:u:l:b:C:R:E:A:K:O:M:X:P:F:c:S:I:p::v::d::t::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				input = optarg;
				break;   
			case 'o':
				object = optarg;
				break;   
			case 's':
				outputSpectraFile = optarg;
                break; 		
			case 'y':		// spectrum type
				spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
				break;                
			case 'w':
				wavelengthCalibration = optarg;
				break;
			case 'V':		// for telluric wl correction
				radialvelocitycorrection = optarg;
				break;
			case 'T':
				telluriccorrection = optarg;
				break;
			case 'm':
				inputFlatFluxCalibration = optarg;
				break;
			case 'u':
				inputWavelengthMaskForUncalContinuum = optarg;
				break;
			case 'l':
				numberOfPointsInUniformSample = atoi(optarg);
				break;
			case 'b':		// normalization binsize
				normalizationBinsize = atoi(optarg);
				break;
                
			case 'C':       // for flux calibration
				fluxCalibration = optarg;
				break;    
			case 'R':       // for flat-response flux calibration
                flatResponse = optarg;
                break;
                
			case 'E':
				exposureTime = atof(optarg);
				break;		
			case 'A':
				AbsoluteCalibration = atoi(optarg)==1;
				break;
                
            case 'K':
				starplusskyInvertSkyFiber = (atoi(optarg)?true:false);
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
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break; 		                
			case 'F':
				spectrumDataFilename = optarg;
				break; 	
			case 'c':
				continuumDataFilename = optarg;
				break; 	
			case 'S':
				scriptfilename = optarg;
				break;  
			case 'I':		// for interactive plots
				interactive = true;
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
		if (input.empty()) {
			throw operaException("operaStarPlusSky: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputSpectraFile.empty()) {
			throw operaException("operaStarPlusSky: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (wavelengthCalibration.empty()) {
			throw operaException("operaStarPlusSky: wcal: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cout << "operaStarPlusSky: input spectrum = " << input << endl; 
			cout << "operaStarPlusSky: object = " << object << endl; 
			cout << "operaStarPlusSky: output spectrum file = " << outputSpectraFile << endl;
			cout << "operaStarPlusSky: spectrum type = " << spectralOrderType << endl;							
			cout << "operaStarPlusSky: wavelength calibration file = " << wavelengthCalibration << endl;
            cout << "operaStarPlusSky: binsize for normalization = " << normalizationBinsize << endl;  
            cout << "operaStarPlusSky: input flux calibration file = " << fluxCalibration << endl; 
            cout << "operaStarPlusSky: input flat response calibration file = " << flatResponse << endl;
			cout << "operaStarPlusSky: exposure time = " << exposureTime << endl;
			cout << "operaStarPlusSky: absolute calibration = " << AbsoluteCalibration << endl;
			cout << "operaStarPlusSky: starplusskyInvertSkyFiber = " << starplusskyInvertSkyFiber << endl;
            cout << "operaStarPlusSky: radialvelocitycorrection = " << radialvelocitycorrection << endl;
            cout << "operaStarPlusSky: telluriccorrection = " << telluriccorrection << endl;
            cout << "operaStarPlusSky: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
            cout << "operaStarPlusSky: inputWavelengthMaskForUncalContinuum = " << inputWavelengthMaskForUncalContinuum << endl;
            cout << "operaStarPlusSky: numberOfPointsInUniformSample = " << numberOfPointsInUniformSample << endl;

            if (ordernumber != NOTPROVIDED) {
                cout << "operaStarPlusSky: ordernumber = " << ordernumber << endl;            
            }             
            
            if (plot) {
                cout << "operaStarPlusSky: plotfilename = " << plotfilename << endl;
                cout << "operaStarPlusSky: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaStarPlusSky: continuumDataFilename = " << continuumDataFilename << endl;                
                cout << "operaStarPlusSky: scriptfilename = " << scriptfilename << endl; 
                if (interactive) {
                    cout << "operaStarPlusSky: interactive = YES" << endl; 
                } else {
                    cout << "operaStarPlusSky: interactive = NO" << endl; 
                }
            }            
		}
        
		/*
		 * Plotting support
		 */
        // DT May 20 2014 - not used -- unsigned numberOfprintouts = 2; // for 3D plotting        
        ofstream *fspecdata = NULL;
        ofstream *fcontinuumdata = NULL;
        
        if (!spectrumDataFilename.empty()) {
            fspecdata = new ofstream();
            fspecdata->open(spectrumDataFilename.c_str());  
        }    
        
        if (!continuumDataFilename.empty()) {
            fcontinuumdata = new ofstream();
            fcontinuumdata->open(continuumDataFilename.c_str());  
        }          

        /*
         * Down to business, read in all the source and calibration data.
         */
        operaSpectralOrderVector spectralOrders(input);
        spectralOrders.ReadSpectralOrders(wavelengthCalibration);
        
        if(!minorderprovided) {
            minorder = spectralOrders.getMinorder();
        }
        if(!maxorderprovided) {
            maxorder = spectralOrders.getMaxorder();
        }
        if(ordernumber != NOTPROVIDED) {
            minorder = ordernumber;
            maxorder = ordernumber;
        }
        
        if (verbose) {
            cout << "operaStarPlusSky: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        }
        
        int minPossibleOrder = 0;
        int maxPossibleOrder = 0;
        
        for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                
                spectralOrder->getSpectralElements()->CreateExtendedvectors(spectralOrder->getSpectralElements()->getnSpectralElements());
                
                // Save the raw flux for later
                spectralOrder->getSpectralElements()->copyTOrawFlux();
                spectralOrder->getSpectralElements()->copyTOnormalizedFlux();
                spectralOrder->getSpectralElements()->copyTOfcalFlux();
                
                operaWavelength *Wavelength = spectralOrder->getWavelength();
                spectralOrder->getSpectralElements()->setwavelengthsFromCalibration(Wavelength);
                spectralOrder->getSpectralElements()->copyTOtell();
                
                if(order < minPossibleOrder || minPossibleOrder==0) {
                    minPossibleOrder = order;
                }
                if(order > maxPossibleOrder) {
                    maxPossibleOrder = order;
                }
            }
        }
        
        if(minPossibleOrder > minorder) {
            minorder = minPossibleOrder;
            if (verbose)
                cout << "operaStarPlusSky: minorder reset to " << minorder << endl;
        }
        if(maxPossibleOrder < maxorder) {
            maxorder = maxPossibleOrder;
            if (verbose)
                cout << "operaStarPlusSky: maxorder reset to " << maxorder << endl;
        }
        
        unsigned NumberofBeams = spectralOrders.getNumberofBeams(minorder, maxorder);
        
        //---------------------------------
        // Load telluric corrected wavelength calibration
		if (!telluriccorrection.empty()) {
            spectralOrders.readTelluricRVINTOExtendendSpectra(telluriccorrection, minorder, maxorder);
		}

        //---------------------------------
        // Load Barycentric RV wavelength correction and also wavelength calibration
        if (!radialvelocitycorrection.empty()) {
            spectralOrders.readRVCorrectionINTOExtendendSpectra(radialvelocitycorrection, wavelengthCalibration, minorder, maxorder);
        }

        bool StarPlusSky = true;
        
        //---------------------------------
        // Correct flat-field
        if (!inputFlatFluxCalibration.empty()) {
            spectralOrders.correctFlatField(inputFlatFluxCalibration, minorder, maxorder, StarPlusSky, starplusskyInvertSkyFiber);
            spectralOrders.saveExtendedRawFlux(minorder, maxorder);
        }

        //---------------------------------
        // Flux Normalization and Flux Calibration Stuff
        
        /*
         * Flux Calibration Stuff ...
         */        
		if (!inputWavelengthMaskForUncalContinuum.empty()) {
            if (!fluxCalibration.empty()) {
                spectralOrders.normalizeAndCalibrateFluxINTOExtendendSpectra(inputWavelengthMaskForUncalContinuum,fluxCalibration, exposureTime, AbsoluteCalibration,numberOfPointsInUniformSample,normalizationBinsize, delta_wl, minorder, maxorder, false, StarPlusSky);
            } else if (fluxCalibration.empty() && !flatResponse.empty()) {
                spectralOrders.normalizeAndApplyFlatResponseINTOExtendendSpectra(inputWavelengthMaskForUncalContinuum,flatResponse,numberOfPointsInUniformSample,normalizationBinsize, delta_wl, minorder, maxorder, false, StarPlusSky);
            } else {
                spectralOrders.normalizeFluxINTOExtendendSpectra(inputWavelengthMaskForUncalContinuum,numberOfPointsInUniformSample,normalizationBinsize, delta_wl, minorder, maxorder, false);
            }
        } else {
            spectralOrders.normalizeOrderbyOrderAndSaveFluxINTOExtendendSpectra(normalizationBinsize, minorder, maxorder, false);
        }

        // output a wavelength calibrated spectrum...
		spectralOrders.setObject(object);
		spectralOrders.WriteSpectralOrders(outputSpectraFile, spectralOrderType);

        if (fspecdata != NULL) {
			fspecdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateExtractionPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(), NumberofBeams, interactive);
            }
        } 
        
	}
	catch (operaException e) {
		cerr << "operaStarPlusSky: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaStarPlusSky: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	" --input=<input.e>"
	" --outputSpectraFile=<FILE_NAME>"	
	" --spectrumtype=<UNS_VALUE>"
	"  -i, --inputUncalibratedSpectrum=<FITS_IMAGE>, .e\n"
	"  -s, --outputCalibratedSpectrum=<FILE_NAME>, Output file name\n"
	"  -y, --spectrumtype=<UNS_VALUE>, Method for extraction\n"
	"  -w, --wavelengthCalibration=<FILE_NAME>, wavelength calibration polynomials\n"
	"  -V, --radialvelocitycorrection=<FILE_NAME>, Barycentric Radial Velocity Correction\n"
	"  -C, --fluxCalibration=<FILE_NAME>, flux calibration file (.fcal) -> this overrides flatResponse\n"
    "  -R, --flatResponse=<FILE_NAME>, flat response calibration file (LE .s)\n"
	"  -b, --normalizationBinsize=<float>,  binsize for normalization\n"
	"  -B, --orderBin=<int>, number or orders to bin for continuum evaluation\n"
	"  -A, --AbsoluteCalibration=<bool>, perform absolute flux calibration\n"
	"  -K, --starplusskyInvertSkyFiber=<bool>, Star+sky: invert sky fiber (default is beam[0]=star and beam[1]=sky). \n"
	"  -l, --usePolynomial=1|0, option to use polynomial for normalization\n"
	"  -r, --orderOfPolynomial=<unsigned>, option to set degree of polynomial for normalization\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -p, --plot,  Turn on plotting\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateExtractionPlot(const char *gnuScriptFileName, const char *outputPlotEPSFileName,const char *dataFileName, unsigned nbeams, bool display)
{
	FILE *fgnu;
	remove(gnuScriptFileName); // delete any existing file with the same name
	
	fgnu = fopen(gnuScriptFileName,"w");
	
	fprintf(fgnu,"unset key\n");
	fprintf(fgnu,"set view 0,0\n");
	fprintf(fgnu,"set iso 100\n");
	fprintf(fgnu,"set samples 100\n");
	fprintf(fgnu,"set pm3d at s\n");
	fprintf(fgnu,"set ticslevel 0\n");   
	
	fprintf(fgnu,"set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
	fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName);
	
	unsigned fluxcol = 6 + 4*nbeams;
	fprintf(fgnu,"splot \"%s\" u 5:1:%u with pm3d\n",dataFileName,fluxcol);
	
	if (display) {
		fprintf(fgnu,"set output\n");
		fprintf(fgnu,"set terminal x11\n");
		fprintf(fgnu,"replot\n");       
		fclose(fgnu);   
		systemf("gnuplot -persist %s",gnuScriptFileName);
	} else {
		fclose(fgnu);  
	}
}

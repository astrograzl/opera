/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaPolarIntensity
 Version: 1.0
 Description: Create polar intensity spectrum.
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
#include "core-espadons/operaPolarIntensity.h"

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

/*! \file operaPolarIntensity.cpp */

using namespace std;

/*! 
 * operaPolarIntensity
 * \author Eder Martioli
 * \brief Create polar intensity spectrum.
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
    
	string input; 
	string object;
	string outputSpectraFile;
	string wavelengthCalibration;
    string inputFlatFluxCalibration;

    string inputWavelengthMaskForUncalContinuum;
    unsigned numberOfPointsInUniformSample = 150;
    
	operaSpectralOrder_t spectralOrderType = CalibratedExtendedBeamSpectrum;
	
    /*
     * Parameters for normalization
     */
	unsigned normalizationBinsize = 300;
    
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
	float exposureTime = 0.0;
    bool AbsoluteCalibration = false;
    
    double delta_wl = 1.0; // wavelength (in nm) range for stiching non-overlapping orders
    /*
     * Apply final wavelength calibration to stitch orders together
     */     
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
		{"outputCalibratedSpectrum",	1, NULL, 's'},  // .s                                                    
		{"spectrumtype",				1, NULL, 'y'},	// spectrum type
		{"wavelengthCalibration",		1, NULL, 'w'},	// wavelength calibration file (.wcal or .auto)
 		{"radialvelocitycorrection",	1, NULL, 'V'},  // Barycentric wavelength correction file (.rvel)
 		{"telluriccorrection",			1, NULL, 'T'},  // Telluric wavelength correction file (.tell)
 		{"inputFlatFluxCalibration",	1, NULL, 'm'},  // flat field spectrum ff_
        {"inputWavelengthMaskForUncalContinuum",           1, NULL, 'u'},
		{"numberOfPointsInUniformSample",   1, NULL, 'l'},

 		{"normalizationBinsize",		1, NULL, 'b'},	// binsize for normalization                                                   
		
        {"fluxCalibration",				1, NULL, 'C'},	// apply flux calibration; file (.fcal)
        {"etime",						1, NULL, 'E'},	// needed for flux calibration
		{"AbsoluteCalibration",         1, NULL, 'A'},  // absolute or relative flux calibration
            
		{"ordernumber",					1, NULL, 'O'},	// just do a particulare order
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
	
	while((opt = getopt_long(argc, argv, "i:o:s:y:w:V:T:m:u:l:b:C:E:A:O:M:X:P:F:c:S:I:p::v::d::t::h", 
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
			case 'E':
				exposureTime = atof(optarg);
				break;		
			case 'A':
				AbsoluteCalibration = atoi(optarg)==1;
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
			throw operaException("operaPolarIntensity: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputSpectraFile.empty()) {
			throw operaException("operaPolarIntensity: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (wavelengthCalibration.empty()) {
			throw operaException("operaPolarIntensity: wcal: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cout << "operaPolarIntensity: input spectrum = " << input << endl; 
			cout << "operaPolarIntensity: object = " << object << endl; 
			cout << "operaPolarIntensity: output spectrum file = " << outputSpectraFile << endl;
			cout << "operaPolarIntensity: spectrum type = " << spectralOrderType << endl;							
			cout << "operaPolarIntensity: wavelength calibration file = " << wavelengthCalibration << endl;
            cout << "operaPolarIntensity: binsize for normalization = " << normalizationBinsize << endl;  
            cout << "operaPolarIntensity: input flux calibration file = " << fluxCalibration << endl;
            cout << "operaPolarIntensity: inputWavelengthMaskForUncalContinuum = " << inputWavelengthMaskForUncalContinuum << endl;
            cout << "operaPolarIntensity: numberOfPointsInUniformSample = " << numberOfPointsInUniformSample << endl;
			cout << "operaPolarIntensity: exposure time = " << exposureTime << endl;
			cout << "operaPolarIntensity: absolute calibration = " << AbsoluteCalibration << endl;
            cout << "operaPolarIntensity: radialvelocitycorrection = " << radialvelocitycorrection << endl;  
			cout << "operaPolarIntensity: telluriccorrection = " << telluriccorrection << endl;
            cout << "operaPolarIntensity: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
          
            if (ordernumber != NOTPROVIDED) {
                cout << "operaPolarIntensity: ordernumber = " << ordernumber << endl;            
            }             
            
            if (plot) {
                cout << "operaPolarIntensity: plotfilename = " << plotfilename << endl;
                cout << "operaPolarIntensity: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaPolarIntensity: continuumDataFilename = " << continuumDataFilename << endl;                
                cout << "operaPolarIntensity: scriptfilename = " << scriptfilename << endl; 
                if (interactive) {
                    cout << "operaPolarIntensity: interactive = YES" << endl; 
                } else {
                    cout << "operaPolarIntensity: interactive = NO" << endl; 
                }
            }            
		}
        
		/*
		 * Plotting support
		 */
        // DT May 20 2014 -- not used unsigned numberOfprintouts = 2; // for 3D plotting
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
			cout << "operaPolarIntensity: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        unsigned NumberofBeams = spectralOrders.getNumberofBeams(minorder, maxorder);

		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			if (spectralOrder->gethasSpectralElements()) {
				spectralOrder->getSpectralElements()->CreateExtendedvectors(spectralOrder->getSpectralElements()->getnSpectralElements());
                // Save the raw flux for later
                spectralOrder->getSpectralElements()->copyTOrawFlux();
                spectralOrder->getSpectralElements()->copyTOnormalizedFlux();
                spectralOrder->getSpectralElements()->copyTOfcalFlux();
                
                if(spectralOrder->gethasWavelength()) {
                    operaWavelength *Wavelength = spectralOrder->getWavelength();
                    spectralOrder->getSpectralElements()->setwavelengthsFromCalibration(Wavelength);
                    spectralOrder->getSpectralElements()->copyTOtell();
                }
			}
		}

        //---------------------------------
        // Load telluric corrected wavelength calibration
		if (!telluriccorrection.empty()) {
            spectralOrders.readTelluricWavelengthINTOExtendendSpectra(telluriccorrection, minorder, maxorder);
		}
        
        //---------------------------------
        // Load Barycentric RV wavelength correction and also wavelength calibration
        if (!radialvelocitycorrection.empty()) {
            spectralOrders.readRVCorrectionINTOExtendendSpectra(radialvelocitycorrection, wavelengthCalibration, minorder, maxorder);
        }
        
        //---------------------------------
        // Correct flat-field
        if (!inputFlatFluxCalibration.empty()) {
            spectralOrders.correctFlatField(inputFlatFluxCalibration, minorder, maxorder, false);
            spectralOrders.saveExtendedRawFlux(minorder, maxorder);
        }

        //---------------------------------
        // Flux Normalization and Flux Calibration Stuff
		if (!fluxCalibration.empty() && !inputWavelengthMaskForUncalContinuum.empty()) {
            spectralOrders.normalizeAndCalibrateFluxINTOExtendendSpectra(inputWavelengthMaskForUncalContinuum,fluxCalibration, exposureTime, AbsoluteCalibration,numberOfPointsInUniformSample,normalizationBinsize, delta_wl, minorder, maxorder, false, false);
        } else if (!inputWavelengthMaskForUncalContinuum.empty()) {
            spectralOrders.normalizeFluxINTOExtendendSpectra(inputWavelengthMaskForUncalContinuum,numberOfPointsInUniformSample,normalizationBinsize, delta_wl, minorder, maxorder, false);
        }

        //---------------------------------
 		// output wavelength/flux calibrated spectrum...
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
		cerr << "operaPolarIntensity: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPolarIntensity: " << operaStrError(errno) << endl;
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
	"  -b, --normalizationBinsize=<float>,  binsize for normalization\n"
	"  -A, --AbsoluteCalibration=<bool>, perform absolute flux calibration\n"    
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

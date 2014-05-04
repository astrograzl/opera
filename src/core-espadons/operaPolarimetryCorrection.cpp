/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaPolarimetryCorrectionCorrection
 Version: 1.0
 Description: Create polar corrected spectrum.
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
#include "core-espadons/operaPolarimetryCorrection.h"

#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaWavelength.h"				// for wavelength polynomial
#include "libraries/operaLib.h"						// systemf
#include "libraries/operaLibCommon.h"               // for anglecoord_t and timecoord_t
#include "libraries/operaHelio.h"					// for sexigesimal conversion

/*! \file operaPolarimetryCorrectionCorrection.cpp */

using namespace std;

/*! 
 * operaPolarimetryCorrectionCorrection
 * \author Doug Teeple
 * \brief Apply wavelength correction based on telluric lines to .p byproduct.
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

/*
 Column #3 gives the unnormalizationd Stokes parameters, in the pn.s and the
 pu.s files.
 
 To get a normalizationd Stokes parameter that goes between -1.0 and +1.0 (or
 -100% and +100%), from a pn.s file or a pu.s file, divide column #3 by
 column #2.
 
 EXAMPLE:
 
 1609679pn.s
 369.0888  2.0876e-01  6.8397e-04 -8.1249e-03  6.7757e-03  4.4399e-03
 369.0911  2.0598e-01  5.0268e-03  4.3682e-03 -1.9970e-03  3.9532e-03
 369.0935  2.0679e-01 -3.9196e-05 -7.5297e-04 -1.5541e-03  4.2366e-03
 
 1609679pu.s
 369.0888  1.5895e+00  5.2079e-03 -6.1864e-02  5.1592e-02  3.3806e-02
 369.0911  1.5683e+00  3.8274e-02  3.3259e-02 -1.5205e-02  3.0099e-02
 369.0935  1.5745e+00 -2.9844e-04 -5.7330e-03 -1.1833e-02  3.2257e-02
 
 If you try the numbers, you get the same results when you divide col #3 by
 col #2.
 
 The 'n' and the 'u' refer to the intensity (if it's normalizationd to 1.0 or
 not). But in all cases, the Stokes given is not normalizationd to anything.
 The user has to do the normalization in both cases.
 
 If the data are reduced WITHOUT the continuum polarization subtracted
 (not recommended, rarely used), Upena gives the normalizationd Stokes
 parameter. So column #3 goes between -1.0 and +1.0 (or -100% and +100%).
 
 */


#define NOTPROVIDED -999

int main(int argc, char *argv[])
{
	int opt;
    
	string polar; 
	string input; 
	string object;
	string outputSpectraFile;
	string wavelengthCalibration;
	string rvel;
	string flatfieldspectrum;
	operaSpectralOrder_t spectralOrderType = LibreEspritsp1Spectrum;
	
    /*
     * Parameters for normalization
     */
    unsigned normalization = 0;
	bool fluxcalibrate = false;
	bool LEFluxCalibrate = false;
    unsigned normalizationBinsize = 100;
    bool usePolynomial = FALSE;
    unsigned orderOfPolynomial = 5;    
    
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
    string LEfluxCalibration;
	float exposureTime = 0.0;
    
    bool AbsoluteCalibration = false;
    int orderBin = 2;
    
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
		{"polar",						1, NULL, 'R'},	// .p
		{"outputCalibratedSpectrum",	1, NULL, 's'},  // .s                                                    
		{"spectrumtype",				1, NULL, 'y'},	// spectrum type
		{"wavelengthCalibration",		1, NULL, 'w'},	// wavelength calibration file (.wcal or .auto)
 		{"radialvelocitycorrection",	1, NULL, 'V'},  // Barycentric wavelength correction file (.rvel)
 		{"telluriccorrection",			1, NULL, 'T'},  // Telluric wavelength correction file (.tell)
 		{"flatfieldspectrum",			1, NULL, 'm'},  // flat field spectrum ff_
        
		{"normalize",					1, NULL, 'N'},	// apply flux normalization                                               
 		{"normalizationBinsize",		1, NULL, 'b'},	// binsize for normalization                                                   
		{"usePolynomial",				1, NULL, 'l'},	// option to use polynomial for normalization 
		{"orderOfPolynomial",			1, NULL, 'r'},	// option to set degree of polynomial for normalization
		
        {"fluxcalibrate",				1, NULL, 'f'},	// apply flux calibration        
        {"fluxCalibration",				1, NULL, 'C'},	// flux calibration file (.fcal)        
        {"LEfluxcalibrate",				1, NULL, 'a'},	// apply flux calibration        
        {"LEfluxCalibration",			1, NULL, 'x'},	// flux calibration; file (.fcal)        
        {"etime",						1, NULL, 'E'},	// needed for flux calibration        
		{"orderBin",                    1, NULL, 'B'},  // needed for flux calibration
		{"AbsoluteCalibration",         1, NULL, 'A'},  // absolute or relative flux calibration

		{"object",						1, NULL, 'o'},	// needed for Libre-Esprit output
		
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
	
	while((opt = getopt_long(argc, argv, "i:o:s:y:w:T:m:V:N:b:l:r:f:a:x:C:E:B:A:O:M:X:P:F:c:S:I:p::v::d::t::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				input = optarg;
				break;   
			case 'R':
				polar = optarg;
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
			case 'T':
				telluriccorrection = optarg;
				break;
			case 'm':
				flatfieldspectrum = optarg;
				break;
			case 'V':		// for telluric wl correction
				radialvelocitycorrection = optarg;
				break;
                
			case 'f':		// for flux calibration
				fluxcalibrate = atoi(optarg)==1;
				break;
			case 'a':		// for flux calibration
				LEFluxCalibrate = atoi(optarg)==1;
				break;
			case 'N':		// for normalization
				normalization = atoi(optarg);
				break;
			case 'b':		// normalization binsize
				normalizationBinsize = atoi(optarg);
				break;
			case 'l':		
				usePolynomial = atoi(optarg)==1;
				break;
			case 'r':
				orderOfPolynomial = atoi(optarg);
				break;  
                
			case 'C':       // for flux calibration
				fluxCalibration = optarg;
				break;    
			case 'x':       // for flux calibration
				LEfluxCalibration = optarg;
				break;    
			case 'E':
				exposureTime = atof(optarg);
				break;		
			case 'B':
				orderBin = atoi(optarg);
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
			throw operaException("operaPolarimetryCorrection: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (outputSpectraFile.empty()) {
			throw operaException("operaPolarimetryCorrection: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		if (wavelengthCalibration.empty()) {
			throw operaException("operaPolarimetryCorrection: wcal: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cout << "operaPolarimetryCorrection: input spectrum = " << input << endl; 
			cout << "operaPolarimetryCorrection: polar = " << polar << endl; 
			cout << "operaPolarimetryCorrection: object = " << object << endl; 
			cout << "operaPolarimetryCorrection: output spectrum file = " << outputSpectraFile << endl;
			cout << "operaPolarimetryCorrection: spectrum type = " << spectralOrderType << endl;							
			cout << "operaPolarimetryCorrection: wavelength calibration file = " << wavelengthCalibration << endl;
            cout << "operaPolarimetryCorrection: binsize for normalization = " << normalizationBinsize << endl;  
            cout << "operaPolarimetryCorrection: input flux calibration file = " << fluxCalibration << endl; 
            cout << "operaPolarimetryCorrection: flux calibrate = " << fluxcalibrate << endl; 
            cout << "operaPolarimetryCorrection: LE flux calibrate = " << LEFluxCalibrate << endl; 
            cout << "operaPolarimetryCorrection: LE input flux calibration file = " << LEfluxCalibration << endl; 
			cout << "operaPolarimetryCorrection: exposure time = " << exposureTime << endl;
			cout << "operaPolarimetryCorrection: order bin = " << orderBin << endl;
			cout << "operaPolarimetryCorrection: absolute calibration = " << AbsoluteCalibration << endl;            
			cout << "operaPolarimetryCorrection: radialvelocitycorrection = " << radialvelocitycorrection << endl;  
			cout << "operaPolarimetryCorrection: telluriccorrection = " << telluriccorrection << endl;
            cout << "operaPolarimetryCorrection: flatfieldspectrum = " << flatfieldspectrum << endl;
          
            if (ordernumber != NOTPROVIDED) {
                cout << "operaPolarimetryCorrection: ordernumber = " << ordernumber << endl;            
            }             
            
            if (plot) {
                cout << "operaPolarimetryCorrection: plotfilename = " << plotfilename << endl;
                cout << "operaPolarimetryCorrection: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaPolarimetryCorrection: continuumDataFilename = " << continuumDataFilename << endl;                
                cout << "operaPolarimetryCorrection: scriptfilename = " << scriptfilename << endl; 
                if (interactive) {
                    cout << "operaPolarimetryCorrection: interactive = YES" << endl; 
                } else {
                    cout << "operaPolarimetryCorrection: interactive = NO" << endl; 
                }
            }            
		}
        
		/*
		 * Plotting support
		 */
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
		operaSpectralOrderVector spectralOrders(polar);
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
		for (int order=minorder; order<=maxorder; order++) {			
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			if (spectralOrder->gethasSpectralElements()) {
				spectralOrder->getSpectralElements()->CreateExtendedvectors(spectralOrder->getSpectralElements()->getnSpectralElements());
			}
		}
		if (spectralOrderType == CalibratedExtendedBeamSpectrum) {
			for (int order=minorder; order<=maxorder; order++) {			
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					spectralOrder->getSpectralElements()->copyTOrawFlux();	// Save the raw flux for later
				}
			}
		}
		if (!telluriccorrection.empty()) {
			spectralOrders.ReadSpectralOrders(telluriccorrection);
			for (int order=minorder; order<=maxorder; order++) {			
				operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
				if (spectralOrder->gethasSpectralElements()) {
					operaWavelength *wavelength = spectralOrder->getWavelength();  
					operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
					spectralElements->setwavelengthsFromCalibration(wavelength);
					spectralElements->copyTOtell();	// Save the tell 
				}
			}
		} 
		spectralOrders.ReadSpectralOrders(wavelengthCalibration);
		
		if (!radialvelocitycorrection.empty()) {        
			spectralOrders.readRadialVelocityCorrection(radialvelocitycorrection);
		}
        /*
         * Flux Calibration Stuff ...
         */
        unsigned NumberofBeams = 0;
        double uncalibratedContinuumFluxForNormalization = 0;
        double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS];
        double spectralBinConstant = 0;
        double BeamSpectralBinConstant[MAXNUMBEROFBEAMS];
        
		if (!fluxCalibration.empty()) {
            spectralOrders.ReadSpectralOrders(fluxCalibration);
            
            for (int order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                if (spectralOrder->gethasSpectralElements()) {
                    NumberofBeams = spectralOrder->getnumberOfBeams();
                    break;
                }
            }
            
            unsigned nsigcut = 3;
            
            spectralOrders.getContinuumFluxesForNormalization(&uncalibratedContinuumFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization,normalizationBinsize, orderBin, nsigcut);
            
            spectralBinConstant = exposureTime;
            for(unsigned beam=0; beam < NumberofBeams; beam++) {
                BeamSpectralBinConstant[beam] = exposureTime;
            }
        }
        
		operaSpectralElements *fluxCalibrationElements = NULL;
		if (LEFluxCalibrate && !LEfluxCalibration.empty()) {
			unsigned length = spectralOrders.getLEElementCount(LEfluxCalibration);
			fluxCalibrationElements = new operaSpectralElements(length);
			spectralOrders.readLEFluxCalibration(LEfluxCalibration, fluxCalibrationElements);
		}
		
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
		
        if (verbose)
			cout << "operaPolarimetryCorrection: minorder ="<< minorder << " maxorder=" << maxorder << endl;    
        
        if(NumberofBeams) {
            NumberofBeams = 0;    // for plotting
        }
        //unsigned numberOfprintouts = 2; // for 3D plotting
        
		for (int order=minorder; order<=maxorder; order++) {			
			
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
			
            if (spectralOrder->gethasPolarimetry() && spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralElements()) {
				
                operaWavelength *wavelength = spectralOrder->getWavelength();  
				operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();  
                
				if (!fluxCalibration.empty() && spectralOrder->gethasSpectralEnergyDistribution()) {
					spectralOrder->applyFluxCalibration(spectralBinConstant, BeamSpectralBinConstant, uncalibratedContinuumFluxForNormalization,uncalibratedContinuumBeamFluxForNormalization, AbsoluteCalibration, fspecdata);
					spectralElements->copyTOfcalFlux();
					spectralElements->copyFROMrawFlux();
				}
				
				//spectralOrder->applyNormalization(normalizationBinsize,orderOfPolynomial,usePolynomial,fspecdata,fcontinuumdata, TRUE, numberOfprintouts);
				spectralElements->copyTOnormalizedFlux();
				spectralElements->copyFROMrawFlux();

				spectralElements->setwavelengthsFromCalibration(wavelength);
				spectralElements->setHasWavelength(true);   
				spectralElements->setHasDistance(false);
				
				if (spectralOrderType != CalibratedExtendedBeamSpectrum) {
					if (normalization) {
						spectralElements->copyFROMnormalizedFlux();
					} else if (!fluxCalibration.empty() && spectralOrder->gethasSpectralEnergyDistribution()) {
						spectralElements->copyFROMfcalFlux();
					}
					if (!telluriccorrection.empty()) {
						spectralElements->copyFROMtell();
					}
				}
				
				if (!normalization && LEFluxCalibrate && !LEfluxCalibration.empty()) {
					spectralOrder->applyFlatResponse(exposureTime, fluxCalibrationElements, fspecdata);
				}
                if (!radialvelocitycorrection.empty()) {
					double radialVelocityInMetersPerSecond = spectralOrders.getBarycentricRadialVelocityCorrection() * 1000.0;
					if (spectralOrderType == CalibratedExtendedBeamSpectrum) {
						unsigned elements = spectralElements->getnSpectralElements();
						while (elements--) {
							double currentElemWavelength = spectralElements->getwavelength(elements);
							double wavelengthCorrection = currentElemWavelength*(radialVelocityInMetersPerSecond/SPEED_OF_LIGHT_M);
							spectralElements->setrvel(wavelengthCorrection, elements);
						}
					} else {
						spectralOrder->applyBarycentricWavelengthCorrection(radialVelocityInMetersPerSecond);
					}
				}
			}
 		}
		// output a wavelength calibrated spectrum...
		spectralOrders.setObject(object);
		spectralOrders.WriteSpectralOrders(outputSpectraFile, spectralOrderType);
        
		if (fluxCalibrationElements != NULL)
			delete fluxCalibrationElements;
		
		if (fspecdata != NULL) {
			fspecdata->close();
            if (!plotfilename.empty() && !scriptfilename.empty()) {
                GenerateExtractionPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(), NumberofBeams, interactive);
            }
        } 
	}
	catch (operaException e) {
		cerr << "operaPolarimetryCorrection: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaPolarimetryCorrection: " << operaStrError(errno) << endl;
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
	"  -N, --normalization=1|0, apply flux normalization\n"
	"  -b, --normalizationBinsize=<float>,  binsize for normalization\n"
	"  -B, --orderBin=<int>, number or orders to bin for continuum evaluation\n"
	"  -A, --AbsoluteCalibration=<bool>, perform absolute flux calibration\n"    
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

/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateFlatFieldFluxCalibration
 Version: 1.0
 Description: Flux Calibration with Standard source 
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
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"						// for operaFitSpline
#include "libraries/ladfit.h"						// for ladfit
#include "libraries/Polynomial.h"
#include "libraries/operaFFT.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaStats.h"					// for operaArrayMedian
#include "libraries/operaSpectralTools.h"			//

#include "core-espadons/operaCreateFlatFieldFluxCalibration.h"

#define NOTPROVIDED -999

/*! \file operaCreateFlatFieldFluxCalibration.cpp */

using namespace std;

double getMedianValueInRange(operaSpectralElements *inputElements, unsigned centerElement, unsigned binsize);

int debug=0, verbose=0, trace=0, plot=0;

/*! 
 * operaCreateFlatFieldFluxCalibration
 * \author Eder Martioli
 * \brief Flux Calibration with Standard source.
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
	
	string inputMasterFlatSpectrum;     // extracted masterflat spectrum (*.e.gz file)
	string outputFluxCalibrationFile;   // fluxCalibration (*.fcal.gz file)
	string wavelengthCalibration;

    double wavelengthForNormalization = 0.0;

	int ordernumber = NOTPROVIDED;
    
    unsigned binsize = 30;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false;            
    
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string spectrumDataFilename;
	string continuumDataFilename;	
	string scriptfilename;	
	
	struct option longopts[] = {
		{"inputMasterFlatSpectrum",         1, NULL, 'i'},
		{"outputFluxCalibrationFile",       1, NULL, 'o'},
		{"wavelengthCalibration",           1, NULL, 'w'},
		{"wavelengthForNormalization",      1, NULL, 'W'},
        {"binsize",                         1, NULL, 'b'},
		{"ordernumber",                     1, NULL, 'O'},
		{"minorder",                        1, NULL, 'M'},
		{"maxorder",                        1, NULL, 'X'},
		{"plotfilename",                    1, NULL, 'P'},
		{"spectrumDataFilename",            1, NULL, 'F'},
		{"continuumDataFilename",           1, NULL, 'C'},        
		{"scriptfilename",                  1, NULL, 'S'},  
		{"interactive",                     0, NULL, 'I'},
		{"plot",				optional_argument, NULL, 'p'},       
		{"verbose",				optional_argument, NULL, 'v'},
		{"debug",				optional_argument, NULL, 'd'},
		{"trace",				optional_argument, NULL, 't'},
		{"help",				no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:o:w:W:b:O:M:X:P:F:C:S:I:p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputMasterFlatSpectrum = optarg;	
				break;    
			case 'o':		// output
				outputFluxCalibrationFile = optarg;
				break;
			case 'w':
				wavelengthCalibration = optarg;
				break;
            case 'W':	
				wavelengthForNormalization = atof(optarg);
				break;
            case 'b':		// binsize
				binsize = atoi(optarg);
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
			case 'C':
				continuumDataFilename = optarg;
				break;                 
			case 'S':
				scriptfilename = optarg;
				break;  
			case 'I':		// for interactive plots
				interactive = true;
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
	
	try {
		// we need an input uncalibrated spectrum...
		if (inputMasterFlatSpectrum.empty()) {
			throw operaException("operaCreateFlatFieldFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an output...
		if (outputFluxCalibrationFile.empty()) {
			throw operaException("operaCreateFlatFieldFluxCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
        if (wavelengthCalibration.empty()) {
			throw operaException("operaCreateFlatFieldFluxCalibration: wcal: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}

		if (verbose) {
			cout << "operaCreateFlatFieldFluxCalibration: inputMasterFlatSpectrum = " << inputMasterFlatSpectrum << endl; 
            cout << "operaCreateFlatFieldFluxCalibration: outputFluxCalibrationFile = " << outputFluxCalibrationFile << endl;
            cout << "operaCreateFlatFieldFluxCalibration: wavelengthCalibration = " << wavelengthCalibration << endl;
            cout << "operaCreateFlatFieldFluxCalibration: wavelengthForNormalization = " << wavelengthForNormalization << endl;
            cout << "operaCreateFlatFieldFluxCalibration: binsize = " << binsize << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaCreateFlatFieldFluxCalibration: ordernumber = " << ordernumber << endl;            
            }   
            if(plot) {
                cout << "operaCreateFlatFieldFluxCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaCreateFlatFieldFluxCalibration: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaCreateFlatFieldFluxCalibration: continuumDataFilename = " << continuumDataFilename << endl;
                cout << "operaCreateFlatFieldFluxCalibration: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cout << "operaCreateFlatFieldFluxCalibration: interactive = YES" << endl; 
                } else {
                    cout << "operaCreateFlatFieldFluxCalibration: interactive = NO" << endl; 
                }
            }            
            
		}
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
        
		operaSpectralOrderVector spectralOrders(inputMasterFlatSpectrum);
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
        
		if (verbose)
			cout << "operaCreateFlatFieldFluxCalibration: minorder ="<< minorder << " maxorder=" << maxorder << endl;
  
        double maxFluxForNormalization = 0.0; //  DT May 08 2014 
        double maxBeamFluxForNormalization[MAXNUMBEROFBEAMS];
        unsigned numberOfBeams = 0;
        unsigned elemIndexPicked = 0;
        int refOrder = -999;
        
        // First loop over orders to search for order with maximum number of
        // elements
        // --> maxNElements
        unsigned maxNElements = 0;
        for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);

            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                operaWavelength *wavelength = spectralOrder->getWavelength();
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                SpectralElements->setwavelengthsFromCalibration(wavelength);
                if(maxNElements < SpectralElements->getnSpectralElements()) {
                    maxNElements = SpectralElements->getnSpectralElements();
                }
            }
        }
  
        // Calculate reference order and maximum flux for normalization
        // --> refOrder, maxFluxForNormalization
        if(wavelengthForNormalization) {
            // If wavelength for normalization is given, then use max flux at wavelength
            // --> refOrder, maxFluxForNormalization
            int *orderWithReferenceFluxForNormalization = new int[MAXORDERS];
            unsigned *elemIndexWithReferenceFluxForNormalization = new unsigned[MAXORDERS];
            unsigned nOrdersPicked = spectralOrders.getElemIndexAndOrdersByWavelength(orderWithReferenceFluxForNormalization,elemIndexWithReferenceFluxForNormalization,wavelengthForNormalization);

            for(unsigned i=0;i<nOrdersPicked;i++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(orderWithReferenceFluxForNormalization[i]);

                if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                    
                    //May 08 2014 -- not used -- operaWavelength *wavelength = spectralOrder->getWavelength();
                    operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                    double tmp_FluxForNormalization = SpectralElements->getFlux(elemIndexWithReferenceFluxForNormalization[i]);
                    
                    if(maxFluxForNormalization < tmp_FluxForNormalization) {
                        maxFluxForNormalization = tmp_FluxForNormalization;
                        refOrder = orderWithReferenceFluxForNormalization[i];
                        elemIndexPicked = elemIndexWithReferenceFluxForNormalization[i];
                        numberOfBeams = spectralOrder->getnumberOfBeams();
                        
                        for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                            operaSpectralElements *BeamElements = spectralOrder->getBeamElements(beam);
                            maxBeamFluxForNormalization[beam] = BeamElements->getFlux(elemIndexPicked);
                        }   
                    }
                }
            }
            
        } else {
            float *flux_tmp = new float[maxNElements];

            for (int order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                
                if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                    operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();

                    unsigned NumberofElementsToBin = binsize; 
                    unsigned NumberOfElementSamples = (unsigned)ceil((float)SpectralElements->getnSpectralElements()/(float)NumberofElementsToBin);
                    
                    for(unsigned k=0;k<NumberOfElementSamples;k++){
                        unsigned firstElement = NumberofElementsToBin*(k);
                        unsigned lastElement =  NumberofElementsToBin*(k+1);
                        if (lastElement > SpectralElements->getnSpectralElements()){
                            lastElement = SpectralElements->getnSpectralElements();   
                        }
                        unsigned np=0;
                        for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
                            flux_tmp[np++] = (float)SpectralElements->getFlux(indexElem);
                        }
                        double medianFlux = (double)operaArrayMedian(np,flux_tmp);
                        
                        if(maxFluxForNormalization < medianFlux) {
                            maxFluxForNormalization = medianFlux;
                            refOrder = order;
                            elemIndexPicked = firstElement + (lastElement - firstElement)/2;
                            numberOfBeams = spectralOrder->getnumberOfBeams();
                            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                                operaSpectralElements *BeamElements = spectralOrder->getBeamElements(beam);
                                np=0;
                                for(unsigned indexElem=firstElement;indexElem < lastElement; indexElem++) {
                                    flux_tmp[np++] = (float)BeamElements->getFlux(indexElem);
                                }
                                maxBeamFluxForNormalization[beam] = (double)operaArrayMedian(np,flux_tmp);
                            }
                        }
                    }
                }
            }
        }
 
        if(debug) {
            cout << "refOrder=" << refOrder << endl;
            cout << "maxFluxForNormalization=" << maxFluxForNormalization << endl;
            for(unsigned beam = 0; beam < numberOfBeams; beam++) {
                cout << "maxBeamFluxForNormalization[" << beam << "]=" << maxBeamFluxForNormalization[beam] << endl;
            }
        }
       
        // Normalize flat-field spectrum by maxFluxForNormalization
        // to generate flat-field flux calibration spectrum
        // --> spectralEnergyDistribution
		for (int order=minorder; order<=maxorder; order++) {
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);

            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                if (verbose) {
					cout << "operaCreateFlatFieldFluxCalibration: Calibrating order number: "<< order << endl;
                }
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                operaWavelength *wavelength = spectralOrder->getWavelength();

                // use info in SpectralElements to create spectralEnergyDistribution
                spectralOrder->createSpectralEnergyDistributionElements(SpectralElements->getnSpectralElements());
                
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();                
                spectralEnergyDistribution->setUncalibratedFluxElements(SpectralElements);
                spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralElements *BeamElements = spectralOrder->getBeamElements(beam);
                    for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
                        BeamElements->setdistd(SpectralElements->getdistd(indexElem), indexElem);
                    }
                    BeamElements->setwavelengthsFromCalibration(wavelength);
                    
                    operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
                    BeamSED->setUncalibratedFluxElements(BeamElements);
                    BeamSED->setwavelengthForNormalization(wavelengthForNormalization);
                }
                
                // normalize spectral Elements
                spectralOrder->normalizeSpectralElementsByConstant(maxFluxForNormalization,maxBeamFluxForNormalization);

                operaSpectralElements *fluxCalibrationElements = spectralEnergyDistribution->getFluxCalibrationElements();
                operaSpectralElements *thruputElements = spectralEnergyDistribution->getThroughputElements();
                
                for(unsigned indexElem=0;indexElem < SpectralElements->getnSpectralElements(); indexElem++) {
                    fluxCalibrationElements->setFlux(SpectralElements->getFlux(indexElem),indexElem);
                    thruputElements->setFlux(SpectralElements->getFlux(indexElem),indexElem);
                    fluxCalibrationElements->setwavelength(SpectralElements->getwavelength(indexElem),indexElem);
                    thruputElements->setwavelength(SpectralElements->getwavelength(indexElem),indexElem);
                    
                    for(unsigned beam=0; beam < numberOfBeams; beam++) {
                        operaSpectralElements *BeamElements = spectralOrder->getBeamElements(beam);
                        operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
                        operaSpectralElements *fluxCalibrationBeamElements = BeamSED->getFluxCalibrationElements();
                        operaSpectralElements *thruputBeamElements = BeamSED->getThroughputElements();

                        fluxCalibrationBeamElements->setFlux(BeamElements->getFlux(indexElem),indexElem);
                        thruputBeamElements->setFlux(BeamElements->getFlux(indexElem),indexElem);
                        fluxCalibrationBeamElements->setwavelength(SpectralElements->getwavelength(indexElem),indexElem);
                        thruputBeamElements->setwavelength(SpectralElements->getwavelength(indexElem),indexElem);
                    }
                }
                
                spectralEnergyDistribution->setHasUncalibratedFlux(true);
                spectralEnergyDistribution->setHasFluxCalibration(true);
                spectralEnergyDistribution->setHasInstrumentThroughput(true);
                for(unsigned beam=0; beam < numberOfBeams; beam++) {
                    operaSpectralEnergyDistribution *BeamSED = spectralOrder->getBeamSED(beam);
                    BeamSED->setHasUncalibratedFlux(true);
                    BeamSED->setHasFluxCalibration(true);
                    BeamSED->setHasInstrumentThroughput(true);
                }
                spectralOrder->sethasSpectralEnergyDistribution(true);
                
                if(debug) {
                    for(unsigned i=0; i<SpectralElements->getnSpectralElements(); i++) {
                        cout << order << " " << i << " "
                        << SpectralElements->getwavelength(i) << " "
                        << SpectralElements->getdistd(i) << " "
                        << SpectralElements->getFlux(i) << " "
                        << spectralEnergyDistribution->getUncalibratedFluxElements()->getFlux(i) << " "
                        << fluxCalibrationElements->getFlux(i) << " "
                        << thruputElements->getFlux(i) << endl;
                    }
               }

			} //if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
            
            if (!spectralOrder->gethasSpectralElements() ||  !spectralOrder->gethasWavelength()) {
				if (verbose)
					cout << "operaCreateFlatFieldFluxCalibration: Skipping order number: "<< order << " no spectral elements." << endl;
			}
        }

        /*
		 * and write output
		 */
		spectralOrders.WriteSpectralOrders(outputFluxCalibrationFile, Fcal);
		
        if (fspecdata != NULL && fcontinuumdata != NULL) {
			fspecdata->close();
            fcontinuumdata->close();
            
            if (!scriptfilename.empty()) {
                GenerateCreateFlatFieldFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(),continuumDataFilename.c_str(), 2, interactive);
            }
        }
        
	}
	catch (operaException e) {
		cerr << "operaCreateFlatFieldFluxCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateFlatFieldFluxCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --inputMasterFlatSpectrum=<SPEC_FILE>"
	" --outputFluxCalibrationFile=<SPEC_FILE>"
	" --wavelengthCalibration=<WCAL_FILE>"
	" --wavelengthForNormalization=<DBL_VALUE>"
	" --binsize=<UNS_VALUE>"
	" --ordernumber=<UNS_VALUE>"
	" --minorder=<UNS_VALUE>"
	" --maxorder=<UNS_VALUE>"
	" --plotContinuum=<BOOL>"
	" --plotfilename=<EPS_FILE>"
	" --spectrumDataFilename=<DATA_FILE>"
	" --continuumDataFilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputMasterFlatSpectrum=/Users/edermartioli/opera//calibrations/FCAL001-13BQ04-Sep12/ff_OLAPAa_pol_Normal.e.gz --wavelengthCalibration=/Users/edermartioli/opera//calibrations/FCAL001-13BQ04-Sep12/OLAPAa_pol_Normal.wcar.gz  --outputFluxCalibrationFile=/Users/edermartioli/opera//calibrations/FCAL001-13BQ04-Sep12/ff_OLAPAa_pol_Normal.fcal.gz  --binsize=50 -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -i, --inputMasterFlatSpectrum=<SPEC_FILE>,  Spectrophotometric standard extracted uncalibrated spectrum \n"
	"  -o, --outputFluxCalibrationFile=<SPEC_FILE>,  Output flux calibration conversion file \n"
	"  -w, --wavelengthCalibration=<WCAL_FILE>,  Input wavelength calibration file \n"
    "  -W, --wavelengthForNormalization=<DBL_VALUE>, Choose a wavlength for which the spectrum should be normalized \n"
	"  -b, --binsize=<UNS_VALUE>, Number of points to bin for continuum estimate \n"
	"  -O, --ordernumber=<UNS_VALUE>, Absolute order number to extract (default=all)\n"
	"  -M, --minorder=<UNS_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<UNS_VALUE>, Define maximum order number\n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -c, --plotContinuum=<BOOL>, Switch to generate plot of continuum or normalized line spectra\n"
	"  -F, --spectrumDataFilename=<DATA_FILE>\n"
	"  -C, --continuumDataFilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateCreateFlatFieldFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display)
{
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str());  // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    *fgnu << "\nset xlabel \"wavelength (nm)\"" << endl;
    *fgnu << "set ylabel \"flux\"" << endl;
    
    *fgnu << "set pointsize 1.0" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
        if (display) {
            *fgnu << "\nset terminal x11" << endl;
            *fgnu << "set output" << endl;
            *fgnu << "replot" << endl;
        } else {
            *fgnu << "\n#set terminal x11" << endl;
            *fgnu << "#set output" << endl;
            *fgnu << "#replot" << endl;
        }
    } else {
        *fgnu << "\nplot \"" << spectrumDataFilename << "\" u 6:7 w d" <<
        ",\"" << continuumDataFilename << "\" u 4:5 w linespoint lw 2.5" << endl;
        
        *fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        *fgnu << "#replot" << endl;
        *fgnu << "#set terminal x11" << endl;
        *fgnu << "#set output" << endl;
    }
    
    fgnu->close();
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());
    }
}

double getMedianValueInRange(operaSpectralElements *inputElements, unsigned centerElement, unsigned binsize) {
    
    unsigned nElements = inputElements->getnSpectralElements();
    
    int i1 = (int)centerElement - (int)binsize/2;
    int i2 = (int)centerElement + (int)binsize/2;
    
    if(i1 < 0) {
        i1 = 0;
        i2 = i1 + (int)binsize;
    }
    if((unsigned)i2 > nElements) {
        i2 = (int)nElements;
        i1 = (int)nElements - (int)binsize;
        if(i1 < 0) {
            i1 = 0;
        }
    }

    float *fluxInRange = new float[nElements];
    unsigned count = 0;
    for(unsigned i=(unsigned)i1;i<(unsigned)i2;i++) {
        fluxInRange[count++] = (float)inputElements->getFlux(i);
    }
    
    double medianFlux = (double)operaArrayMedian(count,fluxInRange);
    delete[] fluxInRange;
    
    return medianFlux;
}

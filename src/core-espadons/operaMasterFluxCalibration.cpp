/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaMasterFluxCalibration
 Version: 1.0
 Description: Create a Master Flux Calibration out of many calibrations. 
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

#include <fstream>
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

#define MAXNUMBEROFFLUXCALIBRATIONFILES 500
#define WAVELENGTH_PRECISION 0.01
#define MINELEMENTS 20

/*! \file operaMasterFluxCalibration.cpp */

using namespace std;

void GenerateMasterFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string outputDataFilename, unsigned NumberofBeams, bool display);
        
/*! 
 * operaMasterFluxCalibration
 * \author Eder Martioli
 * \brief Create a Master Flux Calibration out of many calibrations.
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
	operaArgumentHandler args;
	
	string listofinputfcal;
	double inputconstant = 0.0;
    string outputfcal;
	string inputReferenceSpectrum;      // this is a reference opera spectrum that defines the format for the output
    string inputWaveFile;
    unsigned combineMethod = 1; // Not implemented yet. Eder - Jan/09/2013.
	int ordernumber = NOTPROVIDED;
    int minorder = 22;
    int maxorder = 62;    
    string plotfilename;	
	string spectrumDataFilename;
	string outputDataFilename;
	string scriptfilename;
	bool interactive = false;
	
	args.AddOptionalArgument("inputfcal", listofinputfcal, "", "List of input flux calibration files, separated by spaces");
	args.AddOptionalArgument("inputconstant", inputconstant, 0.0, "This value should be input when there are no fcal files available");
	args.AddRequiredArgument("outputfcal", outputfcal, "Output master flux calibration file");
	args.AddRequiredArgument("inputReferenceSpectrum", inputReferenceSpectrum, "Input reference spectrum file");
	args.AddRequiredArgument("inputWaveFile", inputWaveFile, "Input wavelength calibration file");
	args.AddOptionalArgument("combineMethod", combineMethod, 1, "Method for combining images: 1 = Median, 2 = Mean, 3 = Weighted mean -- WARNING: NOT CURRENTLY IMPLEMENTED"); //NOT IMPLEMENTED YET
	args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
	args.AddOptionalArgument("plotfilename", plotfilename, "", "Output plot eps file name");
	args.AddOptionalArgument("spectrumDataFilename", spectrumDataFilename, "", "Output spectrum data file name");
	args.AddOptionalArgument("outputDataFilename", outputDataFilename, "", "Output data file name");
	args.AddOptionalArgument("scriptfilename", scriptfilename, "", "Output gnuplot script file name");
	args.AddSwitch("interactive", interactive, "For interactive plots");
	
	try {
        args.Parse(argc, argv);
        
        string inputfcal[MAXNUMBEROFFLUXCALIBRATIONFILES];
		unsigned inputFcalIndex = 0;
		SplitStringIntoArray(listofinputfcal, inputfcal, inputFcalIndex, MAXNUMBEROFFLUXCALIBRATIONFILES); // Split list of images into array
        
		// we need at least one input opera flux calibration file..
        if(inputFcalIndex == 0 && inputconstant == 0.0) {
			throw operaException("operaMasterFluxCalibration: no input fcal or constant. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        } else if (inputFcalIndex > MAXNUMBEROFFLUXCALIBRATIONFILES) {
            throw operaException("operaMasterFluxCalibration: number of input fcal exceeds limit. ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }
		// we need an output calibration file name...
		if (outputfcal.empty()) {
			throw operaException("operaMasterFluxCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input reference spectrum...
		if (inputReferenceSpectrum.empty()) {
			throw operaException("operaMasterFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (args.verbose) {
            for(unsigned index=0; index<inputFcalIndex; index++) {
                cout << "operaMasterFluxCalibration: input flux calibration inputfcal["<<index<<"] = " << inputfcal[index] << endl;
            }
            cout << "operaMasterFluxCalibration: inputconstant = " << inputconstant << endl;
			cout << "operaMasterFluxCalibration: output Flux Calibration .fcal file = " << outputfcal << endl;
			cout << "operaMasterFluxCalibration: inputReferenceSpectrum = " << inputReferenceSpectrum << endl;
			cout << "operaMasterFluxCalibration: inputWaveFile = " << inputWaveFile << endl;
			cout << "operaMasterFluxCalibration: combineMethod = " << combineMethod << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaMasterFluxCalibration: ordernumber = " << ordernumber << endl;
            if(args.plot) {
                cout << "operaMasterFluxCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaMasterFluxCalibration: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaMasterFluxCalibration: outputDataFilename = " << spectrumDataFilename << endl;
                cout << "operaMasterFluxCalibration: scriptfilename = " << scriptfilename << endl;
			    cout << "operaMasterFluxCalibration: interactive = " << (interactive ? "YES" : "NO") << endl;
            }
		}
        
        ofstream fspecdata;
        ofstream foutdata;
        if (!spectrumDataFilename.empty()) fspecdata.open(spectrumDataFilename.c_str());
        if (!outputDataFilename.empty()) foutdata.open(outputDataFilename.c_str());
        
		operaSpectralOrderVector spectralOrders(inputReferenceSpectrum);
        spectralOrders.ReadSpectralOrders(inputWaveFile);
        
        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
		if (args.verbose) cout << "operaMasterFluxCalibration: minorder ="<< minorder << " maxorder=" << maxorder << endl;
        
        unsigned NumberofBeams = 0;
        double wavelengthForNormalization = 0;
        
        // initialize output vectors
        for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
			if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                unsigned nElements = SpectralElements->getnSpectralElements();
                
                if(nElements==0) {
                    cerr << "operaMasterFluxCalibration: skipping order " << order << " -> reference spectrum has nElements=0" << endl;
                    continue;
                }
                
                SpectralElements->setwavelengthsFromCalibration(spectralOrder->getWavelength());
                
                if(NumberofBeams == 0) {
                    NumberofBeams = spectralOrder->getnumberOfBeams();
                }
                // create output energy distribution elements for output fcals
                spectralOrder->createSpectralEnergyDistributionElements(nElements);
                operaSpectralEnergyDistribution *SpectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                
                for(unsigned indexElem=0; indexElem<nElements; indexElem++) {
                    if(inputFcalIndex) {
                        SpectralEnergyDistribution->getFluxCalibrationElements()->setFlux(0.0,indexElem);
                        SpectralEnergyDistribution->getThroughputElements()->setFlux(0.0,indexElem);
                    } else {
                        SpectralEnergyDistribution->getFluxCalibrationElements()->setFlux(inputconstant,indexElem);
                        SpectralEnergyDistribution->getThroughputElements()->setFlux(inputconstant,indexElem);
                    }
                    SpectralEnergyDistribution->getFluxCalibrationElements()->setwavelength(SpectralElements->getwavelength(indexElem),indexElem);
                    SpectralEnergyDistribution->getFluxCalibrationElements()->setFluxVariance(0.0,indexElem);
                    SpectralEnergyDistribution->getThroughputElements()->setFluxVariance(0.0,indexElem);
                    
                    for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                        if(inputFcalIndex) {
                            spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->setFlux(0.0,indexElem);
                            spectralOrder->getBeamSED(beam)->getThroughputElements()->setFlux(0.0,indexElem);
                        } else {
                            spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->setFlux(inputconstant,indexElem);
                            spectralOrder->getBeamSED(beam)->getThroughputElements()->setFlux(inputconstant,indexElem);
                        }
                        spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->setFluxVariance(0.0,indexElem);
                        spectralOrder->getBeamSED(beam)->getThroughputElements()->setFluxVariance(0.0,indexElem);
                    }
                    
                }
                spectralOrder->sethasSpectralEnergyDistribution(true);
                SpectralEnergyDistribution->setHasFluxCalibration(true);
                SpectralEnergyDistribution->setHasInstrumentThroughput(true);
            }
        }
        
        float *referenceWavelength = new float[MAXSPECTRALELEMENTSPERORDER];
        float *inputWavelength = new float[MAXSPECTRALELEMENTSPERORDER];
        
        float *inputFluxCal = new float[MAXSPECTRALELEMENTSPERORDER];
        float *inputThroughput = new float[MAXSPECTRALELEMENTSPERORDER];
        float *referenceInputFluxCal = new float[MAXSPECTRALELEMENTSPERORDER];
        float *referenceInputThroughput = new float[MAXSPECTRALELEMENTSPERORDER];
        
        float *inputBeamFluxCal[MAXNUMBEROFBEAMS],*referenceBeamInputFluxCal[MAXNUMBEROFBEAMS];
        float *inputBeamThroughput[MAXNUMBEROFBEAMS],*referenceBeamInputThroughput[MAXNUMBEROFBEAMS];
        
        operaSpectralEnergyDistribution *BeamSED[MAXNUMBEROFBEAMS];
		operaSpectralElements *beamFluxcalibration[MAXNUMBEROFBEAMS];
        operaSpectralElements *beamThroughput[MAXNUMBEROFBEAMS];
        
        for(unsigned beam = 0; beam < NumberofBeams; beam++) {
            inputBeamFluxCal[beam] = new float[MAXSPECTRALELEMENTSPERORDER];
            inputBeamThroughput[beam] = new float[MAXSPECTRALELEMENTSPERORDER];
            referenceBeamInputFluxCal[beam] = new float[MAXSPECTRALELEMENTSPERORDER];
            referenceBeamInputThroughput[beam] = new float[MAXSPECTRALELEMENTSPERORDER];
        }
        
        for(unsigned index=0; index<inputFcalIndex; index++) {
            // read input flux calibration
            operaSpectralOrderVector inputSpectralOrders(inputfcal[index]);
            
            // get wavelengthForNormalization
            if(index==0) {
                for (int order=minorder; order<=maxorder; order++) {
                    operaSpectralOrder *inputSpectralOrder = inputSpectralOrders.GetSpectralOrder(order);
                    if (inputSpectralOrder->gethasSpectralEnergyDistribution()) {
                        operaSpectralEnergyDistribution *inputSpectralEnergyDistribution = inputSpectralOrder->getSpectralEnergyDistribution();
                        if(inputSpectralEnergyDistribution->getwavelengthForNormalization() != 0) {
                            wavelengthForNormalization = inputSpectralEnergyDistribution->getwavelengthForNormalization();
                            break;
                        }
                    }
                }
            }
            
            bool interpolate = FALSE;
            
            for (int order=minorder; order<=maxorder; order++) {
                
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                
                if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralEnergyDistribution()) {
                    
                    if (args.debug) cout << "operaMasterFluxCalibration: processing fcal index=" << index << " order=" << order << endl;
                    
                    operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();
                    operaSpectralEnergyDistribution *SpectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                    
                    SpectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                    
                    unsigned nElements = SpectralElements->getnSpectralElements();
                    
                    // DT May 20 2014 you can't do spline fitting without some reasobale number of elements... if(nElements==0) {
					if(nElements < MINELEMENTS) {
                        cerr << "operaMasterFluxCalibration: skipping order " << order << " -> reference spectrum has nElements=0" << endl;
                        continue;
                    }
                    
                    for(unsigned indexElem=0; indexElem<nElements; indexElem++) {
                        referenceWavelength[indexElem] = (float)SpectralElements->getwavelength(indexElem);
                    }
                    
                    // grab order of input fcal
                    operaSpectralOrder *inputSpectralOrder = inputSpectralOrders.GetSpectralOrder(order);

                    if (inputSpectralOrder->gethasSpectralEnergyDistribution()) {
                        operaSpectralEnergyDistribution *inputSpectralEnergyDistribution = inputSpectralOrder->getSpectralEnergyDistribution();

                        operaSpectralElements *FluxCalibration = inputSpectralEnergyDistribution->getFluxCalibrationElements();
                        operaSpectralElements *InstrumentThroughput = inputSpectralEnergyDistribution->getThroughputElements();
                        unsigned inputnElements = inputSpectralEnergyDistribution->getFluxCalibrationElements()->getnSpectralElements();
                        if(args.verbose) cout << "operaMasterFluxCalibration: processing fcal index=" << index << " order=" << order << " refnElements=" << nElements << " inputnElements=" << inputnElements << endl;
                        
                        if(inputnElements != nElements) {
                            interpolate = TRUE;
                            if(args.verbose) cout << "operaMasterFluxCalibration: using interpolation since inputnElements="<< inputnElements<<" != nElements=" << nElements << endl;
                        }
                        
                        for(unsigned beam = 0; beam < inputSpectralOrder->getnumberOfBeams(); beam++) {
                            BeamSED[beam] = inputSpectralOrder->getBeamSED(beam);
                            beamFluxcalibration[beam] = BeamSED[beam]->getFluxCalibrationElements();
                            beamThroughput[beam] = BeamSED[beam]->getThroughputElements();
                        }

                        // loop below collects data from input fcal
                        for(unsigned indexElem=0; indexElem<inputnElements; indexElem++) {
                            inputWavelength[indexElem] = (float)FluxCalibration->getwavelength(indexElem);

                            if(!interpolate) {
                                if(fabs(inputWavelength[indexElem] - referenceWavelength[indexElem]) > WAVELENGTH_PRECISION) {
                                    //if(inputWavelength[indexElem] != referenceWavelength[indexElem]) {
                                    interpolate = TRUE;
                                    if(args.verbose) cout << "operaMasterFluxCalibration: using interpolation since |inputWavelength - referenceWavelength| = " << fabs(inputWavelength[indexElem] - referenceWavelength[indexElem]) << " > " << WAVELENGTH_PRECISION << endl;
                                }
                            }

                            inputFluxCal[indexElem] = (float)FluxCalibration->getFlux(indexElem);
                            inputThroughput[indexElem] = (float)InstrumentThroughput->getFlux(indexElem);

                            for(unsigned beam = 0; beam < inputSpectralOrder->getnumberOfBeams(); beam++) {
                                inputBeamFluxCal[beam][indexElem] = (float)(beamFluxcalibration[beam]->getFlux(indexElem));
                                inputBeamThroughput[beam][indexElem] = (float)(beamThroughput[beam]->getFlux(indexElem));
                            }

                            // below is for plotting
                            if(fspecdata.is_open()) {
                                fspecdata << index << ' ' << order << ' ' << indexElem << ' '
                                << FluxCalibration->getwavelength(indexElem) << ' '
                                << FluxCalibration->getFlux(indexElem) << ' '
                                << FluxCalibration->getFluxVariance(indexElem) << ' '
                                << InstrumentThroughput->getFlux(indexElem) << ' '
                                << InstrumentThroughput->getFluxVariance(indexElem) << ' ';
                                for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                                    fspecdata << beam << ' '
                                    << beamFluxcalibration[beam]->getFlux(indexElem) << ' '
                                    << beamFluxcalibration[beam]->getFluxVariance(indexElem) << ' '
                                    << beamThroughput[beam]->getFlux(indexElem) << ' '
                                    << beamThroughput[beam]->getFluxVariance(indexElem) << ' ';
                                }
                                fspecdata << endl;
                            }
                        }
                        
                        // below is for plotting
                        if(fspecdata.is_open()) fspecdata << endl;
                                                
                        if(interpolate) {
                            if(args.verbose) cout << "operaMasterFluxCalibration: Starting interpolations" << endl;
                            
                            operaFitSpline(inputnElements,inputWavelength,inputFluxCal,nElements,referenceWavelength,referenceInputFluxCal);
                            operaFitSpline(inputnElements,inputWavelength,inputThroughput,nElements,referenceWavelength,referenceInputThroughput);
                            
                            for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                                operaFitSpline(inputnElements,inputWavelength,inputBeamFluxCal[beam],nElements,referenceWavelength,referenceBeamInputFluxCal[beam]);
                                operaFitSpline(inputnElements,inputWavelength,inputBeamThroughput[beam],nElements,referenceWavelength,referenceBeamInputThroughput[beam]);
                            }
                        } else {
                            for(unsigned indexElem=0; indexElem<inputnElements; indexElem++) {
                                referenceInputFluxCal[indexElem] = inputFluxCal[indexElem];
                                referenceInputThroughput[indexElem] = inputThroughput[indexElem];
                                
                                for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                                    referenceBeamInputFluxCal[beam][indexElem] = inputBeamFluxCal[beam][indexElem];
                                    referenceBeamInputThroughput[beam][indexElem] = inputBeamThroughput[beam][indexElem];
                                }
                            }
                        }

                        if(args.debug) cout << "operaMasterFluxCalibration: adding values to the output fcal" << endl;
                        
                        for(unsigned indexElem=0; indexElem<nElements; indexElem++) {
                            double fluxCal = SpectralEnergyDistribution->getFluxCalibrationElements()->getFlux(indexElem) + (double)referenceInputFluxCal[indexElem]/(double)inputFcalIndex;
                            double thru = SpectralEnergyDistribution->getThroughputElements()->getFlux(indexElem) + (double)referenceInputThroughput[indexElem]/(double)inputFcalIndex;
                            
                            SpectralEnergyDistribution->getFluxCalibrationElements()->setFlux(fluxCal,indexElem);
                            SpectralEnergyDistribution->getThroughputElements()->setFlux(thru,indexElem);
                            
                            for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                                double beamfluxCal = spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFlux(indexElem) + (double)referenceBeamInputFluxCal[beam][indexElem]/(double)inputFcalIndex;
                                double beamthru = spectralOrder->getBeamSED(beam)->getThroughputElements()->getFlux(indexElem) + (double)referenceBeamInputThroughput[beam][indexElem]/(double)inputFcalIndex;
                                
                                spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->setFlux(beamfluxCal,indexElem);
                                spectralOrder->getBeamSED(beam)->getThroughputElements()->setFlux(beamthru,indexElem);
                            }
                            
                            // below is for plotting
                            if(foutdata.is_open() && index==inputFcalIndex-1) {
                                foutdata << index << ' ' << order << ' ' << indexElem << ' '
                                << SpectralElements->getwavelength(indexElem) << ' '
                                << SpectralElements->getFlux(indexElem) << ' '
                                << SpectralElements->getFluxVariance(indexElem) << ' '
                                << SpectralEnergyDistribution->getFluxCalibrationElements()->getFlux(indexElem) << ' '
                                << SpectralEnergyDistribution->getFluxCalibrationElements()->getFluxVariance(indexElem) << ' ';
                                for(unsigned beam = 0; beam < spectralOrder->getnumberOfBeams(); beam++) {
                                    foutdata << beam << ' '
                                    << spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFlux(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getFluxCalibrationElements()->getFluxVariance(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getThroughputElements()->getFlux(indexElem) << ' '
                                    << spectralOrder->getBeamSED(beam)->getThroughputElements()->getFluxVariance(indexElem) << ' ';
                                }
                                foutdata << endl;
                            }
                        }
                        // below is for plotting
                        if(foutdata.is_open() && index==inputFcalIndex-1) {
                            foutdata << endl;
                        }
                        
                    } else if (!inputSpectralOrder->gethasSpectralEnergyDistribution()) {
                        if (args.verbose) cout << "operaMasterFluxCalibration: Skipping order number: "<< order << "for index=" << index << " no input SpectralEnergyDistribution." << endl;
                    }
                } // if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength() && spectralOrder->gethasSpectralEnergyDistribution()) {
            } //for (int order=minorder; order<=maxorder; order++) {
        } //for(unsigned index=0; index<inputFcalIndex; index++) {
        
        /*
         * and write output
         */
        spectralOrders.WriteSpectralOrders(outputfcal, Fcal);
        
        delete[] referenceWavelength;
        delete[] inputWavelength;
        delete[] inputFluxCal;
        delete[] inputThroughput;
        delete[] referenceInputFluxCal;
        delete[] referenceInputThroughput;
        
        for(unsigned beam = 0; beam < NumberofBeams; beam++) {
            delete[] inputBeamFluxCal[beam];
            delete[] inputBeamThroughput[beam];
            delete[] referenceBeamInputFluxCal[beam];
            delete[] referenceBeamInputThroughput[beam];
        }
        
        if (fspecdata.is_open() && foutdata.is_open()) {
            fspecdata.close();
            foutdata.close();
            
            if (!scriptfilename.empty()) {
                GenerateMasterFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(),outputDataFilename.c_str(), NumberofBeams, interactive);
            }
        }
	}
	catch (operaException e) {
		cerr << "operaMasterFluxCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaMasterFluxCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

void GenerateMasterFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string outputDataFilename, unsigned NumberofBeams, bool display)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str());  // delete any existing file with the same name
    ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "\nset xlabel \"wavelength (nm)\"" << endl;
    fgnu << "set ylabel \"flux calibration\"" << endl;
    
    fgnu << "set pointsize 0.5" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nplot \"" << spectrumDataFilename << "\" u 4:5 w d" <<
        ",\"" << outputDataFilename << "\" u 4:7 w l lw 3" << endl;
        
        if (display) {
            fgnu << "\nset terminal x11" << endl;
            fgnu << "set output" << endl;
            fgnu << "replot" << endl;
        } else {
            fgnu << "\n#set terminal x11" << endl;
            fgnu << "#set output" << endl;
            fgnu << "#replot" << endl;
        }
    } else {
        fgnu << "\nplot \"" << spectrumDataFilename << "\" u 4:5 w d" <<
        ",\"" << outputDataFilename << "\" u 4:7 w l lw 3" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if(!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

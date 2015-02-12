/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaCreateFluxCalibration
 Version: 1.0
 Description: Flux Calibration with Standard source 
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
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaLibCommon.h"
#include "libraries/operaFit.h"						// for operaLMFitPolynomial
#include "libraries/Polynomial.h"
#include "libraries/operaFFT.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaSpectralTools.h"			// void calculateUniformSample, getFluxAtWavelength

#include "core-espadons/operaCreateFluxCalibration.h"

#define NOTPROVIDED -999
#define MAXFLUXREFERENCELENGTH 20000
#define MAXNUMBEROFREFWLRANGES 1000

/*! \file operaCreateFluxCalibration.cpp */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*
 * the reference spectrum
 */
static unsigned nPointsInReferenceSpectrum = 0;
static double referenceWavelength[MAXFLUXREFERENCELENGTH];
static double referenceIntensity[MAXFLUXREFERENCELENGTH];
static double referenceNormIntensity[MAXFLUXREFERENCELENGTH];
static double referenceVariance[MAXFLUXREFERENCELENGTH];  

unsigned getReferenceSpectrumRange(double wl0, double wlf, double **wl, double **flux, double **normflux, double **fluxvar);
unsigned readReferenceSpectrum(string reference_spectrum, double *referenceWavelength, double *referenceIntensity, double *referenceVariance);

void normalizeIntensityByMaximum(unsigned np, double *intensity, double *variance);
void normalizeIntensityByReferenceWavelength(unsigned np, double *intensity, double *wavelength, double *outputNormIntensity, double refWavelength);
double getReferenceFlux(unsigned np, double *intensity, double *wavelength, double refWavelength);

double operaArrayMaxValue_d(unsigned np, const double *xarray, const double *yarray, double *maxx);
unsigned getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, float *refContinuumwl,float *refContinuumflux,float *refContinuumNormflux);
unsigned getReferenceSpectrumRange(unsigned nRefContinuum,double *refContinuumwl,double *refContinuumflux,double *refContinuumNormflux,double wl0,double wlf, double **wl, double **flux, double **normflux);


/*! 
 * operaCreateFluxCalibration
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
	
	string inputUncalibratedSpectrum;
	string inputCalibratedSpectrum;
    string inputFlatFluxCalibration;
    
	string outputFluxCalibrationFile;
    string inputWavelengthMaskForRefContinuum;
    string inputWavelengthMaskForUncalContinuum;
    
	string inputaper;
	string inputWaveFile;
	
    double wavelengthForNormalization = 548;
    
    unsigned numberOfPointsInUniformSample = 200;
    unsigned numberOfPointsInUniformRefSample = 70;
    
    double exposureTime = 1;   // in seconds

	struct pixelsize {
		unsigned x, y;
	} pixelsize = {1, 1};
    
	int ordernumber = NOTPROVIDED;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false;            
    
    unsigned binsize = 100;
    
    double delta_wl = 1.0; // wavelength (in nm) range for stiching non-overlapping orders
    
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
    string plotfilename;	
	string spectrumDataFilename;
	string continuumDataFilename;	
	string scriptfilename;	
	
	struct option longopts[] = {
		{"inputUncalibratedSpectrum",       1, NULL, 'i'},
		{"inputCalibratedSpectrum",         1, NULL, 'c'},
		{"inputFlatFluxCalibration",        1, NULL, 'f'},
        {"inputWavelengthMaskForRefContinuum",             1, NULL, 'm'},
        {"inputWavelengthMaskForUncalContinuum",           1, NULL, 'u'},
        {"inputWaveFile",                   1, NULL, 'w'},
		{"outputFluxCalibrationFile",       1, NULL, 'o'},
		{"inputApertureFile",               1, NULL, 'a'},   
		{"wavelengthForNormalization",		1, NULL, 'L'},
		{"exposureTime",                    1, NULL, 'E'},
		{"pixelsize",                       1, NULL, 'D'},
		{"ordernumber",                     1, NULL, 'O'},	
		{"minorder",                        1, NULL, 'M'},
		{"maxorder",                        1, NULL, 'X'},
		{"numberOfPointsInUniformSample",   1, NULL, 'l'},
		{"numberOfPointsInUniformRefSample",1, NULL, 'g'},
		{"binsize",                         1, NULL, 'b'},
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
	
	while((opt = getopt_long(argc, argv, "i:c:f:m:u:w:o:a:L:E:D:O:M:X:l:g:b:P:F:C:S:I:p::v::d::t::h",  longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputUncalibratedSpectrum = optarg;	
				break;    
			case 'c':
				inputCalibratedSpectrum = optarg;	
				break;
			case 'f':
				inputFlatFluxCalibration = optarg;
				break;
			case 'm':
				inputWavelengthMaskForRefContinuum = optarg;
				break;
			case 'u':
				inputWavelengthMaskForUncalContinuum = optarg;
				break;
			case 'w':
				inputWaveFile = optarg;
				break;
			case 'o':		// output
				outputFluxCalibrationFile = optarg;
				break;
			case 'a':
				inputaper = optarg;
				break;
			case 'L':
				wavelengthForNormalization = atof(optarg);
				break;
			case 'E':
				exposureTime = atof(optarg);
				break;
			case 'D':		// pixel size in microns
				if (strlen(optarg))
					sscanf(optarg, "%u %u", &pixelsize.x, &pixelsize.y);
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
			case 'l':
				numberOfPointsInUniformSample = atoi(optarg);
				break;
			case 'g':
				numberOfPointsInUniformRefSample = atoi(optarg);
				break;                
			case 'b':		// binsize
				binsize = atoi(optarg);
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
		if (inputUncalibratedSpectrum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input template spectrum...
		if (inputCalibratedSpectrum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need an input flat flux calibration spectrum...
		if (inputFlatFluxCalibration.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		// we need an output...
		if (outputFluxCalibrationFile.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a wavelength calibration file...
		if (inputWaveFile.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for ref...
		if (inputWavelengthMaskForRefContinuum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need a wavelength mask file for uncal...
		if (inputWavelengthMaskForUncalContinuum.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		// we need a aperture file...
		if (inputaper.empty()) {
			throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		if (verbose) {
			cout << "operaCreateFluxCalibration: input uncalibrated spectrum file = " << inputUncalibratedSpectrum << endl; 
			cout << "operaCreateFluxCalibration: input calibrated spectrum file = " << inputCalibratedSpectrum << endl;
			cout << "operaCreateFluxCalibration: inputFlatFluxCalibration = " << inputFlatFluxCalibration << endl;
			cout << "operaCreateFluxCalibration: inputWavelengthMaskForRefContinuum = " << inputWavelengthMaskForRefContinuum << endl;
			cout << "operaCreateFluxCalibration: inputWavelengthMaskForUncalContinuum = " << inputWavelengthMaskForUncalContinuum << endl;
			cout << "operaCreateFluxCalibration: inputWaveFile = " << inputWaveFile << endl;
            cout << "operaCreateFluxCalibration: output flux calibration file = " << outputFluxCalibrationFile << endl;
			cout << "operaCreateFluxCalibration: input aperture file = " << inputaper << endl;                        
            cout << "operaCreateFluxCalibration: wavelengthForNormalization= " << wavelengthForNormalization << " nm" << endl;
            cout << "operaCreateFluxCalibration: exposure time= " << exposureTime << " seconds" << endl;
			cout << "operaCreateFluxCalibration: numberOfPointsInUniformSample = " << numberOfPointsInUniformSample << endl;
			cout << "operaCreateFluxCalibration: numberOfPointsInUniformRefSample = " << numberOfPointsInUniformRefSample << endl;
			cout << "operaCreateFluxCalibration: pixelsize {x,y} = {" << pixelsize.x << "," << pixelsize.y << "}" << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaCreateFluxCalibration: ordernumber = " << ordernumber << endl;            
            }   
            cout << "operaCreateFluxCalibration: binsize = " << binsize << endl;  
            if(plot) {
                cout << "operaCreateFluxCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaCreateFluxCalibration: spectrumDataFilename = " << spectrumDataFilename << endl;
                cout << "operaCreateFluxCalibration: continuumDataFilename = " << continuumDataFilename << endl;
                cout << "operaCreateFluxCalibration: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cout << "operaCreateFluxCalibration: interactive = YES" << endl; 
                } else {
                    cout << "operaCreateFluxCalibration: interactive = NO" << endl; 
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
        
		operaSpectralOrderVector spectralOrders(inputUncalibratedSpectrum);
        spectralOrders.ReadSpectralOrders(inputaper);
        spectralOrders.ReadSpectralOrders(inputWaveFile);

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
			cout << "operaCreateFluxCalibration: minorder ="<< minorder << " maxorder=" << maxorder << endl;        

		/*
		 * Flux calibration reference file:
         *
		 * Read reference calibrated spectrum
		 *		lambda vs. intensity, intensityVariance (optional)
		 */
        nPointsInReferenceSpectrum = readReferenceSpectrum(inputCalibratedSpectrum, referenceWavelength, referenceIntensity, referenceVariance);        
        normalizeIntensityByReferenceWavelength(nPointsInReferenceSpectrum,referenceIntensity,referenceWavelength,referenceNormIntensity,wavelengthForNormalization);
        
        //---------------------------------
        // Loop over orders to set maximum number of elements, set wavelength and the number of beams
        // --> maxNElements & NumberofBeams
        unsigned NumberofBeams = spectralOrders.getNumberofBeams(minorder, maxorder);
        unsigned maxNElements = spectralOrders.getMaxNumberOfElementsInOrder(minorder, maxorder);
        
        for (int order=minorder; order<=maxorder; order++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            if (spectralOrder->gethasSpectralElements() && spectralOrder->gethasWavelength()) {
                operaWavelength *wavelength = spectralOrder->getWavelength();
                operaSpectralElements *SpectralElements = spectralOrder->getSpectralElements();                
                SpectralElements->setwavelengthsFromCalibration(wavelength);
            }
        }
        if (verbose) {
			cout << "operaCreateFluxCalibration: NumberofBeams = " << NumberofBeams << endl;
        }

        if(NumberofBeams == 0) {
            throw operaException("operaCreateFluxCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
        }

        //---------------------------------
        // below one can define constants for flux calibration
        // DT May 08 2014 -- not used -- double AreaOfOnePixel = (pixelsize.x * pixelsize.y); // for example, the area of a pixel could also be included
        double spectralBinConstant = exposureTime;
        double BeamSpectralBinConstant[MAXNUMBEROFBEAMS];
        for(unsigned beam=0; beam < NumberofBeams; beam++) {
            // uncalibratedContinuumBeamFluxForNormalization[beam] = spectralOrders.GetSpectralOrder(orderpicked)->getBeamSED(beam)->getUncalibratedFluxElements()->getFlux(elemIndexPicked);
            BeamSpectralBinConstant[beam] = exposureTime;
        }
        
        //---------------------------------
        // Correct flat-field
        if (!inputFlatFluxCalibration.empty()) {
            spectralOrders.correctFlatField(inputFlatFluxCalibration, minorder, maxorder, false);
        }

        //---------------------------------
        // Calculate a clean sample of the continuum from the ref spectrum
        float *refContinuumwl = new float[MAXNUMBEROFREFWLRANGES];
        float *refContinuumflux = new float[MAXNUMBEROFREFWLRANGES];
        float *refContinuumNormflux = new float[MAXNUMBEROFREFWLRANGES];
        unsigned nRefContinuum = getContinuumFromInputReferenceSpectrum(inputWavelengthMaskForRefContinuum,refContinuumwl,refContinuumflux,refContinuumNormflux);

        double refFluxForNormalization = (double)getFluxAtWavelength(nRefContinuum,refContinuumwl,refContinuumflux,wavelengthForNormalization);
        
        float *uniformRef_wl = new float[numberOfPointsInUniformRefSample];
        float *uniformRef_flux = new float[numberOfPointsInUniformRefSample];
        
        calculateUniformSample(nRefContinuum,refContinuumwl,refContinuumflux,numberOfPointsInUniformRefSample,uniformRef_wl,uniformRef_flux);

        if(debug) {
            // original sample
            for(unsigned i=0;i<nRefContinuum;i++) {
                cout << refContinuumwl[i] << " "
                << refContinuumflux[i] << " "
                << refContinuumNormflux[i] << endl;
            }
            // uniform sample
            for (unsigned i=0; i<numberOfPointsInUniformRefSample; i++) {
                cout << uniformRef_wl[i] << " " << uniformRef_flux[i] << endl;
            }
        }

        float *uniform_wl = new float[numberOfPointsInUniformSample];
        float *uniform_flux = new float[numberOfPointsInUniformSample];
        float *uniform_Beamflux[MAXNUMBEROFBEAMS];
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            uniform_Beamflux[beam] = new float[numberOfPointsInUniformSample];
        }
        
        spectralOrders.calculateCleanUniformSampleOfContinuum(minorder,maxorder,binsize,delta_wl,inputWavelengthMaskForUncalContinuum,numberOfPointsInUniformSample,uniform_wl,uniform_flux,uniform_Beamflux,TRUE);
        
        double uncalFluxForNormalization = (double)getFluxAtWavelength(numberOfPointsInUniformSample,uniform_wl,uniform_flux,wavelengthForNormalization);
        double uncalBeamFluxForNormalization[MAXNUMBEROFBEAMS];
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            uncalBeamFluxForNormalization[beam] = (double)getFluxAtWavelength(numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam],wavelengthForNormalization);
        }

        //---------------------------------
        // Use clean sample of the continuum from both the reference and uncalibrated spectrum
        // to calculate flux calibration and instrument throughput for all orders
        // Results are saved in the SED classes which will be written out to fcal
        
        float *CalibratedModelFlux = new float[maxNElements];
        float *UncalibratedModelFlux = new float[maxNElements];
        float *UncalibratedModelBeamFlux[MAXNUMBEROFBEAMS];
        for(unsigned beam=0;beam<NumberofBeams;beam++) {
            UncalibratedModelBeamFlux[beam] = new float[maxNElements];
        }
        float *elemWavelength = new float[maxNElements];
        
        operaSpectralEnergyDistribution *BeamSED[MAXNUMBEROFBEAMS];
        operaSpectralElements *BeamFluxCalibrationElements[MAXNUMBEROFBEAMS];
        operaSpectralElements *BeamThroughputElements[MAXNUMBEROFBEAMS];

        for (int order=minorder; order<=maxorder; order++) {
            
			operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            if (spectralOrder->gethasWavelength() &&
                spectralOrder->gethasSpectralElements() &&
                spectralOrder->gethasSpectralEnergyDistribution()) {
				            
				operaSpectralElements *spectralElements = spectralOrder->getSpectralElements();
                operaSpectralEnergyDistribution *spectralEnergyDistribution = spectralOrder->getSpectralEnergyDistribution();
                spectralEnergyDistribution->setwavelengthForNormalization(wavelengthForNormalization);
                
                unsigned nElements = spectralElements->getnSpectralElements();

                spectralEnergyDistribution->setUncalibratedFluxElements(spectralElements);
                spectralEnergyDistribution->setCalibratedFluxElements(spectralElements);
                
                operaSpectralElements *calibratedFluxElements = spectralEnergyDistribution->getCalibratedFluxElements();
                operaSpectralElements *uncalibratedFluxElements = spectralEnergyDistribution->getUncalibratedFluxElements();
                
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    BeamSED[beam] = spectralOrder->getBeamSED(beam);
                    BeamSED[beam]->setUncalibratedFluxElements(spectralOrder->getBeamElements(beam));
                    BeamSED[beam]->setCalibratedFluxElements(spectralOrder->getBeamElements(beam));
                }
                
                for(unsigned i=0;i<nElements;i++) {
                    elemWavelength[i] =  (float)spectralElements->getwavelength(i);
                }
                if(!continuumDataFilename.empty()) {
                    for(unsigned i=0;i<numberOfPointsInUniformSample;i++) {
                       *fcontinuumdata << uniform_wl[i] << ' ' << uniform_flux[i] << endl;
                    }
                }
                operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_flux,nElements,elemWavelength,UncalibratedModelFlux);
                operaFitSpline(numberOfPointsInUniformRefSample,uniformRef_wl,uniformRef_flux,nElements,elemWavelength,CalibratedModelFlux);
                
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    operaFitSpline(numberOfPointsInUniformSample,uniform_wl,uniform_Beamflux[beam],nElements,elemWavelength,UncalibratedModelBeamFlux[beam]);
                }
                
                operaSpectralElements *FluxCalibrationElements = spectralEnergyDistribution->getFluxCalibrationElements();
                operaSpectralElements *ThroughputElements = spectralEnergyDistribution->getThroughputElements();
                for(unsigned beam=0;beam<NumberofBeams;beam++) {
                    BeamFluxCalibrationElements[beam] = BeamSED[beam]->getFluxCalibrationElements();
                    BeamThroughputElements[beam] = BeamSED[beam]->getThroughputElements();
                }
                
                for(unsigned i=0;i<nElements;i++) {
                    calibratedFluxElements->setFlux((double)CalibratedModelFlux[i]/refFluxForNormalization,i);
                    uncalibratedFluxElements->setFlux((double)UncalibratedModelFlux[i]/uncalFluxForNormalization,i);
                    
                    FluxCalibrationElements->setFlux((double)(CalibratedModelFlux[i]/(UncalibratedModelFlux[i]/spectralBinConstant)),i);
                    ThroughputElements->setFlux((double)(CalibratedModelFlux[i]/UncalibratedModelFlux[i])*(uncalFluxForNormalization/refFluxForNormalization),i);
                    
                    for(unsigned beam=0;beam<NumberofBeams;beam++) {
                        BeamFluxCalibrationElements[beam]->setFlux((double)(CalibratedModelFlux[i]/(UncalibratedModelBeamFlux[beam][i]/BeamSpectralBinConstant[beam])),i);
                        BeamThroughputElements[beam]->setFlux((double)(CalibratedModelFlux[i]/UncalibratedModelBeamFlux[beam][i])*(uncalBeamFluxForNormalization[beam]/refFluxForNormalization),i);
                    }
                    
                    //cout << elemWavelength[i] << " " << CalibratedModelFlux[i] << endl;
                }
                spectralEnergyDistribution->setHasCalibratedFlux(true);
                spectralEnergyDistribution->setHasUncalibratedFlux(true);
                
                if(!spectrumDataFilename.empty()) {
                    for(unsigned i=0; i<spectralElements->getnSpectralElements(); i++) {
                        *fspecdata << order << " " << i << " "
                        << spectralElements->getwavelength(i) << " "
                        << spectralElements->getFlux(i) << " "
                        << UncalibratedModelFlux[i] << " "
                        << CalibratedModelFlux[i] << " "
                        << spectralEnergyDistribution->getUncalibratedFluxElements()->getFlux(i) << " "
                        << spectralEnergyDistribution->getCalibratedFluxElements()->getFlux(i) << " "
                        << spectralEnergyDistribution->getFluxCalibrationElements()->getFlux(i) << " "
                        << spectralEnergyDistribution->getThroughputElements()->getFlux(i) << " ";
                        for(unsigned beam=0;beam<NumberofBeams;beam++) {
                            *fspecdata << BeamSED[beam]->getUncalibratedFluxElements()->getFlux(i) << " "
                                << BeamSED[beam]->getFluxCalibrationElements()->getFlux(i) << " "
                                << BeamSED[beam]->getThroughputElements()->getFlux(i) << " ";
                        }
                        *fspecdata << endl;
                    }
                }
            } else {
                spectralOrder->sethasSpectralEnergyDistribution(FALSE);
            }
        }
        
		/*
		 * and write out fcal
		 */
		spectralOrders.WriteSpectralOrders(outputFluxCalibrationFile, Fcal);
		
        if (fspecdata != NULL && fcontinuumdata != NULL) {
			fspecdata->close();
            fcontinuumdata->close();
            
            if (!scriptfilename.empty()) {
                GenerateCreateFluxCalibrationPlot(scriptfilename.c_str(),plotfilename.c_str(),spectrumDataFilename.c_str(),continuumDataFilename.c_str(), NumberofBeams, interactive);
            }
        }
        
	}
	catch (operaException e) {
		cerr << "operaCreateFluxCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaCreateFluxCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
} 

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --inputUncalibratedSpectrum=<SPEC_FILE>"
	" --inputCalibratedSpectrum=<TEMPLATE_FILE>"    
	" --inputWaveFile=<WAVE_FILE>"
	" --outputFluxCalibrationFile=<SPEC_FILE>"
	" --inputApertureFile=<APER_FILE>"
	" --spectrumtype=<UNS_VALUE>"
	" --wavelengthForNormalization=<DOUBLE_VALUE>"
	" --exposureTime=<FLOAT_VALUE>"
	" --pixelsize=<FLOAT_VALUE>"
	" --ordernumber=<UNS_VALUE>"
	" --minorder=<UNS_VALUE>"
	" --maxorder=<UNS_VALUE>"
	" --numberOfPointsInUniformSample=<UNS_VALUE>"
	" --numberOfPointsInUniformRefSample=<UNS_VALUE>"
	" --binsize=<UNS_VALUE>"
	" --usePolynomial=<BOOL>"
	" --orderOfPolynomial=<UNS_VALUE>"
	" --generate3DPlot=<BOOL>"
	" --generateBeamPlot=<BOOL>"
	" --plotContinuum=<BOOL>"
	" --plotfilename=<EPS_FILE>"
	" --spectrumDataFilename=<DATA_FILE>"
	" --continuumDataFilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputCalibratedSpectrum=HR1544_operaFluxCal.dat --inputUncalibratedSpectrum=1515004.e.gz --inputWaveFile=/Users/edermartioli/opera/calibrations/GalileanMoons/OLAPAa_pol_Normal.wcar.gz --outputFluxCalibrationFile=1515004.fcal.gz --inputApertureFile=/Users/edermartioli/opera/calibrations/GalileanMoons/OLAPAa_pol_Normal.aper.gz --normalizeCalibratedSpectrum=1 --binsize=210 --exposureTime=1 --spectrumDataFilename=1515004.spec --continuumDataFilename=1515004.cont --scriptfilename=1515004fcal.gnu -v \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -i, --inputUncalibratedSpectrum=<SPEC_FILE>,  Spectrophotometric standard extracted uncalibrated spectrum \n"
	"  -c, --inputCalibratedSpectrum=<TEMPLATE_FILE>,  Spectrophotometric standard template calibrated spectrum \n"
	"  -w, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file\n"
	"  -o, --outputFluxCalibrationFile=<SPEC_FILE>,  Output flux calibration conversion file \n"
	"  -a, --inputApertureFile=<APER_FILE>, Input aperture calibration file\n"
	"  -L, --wavelengthForNormalization=<DOUBLE_VALUE>, Wavelength (nm) for normalization of reference spectrum\n"
	"  -E, --exposureTime=<WAVE_FILE>, Exposure time of input uncalibrated spectrum\n"
	"  -D, --pixelsize=<WAVE_FILE>, Detector pixel size\n"
	"  -O, --ordernumber=<UNS_VALUE>, Absolute order number to extract (default=all)\n"
	"  -M, --minorder=<UNS_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<UNS_VALUE>, Define maximum order number\n"
	"  -l, --numberOfPointsInUniformSample=<UNS_VALUE>, Define lowest order to consider in the fit across orders\n"
	"  -g, --numberOfPointsInUniformRefSample=<UNS_VALUE>, Define highest order to consider in the fit across orders\n"
	"  -b, --binsize=<UNS_VALUE>, Number of points to bin for continuum estimate \n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -E, --generate3DPlot=<BOOL>, Switch to generate 3D or 2D plot spectra\n"
	"  -B, --generateBeamPlot=<BOOL>, Switch to generate plot of beams or full slit spectra\n"
	"  -c, --plotContinuum=<BOOL>, Switch to generate plot of continuum or normalized line spectra\n"
	"  -F, --spectrumDataFilename=<DATA_FILE>\n"
	"  -C, --continuumDataFilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateCreateFluxCalibrationPlot(string gnuScriptFileName, string outputPlotEPSFileName, string spectrumDataFilename, string continuumDataFilename, unsigned NumberofBeams, bool display)
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

/*
 * Read the the full reference spectrum
 */
unsigned readReferenceSpectrum(string reference_spectrum, double *referenceWavelength, double *referenceIntensity, double *referenceVariance) {
	ifstream astream;
	string dataline;
    
	double tmpwl = -1.0; 
	double tmpi = -1.0; 
	unsigned np = 0;
	
	astream.open(reference_spectrum.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf", &tmpwl, &tmpi);
                    
                    referenceWavelength[np] = tmpwl;
                    referenceIntensity[np] = tmpi;
                    referenceVariance[np] = tmpi;
                    np++;  
                }	// skip comments
            }
		} // while (astream.good())
        
		if (np > 0) {
			if (verbose) {
				printf("          [Reference] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, referenceWavelength[0], referenceWavelength[np/2], referenceWavelength[np-1]);
			}
		} else {
			printf("          [Reference] no points found in flux reference file.\n");
		}
		astream.close();
	}	// if (astream.open())
	return np;
}

/*
 * get a subset of the reference spectrum for this order only, between wl0 and wlf
 */
unsigned getReferenceSpectrumRange(double wl0, double wlf, double **wl, double **flux, double **normflux, double **fluxvar) {
	unsigned firstline = 0;
	unsigned np = 0;
    
	for (np=0; np<nPointsInReferenceSpectrum; np++) {
		if (referenceWavelength[np] >= wl0) {
			if (firstline == 0 && np) {
                    *flux = &referenceIntensity[np-1];
                    *normflux = &referenceNormIntensity[np-1];
                    *wl = &referenceWavelength[np-1];
                    *fluxvar = &referenceVariance[np-1];
                    firstline = np-1;
			}
			if (referenceWavelength[np] > wlf) {
                np++;
				break;
            }
		}
	}
	if (np == nPointsInReferenceSpectrum)
		np--;
	if (np > firstline) {
		return (np-firstline);
	} else {
		return 0;
	}
}

void normalizeIntensityByReferenceWavelength(unsigned np, double *intensity, double *wavelength, double *outputNormIntensity, double refWavelength) {
    double referenceFlux = getReferenceFlux(np,intensity,wavelength,refWavelength);
    
	for(unsigned i=0;i<np;i++) {
        outputNormIntensity[i] = intensity[i]/referenceFlux;
	}
}

double getReferenceFlux(unsigned np, double *intensity, double *wavelength, double refWavelength) {
    
    float *wavelengthData_f = new float[np];
    float *fluxData_f = new float[np];
    
    for(unsigned i=0;i<np;i++) {
        wavelengthData_f[i] = (float)wavelength[i];
        fluxData_f[i] =  (float)intensity[i];
    }
    
    unsigned nElements = 1;
    
    float *referenceFlux = new float[nElements];
    float *referencewl = new float[nElements];
    
    for(unsigned i=0;i<nElements;i++) {
        referencewl[i] = (float)refWavelength;
    }
    operaFitSpline(np,wavelengthData_f,fluxData_f,nElements,referencewl,referenceFlux);
    
    double outputflux = (double)(referenceFlux[0]);
    
    delete[] wavelengthData_f;
    delete[] fluxData_f;
    delete[] referenceFlux;
    delete[] referencewl;
    
	return outputflux;
}

void normalizeIntensityByMaximum(unsigned np, double *intensity, double *variance) {
	double maxIntensity = -3.4e+38;
	
	for(unsigned i=0;i<np;i++) {
		if(intensity[i] > maxIntensity)
			maxIntensity = intensity[i];
	}
	
	for(unsigned i=0;i<np;i++) {
        intensity[i] /= maxIntensity;
        variance[i] /= maxIntensity;
	}
}


double operaArrayMaxValue_d(unsigned np, const double *xarray, const double *yarray, double *maxx) {
	double ymax = -3.4e+38;
	double xmax = 0;
    
	while (np--) {
		if(*yarray > ymax) {
			ymax = *yarray;
            xmax = *xarray;
        }
		yarray++;
        xarray++;
	}
    *maxx = xmax;
	return ymax;
}


unsigned getContinuumFromInputReferenceSpectrum(string inputWavelengthMaskForRefContinuum, float *refContinuumwl,float *refContinuumflux,float *refContinuumNormflux) {
    
    double *wl0_vector = new double[MAXNUMBEROFREFWLRANGES];
    double *wlf_vector = new double[MAXNUMBEROFREFWLRANGES];
    
    unsigned nRangesInWLMask = readContinuumWavelengthMask(inputWavelengthMaskForRefContinuum,wl0_vector,wlf_vector);
    
    unsigned nTotalPoints = 0;
    for(unsigned k=0;k<nRangesInWLMask; k++){
        
        double *refwl = NULL, *refflux = NULL, *refnormflux = NULL, *reffluxvar = NULL;
        unsigned nPointsInReference = getReferenceSpectrumRange(wl0_vector[k],wlf_vector[k],&refwl,&refflux,&refnormflux,&reffluxvar);
        
        double ref_wl=0;
        double ref_maxFlux = 0;
        double ref_maxNormFlux = 0;
        
        ref_maxFlux = operaArrayMaxValue_d(nPointsInReference,refwl,refflux,&ref_wl);
        ref_maxNormFlux = operaArrayMaxValue_d(nPointsInReference,refwl,refnormflux, &ref_wl);
        
        if(debug) {
            cout << k << " "
            << nPointsInReference << " "
            << wl0_vector[k] << " "
            << wlf_vector[k] << " "
            << ref_wl << " "
            << ref_maxFlux << " "
            << ref_maxNormFlux << endl;
        }
        
        if(ref_wl && ref_maxFlux && ref_maxNormFlux) {
            refContinuumwl[nTotalPoints] = (float)ref_wl;
            refContinuumflux[nTotalPoints] = (float)ref_maxFlux;
            refContinuumNormflux[nTotalPoints] = (float)ref_maxNormFlux;
            nTotalPoints++;
        }
    }
    
    delete[] wl0_vector;
    delete[] wlf_vector;
    
    return nTotalPoints;
}


/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaInstrumentProfileCalibration
 Version: 1.0
 Description: Create the Instrument Profile 
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
#include <pthread.h>

#include <fstream>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaInstrumentProfileCalibration.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaEspadonsImage.h"			// for imtype_t
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaInstrumentProfile.h"		// for operaInstrumentProfile
#include "libraries/operaExtractionAperture.h"
#include "libraries/operaSpectralLines.h"

#include "libraries/operaMath.h"                    // for LengthofPolynomial
#include "libraries/operaLibCommon.h"
#include "libraries/operaLib.h"						// systemf
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaFit.h"	
#include "libraries/operaFFT.h"	
#include "libraries/operaParameterAccess.h"
#include "libraries/operaConfigurationAccess.h"

#define NOTPROVIDED -999

/*! \file operaInstrumentProfileCalibration.cpp */

using namespace std;

/*!
 * operaInstrumentProfileCalibration
 * \author Eder Martioli
 * \brief Create the Instrument Profile.
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
 * These variable have been made global for thread support.
 * Teeple Dec 20 2012
 */

int debug=0, verbose=0, trace=0, plot=0;

string geometryfilename;
string outputprof;

string masterbias;
string masterflat;
string mastercomparison;
string badpixelmask;
string masterfabperot;
string gainfilename;

string plotfilename;
string datafilename;
string scriptfilename;

float spectralElementHeight = 1.0;
float referenceLineWidth = 2.5;

float gain = 1.12;
float noise = 3.5;

float LocalMaxFilterWidth = 2.5*referenceLineWidth;
float DetectionThreshold = 0.2;
float MinPeakDepth = 1.5*noise;

unsigned minorder = 22;
bool minorderprovided = false;
unsigned maxorder = 62;
bool maxorderprovided = false;

int ordernumber = NOTPROVIDED;

/* "method" defines how all individual IP measurements from spectral lines should be combined.
 *  It supports the following methods:
 *      1. weighted mean
 *      2. median combine
 *      3. polynomial fit
 *      4. polynomial fit with median binning (use binsize provided)
 */
int method = 1;

unsigned binsize = 80;
float tilt = -3.0;

struct ipDimensions {
    unsigned xsize, ysize;
    unsigned xsampling, ysampling;
} ipDimensions = {0,0,0,0};  // {Xsize, Xsampling, Ysize, Ysampling}

operaFITSImage *fabperot = NULL;
operaFITSImage *badpix = NULL;
operaFITSImage *bias = NULL;
operaFITSImage *flat = NULL;
operaFITSImage *comp = NULL;

unsigned IPxsize = 0;
unsigned IPxsampling = 0;
unsigned IPysize = 0;
unsigned IPysampling = 0;

operaSpectralOrderVector spectralOrders;

unsigned maxthreads = 1;

bool interactive = false;

/*
 * Thread Support to process all orders in parallel
 */

typedef struct thread_args {
	int order;
} thread_args_t;

pthread_t *threads = NULL;
thread_args_t *thread_args = NULL;

void *processOrder(void *argument) {
	thread_args_t *thread_args_s = (thread_args_t *)argument;
	int order = thread_args_s->order;
    
    if (verbose) {
        cout << "operaInstrumentProfileCalibration: Processing order = " << order << endl;
    }
    // create pointer to current spectral order:
    operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
    
    // create a set of spectral elements based on given element height:
    spectralOrder->setSpectralElementsByHeight(spectralElementHeight);
    
    // create vector to store information of instrument profile for current spectral order:
    unsigned MaxNDataPoints = spectralOrder->getSpectralElements()->getnSpectralElements();
    
    unsigned sampleElementForPlot = spectralOrder->getSpectralElements()->getnSpectralElements()/2;
    
    spectralOrder->setInstrumentProfileVector(IPxsize,IPxsampling,IPysize,IPysampling,MaxNDataPoints);
    
    spectralOrder->measureInstrumentProfileAlongRowsInto2DWithGaussian(*flat,*badpix,binsize,referenceLineWidth,tilt,false,sampleElementForPlot,NULL);
    
    if(spectralOrder->gethasInstrumentProfile()) {
        if (fabperot) {
            double MaxContamination = 1.0;
            double amplitudeCutOff = 3*noise;
            unsigned nSigCut = 3;
            
            operaSpectralLines *spectralLines = NULL;
            try {
                spectralOrder->calculateXCorrBetweenIPandImage(*fabperot, *badpix, NULL);
                spectralOrder->setSpectralLines(*fabperot, *badpix, *bias, noise, gain, referenceLineWidth, DetectionThreshold, LocalMaxFilterWidth, MinPeakDepth);
                spectralOrder->sethasSpectralLines(true);
                spectralLines = spectralOrder->getSpectralLines();
                if (verbose)
                    cout << "operaInstrumentProfileCalibration: " << spectralLines->getnLines() << " lines found in order " << order << " of Fabry-Perot." << endl;
            } catch (operaException e) {
                if (verbose)
                    cout << "operaInstrumentProfileCalibration: No lines found in order " << order << " of Fabry-Perot." <<endl;
                spectralOrder->sethasInstrumentProfile(false);
            }
            if (spectralOrder->gethasSpectralLines() && spectralLines && spectralLines->getnLines() > 0) {
                if(method == 1) {
                    spectralOrder->measureInstrumentProfileUsingWeightedMean(*fabperot, *badpix, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot,NULL);
                } else if (method == 2) {
                    spectralOrder->measureInstrumentProfileUsingMedian(*fabperot, *badpix, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot, NULL);
                } else if (method == 3) {
                    spectralOrder->measureInstrumentProfile(*fabperot, *badpix, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot, NULL);
                } else if (method == 4) {
                    spectralOrder->measureInstrumentProfileWithBinning(*fabperot, *badpix, binsize, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot, NULL);
                }
                spectralOrder->recenterOrderPosition();
            }
            
        } else {
            double MaxContamination = 0.05;      // accept up to 1% flux contamination from neighbor line
            double amplitudeCutOff = 2*noise;  // limit lines with amplitude greater than 10 x CCD noise
            unsigned nSigCut = 2;               // limit lines with width up to twice larger than the median deviation
            
            operaSpectralLines *spectralLines = NULL;
            try {
                spectralOrder->calculateXCorrBetweenIPandImage(*comp, *badpix, NULL);
                spectralOrder->setSpectralLines(*comp, *badpix, *bias, noise, gain, referenceLineWidth, DetectionThreshold, LocalMaxFilterWidth, MinPeakDepth);
                spectralOrder->sethasSpectralLines(true);
                spectralLines = spectralOrder->getSpectralLines();
                if (verbose)
                    cout << "operaInstrumentProfileCalibration: " << spectralLines->getnLines() << " lines found in order " << order << " of Comparison." << endl;
            } catch (operaException e) {
                if (verbose)
                    cout << "operaInstrumentProfileCalibration: No lines found in order " << order << " of Comparison." << endl;
                spectralOrder->sethasInstrumentProfile(false);
            }
            if(spectralOrder->gethasSpectralLines() && spectralLines != NULL && spectralLines->getnLines() > 0) {
                if(method == 1) {
                    spectralOrder->measureInstrumentProfileUsingWeightedMean(*comp, *badpix, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot, NULL);
                } else if (method == 2) {
                    spectralOrder->measureInstrumentProfileUsingMedian(*comp, *badpix, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot, NULL);
                } else if (method == 3) {
                    spectralOrder->measureInstrumentProfile(*comp, *badpix, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot, NULL);
                } else if (method == 4) {
                    spectralOrder->measureInstrumentProfileWithBinning(*comp, *badpix, binsize, MaxContamination, amplitudeCutOff, nSigCut,sampleElementForPlot, NULL);
                }
                
                spectralOrder->recenterOrderPosition();
            }
        }
    }
    // clean up so we don't run out of memory
    spectralOrder->getInstrumentProfile()->deleteDataCubes();
    if (spectralOrder->gethasSpectralLines()) {
        delete spectralOrder->getSpectralLines();
        spectralOrder->sethasSpectralLines(false);
    }
    if (spectralOrder->gethasSpectralElements()) {
        delete spectralOrder->getSpectralElements();
        spectralOrder->sethasSpectralElements(false);
    }
	return NULL;
}

static void processSingleOrder(int order) {
    thread_args[0].order = order;
    processOrder((void *) &thread_args[0]);
}

static bool spawnthreads(int order, int maxorder, int count) {
	int j = 0;
    for (int i=order; i<=maxorder; i++) {
		thread_args[i].order = i;
		if (pthread_create(&threads[i], NULL, processOrder, (void *) &thread_args[i]) != 0)
			return false;
        if (++j >= count)
            break;
	}
    return true;
}

static bool waitthreads(int order, int maxorder, int count) {
	int j = 0;
	for (int i=order; i<=maxorder; i++) {
		if (pthread_join(threads[i], NULL) != 0)
			return false;
        if (++j >= count)
            break;
	}
    return true;
}

static bool processOrders(int minorder, int maxorder) {
	for (int order=minorder; order<=maxorder; order+=maxthreads) {
        spawnthreads(order, maxorder, maxthreads);
        waitthreads(order, maxorder, maxthreads);
	}
	return true;
}

int main(int argc, char *argv[])
{
	int opt;
	
	
	struct option longopts[] = {
		{"outputProf",1, NULL, 'o'},         
		{"gainfilename",1, NULL, 'G'},        
		{"geometryfilename",1, NULL, 'g'},
		{"masterbias",1, NULL, 'b'},
		{"masterflat",1, NULL, 'f'},
		{"mastercomparison",1, NULL, 'c'},
		{"badpixelmask",1, NULL, 'm'},        
		{"masterfabperot",1, NULL, 'a'},
		{"method",1, NULL, 'M'},	
		{"spectralElementHeight",1, NULL, 'H'},
		{"referenceLineWidth",1, NULL, 'W'},
		{"LocalMaxFilterWidth",1, NULL, 'L'},	// DT May 19 2014 movded to parameters
		{"DetectionThreshold",1, NULL, 'R'},	// DT May 19 2014 movded to parameters
		{"MinPeakDepth",1, NULL, 'E'},			// DT May 19 2014 movded to parameters
        {"ordernumber",1, NULL, 'O'},
		{"minorder",1, NULL, 'N'},
		{"maxorder",1, NULL, 'X'},         
        {"binsize",1, NULL, 'B'},
        {"tilt",1, NULL, 'T'},
		{"ipDimensions",1, NULL, 'D'},
		{"maxthreads",1, NULL, 'x'},
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},
		{"interactive",0, NULL, 'I'},
		
		{"plot",		optional_argument, NULL, 'p'},       
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "o:G:g:b:f:c:m:a:M:A:H:W:O:N:X:B:L:R:T:E:D:P:x:F:S:I:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'o':
				outputprof = optarg;
				break;                  
			case 'G':		// gain
				gainfilename = optarg;
				break;                
			case 'g':
				geometryfilename = optarg;
				break;     
			case 'b':
				masterbias = optarg;
				break;
			case 'f':		// masterflat
				masterflat = optarg;
				break;
			case 'c':
				mastercomparison = optarg;
				break;				
			case 'm':		// badpixelmask
				badpixelmask = optarg;
				break;
			case 'a':		// masterfabperot
				masterfabperot = optarg;
				break;
			case 'M':		// masterfabperot
				method = atoi(optarg);
				break;                  
			case 'H':		// element height in pixels
				spectralElementHeight = atof(optarg);
				break; 	
			case 'W':		// line width for reference (in pixel units)
				referenceLineWidth = atof(optarg);
				break;                 
			case 'L':
				LocalMaxFilterWidth = atof(optarg);
				break;                 
			case 'R':
				DetectionThreshold = atof(optarg);
				break;                 
			case 'E':
				MinPeakDepth = atof(optarg) * noise;
				break;                 
			case 'T':		// tilt in degree
				tilt = atof(optarg);
				break;                 
			case 'O':
				ordernumber = atoi(optarg);
				break;	
			case 'N':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break;                
			case 'B':
				binsize = atoi(optarg);
				break;                
			case 'D':		// ipDimensions
				if (strlen(optarg))
					sscanf(optarg, "%u %u %u %u", &ipDimensions.xsize, &ipDimensions.xsampling, &ipDimensions.ysize, &ipDimensions.ysampling);
				break;                  
			case 'x':
                maxthreads = atoi(optarg);
                break;
			case 'P':
				plotfilename = optarg;
				plot = 1;
				break; 		                
			case 'F':
				datafilename = optarg;
				break; 	
			case 'S':
				scriptfilename = optarg;
				break;
			case 'I':		// for interactive plots
				interactive = true;
				break;

			case 'p':
				plot = 1;
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
        if(method!=1 && method !=2 && method != 3 && method != 4) {
            method = 1;
            cout <<  "operaInstrumentProfileCalibration: Warning, undefined method. Using method = 1 (weighted mean)" << endl;	
        }
		// we need a geometry file...
		if (geometryfilename.empty()) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a masterbias...
		if (masterbias.empty()) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a masterflat...
		if (masterflat.empty()) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a comparison...
		if (mastercomparison.empty()) {
			throw operaException("operaInstrumentProfileCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}								
		
		if (verbose) {
			cout << "operaInstrumentProfileCalibration: geometryfilename = " << geometryfilename << endl; 
			cout << "operaInstrumentProfileCalibration: outputProf = " << outputprof << endl;
			cout << "operaInstrumentProfileCalibration: masterbias = " << masterbias << endl;            
			cout << "operaInstrumentProfileCalibration: masterflat = " << masterflat << endl; 	
			cout << "operaInstrumentProfileCalibration: mastercomparison = " << mastercomparison << endl; 	
			cout << "operaInstrumentProfileCalibration: badpixelmask = " << badpixelmask << endl;     
			cout << "operaInstrumentProfileCalibration: method = " << method << endl;                 
			cout << "operaInstrumentProfileCalibration: masterfabperot = " << masterfabperot << endl; 			
			cout << "operaInstrumentProfileCalibration: spectralElementHeight = " << spectralElementHeight << endl;            
			cout << "operaInstrumentProfileCalibration: referenceLineWidth = " << referenceLineWidth << endl;
			cout << "operaInstrumentProfileCalibration: maxthreads = " << maxthreads << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaInstrumentProfileCalibration: ordernumber = " << ordernumber << endl;            
            }
            if(plot) {
                cout << "operaInstrumentProfileCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaInstrumentProfileCalibration: datafilename = " << datafilename << endl;
                cout << "operaInstrumentProfileCalibration: scriptfilename = " << scriptfilename << endl;                
            }            
		}
        
        ofstream *fdata = NULL;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }
        
        bias = new operaFITSImage(masterbias, tfloat, READONLY);
		flat = new operaFITSImage(masterflat, tfloat, READONLY);
		comp = new operaFITSImage(mastercomparison, tfloat, READONLY);
		
//		flat -= bias;			// remove bias from masterflat
		
        IPxsize = ipDimensions.xsize;
        IPxsampling = ipDimensions.xsampling;
        IPysize = ipDimensions.ysize;
        IPysampling = ipDimensions.ysampling;
        
		if (!masterfabperot.empty()){
			fabperot = new operaFITSImage(masterfabperot, tfloat, READONLY);
		}
		
		if (!badpixelmask.empty()){              
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat->getnaxis1(),flat->getnaxis2(),tfloat);
            *badpix = 1.0;
        }
        
		spectralOrders.ReadSpectralOrders(geometryfilename);
		spectralOrders.readGainNoise(gainfilename);
		
        unsigned amp = 0;
		gain = spectralOrders.getGainBiasNoise()->getGain(amp);
		noise = spectralOrders.getGainBiasNoise()->getNoise(amp);
		//double biaslevel = spectralOrders.getGainBiasNoise()->getBias(amp);
        
        if(!IPxsize)
            IPxsize = 30;
        if(!IPxsampling)
            IPxsampling = 5;
        if(!IPysize)
            IPysize = 6;
        if(!IPysampling)
            IPysampling = 5;
#if 0
        //
		// -- DT May 19 2014
		// ... moved to parameters
        // -- E. Martioli 14 May 2014
        // Changed values below to increase sensitivity on fainter lines
        // in GRACES comparison spectra. Previous values were the following:
		//
        // LocalMaxFilterWidth = 2.5*referenceLineWidth;
        // DetectionThreshold = 0.2;
        // MinPeakDepth = 1.5*noise;
        // Note: the problem here could be when using FP exposures, where it
        //       could detect ghosts as if they were specral lines.
        
        if(!LocalMaxFilterWidth)
            LocalMaxFilterWidth = 2.5*referenceLineWidth;
        if(!DetectionThreshold)
            DetectionThreshold = 0.1;
        if(!MinPeakDepth)
            MinPeakDepth = 1.25*noise;
#endif        
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
			cout << "operaInstrumentProfileCalibration: binsize = " << binsize << endl;            
			cout << "operaInstrumentProfileCalibration: ipDimensions = [" << IPxsize << " pxl (" << IPxsampling << " ppp)," << IPysize <<" pxl ("<< IPysampling << " ppp)]" << endl; 
			cout << "operaInstrumentProfileCalibration: minorder = " << minorder << " maxorder = " << maxorder << endl;            
		}

        unsigned long nthreads = maxorder+1;
        threads = (pthread_t *)calloc(nthreads, sizeof(pthread_t*));
        thread_args = (thread_args_t *)calloc(nthreads, sizeof(thread_args_t));

        if (maxthreads > 1) {
            processOrders(minorder, maxorder);
        } else {
            for (unsigned order=minorder; order<=maxorder; order++) {
                processSingleOrder(order);
            }
        }
        
        if (fdata != NULL) {
            int minorderWithIP = NOTPROVIDED;
            int maxorderWithIP = NOTPROVIDED;
            
            for (unsigned order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                if(spectralOrder->gethasGeometry() && spectralOrder->gethasInstrumentProfile()) {
                    operaGeometry *Geometry = spectralOrder->getGeometry();
                    operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
                    float distd = (float)Geometry->CalculateDistance(Geometry->getYmin(), (Geometry->getYmax() - Geometry->getYmin())/2);
                    instrumentProfile->printModel(distd,order,fdata);
                    if(minorderWithIP == NOTPROVIDED){
                        minorderWithIP = (int)order;
                    }
                    maxorderWithIP = (int)order;
                }
            }
            
            fdata->close();
            if (!scriptfilename.empty()) {
                GenerateInstrumentProfile3DPlot(scriptfilename,plotfilename,datafilename,minorderWithIP,maxorderWithIP, IPxsize, IPysize, interactive);
            }
        }
		
		//        spectralOrders.WriteSpectralOrders(geometryfilename, Geom);	// THIS IS A SOMEWHAT DANGEROUS BYPRODUCT
		spectralOrders.WriteSpectralOrders(outputprof, Prof);

		bias->operaFITSImageClose();
		flat->operaFITSImageClose();
		comp->operaFITSImageClose();
        
        if(badpix) {
            delete badpix;
        }
		if (fabperot){
            delete fabperot;
		}
		if (bias){
            delete bias;
		}
		if (flat){
            delete flat;
		}
		if (comp){
            delete comp;
		}
	}
	catch (operaException e) {
		cerr << "operaInstrumentProfileCalibration: "  << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaInstrumentProfileCalibration: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaInstrumentProfileCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	" --outputProf=<PROF_FILE>"    
    " --gainfilename=<GAIN_FILE>"    
	" --geometryfilename=<GEOM_FILE>"
	" --masterbias=<FITS_IMAGE>"	
	" --masterflat=<FITS_IMAGE>"			
	" --mastercomparison=<FITS_IMAGE>"	
	" --badpixelmask=<FITS_IMAGE>"	
	" --masterfabperot=<FITS_IMAGE>"
	" --method=<INT_OPTVALUE>"    
	" --spectralElementHeight=<FLT_VALUE>"
	" --referenceLineWidth=<FLT_VALUE>"
	" --ordernumber=<INT_VALUE>" 
	" --minorder=<INT_VALUE>"
	" --maxorder=<INT_VALUE>"         
	" --binsize=<UNS_VALUE>"
	" --ipDimensions=<\"UNS_VALUE,UNS_VALUE,UNS_VALUE,UNS_VALUE\">"
	" --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>" 
	" --interactive=<BOOL>\n\n"     
	
	" Example: "+string(modulename)+" --geometryfilename=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.geom --masterbias=/opera/calibrations/11AQ14-Jul08/masterbias_OLAPAa_sp2_Normal.fits --masterflat=/opera/calibrations/11AQ14-Jul08/masterflat_OLAPAa_sp2_Normal.fits --mastercomparison=/opera/calibrations/11AQ14-Jul08/mastercomparison_OLAPAa_sp2_Normal.fits --badpixelmask=/opera/config/badpix_olapa-a.fits --ipDimensions=\"28 5 10 5\" --gain=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.gain --referenceLineWidth=2.5  --masterfabperot=/opera/calibrations/11AQ14-Jul08/masterfabperot_OLAPAa_sp2_Normal.fits --spectralElementHeight=1.0 --outputProf=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.prof --binsize=90 -O 23 -P test.eps -F test.dat -S test.gnu \n\n"
	" -h, --help  display help message\n"
	" -v, --verbose,  Turn on message sending\n"
	" -d, --debug,  Turn on debug messages\n"
	" -t, --trace,  Turn on trace messages\n"
	" -o, --outputProf=<PROF_FILE>, Output instrument profile file\n"    
    " -G, --gainfilename=<GAIN_FILE>, Input gain/noise file\n"
    " -g, --geometryfilename=<GEOM_FILE>, Input geometry file\n"
	" -b, --masterbias=<FITS_IMAGE>, Input Master Bias FITS image\n"	
	" -f, --masterflat=<FITS_IMAGE>, Input Master Flat-Field FITS image\n"		
	" -c, --mastercomparison=<FITS_IMAGE>, Input Master Comparison (ThAr) FITS image\n"	
	" -m, --badpixelmask=<FITS_IMAGE>, FITS image for badpixel mask\n"
	" -a, --masterfabperot=<FITS_IMAGE>, Input Master Fabry-Perot FITS image\n"	    
	" -M, --method=<UNS_OPTVALUE>, Method to combine IP measurements from spectral lines\n"
    "                              Available options are method = 1, 2, 3 or 4, where: \n"
    "                              1. weighted mean \n"  
    "                              2. median combine \n"  
    "                              3. polynomial fit \n"  
    "                              4. polynomial fit with median binning (use binsize provided) \n"  
    " -H, --spectralElementHeight=<FLT_VALUE>, Height of spectral element in Y-direction in pixel units\n"
    " -W, --referenceLineWidth=<FLT_VALUE>, Spectral line width for reference in pixel units\n"
	" -O, --ordernumber=<INT_VALUE>, Pick order number to plot IP model (default = all)\n"
	" -N, --minorder=<INT_VALUE>, Define minimum order number\n"
	" -X, --maxorder=<INT_VALUE>, Define maximum order number\n"         
    " -B, --binsize=<UNS_VALUE>, Number of points to bin for IP measurements\n"    
	" -D, --ipDimensions=<\"UNS_VALUE,UNS_VALUE,UNS_VALUE,UNS_VALUE\"> \n"
	" -P, --plotfilename=<EPS_FILE>\n"
	" -F, --datafilename=<DATA_FILE>\n"
	" -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	" -I, --interactive=<BOOL>\n\n";
}

void GenerateInstrumentProfile3DPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned minorderWithIP, unsigned maxorderWithIP, unsigned IPxsize, unsigned IPysize, bool display)
{
    ofstream *fgnu = NULL;
    
    if (!gnuScriptFileName.empty()) {
        remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
        fgnu = new ofstream();
        fgnu->open(gnuScriptFileName.c_str());
    } else {
        exit(EXIT_FAILURE);
    }
    
    *fgnu << "reset" << endl;
    *fgnu << "unset key" << endl;
    *fgnu << "set view 0,0" << endl;
    
    *fgnu << "set palette color" << endl;
    *fgnu << "set palette gamma 2.5" << endl;
    *fgnu << "set pm3d map" << endl;
    *fgnu << "unset ztics" << endl;
    *fgnu << "set cblabel \"flux fraction\"" << endl;
    *fgnu << endl;
    *fgnu << "set xrange[-5:(5*" << IPxsize << "+5)]" << endl;
    *fgnu << "set yrange[-5:(("<< (maxorderWithIP - minorderWithIP) <<"*" << IPysize << "/5)+5)]" << endl;
    
    
    for(unsigned order=minorderWithIP;order<=maxorderWithIP;order++) {
        double xlabelpos = ((double)order-((double)minorderWithIP+(double)(5*floor((order-minorderWithIP)/5))))*(double)IPxsize+1;
        double ylabelpos = (double)floor((order-minorderWithIP)/5)*(double)IPysize+1;
        *fgnu << "set label \"" << order << "\" at " << xlabelpos << "," << ylabelpos << " front font \"Helvetica,7\"" << endl;
    }

    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nsplot \"" << dataFileName << "\" u (($1-("<<minorderWithIP <<"+5*int(($1-"<< minorderWithIP <<")/5)))*" << IPxsize << "+" << (float)IPxsize/2 << "+$4):(int(($1-"<<minorderWithIP <<")/5)*" << IPysize << "+(" << (float)IPysize/2 << ")+$5):7 with pm3d" << endl;
        
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
        
        *fgnu << "\nsplot \"" << dataFileName << "\" u (($1-("<<minorderWithIP <<"+5*int(($1-"<< minorderWithIP <<")/5)))*" << IPxsize << "+" << (float)IPxsize/2 << "+$4):(int(($1-"<<minorderWithIP <<")/5)*" << IPysize << "+(" << (float)IPysize/2 << ")+$5):7 with pm3d" << endl;
        
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



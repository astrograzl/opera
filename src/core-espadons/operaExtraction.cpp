/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaExtraction
 Version: 1.0
 Description: Extract spectrum usng various alogorithms.
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

#include <pthread.h>

#include "globaldefines.h"
#include "operaError.h"
#include "core-espadons/operaExtraction.h"

#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaInstrumentProfile.h"		// for operaInstrumentProfile
#include "libraries/GainBiasNoise.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaEspadonsImage.h"			// for imtype_t
#include "libraries/operaFITSProduct.h"	
#include "libraries/operaLib.h"						// for systemf and itos

#include "libraries/operaLibCommon.h"
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaFit.h"	
#include "libraries/operaFFT.h"

#define NOTPROVIDED -999

/*! \file operaExtraction.cpp */

using namespace std;

/*! 
 * operaExtraction
 * \author Eder Martioli
 * \brief Module to extract spectra usng various alogorithms.
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


string inputImage;
string outputSpectraFile;
string spectrumtypename;

string inputgain;
string inputgeom;
string inputprof;
string inputaper;

string masterbias;
string masterflat;
string normalizedflatfile;
string badpixelmask;

int ordernumber = NOTPROVIDED;

int minorder = 22;
bool minorderprovided = false;
int maxorder = 62;
bool maxorderprovided = false;

unsigned backgroundBinsize = 1;
double effectiveApertureFraction = 0.99;

unsigned sigmaclip = 50;
unsigned iterations = 3;
bool onTargetProfile = false;
bool usePolynomialFit = false;
bool removeBackground = false;

bool starplusskymode = false;
bool starplusskyInvertSkyFiber = false;

operaSpectralOrder_t spectralOrderType = RawBeamSpectrum;

unsigned maxthreads = 1;

bool interactive = false;

bool noCrossCorrelation = false;

int debug=0, verbose=0, trace=0, plot=0;

string plotfilename;
string datafilename;
string scriptfilename;

operaSpectralOrderVector spectralOrders;
GainBiasNoise *gainBiasNoise = NULL;

operaFITSImage *bias = NULL;
operaFITSImage *normalizedflat = NULL;
operaFITSImage *badpix = NULL;
operaFITSImage *flat = NULL;
operaFITSImage *object = NULL;

ofstream *fdata = NULL;

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
    
    operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
    
    if (!spectralOrder->gethasGeometry()) {
        if (verbose)
            cout << "operaExtraction: Skipping order number: "<< order << " no geometry." << endl;
    }
    if (!spectralOrder->gethasInstrumentProfile()) {
        if (verbose)
            cout << "operaExtraction: Skipping order number: "<< order << " no instrument profile." << endl;
    }
    if (!spectralOrder->gethasExtractionApertures()) {
        if (verbose)
            cout << "operaExtraction: Skipping order number: "<< order << " no extraction aperture." << endl;
    }
    if (spectralOrder->gethasGeometry() &&
        spectralOrder->gethasInstrumentProfile() &&
        spectralOrder->gethasExtractionApertures()) {
        
        if (verbose)
            cout << "operaExtraction: Processing order number: "<< order << endl;
        
        switch (spectralOrderType) {
            case RawBeamSpectrum:
                spectralOrder->extractRawSpectrum(*object, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction, NULL);
                if(!noCrossCorrelation) {
                    spectralOrder->calculateXCorrBetweenIPandImage(*object,*badpix,NULL);
                }
                break;
            case StandardBeamSpectrum:
                spectralOrder->extractStandardSpectrum(*object, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction, backgroundBinsize, NULL);
                if(!noCrossCorrelation) {
                    spectralOrder->calculateXCorrBetweenIPandImage(*object,*badpix,NULL);
                }
                break;
            case OptimalBeamSpectrum:
                if(noCrossCorrelation) {
                    spectralOrder->extractOptimalSpectrum(*object, *flat, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction, backgroundBinsize,sigmaclip,iterations, onTargetProfile, usePolynomialFit, removeBackground, verbose, FALSE, NULL);
                } else {
                    spectralOrder->extractOptimalSpectrum(*object, *flat, *normalizedflat, *bias, *badpix, *gainBiasNoise, effectiveApertureFraction, backgroundBinsize,sigmaclip,iterations, onTargetProfile, usePolynomialFit, removeBackground, verbose, TRUE, NULL);
                }
                
                break;
            case OperaOptimalBeamSpectrum:
                //spectralOrder->extractOPERAOptimal(object, flat, *bias, *badpix, backgroundBinsize, *gainBiasNoise);
                break;
            default:
                break;
        }
        
        if(starplusskymode) {
            spectralOrder->calculateStarAndSkyElements(starplusskyInvertSkyFiber, NULL);
        }
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
		{"inputImage",1, NULL, 'i'},
		{"outputSpectraFile",1, NULL, 's'},		
		{"inputGainFile",1, NULL, 'e'},        
		{"inputGeometryFile",1, NULL, 'g'},
		{"inputInstrumentProfileFile",1, NULL, 'r'},
		{"inputApertureFile",1, NULL, 'a'},
		{"masterbias",1, NULL, 'b'},	
		{"masterflat",1, NULL, 'f'},		
		{"normalizedflat",1, NULL, 'n'},		
		{"badpixelmask",1, NULL, 'm'},
		{"ordernumber",1, NULL, 'O'},	
		{"minorder",1, NULL, 'M'},
		{"maxorder",1, NULL, 'X'},       
        {"effectiveApertureFraction",1, NULL, 'A'},
        {"backgroundBinsize",1, NULL, 'B'},
		{"sigmaclip",1, NULL, 'C'},
		{"iterations",1, NULL, 'R'},  
		{"onTargetProfile",1, NULL, 'J'},  
        {"usePolynomialFit",1, NULL, 'Y'},
        {"removeBackground",1, NULL, 'G'},
        {"starplusskymode",1, NULL, 'K'},
        {"starplusskyInvertSkyFiber",1, NULL, 'V'},
        {"noCrossCorrelation",1, NULL, 'x'},
		{"spectrumtype",1, NULL, 'T'},	
		{"spectrumtypename",1, NULL, 'N'},	
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},  
		{"interactive",0, NULL, 'I'}, 
		{"maxthreads",1, NULL, 'y'},
		
		{"plot",0, NULL, 'p'},
		{"verbose",0, NULL, 'v'},
		{"debug",0, NULL, 'd'},
		{"trace",0, NULL, 't'},
		{"help",0, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "i:s:e:g:r:a:b:f:n:m:O:M:X:A:B:C:R:J:Y:G:K:V:x:T:N:P:F:S:I:y:v::d::p::t::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'i':
				inputImage = optarg;
				break;   
			case 's':
				outputSpectraFile = optarg;
                break; 				
			case 'e':
				inputgain = optarg;
				break;	                
			case 'g':
				inputgeom = optarg;
				break;				
			case 'r':
				inputprof = optarg;
				break;	
			case 'a':
				inputaper = optarg;
				break;
			case 'b':		// bias
				masterbias = optarg;
				break;  
			case 'f':		// flat
				masterflat = optarg;
				break;  
			case 'n':		// normalized flat
				normalizedflatfile = optarg;
				break;  				
			case 'm':		// badpixelmask
				badpixelmask = optarg;
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
			case 'A':
				effectiveApertureFraction = atof(optarg);
				break;
			case 'B':
				backgroundBinsize = atoi(optarg);
				break;
			case 'C':
				sigmaclip = atoi(optarg);
				break;  
			case 'R':
				iterations = atoi(optarg);
				break;  
			case 'J':
				onTargetProfile = (atoi(optarg)?true:false);
				break;  
			case 'Y':
				usePolynomialFit = (atoi(optarg)?true:false);
				break;
			case 'G':
				removeBackground = (atoi(optarg)?true:false);
				break;
			case 'K':
				starplusskymode = (atoi(optarg)?true:false);
				break;
			case 'V':
				starplusskyInvertSkyFiber = (atoi(optarg)?true:false);
				break;
			case 'x':
				noCrossCorrelation = (atoi(optarg)?true:false);
				break;                                  
			case 'T':		// spectrum type
				spectralOrderType = (operaSpectralOrder_t)atoi(optarg);
				break;
			case 'N':		// spectrum type name for verbose mode
				spectrumtypename = optarg;
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
			case 'y':
                maxthreads = atoi(optarg);
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
	
	/*Start the module here*/
	
	try {
		// we need an image...
		if (inputImage.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a gain file...
		if (inputgain.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}        
		// we need a geometry file...
		if (inputgeom.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		// we need an instrument profile file...
		if (inputprof.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}	        
		// we need a aperture file...
		if (inputaper.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}	
		// we need a master flat...
		if (masterflat.empty()) {
			throw operaException("operaExtraction: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}			
		
		if (verbose) {
			cout << "operaExtraction: inputImage = " << inputImage << endl; 
			cout << "operaExtraction: outputSpectraFile = " << outputSpectraFile << endl;
			cout << "operaExtraction: inputGainFile = " << inputgain << endl;            
			cout << "operaExtraction: inputGeometryFile = " << inputgeom << endl; 			
			cout << "operaExtraction: inputInstrumentProfileFile = " << inputprof << endl;
			cout << "operaExtraction: inputApertureFile = " << inputaper << endl;            
			cout << "operaExtraction: masterbias = " << masterbias << endl; 	
			cout << "operaExtraction: masterflat = " << masterflat << endl; 	
			cout << "operaExtraction: normalizedflatfile = " << normalizedflatfile << endl; 				
			cout << "operaExtraction: badpixelmask = " << badpixelmask << endl;
			cout << "operaExtraction: spectrumtype = " << spectralOrderType << endl;
			cout << "operaExtraction: spectrumtypename = " << spectrumtypename << endl;
            cout << "operaExtraction: effectiveApertureFraction = " << effectiveApertureFraction << endl;
            cout << "operaExtraction: backgroundBinsize = " << backgroundBinsize << endl;
            cout << "operaExtraction: starplusskymode = " << starplusskymode << endl;
            cout << "operaExtraction: starplusskyInvertSkyFiber = " << starplusskyInvertSkyFiber << endl;
            cout << "operaExtraction: sigmaclip = " << sigmaclip << endl;
			cout << "operaExtraction: iterations = " << iterations << endl;
			cout << "operaExtraction: onTargetProfile = " << onTargetProfile << endl;
			cout << "operaExtraction: usePolynomialFit = " << usePolynomialFit << endl;
			cout << "operaExtraction: removeBackground = " << removeBackground << endl;
			cout << "operaExtraction: noCrossCorrelation = " << noCrossCorrelation << endl;
			cout << "operaExtraction: removeBackground = " << removeBackground << endl;
            
            if(ordernumber != NOTPROVIDED) {
                cout << "operaExtraction: ordernumber = " << ordernumber << endl;            
            }   
			cout << "operaExtraction: backgroundBinsize = " << backgroundBinsize << endl;            
            if(plot) {
                cout << "operaExtraction: plotfilename = " << plotfilename << endl;
                cout << "operaExtraction: datafilename = " << datafilename << endl;
                cout << "operaExtraction: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cout << "operaExtraction: interactive = YES" << endl; 
                } else {
                    cout << "operaExtraction: interactive = NO" << endl; 
                }
            }            
            
		}
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }          
        
		flat = new operaFITSImage(masterflat, tfloat, READONLY);
        object = new operaFITSImage(inputImage, tfloat, READONLY);		
		
		if (!masterbias.empty()){              
			bias = new operaFITSImage(masterbias, tfloat, READONLY);
		} else {
            bias = new operaFITSImage(object->getnaxis1(),object->getnaxis2(),tfloat);
            *bias = 0.0;
        }	        
                
		if (!badpixelmask.empty()){              
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(object->getnaxis1(),object->getnaxis2(),tfloat);
            *badpix = 1.0;
        }		
        
		if (!normalizedflatfile.empty()){              
			normalizedflat = new operaFITSImage(normalizedflatfile, tfloat, READONLY);
		} else {
            normalizedflat = new operaFITSImage(object->getnaxis1(),object->getnaxis2(),tfloat);
            *normalizedflat = 1.0;
        }        		

		spectralOrders.ReadSpectralOrders(inputgeom);
        spectralOrders.ReadSpectralOrders(inputaper);
        spectralOrders.ReadSpectralOrders(inputprof);

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
			cout << "operaExtraction: minorder ="<< minorder << " maxorder=" << maxorder << endl;        
        
		/*
		 * Add the order data to the raw spectrum product
		 */
		spectralOrders.readGainNoise(inputgain);
        gainBiasNoise = spectralOrders.getGainBiasNoise();
  
        /* 
         * Uncomment below to introduce an artificial gain/noise
         *   to the second amplifier region. This may be wanted to
         *   simulate the effects of gain/noise in the reduction
         */
        // gainBiasNoise->setGain(1,gainBiasNoise->getGain(0)*50);
        // gainBiasNoise->setNoise(1,gainBiasNoise->getNoise(0));
        
        if (verbose) {
            for(unsigned amp=0;amp<gainBiasNoise->getAmps();amp++) {
                DATASEC_t datasec;
                gainBiasNoise->getDatasec(amp, datasec);
                
                double gain = gainBiasNoise->getGain(amp);
                double noise = gainBiasNoise->getNoise(amp);
                
                cout << "operaExtraction: "
                << " namps=" << gainBiasNoise->getAmps()
                << " amp="   << amp
                << " gain="  << gain
                << " noise=" << noise
                << " DataSec=[" << datasec.x1 << ":" << datasec.x2 <<","<< datasec.y1 << ":" << datasec.y2 << "]"
                << endl;
            }
		}
        
        unsigned long nthreads = maxorder+1;
        threads = (pthread_t *)calloc(nthreads, sizeof(pthread_t*));
        thread_args = (thread_args_t *)calloc(nthreads, sizeof(thread_args_t));
        
        if (maxthreads > 1) {
            processOrders(minorder, maxorder);
        } else {
            for (int order=minorder; order<=maxorder; order++) {
                processSingleOrder(order);
            }
        }
        
        unsigned NumberofBeams = spectralOrders.GetSpectralOrder(minorder)->getnumberOfBeams(); // for plotting
        
        for (int order=minorder; order<=maxorder; order++) {
           if (fdata != NULL) {
               operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
               for(unsigned slitview=0; slitview<2; slitview++) {
                    spectralOrder->printBeamSpectrum(itos(slitview),fdata);
                }
                *fdata << endl;
            }
        }
		// output a spectrum...
		spectralOrders.WriteSpectralOrders(outputSpectraFile, spectralOrderType);
        
		object->operaFITSImageClose();
		
        flat->operaFITSImageClose();
        
        if(bias)
            delete bias;
        if(badpix)
            delete badpix;
        if(normalizedflat)
            delete normalizedflat;
        if(object)
            delete object;
        if(flat)
            delete flat;
        
        if (fdata != NULL) {
            fdata->close();
            if (!scriptfilename.empty()) {
                GenerateExtraction3DSpecPlot(scriptfilename,plotfilename,datafilename, NumberofBeams, interactive);
            }
        }             
	}
	catch (operaException e) {
		cerr << "operaExtraction: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaExtraction: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaExtraction: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -i, --inputImage=<FITS_IMAGE>, Input FITS image to extract spectrum\n"
	"  -s, --outputSpectraFile=<FILE_NAME>, Output file name\n"
	"  -e, --inputGainFile=<GAIN_FILE>, Input noise/gain file\n"
	"  -g, --inputGeometryFile=<GEOM_FILE>, Input geometry file\n"
	"  -p, --inputInstrumentProfileFile=<PROF_FILE>, Input instrument profile file\n"
	"  -a, --inputApertureFile=<APER_FILE>, Input extraction aperture file\n"
	"  -b, --masterbias=<FITS_IMAGE>, FITS image with masterbias\n"
	"  -f, --masterflat=<FITS_IMAGE>, FITS image with masterflat\n"
	"  -n, --normalizedflat=<FITS_IMAGE>, FITS image with normalized flat-field\n"
	"  -m, --badpixelmask=<FITS_IMAGE>, FITS image with badpixel mask\n"
	"  -O, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
	"  -N, --minorder=<INT_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<INT_VALUE>, Define maximum order number\n"
    "  -A, --effectiveApertureFraction=<FLT_VALUE>, Fraction of aperture width to set spectral element size\n"
    "  -B, --backgroundBinsize=<UNS_VALUE>, Number of points to bin for IP measurements\n"
    "  -C, --sigmaclip=<UNS_VALUE>, Variance threshold for optimal extraction\n"
    "  -R, --iterations=<UNS_VALUE>, Number iterations for optimal extraction\n"
    "  -J, --onTargetProfile=<BOOL>, Measure spatial profile on-target instead of using flat-field\n"
    "  -Y, --usePolynomialFit=<BOOL>, Use polynomial instead of median for first measurement of profile\n"
    "  -G, --removeBackground=<BOOL>, Remove background. May be turned Off if target is bright. \n"
    "  -K, --starplusskymode=<BOOL>, Star+sky: main flux is the sum of right beams minus sum of left beams. \n"
    "  -V, --starplusskyInvertSkyFiber=<BOOL>, Star+sky: invert sky fiber (default is beam[0]=star and beam[1]=sky). \n"
    "  -x, --noCrossCorrelation=<BOOL>, Use this flag to turn-off the cross-correlation calculation\n"
	"  -T, --spectrumtype=<UNS_VALUE>, Method for extraction\n"
    "                              Available options are = 5, 6, 7 or 8, where: \n"
    "                              5. Raw Flux Sum (default)\n"
    "                              6. Standard Flux \n"
    "                              7. Optimal Extraction; ref: Horne, K., (1986) & Marsh, T.R., (1989)\n"
    "                              8. OPERA Optimal Extraction\n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputImage=/data/espadons/51Peg-12BQ10-Dec01//1599047o.fits --badpixelmask=/Users/edermartioli/opera-1.0//config/badpix_olapa-a.fits.fz --masterbias=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/masterbias_OLAPAa_pol_Normal.fits.fz --masterflat=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/masterflat_OLAPAa_pol_Normal.fits.fz  --inputInstrumentProfileFile=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/OLAPAa_pol_Normal.prof.gz --inputGeometryFile=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/OLAPAa_pol_Normal.geom.gz --inputApertureFile=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/OLAPAa_pol_Normal.aper.gz --inputGainFile=/Users/edermartioli//opera//calibrations/51Peg-12BQ10-Dec01/OLAPAa_pol_Normal.gain.gz --backgroundBinsize=120 --sigmaclip=50 --onTargetProfile=1 --iterations=3 --ordernumber=-999 --usePolynomialFit=0 --spectrumtype=7 --spectrumtypename=OptimalBeamSpectrum  --outputSpectraFile=1599047.e.gz --removeBackground=0 --scriptfilename=optimalSpec.gnu --datafilename=optimalSpec.dat --plotfilename=optimalSpec.eps --minorder=28 --maxorder=32 -v -t\n\n"
    ;
}

/*
 * Generate 3D plot for spectra of atlas + comparison + identified lines
 */
void GenerateExtraction3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string datafilename, unsigned numberOfBeams, bool display) {
    
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
    
    *fgnu << "set palette gray" << endl;
    *fgnu << "set palette gamma 2.0" << endl;
    *fgnu << "set pm3d map" << endl;
    *fgnu << "unset ztics" << endl;
    *fgnu << "set cblabel \"flux\"" << endl;

    *fgnu << "set xrange[-200:*]" << endl;
   
    *fgnu << "\nset xlabel \"distance (pixels)\"" << endl;
    *fgnu << "set ylabel \"order number\"" << endl;
    
    unsigned fluxColumnForFirstBeam = 13;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;
        
        *fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.4*$1 - 0.8 + 0.25 + 0.125):" << fluxColumnForFirstBeam <<" w pm3d";
        for(unsigned beam=1; beam<numberOfBeams; beam++) {
            unsigned fluxcol = fluxColumnForFirstBeam + 4*beam;
            *fgnu << ",\"\" u 6:($2 + 0.4*$1 - 0.35 + 0.25 + 0.125):" << fluxcol << " w pm3d";
        }
        *fgnu << endl;
    
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
        *fgnu << endl;
        
        *fgnu << "splot \"" << datafilename << "\" u 6:($2 + 0.4*$1 - 0.8 + 0.25 + 0.125):" << fluxColumnForFirstBeam <<" w pm3d";
        for(unsigned beam=1; beam<numberOfBeams; beam++) {
            unsigned fluxcol = fluxColumnForFirstBeam + 4*beam;
            *fgnu << ",\"\" u 6:($2 + 0.4*$1 - 0.35 + 0.25 + 0.125):" << fluxcol << " w pm3d";
        }
        *fgnu << endl;
        
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

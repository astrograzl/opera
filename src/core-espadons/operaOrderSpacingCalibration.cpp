/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaOderSpacingCalibration
 Version: 1.0
 Description: this module performs the measurement of spacing between orders
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
#include "core-espadons/operaOrderSpacingCalibration.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/Polynomial.h"				// for Polynomial
#include "libraries/operaInstrumentProfile.h"	// for operaInstrumentProfile
#include "libraries/operaSpectralElements.h"	// for operaSpectralElement
#include "libraries/operaGeometry.h"			// for operaGeometry
#include "libraries/operaSpectralOrder.h"		// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"	// for operaSpectralOrderVector
#include "libraries/GainBiasNoise.h"

#include "libraries/operaLibCommon.h"
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"					// for MAXORDERS
#include "libraries/operaFit.h"	
#include "libraries/operaFFT.h"
#include "libraries/operaParameterAccess.h"
#include "libraries/operaConfigurationAccess.h"

/*! \file operaOrderSpacingCalibration.cpp */

using namespace std;

/*! 
 * operaOrderSpacingCalibration
 * \author Doug Teeple
 * \brief Order spacing measurements.
 * \arg argc
 * \arg argv
 * \note --orderspacingoutput=...
 * \note --masterbias=...
 * \note --masterflat=...
 * \note --badpixelmask=...
 * \note --inputGainFile=...
 * \note --subformat=...
 * \note --detectionMethod=...
 * \note --FFTfilter=...
 * \note --aperture=...
 * \note --numberOfsamples=...
 * \note --plotfilename=...
 * \note --datafilename=...
 * \note --scriptfilename=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;

	string orderspacingoutput;
	string masterbias; 
	string masterflat; 
	string badpixelmask;
	string inputGainFile;
	struct subformat {
		unsigned x1, x2;
		unsigned y1, y2;
	} subformat = {8, 2040, 3, 4600 };
	bool FFTfilter = false;
	int aperture = 20;	
	int detectionMethod = 1; // 1. Gaussian, 2. IP, 3. Top-hat
	
    unsigned nsamples = 30;
    unsigned sampleCenterPosition = 1;  // This position is with respect to the dispersion direction (rows for Espadons)
    
    double gain = 1;
    double noise = 5;
    
    string plotfilename;
	string datafilename;
	string scriptfilename;
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
	struct option longopts[] = {
		{"orderspacingoutput",1, NULL, 'o'},
		{"masterbias",1, NULL, 'b'},
		{"masterflat",1, NULL, 'f'},
		{"badpixelmask",1, NULL, 'm'},
		{"defaultgain",1, NULL, 'a'},
		{"defaultoise",1, NULL, 'n'},
		{"inputGainFile",1, NULL, 'g'},
		{"subformat",1, NULL, 's'},
		{"detectionMethod",1, NULL, 'M'},	
		{"FFTfilter",1, NULL, 'R'},			
		{"aperture",1, NULL, 'A'},
        {"numberOfsamples",1, NULL, 'N'},
		{"sampleCenterPosition",1, NULL, 'Y'},  // This position is with respect to the dispersion direction (rows for Espadons)
        
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},

		{"interactive",	optional_argument, NULL, 'I'},         
		{"plot",		optional_argument, NULL, 'p'},
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "o:b:f:m:g:a:n:s:M:R:A:N:Y:P:F:S:I::p::v::d::t::h",
							 longopts, NULL))  != -1)
	{
		
		switch(opt) 
		{
			case 'o':		// output
				orderspacingoutput = optarg;
				break;
			case 'b':		// output
				masterbias = optarg;
				break;
			case 'f':		// masterflat
				masterflat = optarg;
				break;
			case 'm':		// badpixelmask
				badpixelmask = optarg;
				break;            
			case 'g':		// gain / noise / bias
                inputGainFile = optarg;
                break;
			case 'a':		// gain
                gain = atof(optarg);
                break;
			case 'n':		// noise
                noise = atof(optarg);
                break;
			case 's':		// subformat
				if (strlen(optarg))
					sscanf(optarg, "%u %u %u %u", &subformat.x1, &subformat.x2, &subformat.y1, &subformat.y2);
				break;                                                    
			case 'M':
				detectionMethod = atoi(optarg); // 1. Gaussian, 2. IP, 3. Top-hat	
				break;            
			case 'R':
				FFTfilter = (atoi(optarg)?true:false); 
				break;				                    
			case 'A':		// aperture in pixels
				aperture = atoi(optarg);
				break;     
			case 'N':
				nsamples = atoi(optarg);
				break;
            case 'Y':
				sampleCenterPosition = atoi(optarg);
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

		// we need a masterflat...
		if (masterflat.empty()) {
			throw operaException("operaOrderSpacingCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need some place to put data...
		if (orderspacingoutput.empty()) {
			throw operaException("operaOrderSpacingCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		if (verbose) {
			cout << "operaOrderSpacingCalibration: orderspacingoutput = " << orderspacingoutput << endl;
			cout << "operaOrderSpacingCalibration: masterbias = " << masterbias << endl;            
			cout << "operaOrderSpacingCalibration: masterflat = " << masterflat << endl; 	
			cout << "operaOrderSpacingCalibration: badpixelmask = " << badpixelmask << endl;
			cout << "operaOrderSpacingCalibration: inputGainFile = " << inputGainFile << endl;            
			cout << "operaOrderSpacingCalibration: subformat = " << subformat.x1 << " " << subformat.x2 << " "<< subformat.y1  << " "<< subformat.y2 << "\n";
			cout << "operaOrderSpacingCalibration: detectionMethod = " << detectionMethod << endl;
			cout << "operaOrderSpacingCalibration: FFTfilter = " << FFTfilter << endl;
			cout << "operaOrderSpacingCalibration: aperture = " << aperture << endl;
			cout << "operaOrderSpacingCalibration: nsamples = " << nsamples << endl;
            if(plot) {
                cout << "operaOrderSpacingCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaOrderSpacingCalibration: datafilename = " << datafilename << endl;
                cout << "operaOrderSpacingCalibration: scriptfilename = " << scriptfilename << endl;
                if(interactive) {
                    cout << "operaOrderSpacingCalibration: interactive = YES" << endl;
                } else {
                    cout << "operaOrderSpacingCalibration: interactive = NO" << endl;
                }
            }            
		}
        
        ofstream *fdata = NULL;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());
        }
        
		/*
		 * open input images and load data into an in-memory image
		 */
		
		unsigned x1 = subformat.x1;
		unsigned y1 = subformat.y1;
		unsigned nx = subformat.x2 - subformat.x1;
		unsigned ny = subformat.y2 - subformat.y1;
                
        operaFITSImage flat(masterflat, tfloat, READONLY);

        operaFITSImage *bias = NULL;
        
		if (!masterbias.empty()){
			bias = new operaFITSImage(masterbias, tfloat, READONLY);
		} else {
            bias = new operaFITSImage(flat.getnaxis1(),flat.getnaxis2(),tfloat);
            *bias = 0.0;
        }
            
		flat -= (*bias);			// remove bias from masterflat
        if(bias) {
            delete bias;
		}
        operaFITSImage *badpix = NULL;
        
		if (!badpixelmask.empty()){
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat.getnaxis1(),flat.getnaxis2(),tfloat);
            *badpix = 1.0;
        }
		
		if (verbose) {
			cout << "operaOrderSpacingCalibration: x1,y1,nx,ny = " << x1 << ' ' << y1 << ' ' << nx  << ' ' << ny << '\n'; 
		}
		
        float slit = (float)aperture;
        
		operaSpectralOrderVector spectralOrders(MAXORDERS, ny, ny, 0);        
		
		/*
		 * read gain and noise from input file
		 */
        unsigned amp = 0;
        if (!inputGainFile.empty()) {
            spectralOrders.readGainNoise(inputGainFile);
            gain = spectralOrders.getGainBiasNoise()->getGain(amp);
            noise = spectralOrders.getGainBiasNoise()->getNoise(amp);
		}
		if (verbose)
			cout << "operaOrderSpacingCalibration: gain="<< gain << " noise=" << noise << endl;
        
        spectralOrders.fitOrderSpacingPolynomial(flat, *badpix, slit, nsamples, sampleCenterPosition, detectionMethod, FFTfilter, (float)gain, (float)noise, subformat.x1, subformat.x2, subformat.y1, subformat.y2, fdata);

        if (fdata != NULL) {
            fdata->close();
            if (!scriptfilename.empty()) {
                GenerateOrderSpacingPlot(scriptfilename,plotfilename,datafilename, interactive);
            }
        }
        
		//
		// now flush out the order spacing output
		//
		if (!orderspacingoutput.empty()) {
			spectralOrders.WriteSpectralOrders(orderspacingoutput, Orderspacing);
		}
        if(badpix)
            delete badpix;
        
		flat.operaFITSImageClose();
	}
	catch (operaException e) {
		cerr << "operaOrderSpacingCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaOrderSpacingCalibration: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaOrderSpacingCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --orderspacingoutput=<ORDP_FILE>"
	" --masterbias=<FITS_IMAGE>"
	" --masterflat=<FITS_IMAGE>"
	" --badpixelmask=<FITS_IMAGE>"
    " --inputGainFile=<GAIN_FILE>"
    " --subformat=<\"UNS_VALUE UNS_VALUE UNS_VALUE UNS_VALUE\">"
    " --detectionMethod=<UNS_VALUE>"
    " --FFTfilter=<BOOL>"
    " --aperture=<UNS_VALUE>"
    " --numberOfsamples=<UNS_VALUE>"
    " --sampleCenterPosition=<UNS_VALUE>"    
    " --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --masterbias=/Users/edermartioli/opera//calibrations/PolarData/masterbias_OLAPAa_pol_Normal.fits.fz --masterflat=/Users/edermartioli//opera//calibrations/PolarData/masterflat_OLAPAa_pol_Normal.fits.fz --badpixelmask=/Users/edermartioli/opera-1.0//config/badpix_olapa-a.fits.fz --inputGainFile=/Users/edermartioli//opera//calibrations/PolarData/OLAPAa_pol_Normal.gain.gz --subformat=\"8 2040 3 4600\" --aperture=26 --detectionMethod=2 --FFTfilter=0 --numberOfsamples=30  --sampleCenterPosition=2300 --orderspacingoutput=OLAPAa_pol_Normal.ordp.gz --datafilename=spacing.dat --scriptfilename=spacing.gnu --plotfilename=spacing.eps -v\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -o, --orderspacingoutput=<FILE_NAME>, Order spacing output file name\n"    
	"  -b, --masterbias=<FITS_IMAGE>, FITS image with masterbias\n"
	"  -f, --masterflat=<FITS_IMAGE>, FITS image with masterflat\n"
	"  -m, --badpixelmask=<FITS_IMAGE>, FITS image with badpixel mask\n"
	"  -g, --inputGainFile=<GAIN_FILE>, Input noise/gain/bias file\n"
    "  -s, --subformat=<\"UNS_VALUE UNS_VALUE UNS_VALUE UNS_VALUE\">, Image subformat to be inspected\n"
	"  -M, --detectionMethod=<UNS_VALUE>, Method for detecting orders\n"
    "                              Available options are = 1, 2, and  3, where: \n"
    "                              1. Gaussian (default)\n"
    "                              2. IP \n"
    "                              3. Top-hat\n"
    "  -R, --FFTfilter=<BOOL>, Activate Fourier smoothing filter\n"
    "  -A, --aperture=<UNS_VALUE>, Aperture size in pixel units\n"
    "  -N, --numberOfsamples=<UNS_VALUE>, Number of row samples for detecting orders\n"
    "  -Y, --sampleCenterPosition=<UNS_VALUE>, Detector position to center samples. Position along the dispersion direction (rows for Espadons)\n"    
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateOrderSpacingPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display)
{
    FILE *fgnu;
    remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	
    fgnu = fopen(gnuScriptFileName.c_str(),"w");
    
    fprintf(fgnu,"reset\n");
    fprintf(fgnu,"unset key\n");
    fprintf(fgnu,"\nset xlabel \"order center (col pixels)\"\n");
    fprintf(fgnu,"set ylabel \"order separation (pixels)\"\n");
    
    if(!outputPlotEPSFileName.empty()) {
        fprintf(fgnu,"\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
        fprintf(fgnu,"set output \"%s\"\n",outputPlotEPSFileName.c_str());
        
        fprintf(fgnu,"\nplot \"%s\" u 1:2:3 with yerr,\"%s\" u 1:4 w l\n",dataFileName.c_str(),dataFileName.c_str());
        
        if (display) {
            fprintf(fgnu,"\nset terminal x11\n");
            fprintf(fgnu,"set output\n");
            fprintf(fgnu,"replot\n");
        } else {
            fprintf(fgnu,"\n#set terminal x11\n");
            fprintf(fgnu,"#set output\n");
            fprintf(fgnu,"#replot\n"); 
        }
    } else {
        fprintf(fgnu,"\nplot \"%s\" u 1:2:3 with yerr,\"%s\" u 1:4 w l\n",dataFileName.c_str(),dataFileName.c_str());
        
        fprintf(fgnu,"\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14\n");
        fprintf(fgnu,"#set output \"outputPlotEPSFileName.eps\"\n");
        fprintf(fgnu,"#replot\n");
        fprintf(fgnu,"#set terminal x11\n");
        fprintf(fgnu,"#set output\n");
    }
    
    fclose(fgnu);
    
    if (display) {
        systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    } else {
        if(!outputPlotEPSFileName.empty())
            systemf("gnuplot %s",gnuScriptFileName.c_str());        
    }
}


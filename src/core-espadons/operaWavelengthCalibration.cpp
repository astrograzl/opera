/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaWavelengthCalibration
 Version: 1.0
 Description: Wavelength Calibration 
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
#include <stdarg.h>
#include <getopt.h>
#include <math.h>									// for pow

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "globaldefines.h"
#include "operaError.h"

#include "libraries/operaException.h"
#include "libraries/operaFITSProduct.h"
#include "libraries/operaSpectralElements.h"
#include "libraries/operaSpectralFeature.h"
#include "libraries/operaSpectralLines.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralOrderVector.h"
#include "libraries/GainBiasNoise.h"
#include "libraries/Gaussian.h"
#include "libraries/Polynomial.h"						// for Polynomial

#include "libraries/operaSpectralTools.h"

#include "libraries/operaLibCommon.h"					// for doubleValue_t
#include "libraries/operaLib.h"							// for itos
#include "libraries/operaMath.h"						// for LengthofPolynomial
#include "libraries/operaCCD.h"							// for MAXORDERS
#include "libraries/operaFFT.h"							// for operaXCorrelation
#include "libraries/gzstream.h"							// for gzstream - read compressed reference spectra

#include "core-espadons/operaWavelengthCalibration.h"

#define DEBUG false


/*! \file operaWavelengthCalibration.cpp */
/*! \ingroup core */

using namespace std;

int debug=0, verbose=0, trace=0, plot=0;

/*
 * the reference atlas spectrum
 */
static unsigned nPointsInAtlasSpectrum = 0;
static double atlasWavelength[MAXFULLREFWAVELENGTHS];
static double atlasIntensity[MAXFULLREFWAVELENGTHS];   
static double atlasVariance[MAXFULLREFWAVELENGTHS];   

/*
 * the reference atlas spectral lines
 */
static unsigned thatlaslines = 0;
static double thAtlasWavelength[MAXREFWAVELENGTHS];
static double thAtlasIntensity[MAXREFWAVELENGTHS]; 

/*
 
 Below it follows in a few words a 1st-pass for the wavelength calibration
 algorithm.
 
 For each spectral order do the following steps:
 
 1. Read ThAr raw spectrum; intensity versus distance in pixel units:
 I(d) vs. d
 
 2. Measure total distance "D" (in pixel units) covered by the order. This
 is given by the line integral of the polynomial that describes the center
 of the order.
 
 3. Read wavelength range covered by the order: wl0,wlf
 
 4. Calculate first order solution:
 wl = f(d), where f(d) as first order is given by
 
 f(d) = wl0 + ((wlf - wl0)/D)*d
 
 assuming f(d=0) = wl0.
 
 5. Read ThAr atlas of spectral lines within the range covered by the order
 [wl0:wlf]. The atlas consists of line wavelength (l_wl), error (l_wlerr),
 and line relative intensity (l_i).
 
 7. Once we have a table of l_i, l_d, and l_derr, then we can calculate the
 maximum cross-correlation between this and the atlas data to identify the
 lines. The identification usually doesn't go one-by-one, so we will end up
 having to do some cleaning for either the undetected or over-detected
 lines.
 
 8. Now one can use the table (l_d+/-l_derr) versus (l_wl+/-l_wlerr) to
 find the wavelength solution by fitting a polynomial to these data.
 
 note that the polynomial should be an update to the first-order solution.
 
 The update is intended for two reasons:
 
 First because I have experienced before that the higher order terms are so
 small when compared to the first order that the fitting routine can get in
 trouble to find good solutions.
 
 Another reason is that we want to use our first solution to exclude
 outliers and then run the fitting again as many times as it gets to give
 the best solution. So, we will have to update our solution as we get our
 dataset cleaner or as we gather more information.
 
 */

/* 
 * operaWavelengthCalibration
 * \author Doug Teeple
 * \brief wavelength calibration.
 * \arg argc
 * \arg argv
 * \note --outputWave=...
 * \note --atlas_lines=...
 * \note --thcal=...
 * \note --geom=...
 * \note --wlcal_initialguess=...
 * \note --binsize=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

#define NOTPROVIDED -999

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cerr <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdthp]" +
	" --atlas_lines=<ATLASLINES_FILE>"
	" --atlas_spectrum=<ATLASSPEC_FILE>"
	" --uncalibrated_lines=<RAWLINES_FILE>"
	" --uncalibrated_spectrum=<RAWSPEC_FILE>"
	" --inputGeomFile=<GEOM_FILE>"
	" --inputWaveFile=<WAVE_FILE>"
	" --outputWaveFile=<WAVE_FILE>"
	" --parseSolution=<BOOL>"
	" --normalizeUncalibratedSpectrum=<BOOL>"
	" --normalizationBinSize=<UNS_VALUE>"
	" --LocalMaxFilterWidth=<DBL_VALUE>"    
    " --ParRangeSizeInPerCent=<FLT_VALUE>"
	" --NpointsPerPar=<UNS_VALUE>"
	" --maxNIter=<UNS_VALUE>"
	" --minNumberOfLines=<UNS_VALUE>"
	" --maxorderofpolynomial=<UNS_VALUE>"
	" --dampingFactor=<FLT_VALUE>"
	" --initialAcceptableMismatch=<FLT_VALUE>"
	" --nsigclip=<FLT_VALUE>"
	" --DetectionThreshold=<FLT_VALUE>"
	" --MinPeakDepth=<FLT_VALUE>"
    " --ordernumber=<INT_VALUE"
	" --minorder=<INT_VALUE>"
	" --maxorder=<INT_VALUE>"
    " --ordersplotfilename=<EPS_FILE>"
    " --specplotfilename=<EPS_FILE>"
    " --ordersscriptfilename=<GNUPLOT_FILE>"
    " --specscriptfilename=<GNUPLOT_FILE>"
    " --ordersdatafilename=<DATA_FILE>"
    " --atlasdatafilename=<DATA_FILE>"
    " --compdatafilename=<DATA_FILE>"
    " --linesdatafilename=<DATA_FILE>"
	" --generate3DPlot=<BOOL>"
    " --subtractCentralWavelength=<BOOL>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --atlas_lines=/Users/edermartioli/opera-1.0//config/thar_MM201006.dat.gz --atlas_spectrum=/Users/edermartioli/opera-1.0//config/LovisPepe_ThArAtlas.dat.gz --uncalibrated_spectrum=/Users/edermartioli//opera//spectra/StarOnlyData/th_OLAPAa_sp2_Normal.e.gz --uncalibrated_linewidth=1.25 --inputGeomFile=/Users/edermartioli//opera//calibrations/StarOnlyData/OLAPAa_sp2_Normal.geom.gz --inputWaveFile=/Users/edermartioli/opera-1.0//config/wcal_ref.dat.gz --ordernumber=-999 --parseSolution=0 --ParRangeSizeInPerCent=0.1 --NpointsPerPar=300 --maxNIter=40 --minNumberOfLines=50 --maxorderofpolynomial=4 --dampingFactor=0.85 --initialAcceptableMismatch=1.0 --nsigclip=3  --outputWaveFile=/Users/edermartioli//opera//calibrations/StarOnlyData/OLAPAa_sp2_Normal.wcar.gz --ordersplotfilename=waveorders.eps --specplotfilename=wavespec.eps --ordersscriptfilename=waveorders.gnu --specscriptfilename=wavespec.gnu --ordersdatafilename=waveorders.dat --atlasdatafilename=wavespecatlas.dat --compdatafilename=wavespeccomp.dat --linesdatafilename=wavespeclines.dat --generate3DPlot=0 --subtractCentralWavelength=0 --minorder=24 --maxorder=46 -v -t\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -a, --atlas_lines=<ATLASLINES_FILE>, File containing the atlas of reference lines\n"
	"  -s, --atlas_spectrum=<ATLASSPEC_FILE>, File containing the spectrum of reference atlas\n"
	"  -l, --uncalibrated_lines=<RAWLINES_FILE>, File containing the uncalibrated raw lines\n"
	"  -u, --uncalibrated_spectrum=<RAWSPEC_FILE>, File containing the uncalibrated raw spectrum\n"
	"  -g, --inputGeomFile=<GEOM_FILE>, Input geometry calibration file\n"
	"  -r, --inputWaveFile=<WAVE_FILE>, Input wavelength calibration file (initial guess)\n"
	"  -o, --outputWaveFile=<WAVE_FILE>, Output wavelength calibration file to store final solution\n"
	"  -e, --parseSolution=<BOOL>, Parse parameters to search solution - use only when initial guess is poor\n"
	"  -n, --normalizeUncalibratedSpectrum=<BOOL>, Normalize uncalibrated spectrum\n"
	"  -B, --normalizationBinSize=<UNS_VALUE>, Binsize to be used for normalization\n"
	"  -B, --LocalMaxFilterWidth=<DBL_VALUE>, Parameter (in units of line width) to set a window filter to guarantee a line is not detected twice\n"
    "  -R, --ParRangeSizeInPerCent=<FLT_VALUE>, Parameter to define the range within which a coefficient will be changed to search the solution\n"
	"  -T, --NpointsPerPar=<UNS_VALUE>, Parameter to define the number of times a coefficient will be changed\n"
	"  -M, --maxNIter=<UNS_VALUE>, Maximum number of iterations for shrinking the acceptable mismatch\n"
	"  -L, --minNumberOfLines=<UNS_VALUE>, Minimum number of lines to keep when shrinking the acceptable mismatch difference between atlas and comparison\n"
	"  -Y, --maxorderofpolynomial=<UNS_VALUE>, Maximum degree of polynomial for wavelength solution\n"
	"  -D, --dampingFactor=<FLT_VALUE>, Damping factor for shrinking the size of the acceptableMismatch on each iteration.  This factor may be set between 0 and 1\n"
	"  -A, --initialAcceptableMismatch=<FLT_VALUE>, Initial acceptable mismatch difference between atlas and comparison lines (in units of line width)\n"
	"  -C, --nsigclip=<FLT_VALUE>, Threshold (in units of rms) for clipping lines.\n"
	"  -y  --DetectionThreshold=<FLT_VALUE>, threshold to regulate the sensitivity of line detection. Must be between 0 and 1."
	"  -z  --MinPeakDepth=<FLT_VALUE>, limit that also regulates the sensitity of line detection in units of noise."
    "  -O, --ordernumber=<INT_VALUE>, Absolute order number to extract (default=all)\n"
	"  -N, --minorder=<INT_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<INT_VALUE>, Define maximum order number\n"
	"  -P, --ordersplotfilename=<EPS_FILE>\n"
	"  -Q, --specplotfilename=<EPS_FILE>\n"
	"  -F, --ordersdatafilename=<DATA_FILE>\n"
	"  -G, --atlasdatafilename=<DATA_FILE>\n"
	"  -H, --compdatafilename=<DATA_FILE>\n"
	"  -J, --linesdatafilename=<DATA_FILE>\n"
	"  -S, --ordersscriptfilename=<GNUPLOT_FILE>\n"
	"  -U, --specscriptfilename=<GNUPLOT_FILE>\n"
	"  -E, --generate3DPlot=<BOOL>, Choose a 3D plot of the spectra instead of 2D\n"
	"  -K, --subtractCentralWavelength=<BOOL>,Choose to subtract order central wavelength for plot\n"
    "  -I, --interactive=<BOOL>\n\n";
}

int main(int argc, char *argv[])
{
	int opt;
	
	string outputWave; 
    
	string atlas_lines; 
    string atlas_spectrum;
	
    string uncalibrated_lines;
    string uncalibrated_spectrum;
    
    double uncalibrated_linewidth = 1.5; // default is 1.5 pixels
    
	string geometryfilename; 
	string wlcal_initialguess; 
    
    string ordersplotfilename;
    string specplotfilename;
	
    string ordersscriptfilename;
    string specscriptfilename;
    
    string atlasdatafilename;
	string compdatafilename;
	string linesdatafilename;
	string ordersdatafilename;
    
    bool generate3DPlot = FALSE;
    bool subtractCentralWavelength = TRUE;
    
    bool parseSolution = false;
    bool normalizeUncalibratedSpectrum = false;
    unsigned normalizationBinSize = 150;
    
	bool interactive = false;    
	
    unsigned minorder = 0;
    bool minorderprovided = false;
    unsigned maxorder = MAXORDERS;
    bool maxorderprovided = false; 
    
    int ordernumber = NOTPROVIDED;	
	
    /*
     * The parameters below we don't know yet whether they will be input
     */
    double DetectionThreshold = 0.05;    // threshold to regulate the sensitivity of line detection. Must be between 0 and 1.
    double LocalMaxFilterWidth = 3.0;    // parameter to set a window filter to guarantee a line is not detected twice. It's in units of line width
    double MinPeakDepth = 0.25;           // limit that also regulates the sensitity of line detection in units of noise.
    
    double ParRangeSizeInPerCent = 0.1;  // define the range within which a coefficient will be changed to calcuate the x-correlation
    unsigned NpointsPerPar = 300;        // define the number of times a coefficient is changed
    unsigned maxpolyorder = 2;           // define the coefficients to be changed in order to search for maximum x-correlation
    unsigned nIterforXcorr = 2;          // number of iterations to search for maximum x-correlation
    
    unsigned maxNIter = 30;                 // maximum number of iterations for shrinking the acceptable mismatch
    unsigned minNumberOfLines = 40;         // minimum number of lines to stop shrinking the acceptable mismatch difference between atlas and comparison
    unsigned maxorderofpolynomial = 4;      // maximum degree of polynomial for wavelength solution
    double dampingFactor = 0.90;            // Damping factor to shrink the size of the quantity acceptableMismatch on each iteration.  This factor may be set between 0 to 1.
    double initialAcceptableMismatch = 1.0; // initial acceptable mismatch difference between atlas and comparison lines. Used for identification. In units of line width.
    double nsigclip = 3.0;						// Threshold (in units of rms) for clipping lines.
    
	struct option longopts[] = {
		{"atlas_lines",1, NULL, 'a'},
		{"atlas_spectrum",1, NULL, 's'},        
		{"uncalibrated_lines",1, NULL, 'l'},					// operaExtractSpactralLines does this for us
		{"uncalibrated_spectrum",1, NULL, 'u'},					// one may use the spectrum     
		{"uncalibrated_linewidth",1, NULL, 'w'},				// need this if using input spectrum  
		{"inputGeomFile",		1, NULL, 'g'},	
		{"inputWaveFile",		1, NULL, 'r'},				// stores polynomials for each order to get a first guess as to the wavelengths
		{"outputWaveFile",		1, NULL, 'o'},
		{"parseSolution",1, NULL, 'e'},                         // need this if running with poor initial guess  
		{"normalizeUncalibratedSpectrum",1, NULL, 'n'},         // normalize uncalibrated spectrum
		{"normalizationBinSize",1, NULL, 'B'},         //  binsize to be used for normalization
		{"LocalMaxFilterWidth",1, NULL, 'W'},         // parameter to set a window filter to guarantee a line is not detected twice. It's in units of line width
		{"ParRangeSizeInPerCent",1, NULL, 'R'},         // define the range within which a coefficient will be changed to calcuate the x-correlation
		{"NpointsPerPar",1, NULL, 'T'},                 // define the number of times a coefficient is changed
		{"maxNIter",1, NULL, 'M'},                      // maximum number of iterations for shrinking the acceptable mismatch
		{"minNumberOfLines",1, NULL, 'L'},              // minimum number of lines to stop shrinking the acceptable mismatch difference between atlas and comparison
		{"maxorderofpolynomial",1, NULL, 'Y'},          // maximum degree of polynomial for wavelength solution
		{"dampingFactor",1, NULL, 'D'},                 // Damping factor to shrink the size of the quantity acceptableMismatch on each iteration.  This factor may be set between 0 to 1.
		{"initialAcceptableMismatch",1, NULL, 'A'},     // initial acceptable mismatch difference between atlas and comparison lines. Used for identification. In units of line width.
		{"nsigclip",1, NULL, 'C'},                      // Threshold (in units of rms) for clipping lines.		{"initialAcceptableMismatch",1, NULL, 'A'},     // initial acceptable mismatch difference between atlas and comparison lines. Used for identification. In units of line width.
		{"DetectionThreshold",1, NULL, 'y'},
		{"MinPeakDepth",1, NULL, 'z'},
        {"ordernumber",1, NULL, 'O'},
		{"minorder",1, NULL, 'N'},
		{"maxorder",1, NULL, 'X'},

        {"ordersplotfilename",1, NULL, 'P'},
        {"specplotfilename",1, NULL, 'Q'},
        {"ordersscriptfilename",1, NULL, 'S'},
        {"specscriptfilename",1, NULL, 'U'},
        {"ordersdatafilename",1, NULL, 'F'},
        {"atlasdatafilename",1, NULL, 'G'},
        {"compdatafilename",1, NULL, 'H'},
        {"linesdatafilename",1, NULL, 'J'},
        {"generate3DPlot",1, NULL, 'E'},
        {"subtractCentralWavelength",1, NULL, 'K'},
        
		{"interactive",0, NULL, 'I'},
		{"plot",		optional_argument, NULL, 'p'},       
		{"verbose",		optional_argument, NULL, 'v'},
		{"debug",		optional_argument, NULL, 'd'},
		{"trace",		optional_argument, NULL, 't'},
		{"help",		no_argument, NULL, 'h'},
		{0,0,0,0}};

	while((opt = getopt_long(argc, argv, "a:s:l:u:w:g:r:o:e:n:B:W:R:T:M:L:Y:D:A:C:y:z:O:N:X:P:Q:S:U:E:K:F:G:H:J:I:p::v::d::t::h",
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{
			case 'a':		// atlas_lines
				atlas_lines = optarg;
				break;
			case 's':		// atlas full spectrum
				atlas_spectrum = optarg;
				break;                
			case 'l':		// uncalibrated_lines
				uncalibrated_lines = optarg;
				break;    
			case 'u':		// uncalibrated_spectrum
				uncalibrated_spectrum = optarg;
				break;          
			case 'w':		// uncalibrated_linewidth: necessary only if input spectrum
				uncalibrated_linewidth = atof(optarg);
				break;                 
			case 'g':		// geometryfilename
				geometryfilename = optarg;
				break; 
			case 'r':		// wavelngth calibration reference polynomials
				wlcal_initialguess = optarg;
				break;  
			case 'o':		// output
				outputWave = optarg;
				break;
			case 'e':		// parseSolution: necessary if using a poor input initial guess
				parseSolution = atoi(optarg) == 1;
				break;
			case 'n':		// normalizeUncalibratedSpectrum: normalize uncalibrated input spectrum
				normalizeUncalibratedSpectrum = atoi(optarg) == 1;
				break;
			case 'B':		// normalizationBinSize: binsize for normalization
				normalizationBinSize = atoi(optarg);
				break;
			case 'W':		// parameter to set a window filter to guarantee a line is not detected twice. It's in units of line width
				LocalMaxFilterWidth = atof(optarg);
				break;
			case 'R':
				ParRangeSizeInPerCent = atof(optarg);
				break;
			case 'T':
				NpointsPerPar = atoi(optarg);
				break;
			case 'M':
				maxNIter = atoi(optarg);
				break;
            case 'L':
				minNumberOfLines = atoi(optarg);
				break;
			case 'Y':
				maxorderofpolynomial = atoi(optarg);
				break;
			case 'D':
				dampingFactor = atof(optarg);
				break;
			case 'A':
				initialAcceptableMismatch = atof(optarg);
				break;
			case 'C':
				nsigclip = atof(optarg);
				break;
			case 'y':
				DetectionThreshold = atof(optarg);
				break;
			case 'z':
				MinPeakDepth = atof(optarg);
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
            case 'P':
				ordersplotfilename = optarg;
				plot = 1;
				break;
            case 'Q':
				specplotfilename = optarg;
				plot = 1;
				break;
			case 'F':
				ordersdatafilename = optarg;
				break; 	
			case 'G':
				atlasdatafilename = optarg;
				break;
			case 'H':
				compdatafilename = optarg;
				break;
			case 'J':
				linesdatafilename = optarg;
				break;
			case 'S':
				ordersscriptfilename = optarg;
				break;
			case 'U':
				specscriptfilename = optarg;
				break;
			case 'E':
                generate3DPlot = (atoi(optarg)?true:false);
				break;
			case 'K':
                subtractCentralWavelength = (atoi(optarg)?true:false);
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
		// we need a atlas_lines lines or spectrum...
		if (atlas_lines.empty() && atlas_spectrum.empty()) {
			throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		
		// we need EITHER thorium uncalibrated lines or spectrum...
		if (uncalibrated_lines.empty() && uncalibrated_spectrum.empty()) {
			throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		// we need a geometryfilename...
		if (geometryfilename.empty()) {
			throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
		// we need a wlcal_initialguess, the initial guess at a polynomial...
		if (wlcal_initialguess.empty()) {
			throw operaException("operaWavelengthCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}
        
		if (verbose) {
			cerr << "operaWavelengthCalibration: atlas_lines = " << atlas_lines << endl;
			cerr << "operaWavelengthCalibration: atlas_spectrum = " << atlas_spectrum << endl;            
			cerr << "operaWavelengthCalibration: uncalibrated_lines = " << uncalibrated_lines << endl;
			cerr << "operaWavelengthCalibration: uncalibrated_spectrum = " << uncalibrated_spectrum << endl;            
			cerr << "operaWavelengthCalibration: geometryfilename = " << geometryfilename << endl;            
			cerr << "operaWavelengthCalibration: wlcal_initialguess = " << wlcal_initialguess << endl; 
			cerr << "operaWavelengthCalibration: outputWave = " << outputWave << endl;            
			cerr << "operaWavelengthCalibration: parseSolution = " << parseSolution << endl;
			cerr << "operaWavelengthCalibration: normalizeUncalibratedSpectrum = " << normalizeUncalibratedSpectrum << endl;
			cerr << "operaWavelengthCalibration: normalizationBinSize = " << normalizationBinSize << endl;
			cerr << "operaWavelengthCalibration: LocalMaxFilterWidth = " << LocalMaxFilterWidth << endl;
			cerr << "operaWavelengthCalibration: ParRangeSizeInPerCent = " << ParRangeSizeInPerCent << endl;
			cerr << "operaWavelengthCalibration: NpointsPerPar = " << NpointsPerPar << endl;
			cerr << "operaWavelengthCalibration: maxNIter = " << maxNIter << endl;
			cerr << "operaWavelengthCalibration: minNumberOfLines = " << minNumberOfLines << endl;
			cerr << "operaWavelengthCalibration: maxorderofpolynomial = " << maxorderofpolynomial << endl;
			cerr << "operaWavelengthCalibration: dampingFactor = " << dampingFactor << endl;
			cerr << "operaWavelengthCalibration: initialAcceptableMismatch = " << initialAcceptableMismatch << endl;
			cerr << "operaWavelengthCalibration: nsigclip = " << nsigclip << endl;
			cerr << "operaWavelengthCalibration: DetectionThreshold = " << DetectionThreshold << endl;
			cerr << "operaWavelengthCalibration: MinPeakDepth = " << MinPeakDepth << endl;
            if(ordernumber != NOTPROVIDED) {
                cerr << "operaWavelengthCalibration: ordernumber = " << ordernumber << endl;            
            }
            if(plot) {                
                cerr << "operaWavelengthCalibration: ordersplotfilename = " << ordersplotfilename << endl;
                cerr << "operaWavelengthCalibration: specplotfilename = " << specplotfilename << endl;
                cerr << "operaWavelengthCalibration: ordersscriptfilename = " << ordersscriptfilename << endl;
                cerr << "operaWavelengthCalibration: specscriptfilename = " << specscriptfilename << endl;
                cerr << "operaWavelengthCalibration: ordersdatafilename = " << ordersdatafilename << endl;
                cerr << "operaWavelengthCalibration: atlasdatafilename = " << atlasdatafilename << endl;
                cerr << "operaWavelengthCalibration: compdatafilename = " << compdatafilename << endl;
                cerr << "operaWavelengthCalibration: linesdatafilename = " << linesdatafilename << endl;
                cerr << "operaWavelengthCalibration: generate3DPlot = " << generate3DPlot << endl;
                cerr << "operaWavelengthCalibration: subtractCentralWavelength = " << subtractCentralWavelength << endl;
            }            
		}
		
		ofstream *fatlasdata = NULL;
		ofstream *fcompdata = NULL;
		ofstream *flinesdata = NULL;
		ofstream *fordersdata = NULL;
        
        if (!atlasdatafilename.empty()) {
            fatlasdata = new ofstream();
            fatlasdata->open(atlasdatafilename.c_str());
        }
        if (!compdatafilename.empty()) {
            fcompdata = new ofstream();
            fcompdata->open(compdatafilename.c_str());
        }
        if (!linesdatafilename.empty()) {
            flinesdata = new ofstream();
            flinesdata->open(linesdatafilename.c_str());
        }
        if (!ordersdatafilename.empty()) {
            fordersdata = new ofstream();
            fordersdata->open(ordersdatafilename.c_str());
        }
        
		operaSpectralOrderVector spectralOrders(geometryfilename);	// get the geometry
		
		spectralOrders.ReadSpectralOrders(wlcal_initialguess); // read wavelength calibration reference first guess
        
        
		if (!uncalibrated_lines.empty()) {     
            spectralOrders.ReadSpectralOrders(uncalibrated_lines); // This merges in the uncalibrated lines information
        }
		
        if (!uncalibrated_spectrum.empty()) { 
            spectralOrders.ReadSpectralOrders(uncalibrated_spectrum); // This merges in the uncalibrated spectrum information
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
		
		if (verbose) {
			cerr << "operaWavelengthCalibration: minorder = " << minorder << " maxorder = " << maxorder << endl;            
		}        
        
		/*
		 * Read wavelength calibration reference first guess covered by the order: wl0,wlf
		 */        
		/*
		 * Read ThAr atlas spectrum
		 *		lambda vs. intensity, intensityVariance
		 */
		if (!atlas_spectrum.empty()) {      
			if (verbose) {
				cerr << "operaWavelengthCalibration: reading atlas spectrum " << atlas_spectrum << endl;            
			}        
            nPointsInAtlasSpectrum = readAtlasSpectrum(atlas_spectrum, atlasWavelength, atlasIntensity, atlasVariance);
        }
		/*
		 * Read ThAr atlas lines
		 *		lambda vs. intensity
		 */        
		if (!atlas_lines.empty()) {         
			if (verbose) {
				cerr << "operaWavelengthCalibration: reading atlas lines " << atlas_lines << endl;            
			}        
            thatlaslines = readThoriumArgonAtlas(atlas_lines, thAtlasWavelength, thAtlasIntensity);        
        }
        
        /*
         * vectors for uncalibrated spectral lines
         */
        unsigned rawlinesinorder=0;
        double rawlinecenter[MAXREFWAVELENGTHSPERORDER];
        double rawlinecenterError[MAXREFWAVELENGTHSPERORDER];        
        double rawlineflux[MAXREFWAVELENGTHSPERORDER];
        double rawlinesigma[MAXREFWAVELENGTHSPERORDER];
        float frawlinesigma[MAXREFWAVELENGTHSPERORDER];
        
        /*
         * vectors for atlas spectral lines
         */        
        unsigned atlaslinesinorder=0;
        double atlasLineswl[MAXREFWAVELENGTHSPERORDER]; 
        double atlasLineswlError[MAXREFWAVELENGTHSPERORDER];         
        double atlasLinesflux[MAXREFWAVELENGTHSPERORDER];        
		
        /*
         * cross-correlation
         */          
        double atlasXcorr[MAXPOINTSINSIMULATEDSPECTRUM];        
        double convolvedAtlas[MAXPOINTSINSIMULATEDSPECTRUM];        
		
        /*
         * Initialize linewidth with input uncalibrated_linewidth plus-minus 20% error. 
         */           
        double rawlinewidth;
        double rawlinewidth_err;
        
        operaSpectralOrder *spectralOrder[MAXORDERS];
        for (unsigned order=minorder; order<=maxorder; order++) {
            spectralOrder[order] = spectralOrders.GetSpectralOrder(order);
        }
        
        int ordershift = 0;
        
        for (unsigned order=minorder; order<=maxorder; order++) {
            unsigned bestnpar = maxorderofpolynomial; 
            
            rawlinewidth = uncalibrated_linewidth;
            rawlinewidth_err = uncalibrated_linewidth*0.2;
            
            spectralOrder[order] = spectralOrders.GetSpectralOrder(order);
             spectralOrder[order+ordershift] = spectralOrders.GetSpectralOrder(order+ordershift);
            
			if (spectralOrder[order]->gethasSpectralElements() && spectralOrder[order]->gethasGeometry() && spectralOrder[order+ordershift]->gethasWavelength()) {
				operaGeometry *geometry = spectralOrder[order]->getGeometry();
				operaWavelength *wavelength = spectralOrder[order+ordershift]->getWavelength();
                
                operaSpectralElements *spectralElements = NULL;
                if(spectralOrder[order]->gethasSpectralElements()) {
					
                    if(normalizeUncalibratedSpectrum) {
                        spectralOrder[order]->applyNormalizationForEmissionSpectrum(normalizationBinSize,0,FALSE,NULL,NULL,TRUE,1);
                    }
                    spectralElements = spectralOrder[order]->getSpectralElements();
                    
                    for(unsigned elemIndex=0;elemIndex<spectralElements->getnSpectralElements();elemIndex++) {
                        double elmewavelength = wavelength->evaluateWavelength(spectralElements->getdistd(elemIndex));
                        spectralElements->setwavelength(elmewavelength,elemIndex);
                        if(debug) // at this point one can check the initial wavelength calibration and the input comparison flux
                            cout << spectralElements->getwavelength(elemIndex) << " " << spectralElements->getFlux(elemIndex) << endl;
                    }
                    spectralElements->setHasWavelength(TRUE);
                }

                double dmin = 0.0;
                double dmax = (double)geometry->CalculateAndSetOrderLength();
                wavelength->setDmin(dmin);                
                wavelength->setDmax(dmax);
                
                if (verbose) {
					printf("operaWavelengthCalibration: Order %d: [geom] ymin = %.2f ymax = %.2f dmin = %.2f dmax = %.2f \n", spectralOrder[order]->getorder(), geometry->getYmin(), geometry->getYmax(), dmin, dmax);
				}   
				
                /*
                 * Below it calculates the initial and final wavelength based on the geometry calibration
                 * Note: this will be used to select the atlas range.  
                 */                     
                double wl0 = wavelength->getinitialWavelength();
                double wlf = wavelength->getfinalWavelength();
                wl0 -= (wavelength->getcentralWavelength())*(ParRangeSizeInPerCent/100);
                wlf += (wavelength->getcentralWavelength())*(ParRangeSizeInPerCent/100);                    
                
                double wl_central = wavelength->getcentralWavelength();
                
				if (verbose) {
					printf("operaWavelengthCalibration: Order %d: [wave] wavelength selected range: wl0 = %.2f wlc = %.2f wlf = %.2f\n",  spectralOrder[order]->getorder(), wl0, wavelength->getcentralWavelength(),wlf);
				}
				
                /*
                 * Read raw distance and flux from uncalibrated_lines file
                 */
                if (!uncalibrated_lines.empty()) {             
					if (verbose) {
						cerr << "operaWavelengthCalibration: reading uncalibrated lines " << uncalibrated_lines << endl;            
					}        
                    rawlinesinorder = getRawLines(spectralOrder[order], rawlinecenter, rawlinecenterError, rawlineflux, rawlinesigma);
                    for (unsigned i=0; i<rawlinesinorder; i++)
                        frawlinesigma[i] = (float)rawlinesigma[i];	
                    // DT Sept 18 -- the Quick version can be used since frawlinesigma is never used again...
					rawlinewidth = (double)operaArrayMedianQuick(rawlinesinorder,frawlinesigma);
                    rawlinewidth_err = (double)operaArrayMedianSigmaQuick(rawlinesinorder,frawlinesigma,(float)rawlinewidth);                
                }
                /*
                 * Below it reads raw distance and flux of lines detected from comparison input spectrum
                 * Note: the step below overrides the lines obtained from uncalibrated_lines above if both are provided. 
                 */            
                if (!uncalibrated_spectrum.empty() && spectralOrder[order]->gethasSpectralElements()) {
					if (verbose) {
						cerr << "operaWavelengthCalibration: reading uncalibrated spectrum " << uncalibrated_spectrum << endl;            
					}        
                    operaSpectralElements *compSpectrum = spectralOrder[order]->getSpectralElements();
                    
                    operaSpectralLines compLines(compSpectrum,rawlinewidth,distance_disp);
					
                    operaFluxVector *compfluxvector = compSpectrum->getFluxVector();
                    if(debug) {
                        for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                            cout << compSpectrum->getdistd(i) << " " << compfluxvector->getflux(i) << " " << compSpectrum->getXCorrelation(i) << endl;
                        }
                    }
                    if(!compSpectrum->getHasXCorrelation()){
						double compSpectrumdistd[MAXPOINTSINSIMULATEDSPECTRUM];         
						double compSpectrumflux[MAXPOINTSINSIMULATEDSPECTRUM]; 
                        double compXcorr[MAXPOINTSINSIMULATEDSPECTRUM];
                        
						for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
							compSpectrumdistd[i] = compSpectrum->getdistd(i);
							compSpectrumflux[i] = compfluxvector->getflux(i);
						}
                        
                        calculateXCorrWithGaussian(compSpectrum->getnSpectralElements(), compSpectrumdistd,compSpectrumflux,compXcorr,rawlinewidth);
						for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
							compSpectrum->setXCorrelation(compXcorr[i], i);
                        }
						compSpectrum->setHasXCorrelation(true); 
                    }
                    double CompLocalMaxFilterWidth = LocalMaxFilterWidth*rawlinewidth;
					
                    double meanVariance = 0;
                    unsigned nvarpoints = 0;
                    for (unsigned i=0; i<compSpectrum->getnSpectralElements(); i++) {
                        if(!isnan(compfluxvector->getvariance(i))) {
                            meanVariance += compfluxvector->getvariance(i);
                            nvarpoints++;
                        }
                    } 
                    double CompMinPeakDepth = MinPeakDepth*sqrt(meanVariance/(double)nvarpoints);
					
                    compLines.detectSpectralFeatures(DetectionThreshold,CompLocalMaxFilterWidth,CompMinPeakDepth);
					
                    unsigned line = 0;
					
                    for(unsigned feature=0;feature<compLines.getNFeatures();feature++) {
                        operaSpectralFeature *currentFeature = compLines.getSpectralFeature(feature);
                        double *center = currentFeature->getGaussianFit()->getCenterVector();
                        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
                        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
                        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector(); 
                        for(unsigned l=0; l<currentFeature->getnLines(); l++) { 
                            rawlinecenter[line] = center[l];
                            rawlinecenterError[line] = centerError[l];
                            rawlineflux[line] = amplitude[l];
                            rawlinesigma[line] = sigma[l];
                            frawlinesigma[line] = (float)rawlinesigma[line];                            
                            if(debug)
                                cout << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << " " << " " << currentFeature->getnLines() << " " << currentFeature->getGaussianFit()->getGaussianChisqr() << endl;
                            line++;
                        }
                    }
                    
					rawlinesinorder = line;
                    if (rawlinesinorder > 0) {
						rawlinewidth = (double)operaArrayMedian(rawlinesinorder,frawlinesigma);
						rawlinewidth_err = (double)operaArrayMedianSigma(rawlinesinorder,frawlinesigma,(float)rawlinewidth);   
						//                    rawlinewidth = (double)operaArrayMean(rawlinesinorder,frawlinesigma);                    
						//                    rawlinewidth_err = (double)operaArraySig(rawlinesinorder,frawlinesigma);
						
                        
                        line = 0;
                        for(unsigned feature=0;feature<compLines.getNFeatures();feature++) {
                            operaSpectralFeature *currentFeature = compLines.getSpectralFeature(feature);
                            double *center = currentFeature->getGaussianFit()->getCenterVector();
                            double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();
                            double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
                            double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector();
                            for(unsigned l=0; l<currentFeature->getnLines(); l++) {
                                if(sigma[l] > rawlinewidth - (double)nsigclip*rawlinewidth_err && sigma[l] < rawlinewidth + (double)nsigclip*rawlinewidth_err) {
                                    rawlinecenter[line] = center[l];
                                    rawlinecenterError[line] = centerError[l];
                                    rawlineflux[line] = amplitude[l];
                                    rawlinesigma[line] = sigma[l];
                                    frawlinesigma[line] = (float)rawlinesigma[line];
                                    if(debug)
                                        cout << center[l] <<  " " << centerError[l] << " " << amplitude[l] << " " << sigma[l] << " " << " " << currentFeature->getnLines() << " " << currentFeature->getGaussianFit()->getGaussianChisqr() << endl;
                                    line++;
                                }
                            }
                        }
                        
                        rawlinesinorder = line;
                        if (rawlinesinorder > 0) {
                            rawlinewidth = (double)operaArrayMedian(rawlinesinorder,frawlinesigma);
                            rawlinewidth_err = (double)operaArrayMedianSigma(rawlinesinorder,frawlinesigma,(float)rawlinewidth);
                            compSpectrum->setHasWavelength(true);
                        }
					}
                } 
                if (rawlinesinorder == 0) {
                    if (verbose) {
						printf("operaWavelengthCalibration: Warning: Order %d: [Comparison] No lines detected from input comparison. Skipping calibration.\n", order); 
                    }
					continue;
                } else {
                    /*
                     * Below it normalizes the raw lines
                     */                                      
                    if (verbose) {
                        printf("operaWavelengthCalibration: Order %d: [Comparison] %d lines in comparison between wl0 = %.2f and wlf = %.2f.\n", order, rawlinesinorder,wavelength->evaluateWavelength(rawlinecenter[0]),wavelength->evaluateWavelength(rawlinecenter[rawlinesinorder-1])); 
                    }
                    wavelength->createComparisonDataVectors(rawlinesinorder,rawlinecenter,rawlinecenterError,rawlineflux);
                }
				
                /*
                 * Read wavelength and flux from atlas_lines file
                 */                                         
				
                if (!atlas_lines.empty()) {   
					if (verbose) {
						cerr << "operaWavelengthCalibration: reading atlas lines " << atlas_lines << endl;            
					}        
                    double *thwl = NULL, *thintensity = NULL;
                    atlaslinesinorder = getThoriumArgonAtlasRange(wl0, wlf, &thwl, &thintensity);  
                    for(unsigned l=0; l<atlaslinesinorder; l++) {
                        atlasLineswl[l] = *thwl++;
                        atlasLinesflux[l] = *thintensity++;
                        atlasLineswlError[l] = wavelength->convertPixelToWavelength(rawlinewidth);
                    }
                } 
				
                /*
                 * Below it reads wavelength and flux from the atlas spectrum
                 */    
                if (!atlas_spectrum.empty()) {    
					if (verbose) {
						cerr << "operaWavelengthCalibration: reading atlas spectrum " << atlas_spectrum << endl;            
					}        
                    double *thwl = NULL, *thintensity = NULL, *thvar = NULL;
                    unsigned npatlasspecinorder = getAtlasSpectrumRange(wl0, wlf, &thwl, &thintensity, &thvar);                


                   if(fatlasdata != NULL && generate3DPlot==TRUE) {
                       double maxatlasflux = -BIG;
                       for (unsigned i=0; i<npatlasspecinorder; i++) {
                           if(thintensity[i] > maxatlasflux)
                               maxatlasflux = thintensity[i];
                       }
                        /*** NOTE ***/
                        // Below it produces data for a 3D plot of spectrum.
                       
                        for(unsigned slitview=0;slitview<2;slitview++){
                            for (unsigned i=0; i<npatlasspecinorder; i++) {
                                *fatlasdata << order << " " << thwl[i] << " " << thintensity[i]/maxatlasflux << " " << slitview << " " << wl_central << endl;
                            }
                            *fatlasdata << endl;
                        }
                        *fatlasdata << endl;
                   } else if (fatlasdata != NULL && generate3DPlot==FALSE) {
                       /*** NOTE ***/
                       // Below it produces data for a 2D plot of spectrum.
                       double maxatlasflux = -BIG;
                       for (unsigned i=0; i<npatlasspecinorder; i++) {
                           if(thintensity[i] > maxatlasflux)
                               maxatlasflux = thintensity[i];
                       }
                       for (unsigned i=0; i<npatlasspecinorder; i++) {
                           *fatlasdata << order << " " << thwl[i] << " " << thintensity[i]/maxatlasflux << " " << wl_central << endl;
                       }
                       *fatlasdata << endl;
                   }
                    
                    /*
                     * Below it calculates the cross-correlation between the atlas spectrum and a gaussian function. 
                     */                                   
					if (verbose) {
						cerr << "operaWavelengthCalibration: calculating cross correlation with Gaussian.." << endl;            
					}        
                    calculateXCorrWithGaussian(npatlasspecinorder,thwl,thintensity,atlasXcorr,wavelength->convertPixelToWavelength(rawlinewidth));                
                    
                    /*
                     * Below it degrades the resolution of the atlas to the resolution of raw lines.
                     * The degradation is done by convolving the spectrum with a gaussian. 
                     */                      
					if (verbose) {
						cerr << "operaWavelengthCalibration: Convolving specctrum with Gaussian.." << endl;            
					}        
                    convolveSpectrumWithGaussian(npatlasspecinorder,thwl,thintensity,convolvedAtlas,wavelength->convertPixelToWavelength(rawlinewidth));                
					
                    /*
                     * Below it reads the atlas spectrum into an operaSpectralElements class
                     */                     
                    operaSpectralElements atlasSpectrum(npatlasspecinorder);
					for (unsigned i=0; i<npatlasspecinorder; i++) {
						atlasSpectrum.setXCorrelation(atlasXcorr[i], i);
					}
                    atlasSpectrum.setHasXCorrelation(true); 
                    //atlasSpectrum.setwavelengthVector(thwl);
					for (unsigned i=0; i<npatlasspecinorder; i++) {
						atlasSpectrum.setwavelength(thwl[i], i);
					}
                    atlasSpectrum.setHasWavelength(true);
                    operaFluxVector atlasfluxvector(convolvedAtlas,thvar,npatlasspecinorder);
                    atlasSpectrum.setFluxVector(&atlasfluxvector);	// copies, does not set local stack address...
                    atlasSpectrum.setHasRawSpectrum(true);
					
                    /*
                     * Below it creates an operaSpectralLines class for the atlas lines
                     */                     
                    operaSpectralLines atlasLines(&atlasSpectrum, wavelength->convertPixelToWavelength(rawlinewidth), wavelength_disp);
                    if(debug) {
                        for (unsigned i=0; i<npatlasspecinorder; i++) {
                            cout << atlasSpectrum.getwavelength(i)  << " " << atlasSpectrum.getFlux(i) << " " << atlasSpectrum.getFluxVariance(i) << " " << atlasSpectrum.getXCorrelation(i) << endl;
                        }
                    }
				
                    /*
                     * Below it set detection thresholds and run the algorithm to detect spectral lines in the atlas
                     */  
					
                    double AtlasLocalMaxFilterWidth = LocalMaxFilterWidth*wavelength->convertPixelToWavelength(uncalibrated_linewidth);
                    double AtlasMinPeakDepth = MinPeakDepth*sqrt(operaArrayMean_d(npatlasspecinorder,thvar));
                    atlasLines.detectSpectralFeatures(DetectionThreshold,AtlasLocalMaxFilterWidth,AtlasMinPeakDepth);
                    
                    /*
                     * Below it reads the atlas lines information
                     */                      
                    atlaslinesinorder = atlasLines.getnLines();
                    
                    unsigned line = 0;
                    
                    for(unsigned feature=0;feature<atlasLines.getNFeatures();feature++) {
                        operaSpectralFeature *currentFeature = atlasLines.getSpectralFeature(feature);
                        double *center = currentFeature->getGaussianFit()->getCenterVector();
                        double *centerError = currentFeature->getGaussianFit()->getCenterErrorVector();                        
                        double *sigma = currentFeature->getGaussianFit()->getSigmaVector();
                        double *amplitude = currentFeature->getGaussianFit()->getAmplitudeVector(); 
                        for(unsigned l=0; l<currentFeature->getnLines(); l++) {
                            if(center[l] > wl0 && center[l] < wlf) {
                                atlasLineswl[line] = center[l];
                                atlasLineswlError[line] = centerError[l];
                                atlasLinesflux[line] = amplitude[l];
                                if(debug)
                                    cout << center[l] <<  " " << amplitude[l] << " " << sigma[l] << endl;
                                line++;
                            }
                        }
                    }
                    atlaslinesinorder = line;
                }
				
                if (atlaslinesinorder == 0) {
                    printf("operaWavelengthCalibration: Warning:  Order %d: [Atlas] No lines detected from input atlas. Skipping calibration.\n", order);
                    continue;
                } else {                             
                    if (verbose) {
                        printf("operaWavelengthCalibration: Order %d: [Atlas] %d lines detected in input atlas between wl0 = %.2f and wlf = %.2f .\n", order, atlaslinesinorder, wl0, wlf);
                    }
                    wavelength->createAtlasDataVectors(atlaslinesinorder,atlasLineswl, atlasLineswlError,atlasLinesflux);
                }
				
                /*
                 * At this point all possible lines either in comparison or in the Atlas have been read.
                 * So, it starts identification of lines and refining wavelength solution.
                 */
                wavelength->createDataVectors(MAXREFWAVELENGTHS);
                doubleValue_t ResolutionElementInPixels = {2*rawlinewidth, rawlinewidth_err};
                wavelength->calculateSpectralResolution(ResolutionElementInPixels);

                if(parseSolution) {
                    /*
                     * The loop below adjust the wavelength solution to the maximum correlation between the two sets of lines.
                     * The adjustment is done by varying each coefficient of the polynomial solution and finding the maximum
                     * x-correlation.
                     * maxpolyorder is the highest coefficient that is changed in the search
                     * NpointsPerPar is the number of trial solutions for a given coefficient
                     * ParRangeSizeInPerCent is a quantity that defined the range for the search. The range is set
                     * relative to the value of the coefficient.
                     */
                    for(unsigned iter = 0; iter < nIterforXcorr; iter++) {
                        wavelength->refineWavelengthSolutionByXCorrelation(NpointsPerPar, ParRangeSizeInPerCent, maxpolyorder);
                    }
                } else {
                    /*
                     *  If choose not to parse parameters then it scans over
                     *  the zeroth order coefficient in order to find the highest
                     *  matching rate. This allows the spectrum to be identified 
                     *  even when there is a shift of the spectra on the detector
                     */
                    wavelength->refineWavelengthSolutionByFindingMaxMatching(NpointsPerPar, ParRangeSizeInPerCent, initialAcceptableMismatch);
                }
                
                double acceptableMismatch = initialAcceptableMismatch; // in units of sigma
                
                wavelength->calculateSpectralResolution(ResolutionElementInPixels);
                wavelength->matchAtlaswithComparisonLines(acceptableMismatch);
                                
				//
				//if wavelength->getnDataPoints() == 0, the result is -1!!!
				//
				Polynomial *wavelengthPolynomial = wavelength->getWavelengthPolynomial();
				double *bestpar = (double *)wavelengthPolynomial->getVector();                
				double minchisqr = wavelengthPolynomial->getChisqr();
				bestnpar = wavelengthPolynomial->getOrderOfPolynomial();
                                
                if (wavelength->getnDataPoints() > 0) {
					if(wavelength->getnDataPoints() <= bestnpar) {
						bestnpar = wavelength->getnDataPoints()-1;
					}                
					wavelength->CalculateWavelengthSolution(bestnpar,false);
					
					unsigned nochangeinChisqr = 0;
					
					if (verbose) {                  
						double ComparisonMatchPercentage =  wavelength->getPerCentageOfComparisonMatch();    
						double AtlasMatchPercentage =  wavelength->getPerCentageOfAtlasMatch();                     
						printf("\noperaWavelengthCalibration: Order %d: Initial Solution:\n", order);
						printf("operaWavelengthCalibration: Order %d: ", order);
                        wavelength->getWavelengthPolynomial()->printEquation(&cout);           
						cout <<  " chisqr=" << minchisqr << endl; 
						printf("operaWavelengthCalibration: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order,wavelength->getnDataPoints(),wavelength->getinitialWavelength(),wavelength->getfinalWavelength());
						printf("operaWavelengthCalibration: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order,AtlasMatchPercentage);
						printf("operaWavelengthCalibration: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order,ComparisonMatchPercentage);                        
					}
					/*
					 * Below it starts the iterations to improve the polynomial that gives the pixel-to-wavelength solution
					 */  
					for(unsigned iter=0; iter < maxNIter; iter++) {	
						
						if (debug) {
							printf("operaWavelengthCalibration: Refinement iteration %d\n", iter);
						}
                        
						/*** NOTE ***/
						// The acceptable mismatch can start considerably big and then shrink down as calibration
                        // gets better. However, it will only shrink either to a minimum value or minimum number of lines. 
						if(wavelength->getnDataPoints() > minNumberOfLines) {                                          
							acceptableMismatch *= dampingFactor; 
						}
						wavelength->calculateSpectralResolution(ResolutionElementInPixels);
                        
						wavelength->matchAtlaswithComparisonLines(acceptableMismatch);
                        
						wavelength->filterDataPointsBySigmaClip((double)nsigclip);
						
						wavelength->RefineWavelengthSolution(bestnpar,false);
						
						if(wavelength->getnDataPoints() <= bestnpar) {
                            bestnpar = wavelength->getnDataPoints()-1;
						}
						
						// Note that witherrors is set to false here, so polynomial errors are all zero
						wavelength->CalculateWavelengthSolution(bestnpar,false);
						
						wavelengthPolynomial = wavelength->getWavelengthPolynomial();
						bestpar = (double *)wavelengthPolynomial->getVector();
						bestnpar = wavelengthPolynomial->getOrderOfPolynomial();
						
						if(wavelengthPolynomial->getChisqr() < minchisqr) {
							minchisqr = wavelengthPolynomial->getChisqr();
							nochangeinChisqr = 0;
						} else if(wavelengthPolynomial->getChisqr() == minchisqr) {
							nochangeinChisqr++;
						}
                    
                        if(debug) {
                            double ComparisonMatchPercentage =  wavelength->getPerCentageOfComparisonMatch();
                            double AtlasMatchPercentage =  wavelength->getPerCentageOfAtlasMatch();
                            printf("\noperaWavelengthCalibration: Order %d: Initial Solution:\n", order);
                            printf("operaWavelengthCalibration: Order %d: ", order);
                            wavelength->getWavelengthPolynomial()->printEquation(&cout);
                            cout <<  " chisqr=" << minchisqr << endl;
                            printf("operaWavelengthCalibration: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order,wavelength->getnDataPoints(),wavelength->getinitialWavelength(),wavelength->getfinalWavelength());
                            printf("operaWavelengthCalibration: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order,AtlasMatchPercentage);
                            printf("operaWavelengthCalibration: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order,ComparisonMatchPercentage);
                        }
                        
						if(nochangeinChisqr > 3) {
							break;
						}
					} // for(unsigned iter=0; iter < nIter; iter++)
				} else {
					printf("operaWavelengthCalibration: Order %d: ZERO points to calculate wavelength solution, skipping order,,.\n", order);
				}
               
                if(wavelength->getnDataPoints() > 0) {
					if(debug) {
						for(unsigned l=0 ; l<wavelength->getnDataPoints(); l++) {
							cout << order << " " << wavelength->getDistance(l) << " " << wavelength->getWavelength(l) << " " << wavelength->evaluateWavelength(wavelength->getDistance(l)) << " " << wavelength->getWavelength(l) - wavelength->evaluateWavelength(wavelength->getDistance(l)) << " " << wavelength->getWavelengthError(l) << endl;
						}
					}
                    // wavelength->getWavelengthPolynomial()->printEquation(&cout);
					if (verbose) {
						double ComparisonMatchPercentage =  wavelength->getPerCentageOfComparisonMatch();    
						double AtlasMatchPercentage =  wavelength->getPerCentageOfAtlasMatch();                     
						wavelengthPolynomial = wavelength->getWavelengthPolynomial();
						bestpar = (double *)wavelengthPolynomial->getVector();
						bestnpar = wavelengthPolynomial->getOrderOfPolynomial();
                        printf("\n");
						printf("operaWavelengthCalibration: Order %d: Wavelength solution after done shrinking:\n", order);
						printf("operaWavelengthCalibration: Order %d: ", order);
                        wavelength->getWavelengthPolynomial()->printEquation(&cout);
						cout <<  " chisqr=" << wavelength->getWavelengthPolynomial()->getChisqr() << endl;
						printf("operaWavelengthCalibration: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order,wavelength->getnDataPoints(),wavelength->getinitialWavelength(),wavelength->getfinalWavelength());
						printf("operaWavelengthCalibration: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order,AtlasMatchPercentage);
						printf("operaWavelengthCalibration: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order,ComparisonMatchPercentage);                        
					}
                    /*
                    cout << order << " ";
                    for(unsigned k=0;k<bestnpar;k++) {
                        cout << bestpar[k] << " ";
                    }
                    cout << endl;
                     */
					/*
					 ***** at this point there should be a good wavelength solution *****
					 */ 
					/*
					 * Below it calculates the radial velocity precision
					 */                 
					wavelength->calculateRadialVelocityPrecision();
					
					/*
					 * Calculate the spectral resolution
					 */                  
					if (verbose){
						printf("operaWavelengthCalibration: -----------------------------------------------------------------\n");
						printf("operaWavelengthCalibration: Order %d: Radial velocity precision = %.2f m/s. Using %d spectral lines.\n", order, wavelength->getRadialVelocityPrecision(),wavelength->getnDataPoints());
						printf("operaWavelengthCalibration: Order %d: Wavelength RMS precision = %.10f nm  Median Precision = %.10f nm.\n", order,wavelength->calculateWavelengthRMSPrecision(),wavelength->calculateWavelengthMedianPrecision());
						printf("operaWavelengthCalibration: Order %d: [Comparison Lines] median sigma = %.2f +/- %.2f.\n", order, rawlinewidth, rawlinewidth_err);
					}
					wavelength->calculateSpectralResolution(ResolutionElementInPixels);
					
					if (debug)
						printf("%d %.2f %.2f %.2f %.2f %.2f\n", order, wavelength->getcentralWavelength(),ResolutionElementInPixels.value, ResolutionElementInPixels.error, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
					
					if (verbose) {
						printf("operaWavelengthCalibration: Order %d: Spectral Resolution = %.2f +/- %.2f.\n", order, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
						printf("operaWavelengthCalibration: -----------------------------------------------------------------\n\n");                
					}
                    if(debug)
                        printf("%d %.2f %d %.2f %.2f\n", order, wavelength->getRadialVelocityPrecision(),wavelength->getnDataPoints(),wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
					
					if (debug) {                
						cout << order << " " 
						<< wavelength->getinitialWavelength() <<  " "
						<< wavelength->getcentralWavelength() << " " 
						<< wavelength->getfinalWavelength() << " "
						<< wavelength->calculateWavelengthRMSPrecision() << " "
						<< wavelength->calculateWavelengthMedianPrecision() << " "
						<< wavelength->getRadialVelocityPrecision() << " "
						<< ResolutionElementInPixels.value << " "
						<< ResolutionElementInPixels.error << " "
						<< wavelength->getSpectralResolution().value << " "
						<< wavelength->getSpectralResolution().error << endl;
					}
				}
                
                wavelength->createComparisonDataVectors(rawlinesinorder,rawlinecenter,rawlinecenterError,rawlineflux);
                wavelength->createAtlasDataVectors(atlaslinesinorder,atlasLineswl, atlasLineswlError,atlasLinesflux);
                wavelength->calculateSpectralResolution(ResolutionElementInPixels);
                wavelength->matchAtlaswithComparisonLines(ResolutionElementInPixels.value/2);
                bestnpar = maxorderofpolynomial;
                
                if (wavelength->getnDataPoints() > 1) {
					if(wavelength->getnDataPoints() <= bestnpar) {
						bestnpar = wavelength->getnDataPoints()-1;
					}
                    wavelength->filterDataPointsBySigmaClip((double)nsigclip/2);
					wavelength->CalculateWavelengthSolution(bestnpar,false);
                    
                    double ComparisonMatchPercentage =  wavelength->getPerCentageOfComparisonMatch();
                    double AtlasMatchPercentage =  wavelength->getPerCentageOfAtlasMatch();
                    
					if (verbose) {
                        printf("operaWavelengthCalibration: *****************************************************************\n");
						printf("operaWavelengthCalibration: Order %d: Final Solution:\n", order);
						printf("operaWavelengthCalibration: Order %d: ", order);
                        wavelength->getWavelengthPolynomial()->printEquation(&cout);
						cout <<  " chisqr=" << wavelength->getWavelengthPolynomial()->getChisqr() << endl;
						printf("operaWavelengthCalibration: Order %d: %u lines matched between wl0 = %.2f  wlf = %.2f.\n", order,wavelength->getnDataPoints(),wavelength->getinitialWavelength(),wavelength->getfinalWavelength());
						printf("operaWavelengthCalibration: Order %d: [Atlas]    matched %.2f %% of detected lines.\n", order,AtlasMatchPercentage);
						printf("operaWavelengthCalibration: Order %d: [Comparison] matched %.2f %% of detected lines.\n", order,ComparisonMatchPercentage);
                        printf("operaWavelengthCalibration: *****************************************************************\n");                        
					}
                    if (fordersdata != NULL) {
						*fordersdata << order << " "
						<< wavelength->getinitialWavelength() <<  " "
						<< wavelength->getcentralWavelength() << " "
						<< wavelength->getfinalWavelength() << " "
						<< wavelength->calculateWavelengthRMSPrecision() << " "
						<< wavelength->calculateWavelengthMedianPrecision() << " "
						<< ResolutionElementInPixels.value << " "
						<< ResolutionElementInPixels.error << " "
						<< wavelength->getnDataPoints() << " "
						<< ComparisonMatchPercentage  << " "                       
						<< AtlasMatchPercentage << " "                        
                        << wavelength->getRadialVelocityPrecision() << " ";
                    }
                    
					wavelength->calculateRadialVelocityPrecision();

					if (verbose){
						printf("operaWavelengthCalibration: -----------------------------------------------------------------\n");
						printf("operaWavelengthCalibration: Order %d: Radial velocity precision = %.2f m/s. Using %d spectral lines.\n", order, wavelength->getRadialVelocityPrecision(),wavelength->getnDataPoints());
						printf("operaWavelengthCalibration: Order %d: Wavelength RMS precision = %.10f nm\n", order,wavelength->calculateWavelengthRMSPrecision());
                        printf("operaWavelengthCalibration: Order %d: Median Precision = %.10f nm.\n", order,wavelength->calculateWavelengthMedianPrecision());
					}
					wavelength->calculateSpectralResolution(ResolutionElementInPixels);
					
					if (debug)
						printf("%d %.2f %.2f %.2f %.2f %.2f\n", order, wavelength->getcentralWavelength(),ResolutionElementInPixels.value, ResolutionElementInPixels.error, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
					
					if (verbose) {
						printf("operaWavelengthCalibration: Order %d: Spectral Resolution = %.2f +/- %.2f.\n", order, wavelength->getSpectralResolution().value, wavelength->getSpectralResolution().error);
						printf("operaWavelengthCalibration: -----------------------------------------------------------------\n\n");
					}

                    if (fordersdata != NULL) {
						*fordersdata << wavelength->getRadialVelocityPrecision() << " "
						<< wavelength->getSpectralResolution().value << " "
						<< wavelength->getSpectralResolution().error/2 << endl;
					}
                }
               
                double maxflux = -BIG;
                if(flinesdata != NULL || fcompdata != NULL) {
                    for (unsigned i=0; i<spectralElements->getnSpectralElements(); i++) {
                        if (spectralElements->getFlux(i) > maxflux) {
                            maxflux = spectralElements->getFlux(i);
                        }
                    }
                }
                
                if(fcompdata != NULL && generate3DPlot==TRUE) {
                    /*** NOTE ***/
                    // Below it produces data for a 3D plot of spectrum.

                    for(unsigned slitview=0;slitview<2;slitview++){
                        unsigned lastline = 0;
                        for (unsigned i=0; i<spectralElements->getnSpectralElements(); i++) {
                            double dist = spectralElements->getdistd(i);
                            double wl = wavelength->evaluateWavelength(dist);
                            double flux = spectralElements->getFlux(i);
                            double matchlinesflux = 0;
                            
                            for(unsigned lines=lastline;lines<wavelength->getnDataPoints();lines++){
                                if(dist > wavelength->getDistance(lines) - rawlinewidth && dist < wavelength->getDistance(lines) + rawlinewidth) {
                                    matchlinesflux = 1;
                                    lastline = lines;
                                    break;
                                } else if (dist > wavelength->getDistance(lines) + rawlinewidth) {
                                    lastline++;
                                    break;
                                }
                            }
                            
                            *fcompdata << order << " " << dist << " " << wl << " " << flux/maxflux << " " << matchlinesflux << " " << slitview  << " " << wl_central << endl;
                        }
                        *fcompdata << endl;
                    }
                    *fcompdata << endl;
                } else if (fcompdata != NULL && generate3DPlot==FALSE){
                    /*** NOTE ***/
                    // Below it produces data for a 2D plot of spectrum.
                    
                    for (unsigned i=0; i<spectralElements->getnSpectralElements(); i++) {
                        double dist = spectralElements->getdistd(i);
                        double wl = wavelength->evaluateWavelength(dist);
                        double flux = spectralElements->getFlux(i);
                        *fcompdata << order << " " << dist << " " << wl << " " << flux/maxflux << " " << wl_central << endl;
                    }
                    *fcompdata << endl;
                }
                
                if(flinesdata != NULL) {
                    for(unsigned lines=0;lines<wavelength->getnDataPoints();lines++){
                        *flinesdata << order << " " << wavelength->getDistance(lines) << " " << wavelength->getWavelength(lines) << " " << ((wavelength->getcomparisonLinesflux(wavelength->getMatchComparisonIndex(lines))/maxflux)/2.0) << " " << wavelength->convertPixelToWavelength(ResolutionElementInPixels.value) << " " << wl_central << endl;
                    }
                    *flinesdata << endl;
                }
                
            } else if (!spectralOrder[order]->gethasWavelength()) {
                if (verbose) {
                    printf("operaWavelengthCalibration: Order %d: has no associated wavelength reference calibration data. Wavelength calibration not possible.\n",  spectralOrder[order]->getorder());
                }
            } else if (!spectralOrder[order]->gethasGeometry()) {
                if (verbose) {
                    printf("operaWavelengthCalibration: Order %d: has no associated geometry data. Wavelength calibration not possible.\n",  spectralOrder[order]->getorder());
                }
                spectralOrder[order]->sethasWavelength(false);
            } else if (!spectralOrder[order]->gethasSpectralElements()) {
                if (verbose) {
                    printf("operaWavelengthCalibration: Order %d: has no associated spectral elements data. Wavelength calibration not possible.\n",  spectralOrder[order]->getorder());
                }
                spectralOrder[order]->sethasWavelength(false);
            } else {
                if (verbose) {
                    printf("operaWavelengthCalibration: Order %d: has neither geometry nor wavelength reference calibration data. Wavelength calibration not possible.\n",  spectralOrder[order]->getorder());
                }
                spectralOrder[order]->sethasWavelength(false);
            }            
		} // for (unsigned order=minorder; order<=maxorder; order++)
       
        /*
         * Wavelength Orders Info Plot: plot spectral resolution, rms in nm, radial velocity precision, data and polynomial solution
         */
        if (fordersdata != NULL) {
            fordersdata->close();
            if (!ordersscriptfilename.empty()) {
               GenerateWavelengthOrdersPlot(ordersscriptfilename, ordersplotfilename, ordersdatafilename, interactive);
            }
        }  

        /*
         * Wavelength Spectrum Plot: plot atlas and comparison spectra and final set of matched lines.
         */
        if(generate3DPlot) {
            if (fatlasdata != NULL && fcompdata != NULL) {
                fatlasdata->close();
                fcompdata->close();
                if (!specscriptfilename.empty()) {
                    GenerateWavelength3DSpecPlot(specscriptfilename, specplotfilename, atlasdatafilename, compdatafilename, subtractCentralWavelength, interactive);
                }
            }
        } else {
            if (fatlasdata != NULL && fcompdata != NULL && flinesdata != NULL) {
                fatlasdata->close();
                fcompdata->close();
                flinesdata->close();
                
                if (!specscriptfilename.empty()) {
                    GenerateWavelengthSpecPlot(specscriptfilename, specplotfilename, atlasdatafilename, compdatafilename, linesdatafilename, subtractCentralWavelength, interactive);
                }
            }
        }
        spectralOrders.WriteSpectralOrders(outputWave, Wave);  
        
        //delete[] convolvedAtlas;
        exit(EXIT_SUCCESS);
	}
	catch (operaException e) {
		cerr << "operaWavelengthCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaWavelengthCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
}


/*
 * Generate multiple plot containing statistical info about wavelength calibration
 */

void GenerateWavelengthOrdersPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName,bool display)
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
    *fgnu << "unset key\n" << endl;
    *fgnu << "NX=2; NY=2" << endl;
    *fgnu << "DX=0.1; DY=0.1; SX=0.42; SY=0.42" << endl;
    *fgnu << "set bmargin DX; set tmargin DX; set lmargin DY; set rmargin DY" << endl;
    *fgnu << "set size SX*NX+DX*2,SY*NY+DY*2" << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "set multiplot" << endl;

        *fgnu << "set size 0.9*SX,0.9*SY" << endl;
        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2tics; set ylabel \"{/Symbol l} precision (nm)\"" << endl;
        *fgnu << "set origin 0.75*DX,DY+SY" << endl;
        *fgnu << "set key" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:5 t \"rms of residuals\" w linespoint lw 1.5, \"\" u 1:6  t \"median of residuals\" w linespoint lw 1.5" << endl;
        
        *fgnu << "unset key" << endl;
        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset ytics; set y2tics mirror" << endl;
        *fgnu << "unset ylabel; set y2label \"spectral resolution {/Symbol l}/{/Symbol Dl}\"" << endl;
        *fgnu << "set origin DX+SX,DY+SY" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:14:15 w yerr pt 7 lw 1.5 lt 1,\"\" u 1:14 w linespoint lt 1" << endl;
        
        *fgnu << "set key" << endl;
        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2label; set ylabel \"Radial velocity precision (m/s)\"" << endl;
        *fgnu << "set origin 0.75*DX,DY" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:13 t \"full set of lines\" w linespoint pt 7 lt 3 lw 2, \"\" u 1:12 t \"clean sample\" w linespoint pt 7 lt 4 lw 2" << endl;
        *fgnu << "unset key" << endl;
        
        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset y2tics; unset y2label" << endl;
        *fgnu << "set origin DX+SX,DY" << endl;
        *fgnu << "set ytics nomirror" << endl;
        *fgnu << "set ylabel \"% matching lines\" offset +1.5,0l" << endl;
        *fgnu << "set key bottom" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:10 t \"% of comparison lines\" w linespoint pt 7 lt 4 lw 1, \"\" u 1:11 t \"% of atlas lines\" w linespoint pt 6 lt 4 lw 1" << endl;
        
        *fgnu << "unset key" << endl;
        *fgnu << "unset ytics; set y2tics" << endl;
        *fgnu << "unset ylabel; set y2label \"Number of matching lines\"" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:9 w linespoint lt 2 lw 3" << endl;

        *fgnu << "unset multiplot" << endl;
        
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
        *fgnu << "set multiplot" << endl;
        
        *fgnu << "set size 0.9*SX,0.9*SY" << endl;
        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2tics; set ylabel \"{/Symbol l} precision (nm)\"" << endl;
        *fgnu << "set origin 0.75*DX,DY+SY" << endl;
        *fgnu << "set key" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:5 t \"rms of residuals\" w linespoint lw 1.5, \"\" u 1:6  t \"median of residuals\" w linespoint lw 1.5" << endl;
        
        *fgnu << "unset key" << endl;
        *fgnu << "unset xlabel" << endl;
        *fgnu << "unset ytics; set y2tics mirror" << endl;
        *fgnu << "unset ylabel; set y2label \"spectral resolution {/Symbol l}/{/Symbol Dl}\"" << endl;
        *fgnu << "set origin DX+SX,DY+SY" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:14:15 w yerr pt 7 lw 1.5 lt 1,\"\" u 1:14 w linespoint lt 1" << endl;
        
        *fgnu << "set key" << endl;
        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset y2tics; set ytics" << endl;
        *fgnu << "unset y2label; set ylabel \"Radial velocity precision (m/s)\"" << endl;
        *fgnu << "set origin 0.75*DX,DY" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:13 t \"full set of lines\" w linespoint pt 7 lt 3 lw 2, \"\" u 1:12 t \"clean sample\" w linespoint pt 7 lt 4 lw 2" << endl;
        *fgnu << "unset key" << endl;
        
        *fgnu << "unset x2label; set xlabel \"order number\"" << endl;
        *fgnu << "unset y2tics; unset y2label" << endl;
        *fgnu << "set origin DX+SX,DY" << endl;
        *fgnu << "set ytics nomirror" << endl;
        *fgnu << "set ylabel \"% matching lines\" offset +1.5,0l" << endl;
        *fgnu << "set key bottom" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:10 t \"% of comparison lines\" w linespoint pt 7 lt 4 lw 1, \"\" u 1:11 t \"% of atlas lines\" w linespoint pt 6 lt 4 lw 1" << endl;
        
        *fgnu << "unset key" << endl;
        *fgnu << "unset ytics; set y2tics" << endl;
        *fgnu << "unset ylabel; set y2label \"Number of matching lines\"" << endl;
        *fgnu << "plot \"" << dataFileName << "\" u 1:9 w linespoint lt 2 lw 3" << endl;
        
        *fgnu << endl;
        
        *fgnu << "unset multiplot" << endl;
        
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
 * Generate 2D plot for spectra of atlas + comparison + identified lines
 */
void GenerateWavelengthSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, string linesdatafilename, bool subtractCentralWavelength, bool display) {
    
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
    if(subtractCentralWavelength) {
        *fgnu << "\nset xlabel \"{/Symbol l} - {/Symbol l}_c (nm)\"" << endl;
    } else {
        *fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
    }
    *fgnu << "set ylabel \"order number + norm flux\"" << endl;
        
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        *fgnu << endl;

        if(subtractCentralWavelength) {
            *fgnu << "plot \"" << atlasdatafilename << "\" u ($2-$4):($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$5):($1 + $4) w l lt 3, ";
            *fgnu << "\"" << linesdatafilename << "\" u ($3-$6):($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
        } else {
            *fgnu << "plot \"" << atlasdatafilename << "\" u 2:($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u 3:($1 + $4) w l lt 3, ";
            *fgnu << "\"" << linesdatafilename << "\" u 3:($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
        }
        
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
        
        if(subtractCentralWavelength) {
            *fgnu << "plot \"" << atlasdatafilename << "\" u ($2-$4):($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$5):($1 + $4) w l lt 3, ";
            *fgnu << "\"" << linesdatafilename << "\" u ($3-$6):($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
        } else {
            *fgnu << "plot \"" << atlasdatafilename << "\" u 2:($1 + $3) w l lt 4, ";
            *fgnu << "\"" << compdatafilename << "\" u 3:($1 + $4) w l lt 3, ";
            *fgnu << "\"" << linesdatafilename << "\" u 3:($1 + $4):5 w xerr pt 7 lw 2 lt 1" << endl;
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

/*
 * Generate 3D plot for spectra of atlas + comparison + identified lines 
 */
void GenerateWavelength3DSpecPlot(string gnuScriptFileName, string outputPlotEPSFileName, string atlasdatafilename, string compdatafilename, bool subtractCentralWavelength, bool display) {
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
    *fgnu << "set cblabel \"normalized flux\"" << endl;
    
    if(subtractCentralWavelength) {
        *fgnu << "\nset xlabel \"{/Symbol l} - {/Symbol l}_c (nm)\"" << endl;
    } else {
        *fgnu << "\nset xlabel \"{/Symbol l} (nm)\"" << endl;
    }
    *fgnu << "set ylabel \"order number\"" << endl;
    
    *fgnu << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced mono solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        if(subtractCentralWavelength) {
            *fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2-$5):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;
        } else {
            *fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;

        }
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
        
        if(subtractCentralWavelength) {
            *fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2-$5):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3-$7):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;
        } else {
            *fgnu << "\nsplot \"" << atlasdatafilename << "\" u ($2):($1 + $4/4 - 0.5 + 0.25 + 0.25 + 0.175):3 w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.25 + 0.125):(1-$5) w pm3d, ";
            *fgnu << "\"" << compdatafilename << "\" u ($3):($1 + $6/4 - 0.5 + 0.075):4 w pm3d" << endl;
			
        }
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
 * Generate plot for wavelength solution - NOT IMPLEMENTED YET
 */
void GenerateWavelengthSolutionPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display)
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
    *fgnu << "\nset xlabel \"image rows (pixels)\"" << endl;
    *fgnu << "set ylabel \"image cols (pixels)\"" << endl;
    
    *fgnu << "set pointsize 0.5" << endl;
    
//    *fgnu << "set xrange[" << row0 << ":" << rowf << "]" << endl;
//    *fgnu << "set yrange[" << col0 << ":" << colf << "]" << endl;
    
    for(unsigned k=0;k<npolynomials;k++) {
        *fgnu << "poly" << k;
        polynomials[k]->printEquation(fgnu);
    }
    
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nplot \"" << dataFileName << "\" u 4:3 w p pt 7";
        
        for(unsigned k=0;k<npolynomials;k++) {
            *fgnu << ", poly" << k << "f(x)";
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
        *fgnu << "\nplot \"" << dataFileName << "\" u 4:3 w p pt 7";
        
        for(unsigned k=0;k<npolynomials;k++) {
            *fgnu << ", poly" << k << "f(x)";
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


/*
 * get the raw lines, between wl0 and wlf
 */
unsigned getRawLines(operaSpectralOrder *spectralOrder, double *rawlinecenter, double *rawlinecenterError, double *rawlineflux, double *rawlinesigma) {
	unsigned rawlinesinorder = 0;
	
	if (spectralOrder->gethasSpectralLines()) {
		operaSpectralLines *spectralLines = spectralOrder->getSpectralLines();
		unsigned nFeatures = spectralLines->getNFeatures();
		for (unsigned featurenumber=0;featurenumber<nFeatures;featurenumber++) {
			operaSpectralFeature *spectralFeature = spectralLines->getSpectralFeature(featurenumber);
			double *center = spectralFeature->getGaussianFit()->getCenterVector();
			double *sigma = spectralFeature->getGaussianFit()->getSigmaVector();
			double *amplitude = spectralFeature->getGaussianFit()->getAmplitudeVector();        
			double *centerError = spectralFeature->getGaussianFit()->getCenterErrorVector();
			//double *sigmaError = spectralFeature->getGaussianFit()->getSigmaErrorVector();
			//double *amplitudeError = spectralFeature->getGaussianFit()->getAmplitudeErrorVector(); 
			
			for (unsigned line=0; line<spectralFeature->getnLines(); line++) {
				rawlinecenter[rawlinesinorder] = *center++;
                rawlinecenterError[rawlinesinorder] = *centerError++;
				rawlineflux[rawlinesinorder] = *amplitude++;
                rawlinesigma[rawlinesinorder] = *sigma++;
				rawlinesinorder++;
			}
		}
	}
	return rawlinesinorder;
}
/*
 * get a subset of the thatlas lines for this order only, between wl0 and wlf
 */ 
unsigned getThoriumArgonAtlasRange(double wl0, double wlf, double **thwl, double **thi) {
	unsigned firstline = 0;
	unsigned line = 0;
	
	for (line=0; line<thatlaslines; line++) {
		if (thAtlasWavelength[line] >= wl0) {
			if (firstline == 0) {
				*thi = &thAtlasIntensity[line];
				*thwl = &thAtlasWavelength[line];
				firstline = line;
			}
			if (thAtlasWavelength[line] > wlf)
				break;
		}
	}
	if (line)
		line--;
	if (line > firstline) {
		return (line-firstline);
	} else {
		return 0;
	}
}

/*
 * Read the entire thorium argon atlas
 * and normalize the results
 */

unsigned readThoriumArgonAtlas(string atlas_lines, double *thAtlasWavelength, double *thAtlasIntensity) {
	igzstream astream;
	string dataline;
	double tmp = -1.0; 
	double tmpi = -1.0; 
	unsigned line = 0;
	
	astream.open(atlas_lines.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					char buff[8];
					sscanf(dataline.c_str(), "%*f %lf %lf %s", &tmp, &tmpi, buff);
                    
					if (strcmp(buff, "Th") && strcmp(buff, "Ar")) {
						if (DEBUG)
							printf("Skipping non-thorium-argon atlas entry %s.\n", buff); 
					} else {
						tmp *= 0.1;
						tmpi = pow(10, tmpi);
                        
                        thAtlasWavelength[line] = tmp;
                        thAtlasIntensity[line] = tmpi;
                        line++;  
                    }
 				}	// skip comments
			}	// if strlen
		} // while (astream.good())
		line--;
		if (line > 0) {
			if (verbose) {
				printf("          [Atlas] %d lines found wl0=%.2f wlc=%.2f wlf=%.2f\n", line, thAtlasWavelength[0], thAtlasWavelength[line/2], thAtlasWavelength[line-1]);
			}
		} else {
			if (verbose) {
				printf("          [Atlas] no lines found in atlas.\n");
			}
		}
		astream.close();
	}	// if (astream.open()
	return line;
}

/*
 * Read the the full atlas spectrum
 */
unsigned readAtlasSpectrum(string atlas_spectrum, double *atlasWavelength, double *atlasIntensity, double *atlasVariance) {
	igzstream astream;
	string dataline;
    
	double tmpwl = -1.0; 
	double tmpi = -1.0; 
	double tmp1 = -1.0; 
	double tmp2 = -1.0; 
	double tmpvar = -1.0; 
	unsigned np = 0;
	
	astream.open(atlas_spectrum.c_str());
	if (astream.is_open()) {
		while (astream.good()) {
			getline(astream, dataline);
			if (strlen(dataline.c_str())) {
				if (dataline.c_str()[0] == '#') {
					// skip comments
				} else {
					sscanf(dataline.c_str(), "%lf %lf %lf %lf %lf", &tmpwl, &tmpi, &tmp1, &tmp2, &tmpvar);
                    
                    atlasWavelength[np] = 0.1*tmpwl;
                    atlasIntensity[np] = tmpi;
                    atlasVariance[np] = tmpvar;
                    np++;  
                }	// skip comments
            }
		} // while (astream.good())
		
		if (np > 0) {
			if (verbose) {
				printf("          [Atlas] %d points found wl0=%.2f wlc=%.2f wlf=%.2f\n", np, atlasWavelength[0], atlasWavelength[np/2], atlasWavelength[np-1]);
			}
		} else {
			printf("          [Atlas] no points found in atlas.\n");
		}
		astream.close();
	}	// if (astream.open()
	return np;
}

/*
 * get a subset of the atlas spectrum for this order only, between wl0 and wlf
 */
unsigned getAtlasSpectrumRange(double wl0, double wlf, double **thwl, double **thi, double **thvar) {
	unsigned firstline = 0;
	unsigned np = 0;
    
	for (np=0; np<nPointsInAtlasSpectrum; np++) {
		if (atlasWavelength[np] >= wl0) {
			if (firstline == 0) {
				*thi = &atlasIntensity[np];
				*thwl = &atlasWavelength[np];
                *thvar = &atlasVariance[np];
				firstline = np;
			}
			if (atlasWavelength[np] > wlf)
				break;
		}
	}
	if (np)
		np--;
	if (np > firstline) {
		return (np-firstline);
	} else {
		return 0;
	}
}


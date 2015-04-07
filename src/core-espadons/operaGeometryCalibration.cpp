/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                    ****
 *******************************************************************
 Module name: operaGeometryCalibration
 Version: 1.0
 Description: Finds the location of orders.
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
#include "core-espadons/operaGeometryCalibration.h"

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
#include "libraries/ladfit.h" // for ladfit_d

#define MINIMUMORDERTOCONSIDER 15

/*! \file operaGeometryCalibration..cpp */

using namespace std;

/*! 
 * operaGeometryCalibration
 * \author Doug Teeple
 * \brief CCD geometry calculations.
 * \arg argc
 * \arg argv
 * \note --outputGeomFile=...
 * \note --masterbias=...
 * \note --masterflat=...
 * \note --badpixelmask=...
 * \note --subformat=...
 * \note --aperture=...
 * \note --detectionMethod=...
 * \note --FFTfilter=...
 * \note --referenceOrderSamplePosition=...
 * \note --minordertouse=...
 * \note --orderOfTracingPolynomial=...
 * \note --colDispersion=...
 * \note --invertOrders=...
 * \note --binsize=...
 * \throws operaException cfitsio error code
 * \throws operaException operaErrorNoInput
 * \throws operaException operaErrorNoOuput
 * \return EXIT_STATUS
 * \ingroup core
 */

int main(int argc, char *argv[])
{
	int opt;
    
	string outputGeomFile;
	string masterbias;
	string masterflat;
	string badpixelmask;
	string inputGainFile;
	string inputOrderSpacing;
    
	struct subformat {
		unsigned x1, x2;
		unsigned y1, y2;
	} subformat = {8, 2040, 3, 4600 };
    
    unsigned referenceOrderSamplePosition = 1;  // This position is with respect to the dispersion direction (rows for Espadons)
    
	unsigned minordertouse = 0;
	unsigned maxorders = MAXORDERS;
	unsigned orderOfTracingPolynomial = 4;
    
    bool recenterIPUsingSliceSymmetry = TRUE; // This flag will allow the spatial profile to be recentered such as to divide half number of slices on each side
    unsigned totalNumberOfSlices = 6; // Number of slices used in routine to recenter spatial profile
    
	dispersiondirection_t colDispersion = up;
	bool invertOrders = true;
	bool FFTfilter = false;
	bool witherrors = false;
    bool graces = false;
	int aperture = 20;
	int detectionMethod = 1; // 1. Gaussian, 2. IP, 3. Top-hat
	unsigned binsize = 1;
    double gain = 1;
    double noise = 5;
    unsigned nsamples = 3;
    string plotfilename;
	string datafilename;
	string scriptfilename;
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
    
	struct option longopts[] = {
		{"outputGeomFile",1, NULL, 'o'},
		{"masterbias",1, NULL, 'b'},
		{"masterflat",1, NULL, 'f'},
		{"badpixelmask",1, NULL, 'm'},
        {"inputGainFile",1, NULL, 'g'},
		{"inputOrderSpacing",1, NULL, 'c'},
        
		{"subformat",1, NULL, 's'},
		{"referenceOrderSamplePosition",1, NULL, 'Y'},  // This position is with respect to the dispersion direction (rows for Espadons)
        {"minordertouse",1, NULL, 'i'},
		{"maxorders",1, NULL, 'X'},
		{"orderOfTracingPolynomial",1, NULL, 'N'},
		{"recenterIPUsingSliceSymmetry",1, NULL, 'y'},
		{"totalNumberOfSlices",1, NULL, 'C'},
        {"detectionMethod",1, NULL, 'M'},
		{"FFTfilter",1, NULL, 'R'},
		{"graces",1, NULL, 'G'},
		{"colDispersion",1, NULL, 'D'},
		{"aperture",1, NULL, 'A'},
		{"invertOrders",1, NULL, 'V'},
		{"binsize",1, NULL, 'B'},
		{"witherrors",0, NULL, 'E'},
        {"nsamples",1, NULL, 'a'},
		
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
	
	while((opt = getopt_long(argc, argv, "o:b:f:m:g:c:s:Y:i:X:N:y:C:M:R:G:D:A:V:B:E:a:P:F:S:I::p::v::d::t::h",
							 longopts, NULL))  != -1)
	{
		
		switch(opt)
		{
			case 'o':		// outputGeomFile
				outputGeomFile = optarg;
				break;
			case 'b':
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
			case 'c':
				inputOrderSpacing = optarg; // input order spacing calibration file
				break;
			case 's':		// image subformat
				if (strlen(optarg))
					sscanf(optarg, "%u %u %u %u", &subformat.x1, &subformat.x2, &subformat.y1, &subformat.y2);
				break;
            case 'Y':
				referenceOrderSamplePosition = atoi(optarg);
				break;
			case 'i':
				minordertouse = atoi(optarg);
				break;
			case 'X':
				maxorders = atoi(optarg);
				break;
            case 'N':
				orderOfTracingPolynomial = atoi(optarg);
				break;
            case 'y':
				recenterIPUsingSliceSymmetry = (atoi(optarg)?true:false);
				break;
            case 'C':
				totalNumberOfSlices = atoi(optarg);
				break;
			case 'M':
				detectionMethod = atoi(optarg); // 1. Gaussian, 2. IP, 3. Top-hat
				break;
			case 'R':
				FFTfilter = (atoi(optarg)?true:false);
				break;
			case 'G':
				graces = (atoi(optarg)?true:false);
				break;
			case 'D':
				colDispersion = (atoi(optarg)?up:down);
				break;
			case 'A':		// aperture in pixels
				aperture = atoi(optarg);
				break;
			case 'V':
				invertOrders = (atoi(optarg)?true:false);
				break;
			case 'B':
				binsize = atoi(optarg);
				break;
			case 'E':
				witherrors = (atoi(optarg)?true:false);
				break;
			case 'a':
				nsamples = atoi(optarg);
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
			throw operaException("operaGeometryCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need some place to put data...
		if (outputGeomFile.empty()) {
			throw operaException("operaGeometryCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
		// we need an input order spacing calibration file
		if (inputOrderSpacing.empty()) {
			throw operaException("operaGeometryCalibration: ", operaErrorNoOutput, __FILE__, __FUNCTION__, __LINE__);
		}
        
		if (verbose) {
			cout << "operaGeometryCalibration: outputGeomFile = " << outputGeomFile << endl;
			cout << "operaGeometryCalibration: masterbias = " << masterbias << endl;
			cout << "operaGeometryCalibration: masterflat = " << masterflat << endl;
			cout << "operaGeometryCalibration: badpixelmask = " << badpixelmask << endl;
			cout << "operaGeometryCalibration: inputGainFile = " << inputGainFile << endl;
			cout << "operaGeometryCalibration: inputOrderSpacing = " << inputOrderSpacing << endl;
			cout << "operaGeometryCalibration: subformat = " << subformat.x1 << " " << subformat.x2 << " "<< subformat.y1  << " "<< subformat.y2 << "\n";
			cout << "operaGeometryCalibration: referenceOrderSamplePosition = " << referenceOrderSamplePosition << endl;
			cout << "operaGeometryCalibration: minordertouse = " << minordertouse << endl;
			cout << "operaGeometryCalibration: maxorders = " << maxorders << endl;
			cout << "operaGeometryCalibration: orderOfTracingPolynomial = " << orderOfTracingPolynomial << endl;
			cout << "operaGeometryCalibration: recenterIPUsingSliceSymmetry = " << recenterIPUsingSliceSymmetry << endl;
			cout << "operaGeometryCalibration: totalNumberOfSlices = " << totalNumberOfSlices << endl;
			cout << "operaGeometryCalibration: detectionMethod = " << detectionMethod << endl;
			cout << "operaGeometryCalibration: FFTfilter = " << FFTfilter << endl;
			cout << "operaGeometryCalibration: graces = " << graces << endl;
			cout << "operaGeometryCalibration: colDispersion = " << colDispersion << endl;
			cout << "operaGeometryCalibration: aperture = " << aperture << endl;
			cout << "operaGeometryCalibration: invertOrders = " << invertOrders << endl;
			cout << "operaGeometryCalibration: binsize = " << binsize << endl;
			cout << "operaGeometryCalibration: witherrors = " << witherrors << endl;
			cout << "operaGeometryCalibration: nsamples = " << nsamples << endl;
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
        
        operaFITSImage *flat = new operaFITSImage(masterflat, tfloat, READONLY);
        
        operaFITSImage *bias = NULL;
        
		if (!masterbias.empty()){
			bias = new operaFITSImage(masterbias, tfloat, READONLY);
		} else {
            bias = new operaFITSImage(flat->getnaxis1(),flat->getnaxis2(),tfloat);
            *bias = 0.0;
        }
        //---
        
        /*
         * Remove bias from masterflat
         */
		(*flat) -= (*bias);
        if(bias) {
            delete bias;
		}
        //---
        
        /*
         * Open badpixel mask
         */
        operaFITSImage *badpix = NULL;
        
		if (!badpixelmask.empty()){
			badpix = new operaFITSImage(badpixelmask, tfloat, READONLY);
		} else {
            badpix = new operaFITSImage(flat->getnaxis1(),flat->getnaxis2(),tfloat);
            *badpix = 1.0;
        }
		//---
        
		if (verbose) {
			cout << "operaGeometryCalibration: x1,y1,nx,ny = " << x1 << ' ' << y1 << ' ' << nx  << ' ' << ny << '\n';
		}
		
        /*
         * Setup slit aperture, and create spectral order vector object.
         */
        double slit = (double)aperture;
        unsigned uslit = (unsigned)aperture;
        
		operaSpectralOrderVector spectralOrders(MAXORDERS, ny, ny, 0);
		//---
        
        /*
         * read order spacing polynomial
         */
        spectralOrders.readOrderSpacingPolynomial(inputOrderSpacing);
        unsigned npars = spectralOrders.getOrderSpacingPolynomial()->getOrderOfPolynomial();
        double *par = spectralOrders.getOrderSpacingPolynomial()->getVector();
        //---
        
		/*
		 * read gain and noise
		 */
        unsigned amp = 0;
        if (!inputGainFile.empty()) {
            spectralOrders.readGainNoise(inputGainFile);
            gain = spectralOrders.getGainBiasNoise()->getGain(amp);
            noise = spectralOrders.getGainBiasNoise()->getNoise(amp);
		}
		if (verbose)
			cout << "operaGeometryCalibration: gain="<< gain << " noise=" << noise << endl;
        //---
        
        /*
         * Set up variables in Geometry class
         */
        for (unsigned order=0; order<MAXORDERS; order++) {
			spectralOrders.GetSpectralOrder(order)->getGeometry()->setapertureWidth(aperture);
			spectralOrders.GetSpectralOrder(order)->getGeometry()->setdispersionDirection(colDispersion);
			spectralOrders.GetSpectralOrder(order)->getGeometry()->setNumberofPointsToBinInYDirection(binsize);
            spectralOrders.GetSpectralOrder(order)->getGeometry()->setYmin(subformat.y1);
            spectralOrders.GetSpectralOrder(order)->getGeometry()->setYmax(subformat.y2);
		}
		//---
        
        /*
         * Measure instrument spatial profile from several samples of rows
         */
        float *ipfunc = new float[uslit];
        float *ipx = new float[uslit];
        float *iperr = new float[uslit];
        
        spectralOrders.measureIPAlongRowsFromSamples(*flat, *badpix, slit, nsamples, FFTfilter,(float)gain, (float)noise, subformat.x1, subformat.x2, subformat.y1, subformat.y2, ipfunc, ipx, iperr);
		//---
        
        /*
         * Apply centering corrections to the measured spatial profile
         */
        bool applyXcenterCorrection = TRUE;
        float xnewcenter = calculatePhotoCenter(uslit, ipx, ipfunc, iperr, applyXcenterCorrection);
        
        if(recenterIPUsingSliceSymmetry) {
            xnewcenter = calculateCenterOfSymmetry(uslit, ipx, ipfunc, iperr, totalNumberOfSlices,applyXcenterCorrection);
#ifdef PRINT_DEBUG
            for(unsigned i=0;i<uslit;i++) {
                cout << ((float)i - slit/2) <<  " " << ipx[i] <<  " " << ipfunc[i] << " " << iperr[i] << endl;
            }
#endif
        }
        
#ifdef PRINT_DEBUG
        for(unsigned i=0;i<uslit;i++) {
            cout << ((float)i - slit/2) <<  " " << ipx[i] <<  " " << ipfunc[i] << " " << iperr[i] << endl;
        }
#endif
        //---
        
        /*
         * Set up memory for variables, binsize and number of samples
         */
		unsigned NumberOfySamples;
		unsigned NumberofPointsToBinInYDirection = binsize;
        unsigned NumberofPointsToBinInYForReference = 2*NumberofPointsToBinInYDirection;

        float *fxref = new float[nx];
        float *fyref = new float[ny];
        float xmean[MAXORDERS],xmeanerr[MAXORDERS],ymean[MAXORDERS];

		NumberOfySamples = (unsigned)ceil((float)(ny-y1)/(float)NumberofPointsToBinInYDirection);
        //---
        
        /*
         * First pass only using rows around reference row to figure out a reference order map.
         */
        // Note: for the reference sample we use twice more points to obain higher SNR
        unsigned firstY = y1;
        unsigned lastY = ny;
        if(referenceOrderSamplePosition >= NumberofPointsToBinInYForReference/2) {
            firstY = referenceOrderSamplePosition - NumberofPointsToBinInYForReference/2;
        }
        if(referenceOrderSamplePosition + NumberofPointsToBinInYForReference/2 < ny) {
            lastY = referenceOrderSamplePosition + NumberofPointsToBinInYForReference/2;
        }
        
        float y = 0;
        unsigned np = getRowBinnedData(flat,x1,nx,nx,firstY,lastY,ny,fxref,fyref,&y,FFTfilter);
        
        if(debug) {
            for (unsigned i=0;i<np;i++) {
                cout << fxref[i] << " " << fyref[i] << endl;
            }
        }
        
        unsigned nords = geometryDetectOrders(np,fxref,fyref,slit,ipfunc,binsize,(float)noise,(float)gain,xmean,ymean,xmeanerr,detectionMethod, witherrors, graces);

        float *xref = new float[MAXORDERS];
        float *yref = new float[MAXORDERS];
        float *xreferr = new float[MAXORDERS];
        int *AbsRefOrdNumber = new int[MAXORDERS];
        
        unsigned nrefs = operaCCDDetectOrderMapBasedOnSpacingPolynomial(np,fxref,fyref,uslit,ipfunc,ipx,slit,(float)noise,(float)gain,npars,par,nords,xmean,ymean,xmeanerr,xref,yref,xreferr,AbsRefOrdNumber,minordertouse,maxorders);
        // -- End of reference map
        if(debug) {
            for (unsigned i=0;i<nrefs;i++) {
                cout << AbsRefOrdNumber[i] << " " << xref[i] << " " << yref[i] << " " << xreferr[i] << endl;
            }
        }
        //---
        
        /*
         * Create memory space to save data obtained from samples
         */
        
        float *fx = new float[nx];
        float *fy = new float[ny];
        
        float *xord_tmp = new float[MAXORDERS];
        float *yord_tmp = new float[MAXORDERS];
        float *xerrord_tmp = new float[MAXORDERS];
        int *AbsOrdNumber_tmp = new int[MAXORDERS];
        
        float *ypos = new float[NumberOfySamples];
        unsigned *newnords = new unsigned[NumberOfySamples];

        float *xords[NumberOfySamples];
        float *yords[NumberOfySamples];
        float *xerrords[NumberOfySamples];
        int *AbsOrdNumbers[NumberOfySamples];
        
        for(unsigned k=0;k<NumberOfySamples;k++){
            xords[k] = new float[MAXORDERS];
            yords[k] = new float[MAXORDERS];
            xerrords[k] = new float[MAXORDERS];
            AbsOrdNumbers[k] = new int[MAXORDERS];
            newnords[k] = 0;
        }
        //---
        
        /*
         * Once first order is figured out, then we save the data for each order. Note that
         * the absolute order number is important to make sure the data are always associated
         * to the correct order number. 
         */
        
        //-- Figure out which bin contains the reference order
        unsigned kref = (unsigned)round(float(referenceOrderSamplePosition - y1)/(float)NumberofPointsToBinInYDirection);
        
        //-- Start detecting orders in samples ABOVE reference row
        for (unsigned i=0;i<nrefs;i++) {
            xord_tmp[i] = xref[i];
            yord_tmp[i] = yref[i];
            xerrord_tmp[i] = xreferr[i];
            AbsOrdNumber_tmp[i] = AbsRefOrdNumber[i];
        }
        for(unsigned k=kref;k<NumberOfySamples;k++){
            unsigned firstY = y1 + NumberofPointsToBinInYDirection*(k);
            unsigned lastY =  y1 + NumberofPointsToBinInYDirection*(k+1);
            if(lastY >= ny) break;
            
            unsigned np = getRowBinnedData(flat,x1,nx,nx,firstY,lastY,ny,fx,fy,&ypos[k],FFTfilter);
            unsigned nords = geometryDetectOrders(np,fx,fy,slit,ipfunc,binsize,(float)noise,(float)gain,xmean,ymean,xmeanerr,detectionMethod, witherrors, graces);
            
            newnords[k] = operaCCDDetectMissingOrdersUsingNearMap(np,fx,fy,uslit,ipfunc,ipx,slit,(float)noise,(float)gain,npars,par,nords,xmean,ymean,xmeanerr,nrefs,xord_tmp,yord_tmp,AbsOrdNumber_tmp,xords[k],yords[k],xerrords[k],AbsOrdNumbers[k]);
            
            if(debug) {
                for (unsigned i=0;i<newnords[k];i++) {
                    cout << xords[k][i] << " " << yords[k][i] << " " << xerrords[k][i]<< endl;
                }
            }

            //-- save current sample in the tmp to be used in the next loop around
            for (unsigned i=0;i<nrefs;i++) {
                xord_tmp[i] = xords[k][i];
                yord_tmp[i] = yords[k][i];
                xerrord_tmp[i] = xerrords[k][i];
                AbsOrdNumber_tmp[i] = AbsOrdNumbers[k][i];
            }
        }
        
        //-- Then detect orders in samples BELOW reference row
        for (unsigned i=0;i<nrefs;i++) {
            xord_tmp[i] = xref[i];
            yord_tmp[i] = yref[i];
            xerrord_tmp[i] = xreferr[i];
            AbsOrdNumber_tmp[i] = AbsRefOrdNumber[i];
        }
        for(unsigned k=kref-1;k>0;k--){
            unsigned firstY = y1 + NumberofPointsToBinInYDirection*(k);
            unsigned lastY =  y1 + NumberofPointsToBinInYDirection*(k+1);
            if(lastY >= ny) break;
            
            unsigned np = getRowBinnedData(flat,x1,nx,nx,firstY,lastY,ny,fx,fy,&ypos[k],FFTfilter);
            unsigned nords = geometryDetectOrders(np,fx,fy,slit,ipfunc,binsize,(float)noise,(float)gain,xmean,ymean,xmeanerr,detectionMethod, witherrors, graces);
        
            newnords[k] = operaCCDDetectMissingOrdersUsingNearMap(np,fx,fy,uslit,ipfunc,ipx,slit,(float)noise,(float)gain,npars,par,nords,xmean,ymean,xmeanerr,nrefs,xord_tmp,yord_tmp,AbsOrdNumber_tmp,xords[k],yords[k],xerrords[k],AbsOrdNumbers[k]);
            
            //-- save current sample in the tmp to be used in the next loop around
            for (unsigned i=0;i<nrefs;i++) {
                xord_tmp[i] = xords[k][i];
                yord_tmp[i] = yords[k][i];
                xerrord_tmp[i] = xerrords[k][i];
                AbsOrdNumber_tmp[i] = AbsOrdNumbers[k][i];
            }
        }
        //---
        
        
        /*
         * Now put data into plot file and geometry class, so it can fit the tracing polynomials
         */
		for(unsigned k=0;k<NumberOfySamples;k++){
            if (fdata != NULL) {
                for(unsigned i=0;i<newnords[k];i++) {
                    *fdata << k << " " <<  AbsOrdNumbers[k][i] << " " << xords[k][i] << " " << ypos[k] << " " << xerrords[k][i] << " " << yords[k][i] << endl;
                }
            }

			unsigned lastorder = 0;
			unsigned order = 0;
			for (unsigned i=0;i<newnords[k];i++) {
				operaSpectralOrder *SpectralOrder = NULL;
				for (order=lastorder; order<MAXORDERS; order++) {
					if (spectralOrders.GetSpectralOrder(order)->getorder() == (unsigned)AbsOrdNumbers[k][i]) {
						SpectralOrder = spectralOrders.GetSpectralOrder(order);
						lastorder = order+1;
						break;
					}
				}
				// SpectralOrder should always be found, but...
				if (SpectralOrder) {
					SpectralOrder->getGeometry()->addOrderCenterValue(xords[k][i], ypos[k], yords[k][i], xerrords[k][i]);
				}
			}
        }
        //---
        
        /*
         * Now calculate the polynomials for each order
         */
        Polynomial *polynomials[MAXORDERS]; // for plot

        int *LastAbsOrdNumbers = new int[MAXORDERS];
        
		unsigned ordercount = 0;
		for (unsigned i=0;i<maxorders;i++) {
			double chisqr;	// NOTE: The chisqr is not used here, but it is saved in the polynomial itself
            
            LastAbsOrdNumbers[i] = minordertouse + (int)i;
            
			if (spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->getGeometry()->getNdatapoints() > orderOfTracingPolynomial) {	// i.e. we did add data to this order
                spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->getGeometry()->traceOrder(orderOfTracingPolynomial, chisqr, witherrors);		// i.e. fit polynomial to the data of a single order
				spectralOrders.setMaxorder(LastAbsOrdNumbers[i]);
				spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->sethasGeometry(true);
				
                if (ordercount == 0) {
					spectralOrders.setMinorder(LastAbsOrdNumbers[i]);
				}
                polynomials[ordercount] = spectralOrders.GetSpectralOrder(LastAbsOrdNumbers[i])->getGeometry()->getCenterPolynomial();
                ordercount++;
			}
        }
                
		spectralOrders.setCount(ordercount);
        //---
        if (fdata != NULL) {
            fdata->close();
            if (!scriptfilename.empty()) {
                GenerateGeometryPlot(scriptfilename, plotfilename, datafilename, ordercount, polynomials, interactive,subformat.x1,subformat.x2,subformat.y1,subformat.y2,LastAbsOrdNumbers);
            }
        }
		
        /*
         * Write out geometry calibration file
         */
		if (!outputGeomFile.empty()) {
			spectralOrders.WriteSpectralOrders(outputGeomFile, Geom);
		}
        //---
        
        if(badpix)
            delete badpix;
		flat->operaFITSImageClose();
		
        delete[] ipfunc;
        delete[] ipx;
        delete[] iperr;  
        delete[] fx;
        delete[] fy;
	}
	catch (operaException e) {
		cerr << "operaGeometryCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (string s) {
		cerr << "operaGeometryCalibration: " << s << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaGeometryCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" +
	" --outputGeomFile=<GEOM_FILE>"
	" --masterbias=<FITS_IMAGE>"
	" --masterflat=<FITS_IMAGE>"
	" --badpixelmask=<FITS_IMAGE>"
    " --inputGainFile=<GAIN_FILE>"
	" --inputOrderSpacing=<ORDP_FILE>"
    " --subformat=<\"UNS_VALUE UNS_VALUE UNS_VALUE UNS_VALUE\">"
    " --referenceOrderSamplePosition=<UNS_VALUE>"
    " --minordertouse=<UNS_VALUE>"
    " --maxorders=<UNS_VALUE>"
    " --orderOfTracingPolynomial=<UNS_VALUE>"
    " --recenterIPUsingSliceSymmetry=<BOOL>"
    " --totalNumberOfSlices=<UNS_VALUE>"
    " --detectionMethod=<UNS_VALUE>"
    " --FFTfilter=<BOOL>"
    " --colDispersion=<UNS_VALUE>"
    " --aperture=<UNS_VALUE>"
    " --invertOrders=<UNS_VALUE>"    
    " --binsize=<UNS_VALUE>"
    " --witherrors=<BOOL>"
    " --nsamples=<UNS_VALUE>"    
    " --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>"
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --masterbias=/Users/edermartioli//opera//calibrations/PolarData/masterbias_OLAPAa_pol_Normal.fits.fz --masterflat=/Users/edermartioli//opera//calibrations/PolarData/masterflat_OLAPAa_pol_Normal.fits.fz --badpixelmask=/Users/edermartioli/opera-1.0//config/badpix_olapa-a.fits.fz --inputGainFile=/Users/edermartioli//opera//calibrations/PolarData/OLAPAa_pol_Normal.gain.gz --inputOrderSpacing=/Users/edermartioli//opera//calibrations/PolarData/OLAPAa_pol_Normal.ordp.gz --subformat=\"8 2040 3 4600\" --aperture=26 --detectionMethod=2 --FFTfilter=0 --referenceOrderSamplePosition=2300 --nsamples=3 --maxorders=44 --minordertouse=17  --orderOfTracingPolynomial=4 --binsize=20 --colDispersion=1 --invertOrders=1  --outputGeomFile=OLAPAa_pol_Normal.geom.gz --scriptfilename=geom.gnu --datafilename=geom.dat --plotfilename=geom.eps\n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -o, --outputGeomFile=<GEOM_FILE>, Geometry output calibration file name\n"
	"  -b, --masterbias=<FITS_IMAGE>, FITS image with masterbias\n"
	"  -f, --masterflat=<FITS_IMAGE>, FITS image with masterflat\n"
	"  -m, --badpixelmask=<FITS_IMAGE>, FITS image with badpixel mask\n"
	"  -g, --inputGainFile=<GAIN_FILE>, Input noise/gain/bias file\n"
	"  -c, --inputOrderSpacing=<ORDP_FILE>, Order spacing input file name\n"
    "  -s, --subformat=<\"UNS_VALUE UNS_VALUE UNS_VALUE UNS_VALUE\">, Image subformat to be inspected\n"    
    "  -Y, --referenceOrderSamplePosition=<UNS_VALUE>, Detector position to pick samples. Position along the dispersion direction (rows for Espadons)\n"
    "  -i, --minordertouse=<UNS_VALUE>, Number of first useful order\n"
    "  -X, --maxorders=<UNS_VALUE>, Maximum number of orders to use\n"
    "  -N, --orderOfTracingPolynomial=<UNS_VALUE>, Degree of polynomial to trace order positions\n"
    "  -y, --recenterIPUsingSliceSymmetry=<BOOL>, Boolean to allow the IP to be recentered placing half number of slices on each side\n"
    "  -C, --totalNumberOfSlices=<UNS_VALUE>, Total number of slices. Used in routine to recenter spatial profile\n"
    "  -M, --detectionMethod=<UNS_VALUE>, Method for detecting orders\n"
    "                              Available options are = 1, 2, and  3, where: \n"
    "                              1. Gaussian (default)\n"
    "                              2. IP \n"
    "                              3. Top-hat\n"
    "  -D, --FFTfilter=<BOOL>, Activate Fourier smoothing filter\n"
    "  -R, --colDispersion=<UNS_VALUE>, Define dispersion direction: vertical (along cols) = 1; horizontal (along rows) = 2\n"
    "  -A, --aperture=<UNS_VALUE>, Aperture size in pixel units\n"
    "  -V, --invertOrders=<BOOL>, Select this option to invert the counting of order numbers\n"
    "  -B, --binsize=<UNS_VALUE>, Number of rows to bin in order to detect order positions\n"
    "  -E, --witherrors=<BOOL>, Use error bars for polynomial fit\n"
    "  -a, --nsamples=<UNS_VALUE>, Number of row samples for detecting orders\n"
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n"
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateGeometryPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display, unsigned col0, unsigned colf, unsigned row0, unsigned rowf, int AbsOrdNumber[])
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
    
    *fgnu << "set xrange[" << row0 << ":" << rowf << "]" << endl;
    *fgnu << "set yrange[" << col0 << ":" << colf << "]" << endl;
    
    
    for(unsigned k=0;k<npolynomials;k++) {
        *fgnu << "set label \"" << AbsOrdNumber[k] << "\" at " << rowf+50 << "," << polynomials[k]->Evaluate((double)rowf) << " rotate by 90 font \"Helvetica,4\"" << endl;
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
    
float calculateCenterOfSymmetry(unsigned np, float *ipx, float *ipfunc, float *iperr, unsigned totalNumberOfSlices, bool applyXcenterCorrection) {
    float xphotocenter = 0;
    
    unsigned sliceSize = (unsigned)ceil((float)(np/totalNumberOfSlices));
   
    float minflux = BIG;
    float minx = 0;
    
    for(unsigned i=(np-sliceSize)/2; i<(np+sliceSize)/2; i++) {
        if(ipfunc[i] < minflux){
            minflux = ipfunc[i];
            minx = ipx[i];
        }
    }

    if(minx) {
        xphotocenter = minx;
    }
    
    if(applyXcenterCorrection) {
        for(unsigned i=0; i<np; i++)
            ipx[i] -= xphotocenter;
    }
    return xphotocenter;
}

float calculatePhotoCenter(unsigned np, float *ipx, float *ipfunc, float *iperr, bool applyXcenterCorrection) {
    float xphotocenter = 0;
    float totalflux = 0;
    for(unsigned i=0; i<np; i++) {
        totalflux += ipfunc[i];
    }
    for(unsigned i=0; i<np; i++) {
        xphotocenter += ipx[i]*ipfunc[i]/totalflux;
    }
    if(applyXcenterCorrection) {
        for(unsigned i=0; i<np; i++)
            ipx[i] -= xphotocenter;
    }
    
    return xphotocenter;
}


unsigned getRowBinnedData(operaFITSImage *flat,unsigned x1,unsigned x2,unsigned nx,unsigned y1,unsigned y2,unsigned ny,float *fx,float *fy,float *yout, bool FFTfilter) {
    
    unsigned ybinsize = (unsigned)fabs(y2 - y1);
    
    float *fytmp = new float[ny];
    float *fysample = new float[ybinsize];

    unsigned np = 0;
    
    for (unsigned xx=x1; xx<x2; xx++) {
        unsigned ns=0;
        float ysample=0.0;
        
        for (unsigned yy=y1; yy<y2 ; yy++) {
            fysample[ns++] = (*flat)[yy][xx];
            ysample += (float)yy + 0.5;
        }
        ysample /= (float)ns;
        fx[np] = (float)xx + 0.5;
        
        if(FFTfilter){
            fytmp[np] = operaArrayMedian(ns,fysample);
        } else {
            fy[np] = operaArrayMedian(ns,fysample);
        }
        
#ifdef PRINT_DEBUG
        if(FFTfilter){
            cout << fx[np]  << " " << fytmp[np] << endl;
        } else {
            cout << fx[np]  << " " << fy[np] << endl;
        }
#endif
        *yout = ysample;
        np++;
    }
    
    if(FFTfilter){
        operaFFTLowPass(np,fytmp,fy,0.1);
    }
    
    delete[] fysample;
    delete[] fytmp;
    
    return np;
}

unsigned geometryDetectOrders(unsigned np,float *fx,float *fy,unsigned uslit,float *ipfunc, unsigned binsize, float noise,float gain,float *xmean,float *ymean,float *xmeanerr,int detectionMethod, bool witherrors, bool graces) {
    
    unsigned nords = 0;
    
    double slit = (double)uslit;
    double sigma = slit/4.0;
    double threshold = DETECTTHRESHOLD;
    
    if(detectionMethod == 1) {
        if (witherrors) {
            nords = operaCCDDetectPeaksWithErrorsUsingGaussian(np,fx,fy,sigma,(float)noise,(float)gain,threshold,xmean,ymean,xmeanerr);
        } else {
            nords = operaCCDDetectPeaksWithGaussian(np,fx,fy,sigma,(float)noise,(float)gain,threshold,xmean,ymean);
        }
    } else if (detectionMethod == 2) {
        if (witherrors) {
            nords = operaCCDDetectPeaksWithErrorsUsingIP(np,fx,fy,uslit,ipfunc,(float)noise,(float)gain,threshold/2,xmean,ymean,xmeanerr);
        } else {
            if(graces) {
                nords = operaCCDDetectPeaksByXCorrWithIP(np,fx,fy,uslit,ipfunc,(float)noise/sqrt((float)binsize),(float)gain,threshold,xmean,ymean);
            } else {
                /*
                 * The function below does not work on GRACES data. The one above should be better but
                 * it hasn't been tested yet for ESPaDOnS@CFHT. E. Martioli May 14 2014.
                 */
                nords = operaCCDDetectPeaksWithIP(np,fx,fy,uslit,ipfunc,(float)noise,(float)gain,threshold/2,xmean,ymean);
            }
        }
    } else if (detectionMethod == 3) {
        if (witherrors) {
            nords = operaCCDDetectPeaksWithErrorsUsingTopHat(np,fx,fy,uslit,(float)noise,(float)gain,threshold,xmean,ymean,xmeanerr);
        } else {
            nords = operaCCDDetectPeaksWithTopHat(np,fx,fy,uslit,(float)noise,(float)gain,threshold,xmean,ymean);
        }
    }

#ifdef PRINT_DEBUG
    for(unsigned i=0;i<nords;i++) {
        cout << xmean[i] << " " << y << " " << xmeanerr[i] << " " << ymean[i] << endl;
    }
#endif
    return nords;
}



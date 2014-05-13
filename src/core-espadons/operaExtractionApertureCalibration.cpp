/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaExtractionApertureCalibration
 Version: 1.0
 Description: Calibrate aperture for extraction
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Apr/2012
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
#include "core-espadons/operaExtractionApertureCalibration.h"

#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"			// for operaSpectralOrder
#include "libraries/operaSpectralOrderVector.h"		// for operaSpectralOrderVector
#include "libraries/operaSpectralElements.h"		// for operaSpectralOrder_t
#include "libraries/operaInstrumentProfile.h"		// for operaInstrumentProfile
#include "libraries/operaFITSImage.h"
#include "libraries/operaFITSSubImage.h"
#include "libraries/operaEspadonsImage.h"			// for imtype_t
#include "libraries/operaFITSProduct.h"	
#include "libraries/operaLib.h"						// systemf

#include "libraries/operaExtractionAperture.h"
#include "libraries/operaGeometricShapes.h"

#include "libraries/operaLibCommon.h"
#include "libraries/operaImage.h"
#include "libraries/operaStats.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/operaFit.h"	
#include "libraries/operaFFT.h"	
#include "libraries/ladfit.h"

#define NOTPROVIDED -999

/*! \file operaExtractionApertureCalibration.cpp */

using namespace std;

/*!
 * operaExtractionApertureCalibration
 * \author Eder Martioli
 * \brief Tool to calibrate the aperture for extraction.
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
	
    string inputorderspacing;
	string inputgeom; 	
	string inputprof;  	
	
	string outputApertureFile;
    
    int minorder = 22;
    bool minorderprovided = false;
    int maxorder = 62;    
    bool maxorderprovided = false;    
    
	int ordernumber = NOTPROVIDED;		
    unsigned pickImageRow = 0;
    unsigned nRowSamples = 1;
    
    unsigned xbin = 10;
    
    unsigned numberOfBeams = 1; // 1 for star-only;  2 for polar/s+s
	
    float gapBetweenBeams = 0;
    
    float apertureWidth = 26.0;
    float apertureHeight = 0.6;
    float backgroundAperture = 2.0;
    
    bool interactive = false;
    
	int debug=0, verbose=0, trace=0, plot=0;
	
    string plotfilename;	
	string datafilename;	
	string scriptfilename;	
    
	struct option longopts[] = {      
		{"outputApertureFile",1, NULL, 'o'},	
		{"inputorderspacing",1, NULL, 's'},          
		{"inputgeom",1, NULL, 'g'},
		{"inputprof",1, NULL, 'i'}, 
		{"ordernumber",1, NULL, 'O'},
		{"pickImageRow",1, NULL, 'R'},
		{"nRowSamples",1, NULL, 'L'},   
		{"xbin",1, NULL, 'x'}, 
		{"minorder",1, NULL, 'M'},
		{"maxorder",1, NULL, 'X'},          
		{"numberOfBeams",1, NULL, 'N'},	
		{"gapBetweenBeams",1, NULL, 'G'},	        
		{"apertureWidth",1, NULL, 'W'},	
		{"apertureHeight",1, NULL, 'H'},	
		{"backgroundAperture",1, NULL, 'B'},			
		{"plotfilename",1, NULL, 'P'},
		{"datafilename",1, NULL, 'F'},
		{"scriptfilename",1, NULL, 'S'},  
		{"interactive",0, NULL, 'I'}, 
		
		{"plot",optional_argument, NULL, 'p'},       
		{"verbose",optional_argument, NULL, 'v'},
		{"debug",optional_argument, NULL, 'd'},
		{"trace",optional_argument, NULL, 't'},
		{"help",no_argument, NULL, 'h'},
		{0,0,0,0}};
	
	while((opt = getopt_long(argc, argv, "o:s:g:i:O:R:L:x:M:X:N:G:W:H:B:P:F:S:I:v::d::t::p::h", 
							 longopts, NULL))  != -1)
	{
		switch(opt) 
		{             
			case 'o':
				outputApertureFile = optarg;
				break; 	                
			case 's':
				inputorderspacing = optarg;
				break;                
			case 'g':
				inputgeom = optarg;
				break;				
			case 'i':
				inputprof = optarg;
				break;					
			case 'O':
				ordernumber = atoi(optarg);
				break;	
			case 'R':
				pickImageRow = atoi(optarg);
				break;  
			case 'L':
				nRowSamples = atoi(optarg);
				break; 
			case 'x':
				xbin = atoi(optarg);
				break;
			case 'M':
				minorder = atoi(optarg);
                minorderprovided = true;
				break;  
			case 'X':
				maxorder = atoi(optarg);
                maxorderprovided = true;
				break;                  
			case 'N':
				numberOfBeams = atoi(optarg);
				break;	
			case 'G':		// width of gap betwenn beams
				gapBetweenBeams = atof(optarg);
				break;                 
			case 'W':		// aperture width for extraction in pixels
				apertureWidth = atof(optarg);
				break; 
			case 'H':		// aperture height for extraction in pixels
				apertureHeight = atof(optarg);
				break; 
			case 'B':		// aperture width for background subtraction in pixels
				backgroundAperture = atof(optarg);
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
	
	/*Start the module here*/
	
	try {
		// we need a geometry file...
		if (inputgeom.empty()) {
			throw operaException("operaExtractionApertureCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		// we need an instrument profile file...
		if (inputprof.empty()) {
			throw operaException("operaExtractionApertureCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		
		if (verbose) {
			cout << "operaExtractionApertureCalibration: outputApertureFile = " << outputApertureFile << endl;            
			cout << "operaExtractionApertureCalibration: inputorderspacing = " << inputorderspacing << endl;             
			cout << "operaExtractionApertureCalibration: inputgeom = " << inputgeom << endl; 
			cout << "operaExtractionApertureCalibration: inputprof = " << inputprof << endl;
            if(ordernumber != NOTPROVIDED) {
                cout << "operaExtractionApertureCalibration: ordernumber = " << ordernumber << endl;            
            }
            if(pickImageRow) {
                cout << "operaExtractionApertureCalibration: pickImageRow = " << pickImageRow << endl;       
            } else {
                cout << "operaExtractionApertureCalibration: nRowSamples = " << nRowSamples << endl;       
            }
            cout << "operaExtractionApertureCalibration: xbin = " << xbin << endl;
            
			cout << "operaExtractionApertureCalibration: numberOfBeams = " << numberOfBeams << endl;       
			cout << "operaExtractionApertureCalibration: apertureWidth = " << apertureWidth << endl;            
			cout << "operaExtractionApertureCalibration: apertureHeight = " << apertureHeight << endl;
			cout << "operaExtractionApertureCalibration: backgroundAperture = " << backgroundAperture << endl;            
            if(plot) {
                cout << "operaExtractionApertureCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaExtractionApertureCalibration: datafilename = " << datafilename << endl;
                cout << "operaExtractionApertureCalibration: scriptfilename = " << scriptfilename << endl; 
                if(interactive) {
                    cout << "operaExtractionApertureCalibration: interactive = YES" << endl; 
                } else {
                    cout << "operaExtractionApertureCalibration: interactive = NO" << endl; 
                }
            }            
		}
        ofstream *fdata = NULL;
        
        if (!datafilename.empty()) {
            fdata = new ofstream();
            fdata->open(datafilename.c_str());  
        }  
        
		operaSpectralOrderVector spectralOrders(inputgeom);
        spectralOrders.ReadSpectralOrders(inputprof);
        
        if(!inputorderspacing.empty()) {
            spectralOrders.ReadSpectralOrders(inputorderspacing);
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
		
        
        /*
         * Set aperture Rectangle aperture
         */   
        float width = apertureWidth;
        float height = apertureHeight;    
        
        float *tiltAngle = new float[MAXORDERS];
        float *FluxFraction = new float[MAXORDERS];     
        
        float *tiltWithinOrder = new float[nRowSamples];
        float *FluxFractionWithinOrder = new float[nRowSamples];        
        
        unsigned orderIndex=0;
        
        float *xphotocenter  = new float[nRowSamples*MAXORDERS];
        float *yphotocenter = new float[nRowSamples*MAXORDERS];
        float *photosum = new float[nRowSamples*MAXORDERS];
        unsigned nTotalPoints = 0;
        
		for (int order=minorder; order<=maxorder; order++) {	
            
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                    
			if (spectralOrder->gethasGeometry() && spectralOrder->gethasInstrumentProfile()) {
                
                operaGeometry *Geometry = spectralOrder->getGeometry();
                operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
                
                unsigned NXPoints = instrumentProfile->getNXPoints();
                unsigned NYPoints = instrumentProfile->getNYPoints();
                
                float ymiddle = (Geometry->getYmax() + Geometry->getYmin())/2;
                float ystep = (Geometry->getYmax() + Geometry->getYmin())/(float)nRowSamples;
                
                if(pickImageRow) {
                    nRowSamples = 1;
                }
        
                float normalizationFactor = 0;
                for (unsigned j=0; j<NYPoints; j++) {
                    for (unsigned i=0; i<NXPoints; i++) {
                        normalizationFactor += (float)instrumentProfile->getipDataFromPolyModel(ymiddle,i,j);
                        if(debug){ // for 3D plot
                            cout <<instrumentProfile->getIPixXCoordinate(i) << " "
                            << instrumentProfile->getIPixYCoordinate(j) << " "
                            << (float)instrumentProfile->getipDataFromPolyModel(ymiddle,i,j) << endl;
                        }
                    }
                    //  cout << endl; // for 3D plot
                }
                
                nTotalPoints = 0;
                
                for(unsigned k=0; k<nRowSamples; k++) {
                    
                    float ycenter = ymiddle - (float)(nRowSamples-1)*ystep/2 + ystep*(float)k + 0.5;
                    
                    if(pickImageRow) {
                        ycenter = pickImageRow + 0.5;
                    }
                    
                    unsigned nlocalpoints = 0;
                    double photosum_tmp = 0;
                    double yphotocenter_tmp = 0;
                    double xphotocenter_tmp = 0;
                    
                    for (unsigned i=xbin; i<NXPoints-xbin+1; i++) {
                        
                        float xcoord = instrumentProfile->getIPixXCoordinate(i);
                       
                        float localphotosum = 0;
                        for (unsigned j=0; j<NYPoints; j++) {
                            float ycoord = instrumentProfile->getIPixYCoordinate(j);
                            yphotocenter_tmp += ycoord*(float)instrumentProfile->getipDataFromPolyModel(ycenter,i,j);
                            localphotosum += (float)instrumentProfile->getipDataFromPolyModel(ycenter,i,j);
                        }
                                                
                        xphotocenter_tmp += localphotosum*xcoord;
                        
                        photosum_tmp += localphotosum;
                        
                        if(nlocalpoints == xbin) {
                            if(debug)
                                cout << order << " " << xphotocenter_tmp/photosum_tmp << " " << yphotocenter_tmp/photosum_tmp << " " << photosum_tmp << endl;
                            
                            xphotocenter[nTotalPoints] = xphotocenter_tmp;
                            yphotocenter[nTotalPoints] = yphotocenter_tmp;
                            photosum[nTotalPoints] = photosum_tmp;
                            nTotalPoints++;
                            
                            xphotocenter_tmp = 0;
                            yphotocenter_tmp = 0;
                            photosum_tmp = 0;
                            nlocalpoints = 0;
                        }
                        nlocalpoints++;
                    }
                }
                
                float am,bm,abdevm;
                
                ladfit(xphotocenter,yphotocenter,nTotalPoints,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
                
                tiltAngle[orderIndex]=180*atan(bm)/(M_PI);
                                
                // Set up extraction rectangle in the IP subpixel coordinate system
                Line extractionLine(bm,am,height,width);
                
                operaExtractionAperture lineAperture(&extractionLine,instrumentProfile,ymiddle);
                PixelSet *aperturePixels = lineAperture.getSubpixels();
                
                float FluxFractionCollected = 0;
                for(unsigned i=0; i<aperturePixels->getNPixels(); i++){
                    if(aperturePixels->getiIndex(i) >= 0 && aperturePixels->getiIndex(i) < (int)NXPoints &&
                       aperturePixels->getjIndex(i) >= 0 && aperturePixels->getjIndex(i) < (int)NYPoints ) {
                        FluxFractionCollected += (float)instrumentProfile->getipDataFromPolyModel(ymiddle,(unsigned)aperturePixels->getiIndex(i),(unsigned)aperturePixels->getjIndex(i));
                    }
                }
                
                FluxFraction[orderIndex] = FluxFractionCollected/normalizationFactor;
                               
                if(debug)
                    cout << order << " " << tiltAngle[orderIndex] << " " << 180*(abdevm/0.674433)/(M_PI) << " " << FluxFraction[orderIndex] << endl;
                
                if(verbose)
                    cout << "operaExtractionApertureCalibration: # order = " << order << " tilt = " << tiltAngle[orderIndex] << " +/- " << 180*(abdevm/0.674433)/(M_PI) << " FluxFraction=" << FluxFraction[orderIndex] << endl;

                spectralOrder->setTiltInDegrees(tiltAngle[orderIndex],180*(abdevm/0.674433)/(M_PI));
       
                orderIndex++;
       
            }
        }
        
        float tilt = operaArrayMedian(orderIndex,tiltAngle);
        float tilterror = operaArrayMedianSigma(orderIndex,tiltAngle,tilt);
        
        if(verbose)
            cout << "operaExtractionApertureCalibration: # Final tilt = " << tilt << " +/- " << tilterror << endl;
        
        float beamFluxFraction[MAXNUMBEROFBEAMS];
        unsigned outBoundPoints[MAXNUMBEROFBEAMS];
        operaPoint beamMidPoint(0,0);
        
        operaExtractionAperture *aperture[MAXNUMBEROFBEAMS];		// DO NOT DELETE
        operaExtractionAperture *leftBackgroundAperture = NULL;   	// DO NOT DELETE     
        operaExtractionAperture *rightBackgroundAperture = NULL;	// DO NOT DELETE
        float leftBackgroundFluxFraction;
        float rightBackgroundFluxFraction;        
        
        /*
         * Loop below will calculate apertures and save the information
         */ 
        unsigned ystackIndex = 0;
        orderIndex = 0;
        for (int order=minorder; order<=maxorder; order++) {
            
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
            
            if (spectralOrder->gethasGeometry() && spectralOrder->gethasInstrumentProfile()) {
                
                if(verbose)
                    cout << "operaExtractionApertureCalibration: Calculating the aperture of order " << order << endl;
                
                operaGeometry *Geometry = spectralOrder->getGeometry();
                operaInstrumentProfile *instrumentProfile = spectralOrder->getInstrumentProfile();
                
                spectralOrder->setnumberOfBeams(numberOfBeams);
                
                unsigned NXPoints = instrumentProfile->getNXPoints();		
                unsigned NYPoints = instrumentProfile->getNYPoints();
                
                /*
                 * Set up FITS image to store IP data - consider IP at center of order
                 */
                CMatrix ipImage = newCMatrix(NXPoints,NYPoints);        
                
                float ycenter;
                if(pickImageRow) {
                    ycenter = (float)pickImageRow + 0.5;
                } else {
                    ycenter = (Geometry->getYmax() + Geometry->getYmin())/2 + 0.5;
                }
                
                if((order - minorder) >= 5) {
                    ystackIndex++;
                    minorder = order;                        
                }
                
                float shiftX = (float)((order - minorder)*instrumentProfile->getxsize());
                float shiftY = (float)(ystackIndex*instrumentProfile->getysize());
                
                if(ordernumber != NOTPROVIDED) {
                    shiftX = 0;
                    shiftY = 0;
                }                        
                
                float normalizationFactor = 0;
                for (unsigned j=0; j<NYPoints; j++) {	 
                    for (unsigned i=0; i<NXPoints; i++) {
                        ipImage[j][i] = (float)instrumentProfile->getipDataFromPolyModel(ycenter,i,j);
                        normalizationFactor += (float)instrumentProfile->getipDataFromPolyModel(ycenter,i,j);
                    }
                }
                
                operaPoint MidPoint(0,0);

                Line extractionLine(tan(tiltAngle[orderIndex]*M_PI/180),height,width,MidPoint);
                
                float beamWidth = width/(float)numberOfBeams - gapBetweenBeams;
                float widthX = width*cos(tiltAngle[orderIndex]*M_PI/180);
                float aperXsize = widthX/(float)numberOfBeams;
                for(unsigned k=0;k<numberOfBeams;k++) {
                    beamFluxFraction[k]=0;
                    aperture[k] = NULL;
                }
                
                for(unsigned k=0;k<numberOfBeams;k++) {
                    
                    if(verbose)
                        cout << "operaExtractionApertureCalibration: beam " << k << endl;
                    
                    outBoundPoints[k] = 0;
                    
                    float beamXcenter = MidPoint.getXcoord() - (float)(numberOfBeams-1)*aperXsize/2 + aperXsize*(float)k;
                    float beamYcenter = extractionLine.getYcoord(beamXcenter);            
                    
                    beamMidPoint.setPoint(beamXcenter,beamYcenter);
                    Line beamExtractionLine((float)tan(tiltAngle[orderIndex]*M_PI/180),height,beamWidth,beamMidPoint);
                    
                    // NOTE: This should not be deleted as it is set into the spectralorder a a pointer
                    aperture[k] = new operaExtractionAperture(&beamExtractionLine,instrumentProfile,ycenter);                          
                    
                    PixelSet *aperturePixels = aperture[k]->getSubpixels();   
                    
                    for(unsigned i=0; i<aperturePixels->getNPixels(); i++){
                        if(aperturePixels->getiIndex(i) >= 0 && aperturePixels->getiIndex(i) < (int)NXPoints &&
                           aperturePixels->getjIndex(i) >= 0 && aperturePixels->getjIndex(i) < (int)NYPoints) {
                                
                                ipImage[(unsigned)aperturePixels->getjIndex(i)][(unsigned)aperturePixels->getiIndex(i)] = 0;
                                beamFluxFraction[k] += (float)instrumentProfile->getipDataFromPolyModel(ycenter,(unsigned)aperturePixels->getiIndex(i),(unsigned)aperturePixels->getjIndex(i))/normalizationFactor;
                            
                        } else {
                            outBoundPoints[k]++;
                        }
                    }  
                    
                    if(verbose)
                        cout << "operaExtractionApertureCalibration: # order = " << order << " beam = " << k << " FluxFraction = " << beamFluxFraction[k] << " # out-of-bound points = " << outBoundPoints[k] << endl;
                    
                    aperture[k]->setFluxFraction(beamFluxFraction[k]);
                    
					// DT May 14 2014, if we have no beam flux, we have no aperture...
					if (beamFluxFraction[k] > 0) {
						spectralOrder->setExtractionApertures(k,aperture[k]);
						spectralOrder->sethasExtractionApertures(true);           
					} else {
						spectralOrder->sethasExtractionApertures(false);           
					}
                }
                
                
                if(verbose)
                    cout << "operaExtractionApertureCalibration: Setting up background aperture " << endl;
                /*
                 * set up background aperture:
                 */
                float leftbkgXcenter = MidPoint.getXcoord() - (float)numberOfBeams*aperXsize/2 - backgroundAperture/2;
                float leftbkgYcenter = extractionLine.getYcoord(leftbkgXcenter);            
                
                operaPoint leftBackgroundMidPoint(leftbkgXcenter,leftbkgYcenter);
                Line leftBackgroundExtractionLine(tan(tiltAngle[orderIndex]*M_PI/180),height,backgroundAperture,leftBackgroundMidPoint);
                
                // NOTE: This should not be deleted as the pointer is set into the spectralorder
                leftBackgroundAperture = new operaExtractionAperture(&leftBackgroundExtractionLine,instrumentProfile,ycenter);                          
                PixelSet *leftaperturePixels = leftBackgroundAperture->getSubpixels();   
                
                leftBackgroundFluxFraction = 0;
                unsigned leftBackgroundoutBoundPoints = 0;            
                for(unsigned i=0; i<leftaperturePixels->getNPixels(); i++) {
                    if(leftaperturePixels->getiIndex(i) >= 0 && leftaperturePixels->getiIndex(i) < (int)NXPoints &&
                       leftaperturePixels->getjIndex(i) >= 0 && leftaperturePixels->getjIndex(i) < (int)NYPoints ) {                
                        ipImage[(unsigned)leftaperturePixels->getjIndex(i)][(unsigned)leftaperturePixels->getiIndex(i)] = 0;
                        leftBackgroundFluxFraction += (float)instrumentProfile->getipDataFromPolyModel(ycenter,leftaperturePixels->getiIndex(i),leftaperturePixels->getjIndex(i))/normalizationFactor;                
                    } else {
                        leftBackgroundoutBoundPoints++;
                    }
                }              
                if(verbose)
                    cout << "operaExtractionApertureCalibration: # order = " << order << " Left Background Flux Fraction = "<< leftBackgroundFluxFraction << " # out-of-bound points = " << leftBackgroundoutBoundPoints << endl;
                
                leftBackgroundAperture->setFluxFraction(leftBackgroundFluxFraction);
                spectralOrder->setBackgroundApertures(0,leftBackgroundAperture);
                
                float rightbkgXcenter = MidPoint.getXcoord() + (float)numberOfBeams*aperXsize/2 + backgroundAperture/2;
                float rightbkgYcenter = extractionLine.getYcoord(rightbkgXcenter);            
                
                operaPoint rightBackgroundMidPoint(rightbkgXcenter,rightbkgYcenter);
                Line rightBackgroundExtractionLine(tan(tiltAngle[orderIndex]*M_PI/180),height,backgroundAperture,rightBackgroundMidPoint);
                
                // NOTE: This should not be deleted as the pointer is set into the spectralorder
                rightBackgroundAperture = new operaExtractionAperture(&rightBackgroundExtractionLine,instrumentProfile,ycenter);                          
                PixelSet *rightaperturePixels = rightBackgroundAperture->getSubpixels();   
                
                rightBackgroundFluxFraction=0;
                unsigned rightBackgroundoutBoundPoints = 0;
                for(unsigned i=0; i<rightaperturePixels->getNPixels(); i++){
                    if(rightaperturePixels->getiIndex(i) >= 0 && rightaperturePixels->getiIndex(i) < (int)NXPoints &&
                       rightaperturePixels->getjIndex(i) >= 0 && rightaperturePixels->getjIndex(i) < (int)NYPoints ) {                
                        ipImage[(unsigned)rightaperturePixels->getjIndex(i)][(unsigned)rightaperturePixels->getiIndex(i)] = 0;
                        rightBackgroundFluxFraction += (float)instrumentProfile->getipDataFromPolyModel(ycenter,(unsigned)rightaperturePixels->getiIndex(i),(unsigned)rightaperturePixels->getjIndex(i))/normalizationFactor;                
                    } else {
                        rightBackgroundoutBoundPoints++;
                    }
                }              
                rightBackgroundAperture->setFluxFraction(rightBackgroundFluxFraction);
                spectralOrder->setBackgroundApertures(1,rightBackgroundAperture);
                
                if(verbose)
                    cout << "operaExtractionApertureCalibration: # order = " << order << " Right Background Flux Fraction = "<< rightBackgroundFluxFraction << " # out-of-bound points = " << rightBackgroundoutBoundPoints << endl << endl;            
                
                for (unsigned j=0; j<NYPoints; j++) {	 
                    for (unsigned i=0; i<NXPoints; i++) {
                        if (fdata != NULL) {
                            *fdata << order << " " 
                            << instrumentProfile->getIPixXCoordinate(i) + shiftX << " " 
                            << instrumentProfile->getIPixYCoordinate(j) + shiftY << " " 
                            << ipImage[j][i]/normalizationFactor << endl;
                        } else if (debug) {
                            cout << order << " " 
                            << instrumentProfile->getIPixXCoordinate(i) + shiftX << " " 
                            << instrumentProfile->getIPixYCoordinate(j) + shiftY << " " 
                            << ipImage[j][i]/normalizationFactor << endl;
                        }
                        
                    }
                    if (fdata != NULL) {
                        *fdata << endl;
                    } else {
                        if (debug)
                            cout << endl;
                    }                
                }
                if (fdata != NULL) {
                    *fdata << endl;
                } else {
                    if (debug)
                        cout << endl;
                }
                orderIndex++;
                deleteCMatrix(ipImage);
            } else {
                if(verbose)
                    cout << "operaExtractionApertureCalibration:aperture of order " << order << " skipped." << endl;
            }
        }        
        
        // output aperture to file...
        spectralOrders.WriteSpectralOrders(outputApertureFile,Aperture);		
        
        if (fdata != NULL) {
            fdata->close();
            if (!scriptfilename.empty()) {
                GenerateExtractionAperturePlot(scriptfilename,plotfilename,datafilename,interactive);
            }
            
        } 
 
        delete[] xphotocenter;
        delete[] yphotocenter;
        delete[] photosum;
        
        delete[] tiltAngle;
        delete[] FluxFraction;
        delete[] tiltWithinOrder;
        delete[] FluxFractionWithinOrder;
	}
        
        
	catch (operaException e) {
		cerr << "operaExtractionApertureCalibration: " << e.getFormattedMessage() << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cerr << "operaExtractionApertureCalibration: " << operaStrError(errno) << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/* Print out the proper program usage syntax */
static void printUsageSyntax(char * modulename) {
	cout <<
	"\n"
	" Usage: "+string(modulename)+"  [-vdth]" + 
	" --outputApertureFile=<APER_FILE>"	
	" --inputgeom=<GEOM_FILE>"
	" --inputprof=<PROF_FILE>"
	" --inputorderspacing=<ORDS_FILE>"
	" --ordernumber=<INT_VALUE>"
	" --pickImageRow=<UNS_VALUE>"   
	" --nRowSamples=<UNS_VALUE>"   
	" --xbin=<UNS_VALUE>"
	" --minorder=<INT_VALUE>"
	" --maxorder=<INT_VALUE>"
	" --numberOfBeams=<UNS_VALUE>"
    " --gapBetweenBeams=<FLT_VALUE>"    
	" --apertureWidth=<FLT_VALUE>"
	" --apertureHeight=<FLT_VALUE>"
	" --backgroundAperture=<FLT_VALUE>"
	" --plotfilename=<EPS_FILE>"
	" --datafilename=<DATA_FILE>"
	" --scriptfilename=<GNUPLOT_FILE>" 
	" --interactive=<BOOL>\n\n"
	" Example: "+string(modulename)+" --inputgeom=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.geom --inputprof=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.prof --outputApertureFile=/opera/calibrations/11AQ14-Jul08/OLAPAa_sp2_Normal.aper --numberOfBeams=1 --apertureWidth=24 --apertureHeight=0.6 --backgroundAperture=2 -P /opera/visuals/11AQ14-Jul08/OLAPAa_sp2_Normal.eps -F testApe.dat -S testApe.gnu --minorder=22 --maxorder=61 --gapBetweenBeams=0 --nRowSamples=10 --MinTiltAngle=-5 --MaxTiltAngle=0 --tiltAnglePrecision=0.1 \n\n"
	"  -h, --help  display help message\n"
	"  -v, --verbose,  Turn on message sending\n"
	"  -d, --debug,  Turn on debug messages\n"
	"  -t, --trace,  Turn on trace messages\n"
	"  -o, --outputApertureFile=<APER_FILE>, Output aperture file name\n"
	"  -g, --inputgeom=<GEOM_FILE>, Input geometry file\n"
	"  -p, --inputprof=<PROF_FILE>, Input instrument profile file\n"
	"  -s, --inputorderspacing=<ORDS_FILE>, Input order spacing file\n"
	"  -O, --ordernumber=<INT_VALUE>, Pick order number (default = all)\n"
	"  -R, --pickImageRow=<UNS_VALUE>, Pick row number to plot IP model (default = naxis2/2)\n"  
	"  -L, --nRowSamples=<UNS_VALUE>, Number equally spaced rows to sample IP (default = 1)\n"    
	"  -x, --xbin=<UNS_VALUE>, Number of IP points to bin in x-direction\n" 
	"  -M, --minorder=<INT_VALUE>, Define minimum order number\n"
	"  -X, --maxorder=<INT_VALUE>, Define maximum order number\n"
	"  -N, --numberOfBeams=<UNS_VALUE>, Number of beams to split aperture\n"
	"  -G, --gapBetweenBeams=<FLT_VALUE>, Gap between beams in pixel units\n"   
	"  -W, --apertureWidth=<FLT_VALUE>, Aperture width in pixel units\n"
	"  -H, --apertureHeight=<FLT_VALUE>, Aperture height in pixel units\n"
	"  -B, --backgroundAperture=<FLT_VALUE>, Aperture width for background \n"    
	"  -P, --plotfilename=<EPS_FILE>\n"
	"  -F, --datafilename=<DATA_FILE>\n"
	"  -S, --scriptfilename=<GNUPLOT_FILE>\n" 
	"  -I, --interactive=<BOOL>\n\n";
}

void GenerateExtractionAperturePlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display)
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
	 
    if(!outputPlotEPSFileName.empty()) {
        *fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        *fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        *fgnu << "\nsplot \"" << dataFileName << "\" u 2:3:4 with pm3d" << endl;
        
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
        
        *fgnu << "\nsplot \"" << dataFileName << "\" u 2:3:4 with pm3d" << endl;
        
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

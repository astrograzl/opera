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
#include "libraries/operaIOFormats.h"
#include "libraries/operaCCD.h"						// for MAXORDERS
#include "libraries/ladfit.h"
#include "libraries/operaArgumentHandler.h"
#include "libraries/operaCommonModuleElements.h"

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

void GenerateExtractionAperturePlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display);

void GenerateExtractionApertureTiltPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName1, string dataFileName2, float tiltAngle, float tiltAngleError, bool display);

int main(int argc, char *argv[])
{
	operaArgumentHandler args;
		
    string outputApertureFile;
	string inputgeom;
	string inputprof;
	string inputorderspacing;
    unsigned nRowSamples = 1;
	unsigned pickImageRow = 0;
    unsigned xbin = 10;
    unsigned numberOfBeams = 1; // 1 for star-only;  2 for polar/s+s
    double gapBetweenBeams = 0;
    double apertureWidth = 26.0;
    double apertureHeight = 0.6;
    double backgroundAperture = 2.0;
    bool constantTilt = false;
	int ordernumber = NOTPROVIDED;
	int minorder = NOTPROVIDED;
    int maxorder = NOTPROVIDED;
    string plotfilename;	
	string datafilename;	
	string scriptfilename;
	bool interactive = false;
    string tiltplotfilename;
    string tiltdata1filename;
    string tiltdata2filename;
    string tiltscriptfilename;
    
    args.AddRequiredArgument("outputApertureFile", outputApertureFile, "Output aperture file name");
    args.AddRequiredArgument("inputgeom", inputgeom, "Input geometry file");
    args.AddRequiredArgument("inputprof", inputprof, "Input instrument profile file");
    args.AddOptionalArgument("inputorderspacing", inputorderspacing, "", "Input order spacing file");
    args.AddRequiredArgument("nRowSamples", nRowSamples, "Number of equally spaced rows to sample IP");
    args.AddOptionalArgument("pickImageRow", pickImageRow, 0, "Specific single row to use for IP model (overrides nRowSamples)");
    args.AddRequiredArgument("xbin", xbin, "Number of IP points to bin in x-direction");
    args.AddRequiredArgument("numberOfBeams", numberOfBeams, "Number of beams to split aperture");
    args.AddRequiredArgument("gapBetweenBeams", gapBetweenBeams, "Gap between beams in pixel units");
    args.AddRequiredArgument("apertureWidth", apertureWidth, "Aperture width in pixel units");
    args.AddRequiredArgument("apertureHeight", apertureHeight, "Aperture height in pixel units");
    args.AddRequiredArgument("backgroundAperture", backgroundAperture, "Aperture width for background in pixel units");
    args.AddOptionalArgument("constantTilt", constantTilt, false, "Set median tilt angle to all orders");
    args.AddOrderLimitArguments(ordernumber, minorder, maxorder, NOTPROVIDED);
    args.AddPlotFileArguments(plotfilename, datafilename, scriptfilename, interactive);
    args.AddOptionalArgument("tiltplotfilename", tiltplotfilename, "", "Tilt plot eps file name");
    args.AddOptionalArgument("tiltdata1filename", tiltdata1filename, "", "Tilt first data file name");  
    args.AddOptionalArgument("tiltdata2filename", tiltdata2filename, "", "Tilt second data file name");  
    args.AddOptionalArgument("tiltscriptfilename", tiltscriptfilename, "", "Tilt gnuplot script file name");  
	
	try {
		args.Parse(argc, argv);
		
		// we need a geometry file...
		if (inputgeom.empty()) {
			throw operaException("operaExtractionApertureCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		// we need an instrument profile file...
		if (inputprof.empty()) {
			throw operaException("operaExtractionApertureCalibration: ", operaErrorNoInput, __FILE__, __FUNCTION__, __LINE__);	
		}		
		
		if (args.verbose) {
			cout << "operaExtractionApertureCalibration: outputApertureFile = " << outputApertureFile << endl;            
			cout << "operaExtractionApertureCalibration: inputorderspacing = " << inputorderspacing << endl;             
			cout << "operaExtractionApertureCalibration: inputgeom = " << inputgeom << endl; 
			cout << "operaExtractionApertureCalibration: inputprof = " << inputprof << endl;
            if(ordernumber != NOTPROVIDED) cout << "operaExtractionApertureCalibration: ordernumber = " << ordernumber << endl;
            if(pickImageRow) cout << "operaExtractionApertureCalibration: pickImageRow = " << pickImageRow << endl;       
            else cout << "operaExtractionApertureCalibration: nRowSamples = " << nRowSamples << endl;       
            cout << "operaExtractionApertureCalibration: xbin = " << xbin << endl;
			cout << "operaExtractionApertureCalibration: numberOfBeams = " << numberOfBeams << endl;       
			cout << "operaExtractionApertureCalibration: apertureWidth = " << apertureWidth << endl;            
			cout << "operaExtractionApertureCalibration: apertureHeight = " << apertureHeight << endl;
			cout << "operaExtractionApertureCalibration: backgroundAperture = " << backgroundAperture << endl;
			cout << "operaExtractionApertureCalibration: constantTilt = " << constantTilt << endl;
            if(args.plot) {
                cout << "operaExtractionApertureCalibration: plotfilename = " << plotfilename << endl;
                cout << "operaExtractionApertureCalibration: datafilename = " << datafilename << endl;
                cout << "operaExtractionApertureCalibration: scriptfilename = " << scriptfilename << endl; 
                cout << "operaExtractionApertureCalibration: interactive = " << (interactive ? "YES" : "NO") << endl;
                cout << "operaExtractionApertureCalibration: tiltplotfilename = " << tiltplotfilename << endl;
                cout << "operaExtractionApertureCalibration: tiltdata1filename = " << tiltdata1filename << endl;
                cout << "operaExtractionApertureCalibration: tiltdata2filename = " << tiltdata2filename << endl;
                cout << "operaExtractionApertureCalibration: tiltscriptfilename = " << tiltscriptfilename << endl;
            }
		}
        
        ofstream fdata;
        ofstream ftiltdata1;
        ofstream ftiltdata2;
        if (!datafilename.empty()) fdata.open(datafilename.c_str());
        if (!tiltdata1filename.empty()) ftiltdata1.open(tiltdata1filename.c_str());
        if (!tiltdata2filename.empty()) ftiltdata2.open(tiltdata2filename.c_str());
        
		operaSpectralOrderVector spectralOrders;
		operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputgeom);
        operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputprof);
        if(!inputorderspacing.empty()) operaIOFormats::ReadIntoSpectralOrders(spectralOrders, inputorderspacing);
        
        UpdateOrderLimits(ordernumber, minorder, maxorder, spectralOrders);
        
        /*
         * Set aperture Rectangle aperture
         */   
        float width = float(apertureWidth);
        float height = float(apertureHeight);
        
        float *tiltAngle = new float[MAXORDERS];
        float *tiltAngleError = new float[MAXORDERS];
        float *FluxFraction = new float[MAXORDERS];
        
        unsigned orderIndex=0;
        
        float *xphotocenter  = new float[nRowSamples*MAXORDERS];
        float *yphotocenter = new float[nRowSamples*MAXORDERS];
        float *photosum = new float[nRowSamples*MAXORDERS];
        unsigned nTotalPoints = 0;
        
        // For plot tilt datafile
        float *xphotocenterBin  = new float[nRowSamples*MAXORDERS];
        float *yphotocenterBin = new float[nRowSamples*MAXORDERS];
        unsigned nTotalPointsBin = 0;
        
        int *skippedOrderIndex = new int[MAXORDERS];
        unsigned Nskipped = 0;
        
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
                        if(args.debug){ // for 3D plot
                            cout <<instrumentProfile->getIPixXCoordinate(i) << " "
                            << instrumentProfile->getIPixYCoordinate(j) << " "
                            << (float)instrumentProfile->getipDataFromPolyModel(ymiddle,i,j) << endl;
                        }
                    }
                    //cout << endl; // for 3D plot
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
                    unsigned npoints = 0;
                    
                    float *xcoords = new float[NXPoints];
                    float *ycoords = new float[NXPoints*NYPoints];
                    unsigned nyp = 0;

                    //For tilt plot data file
                    nTotalPointsBin = 0;
                    
                    for (unsigned i=xbin; i<NXPoints-xbin; i++) {
                        
                        float xcoord = instrumentProfile->getIPixXCoordinate(i);
                        xcoords[nlocalpoints] = xcoord;
                        
                        float localphotosum = 0;
                        ycoords[nyp]=0;
                        for (unsigned j=0; j<NYPoints; j++) {
                            float ip = (float)instrumentProfile->getipDataFromPolyModel(ycenter,i,j);
                            if(ip>0 && !isnan(ip)) {
                                float ycoord = instrumentProfile->getIPixYCoordinate(j);
                                ycoords[nyp] += ycoord*ip;
                                yphotocenter_tmp += ycoord*ip;
                                localphotosum += ip;
                                npoints++;
                            }
                        }
                        if(localphotosum) {
                            ycoords[nyp] /= localphotosum;
                            nyp++;
                            xphotocenter_tmp += localphotosum*xcoord;
                            photosum_tmp += localphotosum;
                        }
                        
                        nlocalpoints++;
                        
                        if(nlocalpoints == xbin) {
                            
                            xphotocenter[nTotalPoints] = xphotocenter_tmp/photosum_tmp;
                            yphotocenter[nTotalPoints] = yphotocenter_tmp/photosum_tmp;
                            photosum[nTotalPoints] = photosum_tmp;
                            
                            // For plot tilt data file:
                            xphotocenterBin[nTotalPointsBin] = xphotocenter_tmp/photosum_tmp;
                            yphotocenterBin[nTotalPointsBin] = yphotocenter_tmp/photosum_tmp;
                            nTotalPointsBin++;
                            
                            float xdev = operaArraySigma(nlocalpoints,xcoords);
                            float ydev = operaArraySigma(nyp,ycoords);
                            
                            if(args.debug)
                                cout << order << " " << xphotocenter[nTotalPoints] << " " << yphotocenter[nTotalPoints] << " " <<  xdev << " " << ydev << " " << photosum_tmp/float(npoints) << endl;

                            nTotalPoints++;

                            nyp = 0;
                            xphotocenter_tmp = 0;
                            yphotocenter_tmp = 0;
                            photosum_tmp = 0;
                            nlocalpoints = 0;
                            npoints = 0;
                        }
                    }

                    if (!tiltdata1filename.empty()) { // Variable tilt debug
                        /*
                         * This debug section is useful to produce a plot that indicates
                         * any variation of the tilt angle along the detector rows or
                         * along orders.
                         */
                        float amBin,bmBin,abdevmBin;
                        ladfit(xphotocenterBin,yphotocenterBin,nTotalPointsBin,&amBin,&bmBin,&abdevmBin);
                        float tiltAngleBin = 180*atan(bmBin)/(M_PI);
                        float tiltAngleBinError = abdevmBin/0.674433;
                        ftiltdata1 << order << " " << ycenter << " " << tiltAngleBin << " " << tiltAngleBinError << endl;
                    }

                }
                
                float am,bm,abdevm;
                
                ladfit(xphotocenter,yphotocenter,nTotalPoints,&am,&bm,&abdevm); /* robust linear fit: f(x) =  a + b*x */
                
                if(args.debug) {
                    // This debug section is useful to produce plots containing
                    // individual x and y-centroids in the IP as well as
                    // the model obtained from measurements above
                    for (unsigned i=0; i<NXPoints; i++) {
                        float ipXcoord = instrumentProfile->getIPixXCoordinate(i);
                        float ipsum = 0;
                        float ycentroid = 0;
                        unsigned nnp = 0;
                        for (unsigned j=0; j<NYPoints; j++) {
                            float ipYcoord = instrumentProfile->getIPixYCoordinate(j);
                            float ip = (float)instrumentProfile->getipDataFromPolyModel(ymiddle,i,j);
                            if(ip>0) {
                                ipsum += ip;
                                ycentroid += ipYcoord*ip;
                                nnp++;
                            }
                        }
                        float ipAvg = ipsum/(float)nnp;
                        ycentroid /= ipsum;
                        float ipYmodel = am + bm*ipXcoord;
                        cout << ipXcoord << " " << ipYmodel << " " << ycentroid << " " << ipAvg << endl;
                    }
                } // end debug
                
                if (!isnan(bm) && bm) {
        
                    tiltAngle[orderIndex]=180*atan(bm)/(M_PI);
                    tiltAngleError[orderIndex] = abdevm/0.674433;
                    
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
                    
                    if(!tiltdata2filename.empty()) {
                        ftiltdata2 << order << " " << tiltAngle[orderIndex] << " " << tiltAngleError[orderIndex]  << " " << FluxFraction[orderIndex] << endl;
                    }
                    
                    if(args.verbose)
                        cout << "operaExtractionApertureCalibration: # order = " << order << " tilt = " << tiltAngle[orderIndex] << " +/- " << 180*(abdevm/0.674433)/(M_PI) << " FluxFraction=" << FluxFraction[orderIndex] << endl;
                    
                    spectralOrder->setTiltInDegrees(tiltAngle[orderIndex],tiltAngleError[orderIndex]);
                    
                    orderIndex++;
                } else { 
                    skippedOrderIndex[Nskipped++] = order;
                    cout << "operaExtractionApertureCalibration: WARNING can't calculate tilt for order " << order << endl;
                }
            }
        }
        
        float tilt = operaArrayMedian(orderIndex,tiltAngle);
        float tilterror = operaArrayMedianSigma(orderIndex,tiltAngle,tilt);
        
        if(args.verbose) cout << "operaExtractionApertureCalibration: # Final tilt = " << tilt << " +/- " << tilterror << endl;
        
        // Below it assigns the median final tilt to orders that have been skipped
        for (unsigned i=0; i<Nskipped; i++) {
            operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(skippedOrderIndex[i]);
            spectralOrder->setTiltInDegrees(tilt,tilterror);
        }
        
        // Below it assign the median final tilt to all orders
        if(constantTilt) {
            for (int order=minorder; order<=maxorder; order++) {
                operaSpectralOrder *spectralOrder = spectralOrders.GetSpectralOrder(order);
                spectralOrder->setTiltInDegrees(tilt,tilterror);
            }
        }
        
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
                
                if(args.verbose) cout << "operaExtractionApertureCalibration: Calculating the aperture of order " << order << endl;
                
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
                if(pickImageRow) ycenter = (float)pickImageRow + 0.5;
                else ycenter = (Geometry->getYmax() + Geometry->getYmin())/2 + 0.5;
                
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
                
                float beamWidth = width/(float)numberOfBeams - float(gapBetweenBeams);
                float widthX = width*cos(tiltAngle[orderIndex]*M_PI/180);
                float aperXsize = widthX/(float)numberOfBeams;
                for(unsigned k=0;k<numberOfBeams;k++) {
                    beamFluxFraction[k]=0;
                    aperture[k] = NULL;
                }
                
                for(unsigned k=0;k<numberOfBeams;k++) {
                    if(args.verbose) cout << "operaExtractionApertureCalibration: beam " << k << endl;
                    
                    outBoundPoints[k] = 0;
                    
                    float beamXcenter = MidPoint.getXcoord() - (float)(numberOfBeams-1)*aperXsize/2 + aperXsize*(float)k;
                    float beamYcenter = extractionLine.getYcoord(beamXcenter);            
                    
                    beamMidPoint.setPoint(beamXcenter,beamYcenter);
					// DT May 28 2014 this is set in to the aperture and must be a pointer and must not be deleted
                    Line *beamExtractionLine = new Line((float)tan(tiltAngle[orderIndex]*M_PI/180),height,beamWidth,beamMidPoint);
                    
                    // NOTE: This should not be deleted as it is set into the spectralorder a a pointer
                    aperture[k] = new operaExtractionAperture(beamExtractionLine,instrumentProfile,ycenter);                          
                    
                    PixelSet *aperturePixels = aperture[k]->getSubpixels();   
                    
                    for(unsigned i=0; i<aperturePixels->getNPixels(); i++){
                        if(aperturePixels->getiIndex(i) >= 0 && aperturePixels->getiIndex(i) < (int)NXPoints && aperturePixels->getjIndex(i) >= 0 && aperturePixels->getjIndex(i) < (int)NYPoints) {
							ipImage[(unsigned)aperturePixels->getjIndex(i)][(unsigned)aperturePixels->getiIndex(i)] = 0;
							beamFluxFraction[k] += (float)instrumentProfile->getipDataFromPolyModel(ycenter,(unsigned)aperturePixels->getiIndex(i),(unsigned)aperturePixels->getjIndex(i))/normalizationFactor;
                        } else {
                            outBoundPoints[k]++;
                        }
                    }  
                    
                    if(args.verbose)
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
                
                
                if(args.verbose) cout << "operaExtractionApertureCalibration: Setting up background aperture " << endl;
                /*
                 * set up background aperture:
                 */
                float leftbkgXcenter = MidPoint.getXcoord() - (float)numberOfBeams*aperXsize/2 - float(backgroundAperture)/2;
                float leftbkgYcenter = extractionLine.getYcoord(leftbkgXcenter);            
                
                operaPoint leftBackgroundMidPoint(leftbkgXcenter,leftbkgYcenter);
                Line leftBackgroundExtractionLine(tan(tiltAngle[orderIndex]*M_PI/180),height,float(backgroundAperture),leftBackgroundMidPoint);
                
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
                if(args.verbose)
                    cout << "operaExtractionApertureCalibration: # order = " << order << " Left Background Flux Fraction = "<< leftBackgroundFluxFraction << " # out-of-bound points = " << leftBackgroundoutBoundPoints << endl;
                
                leftBackgroundAperture->setFluxFraction(leftBackgroundFluxFraction);
                spectralOrder->setBackgroundApertures(0,leftBackgroundAperture);
                
                float rightbkgXcenter = MidPoint.getXcoord() + (float)numberOfBeams*aperXsize/2 + float(backgroundAperture)/2;
                float rightbkgYcenter = extractionLine.getYcoord(rightbkgXcenter);            
                
                operaPoint rightBackgroundMidPoint(rightbkgXcenter,rightbkgYcenter);
                Line rightBackgroundExtractionLine(tan(tiltAngle[orderIndex]*M_PI/180),height,float(backgroundAperture),rightBackgroundMidPoint);
                
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
                
                if(args.verbose)
                    cout << "operaExtractionApertureCalibration: # order = " << order << " Right Background Flux Fraction = "<< rightBackgroundFluxFraction << " # out-of-bound points = " << rightBackgroundoutBoundPoints << endl << endl;            
                
                for (unsigned j=0; j<NYPoints; j++) {	 
                    for (unsigned i=0; i<NXPoints; i++) {
                        if (fdata.is_open()) {
                            fdata << order << " " 
                            << instrumentProfile->getIPixXCoordinate(i) + shiftX << " " 
                            << instrumentProfile->getIPixYCoordinate(j) + shiftY << " " 
                            << ipImage[j][i]/normalizationFactor << endl;
                        } else if (args.debug) {
                            cout << order << " " 
                            << instrumentProfile->getIPixXCoordinate(i) + shiftX << " " 
                            << instrumentProfile->getIPixYCoordinate(j) + shiftY << " " 
                            << ipImage[j][i]/normalizationFactor << endl;
                        }
                        
                    }
                    if (fdata.is_open()) fdata << endl;
                    else if (args.debug) cout << endl;
                }
                if (fdata.is_open()) fdata << endl;
                else if (args.debug) cout << endl;
                orderIndex++;
                deleteCMatrix(ipImage);
            } else if(args.verbose) {
				cout << "operaExtractionApertureCalibration:aperture of order " << order << " skipped." << endl;
            }
        }        
        
        // output aperture to file...
        operaIOFormats::WriteFromSpectralOrders(spectralOrders, outputApertureFile,Aperture);
        
        if (ftiltdata1.is_open()) ftiltdata1.close();
        if (ftiltdata2.is_open()) ftiltdata2.close();

        if (!tiltdata1filename.empty() && !tiltdata2filename.empty() && !scriptfilename.empty()) {
            GenerateExtractionApertureTiltPlot(tiltscriptfilename, tiltplotfilename, tiltdata1filename, tiltdata2filename, tilt, tilterror, interactive);
        }
 
        if (fdata.is_open()) {
            fdata.close();
            if (!scriptfilename.empty()) {
                GenerateExtractionAperturePlot(scriptfilename,plotfilename,datafilename,interactive);
            }
        }
        
        delete[] xphotocenter;
        delete[] yphotocenter;
        delete[] photosum;
        
        delete[] tiltAngle;
        delete[] FluxFraction;
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

void GenerateExtractionAperturePlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, bool display)
{
	if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());
    
    fgnu << "reset" << endl;
    fgnu << "unset key" << endl;
    fgnu << "set view 0,0" << endl;
    
    fgnu << "set palette color" << endl;
    fgnu << "set palette gamma 2.5" << endl;
    fgnu << "set pm3d map" << endl;
    fgnu << "unset ztics" << endl;
    fgnu << "set cblabel \"flux fraction\"" << endl;
    fgnu << endl;
	 
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nsplot \"" << dataFileName << "\" u 2:3:4 with pm3d" << endl;
        
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
        fgnu << "\nsplot \"" << dataFileName << "\" u 2:3:4 with pm3d" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

void GenerateExtractionApertureTiltPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName1, string dataFileName2, float tiltAngle, float tiltAngleError, bool display)
{
    if (gnuScriptFileName.empty()) exit(EXIT_FAILURE);
	remove(gnuScriptFileName.c_str()); // delete any existing file with the same name
	ofstream fgnu(gnuScriptFileName.c_str());

    fgnu << "reset" << endl;

    fgnu << "set xlabel \"order number\"" << endl;
    fgnu << "set ylabel \"tilt angle (deg)\"" << endl;
    
    fgnu << "set yrange[" << tiltAngle-tiltAngleError*6 << ":" << tiltAngle+tiltAngleError*6 << "]" << endl;
    
    fgnu << "tilt(x) = " << tiltAngle << endl;
    fgnu << "maxtilt(x) = " << tiltAngle+tiltAngleError << endl;
    fgnu << "mintilt(x) = " << tiltAngle-tiltAngleError << endl;
    
    fgnu << endl;
    
    if(!outputPlotEPSFileName.empty()) {
        fgnu << "\nset terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "set output \"" << outputPlotEPSFileName << "\"" << endl;
        
        fgnu << "\nplot \"" << dataFileName1 << "\" u 1:3:4 t \"individual samples\" w yerr pt 6 ps 1 lw 1 "
        << ",\"" << dataFileName2 << "\" u 1:2 notitle w l lw 2 lt 3"
        << ",\"" << dataFileName2 << "\" u 1:2:3 t \"order median tilt\" w yerr pt 7 lt 3 ps 2 lw 2"
        << ",tilt(x) t \"final median tilt\" w l lt -1 lw 2.5, maxtilt(x) notitle w l lt -1, mintilt(x) notitle w l lt -1" << endl;

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
        fgnu << "\nplot \"" << dataFileName1 << "\" u 1:3:4 t \"individual samples\" w yerr pt 6 ps 1 lw 1 "
        << ",\"" << dataFileName2 << "\" u 1:2 notitle w l lw 2 lt 3"
        << ",\"" << dataFileName2 << "\" u 1:2:3 t \"order median tilt\" w yerr pt 7 lt 3 ps 2 lw 2"
        << ",tilt(x) t \"final median tilt\" w l lt -1 lw 2.5, maxtilt(x) notitle w l lt -1, mintilt(x) notitle w l lt -1" << endl;
        
        fgnu << "\n#set terminal postscript enhanced color solid lw 1.5 \"Helvetica\" 14" << endl;
        fgnu << "#set output \"outputPlotEPSFileName.eps\"" << endl;
        fgnu << "#replot" << endl;
        fgnu << "#set terminal x11" << endl;
        fgnu << "#set output" << endl;
    }
    
    fgnu.close();
    
    if (display) systemf("gnuplot -persist %s",gnuScriptFileName.c_str());
    else if (!outputPlotEPSFileName.empty()) systemf("gnuplot %s",gnuScriptFileName.c_str());
}

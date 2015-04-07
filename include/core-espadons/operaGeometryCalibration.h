#ifndef OPERAGEOMETRYCALIBRATION_H
#define OPERAGEOMETRYCALIBRATION_H
/*******************************************************************
 ****                  MODULE FOR OPERA v1.0                     ***
 *******************************************************************
 Module name: operaGeometryCalibration
 Version: 1.0
 Description: Perform a geometric calibration, finding the order polynomials. 
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

/*! \brief  Perform a geometric calibration, finding the order polynomials. */
/*! \file operaGeometryCalibration.h */
/*! \ingroup core */

#include <string>
#include "libraries/Polynomial.h"
#include "libraries/operaFITSImage.h"

/* prototypes */

static void printUsageSyntax(char * prgname);

void GenerateGeometryPlot(string gnuScriptFileName, string outputPlotEPSFileName, string dataFileName, unsigned npolynomials, Polynomial *polynomials[], bool display, unsigned col0, unsigned colf, unsigned row0, unsigned rowf, int AbsOrdNumber[]);

float calculateCenterOfSymmetry(unsigned np, float *ipx, float *ipfunc, float *iperr, unsigned totalNumberOfSlices, bool applyXcenterCorrection);

float calculatePhotoCenter(unsigned np, float *ipx, float *ipfunc, float *iperr, bool applyXcenterCorrection);

unsigned getRowBinnedData(operaFITSImage *flat,unsigned x1,unsigned x2,unsigned nx,unsigned y1,unsigned y2,unsigned ny,float *fx,float *fy,float *yout, bool FFTfilter);

unsigned geometryDetectOrders(unsigned np,float *fx,float *fy,unsigned uslit,float *ipfunc, unsigned binsize, float noise,float gain,float *xmean,float *ymean,float *xmeanerr,int detectionMethod, bool witherrors, bool graces);


#endif

#ifndef OPERALIBCOMMON_H
#define OPERALIBCOMMON_H

/*******************************************************************
****               		OPERA PIPELINE v1.0                     ****
********************************************************************
Library name: operaLibCommon
Version: 1.0
Author(s): CFHT OPERA team
Affiliation: Canada France Hawaii Telescope 
Location: Hawaii USA
Date: Aug/2011

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

/*! 
 * operaLibCommon
 * \author Doug Teeple / Eder Martioli
 * \brief Common C language functions and constants.
 * \file operaLibCommon.h
 * \ingroup libraries
 */

#ifdef __cplusplus
extern "C" {
#endif
	
#include <stdlib.h>
#include <string.h>

	enum operaSpectralOrder_t {
		None = 0, 
		RawSpectrum, 				// outputs of operaStarOnly and operaStarPlusSky...
		StandardSpectrum, 			// uncalibrated / standard vectors of x,y coordinates and flux
		OptimalSpectrum, 			// uncalibrated / optimal vectors of x,y coordinates and flux
		OperaOptimalSpectrum, 		// uncalibrated / opera optimal vectors of x,y coordinates and flux
		RawBeamSpectrum, 			// output of operaExtraction distance-based
		StandardBeamSpectrum, 		// output of operaExtraction distance-based
		OptimalBeamSpectrum, 		// output of operaExtraction distance-based
		OperaOptimalBeamSpectrum, 	// output of operaExtraction distance-based
		CalibratedRawSpectrum, 		// Calibrated to wavelength.
		CalibratedStandardSpectrum, // Calibrated to wavelength
		CalibratedOptimalSpectrum, 	// Calibrated to wavelength
		CalibratedOperaOptimalSpectrum, // Calibrated to wavelength
		CalibratedRawBeamSpectrum, 	// Calibrated to wavelength
		CalibratedStandardBeamSpectrum, // Calibrated to wavelength
		CalibratedOptimalBeamSpectrum, 	// Calibrated to wavelength
		CalibratedOperaOptimalBeamSpectrum, // Calibrated to wavelength
		CalibratedExtendedBeamSpectrum, // multiple wl and flux
		ExtendedPolarimetry,		 // multiple wl and flux
		LibreEspritSpectrum,		// for backwards compatibility calibrated spectrum
		LibreEspritsp1Spectrum,		// for backwards compatibility calibrated spectrum
		LibreEspritsp2Spectrum,		// for backwards compatibility calibrated spectrum
		LibreEspritpolSpectrum,		// for backwards compatibility calibrated spectrum
		LibreEspritpolarimetry,		// for backwards compatibility calibrated spectrum
		LibreEspritSNR,				// for backwards compatibility calibrated SNR
		ReferenceSpectrum,			// solar reference spectra
		Geom,						// geometry vectors / polynomial
		Wave, 						// wavelength vectors / polynomial
		Spec, 						// spectral element vectors / polynomial
		Prof, 						// instrument profile vectors / polynomial
		Disp, 						// Dispersion polynomial
		SNR,						// SNR for each order
		Polarimetry,				// polarimetry
		Orderspacing,				// the geometry order spacing polynomial
		Lines,						// spectral lines
		GainNoise,					// gain noise calibration data
		Aperture,					// aperture data
		Fcal,                       // flux calibration data
		RVel,                       // Radial Velocity Correction
		PRVel,                      // Polarimetry Radial Velocity Correction
		Tell,                       // Telluric wave Correction
		PTell,                      // Polarimetry Autowave Correction
		EETC,                       // Espadons Exposure Time Calculator file
		CSV,                        // Comma-Separated Values
		LastSpectralOrder_t			// the end
	};
	
    enum InstrumentEnvironment_t {
		InstrumentEnvironmentUnkownFormat = 0,
        ObsCond,                    // Observing conditions format
		SkyObj,                     // Object in the sky format
        spectrograph,               // Spectrograph configuration format
		telescope                  // Telescope configuration format
    };
    
#ifndef FOURSIDES
#define FOURSIDES 4
#endif
#ifndef MAXNPOLYGONSIDES
#define MAXNPOLYGONSIDES 300
#endif
#ifndef MAXPOLYNOMIAL
#define MAXPOLYNOMIAL 12
#endif
#ifndef MAXGAUSSCOEFFS
#define MAXGAUSSCOEFFS 1000
#endif
#ifndef MAXESPADONSY
#define MAXESPADONSY 4640
#endif
#ifndef MAXESPADONSX
#define MAXESPADONSX 2080
#endif
#ifndef MAXNUMBEROFBEAMS    
#define MAXNUMBEROFBEAMS 20    
#endif    
#ifndef LEFTANDRIGHT     
#define LEFTANDRIGHT 2   
#endif   
#ifndef MAXPOINTSINSIMULATEDSPECTRUM     
#define MAXPOINTSINSIMULATEDSPECTRUM 100000  
#endif      
#ifndef EPS
#define EPS 1.0e-6
#endif
#ifndef DBL_EPSILON
#define DBL_EPSILON 1.11e-16
#endif
#ifndef SWAP
#define SWAP(a,b) { float temp=(a); (a)=(b); (b)=temp; }
#endif
#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif
#ifndef MAX_ESPADONS_VECTOR_LENGTH
#define MAX_ESPADONS_VECTOR_LENGTH 51200
#endif
#ifndef MAXREFWAVELENGTHSPERORDER
#define MAXREFWAVELENGTHSPERORDER 5000
#endif
#ifndef BIG
#define BIG 1.0e30
#endif
#ifndef true
#define true 1
#endif
#ifndef false
#define false 0
#endif
#ifndef MIN
#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x)   (x)*(x)
#endif
#ifndef PI 
#define PI 3.141592654
#endif
#ifndef TWOPI 
#define  TWOPI             6.28318530717959
#endif
#ifndef PI_OVER_2 
#define  PI_OVER_2         1.57079632679490  /* From Abramowitz & Stegun */   
#endif
#ifndef ARCSEC_IN_RADIAN 
#define  ARCSEC_IN_RADIAN  206264.8062471
#endif
#ifndef DEG_IN_RADIAN 
#define  DEG_IN_RADIAN     57.2957795130823
#endif
#ifndef HRS_IN_RADIAN 
#define  HRS_IN_RADIAN     3.819718634205
#endif
#ifndef SPEED_OF_LIGHT_M 
#define  SPEED_OF_LIGHT_M     	299792458   /* m/s */
#endif
#ifndef SPEED_OF_LIGHT_KMS
#define  SPEED_OF_LIGHT_KMS     	299792.458   /* km/s */
#endif
#ifndef BOLTZMANN_CONSTANT 
#define  BOLTZMANN_CONSTANT     1.38065e-23 /* J/K */
#endif
#ifndef PLANCK_CONSTANT 
#define  PLANCK_CONSTANT     	6.626069e-34    /* J*s */
#endif
#ifndef VEGA_PHOTONFLUX_AT_548NM
#define  VEGA_PHOTONFLUX_AT_548NM     	9.876e16    /* Photon flux for Vega @ 548nm: 9.876e16 ph/s/m2/m */
#endif
#ifndef VEGA_ENERGYFLUX_AT_548NM
#define  VEGA_ENERGYFLUX_AT_548NM     	0.0358   /* Energy flux for Vega @ 548nm: 3650 Jy +- 2% = 1e-26 W/m^2/Hz = 3.58e-8 W/m^2/microns = 0.0358 W/m^2/m */
#endif
#ifndef VEGA_EFF_TEMPERATURE_K
#define  VEGA_EFF_TEMPERATURE_K  9602
#endif
#ifndef SUN_EFF_TEMPERATURE_K
#define  SUN_EFF_TEMPERATURE_K  5780   
#endif
#ifndef VEGA_VBAND_MAGNITUDE
#define  VEGA_VBAND_MAGNITUDE   0.03
#endif

    
/*!
 * startsWith
 * \author Doug Teeple
 * \brief returns > 0 if astring startswith substring.
 * \arg aString - the string to search
 * \arg substring - the snippet
 * \return EXIT_STATUS
 */
STATIC inline unsigned startsWith(const char *aString, const char *substring) {
	unsigned i;
	if (aString == NULL || substring == NULL)
		return 0;
	else if (strlen(aString) < strlen(substring))
		return 0;
	else
		for (i=0; i < strlen(substring); i++)
			if (aString[i] != substring[i])
				return 0;
	return i;
}
/* 
 * float matrix
 */

typedef struct CMatrixStruct {
	unsigned rows;
	unsigned cols;
	unsigned originrow;
	unsigned origincol;
	unsigned originrow2;	// for ||gram
	unsigned origincol2;	// for ||gram
	unsigned pixelationx;
	unsigned pixelationy;
	float **matrix;
	float *data;
} CMatrixStruct_t;

typedef  float** CMatrix;

/*! 
 * \brief newCMatrix(unsigned cols, unsigned rows)
 * \brief creates a CMatrix, note that the float ** is returned
 * \brief for matrix indexing. A guard pointer of NULL is stored
 * \brief so the base address can be found (one past the gaurd).
 * \param cols
 * \param rows
 * \return CMatrix
 */
STATIC inline CMatrix newCMatrix(unsigned Cols, unsigned Rows) {
	CMatrixStruct_t *cmatrix = (CMatrixStruct_t *)malloc(sizeof(CMatrixStruct_t));
	cmatrix->matrix = (float **)malloc(sizeof(float *)*(Rows+2));
	if (!cmatrix->matrix) {
		return NULL;
	}
	cmatrix->originrow = 0;
	cmatrix->origincol = 0;
	cmatrix->originrow2 = 0;
	cmatrix->origincol2 = 0;
	cmatrix->pixelationx = 1;
	cmatrix->pixelationy = 1;
	cmatrix->rows = Rows;
	cmatrix->cols = Cols;
	cmatrix->data = (float *)malloc(Cols*Rows*sizeof(float));
	if (!cmatrix->data) {
		return NULL;
	}
	memset(cmatrix->data, 0, Cols*Rows*sizeof(float));
	unsigned row;
	CMatrix rv = (CMatrix)cmatrix->matrix;
	for (row=0; row<Rows; row++) {
		rv[row] = cmatrix->data + Cols*row;
	}
	rv[row++] = NULL;			// the guard entry
	rv[row] = (float *)cmatrix;	// store the struct address so it an be deleted in entirety
	return cmatrix->matrix;		// return the float **
}

STATIC inline CMatrixStruct_t *getCMatrixBase(CMatrix cmatrix) {
	unsigned row = 0;
	if (cmatrix) {
		// find the guard row
		while (cmatrix[row]) {
			row++;
		}
		// one past is the base address
		row++;
		return (CMatrixStruct_t *)cmatrix[row];		
	}
	return NULL;
}

STATIC inline float *getCMatrixData(CMatrix cmatrix) {
	if (cmatrix) {
		return (float *)getCMatrixBase(cmatrix)->data;		
	}
	return NULL;
}
	
STATIC inline unsigned getCMatrixPixelationx(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->pixelationx;
}

STATIC inline void setCMatrixPixelationx(CMatrix cmatrix, unsigned pixelation) {
	getCMatrixBase(cmatrix)->pixelationx = pixelation;
}

STATIC inline unsigned getCMatrixPixelationy(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->pixelationy;
}

STATIC inline void setCMatrixPixelationy(CMatrix cmatrix, unsigned pixelation) {
	getCMatrixBase(cmatrix)->pixelationy = pixelation;
}

STATIC inline unsigned getCMatrixCols(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->cols;
}

STATIC inline unsigned getCMatrixRows(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->rows;
}

STATIC inline unsigned getCMatrixOriginCol(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->origincol;
}

STATIC inline void setCMatrixOriginCol2(CMatrix cmatrix, unsigned Col) {
	getCMatrixBase(cmatrix)->origincol2 = Col;
}

STATIC inline unsigned getCMatrixOriginCol2(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->origincol2;
}

STATIC inline void setCMatrixOriginCol(CMatrix cmatrix, unsigned Col) {
	getCMatrixBase(cmatrix)->origincol = Col;
}

STATIC inline unsigned getCMatrixOriginRow(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->originrow;
}

STATIC inline void setCMatrixOriginRow2(CMatrix cmatrix, unsigned Row) {
	getCMatrixBase(cmatrix)->originrow2 = Row;
}

STATIC inline unsigned getCMatrixOriginRow2(CMatrix cmatrix) {
	return getCMatrixBase(cmatrix)->originrow2;
}

STATIC inline void setCMatrixOriginRow(CMatrix cmatrix, unsigned Row) {
	getCMatrixBase(cmatrix)->originrow = Row;
}

STATIC inline void deleteCMatrix(CMatrix cmatrix) {
	if (!cmatrix) {
		return;
	}
	CMatrixStruct_t *cm = getCMatrixBase(cmatrix);
	if (cm) {
		if(cm->matrix)
			free(cm->matrix);
		if(cm->data)
			free(cm->data);
		free(cm);
	}
}

/*
 * The Cube
 */
typedef  float*** CCube;

/*! 
 * \brief STATIC inline CCube newCCube(unsigned Cols, unsigned Rows, unsigned Depth)
 * \brief creates a CCube, note that the float *** is returned
 * \brief for matrix indexing. A guard pointer of NULL is stored
 * \param cols
 * \param rows
 * \param Depth
 * \return CMatrix*
 */
STATIC inline CCube newCCube(unsigned Cols, unsigned Rows, unsigned Depth) {
	CMatrix *depth = (CMatrix *)malloc(sizeof(CMatrix *)*(Depth+4));
	if (!depth) {
		return NULL;
	}
	CMatrix *dp = depth;
	for (unsigned d=0; d<Depth; d++) {
		*dp = newCMatrix(Cols, Rows);
		if (!dp) {
			return NULL;
		}
		dp++;
	}
	*dp++ = NULL;				// the guard entry
#ifdef __LP64__
	*dp++ = (CMatrix)(unsigned long)Depth;		// keep the depth for debugging
	*dp++ = (CMatrix)(unsigned long)Cols;		// keep the cols for debugging
	*dp++ = (CMatrix)(unsigned long)Rows;		// keep the rows for debugging
#else
	*dp++ = (CMatrix)Depth;		// keep the depth for debugging
	*dp++ = (CMatrix)Cols;		// keep the cols for debugging
	*dp++ = (CMatrix)Rows;		// keep the rows for debugging
#endif
	return depth;				// return the float ***
}

STATIC inline void deleteCCube(CCube cCube) {
	if (!cCube) {
		return;
	}
	CMatrix *dp = cCube;
	while (*dp) {
		deleteCMatrix(*dp++);
	}
	free(cCube);
}

/*
 * Polynomials
 */

typedef struct PolynomialCoeffs {
	int orderofPolynomial;
	float p[MAXPOLYNOMIAL];
	float e[MAXPOLYNOMIAL];
	float polychisqr;
} PolynomialCoeffs_t;

typedef struct doublePolynomialCoeffs {
	int orderofPolynomial;
	double p[MAXPOLYNOMIAL];
	double e[MAXPOLYNOMIAL];
	double polychisqr;
} doublePolynomialCoeffs_t;

STATIC inline void PolynomialCoeffsToFloat(PolynomialCoeffs_t *f, doublePolynomialCoeffs_t *d) {
	f->orderofPolynomial = d->orderofPolynomial;
	for (unsigned i=0; i<MAXPOLYNOMIAL; i++) {
		f->p[i] = (float)d->p[i];
		f->e[i] = (float)d->e[i];
	}
	f->polychisqr = (float)d->polychisqr;
}
/*
 * The Polynomial Matrix
 */

typedef  PolynomialCoeffs_t*** PolynomialMatrix;

typedef struct PolynomialMatrixStruct {
	unsigned rows;
	unsigned cols;
	PolynomialCoeffs_t ***matrix;		// the index to reach a row of data
	PolynomialCoeffs_t **data;			// the pointers to polynomials
} PolynomialMatrixStruct_t;

/*! 
 * \brief newPolynomialMatrix(unsigned cols, unsigned rows)
 * \brief creates a CMatrix, note that the float ** is returned
 * \brief for matrix indexing. A guard pointer of NULL is stored
 * \brief so the base address can ebe found (one past the gaurd).
 * \param cols
 * \param rows
 * \return CMatrix
 */
STATIC inline PolynomialMatrix newPolynomialMatrix(unsigned Cols, unsigned Rows) {
	PolynomialMatrixStruct_t *pmatrix = (PolynomialMatrixStruct_t *)malloc(sizeof(PolynomialMatrixStruct_t));
	pmatrix->matrix = (PolynomialCoeffs_t ***)malloc(sizeof(PolynomialCoeffs_t ***)*(Rows+2));
	pmatrix->rows = Rows;
	pmatrix->cols = Cols;
	pmatrix->data = (PolynomialCoeffs_t **)malloc(Cols*Rows*sizeof(PolynomialCoeffs_t **));
	PolynomialCoeffs_t **pp = (PolynomialCoeffs_t **)pmatrix->data;
	for (unsigned y=0; y<Rows; y++) {
		for (unsigned x=0; x<Cols; x++) {
			*pp = (PolynomialCoeffs_t *)malloc(sizeof(PolynomialCoeffs_t));
			if (!pp) {
				return NULL;
			}
			memset(*pp, 0, sizeof(PolynomialCoeffs_t));
			pp++;
		}
	}
	unsigned row;
	PolynomialMatrix rv = pmatrix->matrix;
	for (row=0; row<Rows; row++) {
		rv[row] = pmatrix->data + Cols*row;
	}
	rv[row++] = NULL;			// the guard entry
	rv[row] = (PolynomialCoeffs_t **)pmatrix;	// store the struct address so it an be deleted in entirety
	return pmatrix->matrix;		// return the PolynomialCoeffs_t **
}

STATIC inline PolynomialMatrixStruct_t *getPolynomialMatrixBase(PolynomialMatrix pmatrix) {
	unsigned row = 0;
	// find the guard row
	while (pmatrix[row]) {
		row++;
	}
	// one past is the base address
	row++;
	return (PolynomialMatrixStruct_t *)pmatrix[row];
}


STATIC inline unsigned getPolynomialMatrixCols(PolynomialMatrix pmatrix) {
	return getPolynomialMatrixBase(pmatrix)->cols;
}

STATIC inline unsigned getPolynomialMatrixRows(PolynomialMatrix pmatrix) {
	return getPolynomialMatrixBase(pmatrix)->rows;
}

STATIC inline void deletePolynomialMatrix(PolynomialMatrix pmatrix) {
	PolynomialMatrixStruct_t *cm = getPolynomialMatrixBase(pmatrix);
	if (cm) {
		if (cm->matrix)
			free(cm->matrix);
		if (cm->data)
			free(cm->data);
		free(cm);
	}
}

typedef struct doubleValueVector {
	size_t vectorlength;
	size_t maxlength;
	double *xs;
	double *ys;
	double *values;
	double *errors;
} doubleValueVector_t;

typedef struct WavelengthFluxVector {
	size_t vectorlength;
	size_t maxlength;
	double *ws;
	double *values;
	double *errors;
} WavelengthFluxVector_t;

typedef struct doubleValue {
        double value;
        double error;
} doubleValue_t;  

typedef struct skycoord {
	int ra_h;
	int ra_m;
	double ra_s;
	int dec_d;
	int dec_m;
	double dec_s;
} skycoord_t;    
    
typedef struct geocoord {
	int latitude_d;
	int latitude_m;
	double latitude_s;
	int longitude_d;
	int longitude_m;
	double longitude_s;    
} geocoord_t;
    
#ifdef __cplusplus
}
#endif

#endif


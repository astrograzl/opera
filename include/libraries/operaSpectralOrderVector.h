#ifndef OPERASPECTRALORDERVECTOR_H
#define OPERASPECTRALORDERVECTOR_H
/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaSpectralOrderVector
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

#include "operaError.h"
#include "libraries/operaSpectralElements.h"		// for operaSpectralElements
#include "libraries/operaInstrumentProfile.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/GainBiasNoise.h"
#include "libraries/Polynomial.h"	
#include "libraries/LaurentPolynomial.h"	
#include "libraries/operaWavelength.h" // for MAXORDEROFWAVELENGTHPOLYNOMIAL
#include "libraries/operaEspadonsImage.h" // for instrumentmode_t

#define NEWLINES_BETWEEN_ORDERS false

using namespace std;

/*! 
 * \sa class operaSpectralOrderVector
 * \brief Stores an ordered vector of SpectralOrder class instances. 
 * \return none
 * \file operaSpectralOrderVector.h
 * \ingroup libraries
 */
class operaSpectralOrderVector {
	
private:
	operaSpectralOrder **vector;	// a vector of spectral order pointers
	Polynomial *orderSpacingPolynomial; // captures the spacing between orders
	
	string object;					// for Libre-Esprit output
    
    unsigned numberOfDispersionPolynomials;
	LaurentPolynomial *dispersionPolynomial[MAXORDEROFWAVELENGTHPOLYNOMIAL]; // captures the dispersion polynomials
	
    unsigned length;				// total length of vector
	unsigned minorder;				// lowest order number with a value
	unsigned maxorder;				// highest order number with a value
	unsigned sequence;				// for master fluxcalibrations to store sequence #
	instrumentmode_t instrumentmode;	// instrumentmode
	unsigned count;					// Note that this count is the count of actual orders, 
									// not the length of the vector (which will be MAXORDERS)
	GainBiasNoise *gainBiasNoise;	// flat stats
	double BarycentricRadialVelocityCorrection;

public:
	
	operaSpectralOrder **first;		// pointer to the first spectral order with a value
	operaSpectralOrder **last;		// pointer to the last spectralorder with a value
	/*
	 * Constructors
	 */
	
	/*! 
	 * \sa class operaSpectralOrderVector();
	 * \brief Base constructor.
	 * \return void
	 */
	operaSpectralOrderVector();
	/*! 
	 * \sa class operaSpectralOrderVector(string Filename);
	 * \details Base constructor, read a spectral order vector from a filename
	 * \details Filename can contain either a .geom file or a .wcal file as given by the format
	 * \details both of which are vectors of polynomial coefficients.
	 * \param Filename - string Filename to save to
	 * \return void
	 */
	operaSpectralOrderVector(string Filename);
	/*! 
	 * \sa class operaSpectralOrderVector(unsigned length, operaSpectralOrder_t OrderType, unsigned maxdatapoints, unsigned maxValues, unsigned nElements);
	 * \brief Create a SpectralOrderVector of a given type.
	 * \return void
	 */
	operaSpectralOrderVector(unsigned length, unsigned maxdatapoints, unsigned maxValues, unsigned nElements);
	/*
	 * Destructor
	 */
	~operaSpectralOrderVector();

	/*
	 * Methods
	 */
	
	/*!
	 * \sa method unsigned getnumberOfDispersionPolynomials(void);
	 * \brief returns the number of dispersion polynomials.
	 * \return unsigned - numberOfDispersionPolynomials.
	 */
	unsigned getnumberOfDispersionPolynomials(void);
	
	/*!
	 * \sa method void setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials);
	 * \brief sets the number of dispersion polynomials.
	 * \return none.
	 */
	void setnumberOfDispersionPolynomials(unsigned NumberOfDispersionPolynomials);
    
	/*! 
	 * \sa method void freeoperaSpectralOrderVector();
	 * \brief frees the vector contents and the vector.
	 * \return none.
	 */
	void freeSpectralOrderVector();
	/*! 
	 * \sa method unsigned getGainBiasNoise();
	 * \brief returns a pointer to the GainBiasNoise class instance.
	 * \return GainBiasNoise pointer.
	 */
	GainBiasNoise *getGainBiasNoise();
	/*! 
	 * unsigned getBarycentricRadialVelocityCorrection();
	 * \brief returns a double BarycentricRadialVelocityCorrection.
	 */
	double getBarycentricRadialVelocityCorrection();	
	/*!
	 * void getBarycentricRadialVelocityCorrection(double BarycentricRadialVelocityCorrection);
	 * \brief sets the double BarycentricRadialVelocityCorrection.
	 */
	void setBarycentricRadialVelocityCorrection(double BarycentricRadialVelocityCorrection);
	
	/*! 
	 * \sa method unsigned getCount();
	 * \brief returns the count of spectral orders that have content.
	 * \return unsigned - count.
	 */
	unsigned getCount();
	
	/*! 
	 * \sa method unsigned getMinorder();
	 * \brief returns the least order number in the vector
	 * \return unsigned - minorder.
	 */
	unsigned getMinorder();
	
	/*! 
	 * \sa method unsigned getMaxorder();
	 * \brief returns the maximal order number in the vector
	 * \return unsigned - minorder.
	 */
	unsigned getMaxorder();
	
	/*! 
	 * \sa method unsigned setCount();
	 * \brief sets the count of spectral orders that have content.
	 * \return none.
	 */
	void setCount(unsigned Count);
	
	/*! 
	 * \sa void setObject(string object)
	 * \brief sets the object name
	 * \return none.
	 */
	void setObject(string object);
	
	/*! 
	 * \sa string getObject(void);
	 * \brief get the object name
	 * \return none.
	 */
	string getObject(void);
	
	/*! 
	 * \sa void setInstrumentmode(instrumentmode_t Instrumentmode)
	 * \brief sets the Instrumentmode
	 * \return none.
	 */
	void setInstrumentmode(instrumentmode_t Instrumentmode) { instrumentmode = Instrumentmode; }
	
	/*! 
	 * \sa string getInstrumentmode(void);
	 * \brief get the getInstrumentmode
	 */
	instrumentmode_t getInstrumentmode(void) { return instrumentmode; }
	
	/*! 
	 * \sa setSequence(unsigned sequence)
	 * \brief sets the sequence number
	 * \return none.
	 */
	void setSequence(unsigned Sequence);
	
	/*! 
	 * \sa unsigned getSequence(void);
	 * \brief get the sequence number
	 * \return none.
	 */
	unsigned getSequence(void);
	
	/*! 
	 * \sa method unsigned setMinorder();
	 * \brief sets the least order number in the vector
	 * \return none.
	 */
	void setMinorder(unsigned Minorder);
	
	/*! 
	 * \sa method void getMaxorder(unsigned Maxorder);
	 * \brief sets the maximal order number in the vector
	 * \return none.
	 */
	void setMaxorder(unsigned Maxorder);
	
	/*! 
	 * \sa method Polynomial *getOrderSpacingPolynomial(void);
	 * \brief gets the order spacing polynomial
	 * \return Polynomial *.
	 */
	Polynomial *getOrderSpacingPolynomial(void);	
	/*! 
	 * \sa method setOrderSpacingPolynomial(PolynomialCoeffs_t *pc);
	 * \brief sets the order spacing polynomial
	 * \return void.
	 */
	void setOrderSpacingPolynomial(PolynomialCoeffs_t *pc);
	/* 
	 * \sa method Polynomial *getDispersionPolynomial(unsigned index);
	 * \brief gets the dispersion polynomial
	 */
	LaurentPolynomial *getDispersionPolynomial(unsigned index);
	/* 
	 * setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc);
	 * \brief sets the dispersion polynomial
	 */
    void setDispersionPolynomial(unsigned index, const int MinorderOfLaurentPolynomial,const int MaxorderOfLaurentPolynomial, PolynomialCoeffs_t *pc);
	
	/*! 
	 * \sa method void ReadSpectralOrders(string Filename);
	 * \brief augment an existing vector with information from a file
	 * \param Filename - string.
	 * \return none.
	 */
	void ReadSpectralOrders(string Filename);
	/*! 
	 * \sa method operaSpectralOrder* operaSpectralOrderVector::GetSpectralOrder(unsigned order);
	 * \brief Gets an operaSpectralOrder* to a given order, else NULL
	 * \return operaSpectralOrder - pointer to the operaSpectralOrder.
	 */
	operaSpectralOrder* GetSpectralOrder(unsigned order);
	
	/* 
	 * void operaSpectralOrderVector::ReadSpectralOrders(string Filename, operaSpectralOrder_t Format)
	 * \brief Reads from an m.fits product to create spectralorders
	 */
	bool ReadSpectralOrders(string Filename, operaSpectralOrder_t Format);
	/*! 
	 * \sa method void WriteSpectralOrder(string Filename, operaSpectralOrder_t Format, unsigned order=0, unsigned min=0);
	 * \details Writes a SpectralOrder to a File
	 * \details in the right place such that the vector remains ordered.
	 * \details the optional order argument permits incremental addition to the output file, where zero means write all.
	 * \return operaSpectralOrderVector* - pointer to the updated vector.
	 */
	void WriteSpectralOrders(string Filename, operaSpectralOrder_t Format, unsigned order=0, unsigned min=0);
	/*
	 * Read support routines...
	 */
	void readOrdersFromSNR(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);

	void readOrdersFromGeometry(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromWavelength(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromProfile(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromLines(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromCalibratedSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromCalibratedBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder);

	void readOrdersFromCalibratedExtendedBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder);

	void readOrdersFromFluxCalibrationBeamSpectrum(string filename, operaSpectralOrder_t format, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromLibreEspritSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
    
    void readOrdersFromLibreEspritPolarimetry(string filename, stokes_parameter_t StokesParameter, unsigned &count, unsigned &maxorder);
	
	void readOrdersFromReferenceSpectrum(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readDispersionPolynomial(string filename);
	
	void readOrderSpacingPolynomial(string filename);
	
	void readRadialVelocityCorrection(string filename);
	
	void readOrdersFromAperture(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
	
	void readOrdersFromPolar(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);
    
	void readOrdersFromExtendedPolarimetry(string filename, unsigned &count, unsigned &minorder, unsigned &maxorder);

	void readOrdersFromCSV(string filename);
	
	void readGainNoise(string filename);
	
    void fitOrderSpacingPolynomial(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, unsigned sampleCenterPosition, int detectionMethod, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2, ostream *pout);

    void measureIPAlongRowsFromSamples(operaFITSImage &masterFlatImage, operaFITSImage &badpixImage, float slit, unsigned nsamples, bool FFTfilter, float gain, float noise, unsigned x1, unsigned x2, unsigned y1, unsigned y2,float *ipfunc, float *ipx, float *iperr);
    
    unsigned getElemIndexAndOrdersByWavelength(int *orderWithReferenceFluxForNormalization, unsigned *elemIndexWithReferenceFluxForNormalization, double wavelength);

    void measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut);

    void measureContinuumAcrossOrders(unsigned binsize, int orderBin, unsigned nsigcut, unsigned nOrdersPicked, int *orderForWavelength);
    
    void FitFluxCalibrationAcrossOrders(int lowOrderToClip, int highOrderToClip, int orderBin, bool throughput);
    
    void getContinuumFluxesForNormalization(double *uncalibratedContinuumFluxForNormalization, double uncalibratedContinuumBeamFluxForNormalization[MAXNUMBEROFBEAMS],unsigned binsize, int orderBin, unsigned nsigcut);

	unsigned getLEElementCount(string LEfluxCalibration);
	
	void readLEFluxCalibration(string LEfluxCalibration, operaSpectralElements *fluxCalibrationElements);
};

/*
 * now define spectralordervector iterators, only supports iterating in a positive direction
 * from low to higher order
 */
#include <iterator>

class operaSpectralOrderVectorIterator : public iterator<input_iterator_tag, operaSpectralOrderVector> {
	operaSpectralOrderVector *ovp;
public:
	operaSpectralOrderVectorIterator(operaSpectralOrderVector *o) : ovp(o) {}
	operaSpectralOrderVectorIterator(const operaSpectralOrderVectorIterator &it) : ovp(it.ovp) {}
	operaSpectralOrderVectorIterator& operator ++() {ovp++; return *this;}
	bool operator ==(const operaSpectralOrderVectorIterator & rhs) {return ovp == rhs.ovp;}
	bool operator !=(const operaSpectralOrderVectorIterator & rhs) {return ovp != rhs.ovp;}
	operaSpectralOrderVector& operator *() {return *ovp;}
};

#endif

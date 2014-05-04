#ifndef OPERAFITSPRODUCT_H
#define OPERAFITSPRODUCT_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaFITSProduct
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Jul/2011
 
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

class operaMEFFITSProduct;

#include <string>
#include <fitsio.h>
#include "libraries/operaFITSImage.h"
#include "libraries/operaEspadonsImage.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaSpectralOrderVector.h"

#define MODE_STAR_ONLY_COLS 12
#define MODE_POLAR_COLS 24
#define MODE_STAR_PLUS_SKY_COLS 28
#define MAXPRODUCTCOLS 28

/*! 
 * operaFITSProduct
 * \author Doug Teeple
 * \brief This class extands the FITS image to an opera FITS Product image.
 * \file operaFITSProduct.h
 * \ingroup libraries
 */

class operaFITSProduct : public operaFITSImage {
	
private:
	typedef operaFITSImage& super;	// a way of referring to the super class
	instrumentmode_t instrumentmode;
	string colnames[MAXPRODUCTCOLS];
	unsigned columns, datapoints, next_row;
	
public:
	/*! 
	 * \sa class operaFITSProduct()
	 * \brief Basic operaFITSProduct constructor.
	 * \note extends operaFITSImage
	 * \return none
	 */
	operaFITSProduct(void);
	
	/*! 
	 * \sa class operaFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector)
	 * \brief operaFITSProduct constructconsturcted from a spectral order vector.
	 * \note extends operaFITSImage
	 * \param Filename
	 * \param spectralOrerVector
	 * \return none
	 */
	operaFITSProduct(string Filename, operaSpectralOrderVector *spectralOrerVector);
	
	/* 
	 * \sa class operaFITSProduct(operaMEFFITSProduct &mefproduct)
	 * \brief operaFITSProduct constructconsturcted from a spectral order vector.
	 * \note extends operaFITSImage
	 * \param Filename
	 * \param spectralOrerVector
	 * \return none
	 */
	operaFITSProduct(operaMEFFITSProduct &mefproduct);
	
	/*! 
	 * \sa class operaFITSProduct(unsigned Columns, unsigned Rows)
	 * \brief Basic unnames operaFITSProduct constructor with a size.
	 * \note extends operaFITSImage
	 * \param Columns
	 * \param Rows
	 * \return none
	 */
	operaFITSProduct(unsigned Columns, unsigned Rows);
	
	/*! 
	 * \sa class operaFITSProduct(string Filename, int mode=READWRITE|READONLY)
	 * \brief Constructor for readng a FITS file and creating the corresponding object.
	 * \param Filename
	 * \param mode
	 * \throws operaException operaErrorHeaderProblem
	 * \return none
	 */
	operaFITSProduct(string Filename, int mode=READWRITE/*READONLY*/, unsigned Compression = 0);
	/*! 
	 * \sa class operaFITSProduct(string Filename, instrumentmode_t, Instrumentmode, int mode=READWRITE|READONLY);
	 * \brief Constructor for creating a FITS product file with the correct column width.
	 * \param Filename
	 * \param mode
	 * \param Instrumentmode
	 * \throws operaException operaErrorHeaderProblem
	 * \return none
	 */
	operaFITSProduct(string Filename, instrumentmode_t Instrumentmode, unsigned Compression = 0);
	/*! 
	 * \sa class operaFITSProduct(string Filename, string baseOnFilename, unsigned Columns, unsigned Rows, unsigned Compression = 0)
	 * \brief Constructor for readng a FITS object file and creating the corresponding product.
	 * \param Filename the product file to create
	 * \param baseOnFilename the object file from which to get the headers
	 * \throws operaException operaErrorHeaderProblem
	 * \return none
	 */
	operaFITSProduct(string Filename, string baseOnFilename, instrumentmode_t Instrumentmode, unsigned Columns, unsigned Rows, unsigned Compression = 0);
	/*! 
	 * \sa class operaFITSProduct(string Filename, int mode=READWRITE|READONLY)
	 * \brief Constructor for readng a FITS object file and creating the corresponding product.
	 * \param Filename the product file to create
	 * \throws operaException operaErrorHeaderProblem
	 * \return none
	 */
	operaFITSProduct(string Filename, instrumentmode_t Instrumentmode, unsigned Columns, unsigned Rows, unsigned Compression = 0);
	/*! 
	 * \sa class operaFITSProduct()
	 * \brief destructor
	 * \return void
	 */
	~operaFITSProduct(void);
	
	/*
	 * Helper functions
	 */
	/*
	 * Set the default column names
	 */
	void updateColumnNames(void);
	/*
	 * Set the column names into the fits header
	 */
	void setHeaderColumnNames(void);
	
	/*
	 * Interface to spectra in text files compatible with Libre-Esprit
	 */
	
	/*! 
	 * \sa class operaFITSProduct()
	 * \brief addRow - add a row of data to a FITS product
	 * \param rowdata - float *
	 * \return void
	 */
	void addRow(float *rowdata, unsigned Columns, unsigned startColumn);
	/*! 
	 * \sa class operaFITSProduct()
	 * \brief getNextRow - get the current row
	 * \return void
	 */
	unsigned getNextRow(void);	
	/*! 
	 * \sa class operaFITSProduct()
	 * \brief setNextRow - set the current row
	 * \param Row - unsigned
	 * \return void
	 */
	void setNextRow(unsigned Row);
	/*! 
	 * \sa class operaFITSProduct()
	 * \brief readText - read a .s spectrum
	 * \param Filename - string
	 * \return void
	 */
	
	void readText(string Filename);
	
	/*! 
	 * \sa class operaFITSProduct()
	 * \brief writeText - write a .s spectrum
	 * \param Filename - string
	 * \return void
	 */
	void writeText(string Filename);
	
	/*
	 * getters / setters
	 */
	
	/*! 
	 * string getcolname(unsigned col) 
	 * \brief returns the string column name of col.
	 * \param col
	 * \return string
	 */
	string getcolname(unsigned col);
	
	/*! 
	 * unsigned getcols() 
	 * \brief returns number of columns.
	 * \return string
	 */
	unsigned getcols();
	
	/*! 
	 * unsigned getrows() 
	 * \brief returns number of rows.
	 * \return string
	 */
	unsigned getrows();
	/*! 
	 * unsigned getNeextRow() 
	 * \brief returns current row.
	 * \return string
	 */
	instrumentmode_t getmode();
	/*!
	 * instrumentmode_t getmode() 
	 * \brief returns intrument mode.
	 * \return string
	 */
	void setmode(instrumentmode_t Mode);
};
#endif

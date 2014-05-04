/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: Polynomial
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

/*!
 * Polynomial
 * \author Doug Teeple / Eder Martioli
 * \brief This class encapsulates the polynomial object.
 * \file Polynomial.cpp
 * \ingroup libraries
 */

#include <math.h>						// for pow

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/Polynomial.h"
#include "libraries/operaLibCommon.h"

/*
 * Constructors / Destructors
 */
Polynomial::Polynomial() :
orderOfPolynomial(MAXPOLYNOMIAL),
polychisqr(0.0)
{
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		polynomialVector[i] = 0.0;
		polynomialErrors[i] = 0.0;
	}
}

Polynomial::Polynomial(const unsigned OrderOfPolynomial) : 
polychisqr(0.0)
{
	orderOfPolynomial = OrderOfPolynomial;
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		polynomialVector[i] = 0.0;
		polynomialErrors[i] = 0.0;
	}
}

Polynomial::Polynomial(const unsigned OrderOfPolynomial, const double* CoefficientVector) :
polychisqr(0.0)
{
	orderOfPolynomial = OrderOfPolynomial;
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		polynomialVector[i] = CoefficientVector[i];
	}
}

Polynomial::Polynomial(const unsigned OrderOfPolynomial, const double* CoefficientVector, const double* CoefficientErrorVector) :
polychisqr(0.0)
{
	orderOfPolynomial = OrderOfPolynomial;
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		polynomialVector[i] = CoefficientVector[i];
		polynomialErrors[i] = CoefficientErrorVector[i];
	}
}

Polynomial::Polynomial(const PolynomialCoeffs_t* Coefficients) :
polychisqr(0.0)
{
	orderOfPolynomial = Coefficients->orderofPolynomial;
	polychisqr = Coefficients->polychisqr;
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		polynomialVector[i] = Coefficients->p[i];
		polynomialErrors[i] = Coefficients->e[i];
	}
}

Polynomial::Polynomial(const doublePolynomialCoeffs_t* Coefficients) :
polychisqr(0.0)
{
	orderOfPolynomial = Coefficients->orderofPolynomial;
	polychisqr = Coefficients->polychisqr;
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		polynomialVector[i] = Coefficients->p[i];
		polynomialErrors[i] = Coefficients->e[i];
	}
}

Polynomial::~Polynomial() {
}

/*
 * double Get(const unsigned index)
 * \brief This function gets the polynomial value at index.
 * \brief usage: float p3 = Get<float>(3);
 * \brief usage: double p3 = Get3);
 * \param index is a const unsigned index of the coefficient to get
 * \return void
 */
double Polynomial::Get(const unsigned index) {
	if (index > MAXPOLYNOMIAL) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (index > orderOfPolynomial) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return polynomialVector[index];
}

/*
 * void Set(const double x, const unsigned index)
 * \brief This function sets the polynomial value at index.
 * \brief usage: Set<float>((float)x, 3);
 * \brief usage: Set((double)x, 3);
 * \param x is a double value of the polynomial coefficient at index
 * \return void
 */
void Polynomial::Set(const double x, const unsigned index) {
	if (index > MAXPOLYNOMIAL) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (index > orderOfPolynomial) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	polynomialVector[index] = x;
}

/*
 * SimpletType Evaluate(double x)
 * \brief This function returns the value of a given polynomial function.
 * \brief usage: float value = Evaluate<float>((float)x);
 * \brief usage: double value = Evaluate((double)x);
 * \param x is a double input value for which the given polynomial is evaluated
 * \return double value of evaluation of the polynomial
 */
double Polynomial::Evaluate(const double x) {
	double fpoly = polynomialVector[0];
	for (unsigned i=1; i<orderOfPolynomial; i++) {
		fpoly += polynomialVector[i]*pow(x, i);
	}
	return fpoly;
}

/*
 * SimpletType double** getVector();
 * \brief This function returns the double *polynomialVector.
 * \brief usage: float *vec = getvector<float>();
 * \brief usage: double *vec = getvector();
 * \return double * - the vector of polynomial order coefficients
 */
double* Polynomial::getVector() {
	return polynomialVector;
}

/*
 * double double* getErrorVector();
 * \brief This function returns the double *polynomialErrorVector.
 * \brief usage: double *errs = getvector();
 * \return double * - the vector of polynomial coefficient errors
 */
double* Polynomial::getErrorVector() {
	return polynomialErrors;
}

/*
 * unsigned getOrderOfPolynomial();
 * \brief This function returns the unisgned polynomial Order.
 * \brief usage: unsigned npar = getOrderOfPolynomial<float>();
 * \brief usage: unsigned npar = getOrderOfPolynomial();
 * \return double * - the vector of polynomial order coefficients
 */
unsigned Polynomial::getOrderOfPolynomial() {
	return orderOfPolynomial;
}

/*
 * void setOrderOfPolynomial(unsigned Order);
 * \brief This function sets the unsigned polynomial Order.
 * \brief usage: setOrderOfPolynomial(3);
 * \return void
 */
void Polynomial::setOrderOfPolynomial(unsigned Order) {
	orderOfPolynomial = Order;
}

/*
 * void double getCoefficient(unsigned index);
 * \brief This function sets a coefficent at the index.
 * \brief usage: float coeff3 = getCoefficient(3);
 * \return void
 */
double Polynomial::getCoefficient(unsigned index) {
	if (index > MAXPOLYNOMIAL) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (index > orderOfPolynomial) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return polynomialVector[index];
}

/*
 * void setCoefficient(unsigned index, double value);
 * \brief This function sets a coefficent at the index.
 * \brief usage: setOrderOfPolynomial(3);
 * \return void
 */
void Polynomial::setCoefficient(unsigned index, double value) {
	if (index > MAXPOLYNOMIAL) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (index > orderOfPolynomial) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	polynomialVector[index] = value;
}

/*
 * void double getCoefficientError(unsigned index);
 * \brief This function gets a coefficent error at the index.
 * \brief usage: double err = getCoefficientError(3);
 * \return double
 */
double Polynomial::getCoefficientError(unsigned index) {
	if (index > MAXPOLYNOMIAL) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (index > orderOfPolynomial) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	return polynomialErrors[index];
}

/*
 * void setCoefficientError(unsigned index, double value);
 * \brief This function sets a coefficent error at the index.
 * \brief usage: setCoefficientError(3, 0.001);
 * \return void
 */
void Polynomial::setCoefficientError(unsigned index, double value) {
	if (index > MAXPOLYNOMIAL) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (index > orderOfPolynomial) {
		throw operaException("Polynomial: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	polynomialErrors[index] = value;
}

/*
 * struct PolynomialCoeffs_t* getPolynomialCoeffs();
 * \brief This function returns a PolynomialCoeffs_t struct.
 * \brief usage: PolynomialCoeffs_t *p = getPolynomialCoeffs<float>();
 * \note allocates storage that must be freed
 * \return PolynomialCoeffs_t  * - the PolynomialCoeffs_t struct *
 */
PolynomialCoeffs_t* Polynomial::getPolynomialCoeffs() {
	PolynomialCoeffs_t *pcoefficients = (PolynomialCoeffs_t *)malloc(sizeof(PolynomialCoeffs_t));
	pcoefficients->orderofPolynomial = orderOfPolynomial;
	pcoefficients->polychisqr = polychisqr;
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		pcoefficients->p[i] = polynomialVector[i];
		pcoefficients->e[i] = polynomialErrors[i];
	}
	return pcoefficients;
}	

/*
 * void setPolynomialCoeffs();
 * \brief This function sets polynomial values from a PolynomialCoeffs_t struct.
 * \brief usage: setPolynomialCoeffs<float>();
 * \return void
 */
void  Polynomial::setPolynomialCoeffs(PolynomialCoeffs_t* pcoefficients) {
	orderOfPolynomial = pcoefficients->orderofPolynomial;
	polychisqr = pcoefficients->polychisqr;
	for (unsigned i=0; i<orderOfPolynomial; i++) {
		polynomialVector[i] = pcoefficients->p[i];
		polynomialErrors[i] = pcoefficients->e[i];
	}
}	

/*
 * void setChisqr(double Chisqr)
 * \brief This function sets chisqr.
 * \brief usage: setChisqr(0.98);
 * \return void
 */

void Polynomial::setChisqr(double Chisqr) {
	polychisqr = Chisqr;
}

/*
 * double getChisqr(void)
 * \brief This function gets chisqr.
 * \brief usage: double c = getChisqr();
 * \return double
 */
double Polynomial::getChisqr(void) {
	return polychisqr;
}

/*!
 * void printEquation(ostream *pout)
 * \brief This function prints the polynomial equation in Gnuplot format.
 * \note usage: printEquation(&cout);
 * \return void
 */
void Polynomial::printEquation(ostream *pout) {
    if (pout != NULL) {
        *pout << "f(x) =";
        for(unsigned i=0;i<orderOfPolynomial;i++) {
            if(i==0) {
                *pout << " " << polynomialVector[i];
            } else if (i==1) {
                *pout << " + " << polynomialVector[i] << "*x";
            } else {
                *pout << " + " << polynomialVector[i] << "*x**" << i;
            }
        }
        *pout << endl;
    }
}



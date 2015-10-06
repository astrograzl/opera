#ifndef OPERAFLUXVECTOR_H
#define OPERAFLUXVECTOR_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaFluxVector
 Version: 1.0
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Mar/2012
 
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

#include <cmath>	// for sqrt
#include <utility>	// for pair

#include "operaError.h"
#include "libraries/operaException.h"

/*!
 * \file operaFluxVector.h
 */

using namespace std;

/*!
 * \brief Definition of the value of each "tends towards" optional field.
 */
typedef enum {ToDefault, ToINF, ToNAN, ToZero, ToOne} TendsTowards_t;

/*!
 * \author Doug Teeple
 * \author Andre Venne
 * \brief This class encapsulates the flux vector.
 * \ingroup libraries
 * \details
 * 
 * The flux vector stores fluxes and variances. The operators are defined to
 * operate on operaFluxVectors and propagate the variances.
 * 
 * The variances are calculated as followed :
 * 
 * F = F(a,b)
 * DF = Pow(dF/da,2) * Da + Pow(dF/db,2) *Db
 * 
 * where DF is the resulting variance, Da and Db are the variance of the fluxes a and b, dF/da and dF/db are the partial derivatives of F.
 * The fluxes are supposed uncorrelated.
 *
 * The flux vector also has an optional "tends towards" field which allows control of INF / INF situations,
 * where the results can tend towards INF, NaN, 0.0, 1.0 or default result.
 */
class operaFluxVector {
	
private:
	unsigned length;
    TendsTowards_t towards;
	double *fluxes;
	double *variances;
	
public:
	/*!
     * \brief operaFluxVector constructor.
     * \details This constructor creates an operaFluxVector.
     * \param Length Number of elements in the operaFluxVector
     * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(unsigned Length, TendsTowards_t Towards=ToDefault);
    
    /*!
     * \brief operaFluxVector constructor from flux and variance arrays.
     * \details This constructor creates an operaFluxVector by copying the contents of two equal sized arrays holding the flux and variance.
     * \param Fluxes An array with of the given length
     * \param Variances An array with of the given length
     * \param Length Number of elements in the operaFluxVector
     * \param Towards An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(double *Fluxes, double *Variances, unsigned Length, TendsTowards_t Towards=ToDefault);
    
    /*!
     * \brief operaFluxVector constructor from an operaFluxVector.
     * \details This constructor creates an operaFluxVector by copying the contents of an existing operaFluxVector.
     * \param b An operaFluxVector to be copied
     * \param Towards An optional TendsTowards_t value defaults to ToDefault
     */
	operaFluxVector(const operaFluxVector &b, TendsTowards_t Towards=ToDefault);
	
    /*!
     * \brief Destructor.
     */
	~operaFluxVector(void);
	
    /*!
     * \brief Sets the flux vector.
     * \details Copies the contents of an existing flux array to the flux vector.
     * \param Fluxes An array with the same length as the existing flux vector
     * \return void
     */
	void setFluxVector(double *Fluxes);
	
    /*!
     * \brief Sets the variance vector.
     * \details Copies the contents of an existing variance array to the variance vector.
     * \param Variances An array with the same length as the existing variance vector
     * \return void
     */
	void setVarianceVector(double *Variances);
	
    /*!
     * \brief Sets the flux and variance vectors.
     * \details Copies the contents of existing flux and variance arrays to the variance vector.
     * \param Fluxes An array with the same length as the existing flux vector
     * \param Variances An array with the same length as the existing variance vector
     * \return void
     */
	void setVectors(double *Fluxes, double *Variances);
	
    /*!
     * \brief Sets the flux and variance vectors from an operaFluxVector.
     * \details Copies the contents of flux and variance vectors from an existing operaFluxVector.
	 * \param fluxvector An operaFluxVector with the same length as the existing vectors
     * \return void
     */
	void setVector(const operaFluxVector &fluxvector);
	
	/*!
     * \brief Resizes the operaFluxVector.
     * \details Resizes the flux and variance vectors, and sets length to the new size. Elements are preserved up to min(length, newlength).
	 * \param newlength The new length to resize to
     * \return void
     */
	void resize(unsigned newlength);
	
    /*!
     * \brief Gets the flux vector.
     * \return A pointer to the flux vector
     */
	double* getfluxes() { return fluxes; }
	
    /*!
     * \brief Gets the variance vector.
     * \return A pointer to the variance vector
     */
	double* getvariances() { return variances; }
	
    /*!
     * \brief Gets a flux vector element.
     * \param index The index to get
     * \return The flux at index
     */
	double getflux(unsigned index) const { return fluxes[index]; }
	
    /*!
     * \brief Gets a variance vector element.
     * \param index The index to get
     * \return The variance at index
     */
	double getvariance(unsigned index) const { return variances[index]; }
	
    /*!
     * \brief Sets a flux vector element.
     * \param Flux The new flux value
     * \param index The index to set
     * \return void
     */
	void setflux(double Flux, unsigned index) { fluxes[index] = Flux; }
    
    /*!
     * \brief Sets a variance vector element.
     * \param Variance The new variance value
     * \param index The index to set
     * \return void
     */
	void setvariance(double Variance, unsigned index) { variances[index] = Variance; }
    
    /*!
     * \brief Gets the number of elements in the operaFluxVector.
     * \return The length
     */
	unsigned getlength() const { return length; }
    
    /*!
     * \brief Gets the TendsTowards_t value of operaFluxVector.
     * \return The TendsTowards_t value to which the operaFluxVector will tend to in case of infinity operations.
     */
    TendsTowards_t gettowards(void) const { return towards; }
	
    /*!
     * \brief Gets the error of a flux element.
     * \details Gets the square root of a variance vector element, which is the error of the flux element.
     * \param index The index to get
     * \return The error at index
     */
	double geterror(unsigned index) const { return sqrt(variances[index]); }
	
	/*! 
	 * \brief Indexing operator.
	 * \param i The index to get
	 * \note Usage: pair<double,double> p = FluxVector[i];
     * \note To print : cout << p.first << p.second;
	 * \return A pair of doubles containing the flux and variance at the given index
	 */
	pair<double,double> operator[](unsigned index) const {return pair<double,double>(fluxes[index], variances[index]); }
	
    /*! 
	 * \brief Assignment operator.
     * \details The operator copies the fluxes and variances from the right side.
     * \details If the lengths are not equal, the vectors will be reallocated using the length of the right hand side.
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator=(const operaFluxVector& b);
    
    /*! 
	 * \brief Assignment operator.
     * \details The operator every element of the flux vector to the given double and every element of the variance vector to 0.
     * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator=(double d);
	
    /*!
	 * \brief Addition/assignment operator.
     * \details The operator adds and copies the elements on the right side of the operator to the corresponding elements on the left side.
     * \details The variances are updated accordingly.
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */	
	operaFluxVector& operator+=(const operaFluxVector& b);
	
	/*! 
	 * \brief Addition/assignment operator.
     * \details The operator adds and copies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator+=(double d);
	
	/*!
	 * \brief Subtraction/assignment operator.
     * \details The operator subtracts and copies the elements on the right side of the operator to the corresponding elements on the left side.
     * \details The variances are updated accordingly.
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator-=(const operaFluxVector& b);
	
    /*! 
	 * \brief Subtraction/assignment operator.
     * \details The operator subtracts and copies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator-=(double d);
	
    /*!
	 * \brief Multiplication/assignment operator.
     * \details The operator multiplies and copies the elements on the right side of the operator to the corresponding elements on the left side.
     * \details The variances are updated accordingly.
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator*=(const operaFluxVector& b);
	
    /*! 
	 * \brief Multiplication/assignment operator.
     * \details The operator multiplies and copies the double from the right side of the operator to every value of the flux vector on the left side.
     * \details The variances are updated accordingly.
	 * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator*=(double d);
	
    /*!
	 * \brief Division/assignment operator.
     * \details The operator divides the elements on the left side of the operator by the corresponding elements on the right side and copies them to the left side.
     * \details The variances are updated accordingly.
	 * \param b The rhs operaFluxVector
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator/=(const operaFluxVector& b);
	
    /*! 
	 * \brief Division/assignment operator.
     * \details The operator divides every value of the flux vector on the left side of the operator by the double on the right side and copies them to the left side.
     * \details The variances are updated accordingly.
	 * \param d The flux value
	 * \return A reference to the operaFluxVector
	 */
	operaFluxVector& operator/=(double d);
};

/*!
 * \brief Addition operator.
 * \details Adds two operaFluxVectors without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting sum operaFluxVector
 */
inline operaFluxVector operator+(operaFluxVector a, const operaFluxVector& b) { return a += b; }

/*! 
 * \brief Addition operator.
 * \details Adds a flux value to an operaFluxVector without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting sum operaFluxVector
 */
inline operaFluxVector operator+(operaFluxVector a, double d) { return a += d; }

 /*!
 * \brief Subtraction operator.
 * \details Subtracts one operaFluxVector from another without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting difference operaFluxVector
 */
inline operaFluxVector operator-(operaFluxVector a, const operaFluxVector& b) { return a -= b; }

/*! 
 * \brief Subtraction operator.
 * \details Subtracts a flux value from an operaFluxVector without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting difference operaFluxVector
 */
inline operaFluxVector operator-(operaFluxVector a, double d) { return a -= d; }

/*!
 * \brief Multiplication operator.
 * \details Multiplies two operaFluxVectors without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting product operaFluxVector
 */
inline operaFluxVector operator*(operaFluxVector a, const operaFluxVector& b) { return a *= b; }

/*! 
 * \brief Multiplication operator.
 * \details Multilpies an operaFluxVector by a flux value without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting product operaFluxVector
 */
inline operaFluxVector operator*(operaFluxVector a, double d) { return a*=d; }

/*!
 * \brief Division operator.
 * \details Divides one operaFluxVector by another without modifying either operaFluxVector.
 * \param a An operaFluxVector
 * \param b An operaFluxVector
 * \return The resulting quotient operaFluxVector
 */
inline operaFluxVector operator/(operaFluxVector a, const operaFluxVector& b) { return a /= b; }

/*! 
 * \brief Division operator.
 * \details Divides an operaFluxVector by a flux value without modifying it.
 * \param a An operaFluxVector
 * \param d The flux value
 * \return The resulting quotient operaFluxVector
 */
inline operaFluxVector operator/(operaFluxVector a, double d) { return a/=d; }

/*!
 * \brief Square root function.
 * \details The function applies a square root to every element in the flux vector.
 * \param b An operaFluxVector
 * \return The resulting square root operaFluxVector
 */
operaFluxVector Sqrt(operaFluxVector b);

/*!
 * \brief Power function.
 * \details The function raises every element in the flux vector to the specified power.
 * \param b An operaFluxVector
 * \param d The power
 * \return The resulting power operaFluxVector
 */
operaFluxVector Pow(operaFluxVector b, double d);

/*!
 * \brief Sum function.
 * \details Calculuates the sum of a flux vector.
 * \param b An operaFluxVector
 * \return A pair of doubles containing the sums for the fluxes and variances respspectively
 */
pair<double,double> Sum(const operaFluxVector& b);

/*!
 * \brief Mean function.
 * \details Calculuates the mean of a flux vector.
 * \param b An operaFluxVector
 * \return A pair of doubles containing the means of the fluxes and variances respspectively
 */
pair<double,double> Mean(const operaFluxVector& b);

/*!
 * \brief Resizes a vector.
 * \details Resizes a dynamic array, copying over existing elements that fit in the new range.
 * \param vector The existing vector of doubles to be resized
 * \param oldlength The previous length of the vector
 * \param newlength The new length to resize to
 * \return void
 */
void resizeVector(double *&vector, unsigned oldlength, unsigned newlength);

#endif

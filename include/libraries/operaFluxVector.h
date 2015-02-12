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
	bool istemp;
    TendsTowards_t towards;
	double *fluxes;
	double *variances;
	
public:
	/*
	 * Constructors / Destructors
	 */
    
    /*!
     * \brief operaFluxVector constructor.
     * \details This constructor creates an operaFluxVector with vectors of content 0.0.
     * \param Length An unsigned number of elements in the operaFluxVector
     * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
     * \param Istemp An optional bool defaults to false
     * \return void
     */
	operaFluxVector(unsigned Length, TendsTowards_t Towards=ToDefault, bool Istemp=false);
    
    /*!
     * \brief operaFluxVector constructor from a flux pointer and his variance pointer.
     * \details This constructor creates an operaFluxVector with vectors of content given by a flux double pointer and his variance double pointer.
     * \param Fluxes A double pointer
     * \param Variances A double pointer
     * \param Length An unsigned number of elements in the operaFluxVector
     * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
     * \param Istemp An optional bool defaults to false
     * \return void
     */
	operaFluxVector(double *Fluxes, double *Variances, unsigned Length, TendsTowards_t Towards=ToDefault, bool Istemp=false);
    
    /*!
     * \brief operaFluxVector constructor from an operaFluxVector.
     * \details This constructor creates an operaFluxVector with vectors of content given by an operaFluxVector.
     * \param b An operaFluxVector address
     * \param Length An unsigned number of elements in the operaFluxVector
     * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
     * \param Istemp An optional bool defaults to false
     * \return void
     */
	operaFluxVector(operaFluxVector &b, TendsTowards_t Towards=ToDefault, bool Istemp=false);
	
    /*!
     * \brief operaFluxVector constructor from an operaFluxVector.
     * \details This constructor creates an operaFluxVector with vectors of content given by an operaFluxVector.
     * \param b An operaFluxVector pointer
     * \param Length An unsigned number of elements in the operaFluxVector
     * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
     * \param Istemp An optional bool defaults to false
     * \return void
     */
	operaFluxVector(operaFluxVector *b, TendsTowards_t Towards=ToDefault, bool Istemp=false);
	
    /*!
     * \brief Destructor releases the pointers memory.
     * \return void
     */
	~operaFluxVector(void);
	
    /*
	 * Getters/Setters
	 */
    
    /*!
     * \brief Sets the variance vector.
     * \details A function that sets the variance vector content to the value of a given double pointer.
     * \param Variances A double pointer
     * \return void
     */
	void setVarianceVector(double *Variances);
	
    /*!
     * \brief Sets the flux vector.
     * \details A function that sets the flux vector content to the value of a given double pointer.
     * \param Fluxes A double pointer
     * \return void
     */
	void setVector(double *Fluxes);
	
    /*!
     * \brief Sets the flux vector.
     * \details A function that sets the flux and variance vector content to the value of given double pointers.
	 * \param FluxVector An operaFluxVector address
     * \return void
     */
	void setVector(operaFluxVector &fluxvector);
	
    /*!
     * \brief Sets the flux and variance vector.
     * \details A function that sets the flux and variance vector content to the value of given double pointers.
     * \param Fluxes A double pointer
     * \param Variances A double pointer
     * \return void
     */
	void setVectors(double *Fluxes, double *Variances);
	
    /*!
     * \brief Gets the variance vector.
     * \details A function that gets the variance vector.
     * \return A double pointer
     */
	double* getvariances() { return variances; };
	
    /*!
     * \brief Gets the flux vector.
     * \details A function that gets the flux vector.
     * \return A double pointer
     */
	double* getfluxes() { return fluxes; };
	
    /*!
     * \brief Gets a variance vector element.
     * \details A function that gets a variance vector element.
     * \param index An unsigned value
     * \return A double value
     */
	double getvariance(unsigned index) { return variances[index]; };
	
    /*!
     * \brief Gets a flux vector element.
     * \details A function that gets a flux vector element.
     * \param index An unsigned value
     * \return A double value
     */
	double getflux(unsigned index) { return fluxes[index]; };    
    
    /*!
     * \brief Sets a variance vector element.
     * \details A function that sets a variance vector element.
     * \param Variance A double value
     * \param index An unsigned value
     * \return void
     */
	void setvariance(double Variance, unsigned index) {variances[index] = Variance; };
	
    /*!
     * \brief Sets a flux vector element.
     * \details A function that sets a flux vector element.
     * \param Flux A double value
     * \param index An unsigned value
     * \return void
     */
	void setflux(double Flux, unsigned index) {fluxes[index] = Flux; };     
    
    /*!
     * \brief Gets the number of elements in the operaFluxVector.
     * \details A function that gets the length of the vectors, which is the number of elements in the operaFluxVector.
     * \return An unsigned value
     */
	unsigned getlength() { return length; };
    
    /*!
     * \brief Gets the boolean value of operaFluxVector.
     * \details A function that gets the boolean value giving the temporary state of the operaFluxVector.
     * \return A boolean value
     */
    bool getistemp() { return istemp; };
    
    /*!
     * \brief Gets the TendsTowards_t value of operaFluxVector.
     * \details A function that gets the TendsTowards_t value to which the operaFluxVector will tend to in case of infinity operations.
     * \param index An unsigned value
     * \return A TendsTowards_t value
     */
    TendsTowards_t gettowards(void) { return towards; };
	
    /*!
     * \brief Gets the error of a flux element.
     * \details A function that gets the square root of a variance vector element, which is the error of the flux element.
     * \return A double value
     */
	double geterror(unsigned index) { return sqrt(variances[index]); };
	
	/*! 
	 * \brief Indexing operator.
     * \details The operator returns a flux / variance pair pointer of doubles of the vectors element.
	 * \param i An unsigned value
	 * \note Usage: pair<double,double>*p = FluxVector[i];
     * \note To print : printf("%.2f %.2f\n", p->first, p->second);
	 * \return A pair pointer of doubles
	 */
	pair<double,double>* operator[](unsigned i) {return new pair<double,double>(fluxes[i], variances[i]); };
    
	/*! 
	 * \brief Assignment operator.
     * \details The operator copies the fluxes and variances from the right side of the operator to the left side.
     * \details If the operaFluxVector on the right side is temporary, it is deleted.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector a = operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator=(operaFluxVector* b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b->variances; 
		if (length < b->length) {
			throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
		unsigned n = length; 
		while (n--) { *afluxes++ = *bfluxes++; *avariances++ = *bvariances++; }
		if (b->istemp) delete b;
		return *this;
	};
	
    /*! 
	 * \brief Assignment operator.
     * \details The operator copies the fluxes and variances from the right side of the operator to the left side.
     * \details If the operaFluxVector on the right side is temporary, it is deleted.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector a = operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator=(operaFluxVector& b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b.variances; 
		if (length < b.length) {
			throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
		}
		unsigned n = length; 
		while (n--) { *afluxes++ = *bfluxes++; *avariances++ = *bvariances++; }
		if (b.istemp) delete &b;
		return *this;
	};
    
    /*! 
	 * \brief Assignment operator.
     * \details The operator copies the double from the right side of the operator to every value of the flux vector on the left side.
     * \details It assigns the value 0 to the variances.
	 * \param d A double value
	 * \note Usage: operaFluxVector a = double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator=(double d) {
		double *afluxes = (double *)fluxes; 
		double *avariances = (double *)variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ = d; *avariances++ = 0.0; }
		return *this;
	};
	
    /*!
	 * \brief Addition/assignment operator.
     * \details The operator adds and copies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector a += operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator+=(operaFluxVector* b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ += *bfluxes++; *avariances++ += *bvariances++; }
		if (b->istemp) delete b;
		return *this;
	};
	
    /*!
	 * \brief Addition/assignment operator.
     * \details The operator adds and copies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector a += operaFluxVector b;
	 * \return An operaFluxVector address
	 */	
	operaFluxVector& operator+=(operaFluxVector& b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b.variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ += *bfluxes++; *avariances++ += *bvariances++; }
		if (b.istemp) delete &b;
		return *this;
	};
	
	/*! 
	 * \brief Addition/assignment operator.
     * \details The operator adds and copies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d A double value
	 * \note Usage: operaFluxVector a += double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator+=(double d) {
		double *afluxes = (double *)fluxes; 
		double *avariances = (double *)variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ += d; *avariances++ += 0.0; }
		return *this;
	};
	
    /*!
	 * \brief Subtraction/assignment operator.
     * \details The operator subtracts and copies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector a -= operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator-=(operaFluxVector* b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ -= *bfluxes++; *avariances++ += *bvariances++; }
		if (b->istemp) delete b;
		return *this;
	};
	
	/*!
	 * \brief Subtraction/assignment operator.
     * \details The operator subtracts and copies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector a -= operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator-=(operaFluxVector& b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b.variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ -= *bfluxes++; *avariances++ += *bvariances++; }
		if (b.istemp) delete &b;
		return *this;
	};
	
    /*! 
	 * \brief Subtraction/assignment operator.
     * \details The operator subtracts and copies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d A double value
	 * \note Usage: operaFluxVector a -= double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator-=(double d) {
		double *afluxes = (double *)fluxes; 
		double *avariances = (double *)variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ -= d; *avariances++ += 0.0; }
		return *this;
	};
	
    /*!
	 * \brief Multiplication/assignment operator.
     * \details The operator multiplies and copies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector a *= operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator*=(operaFluxVector* b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
		while (n--) { *afluxes *= *bfluxes; double td = *avariances; *avariances++ = ( pow(*afluxes,2) * *bvariances++ ) + ( pow(*bfluxes,2) * td ); afluxes++; bfluxes++; }
		if (b->istemp) delete b;
		return *this;
	};
	
    /*!
	 * \brief Multiplication/assignment operator.
     * \details The operator multiplies and copies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector a *= operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator*=(operaFluxVector& b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b.variances; 
		unsigned n = length; 
		while (n--) { *afluxes *= *bfluxes; double td = *avariances; *avariances++ = ( pow(*afluxes,2) * *bvariances++ ) + ( pow(*bfluxes,2) * td ); afluxes++; bfluxes++; }
		if (b.istemp) delete &b;
		return *this;
	};
	
    /*! 
	 * \brief Multiplication/assignment operator.
     * \details The operator multiplies and copies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d A double value
	 * \note Usage: operaFluxVector a *= double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator*=(double d) {
		double *afluxes = (double *)fluxes; 
		double *avariances = (double *)variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ *= d; *avariances++ *= pow(d,2); }
		return *this;
	};
	
    /*!
	 * \brief Division/assignment operator.
     * \details The operator divides the elements on the left side of the operator by the corresponding elements on the right side and copies them to the left side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector a /= operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator/=(operaFluxVector* b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
        switch (this->towards) {
            case ToDefault:
                while (n--) { *afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; }
                break;
            case ToINF:
				while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = FP_INFINITE; *avariances++ = 0.0; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToNAN:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = FP_NAN; *avariances++ = FP_NAN; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToZero:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = 0.0; *avariances++ = 0.0; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToOne:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = 1.0; *avariances++ = 0.0; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            default:
                break;
        }
		if (b->istemp) delete b;
		return *this;
	};
	
    /*!
	 * \brief Division/assignment operator.
     * \details The operator divides the elements on the left side of the operator by the corresponding elements on the right side and copies them to the left side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector a /= operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator/=(operaFluxVector& b) {
		double *afluxes = (double *)fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)variances; 
		double *bvariances = (double *)b.variances; 
		unsigned n = length; 
        switch (this->towards) {
            case ToDefault:
                while (n--) { *afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; }
                break;
            case ToINF:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = FP_INFINITE; *avariances++ = 0.0; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToNAN:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = FP_NAN; *avariances++ = FP_NAN; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToZero:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = 0.0; *avariances++ = 0.0; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToOne:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*afluxes++ = 1.0; *avariances++ = 0.0; bfluxes++; bvariances++;} else {*afluxes /= *bfluxes; double td = *avariances; *avariances++ = ( pow(*bfluxes,-2) * td) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
                
            default:
                break;
        }
		if (b.istemp) delete &b;
		return *this;
	};
	
    /*! 
	 * \brief Division/assignment operator.
     * \details The operator divides every value of the flux vector on the left side of the operator by the double on the right side and copies them to the left side.
	 * \param d A double value
	 * \note Usage: operaFluxVector a /= double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator/=(double d) {
		double *afluxes = (double *)fluxes; 
		double *avariances = (double *)variances; 
		unsigned n = length; 
		while (n--) { *afluxes++ /= d; *avariances++ /= pow(d,2); }
		return *this;
	};
	
    /*!
	 * \brief Multiplication operator.
     * \details The operator multiplies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector t = operaFluxVector a * operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator*(operaFluxVector* b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes; 
			tvariances = (double *)this->variances; 
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes; 
			tvariances = (double *)t->variances; 
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)this->variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
		while (n--) { *tfluxes++ = *afluxes * *bfluxes; *tvariances++ = ( pow(*afluxes,2) * *bvariances++ ) + ( pow(*bfluxes,2) * *avariances++ ); afluxes++; bfluxes++; }
		if (b->istemp) delete b;
		return *t;
	};
	
    /*!
	 * \brief Multiplication operator.
     * \details The operator multiplies the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector t = operaFluxVector a * operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator*(operaFluxVector& b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes; 
			tvariances = (double *)this->variances; 
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes; 
			tvariances = (double *)t->variances; 
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)this->variances; 
		double *bvariances = (double *)b.variances; 
		unsigned n = length; 
		while (n--) { *tfluxes++ = *afluxes * *bfluxes; *tvariances++ = ( pow(*afluxes,2) * *bvariances++ ) + ( pow(*bfluxes,2) * *avariances++ ); afluxes++; bfluxes++; }
		if (b.istemp) delete &b;
		return *t;
	};
	
    /*! 
	 * \brief Multiplication operator.
     * \details The operator multiplies the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d A double value
	 * \note Usage: operaFluxVector t = operaFluxVector a * double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator*(double d) {
        double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes;
			tvariances = (double *)this->variances;
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes;
			tvariances = (double *)t->variances;
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes;
		double *avariances = (double *)this->variances;
		unsigned n = length;
		while (n--) { *tfluxes++ = *afluxes++ * d; *tvariances++ = *avariances++ * pow(d,2); }
		return *t;
	};
	
    /*!
	 * \brief Division operator.
     * \details The operator divides the elements on the left side of the operator by the corresponding elements on the right side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector t = operaFluxVector a / operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator/(operaFluxVector* b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes; 
			tvariances = (double *)this->variances; 
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes; 
			tvariances = (double *)t->variances; 
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)this->variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
        switch (this->towards) {
            case ToDefault:
                while (n--) { *tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++;  }
                break;
            case ToINF:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = FP_INFINITE; *tvariances++ = 0.0; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToNAN:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = FP_NAN; *tvariances++ = FP_NAN; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToZero:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = 0.0; *tvariances++ = 0.0; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToOne:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = 1.0; *tvariances++ = 0.0; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
                
            default:
                break;
        }
		if (b->istemp) delete b;
		return *t;
	};
	
    /*!
	 * \brief Division operator.
     * \details The operator divides the elements on the left side of the operator by the corresponding elements on the right side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector t = operaFluxVector a / operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator/(operaFluxVector& b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes; 
			tvariances = (double *)this->variances; 
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes; 
			tvariances = (double *)t->variances; 
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)this->variances; 
		double *bvariances = (double *)b.variances; 
		unsigned n = length; 
        switch (this->towards) {
            case ToDefault:
                while (n--) { *tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++;  }
                break;
            case ToINF:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = FP_INFINITE; *tvariances++ = 0.0; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToNAN:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = FP_NAN; *tvariances++ = FP_NAN; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToZero:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = 0.0; *tvariances++ = 0.0; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
            case ToOne:
                while (n--) { if (isinf(*afluxes)&&isinf(*bfluxes)) {*tfluxes++ = 1.0; *tvariances++ = 0.0; afluxes++; bfluxes++; avariances++; bvariances++;} else {*tfluxes++ = *afluxes / *bfluxes; *tvariances++ = ( pow(*bfluxes,-2) * *avariances++ ) + ( pow(*afluxes / pow(*bfluxes,2),2) * *bvariances++ ); afluxes++; bfluxes++; } }
                break;
                
            default:
                break;
        }
		if (b.istemp) delete &b;
		return *t;
	};
	
    /*! 
	 * \brief Division operator.
     * \details The operator divides every value of the flux vector on the left side of the operator by the double on the right side.
	 * \param d A double value
	 * \note Usage: operaFluxVector t = operaFluxVector a / double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator/(double d) {
        double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes;
			tvariances = (double *)this->variances;
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes;
			tvariances = (double *)t->variances;
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes;
		double *avariances = (double *)this->variances;
		unsigned n = length;
		while (n--) { *tfluxes++ = *afluxes++ / d; *tvariances++ = *avariances++ / pow(d,2); }
		return *t;
	};
	
    /*!
	 * \brief Addition operator.
     * \details The operator adds the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector t = operaFluxVector a + operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator+(operaFluxVector* b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes; 
			tvariances = (double *)this->variances; 
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes; 
			tvariances = (double *)t->variances; 
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)this->variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
		while (n--) { *tfluxes++ = *afluxes++ + *bfluxes++; *tvariances++ = *avariances++ + *bvariances++; }
		if (b->istemp) delete b;
		return *t;
	};
	
    /*!
	 * \brief Addition operator.
     * \details The operator adds the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector t = operaFluxVector a + operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator+(operaFluxVector& b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes;
			tvariances = (double *)this->variances;
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes;
			tvariances = (double *)t->variances;
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes;
		double *bfluxes = (double *)b.fluxes;
		double *avariances = (double *)this->variances;
		double *bvariances = (double *)b.variances;
		unsigned n = length;
		while (n--) { *tfluxes++ = *afluxes++ + *bfluxes++; *tvariances++ = *avariances++ + *bvariances++; }
		if (b.istemp) delete &b;
		return *t;
	};
	
    /*! 
	 * \brief Addition operator.
     * \details The operator adds the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d A double value
	 * \note Usage: operaFluxVector t = operaFluxVector a + double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator+(double d) {
        double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes;
			tvariances = (double *)this->variances;
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes;
			tvariances = (double *)t->variances;
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes;
		double *avariances = (double *)this->variances;
		unsigned n = length;
		while (n--) { *tfluxes++ = *afluxes++ + d; *tvariances++ = *avariances++ + 0.0; }
		return *t;
	};
	
    /*!
	 * \brief Subtraction operator.
     * \details The operator subtracts the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector pointer
	 * \note Usage: operaFluxVector t = operaFluxVector a - operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator-(operaFluxVector* b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes; 
			tvariances = (double *)this->variances; 
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes; 
			tvariances = (double *)t->variances; 
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes; 
		double *bfluxes = (double *)b->fluxes; 
		double *avariances = (double *)this->variances; 
		double *bvariances = (double *)b->variances; 
		unsigned n = length; 
		while (n--) { *tfluxes++ = *afluxes++ - *bfluxes++; *tvariances++ = *avariances++ + *bvariances++; }
		if (b->istemp) delete b;
		return *t;
	};
	
    /*!
	 * \brief Subtraction operator.
     * \details The operator subtracts the elements on the right side of the operator to the corresponding elements on the left side.
	 * \param b An operaFluxVector address
	 * \note Usage: operaFluxVector t = operaFluxVector a - operaFluxVector b;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator-(operaFluxVector& b) {
		double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes; 
			tvariances = (double *)this->variances; 
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes; 
			tvariances = (double *)t->variances; 
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes; 
		double *bfluxes = (double *)b.fluxes; 
		double *avariances = (double *)this->variances; 
		double *bvariances = (double *)b.variances; 
		unsigned n = length; 
		while (n--) { *tfluxes++ = *afluxes++ - *bfluxes++; *tvariances++ = *avariances++ + *bvariances++; }
		if (b.istemp) delete &b;
		return *t;
	};
	
    /*! 
	 * \brief Subtraction operator.
     * \details The operator subtracts the double from the right side of the operator to every value of the flux vector on the left side.
	 * \param d A double value
	 * \note Usage: operaFluxVector t = operaFluxVector a - double d;
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator-(double d) {
        double *tfluxes, *tvariances;
		operaFluxVector *t = NULL;
		if (this->istemp) {
			t = this;
			tfluxes = (double *)this->fluxes;
			tvariances = (double *)this->variances;
		} else {
			t = new operaFluxVector(*this, towards, true);
			tfluxes = (double *)t->fluxes;
			tvariances = (double *)t->variances;
		}
		t->towards = this->towards;
		double *afluxes = (double *)this->fluxes;
		double *avariances = (double *)this->variances;
		unsigned n = length;
		while (n--) { *tfluxes++ = *afluxes++ - d; *tvariances++ = *avariances++ + 0.0; }
		return *t;
	};
};

/*!
 * \brief Square root function.
 * \details The function applies a square root to every element in the flux vector.
 * \param b An operaFluxVector pointer
 * \note Usage: operaFluxVector t = Sqrt( operaFluxVector b );
 * \return An operaFluxVector address
 */
static inline operaFluxVector& Sqrt(operaFluxVector* b) {
	operaFluxVector *t = new operaFluxVector(b->getlength(), b->gettowards(), true);
	unsigned n = b->getlength();
	double *bfluxes = b->getfluxes();
	double *tfluxes = t->getfluxes();
	double *bvariances = b->getvariances();
	double *tvariances = t->getvariances();
	while (n--) {
		*tfluxes++ = sqrt(*bfluxes);
		*tvariances++ = *bvariances++ / ( 4.0 * *bfluxes++ );
	}
	if (b->getistemp()) delete b;
	return *t;
}

/*!
 * \brief Square root function.
 * \details The function applies a square root to every element in the flux vector.
 * \param b An operaFluxVector address
 * \note Usage: operaFluxVector t = Sqrt( operaFluxVector b );
 * \return An operaFluxVector address
 */
static inline operaFluxVector& Sqrt(operaFluxVector& b) {
	operaFluxVector *t = new operaFluxVector(b.getlength(), b.gettowards(), true);
	unsigned n = b.getlength();
	double *bfluxes = b.getfluxes();
	double *tfluxes = t->getfluxes();
	double *bvariances = b.getvariances();
	double *tvariances = t->getvariances();
	while (n--) {
		*tfluxes++ = sqrt(*bfluxes);
		*tvariances++ = *bvariances++ / ( 4.0 * *bfluxes++ );
	}
	if (b.getistemp()) delete &b;
	return *t;
}

/*!
 * \brief Power function.
 * \details The function raises to specified power every element in the flux vector.
 * \param b An operaFluxVector pointer
 * \param d A double value
 * \note Usage: operaFluxVector t = Pow( operaFluxVector b , double d );
 * \return An operaFluxVector address
 */
static inline operaFluxVector& Pow(operaFluxVector* b, double d) {
	operaFluxVector *t = new operaFluxVector(b->getlength(), b->gettowards(), true);
	unsigned n = b->getlength();
	double *bfluxes = b->getfluxes();
	double *tfluxes = t->getfluxes();
	double *bvariances = b->getvariances();
	double *tvariances = t->getvariances();
	while (n--) {
		*tfluxes++ = pow(*bfluxes,d);
		*tvariances++ = pow(d * pow(*bfluxes++,d-1.0),2) * *bvariances++;
	}
	if (b->getistemp()) delete b;
	return *t;
}

/*!
 * \brief Power function.
 * \details The function raises to specified power every element in the flux vector.
 * \param b An operaFluxVector address
 * \param d A double value
 * \note Usage: operaFluxVector t = Pow( operaFluxVector b , double d );
 * \return An operaFluxVector address
 */
static inline operaFluxVector& Pow(operaFluxVector& b, double d) {
	operaFluxVector *t = new operaFluxVector(b.getlength(), b.gettowards(), true);
	unsigned n = b.getlength();
	double *bfluxes = b.getfluxes();
	double *tfluxes = t->getfluxes();
	double *bvariances = b.getvariances();
	double *tvariances = t->getvariances();
	while (n--) {
		*tfluxes++ = pow(*bfluxes,d);
		*tvariances++ = pow(d * pow(*bfluxes++,d-1.0),2) * *bvariances++;
	}
	if (b.getistemp()) delete &b;
	return *t;
}

/*!
 * \brief Sum function.
 * \details The function sums every element in the flux vector.
 * \param b An operaFluxVector pointer
 * \note Usage: double sum = Sum( operaFluxVector& b );
 * \note This function allocated storage which must be disposed of by the caller.
 * \return pointer to a pair of doubles
 */
static inline pair<double,double>Sum(operaFluxVector& b) {
	unsigned n = b.getlength();
	double sum = 0.0;
	double var = 0.0;
	double *pfluxes = b.getfluxes();
	double *pvariances = b.getvariances();
	while (n--) {
		sum += *pfluxes++;
		var += *pvariances++;
	}
	if (b.getistemp()) delete &b;
	return pair<double,double>(sum, var);
}

/*!
 * \brief Sum function.
 * \details The function sums every element in the flux vector.
 * \param b An operaFluxVector pointer
 * \note Usage: double sum = Sum( operaFluxVector& b ).first;
 * \note This function allocated storage which must be disposed of by the caller.
 * \return pointer to a pair of doubles
 */
static inline pair<double,double>Sum(operaFluxVector* b) {
	unsigned n = b->getlength();
	double sum = 0.0;
	double var = 0.0;
	double *pfluxes = b->getfluxes();
	double *pvariances = b->getvariances();
	while (n--) {
		sum += *pfluxes++;
		var += *pvariances++;
	}
	if (b->getistemp()) delete b;
	return pair<double,double>(sum, var);
}

/*!
 * \brief Mean function.
 * \details The function returns the mean of the flux vector.
 * \param b An operaFluxVector pointer
 * \note Usage: double mean = Mean( operaFluxVector& b ).first;
 * \note This function allocated storage which must be disposed of by the caller.
 * \return pointer to a pair of doubles
 */
static inline pair<double,double>Mean(operaFluxVector& b) {
	unsigned n = b.getlength();
	double sum = 0.0;
	double var = 0.0;
	double *pfluxes = b.getfluxes();
	double *pvariances = b.getvariances();
	while (n--) {
		sum += *pfluxes++;
		var += *pvariances++;
	}
	if (b.getistemp()) delete &b;
	sum /= b.getlength();
	var /= b.getlength();
	return pair<double,double>(sum, var);
}

/*!
 * \brief Mean function.
 * \details The function returns the mean of the flux vector.
 * \param b An operaFluxVector pointer
 * \note Usage: double mean = Mean( operaFluxVector& b ).first;
 * \note This function allocated storage which must be disposed of by the caller.
 * \return pointer to a pair of doubles
 */
static inline pair<double,double>Mean(operaFluxVector* b) {
	unsigned n = b->getlength();
	double sum = 0.0;
	double var = 0.0;
	double *pfluxes = b->getfluxes();
	double *pvariances = b->getvariances();
	while (n--) {
		sum += *pfluxes++;
		var += *pvariances++;
	}
	if (b->getistemp()) delete b;
	sum /= b->getlength();
	var /= b->getlength();
	return pair<double,double>(sum, var);
}
	
#endif

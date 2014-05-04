/*******************************************************************
 ****               		OPERA PIPELINE v1.0                     ****
 ********************************************************************
 Library name: operaStokesVector
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

#ifndef OPERASTOKESVECTOR_H
#define OPERASTOKESVECTOR_H

#include "globaldefines.h"
#include "libraries/operaFluxVector.h"			// for operaFluxVector

/*!
 * \file operaStokesVector.h
 * \brief This file holds the declaration of the class operaStokesVector.
 * \ingroup libraries
 */

using namespace std;

/*!
 * \brief Definition of the value of each Stokes parameter.
 */
typedef enum { StokesI=0, StokesQ, StokesU, StokesV} stokes_parameter_t;

/*!
 * \brief Convert a stokes_parameter_t value to the corresponding string name.
 * \details This function converts a stokes_parameter_t value to the corresponding string name.
 * \param StokesParameter A stokes_parameter_t value
 * \return A string
 */
string StokesName(stokes_parameter_t StokesParameter);

/*!
 * \author Andre Venne
 * \brief This class encapsulates the Stokes vector.
 * \sa class operaFluxVector, class operaMuellerMatrix, class operaPolarimetry
 * 
 * This class contains an array of 4 operaFluxVector to store each Stokes parameter or
 * any results with the same format, such as the degree of polarization or
 * the 2 null polarization spectrums.
 */
class operaStokesVector {
    
friend class operaMuellerMatrix;
    
private:

protected:
    unsigned length;
    bool istemp;
    operaFluxVector *stokesVector[4];	// The 4 Stokes parameters I, Q, U, V
	
public:
	/*
	 * Constructors / Destructors
	 */
    
    /*!
     * \brief Basic operaStokesVector constructor.
     * \param Istemp An optional bool defaults to false
     * \return void
     */
    operaStokesVector(bool Istemp=false);
    
    /*!
     * \brief operaStokesVector constructor.
     * \param Length An unsigned number of elements in each operaFluxVector
     * \param Istemp An optional bool defaults to false
     * \return void
     */
	operaStokesVector(unsigned Length, bool Istemp=false);
    
    /*!
     * \brief operaStokesVector constructor from an array of 4 operaFluxVector pointers.
     * \param StokesVector An array of 4 operaFluxVector pointers (StokesI, StokesQ, StokesU, StokesV)
     * \param Istemp An optional bool defaults to false
     * \return void
     */
    operaStokesVector(operaFluxVector *StokesVector[4], bool Istemp=false);
    
    /*!
     * \brief operaStokesVector constructor from 4 operaFluxVector pointers.
     * \param StokesVectorI Stokes I operaFluxVector pointer
     * \param StokesVectorQ Stokes Q operaFluxVector pointer
     * \param StokesVectorU Stokes U operaFluxVector pointer
     * \param StokesVectorV Stokes V operaFluxVector pointer
     * \param Istemp An optional bool defaults to false
     * \return void
     */
    operaStokesVector(operaFluxVector *StokesVectorI, operaFluxVector *StokesVectorQ, operaFluxVector *StokesVectorU, operaFluxVector *StokesVectorV, bool Istemp=false);
	
    /*!
     * \brief Destructor releases the pointers memory.
     * \return void
     */
	~operaStokesVector(); 
    
    /*
	 * Getters/Setters
	 */
    
    /*!
     * \brief Sets the length of the vectors.
     * \details A function that sets the value of the variable holding the number of elements in the operaFluxVector.
     * \param Length An unsigned number of elements
     * \return void
     */
    void setLength(unsigned Length);
	
    /*!
     * \brief Gets the length of the vectors.
     * \details A function that returns the value of the variable holding the number of elements in the operaFluxVector.
     * \return An unsigned value
     */
	unsigned getLength(void);
    
    /*!
     * \brief Sets Stokes parameters.
     * \details A function that sets the operaFluxVector of each Stokes parameter.
     * \param StokesVector An operaStokesVector address
     * \return void
     */
    void setStokesParameters(operaStokesVector &StokesVector);
    
    /*!
     * \brief Sets Stokes parameters.
     * \details A function that sets the operaFluxVector of each Stokes parameter.
     * \param StokesIFluxVector An operaFluxVector address
     * \param StokesQFluxVector An operaFluxVector address
     * \param StokesUFluxVector An operaFluxVector address
     * \param StokesVFluxVector An operaFluxVector address
     * \return void
     */
    void setStokesParameters(operaFluxVector &StokesIFluxVector, operaFluxVector &StokesQFluxVector, operaFluxVector &StokesUFluxVector, operaFluxVector &StokesVFluxVector);
    
    /*!
     * \brief Sets a Stokes parameter.
     * \details A function that sets the operaFluxVector of a Stokes parameter.
     * \param StokesIndex A stokes_parameter_t value
     * \param FluxVector An operaFluxVector pointer
     * \return void
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector);
    
    /*!
     * \brief Sets a Stokes parameter element.
     * \details A function that sets the value and the variance of a Stokes parameter element.
     * \param StokesIndex A stokes_parameter_t value
     * \param StokesValue A double value
     * \param Variance A double value
     * \param index An unsigned index to the element
     * \return void
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets a Stokes parameter.
     * \details A function that gets the operaFluxVector of a Stokes parameter.
     * \param StokesIndex A stokes_parameter_t value
     * \return An operaFluxVector pointer
     */
    operaFluxVector* getStokesParameter(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Gets a Stokes parameter value vector.
     * \details A function that gets the value vector of a Stokes parameter.
     * \param StokesIndex A stokes_parameter_t value
     * \return A double pointer
     */
    double* getStokesParameterFluxes(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Gets a Stokes parameter variance vector.
     * \details A function that gets the variance vector of a Stokes parameter.
     * \param StokesIndex A stokes_parameter_t value
     * \return A double pointer
     */
    double* getStokesParameterVariances(stokes_parameter_t StokesIndex);
	
    /*!
     * \brief Gets a Stokes parameter element value.
     * \details A function that gets the value of a Stokes parameter element.
     * \param StokesIndex A stokes_parameter_t value
     * \param index An unsigned index to the element
     * \return A double value
     */
    double getStokesParameterFlux(stokes_parameter_t StokesIndex, unsigned index);
    
    /*!
     * \brief Gets a Stokes parameter element variance.
     * \details A function that gets the variance of a Stokes parameter element.
     * \param StokesIndex A stokes_parameter_t value
     * \param index An unsigned index to the element
     * \return A double value
     */
    double getStokesParameterVariance(stokes_parameter_t StokesIndex, unsigned index);
	
    /*
	 * Operators
	 */
    
    /*! 
	 * \brief Indexing operator.
	 * \details The operator returns the Stokes parameter operaFluxVector.
	 * \param StokesIndex A stokes_parameter_t value
	 * \note Usage: operaFluxVector = StokesVector[StokesIndex];
	 * \return An operaFluxVector address
	 */
	operaFluxVector& operator[](stokes_parameter_t StokesIndex) {
#ifdef RANGE_CHECK
        if (StokesIndex < 0 || StokesIndex >= 4) {
            throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
        }
#endif
        return *(stokesVector[StokesIndex]);
    };
};

#endif

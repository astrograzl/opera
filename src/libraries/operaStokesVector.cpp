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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaStokesVector.h"

/*!
 * \file operaStokesVector.cpp
 * \brief This file holds the implementation of the class operaStokesVector.
 */

/*
 * \author Andre Venne
 * \brief This class encapsulates the Stokes vector.
 * \sa class operaFluxVector, class operaMuellerMatrix, class operaPolarimetry
 * 
 * This class contains an array of 4 operaFluxVector to store each Stokes parameter or
 * any results with the same format, such as the degree of polarization or
 * the 2 null polarization spectrums.
 */

/*
 * \brief Convert a stokes_parameter_t value to the corresponding string name.
 * \details This function converts a stokes_parameter_t value to the corresponding string name.
 * \param StokesParameter A stokes_parameter_t value
 * \return A string
 */
string StokesName(stokes_parameter_t StokesParameter)
{
    string StokesParameterName;
    
    switch (StokesParameter) {
        case StokesI:
            StokesParameterName = "Stokes I";
            break;
        case StokesQ:
            StokesParameterName = "Stokes Q";
            break;
        case StokesU:
            StokesParameterName = "Stokes U";
            break;
        case StokesV:
            StokesParameterName = "Stokes V";
            break;
        default:
            break;
    }
    
    return StokesParameterName;
}

/*
 * Constructors / Destructors
 */

/*
 * \brief Basic operaStokesVector constructor.
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaStokesVector::operaStokesVector(bool Istemp) :
length(0)
{
    istemp = Istemp;
	stokesVector[StokesI] = NULL;
	stokesVector[StokesQ] = NULL;
	stokesVector[StokesU] = NULL;
	stokesVector[StokesV] = NULL;
}

/*
 * \brief operaStokesVector constructor.
 * \param Length An unsigned number of elements in each operaFluxVector
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaStokesVector::operaStokesVector(unsigned Length, bool Istemp) :
length(0)
{
    length = Length;
    istemp = Istemp;
    stokesVector[StokesI] = new operaFluxVector(length);
    stokesVector[StokesQ] = new operaFluxVector(length);
    stokesVector[StokesU] = new operaFluxVector(length);
    stokesVector[StokesV] = new operaFluxVector(length);
}

/*
 * \brief operaStokesVector constructor from an array of 4 operaFluxVector pointers.
 * \param StokesVector An array of 4 operaFluxVector pointers (StokesI, StokesQ, StokesU, StokesV)
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaStokesVector::operaStokesVector(operaFluxVector *StokesVector[4], bool Istemp)
{
    length = StokesVector[StokesI]->getlength();
    istemp = Istemp;
    stokesVector[StokesI] = new operaFluxVector(*StokesVector[StokesI]);
    stokesVector[StokesQ] = new operaFluxVector(*StokesVector[StokesQ]);
    stokesVector[StokesU] = new operaFluxVector(*StokesVector[StokesU]);
    stokesVector[StokesV] = new operaFluxVector(*StokesVector[StokesV]);
}

/*
 * \brief operaStokesVector constructor from 4 operaFluxVector pointers.
 * \param StokesVectorI Stokes I operaFluxVector pointer
 * \param StokesVectorQ Stokes Q operaFluxVector pointer
 * \param StokesVectorU Stokes U operaFluxVector pointer
 * \param StokesVectorV Stokes V operaFluxVector pointer
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaStokesVector::operaStokesVector(operaFluxVector *StokesVectorI, operaFluxVector *StokesVectorQ, operaFluxVector *StokesVectorU, operaFluxVector *StokesVectorV, bool Istemp)
{
    length = StokesVectorI->getlength();
    istemp = Istemp;
    stokesVector[StokesI] = new operaFluxVector(*StokesVectorI);
    stokesVector[StokesQ] = new operaFluxVector(*StokesVectorQ);
    stokesVector[StokesU] = new operaFluxVector(*StokesVectorU);
    stokesVector[StokesV] = new operaFluxVector(*StokesVectorV);
}

/*
 * \brief Destructor releases the pointers memory.
 * \return void
 */
operaStokesVector::~operaStokesVector()
{
	if(stokesVector[StokesI])
        delete stokesVector[StokesI];
    stokesVector[StokesI] = NULL;
    if(stokesVector[StokesQ])
        delete stokesVector[StokesQ];
    stokesVector[StokesQ] = NULL;
    if(stokesVector[StokesU])
        delete stokesVector[StokesU];
    stokesVector[StokesU] = NULL;
    if(stokesVector[StokesV])
        delete stokesVector[StokesV];
    stokesVector[StokesV] = NULL;
}

/*
 * Getters/Setters
 */

/*
 * \brief Resizes the operaStokesVector.
 * \details Resizes each of the 4 operaFluxVectors.
 * \param Length The new length to resize to
 * \return void
 */
void operaStokesVector::resize(unsigned Length)
{
	stokesVector[StokesI]->resize(Length);
    stokesVector[StokesQ]->resize(Length);
    stokesVector[StokesU]->resize(Length);
    stokesVector[StokesV]->resize(Length);
    length = Length;
}

/*
 * \brief Sets the length of the vectors.
 * \details A function that sets the value of the variable holding the number of elements in the operaFluxVector.
 * \param Length An unsigned number of elements
 * \return void
 */
void operaStokesVector::setLength(unsigned Length)
{
    length = Length;
}

/*
 * \brief Gets the length of the vectors.
 * \details A function that returns the value of the variable holding the number of elements in the operaFluxVector.
 * \return An unsigned value
 */
unsigned operaStokesVector::getLength(void) const
{
    return length;
}

/*
 * \brief Sets Stokes parameters.
 * \details A function that sets the operaFluxVector of each Stokes parameter.
 * \param StokesVector An operaStokesVector address
 * \return void
 */
void operaStokesVector::setStokesParameters(operaStokesVector &StokesVector)
{
    stokesVector[StokesI]->setVector(*StokesVector.stokesVector[StokesI]);
    stokesVector[StokesQ]->setVector(*StokesVector.stokesVector[StokesQ]);
    stokesVector[StokesU]->setVector(*StokesVector.stokesVector[StokesU]);
    stokesVector[StokesV]->setVector(*StokesVector.stokesVector[StokesV]);
}

/*
 * \brief Sets Stokes parameters.
 * \details A function that sets the operaFluxVector of each Stokes parameter.
 * \param StokesIFluxVector An operaFluxVector address
 * \param StokesQFluxVector An operaFluxVector address
 * \param StokesUFluxVector An operaFluxVector address
 * \param StokesVFluxVector An operaFluxVector address
 * \return void
 */
void operaStokesVector::setStokesParameters(operaFluxVector &StokesIFluxVector, operaFluxVector &StokesQFluxVector, operaFluxVector &StokesUFluxVector, operaFluxVector &StokesVFluxVector)
{
    stokesVector[StokesI]->setVector(StokesIFluxVector);
    stokesVector[StokesQ]->setVector(StokesQFluxVector);
    stokesVector[StokesU]->setVector(StokesUFluxVector);
    stokesVector[StokesV]->setVector(StokesVFluxVector);
}

/*
 * \brief Sets a Stokes parameter.
 * \details A function that sets the operaFluxVector of a Stokes parameter.
 * \param StokesIndex A stokes_parameter_t value
 * \param FluxVector An operaFluxVector pointer
 * \return void
 */
void operaStokesVector::setStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector)
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    stokesVector[StokesIndex]->setVector(*FluxVector);
}

/*
 * \brief Sets a Stokes parameter element.
 * \details A function that sets the value and the variance of a Stokes parameter element.
 * \param StokesIndex A stokes_parameter_t value
 * \param StokesValue A double value
 * \param Variance A double value
 * \param index An unsigned index to the element
 * \return void
 */
void operaStokesVector::setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index)
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    stokesVector[StokesIndex]->setflux(StokesValue, index);
    stokesVector[StokesIndex]->setvariance(Variance, index);
}

/*
 * \brief Gets a Stokes parameter.
 * \details A function that gets the operaFluxVector of a Stokes parameter.
 * \param StokesIndex A stokes_parameter_t value
 * \return An operaFluxVector pointer
 */
operaFluxVector* operaStokesVector::getStokesParameter(stokes_parameter_t StokesIndex)
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex];
}

/*
 * \brief Gets a Stokes parameter value vector.
 * \details A function that gets the value vector of a Stokes parameter.
 * \param StokesIndex A stokes_parameter_t value
 * \return A double pointer
 */
double* operaStokesVector::getStokesParameterFluxes(stokes_parameter_t StokesIndex)
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex]->getfluxes();
}

/*
 * \brief Gets a Stokes parameter variance vector.
 * \details A function that gets the variance vector of a Stokes parameter.
 * \param StokesIndex A stokes_parameter_t value
 * \return A double pointer
 */
double* operaStokesVector::getStokesParameterVariances(stokes_parameter_t StokesIndex)
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex]->getvariances();
}

/*
 * \brief Gets a Stokes parameter element value.
 * \details A function that gets the value of a Stokes parameter element.
 * \param StokesIndex A stokes_parameter_t value
 * \param index An unsigned index to the element
 * \return A double value
 */
double operaStokesVector::getStokesParameterFlux(stokes_parameter_t StokesIndex, unsigned index) const
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex]->getflux(index);
}

/*
 * \brief Gets a Stokes parameter element variance.
 * \details A function that gets the variance of a Stokes parameter element.
 * \param StokesIndex A stokes_parameter_t value
 * \param index An unsigned index to the element
 * \return A double value
 */
double operaStokesVector::getStokesParameterVariance(stokes_parameter_t StokesIndex, unsigned index) const
{
#ifdef RANGE_CHECK
    if (StokesIndex < 0 || StokesIndex >= 4) {
        throw operaException("operaStokesVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector[StokesIndex]->getvariance(index);
}

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaPolarimetry
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

#include <iostream>
#include <string>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <utility>	// for pair

#include "globaldefines.h"
#include "operaError.h"
#include "libraries/operaException.h"
#include "libraries/operaSpectralOrder.h"
#include "libraries/operaFluxVector.h"
#include "libraries/operaStokesVector.h"
#include "libraries/operaPolarimetry.h"

#include "libraries/operaFit.h"
#include "libraries/operaMath.h"

/*!
 * \file operaPolarimetry.cpp
 */

/*
 * \author Andre Venne
 * \brief This class encapsulates the polarimetry results.
 * \sa class operaStokesVector
 * 
 * This class holds in operaStokesVector classes the 4 Stokes parameters, the 4 associated degrees of polarization
 * and the 2 null polarization spectra for each Stokes parameter.
 */

/*
 * Constructors / Destructors
 */

/*
 * \brief Basic operaPolarimetry constructor.
 * \return void
 */
operaPolarimetry::operaPolarimetry():
length(0),
hasStokesI(false),
hasStokesQ(false),
hasStokesU(false),
hasStokesV(false),
hasDegreeOfStokesI(false),
hasDegreeOfStokesQ(false),
hasDegreeOfStokesU(false),
hasDegreeOfStokesV(false),
hasFirstNullPolarization(false),
hasSecondNullPolarization(false),
hasWavelength(false)
{
    stokesVector = NULL;
    degreeOfPolarization = NULL;
    firstNullPolarization = NULL;
    secondNullPolarization = NULL;  
	wavelength = NULL;	   
}

/*
 * \brief Basic operaPolarimetry constructor.
 * \param Length An unsigned number of elements in each operaStokesVector
 * \return void
 */
operaPolarimetry::operaPolarimetry(unsigned Length):
hasStokesI(false),
hasStokesQ(false),
hasStokesU(false),
hasStokesV(false),
hasDegreeOfStokesI(false),
hasDegreeOfStokesQ(false),
hasDegreeOfStokesU(false),
hasDegreeOfStokesV(false),
hasFirstNullPolarization(false),
hasSecondNullPolarization(false),
hasWavelength(false)
{
	if (Length == 0) {
		throw operaException("operaPolarimetry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
    length = Length;
    stokesVector = new operaStokesVector(length);
    degreeOfPolarization = new operaStokesVector(length);
    firstNullPolarization = new operaStokesVector(length);
    secondNullPolarization = new operaStokesVector(length); 
	wavelength = (double *)malloc(length*sizeof(double));	   
}


/*
 * \brief Basic operaPolarimetry destructor.
 * \return void
 */
operaPolarimetry::~operaPolarimetry()
{
    if (stokesVector)
	   delete stokesVector;
    stokesVector = NULL;
    if (degreeOfPolarization)
        delete degreeOfPolarization;
    degreeOfPolarization = NULL;
    if (firstNullPolarization)
        delete firstNullPolarization;
    firstNullPolarization = NULL;
    if (secondNullPolarization)
        delete secondNullPolarization;
    secondNullPolarization = NULL;
	if (wavelength && hasWavelength)
		free(wavelength);
	wavelength = NULL;
}

/*
 * Getters/Setters
 */


void operaPolarimetry::resize(unsigned Length)
{
	if (Length == 0) {
		throw operaException("operaPolarimetry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (length == 0) {
		throw operaException("operaPolarimetry: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if(Length == length) return;
	stokesVector->resize(Length);
    degreeOfPolarization->resize(Length);
    firstNullPolarization->resize(Length);
    secondNullPolarization->resize(Length);
    resizeVector(wavelength, length, Length);
	length = Length;
}

/*
 * \brief Sets the length of the Stokes vectors.
 * \details A function that sets the value of the variable holding the number of elements in the operaStokesVector.
 * \param Length An unsigned number of elements
 * \return void
 */
void operaPolarimetry::setLength(unsigned Length)
{
	if (Length > length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
    length = Length;
}

/*
 * \brief Gets the length of the Stokes vectors.
 * \details A function that gets the value of the variable holding the number of elements in the operaStokesVector.
 * \return An unsigned value
 */
unsigned operaPolarimetry::getLength(void) const
{
    return length;
}

/*
 * \brief Gets a wavelength value at indexElem.
 * \param indexEleme the index
 * \return Wavelength a wavelength value
 */
double operaPolarimetry::getwavelength(unsigned indexElem) const {
#ifdef RANGE_CHECK
    if (indexElem >= length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return wavelength[indexElem];
}
/*
 * \brief Gets a wavelength vector address.
 * \return Wavelength a wavelength vector address
 */
double *operaPolarimetry::getwavelength(void) {
	return wavelength;
}
/*
 * \brief Sets a wavelength value at indexElem.
 * \param Wavelength a wavelength value
 * \param indexEleme the index
 * \return void
 */
void operaPolarimetry::setwavelength(double Wavelength, unsigned indexElem) {
#ifdef RANGE_CHECK
    if (indexElem >= length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	wavelength[indexElem] = Wavelength;
}

/*
 * \brief Copies wavelength values.
 * \param double *Wavelength a wavelength vector
 * \param length the number of elements
 * \return void
 */
void operaPolarimetry::setwavelength(double *Wavelength, unsigned Length) {
#ifdef RANGE_CHECK
    if (Length >= length) {
		throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	length = Length;
	unsigned i = length;
	while (i--) {
		wavelength[i] = Wavelength[i];
	}
	setHasWavelength(true);
}

/*
 * \brief Sets a Stokes parameter.
 * \details A function that sets the operaFluxVector of a Stokes parameter.
 * \param StokesIndex A stokes_parameter_t value
 * \param FluxVector An operaFluxVector pointer
 * \return void
 */
void operaPolarimetry::setStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    stokesVector->setStokesParameter(StokesIndex, FluxVector);
    
    setHasStokes(StokesIndex,true);
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
void operaPolarimetry::setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    stokesVector->getStokesParameter(StokesIndex)->setflux(StokesValue, index);
    stokesVector->getStokesParameter(StokesIndex)->setvariance(Variance, index);
    
    switch (StokesIndex) {
        case StokesI:
            hasStokesI = true;
            break;
        case StokesQ:
            hasStokesQ = true;
            break;
        case StokesU:
            hasStokesU = true;
            break;
        case StokesV:
            hasStokesV = true;
            break;
    }
}

/*
 * \brief Gets the Stokes vector.
 * \details A function that gets the operaStokesVector of the Stokes vector.
 * \return An operaStokesVector pointer
 */
operaStokesVector* operaPolarimetry::getStokesVector(void)
{
    return stokesVector;
}
const operaStokesVector* operaPolarimetry::getStokesVector(void) const
{
    return stokesVector;
}

/*
 * \brief Gets a Stokes parameter.
 * \details A function that gets the operaFluxVector of a Stokes parameter.
 * \param StokesIndex A stokes_parameter_t value
 * \return An operaFluxVector pointer
 */
operaFluxVector* operaPolarimetry::getStokesParameter(stokes_parameter_t StokesIndex)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector->getStokesParameter(StokesIndex);
}
const operaFluxVector* operaPolarimetry::getStokesParameter(stokes_parameter_t StokesIndex) const
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return stokesVector->getStokesParameter(StokesIndex);
}

/*
 * \brief Sets a Stokes parameter of the degree of polarization vector.
 * \details A function that sets the operaFluxVector of a Stokes parameter of the degree of polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \param FluxVector An operaFluxVector pointer
 * \return void
 */
void operaPolarimetry::setDegreeOfPolarization(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    degreeOfPolarization->setStokesParameter(StokesIndex, FluxVector);
    setHasDegreeOfStokes(StokesIndex,true);
}

/*
 * \brief Sets a Stokes parameter element of the degree of polarization vector.
 * \details A function that sets the value and the variance of a Stokes parameter element of the degree of polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \param DegreeOfPolarizationValue A double value
 * \param Variance A double value
 * \param index An unsigned index to the element
 * \return void
 */
void operaPolarimetry::setDegreeOfPolarization(stokes_parameter_t StokesIndex, double DegreeOfPolarizationValue, double Variance, unsigned index)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    degreeOfPolarization->getStokesParameter(StokesIndex)->setflux(DegreeOfPolarizationValue, index);
    degreeOfPolarization->getStokesParameter(StokesIndex)->setvariance(Variance, index);
    
    switch (StokesIndex) {
        case StokesI:
            hasDegreeOfStokesI = true;
            break;
        case StokesQ:
            hasDegreeOfStokesQ = true;
            break;
        case StokesU:
            hasDegreeOfStokesU = true;
            break;
        case StokesV:
            hasDegreeOfStokesV = true;
            break;
    }
}

/*
 * \brief Gets the degree of polarization vector.
 * \details A function that gets the operaStokesVector of the degree of polarization vector.
 * \return An operaStokesVector pointer
 */
operaStokesVector* operaPolarimetry::getDegreeOfPolarization(void)
{
    return degreeOfPolarization;
}
const operaStokesVector* operaPolarimetry::getDegreeOfPolarization(void) const
{
    return degreeOfPolarization;
}

/*
 * \brief Gets a Stokes parameter of the degree of polarization vector.
 * \details A function that gets the operaFluxVector of a Stokes parameter of the degree of polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \return An operaFluxVector pointer
 */
operaFluxVector* operaPolarimetry::getDegreeOfPolarization(stokes_parameter_t StokesIndex)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return degreeOfPolarization->getStokesParameter(StokesIndex);
}
const operaFluxVector* operaPolarimetry::getDegreeOfPolarization(stokes_parameter_t StokesIndex) const
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return degreeOfPolarization->getStokesParameter(StokesIndex);
}

/*
 * \brief Sets a Stokes parameter of the first null polarization vector.
 * \details A function that sets the operaFluxVector of a Stokes parameter of the first null polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \param FluxVector An operaFluxVector pointer
 * \return void
 */
void operaPolarimetry::setFirstNullPolarization(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    firstNullPolarization->setStokesParameter(StokesIndex, FluxVector);
    
    hasFirstNullPolarization = true;
}

/*
 * \brief Sets a Stokes parameter element of the first null polarization vector.
 * \details A function that sets the value and the variance of a Stokes parameter element of the first null polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \param FirstNullPolarizationValue A double value
 * \param Variance A double value
 * \param index An unsigned index to the element
 * \return void
 */
void operaPolarimetry::setFirstNullPolarization(stokes_parameter_t StokesIndex, double FirstNullPolarizationValue, double Variance, unsigned index)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    firstNullPolarization->getStokesParameter(StokesIndex)->setflux(FirstNullPolarizationValue, index);
    firstNullPolarization->getStokesParameter(StokesIndex)->setvariance(Variance, index);
    
    hasFirstNullPolarization = true;
}

/*
 * \brief Gets the first null polarization vector.
 * \details A function that gets the operaStokesVector of the first null polarization vector.
 * \return An operaStokesVector pointer
 */
operaStokesVector* operaPolarimetry::getFirstNullPolarization(void)
{
    return firstNullPolarization;
}
const operaStokesVector* operaPolarimetry::getFirstNullPolarization(void) const
{
    return firstNullPolarization;
}

/*
 * \brief Gets a Stokes parameter of the first null polarization vector.
 * \details A function that gets the operaFluxVector of a Stokes parameter of the first null polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \return An operaFluxVector pointer
 */
operaFluxVector* operaPolarimetry::getFirstNullPolarization(stokes_parameter_t StokesIndex)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return firstNullPolarization->getStokesParameter(StokesIndex);
}
const operaFluxVector* operaPolarimetry::getFirstNullPolarization(stokes_parameter_t StokesIndex) const
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return firstNullPolarization->getStokesParameter(StokesIndex);
}

/*
 * \brief Sets a Stokes parameter of the second null polarization vector.
 * \details A function that sets the operaFluxVector of a Stokes parameter of the second null polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \param FluxVector An operaFluxVector pointer
 * \return void
 */
void operaPolarimetry::setSecondNullPolarization(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    secondNullPolarization->setStokesParameter(StokesIndex, FluxVector);
    
    hasSecondNullPolarization = true;
}

/*
 * \brief Sets a Stokes parameter element of the second null polarization vector.
 * \details A function that sets the value and the variance of a Stokes parameter element of the second null polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \param SecondNullPolarizationValue A double value
 * \param Variance A double value
 * \param index An unsigned index to the element
 * \return void
 */
void operaPolarimetry::setSecondNullPolarization(stokes_parameter_t StokesIndex, double SecondNullPolarizationValue, double Variance, unsigned index)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    secondNullPolarization->getStokesParameter(StokesIndex)->setflux(SecondNullPolarizationValue, index);
    secondNullPolarization->getStokesParameter(StokesIndex)->setvariance(Variance, index);
    
    hasSecondNullPolarization = true;
}

/*
 * \brief Gets the second null polarization vector.
 * \details A function that gets the operaStokesVector of the second null polarization vector.
 * \return An operaStokesVector pointer
 */
operaStokesVector* operaPolarimetry::getSecondNullPolarization(void)
{
    return secondNullPolarization;
}
const operaStokesVector* operaPolarimetry::getSecondNullPolarization(void) const
{
    return secondNullPolarization;
}

/*
 * \brief Gets a Stokes parameter of the second null polarization vector.
 * \details A function that gets the operaFluxVector of a Stokes parameter of the second null polarization vector.
 * \param StokesIndex A stokes_parameter_t value
 * \return An operaFluxVector pointer
 */
operaFluxVector* operaPolarimetry::getSecondNullPolarization(stokes_parameter_t StokesIndex)
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return secondNullPolarization->getStokesParameter(StokesIndex);
}
const operaFluxVector* operaPolarimetry::getSecondNullPolarization(stokes_parameter_t StokesIndex) const
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);    
    }
#endif
    return secondNullPolarization->getStokesParameter(StokesIndex);
}

/*
 * \brief Gets the boolean value of Stokes I.
 * \details A function that gets the boolean value giving the state of the Stokes I parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasStokesI(void) const
{
    return hasStokesI;
}

/*
 * \brief Gets the boolean value of Stokes Q.
 * \details A function that gets the boolean value giving the state of the Stokes Q parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasStokesQ(void) const
{
    return hasStokesQ;
}

/*
 * \brief Gets the boolean value of Stokes U.
 * \details A function that gets the boolean value giving the state of the Stokes U parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasStokesU(void) const
{
    return hasStokesU;
}

/*
 * \brief Gets the boolean value of Stokes V.
 * \details A function that gets the boolean value giving the state of the Stokes V parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasStokesV(void) const
{
    return hasStokesV;
}

/*!
 * \brief Get the boolean value of a given Stokes.
 * \details A function that gets the boolean value giving the state of a given Stokes parameter in the Stokes vector.
 * \param StokesIndex choice of Stokes parameter
 * \return A boolean value
 */
bool operaPolarimetry::getHasStokes(stokes_parameter_t StokesIndex) const
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
#endif
    
    switch (StokesIndex) {
        case StokesI:
            return hasStokesI;
            break;
        case StokesQ:
            return hasStokesQ;
            break;
        case StokesU:
            return hasStokesU;
            break;
        case StokesV:
            return hasStokesV;
            break;
    }
	return false;
}


/*
 * \brief Sets the boolean value of Stokes I.
 * \details A function that sets the boolean value giving the state of the Stokes I parameter in the Stokes vector.
 * \param HasStokesI A boolean value
 * \return void
 */
void operaPolarimetry::setHasStokesI(bool HasStokesI)
{
    hasStokesI = HasStokesI;
}

/*
 * \brief Sets the boolean value of Stokes Q.
 * \details A function that sets the boolean value giving the state of the Stokes Q parameter in the Stokes vector.
 * \param HasStokesQ A boolean value
 * \return void
 */
void operaPolarimetry::setHasStokesQ(bool HasStokesQ)
{
    hasStokesQ = HasStokesQ;
}

/*
 * \brief Sets the boolean value of Stokes U.
 * \details A function that sets the boolean value giving the state of the Stokes U parameter in the Stokes vector.
 * \param HasStokesU A boolean value
 * \return void
 */
void operaPolarimetry::setHasStokesU(bool HasStokesU)
{
    hasStokesU = HasStokesU;
}

/*
 * \brief Sets the boolean value of Stokes V.
 * \details A function that sets the boolean value giving the state of the Stokes V parameter in the Stokes vector.
 * \param HasStokesV A boolean value
 * \return void
 */
void operaPolarimetry::setHasStokesV(bool HasStokesV)
{
    hasStokesV = HasStokesV;
}

/*!
 * \brief Set the boolean value of a given Stokes.
 * \details A function that sets the boolean value giving the state of a given Stokes parameter in the Stokes vector.
 * \param HasStokes A boolean value
 * \param StokesIndex choice of Stokes parameter
 * \return void
 */
void operaPolarimetry::setHasStokes(stokes_parameter_t StokesIndex, bool HasStokes) {
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
#endif
    
    switch (StokesIndex) {
        case StokesI:
            hasStokesI = HasStokes;
            break;
        case StokesQ:
            hasStokesQ = HasStokes;
            break;
        case StokesU:
            hasStokesU = HasStokes;
            break;
        case StokesV:
            hasStokesV = HasStokes;
            break;
    }
}

/*
 * \brief Gets the boolean value of the degree of polarization of Stokes I.
 * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes I parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasDegreeOfStokesI(void) const
{
    return hasDegreeOfStokesI;
}

/*
 * \brief Gets the boolean value of the degree of polarization of Stokes Q.
 * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes Q parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasDegreeOfStokesQ(void) const
{
    return hasDegreeOfStokesQ;
}

/*
 * \brief Gets the boolean value of the degree of polarization of Stokes U.
 * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes U parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasDegreeOfStokesU(void) const
{
    return hasDegreeOfStokesU;
}

/*
 * \brief Gets the boolean value of the degree of polarization of Stokes V.
 * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes V parameter in the Stokes vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasDegreeOfStokesV(void) const
{
    return hasDegreeOfStokesV;
}

/*!
 * \brief Get the boolean value of the degree of polarization of a given Stokes.
 * \details A function that gets the boolean value giving the state of the degree of polarization of a given Stokes parameter in the Stokes vector.
 * \param StokesIndex choice of Stokes parameter
 * \return A boolean value
 */
bool operaPolarimetry::getHasDegreeOfStokes(stokes_parameter_t StokesIndex) const
{
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
#endif
    
    switch (StokesIndex) {
        case StokesI:
            return hasDegreeOfStokesI;
            break;
        case StokesQ:
            return hasDegreeOfStokesQ;
            break;
        case StokesU:
            return hasDegreeOfStokesU;
            break;
        case StokesV:
            return hasDegreeOfStokesV;
            break;
    }
	return false;
}

/*
 * \brief Sets the boolean value of the degree of polarization of Stokes I.
 * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes I parameter in the Stokes vector.
 * \param HasDegreeOfStokesI A boolean value
 * \return void
 */
void operaPolarimetry::setHasDegreeOfStokesI(bool HasDegreeOfStokesI)
{
    hasDegreeOfStokesI = HasDegreeOfStokesI;
}

/*
 * \brief Sets the boolean value of the degree of polarization of Stokes Q.
 * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes Q parameter in the Stokes vector.
 * \param HasDegreeOfStokesQ A boolean value
 * \return void
 */
void operaPolarimetry::setHasDegreeOfStokesQ(bool HasDegreeOfStokesQ)
{
    hasDegreeOfStokesQ = HasDegreeOfStokesQ;
}

/*
 * \brief Sets the boolean value of the degree of polarization of Stokes U.
 * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes U parameter in the Stokes vector.
 * \param HasDegreeOfStokesU A boolean value
 * \return void
 */
void operaPolarimetry::setHasDegreeOfStokesU(bool HasDegreeOfStokesU)
{
    hasDegreeOfStokesU = HasDegreeOfStokesU;
}

/*
 * \brief Sets the boolean value of the degree of polarization of Stokes V.
 * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes V parameter in the Stokes vector.
 * \param HasDegreeOfStokesV A boolean value
 * \return void
 */
void operaPolarimetry::setHasDegreeOfStokesV(bool HasDegreeOfStokesV)
{
    hasDegreeOfStokesV = HasDegreeOfStokesV;
}

/*!
 * \brief Set the boolean value of the degree of polarization of a given Stokes.
 * \details A function that sets the boolean value giving the state of the degree of polarization of a given Stokes parameter in the Stokes vector.
 * \param HasDegreeOfStokes A boolean value
 * \param StokesIndex choice of Stokes parameter
 * \return void
 */
void operaPolarimetry::setHasDegreeOfStokes(stokes_parameter_t StokesIndex, bool HasDegreeOfStokes) {
#ifdef RANGE_CHECK
    if (StokesIndex >= 4) {
        throw operaException("operaPolarimetry: unrecognized Stokes parameter. ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
#endif
    
    switch (StokesIndex) {
        case StokesI:
            hasDegreeOfStokesI = HasDegreeOfStokes;
            break;
        case StokesQ:
            hasDegreeOfStokesQ = HasDegreeOfStokes;
            break;
        case StokesU:
            hasDegreeOfStokesU = HasDegreeOfStokes;
            break;
        case StokesV:
            hasDegreeOfStokesV = HasDegreeOfStokes;
            break;
    }
}

/*
 * \brief Gets the boolean value of the first null polarization vector.
 * \details A function that gets the boolean value of the first null polarization vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasFirstNullPolarization(void) const
{
    return hasFirstNullPolarization;
}

/*
 * \brief Gets the boolean value of the second null polarization vector.
 * \details A function that gets the boolean value of the second null polarization vector.
 * \return A boolean value
 */
bool operaPolarimetry::getHasSecondNullPolarization(void) const
{
    return hasSecondNullPolarization;
}

/*
 * \brief Sets the boolean value of the first null polarization vector.
 * \details A function that sets the boolean value of the first null polarization vector.
 * \param HasFirstNullPolarization A boolean value
 * \return void
 */
void operaPolarimetry::setHasFirstNullPolarization(bool HasFirstNullPolarization)
{
    hasFirstNullPolarization = HasFirstNullPolarization;
}

/*
 * \brief Sets the boolean value of the second null polarization vector.
 * \details A function that sets the boolean value of the second null polarization vector.
 * \param HasSecondNullPolarization A boolean value
 * \return void
 */
void operaPolarimetry::setHasSecondNullPolarization(bool HasSecondNullPolarization)
{
    hasSecondNullPolarization = HasSecondNullPolarization;
}

/*
 * Methods
 */

/*
 * \brief Calculates the polarization.
 * \details A function that calculates the polarization from the degree of polarization and stores it in the Stokes vector.
 * \return void
 */
void operaPolarimetry::calculatePolarization(void)
{
    if (hasStokesI) {
        if (hasDegreeOfStokesQ) {
			operaFluxVector fluxdiv = *degreeOfPolarization->getStokesParameter(StokesQ) * *stokesVector->getStokesParameter(StokesI);
            stokesVector->setStokesParameter(StokesQ, &fluxdiv);
            hasStokesQ = true;
        }
        
        if (hasDegreeOfStokesU) {
			operaFluxVector fluxdiv = *degreeOfPolarization->getStokesParameter(StokesU) * *stokesVector->getStokesParameter(StokesI);
            stokesVector->setStokesParameter(StokesU, &fluxdiv);
            hasStokesU = true;
        }
        
        if (hasDegreeOfStokesV) {
			operaFluxVector fluxdiv = *(degreeOfPolarization->getStokesParameter(StokesV)) * *(stokesVector->getStokesParameter(StokesI));
            stokesVector->setStokesParameter(StokesV, &fluxdiv);
            hasStokesV = true;
        }
    } else {
        throw operaException("operaPolarimetry: missing Stokes I.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
}

/*
 * \brief Calculates the polarization.
 * \details A function that calculates the polarization from the degree of polarization and stores it in the Stokes vector.
 * \return void
 */
void operaPolarimetry::calculatePolarization(stokes_parameter_t StokesIndex)
{
    if (hasStokesI) {
        if (getHasDegreeOfStokes(StokesIndex)) {	// StokesI will not be done...
			operaFluxVector fluxdiv = *(degreeOfPolarization->getStokesParameter(StokesIndex)) / *(stokesVector->getStokesParameter(StokesI));
            stokesVector->setStokesParameter(StokesIndex, &fluxdiv);
            setHasStokes(StokesIndex,true);
		}
    } else {
        throw operaException("operaPolarimetry: missing Stokes I.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
}


/*
 * \brief Calculates the degree of polarization.
 * \details A function that calculates the degree of polarization from the polarization and stores it in the degree of polarization vector.
 * \return void
 */
void operaPolarimetry::calculateDegreeOfPolarization(void) {
    if (hasStokesI) {
		operaFluxVector fluxdiv = *stokesVector->getStokesParameter(StokesI) / *stokesVector->getStokesParameter(StokesI);
        degreeOfPolarization->setStokesParameter(StokesI, &fluxdiv);
        hasDegreeOfStokesI = true;
        
        if (hasStokesQ) {
			fluxdiv = *stokesVector->getStokesParameter(StokesQ) / *stokesVector->getStokesParameter(StokesI);
            degreeOfPolarization->setStokesParameter(StokesQ, &fluxdiv);
            hasDegreeOfStokesQ = true;
        }
        
        if (hasStokesU) {
			fluxdiv = *stokesVector->getStokesParameter(StokesU) / *stokesVector->getStokesParameter(StokesI);
            degreeOfPolarization->setStokesParameter(StokesU, &fluxdiv);
            hasDegreeOfStokesU = true;
        }
        
        if (hasStokesV) {
			fluxdiv = *stokesVector->getStokesParameter(StokesV) / *stokesVector->getStokesParameter(StokesI);
            degreeOfPolarization->setStokesParameter(StokesV, &fluxdiv);
            hasDegreeOfStokesV = true;
        }
    } else {
        throw operaException("operaPolarimetry: missing Stokes I.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
}

/*
 * \brief Calculates the degree of polarization.
 * \details A function that calculates the degree of polarization from the observed ordinary and extraordinary beam fluxes.
 * \details This function accepts 2 or 4 input pairs of fluxes (polarimetric exposures).
 * \param StokesIndex is the selected Stokes parameter
 * \param *iE[4] is a set of flux vectors for the input beams with a given state of polarization (ordinary beams)
 * \param *iA[4] is a set of flux vectors for the input beams with a given orthogonal state of polarization (extra-ordinary beams)
 * \param NumberOfExposures is the number of input exposures (accepts only 2 or 4)
 * \return void
 */
void operaPolarimetry::calculateDegreeOfPolarization(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures) {
#ifdef DOUG
    if(NumberOfExposures != 2 && NumberOfExposures != 4) {
        throw operaException("operaPolarimetry: NumberOfExposures not valid.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double PairOfExposures = (double)NumberOfExposures / 2.0;
    operaFluxVector PoverI(length);
    operaFluxVector N1(length),N2(length);
    
    if (StokesIndex != StokesQ && StokesIndex != StokesU && StokesIndex != StokesV) {
        throw operaException("operaPolarimetry: unrecognized Stokes parameter.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }

    if(method == Difference || method == DifferenceWithBeamSwapped) {
        operaFluxVector G1(length), G2(length), G3(length), G4(length);
        operaFluxVector D1(length);
        
		if(method == Difference) {
		
			G1 = (*(iE[0]) - *(iA[0]))/(*(iE[0]) + *(iA[0]));
			G2 = (*(iE[1]) - *(iA[1]))/(*(iE[1]) + *(iA[1]));
			G3 = (*(iE[2]) - *(iA[2]))/(*(iE[2]) + *(iA[2]));
			G4 = (*(iE[3]) - *(iA[3]))/(*(iE[3]) + *(iA[3]));
		
		} else if (method == DifferenceWithBeamSwapped) {
			
			for(unsigned i=0;i<NumberOfExposures;i++) {
				switch (i) {
					case 0:
						G1 = *(iE[i]) - *(iE[i+1]) / *(iE[i]) + *(iE[i+1]);
						break;
					case 2:
						G3 = *(iE[i]) - *(iE[i+1]) / *(iE[i]) + *(iE[i+1]);
						break;
					case 1:
						G2 = *(iA[i-1]) - *(iA[i]) / *(iA[i-1]) + *(iA[i]);
						break;
					case 3:
						G4 = *(iA[i-1]) - *(iA[i]) / *(iA[i-1]) + *(iA[i]);
						break;
					default:
						break;
				}
			}
		}
         D1 = G1 - G2;
        
        if(NumberOfExposures==2) {
        
            PoverI = D1 / (2.0 * PairOfExposures);
        
        } else if (NumberOfExposures==4) {
            
            operaFluxVector D2(length);
            operaFluxVector D1s(length), D2s(length);
            
            D2 = G3 - G4;
            
            PoverI = (D1 + D2) / (2.0 * PairOfExposures);
            
            D1s = G1 - G4;
            D2s = G3 - G2;
            
            N1 = (D1 - D2) / (2.0 * PairOfExposures);
            N2 = (D1s - D2s) / (2.0 * PairOfExposures);
        }

    } else if (method == Ratio) {
        operaFluxVector r1(length), r2(length), r3(length), r4(length);
        operaFluxVector R1(length);
        operaFluxVector R(length);
        
		r1 = *(iE[0]) / *(iA[0]);
		r2 = *(iE[1]) / *(iA[1]);
		r3 = *(iE[2]) / *(iA[2]);
		r4 = *(iE[3]) / *(iA[3]);
    
        R1 = r1 / r2;
        
        if(NumberOfExposures==2) {
            
            R = Pow( R1 , 1.0/(2.0 * PairOfExposures) );
            PoverI = (R - 1.0) / (R + 1.0);
            
        } else if (NumberOfExposures==4) {
            
            operaFluxVector R2(length);
            operaFluxVector R1s(length), R2s(length);

            operaFluxVector RN1(length), RN2(length);
            
            R2 = r3 / r4;
            
            R1s = r1 / r4;
            R2s = r3 / r2;
            
            R = Pow( R1 * R2, 1.0 / (2.0 * PairOfExposures) );
            
            PoverI = (R - 1.0) / (R + 1.0);
            
            RN1 = Pow( R1 / R2, 1.0 / (2.0 * PairOfExposures) );
            N1 = (RN1 - 1.0) / (RN1 + 1.0);
            
            RN2 = Pow( R1s / R2s, 1.0 / (2.0 * PairOfExposures) );
            N2 = (RN2 - 1.0) / (RN2 + 1.0);
        }
    } else if (method == NewMethod) {
        PoverI = 0.0;
        N1 = 0.0;
        N2 = 0.0;
    } else {
        throw operaException("operaPolarimetry: unrecognized method to calculate polarization.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);        
    }
#else
    if(NumberOfExposures != 2 && NumberOfExposures != 4) {
        throw operaException("operaPolarimetry: NumberOfExposures not valid.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    double PairOfExposures = (double)NumberOfExposures / 2.0;
    operaFluxVector PoverI(length);
    operaFluxVector N1(length),N2(length);
    
    if (StokesIndex != StokesQ && StokesIndex != StokesU && StokesIndex != StokesV) {
        throw operaException("operaPolarimetry: unrecognized Stokes parameter.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
	
    if(method == Difference || method == DifferenceWithBeamSwapped) {
        operaFluxVector *G[4];
        operaFluxVector D1(length);
        
        for(unsigned i=0;i<NumberOfExposures;i++) {
            G[i] = new operaFluxVector(length);
            
            if(method == Difference) {
				
                *G[i] = (*iE[i] - *iA[i])/(*iE[i] + *iA[i]);
				
            } else if (method == DifferenceWithBeamSwapped) {
                
                if(i==0 || i==2) {
                    *G[i] = (*iE[i] - *iE[i+1])/(*iE[i] + *iE[i+1]);
                } else if (i==1 || i==3) {
                    *G[i] = (*iA[i-1] - *iA[i])/(*iA[i-1] + *iA[i]);
                }
                
            }
        }
		
        D1 = (*G[0]) - (*G[1]);
        
        if(NumberOfExposures==2) {
			
            PoverI = D1 / (2.0 * PairOfExposures);
			
        } else if (NumberOfExposures==4) {
            
            operaFluxVector D2(length);
            operaFluxVector D1s(length), D2s(length);
            
            D2 = (*G[2]) - (*G[3]);
            
            PoverI = (D1 + D2) / (2.0 * PairOfExposures);
            
            D1s = (*G[0]) - (*G[3]);
            D2s = (*G[2]) - (*G[1]);
            
            N1 = (D1 - D2) / (2.0 * PairOfExposures);
            N2 = (D1s - D2s) / (2.0 * PairOfExposures);
        }
		
        for(unsigned i=0;i<NumberOfExposures;i++) {
            delete G[i];
        }
    } else if (method == Ratio) {
        operaFluxVector *r[4];
        operaFluxVector R1(length);
        operaFluxVector R(length);
        
        for(unsigned i=0;i<NumberOfExposures;i++) {
            r[i] = new operaFluxVector(length);
            *r[i] = (*iE[i]) / (*iA[i]);
        }
		
        R1 = (*r[0]) / (*r[1]);
        
        if(NumberOfExposures==2) {
            
            R = Pow( R1 , 1.0/(2.0 * PairOfExposures) );
            PoverI = (R - 1.0) / (R + 1.0);
            
        } else if (NumberOfExposures==4) {
            
            operaFluxVector R2(length);
            operaFluxVector R1s(length),R2s(length);
			
            operaFluxVector RN1(length),RN2(length);
            
            R2 = (*r[2]) / (*r[3]);	// changed from r[0] / r[1] Apr 15 2013 DT
            
            R1s = (*r[0]) / (*r[3]);
            R2s = (*r[2]) / (*r[1]);
            
            R = Pow( R1 * R2 , 1.0/(2.0 * PairOfExposures) );
            
            PoverI = (R - 1.0) / (R + 1.0);
            
            RN1 = Pow( R1 / R2 , 1.0/(2.0 * PairOfExposures) );
            N1 = (RN1 - 1.0) / (RN1 + 1.0);
            
            RN2 = Pow( R1s / R2s , 1.0/(2.0 * PairOfExposures) );
            N2 = (RN2 - 1.0) / (RN2 + 1.0);
        }
        for(unsigned i=0;i<NumberOfExposures;i++) {
            delete r[i];
        }
    } else if (method == NewMethod) {
        PoverI = 0.0;
        N1 = 0.0;
        N2 = 0.0;
    } else {
        throw operaException("operaPolarimetry: unrecognized method to calculate polarization.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);        
    }
#endif
    
    // The fix below is necessary since the angles of the analyzer with respect to the
    // reference system for ESPaDOnS don't follow the same order as described in the literature.
    // added on Dec 04 2013 EM
    if(StokesIndex == StokesQ) {
        PoverI -= (PoverI * 2.0);
    }
    
    setDegreeOfPolarization(StokesIndex, &PoverI);
    
    if (NumberOfExposures==4) {
        setFirstNullPolarization(StokesIndex, &N1);
        setSecondNullPolarization(StokesIndex, &N2);
    }
}

/*
 * \brief Calculates Stokes I and another given Stokes parameter.
 * \details A function that calculates the polarized flux for a given Stokes and the
 * \details total flux for Stokes I given the observed ordinary and extra-ordinary beam fluxes.
 * \details This function accepts either 2 or 4 input pairs of fluxes (polarimetric exposures).
 * \param StokesIndex is the selected Stokes parameter
 * \param *iE[4] is a set of flux vectors for the input beams with a given state of polarization (ordinary beams)
 * \param *iA[4] is a set of flux vectors for the input beams with a given orthogonal state of polarization (extra-ordinary beams)
 * \param NumberOfExposures is the number of input exposures (accepts only 2 or 4)
 * \return void
 */
void operaPolarimetry::calculateStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures) {
    
    if(NumberOfExposures != 2 && NumberOfExposures != 4) {
        throw operaException("operaPolarimetry: NumberOfExposures not valid.", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
    }
    
    operaFluxVector Intensity(length);
    
    Intensity = 0.0;
    for(unsigned i=0;i<NumberOfExposures;i++) {
        Intensity += (*(iE[i]) + *(iA[i])) / ((double)NumberOfExposures*2.0); // DT Apr 26 2013 added *2.0
    }
    setStokesParameter(StokesI, &Intensity);
    
    if (StokesIndex == StokesQ || StokesIndex == StokesU || StokesIndex == StokesV) {
        if(!getHasDegreeOfStokes(StokesIndex)) {
            calculateDegreeOfPolarization(StokesIndex,iE,iA,NumberOfExposures);
        }
    }
    calculatePolarization(StokesIndex);
}


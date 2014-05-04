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

#ifndef OPERAPOLARIMETRY_H
#define OPERAPOLARIMETRY_H

#include "libraries/operaLibCommon.h"			// for CMatrix
#include "libraries/operaFluxVector.h"			// for FluxVector
#include "libraries/Polynomial.h"	
#include "libraries/operaMuellerMatrix.h"
#include "libraries/operaStokesVector.h"

/*!
 * \file operaPolarimetry.h
 * \brief This file holds the declaration of the class operaPolarimetry.
 * \ingroup libraries
 */

using namespace std;

/*!
 * \brief Definition of the method used to calculate polarization.
 */
typedef enum { Difference=1, Ratio, DifferenceWithBeamSwapped, NewMethod} method_t;

/*!
 * \author Andre Venne
 * \brief This class encapsulates the polarimetry results.
 * \sa class operaStokesVector
 * 
 * This class holds in operaStokesVector classes the 4 Stokes parameters, the 4 associated degrees of polarization
 * and the 2 null polarization spectra for each Stokes parameter.
 */
class operaPolarimetry {
	
private:
	unsigned length;
    method_t method;
    
	operaStokesVector *stokesVector;            // The 4 Stokes parameters I, Q, U, V
    operaStokesVector *degreeOfPolarization;	// The 4 degrees of polarization for each Stokes parameters I, Q, U, V
    operaStokesVector *firstNullPolarization;	// The 4 null polarization for the first null spectrum
    operaStokesVector *secondNullPolarization;	// The 4 null polarization for the second null spectrum
    
	double *wavelength;							// wavelength in nm 

    bool hasStokesI;
    bool hasStokesQ;
    bool hasStokesU;
    bool hasStokesV;
    
    bool hasDegreeOfStokesI;
    bool hasDegreeOfStokesQ;
    bool hasDegreeOfStokesU;
    bool hasDegreeOfStokesV;
    
    bool hasFirstNullPolarization;
    bool hasSecondNullPolarization;
    
	bool hasWavelength;    

public:
	/*
	 * Constructors / Destructors
	 */
    
    /*!
     * \brief Basic operaPolarimetry constructor.
     * \return void
     */
	operaPolarimetry();
    
    /*!
     * \brief Basic operaPolarimetry constructor.
     * \param Length An unsigned number of elements in each operaStokesVector
     * \return void
     */
	operaPolarimetry(unsigned Length);
    
    /*!
     * \brief Basic operaPolarimetry destructor.
     * \return void
     */
	~operaPolarimetry();

 	/*
	 * Getters/Setters
	 */
	
    /*!
     * \brief Sets the length of the Stokes vectors.
     * \details A function that sets the value of the variable holding the number of elements in the operaStokesVector.
     * \param Length An unsigned number of elements
     * \return void
     */
	void setLength(unsigned Length);
	
    /*!
     * \brief Gets the length of the Stokes vectors.
     * \details A function that gets the value of the variable holding the number of elements in the operaStokesVector.
     * \return An unsigned value
     */
	unsigned getLength(void);

    /*!
     * \brief Set the method used to calculate polarization.
     * \details Supported methods are Difference=1, Ratio=2, DifferenceWithBeamSwapped=3, NewMethod=4.
     * \param method_t
     * \return void
     */
	void setmethod(method_t Method) {method = Method;};
	
    /*!
     * \brief Get the method used to calculate polarization.
     * \details Supported methods are Difference=1, Ratio=2, DifferenceWithBeamSwapped=3, NewMethod=4.
     * \return A typedef method_t value
     */
	method_t getmethod(void) {return method;};
    
	/*!
	 * \brief Gets a wavelength vector address.
	 * \return Wavelength a wavelength vector address
	 */
	double *getwavelength(void);
	/*!
	 * \brief Gets a wavelength value at indexElem.
	 * \param indexEleme the index
	 * \return Wavelength a wavelength value
	 */
	double getwavelength(unsigned indexElem);
	/*!
	 * \brief Sets a wavelength value at indexElem.
	 * \param Wavelength a wavelength value
	 * \param indexEleme the index
	 * \return void
	 */
	void setwavelength(double Wavelength, unsigned indexElem);
	/*!
	 * \brief Copies wavelength values.
	 * \param double *Wavelength a wavelength vector
	 * \param length the number of elements
	 * \return void
	 */
	void setwavelength(double *Wavelength, unsigned length);
	/*!
	 * \brief Determines whether wavelength information is avilable.
	 * \return bool
	 */
	bool getHasWavelength() { return hasWavelength; };
	/*!
	 * \brief Sets whether wavelength information is avilable.
	 * \param bool HasWavelength
	 * \return void
	 */
	void setHasWavelength(bool HasWavelength) { hasWavelength = HasWavelength; };
    
    /*!
     * \brief Sets a Stokes parameter of the Stokes vector.
     * \details A function that sets the operaFluxVector of a Stokes parameter of the Stokes vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param FluxVector An operaFluxVector pointer
     * \return void
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector);
    
    /*!
     * \brief Sets a Stokes parameter element of the Stokes vector.
     * \details A function that sets the value and the variance of a Stokes parameter element of the Stokes vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param StokesValue A double value
     * \param Variance A double value
     * \param index An unsigned index to the element
     * \return void
     */
    void setStokesParameter(stokes_parameter_t StokesIndex, double StokesValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets the Stokes vector.
     * \details A function that gets the operaStokesVector of the Stokes vector.
     * \return An operaStokesVector pointer
     */
    operaStokesVector* getStokesVector(void);
    
    /*!
     * \brief Gets a Stokes parameter of the Stokes vector.
     * \details A function that gets the operaFluxVector of a Stokes parameter of the Stokes vector.
     * \param StokesIndex A stokes_parameter_t value
     * \return An operaFluxVector pointer
     */
    operaFluxVector* getStokesParameter(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Sets a Stokes parameter of the degree of polarization vector.
     * \details A function that sets the operaFluxVector of a Stokes parameter of the degree of polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param FluxVector An operaFluxVector pointer
     * \return void
     */
    void setDegreeOfPolarization(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector);
    
    /*!
     * \brief Sets a Stokes parameter element of the degree of polarization vector.
     * \details A function that sets the value and the variance of a Stokes parameter element of the degree of polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param DegreeOfPolarizationValue A double value
     * \param Variance A double value
     * \param index An unsigned index to the element
     * \return void
     */
    void setDegreeOfPolarization(stokes_parameter_t StokesIndex, double DegreeOfPolarizationValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets the degree of polarization vector.
     * \details A function that gets the operaStokesVector of the degree of polarization vector.
     * \return An operaStokesVector pointer
     */
    operaStokesVector* getDegreeOfPolarization(void);
    
    /*!
     * \brief Gets a Stokes parameter of the degree of polarization vector.
     * \details A function that gets the operaFluxVector of a Stokes parameter of the degree of polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \return An operaFluxVector pointer
     */
    operaFluxVector* getDegreeOfPolarization(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Sets a Stokes parameter of the first null polarization vector.
     * \details A function that sets the operaFluxVector of a Stokes parameter of the first null polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param FluxVector An operaFluxVector pointer
     * \return void
     */
    void setFirstNullPolarization(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector);
    
    /*!
     * \brief Sets a Stokes parameter element of the first null polarization vector.
     * \details A function that sets the value and the variance of a Stokes parameter element of the first null polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param FirstNullPolarizationValue A double value
     * \param Variance A double value
     * \param index An unsigned index to the element
     * \return void
     */
    void setFirstNullPolarization(stokes_parameter_t StokesIndex, double FirstNullPolarizationValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets the first null polarization vector.
     * \details A function that gets the operaStokesVector of the first null polarization vector.
     * \return An operaStokesVector pointer
     */
    operaStokesVector* getFirstNullPolarization(void);
    
    /*!
     * \brief Gets a Stokes parameter of the first null polarization vector.
     * \details A function that gets the operaFluxVector of a Stokes parameter of the first null polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \return An operaFluxVector pointer
     */
    operaFluxVector* getFirstNullPolarization(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Sets a Stokes parameter of the second null polarization vector.
     * \details A function that sets the operaFluxVector of a Stokes parameter of the second null polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param FluxVector An operaFluxVector pointer
     * \return void
     */
    void setSecondNullPolarization(stokes_parameter_t StokesIndex, operaFluxVector *FluxVector);
    
    /*!
     * \brief Sets a Stokes parameter element of the second null polarization vector.
     * \details A function that sets the value and the variance of a Stokes parameter element of the second null polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \param SecondNullPolarizationValue A double value
     * \param Variance A double value
     * \param index An unsigned index to the element
     * \return void
     */
    void setSecondNullPolarization(stokes_parameter_t StokesIndex, double SecondNullPolarizationValue, double Variance, unsigned index);
    
    /*!
     * \brief Gets the second null polarization vector.
     * \details A function that gets the operaStokesVector of the second null polarization vector.
     * \return An operaStokesVector pointer
     */
    operaStokesVector* getSecondNullPolarization(void);
    
    /*!
     * \brief Gets a Stokes parameter of the second null polarization vector.
     * \details A function that gets the operaFluxVector of a Stokes parameter of the second null polarization vector.
     * \param StokesIndex A stokes_parameter_t value
     * \return An operaFluxVector pointer
     */
    operaFluxVector* getSecondNullPolarization(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Gets the boolean value of Stokes I.
     * \details A function that gets the boolean value giving the state of the Stokes I parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasStokesI(void);
    
    /*!
     * \brief Gets the boolean value of Stokes Q.
     * \details A function that gets the boolean value giving the state of the Stokes Q parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasStokesQ(void);
    
    /*!
     * \brief Gets the boolean value of Stokes U.
     * \details A function that gets the boolean value giving the state of the Stokes U parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasStokesU(void);
    
    /*!
     * \brief Gets the boolean value of Stokes V.
     * \details A function that gets the boolean value giving the state of the Stokes V parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasStokesV(void);
    
    /*!
     * \brief Get the boolean value of a given Stokes.
     * \details A function that gets the boolean value giving the state of a given Stokes parameter in the Stokes vector.
     * \param StokesIndex choice of Stokes parameter
     * \return A boolean value
     */
    bool getHasStokes(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Sets the boolean value of Stokes I.
     * \details A function that sets the boolean value giving the state of the Stokes I parameter in the Stokes vector.
     * \param HasStokesI A boolean value
     * \return void
     */
    void setHasStokesI(bool HasStokesI);
    
    /*!
     * \brief Sets the boolean value of Stokes Q.
     * \details A function that sets the boolean value giving the state of the Stokes Q parameter in the Stokes vector.
     * \param HasStokesQ A boolean value
     * \return void
     */
    void setHasStokesQ(bool HasStokesQ);
    
    /*!
     * \brief Sets the boolean value of Stokes U.
     * \details A function that sets the boolean value giving the state of the Stokes U parameter in the Stokes vector.
     * \param HasStokesU A boolean value
     * \return void
     */
    void setHasStokesU(bool HasStokesU);
    
    /*!
     * \brief Sets the boolean value of Stokes V.
     * \details A function that sets the boolean value giving the state of the Stokes V parameter in the Stokes vector.
     * \param HasStokesV A boolean value
     * \return void
     */
    void setHasStokesV(bool HasStokesV);

    /*!
     * \brief Set the boolean value of a given Stokes.
     * \details A function that sets the boolean value giving the state of a given Stokes parameter in the Stokes vector.
     * \param HasStokes A boolean value
     * \param StokesIndex choice of Stokes parameter
     * \return void
     */
    void setHasStokes(stokes_parameter_t StokesIndex, bool HasStokes);

    /*!
     * \brief Gets the boolean value of the degree of polarization of Stokes I.
     * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes I parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasDegreeOfStokesI(void);
    
    /*!
     * \brief Gets the boolean value of the degree of polarization of Stokes Q.
     * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes Q parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasDegreeOfStokesQ(void);
    
    /*!
     * \brief Gets the boolean value of the degree of polarization of Stokes U.
     * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes U parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasDegreeOfStokesU(void);
    
    /*!
     * \brief Gets the boolean value of the degree of polarization of Stokes V.
     * \details A function that gets the boolean value giving the state of the degree of polarization of the Stokes V parameter in the Stokes vector.
     * \return A boolean value
     */
    bool getHasDegreeOfStokesV(void);
    
    /*!
     * \brief Get the boolean value of the degree of polarization of a given Stokes.
     * \details A function that gets the boolean value giving the state of the degree of polarization of a given Stokes parameter in the Stokes vector.
     * \param StokesIndex choice of Stokes parameter
     * \return A boolean value
     */
    bool getHasDegreeOfStokes(stokes_parameter_t StokesIndex);

    /*!
     * \brief Sets the boolean value of the degree of polarization of Stokes I.
     * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes I parameter in the Stokes vector.
     * \param HasDegreeOfStokesI A boolean value
     * \return void
     */
    void setHasDegreeOfStokesI(bool HasDegreeOfStokesI);
    
    /*!
     * \brief Sets the boolean value of the degree of polarization of Stokes Q.
     * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes Q parameter in the Stokes vector.
     * \param HasDegreeOfStokesQ A boolean value
     * \return void
     */
    void setHasDegreeOfStokesQ(bool HasDegreeOfStokesQ);
    
    /*!
     * \brief Sets the boolean value of the degree of polarization of Stokes U.
     * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes U parameter in the Stokes vector.
     * \param HasDegreeOfStokesU A boolean value
     * \return void
     */
    void setHasDegreeOfStokesU(bool HasDegreeOfStokesU);
    
    /*!
     * \brief Sets the boolean value of the degree of polarization of Stokes V.
     * \details A function that sets the boolean value giving the state of the degree of polarization of the Stokes V parameter in the Stokes vector.
     * \param HasDegreeOfStokesV A boolean value
     * \return void
     */
    void setHasDegreeOfStokesV(bool HasDegreeOfStokesV);
    
    /*!
     * \brief Set the boolean value of the degree of polarization of a given Stokes.
     * \details A function that sets the boolean value giving the state of the degree of polarization of a given Stokes parameter in the Stokes vector.
     * \param HasDegreeOfStokes A boolean value
     * \param StokesIndex choice of Stokes parameter     
     * \return void
     */
    void setHasDegreeOfStokes(stokes_parameter_t StokesIndex, bool HasDegreeOfStokes);
    
    /*!
     * \brief Gets the boolean value of the first null polarization vector.
     * \details A function that gets the boolean value of the first null polarization vector.
     * \return A boolean value
     */
    bool getHasFirstNullPolarization(void);
    
    /*!
     * \brief Gets the boolean value of the second null polarization vector.
     * \details A function that gets the boolean value of the second null polarization vector.
     * \return A boolean value
     */
    bool getHasSecondNullPolarization(void);
    
    /*!
     * \brief Sets the boolean value of the first null polarization vector.
     * \details A function that sets the boolean value of the first null polarization vector.
     * \param HasFirstNullPolarization A boolean value
     * \return void
     */
    void setHasFirstNullPolarization(bool HasFirstNullPolarization);
    
    /*!
     * \brief Sets the boolean value of the second null polarization vector.
     * \details A function that sets the boolean value of the second null polarization vector.
     * \param HasSecondNullPolarization A boolean value
     * \return void
     */
    void setHasSecondNullPolarization(bool HasSecondNullPolarization);
    
	/*
	 * Methods
	 */
    
    /*!
     * \brief Calculates the polarization.
     * \details A function that calculates the polarization from the degree of polarization and stores it in the Stokes vector.
     * \return void
     */
    void calculatePolarization(void);
    
    void calculatePolarization(stokes_parameter_t StokesIndex);
    
    /*!
     * \brief Calculates the degree of polarization.
     * \details A function that calculates the degree of polarization from the polarization and stores it in the degree of polarization vector.
     * \return void
     */
    void calculateDegreeOfPolarization(void);
    
    /*!
     * \brief Calculates the degree of polarization.
     * \details A function that calculates the degree of polarization from the observed ordinary and extraordinary beam fluxes.
     * \details This function accepts 2 or 4 input pairs of fluxes (polarimetric exposures).
     * \param StokesIndex is the selected Stokes parameter
     * \param *iE[4] is a set of flux vectors for the input beams with a given state of polarization (ordinary beams)
     * \param *iA[4] is a set of flux vectors for the input beams with a given orthogonal state of polarization (extra-ordinary beams)
     * \param NumberOfExposures is the number of input exposures (accepts only 2 or 4)
     * \return void
     */
    void calculateDegreeOfPolarization(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures);

    /*!
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
    void calculateStokesParameter(stokes_parameter_t StokesIndex, operaFluxVector *iE[4], operaFluxVector *iA[4], unsigned NumberOfExposures);
    
    
};
#endif

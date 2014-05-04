/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: oepraFluxvector
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

#include <cstring>		// for malloc
#include <stdlib.h>		// for memset

#include "globaldefines.h"
#include "operaError.h"
#include <libraries/operaException.h>
#include <libraries/operaFluxVector.h>

/*!
 * \brief This file holds the implementation of the class operaFluxVector.
 * \file operaFluxVector.cpp
 * \ingroup libraries
 */

/*
 * \author Doug Teeple
 * \author Andre Venne
 * \brief This class encapsulates the flux vector.
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

/*
 * Constructors / Destructors
 */

/*
 * \brief operaFluxVector constructor.
 * \details This constructor creates an operaFluxVector with vectors of content 0.0.
 * \param Length An unsigned number of elements in the operaFluxVector
 * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaFluxVector::operaFluxVector(unsigned Length, TendsTowards_t Towards, bool Istemp) :
length(0),
towards(ToDefault),
fluxes(NULL),
variances(NULL)
{		
	if (Length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	fluxes = new double[Length]; 
	variances = new double[Length];
	length = Length;
	istemp = Istemp;
	towards = Towards;
}

/*
 * \brief operaFluxVector constructor from a flux pointer and his variance pointer.
 * \details This constructor creates an operaFluxVector with vectors of content given by a flux double pointer and his variance double pointer.
 * \param Fluxes A double pointer
 * \param Variances A double pointer
 * \param Length An unsigned number of elements in the operaFluxVector
 * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaFluxVector::operaFluxVector(double *Fluxes, double *Variances, unsigned Length, TendsTowards_t Towards, bool Istemp) :
length(0),
towards(ToDefault),
fluxes(NULL),
variances(NULL)
{		
	if (Length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	fluxes = new double[Length]; 
	variances = new double[Length];
	length = Length;
	istemp = Istemp;
	towards = Towards;
	memcpy(fluxes, Fluxes, sizeof(double)*Length);
	memcpy(variances, Variances, sizeof(double)*Length);
}

/*
 * \brief operaFluxVector constructor from an operaFluxVector.
 * \details This constructor creates an operaFluxVector with vectors of content given by an operaFluxVector.
 * \param b An operaFluxVector address
 * \param Length An unsigned number of elements in the operaFluxVector
 * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaFluxVector::operaFluxVector(operaFluxVector &b, TendsTowards_t Towards, bool Istemp) :
length(0),
towards(ToDefault),
fluxes(NULL),
variances(NULL)
{		
	if (b.length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	fluxes = new double[b.length]; 
	variances = new double[b.length];
	length = b.length;
	istemp = Istemp;
	towards = Towards;
	memcpy(fluxes, b.fluxes, sizeof(double)*length);
	memcpy(variances, b.variances, sizeof(double)*length);
}

/*
 * \brief operaFluxVector constructor from an operaFluxVector.
 * \details This constructor creates an operaFluxVector with vectors of content given by an operaFluxVector.
 * \param b An operaFluxVector pointer
 * \param Length An unsigned number of elements in the operaFluxVector
 * \param TendsTowards_t An optional TendsTowards_t value defaults to ToDefault
 * \param Istemp An optional bool defaults to false
 * \return void
 */
operaFluxVector::operaFluxVector(operaFluxVector *b, TendsTowards_t Towards, bool Istemp) :
length(0),
towards(ToDefault),
fluxes(NULL),
variances(NULL)
{		
	if (b->length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	fluxes = new double[b->length]; 
	variances = new double[b->length];
	length = b->length;
	istemp = Istemp;
	towards = Towards;
	memcpy(fluxes, b->fluxes, sizeof(double)*length);
	memcpy(variances, b->variances, sizeof(double)*length);
}

/*
 * \brief Destructor releases the pointers memory.
 * \return void
 */
operaFluxVector::~operaFluxVector(void) {
	if (fluxes) {
		delete[] fluxes;
	}
	fluxes = NULL;
	if (variances) {
		delete[] variances;
	}
	variances = NULL;
	length = 0;
}

/*
 * Getters/Setters
 */

/*
 * \brief Sets the variance vector.
 * \details A function that sets the variance vector content to the value of a given double pointer.
 * \param Variances A double pointer
 * \return void
 */
void operaFluxVector::setVarianceVector(double *Variances)
{
	if (length == 0) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	memcpy(variances, Variances, sizeof(double)*length);
}

/*
 * \brief Sets the flux vector.
 * \details A function that sets the flux vector content to the value of a given double pointer.
 * \param Fluxes A double pointer
 * \return void
 */
void operaFluxVector::setVector(double *Fluxes)
{
	if (length == 0) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	memcpy(fluxes, Fluxes, sizeof(double)*length);
}

/* 
 * \brief Sets the flux vector.
 * \param FluxVector An operaFluxVector address
 * \param Fluxes A double pointer
 * \param Variances A double pointer
 * \return void
 */
void operaFluxVector::setVector(operaFluxVector &Fluxvector) {
	if (length == 0 || Fluxvector.length == 0) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (length < Fluxvector.length) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	length = Fluxvector.length;
	memcpy(fluxes, Fluxvector.fluxes, sizeof(double)*length);
	memcpy(variances, Fluxvector.variances, sizeof(double)*length);
}

/*
 * \brief Sets the flux and variance vector.
 * \details A function that sets the flux and variance vector content to the value of given double pointers.
 * \param Fluxes A double pointer
 * \param Variances A double pointer
 * \return void
 */
void operaFluxVector::setVectors(double *Fluxes, double *Variances)
{
	if (length == 0) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	memcpy(fluxes, Fluxes, sizeof(double)*length);
	memcpy(variances, Variances, sizeof(double)*length);
}



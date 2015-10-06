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

operaFluxVector::operaFluxVector(unsigned Length, TendsTowards_t Towards) :
length(Length), towards(Towards)
{
	if (length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	fluxes = new double[length]; 
	variances = new double[length];
}

operaFluxVector::operaFluxVector(double *Fluxes, double *Variances, unsigned Length, TendsTowards_t Towards) :
length(Length), towards(Towards)
{
	if (length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	fluxes = new double[length]; 
	variances = new double[length];
	setVectors(Fluxes, Variances);
}

operaFluxVector::operaFluxVector(const operaFluxVector &b, TendsTowards_t Towards) :
length(b.length), towards(Towards)
{
	if (length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	fluxes = new double[length]; 
	variances = new double[length];
	setVectors(b.fluxes, b.variances);
}

operaFluxVector::~operaFluxVector(void)
{
	if (fluxes) delete[] fluxes;
	if (variances) delete[] variances;
	fluxes = NULL;
	variances = NULL;
	length = 0;
}

/*
 * Getters/Setters
 */

void operaFluxVector::setFluxVector(double *Fluxes)
{
	copy(Fluxes, Fluxes+length, fluxes);
}

void operaFluxVector::setVarianceVector(double *Variances)
{
	copy(Variances, Variances+length, variances);
}

void operaFluxVector::setVectors(double *Fluxes, double *Variances)
{
	setFluxVector(Fluxes);
	setVarianceVector(Variances);
}

void operaFluxVector::setVector(const operaFluxVector &Fluxvector) {
	if (length != Fluxvector.length) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);	
	}
	setVectors(Fluxvector.fluxes, Fluxvector.variances);
}

void operaFluxVector::resize(unsigned newlength) {
	if (newlength == length) return;
	double* oldflux = fluxes;
	double* oldvar = variances;
	fluxes = new double[newlength];
	variances = new double[newlength];
	if (newlength > length) {
		setVectors(oldflux, oldvar);
		length = newlength;
	} else {
		length = newlength;
		setVectors(oldflux, oldvar);
	}
	if(oldflux) delete[] oldflux;
	if(oldvar) delete[] oldvar;
}

/*
 * Operators
 */

operaFluxVector& operaFluxVector::operator=(const operaFluxVector& b) {
	if (b.length == 0) {
		throw operaException("operaFluxVector: ", operaErrorZeroLength, __FILE__, __FUNCTION__, __LINE__);	
	}
	if (length != b.length) {
		if (fluxes) delete[] fluxes;
		if (variances) delete[] variances;
		length = b.length;
		fluxes = new double[length]; 
		variances = new double[length];
	}
	towards = b.towards;
	setVector(b);
	return *this;
}

operaFluxVector& operaFluxVector::operator+=(const operaFluxVector& b) {
	if (length != b.length) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	for(unsigned i = 0; i < length; i++) {
		fluxes[i] += b.fluxes[i];
		variances[i] += b.variances[i];
	}
	return *this;
}

operaFluxVector& operaFluxVector::operator-=(const operaFluxVector& b) {
	if (length != b.length) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	for(unsigned i = 0; i < length; i++) {
		fluxes[i] -= b.fluxes[i];
		variances[i] += b.variances[i];
	}
	return *this;
}

operaFluxVector& operaFluxVector::operator*=(const operaFluxVector& b) {
	if (length != b.length) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	for(unsigned i = 0; i < length; i++) {
		variances[i] = fluxes[i]*fluxes[i]*b.variances[i] + b.fluxes[i]*b.fluxes[i]*variances[i];
		fluxes[i] *= b.fluxes[i];
	}
	return *this;
}

operaFluxVector& operaFluxVector::operator/=(const operaFluxVector& b) {
	if (length != b.length) {
		throw operaException("operaFluxVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	}
	for(unsigned i = 0; i < length; i++) {
		if (towards == ToDefault || !isinf(fluxes[i]) || !isinf(b.fluxes[i])) {
			const double b2 = b.fluxes[i]*b.fluxes[i]; //small optimization
			const double aoverb2 = fluxes[i]/b2; //small optimization
			variances[i] = variances[i]/b2 + aoverb2 * aoverb2 * b.variances[i]; // = (v(a)*b^2 + v(b)*a^2)/b^4
			fluxes[i] /= b.fluxes[i];
		} else if (towards == ToINF) {
			fluxes[i] = FP_INFINITE;
			variances[i] = 0.0;
		} else if (towards == ToNAN) {
			fluxes[i] = FP_NAN;
			variances[i] = FP_NAN;
		} else if (towards == ToZero) {
			fluxes[i] = 0.0;
			variances[i] = 0.0;
		} else if (towards == ToOne) {
			fluxes[i] = 1.0;
			variances[i] = 0.0;
		}
	}
	return *this;
}

operaFluxVector& operaFluxVector::operator=(double d) {
	fill_n(fluxes, length, d);
	fill_n(variances, length, d);
	return *this;
}

operaFluxVector& operaFluxVector::operator+=(double d) {
	for(unsigned i = 0; i < length; i++) fluxes[i] += d;
	return *this;
}

operaFluxVector& operaFluxVector::operator-=(double d) {
	for(unsigned i = 0; i < length; i++) fluxes[i] -= d;
	return *this;
}

operaFluxVector& operaFluxVector::operator*=(double d) {
	for(unsigned i = 0; i < length; i++) {
		fluxes[i] *= d;
		variances[i] *= d*d;
	}
	return *this;
}

operaFluxVector& operaFluxVector::operator/=(double d) {
	for(unsigned i = 0; i < length; i++) {
		fluxes[i] /= d;
		variances[i] /= d*d;
	}
	return *this;
}

operaFluxVector Sqrt(operaFluxVector b) {
	for(unsigned i = 0; i < b.getlength(); i++) {
		b.setvariance(b.getvariance(i) / (b.getflux(i)*4.0), i);
		b.setflux(sqrt(b.getflux(i)), i);
	}
	return b;
}

operaFluxVector Pow(operaFluxVector b, double d) {
	for(unsigned i = 0; i < b.getlength(); i++) {
		b.setvariance(pow(pow(b.getflux(i), d-1.0)*d, 2) * b.getvariance(i), i);
		b.setflux(pow(b.getflux(i), d), i);
	}
	return b;
}

pair<double,double> Sum(const operaFluxVector& b) {
	double sum = 0.0, var = 0.0;
	for(unsigned i = 0; i < b.getlength(); i++) {
		sum += b.getflux(i);
		var += b.getvariance(i);
	}
	return pair<double,double>(sum, var);
}

pair<double,double> Mean(const operaFluxVector& b) {
	double sum = 0.0, var = 0.0;
	for(unsigned i = 0; i < b.getlength(); i++) {
		sum += b.getflux(i);
		var += b.getvariance(i);
	}
	sum /= b.getlength();
	var /= b.getlength();
	return pair<double,double>(sum, var);
}

void resizeVector(double *&vector, unsigned oldlength, unsigned newlength) {
	if (newlength == oldlength) return;
	double *newvector = new double[newlength];
	if(newlength < oldlength) copy(vector, vector+newlength, newvector);
	else copy(vector, vector+oldlength, newvector);
	if(vector) delete[] vector;
	vector = newvector;
}

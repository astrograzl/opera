#ifndef OPERAVECTOROPERATIONS_H
#define OPERAVECTOROPERATIONS_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: operaVectorOperations
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

#include "libraries/operaVector.h"
#include <functional>

/*!
 * \brief Performs a unary operation on all elements of an operaVector, storing the results in that vector.
 * \param funct The operation to apply.
 * \param a The vector to be apply the operation to.
 * \return A reference to the vector post-operation.
 */
template <typename T> operaVector& ApplyOperation(T funct, operaVector& a);

/*!
 * \brief Performs a binary operation on all elements of two operaVectors, storing the results in the first vector.
 * \param funct The operation to apply.
 * \param a The vector to be apply the operation to.
 * \param b The second operand for the operation.
 * \return A reference to the vector post-operation.
 */
template <typename T> operaVector& ApplyOperation(T funct, operaVector& a, const operaVector& b);

/*!
 * \brief Performs a binary operation with a constant value on all elements of an operaVector, storing the results in that vector.
 * \param funct The operation to apply.
 * \param a The vector to be apply the operation to.
 * \param d The second operand for the operation.
 * \return A reference to the vector post-operation.
 */
template <typename T> operaVector& ApplyOperation(T funct, operaVector& a, double d);

/*!
 * \brief Performs an operation that returns a single value from all elements of an operaVector.
 * \param funct The binary operation to apply between each element and the current result value.
 * \param a The vector of elements to use.
 * \param initalValue The value to start with for the result.
 * \return The value that results from applying funct with each value of the vector.
 * \details For example, the product of elements of a vector could be calculated using muliplication as funct, and an initial value of 1.
 */
template <typename T> double ScalarOperation(T funct, const operaVector& a, const double initalValue);

/*!
 * \brief Performs an operation that returns a single value from all elements two operaVectors. Equivalent to applying a scalar operation to the result of a binary operation.
 * \param funct The operation to use for the scalar operation.
 * \param pairwiseOp The operation to use for the binary operation.
 * \param a The first vector of elements to use.
 * \param a The second vector of elements to use.
 * \param initalValue The initial value for the scalar operation.
 * \return The value that results from applying funct with each value of the vector that results from applying pairwiseOp on the two input vectors.
 * \details For example, the dot product of two vectors could be calculated using additon as funct, multiplication as pairwiseOp, and an initial value of 0.
 */
template <typename T, typename U> double ScalarOperation(T funct, U pairwiseOp, const operaVector& a, const operaVector& b, const double initalValue);

/*!
 * \brief Removes the type ambiguity from a function pointer for a unary double operation by converting it into a functor.
 * \param funct The function pointer.
 * \return A function wrapper object.
 */
std::pointer_to_unary_function<double, double> ToFunctor(double (*funct)(double));

/*!
 * \brief Removes the type ambiguity from a function pointer for a binary double operation by converting it into a functor.
 * \param funct The function pointer.
 * \return A function wrapper object.
 */	
std::pointer_to_binary_function<double, double, double> ToFunctor(double (*funct)(double, double));

/*!
 * \brief Performs an element-wise addition of one operaVector to another and assigns the result to the vector.
 * \param a The vector to be added to.
 * \param b The vector to be added.
 * \return A reference to the vector post-operation.
 */
operaVector& operator+=(operaVector& a, const operaVector& b);

/*!
 * \brief Adds a value to each element of a vector and assigns the result to the vector.
 * \param a The vector to be added to.
 * \param d The value to be added.
 * \return A reference to the vector post-operation.
 */
operaVector& operator+=(operaVector& a, double d);

/*!
 * \brief Performs an element-wise subtraction of one operaVector from another and assigns the result to the vector.
 * \param a The vector to be subtracted from.
 * \param b The vector to be subtracted.
 * \return A reference to the vector post-operation.
 */
operaVector& operator-=(operaVector& a, const operaVector& b);

/*!
 * \brief Subtracts a value from each element of the vector and assigns the result to the vector.
 * \param a The vector to be subtracted from.
 * \param d The value to be subtracted.
 * \return A reference to the vector post-operation.
 */
operaVector& operator-=(operaVector& a, double d);

/*!
 * \brief Performs an element-wise multiplication of one operaVector by another and assigns the result to the vector.
 * \param a The vector to be multiplied.
 * \param b The vector to be multiplied by.
 * \return A reference to the vector post-operation.
 */
operaVector& operator*=(operaVector& a, const operaVector& b);

/*!
 * \brief Mulitiplies each element of the vector by a value and assigns the result to the vector.
 * \param a The vector to be multiplied.
 * \param d The value to be multiplied by.
 * \return A reference to the vector post-operation.
 */
operaVector& operator*=(operaVector& a, double d);

/*!
 * \brief Performs an element-wise division of one operaVector by another and assigns the result to the vector.
 * \param a The vector to be divided.
 * \param b The vector to be divided by.
 * \return A reference to the vector post-operation.
 */
operaVector& operator/=(operaVector& a, const operaVector& b);

/*!
 * \brief Divides each element of the vector by a value and assigns the result to the vector.
 * \param a The vector to be divided.
 * \param d The value to be divided by.
 * \return A reference to the vector post-operation.
 */
operaVector& operator/=(operaVector& a, double d);

/*!
 * \brief Calculates the minimum element of the vector.
 * \param a The vector.
 * \return The minimum.
 */
double Min(const operaVector& a);

/*!
 * \brief Calculates the minimum element of the vector.
 * \param a The vector.
 * \return The maximum.
 */
double Max(const operaVector& a);

/*!
 * \brief Calculates the median of the elements in the vector.
 * \param a The vector.
 * \return The median.
 */
double Median(operaVector a);

/*!
 * \brief Gets the absolute value of each element of an operaVector.
 * \param a The vector to take the absolute value of.
 * \return The resultant vector.
 */
operaVector Abs(operaVector a);

/*!
 * \brief Gets the square root of each element of an operaVector.
 * \param a The vector to take the square root of.
 * \return The resultant vector.
 */
operaVector Sqrt(operaVector a);

/*!
 * \brief Raises each element of an operaVector to a power.
 * \param a The vector to raise to a power.
 * \param d The power to raise each element to.
 * \return The resultant vector.
 */
operaVector Pow(operaVector a, double d);

/*!
 * \brief Calculates the inner-product (dot product) of one vector with another vector.
 * \param b An operaVector.
 * \param b An operaVector of the same size.
 * \return The calculated inner-product.
 */
double InnerProduct(const operaVector& a, const operaVector& b);

/*!
 * \brief Calculates the magnitude/norm of the vector.
 * \param a The vector.
 * \return The magnitude.
 */
double Magnitude(const operaVector& a);

/*!
 * \brief Calculates the sum of the elements in the vector.
 * \param a The vector.
 * \return The sum.
 */
double Sum(const operaVector& a);

/*!
 * \brief Calculates the mean of the elements in the vector.
 * \param a The vector.
 * \return The mean.
 */
double Mean(const operaVector& a);

/*!
 * \brief Calculates the variance of the elements in the vector.
 * \param a The vector.
 * \return The variance.
 */
double Variance(const operaVector& a);

/*!
 * \brief Calculates the variance of the elements in the vector from a given mean.
 * \param a The vector.
 * \param mean The already calculated mean.
 * \return The variance.
 */
double Variance(operaVector a, double mean);

/*!
 * \brief Calculates the standard deviation of the elements in the vector.
 * \param a The vector.
 * \return The standard deviation.
 */
double StdDev(const operaVector& a);

/*!
 * \brief Calculates the standard deviation of the elements in the vector from a given mean.
 * \param a The vector.
 * \param mean The already calculated mean.
 * \return The standard deviation.
 */
double StdDev(const operaVector& a, double mean);

/*!
 * \brief Calculates the median absolute deviation of the elements in the vector.
 * \param a The vector.
 * \return The median absolute deviation.
 */
double MedianAbsDev(const operaVector& a);

/*!
 * \brief Calculates the median absolute deviation of the elements in the vector from a given median.
 * \param a The vector.
 * \param median The already calculated median.
 * \return The median absolute deviation.
 */
double MedianAbsDev(operaVector a, double median);

/*!
 * \brief Calculates the median standard deviation of the elements in the vector.
 * \param a The vector.
 * \details Note: divides by 0.674433 to convert from absolute deviation to standard deviation.
 * \return The median standard deviation.
 */
double MedianStdDev(const operaVector& a);

/*!
 * \brief Calculates the median standard deviation of the elements in the vector.
 * \param a The vector.
 * \param median The already calculated median.
 * \details Note: divides by 0.674433 to convert from absolute deviation to standard deviation.
 * \return The median standard deviation.
 */
double MedianStdDev(const operaVector& a, double median);

template <typename T> operaVector& ApplyOperation(T funct, operaVector& a, const operaVector& b);

/*!
 * \brief Performs a unary operation on all elements of an operaVector, returning the result.
 * \param funct The operation to apply.
 * \param a The operand for the operation.
 * \return A vector containing the result of the operation.
 */
template <typename T> operaVector Operation(T funct, operaVector a) { return ApplyOperation(funct, a); }

/*!
 * \brief Performs a binary operation on all elements of two operaVectors, returning the result.
 * \param funct The operation to apply.
 * \param a The first operand for the operation.
 * \param b The second operand for the operation.
 * \return A vector containing the result of the operation.
 */
template <typename T> operaVector Operation(T funct, operaVector a, const operaVector& b) { return ApplyOperation(funct, a, b); }


/*!
 * \brief Performs a binary operation with a constant value on all elements of an operaVector, returning the result.
 * \param funct The operation to apply.
 * \param a The first operand for the operation.
 * \param d The second operand for the operation.
 * \return A vector containing the result of the operation.
 */
template <typename T> operaVector Operation(T funct, operaVector a, double d) { return ApplyOperation(funct, a, d); }

/*!
 * \brief Performs an element-wise addition of two operaVectors without modifying either operaVector.
 * \param a The first vector to be added.
 * \param b The second vector to be added.
 * \return The result of the addition.
 */
operaVector operator+(operaVector a, const operaVector& b) { return a += b; }

/*!
 * \brief Performs an element-wise subtraction of two operaVectors without modifying either operaVector.
 * \param a The vector of elements to subtract from.
 * \param b The vector of elements to subtract.
 * \return The result of the subtraction.
 */
operaVector operator-(operaVector a, const operaVector& b) { return a -= b; }

/*!
 * \brief Performs an element-wise multiplication of two operaVectors without modifying either operaVector.
 * \param a The first vector to be multiplied.
 * \param b The second vector to be multiplied.
 * \return The result of the multiplication.
 */
operaVector operator*(operaVector a, const operaVector& b) { return a *= b; }

/*!
 * \brief Performs an element-wise division of two operaVectors without modifying either operaVector.
 * \param a The vector of elements to be divided.
 * \param b The vector of elements to divide by.
 * \return The result of the division.
 */
operaVector operator/(operaVector a, const operaVector& b) { return a /= b; }

/*!
 * \brief Adds a value to each element of an operaVectorwithout modifying the original operaVector.
 * \param a The vector to be added.
 * \param d The value to be added.
 * \return The result of the addition.
 */
operaVector operator+(operaVector a, double d) { return a += d; }

/*!
 * \brief Subtracts a value from each element of an operaVector without modifying the original operaVector.
 * \param a The vector to be subtracted from.
 * \param d The value to be subtracted.
 * \return The result of the subtraction.
 */
operaVector operator-(operaVector a, double d) { return a -= d; }

/*!
 * \brief Multiplies each element of an operaVector by a value without modifying the original operaVector.
 * \param a The vector to be multiplied.
 * \param d The value to be multiplied by.
 * \return The result of the multiplication.
 */
operaVector operator*(operaVector a, double d) { return a *= d; }

/*!
 * \brief Divides each element of an operaVector by a value without modifying the original operaVector.
 * \param a The vector to be divided.
 * \param d The value to be divided by.
 * \return The result of the division.
 */
operaVector operator/(operaVector a, double d) { return a /= d; }

#include "libraries/operaException.h"
#include <algorithm>
#include <numeric>

template <typename T>
operaVector& ApplyOperation(T funct, operaVector& a) {
	std::transform(a.begin(), a.end(), a.begin(), funct);
	return a;
}

template <typename T>
operaVector& ApplyOperation(T funct, operaVector& a, const operaVector& b) {
	if(a.size() != b.size()) throw operaException("operaVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	std::transform(a.begin(), a.end(), b.begin(), a.begin(), funct);
	return a;
}

template <typename T>
operaVector& ApplyOperation(T funct, operaVector& a, double d) {
	std::transform(a.begin(), a.end(), a.begin(), std::bind2nd(funct, d));
	return a;
}

template <typename T>
double ScalarOperation(T funct, const operaVector& a, double initalValue) {
	return std::accumulate(a.begin(), a.end(), initalValue, funct);
}

template <typename T, typename U>
double ScalarOperation(T funct, U pairwiseOp, const operaVector& a, const operaVector& b, double initalValue) {
	if(a.size() != b.size()) throw operaException("operaVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	return std::inner_product(a.begin(), a.end(), b.begin(), initalValue, funct, pairwiseOp);
}

#endif

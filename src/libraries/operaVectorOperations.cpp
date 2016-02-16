#include "globaldefines.h"
#include "libraries/operaVectorOperations.h"
#include "libraries/operaException.h"
#include <functional>
#include <cmath>
#include <numeric>

operaVector& operator+=(operaVector& a, const operaVector& b) {
	return ApplyOperation(std::plus<double>(), a, b);
}

operaVector& operator+=(operaVector& a, double d) {
	return ApplyOperation(std::plus<double>(), a, d);
}

operaVector& operator-=(operaVector& a, const operaVector& b) {
	return ApplyOperation(std::minus<double>(), a, b);
}

operaVector& operator-=(operaVector& a, double d) {
	return ApplyOperation(std::minus<double>(), a, d);
}

operaVector& operator*=(operaVector& a, const operaVector& b) {
	return ApplyOperation(std::multiplies<double>(), a, b);
}

operaVector& operator*=(operaVector& a, double d) {
	return ApplyOperation(std::multiplies<double>(), a, d);
}

operaVector& operator/=(operaVector& a, const operaVector& b) {
	return ApplyOperation(std::divides<double>(), a, b);
}

operaVector& operator/=(operaVector& a, double d) {
	return ApplyOperation(std::divides<double>(), a, d);
}

operaVector Abs(operaVector a) {
	return ApplyOperation(ToFunctor(std::abs), a);
}

operaVector Sqrt(operaVector a) {
	return ApplyOperation(ToFunctor(std::sqrt), a);
}

operaVector Pow(operaVector a, double d) {
	return ApplyOperation(ToFunctor(std::pow), a, d);
}

double InnerProduct(const operaVector& a, const operaVector& b) {
	return ScalarOperation(std::plus<double>(), std::multiplies<double>(), a, b, 0);
}

double Magnitude(const operaVector& a) {
	return std::sqrt(InnerProduct(a, a));
}

double Sum(const operaVector& a) {
	return ScalarOperation(std::plus<double>(), a, 0);
}

double Mean(const operaVector& a) {
	return Sum(a) / a.size();
}

double Variance(const operaVector& a) {
	return Variance(a, Mean(a));
}

double Variance(operaVector a, double mean) {
	a -= mean;
	return InnerProduct(a, a) / a.size();
}

double StdDev(const operaVector& a) {
	return StdDev(a, Mean(a));
}

double StdDev(const operaVector& a, double mean) {
	return std::sqrt(Variance(a, mean));
}

double Min(const operaVector& a) {
	return *(std::min_element(a.begin(), a.end()));
}

double Max(const operaVector& a) {
	return *(std::max_element(a.begin(), a.end()));
}

double Median(operaVector a) {
	unsigned middle = a.size() / 2; //Middle index if size is odd, index after the middle if size is even
	std::nth_element(a.begin(), a.begin() + middle, a.end()); //Partial sort putting element at middle in correct position, all elements less before it, all elements greater after it
	if(a.size() % 2 == 1) return a[middle];
	double prev = *(std::max_element(a.begin(), a.begin() + middle)); //Find the largest element that is less than middle
	return (prev + a[middle]) / 2.0;
}

double MedianAbsDev(const operaVector& a) {
	return MedianAbsDev(a, Median(a));
}

double MedianAbsDev(operaVector a, double median) {
	a -= median;
	ApplyOperation(ToFunctor(std::fabs), a);
	return Median(a);
}

double MedianStdDev(const operaVector& a) {
	return MedianStdDev(a, Median(a));
}

double MedianStdDev(const operaVector& a, double median) {
	return MedianAbsDev(a, median) / 0.674433; //magic number to convert from absolute deviation to standard deviation, assuming a normal distribution
}

std::pointer_to_unary_function<double, double> ToFunctor(double (*funct)(double)) {
	return std::ptr_fun(funct);
}

std::pointer_to_binary_function<double, double, double> ToFunctor(double (*funct)(double, double)) {
	return std::ptr_fun(funct);
}

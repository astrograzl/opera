#include "globaldefines.h"
#include "libraries/operaVector.h"
#include "libraries/operaException.h"
#include <algorithm>

class operaIndexComparator {
private:
    const std::vector<double> &data;
public:
    operaIndexComparator(const std::vector<double> &data) : data(data) {}
    bool operator()(unsigned i, unsigned j) { return data[i] < data[j]; }
};

operaIndexMap::operaIndexMap(unsigned size) : index(size) {
	for(unsigned i = 0; i < size; i++) index[i] = i;
}

operaIndexRange::operaIndexRange(unsigned startindex, unsigned endindex) : start(startindex), end(endindex) { }

unsigned operaIndexRange::size() {
	return end - start;
}

operaVector::operaVector() { }

operaVector::operaVector(unsigned length) : data(length) { }

operaVector::operaVector(double* dataarray, unsigned length) : data(dataarray, dataarray + length) { }

double& operaVector::operator[](unsigned i) {
#ifdef RANGE_CHECK
    if (i >= data.size()) {
		throw operaException("operaVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data[i];
}

const double& operaVector::operator[](unsigned i) const {
#ifdef RANGE_CHECK
    if (i >= data.size()) {
		throw operaException("operaVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data[i];
}

std::vector<double>::iterator operaVector::begin() {
	return data.begin();
}
	
std::vector<double>::const_iterator operaVector::begin() const {
	return data.begin();
}
	
std::vector<double>::iterator operaVector::end() {
	return data.end();
}
	
std::vector<double>::const_iterator operaVector::end() const {
	return data.end();
}

double operaVector::first() const {
#ifdef RANGE_CHECK
    if (data.empty()) {
		throw operaException("operaVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data.front();
}

double operaVector::last() const {
#ifdef RANGE_CHECK
    if (data.empty()) {
		throw operaException("operaVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	return data.back();
}

unsigned operaVector::size() const {
	return data.size();
}

bool operaVector::empty() const {
	return data.empty();
}

void operaVector::clear() {
	data.clear();
}

void operaVector::resize(unsigned newsize) {
	data.resize(newsize);
}

operaIndexRange operaVector::subrange(double min, double max) const {
	return operaIndexRange(std::lower_bound(data.begin(), data.end(), min) - data.begin(), std::upper_bound(data.begin(), data.end(), max) - data.begin());
}

void operaVector::trim(operaIndexRange range) {
#ifdef RANGE_CHECK
    if (range.start >= data.size() || range.end > data.size() || range.start > range.end) {
		throw operaException("operaVector: ", operaErrorIndexOutOfRange, __FILE__, __FUNCTION__, __LINE__);	
	}
#endif
	std::vector<double> newdata(data.begin() + range.start, data.begin() + range.end);
	data.swap(newdata);
}

void operaVector::fill(double value) {
	std::fill(data.begin(), data.end(), value);
}

void operaVector::insert(double newdata) {
	data.push_back(newdata);
}

void operaVector::reverse() {
	std::reverse(data.begin(), data.end());
}

void operaVector::sort() {
	std::stable_sort(data.begin(), data.end());
}

operaIndexMap operaVector::indexsort() const {
	operaIndexMap indexmap(data.size());
	std::sort(indexmap.index.begin(), indexmap.index.end(), operaIndexComparator(data));
	return indexmap;
}

void operaVector::reorder(const operaIndexMap &indexmap) {
	if(indexmap.index.size() != data.size()) throw operaException("operaVector: ", operaErrorLengthMismatch, __FILE__, __FUNCTION__, __LINE__);
	std::vector<double> newdata(data.size());
	for(unsigned i = 0; i < newdata.size(); i++) newdata[i] = data[indexmap.index[i]];
	data.swap(newdata);
}

void operaVector::copyfrom(double* dataarray) {
	std::copy(dataarray, dataarray + data.size(), data.begin());
}

double* operaVector::datapointer() {
	return &data[0];
}

const double* operaVector::datapointer() const {
	return &data[0];
}

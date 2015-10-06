#ifndef OPERAIOFORMATS_H
#define OPERAIOFORMATS_H

#include "libraries/operaSpectralOrderVector.h"

namespace operaIOFormats {
	/*! 
	 * \sa method void WriteSpectralOrder(operaSpectralOrderVector& orders, string filename, operaSpectralOrder_t format);
	 * \details Writes the specified format information from a spectral order vector to a file
	 * \details the optional order argument permits incremental addition to the output file, where zero means write all.
	 * \return none.
	 */
	void WriteFromSpectralOrders(const operaSpectralOrderVector& orders, string filename, operaSpectralOrder_t format);
	
	/*! 
	 * \sa method void ReadIntoSpectralOrders(operaSpectralOrderVector& orders, string filename);
	 * \brief augment an existing spectral order vector with information from a file
	 * \param filename - string.
	 * \return none.
	 */
	void ReadIntoSpectralOrders(operaSpectralOrderVector& orders, string filename);
}

#endif

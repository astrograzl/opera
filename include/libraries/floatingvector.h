#ifndef FLOATINGVECTOR_H
#define FLOATINGVECTOR_H

/*******************************************************************
 ****               		OPERA PIPELINE v1.0                 ****
 *******************************************************************
 Library name: floatingvector
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

// $Date$
// $Id$
// $Revision$
// $Locker$
// $Log$

/*! 
 * \brief float and double vectors in C.
 * \file floatingvector.h
 * \ingroup libraries
 */
/*
 * \details These vectors are guaranteed to be null terminated so
 * iterating is fast.
 */

#ifdef __cplusplus
extern "C" {
#endif
	
	typedef float **floatvector;
	typedef double **doublevector;
	
	
	floatvector newfloatvector(unsigned Length) {
		floatvector vector = (floatvector)malloc((Length+1)*sizeof(float *));
		memset(vector, 0, (Length+1)*sizeof(float *));
		return vector;
	}
	void deletefloatvector(floatvector vector) {
		free(vector);
	}
	
	doublevector newdoublevector(unsigned Length) {
		doublevector vector = (doublevector)malloc((Length+1)*sizeof(double *));
		memset(vector, 0, (Length+1)*sizeof(double *));
		return vector;
	}
	void deletedoublevector(doublevector vector) {
		free(vector);
	}
	
#ifdef __cplusplus
}
#endif

#endif



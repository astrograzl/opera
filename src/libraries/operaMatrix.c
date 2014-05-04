/*******************************************************************
 ****                LIBRARY FOR OPERA v1.0                     ****
 *******************************************************************
 Library name: operaMatrix
 Version: 1.0
 Description: This C library implements matrix routines..
 Author(s): CFHT OPERA team
 Affiliation: Canada France Hawaii Telescope 
 Location: Hawaii USA
 Date: Aug/2011
 Contact: eder@cfht.hawaii.edu
 
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
#include "libraries/operaMatrix.h"

/* 
 * operaMatrix
 * \author Eder Martioli
 * \brief operaMatrix library.
 * \details {This library contains the basic routines for matrix tools.}
 * \file operaMatrix.c
 * \ingroup libraries
 */


/**********************************************************************************/
/**********************************************************************************/
/**** NOTE WELL:                                                               ****/
/**** Theses functions return pointers to CMatrices created inside the         ****/
/**** functions. The caller is responsible for calling deleteCMatrix() when    ****/
/**** finished with the returned CMtarix.                                      ****/
/**********************************************************************************/
/**********************************************************************************/


/* 
 * void printMatrix(CMatrix inputmatrix)
 * \brief  This function prints the values of a matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return void
 */
void printMatrix(CMatrix matrix){
	for(unsigned j=0; j<getCMatrixRows(matrix); j++) {		
		for(unsigned i=0; i<getCMatrixCols(matrix); i++) {
			printf("%6.3f\t",matrix[j][i]); 
		}
		printf("\n");
	}	
}

/* 
 * float MatrixDeterminant(CMatrix inputmatrix)
 * \brief  This function calculates the value of determinant of a square matrix 
 * \brief  it applies the recursive definition of determinate using expansion by minors.
 * \param  inputmatrix is a CMatrix type  
 * \return float value for the determinant  
 * \notes Source: http://paulbourke.net/miscellaneous/determinant/
 */
float MatrixDeterminant(CMatrix inputmatrix)
{
	
	float det;
	
	unsigned NXPoints = getCMatrixCols(inputmatrix);	
	unsigned NYPoints = getCMatrixRows(inputmatrix);
	
	if(NXPoints != NYPoints) {
		operaPError("operaMatrix:MatrixDeterminant ", MatrixNotSquare);
		return FP_NAN;
	}
	
	unsigned n = NXPoints;

	unsigned i,j,j1,j2;
	
	det = 0;
	
	if (n == 1) {
		
		det = inputmatrix[0][0];
		
	} else if (n == 2)  {
		
		det = (inputmatrix[0][0] * inputmatrix[1][1]) - (inputmatrix[0][1] * inputmatrix[1][0]);		

	} else {
		
		det = 0;		
		
		for (j1 = 0 ; j1 < n ; j1++) {
			
			CMatrix minormatrix = newCMatrix((n-1),(n-1));
			
			for (i = 1 ; i < n ; i++) {
				j2 = 0 ;              
				for (j = 0 ; j < n ; j++) {
					if (j == j1) continue;
					minormatrix[j2][i-1] = inputmatrix[j][i];
					j2++;
				}
			}
			
			det += pow(-1.0,1.0 + j1 + 1.0) * inputmatrix[j1][0] * MatrixDeterminant(minormatrix);
			
			deleteCMatrix(minormatrix);
		}
	}
	
	return(det);
}

 /* 
 * CMatrix MatrixTranspose(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the transpose of a given matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for transpose matrix
 */ 
CMatrix MatrixTranspose(CMatrix inputmatrix, CMatrix outputmatrix) {
	
	if(getCMatrixCols(inputmatrix) != getCMatrixRows(outputmatrix)) {
		operaPError("operaMatrix:MatrixTranspose ", MatrixInvalidDimensions);
		return NULL;
	}
	if(getCMatrixCols(outputmatrix) != getCMatrixRows(inputmatrix)) {
		operaPError("operaMatrix:MatrixTranspose ", MatrixInvalidDimensions);
		return NULL;
	}
	unsigned NXPoints = getCMatrixCols(outputmatrix);	
	unsigned NYPoints = getCMatrixRows(outputmatrix);

	for (unsigned j=0;j<NYPoints;j++) {
		for (unsigned i=0;i<NXPoints;i++) {
			outputmatrix[i][j] = inputmatrix[j][i];
		}
	}		
	
	return outputmatrix;
}

/* 
 * float MatrixTrace(CMatrix inputmatrix)
 * \brief  This function calculates the trace of a square matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return float value for the trace  
 * \notes Input matrix must be square 
 */
float MatrixTrace(CMatrix inputmatrix)
{
	float tracevalue;
	
	unsigned NXPoints = getCMatrixCols(inputmatrix);	
	unsigned NYPoints = getCMatrixRows(inputmatrix);
	
	if(NXPoints != NYPoints) {
		operaPError("operaMatrix:MatrixTrace ", MatrixNotSquare);
		return FP_NAN;
	}
	
	tracevalue=0;
	
	for (unsigned j=0;j<NYPoints;j++) {
		for (unsigned i=0;i<NXPoints;i++) {
			if(i==j) {
				tracevalue += inputmatrix[j][i];
			}
		}
	}
	
	return(tracevalue);
}


/* 
 * CMatrix MatrixCofactor(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the cofactor matrix
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for cofactor matrix
 * \notes Input matrix must be square  
 */ 

CMatrix MatrixCofactor(CMatrix inputmatrix, CMatrix outputmatrix) {
	
	unsigned NXPoints = getCMatrixCols(inputmatrix);	
	unsigned NYPoints = getCMatrixRows(inputmatrix);	

	if(NXPoints != NYPoints) {
		operaPError("operaMatrix:MatrixCofactor ", MatrixNotSquare);
		return NULL;
	}
	
	for (unsigned j=0;j<NYPoints;j++) {
		for (unsigned i=0;i<NXPoints;i++) {
			
			CMatrix minormatrix = newCMatrix(NXPoints-1,NYPoints-1);
			if (minormatrix) {
				/* Form the adjoint a_ij */
				unsigned i1 = 0;
				for (unsigned ii=0;ii<NXPoints;ii++) {
					if (ii == i)
						continue;
					unsigned j1 = 0;
					for (unsigned jj=0;jj<NYPoints;jj++) {
						if (jj == j)
							continue;
						minormatrix[j1][i1] = inputmatrix[jj][ii];
						j1++;
					}
					i1++;
				}
				
				/* Calculate the determinate */
				float det = MatrixDeterminant(minormatrix);
				
				/* Fill in the elements of the cofactor */
				outputmatrix[j][i] = pow(-1.0,i+j+2.0) * det;
				
				deleteCMatrix(minormatrix);					
			}
		}
	}
	return outputmatrix;
}


/* 
 * CMatrix MatrixAdjoint(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the adjoint matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for adjoint matrix
 */ 
CMatrix MatrixAdjoint(CMatrix inputmatrix, CMatrix outputmatrix) {
	return MatrixTranspose(MatrixCofactor(inputmatrix, outputmatrix), outputmatrix);
}

/* 
 * CMatrix MatrixAdjoint(CMatrix inputmatrix, CMatrix outputmatrix)
 * \brief  This function calculates the adjoint matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for adjoint matrix
 */ 
CMatrix MatrixMultiplication(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix) {
	
	unsigned NCOLS1 = getCMatrixCols(inputmatrix1);	
	unsigned NROWS1 = getCMatrixRows(inputmatrix1);	

	unsigned NCOLS2 = getCMatrixCols(inputmatrix2);	
	unsigned NROWS2 = getCMatrixRows(inputmatrix2);		
	
	unsigned NCOLS3 = getCMatrixCols(outputmatrix);	
	unsigned NROWS3 = getCMatrixRows(outputmatrix);		
	
	if(NCOLS1 != NROWS2) {
		operaPError("operaMatrix:MatrixMultiplication ", MatrixInvalidDimensions);
		return NULL;
	}	
	unsigned NCOMMON = NCOLS1;
	
	unsigned NCOLSRES = NCOLS2;	
	unsigned NROWSRES = NROWS1;	
	
	if(NCOLSRES != NCOLS3 || NROWSRES != NROWS3) {
		operaPError("operaMatrix:MatrixMultiplication ", MatrixInvalidDimensions);
		return NULL;
	}	
	for (unsigned j=0;j<NROWSRES;j++) {
		for (unsigned i=0;i<NCOLSRES;i++) {
			outputmatrix[j][i] = 0;
			for (unsigned k=0;k<NCOMMON;k++) {
				outputmatrix[j][i] += inputmatrix1[j][k]*inputmatrix2[k][i];  
			}
		}
	}				
	return outputmatrix;
}

/* 
 * CMatrix MatrixMultiplicationbyConstant(CMatrix inputmatrix, float constantValue, CMatrix outputmatrix)
 * \brief  This function multiply a matrix by a constant 
 * \param  inputmatrix is a CMatrix type 
 * \param  float constantValue  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixMultiplicationbyConstant(CMatrix inputmatrix, float constantValue, CMatrix outputmatrix) {
	
	unsigned NCOLS = getCMatrixCols(inputmatrix);	
	unsigned NROWS = getCMatrixRows(inputmatrix);	
	
	for (unsigned j=0;j<NROWS;j++) {
		for (unsigned i=0;i<NCOLS;i++) {
			outputmatrix[j][i] = inputmatrix[j][i]*constantValue;
		}
	}				
	return outputmatrix;
}

/* 
 * CMatrix MatrixAddition(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix)
 * \brief  This function calculates the sum between two matrices 
 * \param  inputmatrix1 and inputmatrix2 are of CMatrix type  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixAddition(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix) {
	
	unsigned NCOLS1 = getCMatrixCols(inputmatrix1);	
	unsigned NROWS1 = getCMatrixRows(inputmatrix1);	
	
	unsigned NCOLS2 = getCMatrixCols(inputmatrix2);	
	unsigned NROWS2 = getCMatrixRows(inputmatrix2);		
	
	if(NCOLS1 != NCOLS2 || NROWS1 != NROWS2) {
		operaPError("operaMatrix:MatrixAddition ", MatrixInvalidDimensions);
		return NULL;
	}	
	
	unsigned NCOLSRES = NCOLS1;	
	unsigned NROWSRES = NROWS1;	
	
	for (unsigned j=0;j<NROWSRES;j++) {
		for (unsigned i=0;i<NCOLSRES;i++) {
			outputmatrix[j][i] = inputmatrix1[j][i] + inputmatrix2[j][i];
		}
	}	
	return outputmatrix;
}

/* 
 * CMatrix MatrixSubtraction(CMatrix inputmatrix1, CMatrix inputmatrix2)
 * \brief  This function calculates the subtraction between two matrices 
 * \param  inputmatrix1 and inputmatrix2 are of CMatrix type  
 * \return CMatrix for result matrix
 */ 
CMatrix MatrixSubtraction(CMatrix inputmatrix1, CMatrix inputmatrix2, CMatrix outputmatrix) {
	
	unsigned NCOLS1 = getCMatrixCols(inputmatrix1);	
	unsigned NROWS1 = getCMatrixRows(inputmatrix1);	
	
	unsigned NCOLS2 = getCMatrixCols(inputmatrix2);	
	unsigned NROWS2 = getCMatrixRows(inputmatrix2);		
	
	unsigned NCOLS3 = getCMatrixCols(outputmatrix);	
	unsigned NROWS3 = getCMatrixRows(outputmatrix);		
	
	if(NCOLS1 != NCOLS2 || NROWS1 != NROWS2 || NROWS1 != NROWS3 || NCOLS1 != NCOLS3) {
		operaPError("operaMatrix:MatrixSubtraction ", MatrixInvalidDimensions);
		return NULL;
	}	
	
	for (unsigned j=0;j<NROWS3;j++) {
		for (unsigned i=0;i<NCOLS3;i++) {
			outputmatrix[j][i] = inputmatrix1[j][i] - inputmatrix2[j][i];
		}
	}		
	return outputmatrix;
}

/* 
 * CMatrix MatrixInverse(CMatrix inputmatrix)
 * \brief  This function calculates the inverse matrix 
 * \param  inputmatrix is a CMatrix type  
 * \return CMatrix for output inverse matrix
 */ 
CMatrix MatrixInverse(CMatrix inputmatrix, CMatrix outputmatrix) {
	
	float det = MatrixDeterminant(inputmatrix);
	
	if(det == 0) {
		operaPError("operaMatrix:MatrixInverse ", MatrixZeroDeterminant);
		return NULL;
	}	
	
	return MatrixMultiplicationbyConstant(MatrixAdjoint(inputmatrix, outputmatrix), (1.0/det), outputmatrix);
}


/* 
 * CMatrix RotationMatrix2D(float angleInDegrees, CMatrix outputmatrix) {
 * \brief  This function produces an operator matrix to rotate 
 * \brief  points in the xy-Cartesian plane counterclockwise through 
 * \brief  a given angle about the origin of the Cartesian coordinate system. 
 * \param  float angleInDegrees is the rotation angle in degrees. 
 * \return CMatrix (2X2) for rotation matrix
 */ 
CMatrix RotationMatrix2D(float angleInDegrees, CMatrix outputmatrix) {
    
    unsigned NCOLS = 2;
    unsigned NROWS = 2;
    
	unsigned NCOLS2 = getCMatrixCols(outputmatrix);	
	unsigned NROWS2 = getCMatrixRows(outputmatrix);		

	if(NCOLS != NCOLS2 || NROWS != NROWS2) {
		operaPError("operaMatrix:MatrixSubtraction ", MatrixInvalidDimensions);
		return NULL;
	}	
    float angleInRadians = angleInDegrees*M_PI/180.0;
    
    outputmatrix[0][0] = cos(angleInRadians);
    outputmatrix[0][1] = -sin(angleInRadians);
    outputmatrix[1][0] = sin(angleInRadians);
    outputmatrix[1][1] = cos(angleInRadians); 
	
    return outputmatrix;
}


/*
printMatrix
MatrixDeterminant
MatrixTranspose 	
MatrixTrace
MatrixCofactor
MatrixAdjoint
MatrixMultiplication
MatrixMultiplicationbyConstant
MatrixAddition
MatrixSubtraction
MatrixInverse 	
 
CholeskyDecomposition
Moore-PenroseInverse
QRDecomposition

LUDecomposition 	
SingularValueDecomposition 	
SystemsofLinearEquations
MatrixRank 	
LQDecomposition
Eigenvalues
Eigenvectors
*/

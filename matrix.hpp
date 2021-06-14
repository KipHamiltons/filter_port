// Copyright (c) 2014, Freescale Semiconductor, Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Freescale Semiconductor, Inc. nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL FREESCALE SEMICONDUCTOR, INC. BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// This file contains matrix manipulation functions.
//

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>

#include "build.hpp"

// compile time constants that are private to this file
static constexpr double CORRUPTMATRIX = 0.001;  // column vector modulus limit for rotation matrix

namespace filter::matrix {
    // function sets the 3x3 matrix A to the identity matrix
    template <typename Scalar>
    inline void f3x3matrixAeqI(Scalar A[][3]) {
        Scalar* pAij;  // pointer to A[i][j]
        int i, j;      // loop counters

        for (i = 0; i < 3; i++) {
            // set pAij to &A[i][j=0]
            pAij = A[i];
            for (j = 0; j < 3; j++) {
                *(pAij++) = 0.0F;
            }
            A[i][i] = 1.0F;
        }
        return;
    }

    // function sets the matrix A to the identity matrix
    template <typename Scalar>
    inline void fmatrixAeqI(Scalar* A[], int rc) {
        // rc = rows and columns in A

        Scalar* pAij;  // pointer to A[i][j]
        int i, j;      // loop counters

        for (i = 0; i < rc; i++) {
            // set pAij to &A[i][j=0]
            pAij = A[i];
            for (j = 0; j < rc; j++) {
                *(pAij++) = 0.0F;
            }
            A[i][i] = 1.0F;
        }
        return;
    }

    // function sets every entry in the 3x3 matrix A to a constant scalar
    template <typename Scalar>
    inline void f3x3matrixAeqScalar(Scalar A[][3], Scalar scalar) {
        Scalar* pAij;  // pointer to A[i][j]
        int i, j;      // counters

        for (i = 0; i < 3; i++) {
            // set pAij to &A[i][j=0]
            pAij = A[i];
            for (j = 0; j < 3; j++) {
                *(pAij++) = scalar;
            }
        }
        return;
    }

    // function uses Gauss-Jordan elimination to compute the inverse of matrix A in situ
    // on exit, A is replaced with its inverse
    template <typename Scalar>
    inline void fmatrixAeqInvA(Scalar* A[], int iColInd[], int iRowInd[], int iPivot[], int isize) {
        Scalar largest;            // largest element used for pivoting
        Scalar scaling;            // scaling factor in pivoting
        Scalar recippiv;           // reciprocal of pivot element
        Scalar ftmp;               // temporary variable used in swaps
        int i, j, k, l, m;         // index counters
        int iPivotRow, iPivotCol;  // row and column of pivot element

        // to ainline void compiler warnings
        iPivotRow = iPivotCol = 0;

        // initialize the pivot array to 0
        for (j = 0; j < isize; j++) {
            iPivot[j] = 0;
        }

        // main loop i over the dimensions of the square matrix A
        for (i = 0; i < isize; i++) {
            // zero the largest element found for pivoting
            largest = 0.0F;
            // loop over candidate rows j
            for (j = 0; j < isize; j++) {
                // check if row j has been previously pivoted
                if (iPivot[j] != 1) {
                    // loop over candidate columns k
                    for (k = 0; k < isize; k++) {
                        // check if column k has previously been pivoted
                        if (iPivot[k] == 0) {
                            // check if the pivot element is the largest found so far
                            if (std::fabs(A[j][k]) >= largest) {
                                // and store this location as the current best candidate for pivoting
                                iPivotRow = j;
                                iPivotCol = k;
                                largest   = (Scalar) std::fabs(A[iPivotRow][iPivotCol]);
                            }
                        }
                        else if (iPivot[k] > 1) {
                            // zero determinant situation: exit with identity matrix
                            fmatrixAeqI(A, isize);
                            return;
                        }
                    }
                }
            }
            // increment the entry in iPivot to denote it has been selected for pivoting
            iPivot[iPivotCol]++;

            // check the pivot rows iPivotRow and iPivotCol are not the same before swapping
            if (iPivotRow != iPivotCol) {
                // loop over columns l
                for (l = 0; l < isize; l++) {
                    // and swap all elements of rows iPivotRow and iPivotCol
                    ftmp            = A[iPivotRow][l];
                    A[iPivotRow][l] = A[iPivotCol][l];
                    A[iPivotCol][l] = ftmp;
                }
            }

            // record that on the i-th iteration rows iPivotRow and iPivotCol were swapped
            iRowInd[i] = iPivotRow;
            iColInd[i] = iPivotCol;

            // check for zero on-diagonal element (singular matrix) and return with identity matrix if detected
            if (A[iPivotCol][iPivotCol] == 0.0F) {
                // zero determinant situation: exit with identity matrix
                fmatrixAeqI(A, isize);
                return;
            }

            // calculate the reciprocal of the pivot element knowing it's non-zero
            recippiv = 1.0F / A[iPivotCol][iPivotCol];
            // by definition, the diagonal element normalizes to 1
            A[iPivotCol][iPivotCol] = 1.0F;
            // multiply all of row iPivotCol by the reciprocal of the pivot element including the diagonal element
            // the diagonal element A[iPivotCol][iPivotCol] now has value equal to the reciprocal of its previous value
            for (l = 0; l < isize; l++) {
                A[iPivotCol][l] *= recippiv;
            }
            // loop over all rows m of A
            for (m = 0; m < isize; m++) {
                if (m != iPivotCol) {
                    // scaling factor for this row m is in column iPivotCol
                    scaling = A[m][iPivotCol];
                    // zero this element
                    A[m][iPivotCol] = 0.0F;
                    // loop over all columns l of A and perform elimination
                    for (l = 0; l < isize; l++) {
                        A[m][l] -= A[iPivotCol][l] * scaling;
                    }
                }
            }
        }  // end of loop i over the matrix dimensions

        // finally, loop in inverse order to apply the missing column swaps
        for (l = isize - 1; l >= 0; l--) {
            // set i and j to the two columns to be swapped
            i = iRowInd[l];
            j = iColInd[l];

            // check that the two columns i and j to be swapped are not the same
            if (i != j) {
                // loop over all rows k to swap columns i and j of A
                for (k = 0; k < isize; k++) {
                    ftmp    = A[k][i];
                    A[k][i] = A[k][j];
                    A[k][j] = ftmp;
                }
            }
        }

        return;
    }

}  // namespace filter::matrix

#endif  // MATRIX_HPP
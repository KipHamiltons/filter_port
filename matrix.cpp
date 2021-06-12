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
#include "matrix.hpp"

#include <cmath>

#include "build.hpp"
#include "math.h"
// #include "stdio.h"
// #include "stdlib.h"
// #include "string.h"
// #include "time.h"

// compile time constants that are private to this file
#define CORRUPTMATRIX 0.001F  // column vector modulus limit for rotation matrix

namespace filter::matrix {

    // function uses Gauss-Jordan elimination to compute the inverse of matrix A in situ
    // on exit, A is replaced with its inverse
    void fmatrixAeqInvA(double* A[], int iColInd[], int iRowInd[], int iPivot[], int isize) {
        double largest;            // largest element used for pivoting
        double scaling;            // scaling factor in pivoting
        double recippiv;           // reciprocal of pivot element
        double ftmp;               // temporary variable used in swaps
        int i, j, k, l, m;         // index counters
        int iPivotRow, iPivotCol;  // row and column of pivot element

        // to avoid compiler warnings
        iPivotRow = iPivotCol = 0;

        // initialize the pivot array to 0
        for (j = 0; j < isize; j++) {
            iPivot[j] = 0;
        }

        // main loop i over the dimensions of the square matrix A
        for (i = 0; i < isize; i++) {
            // zero the largest element found for pivoting
            largest = 0.0;
            // loop over candidate rows j
            for (j = 0; j < isize; j++) {
                // check if row j has been previously pivoted
                if (iPivot[j] != 1) {
                    // loop over candidate columns k
                    for (k = 0; k < isize; k++) {
                        // check if column k has previously been pivoted
                        if (iPivot[k] == 0) {
                            // check if the pivot element is the largest found so far
                            if (fabs(A[j][k]) >= largest) {
                                // and store this location as the current best candidate for pivoting
                                iPivotRow = j;
                                iPivotCol = k;
                                largest   = std::fabs(A[iPivotRow][iPivotCol]);
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
            if (A[iPivotCol][iPivotCol] == 0.0) {
                // zero determinant situation: exit with identity matrix
                fmatrixAeqI(A, isize);
                return;
            }

            // calculate the reciprocal of the pivot element knowing it's non-zero
            recippiv = 1.0 / A[iPivotCol][iPivotCol];
            // by definition, the diagonal element normalizes to 1
            A[iPivotCol][iPivotCol] = 1.0;
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
                    A[m][iPivotCol] = 0.0;
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


    // function sets the matrix A to the identity matrix
    void fmatrixAeqI(double* A[], int rc) {
        // rc = rows and columns in A

        double* pAij;  // pointer to A[i][j]
        int i, j;      // loop counters

        for (i = 0; i < rc; i++) {
            // set pAij to &A[i][j=0]
            pAij = A[i];
            for (j = 0; j < rc; j++) {
                *(pAij++) = 0.0;
            }
            A[i][i] = 1.0;
        }
        return;
    }

    // function sets the 3x3 matrix A to the identity matrix
    void f3x3matrixAeqI(Eigen::Matrix<double, 3, 3>& A) {
        A = Eigen::Matrix<double, 3, 3>::Identity();
    }

    // function multiplies all elements of 3x3 matrix A by the specified scalar
    void f3x3matrixAeqAxScalar(Eigen::Matrix<double, 3, 3>& A, double Scalar) {
        A = Scalar * A;
    }

    // function negates all elements of 3x3 matrix A
    void f3x3matrixAeqMinusA(Eigen::Matrix<double, 3, 3>& A) {
        A = -A;
    }

    // function directly calculates the symmetric inverse of a symmetric 3x3 matrix
    // only the on and above diagonal terms in B are used and need to be specified
    void f3x3matrixAeqInvSymB(Eigen::Matrix<double, 3, 3>& A, Eigen::Matrix<double, 3, 3>& B) {
        double fB11B22mB12B12 = 0.0;  // B(1,1) * B(2,2) - B(1,2) * B(1,2)
        double fB12B02mB01B22 = 0.0;  // B(1,2) * B(0,2) - B(0,1) * B(2,2)
        double fB01B12mB11B02 = 0.0;  // B(0,1) * B(1,2) - B(1,1) * B(0,2)
        double ftmp           = 0.0;  // determinant and then reciprocal

        // calculate useful products
        fB11B22mB12B12 = B(1, 1) * B(2, 2) - B(1, 2) * B(1, 2);
        fB12B02mB01B22 = B(1, 2) * B(0, 2) - B(0, 1) * B(2, 2);
        fB01B12mB11B02 = B(0, 1) * B(1, 2) - B(1, 1) * B(0, 2);

        // set ftmp to the determinant of the input matrix B
        ftmp = B(0, 0) * fB11B22mB12B12 + B(0, 1) * fB12B02mB01B22 + B(0, 2) * fB01B12mB11B02;

        // set A to the inverse of B for any determinant except zero
        if (ftmp != 0.0F) {
            ftmp    = 1.0F / ftmp;
            A(0, 0) = fB11B22mB12B12 * ftmp;
            A(1, 0) = A(0, 1) = fB12B02mB01B22 * ftmp;
            A(2, 0) = A(0, 2) = fB01B12mB11B02 * ftmp;
            A(1, 1)           = (B(0, 0) * B(2, 2) - B(0, 2) * B(0, 2)) * ftmp;
            A(2, 1) = A(1, 2) = (B(0, 2) * B(0, 1) - B(0, 0) * B(1, 2)) * ftmp;
            A(2, 2)           = (B(0, 0) * B(1, 1) - B(0, 1) * B(0, 1)) * ftmp;
        }
        else {
            // provide the identity matrix if the determinant is zero
            f3x3matrixAeqI(A);
        }
    }

    // function re-orthonormalizes a 3x3 rotation matrix
    void fmatrixAeqRenormRotA(Eigen::Matrix<double, 3, 3>& A) {
        double ftmp;  // scratch variable

        // normalize the X column of the low pass filtered orientation matrix
        ftmp = std::sqrt(A(0, 0) * A(0, 0) + A(1, 0) * A(1, 0) + A(2, 0) * A(2, 0));
        if (ftmp > CORRUPTMATRIX) {
            // normalize the x column vector
            ftmp = 1.0F / ftmp;
            A(0, 0) *= ftmp;
            A(1, 0) *= ftmp;
            A(2, 0) *= ftmp;
        }
        else {
            // set x column vector to {1, 0, 0}
            A(0, 0) = 1.0F;
            A(1, 0) = A(2, 0) = 0.0F;
        }

        // force the y column vector to be orthogonal to x using y = y-(x.y)x
        ftmp = A(0, 0) * A(0, 1) + A(1, 0) * A(1, 1) + A(2, 0) * A(2, 1);
        A(0, 1) -= ftmp * A(0, 0);
        A(1, 1) -= ftmp * A(1, 0);
        A(2, 1) -= ftmp * A(2, 0);

        // normalize the y column vector
        ftmp = std::sqrt(A(0, 1) * A(0, 1) + A(1, 1) * A(1, 1) + A(2, 1) * A(2, 1));
        if (ftmp > CORRUPTMATRIX) {
            // normalize the y column vector
            ftmp = 1.0F / ftmp;
            A(0, 1) *= ftmp;
            A(1, 1) *= ftmp;
            A(2, 1) *= ftmp;
        }
        else {
            // set y column vector to {0, 1, 0}
            A(1, 1) = 1.0F;
            A(0, 1) = A(2, 1) = 0.0F;
        }

        // finally set the z column vector to x vector cross y vector (automatically normalized)
        A(0, 2) = A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1);
        A(1, 2) = A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1);
        A(2, 2) = A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1);

        return;
    }
}  // namespace filter::matrix
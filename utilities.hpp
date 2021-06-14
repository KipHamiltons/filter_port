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
// This file contains functions designed to operate on, or compute, orientations.
// These may be in rotation matrix form, quaternion form, or Euler angles.
// It also includes functions designed to operate with specify reference frames
// (Android, Windows 8, NED).
//
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>

#include "build.hpp"

// compile time constants that are private to this file
static constexpr double SMALLQ0      = 0.01;   // limit of quaternion scalar component requiring special algorithm
static constexpr double CORRUPTQUAT  = 0.001;  // threshold for deciding rotation quaternion is corrupt
static constexpr double SMALLMODULUS = 0.01;   // limit where rounding errors may appear

namespace filter::utilities {

    template <typename Scalar>
    [[nodiscard]] constexpr Scalar fasin_deg(Scalar x) {
        const Scalar asin_rad = std::asin(x);
        const Scalar asin_deg = asin_rad * 180.0f * M_1_PI;
        return asin_deg;
    }


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

    // function normalizes a rotation quaternion and ensures q0 is non-negative
    template <typename Scalar>
    inline void fqAeqNormqA(struct fquaternion* pqA) {
        Scalar fNorm;  // quaternion Norm

        // calculate the quaternion Norm
        fNorm = std::sqrt(pqA->q0 * pqA->q0 + pqA->q1 * pqA->q1 + pqA->q2 * pqA->q2 + pqA->q3 * pqA->q3);
        if (fNorm > CORRUPTQUAT) {
            // general case
            fNorm = 1.0 / fNorm;
            pqA->q0 *= fNorm;
            pqA->q1 *= fNorm;
            pqA->q2 *= fNorm;
            pqA->q3 *= fNorm;
        }
        else {
            // return with identity quaternion since the quaternion is corrupted
            pqA->q0 = 1.0;
            pqA->q1 = pqA->q2 = pqA->q3 = 0.0;
        }

        // correct a negative scalar component if the function was called with negative q0
        if (pqA->q0 < 0.0) {
            pqA->q0 = -pqA->q0;
            pqA->q1 = -pqA->q1;
            pqA->q2 = -pqA->q2;
            pqA->q3 = -pqA->q3;
        }

        return;
    }

    // set a quaternion to the unit quaternion
    inline void fqAeq1(struct fquaternion* pqA) {
        pqA->q0 = 1.0;
        pqA->q1 = pqA->q2 = pqA->q3 = 0.0;

        return;
    }

    // function compute the quaternion product qA * qB
    inline void qAeqBxC(struct fquaternion* pqA, const struct fquaternion* pqB, const struct fquaternion* pqC) {
        pqA->q0 = pqB->q0 * pqC->q0 - pqB->q1 * pqC->q1 - pqB->q2 * pqC->q2 - pqB->q3 * pqC->q3;
        pqA->q1 = pqB->q0 * pqC->q1 + pqB->q1 * pqC->q0 + pqB->q2 * pqC->q3 - pqB->q3 * pqC->q2;
        pqA->q2 = pqB->q0 * pqC->q2 - pqB->q1 * pqC->q3 + pqB->q2 * pqC->q0 + pqB->q3 * pqC->q1;
        pqA->q3 = pqB->q0 * pqC->q3 + pqB->q1 * pqC->q2 - pqB->q2 * pqC->q1 + pqB->q3 * pqC->q0;

        return;
    }

    // function compute the quaternion product qA = qA * qB
    inline void qAeqAxB(struct fquaternion* pqA, const struct fquaternion* pqB) {
        struct fquaternion qProd;

        // perform the quaternion product
        qProd.q0 = pqA->q0 * pqB->q0 - pqA->q1 * pqB->q1 - pqA->q2 * pqB->q2 - pqA->q3 * pqB->q3;
        qProd.q1 = pqA->q0 * pqB->q1 + pqA->q1 * pqB->q0 + pqA->q2 * pqB->q3 - pqA->q3 * pqB->q2;
        qProd.q2 = pqA->q0 * pqB->q2 - pqA->q1 * pqB->q3 + pqA->q2 * pqB->q0 + pqA->q3 * pqB->q1;
        qProd.q3 = pqA->q0 * pqB->q3 + pqA->q1 * pqB->q2 - pqA->q2 * pqB->q1 + pqA->q3 * pqB->q0;

        // copy the result back into qA
        *pqA = qProd;

        return;
    }

    // function compute the quaternion product conjg(qA) * qB
    inline struct fquaternion qconjgAxB(const struct fquaternion* pqA, const struct fquaternion* pqB) {
        struct fquaternion qProd;

        qProd.q0 = pqA->q0 * pqB->q0 + pqA->q1 * pqB->q1 + pqA->q2 * pqB->q2 + pqA->q3 * pqB->q3;
        qProd.q1 = pqA->q0 * pqB->q1 - pqA->q1 * pqB->q0 - pqA->q2 * pqB->q3 + pqA->q3 * pqB->q2;
        qProd.q2 = pqA->q0 * pqB->q2 + pqA->q1 * pqB->q3 - pqA->q2 * pqB->q0 - pqA->q3 * pqB->q1;
        qProd.q3 = pqA->q0 * pqB->q3 - pqA->q1 * pqB->q2 + pqA->q2 * pqB->q1 - pqA->q3 * pqB->q0;

        return qProd;
    }

    // Aerospace NED accelerometer 3DOF tilt function computing rotation matrix fR
    template <typename Scalar>
    inline void f3DOFTiltNED(Scalar fR[][3], Scalar fGp[]) {
        // the NED self-consistency twist occurs at 90 deg pitch

        // local variables
        int i;                 // counter
        Scalar fmodGxyz;       // modulus of the x, y, z accelerometer readings
        Scalar fmodGyz;        // modulus of the y, z accelerometer readings
        Scalar frecipmodGxyz;  // reciprocal of modulus
        Scalar ftmp;           // scratch variable

        // compute the accelerometer squared magnitudes
        fmodGyz  = fGp[Y] * fGp[Y] + fGp[Z] * fGp[Z];
        fmodGxyz = fmodGyz + fGp[X] * fGp[X];

        // check for freefall special case where no solution is possible
        if (fmodGxyz == 0.0) {
            f3x3matrixAeqI(fR);
            return;
        }

        // check for vertical up or down gimbal lock case
        if (fmodGyz == 0.0) {
            f3x3matrixAeqScalar(fR, 0.0);
            fR[Y][Y] = 1.0;
            if (fGp[X] >= 0.0) {
                fR[X][Z] = 1.0;
                fR[Z][X] = -1.0;
            }
            else {
                fR[X][Z] = -1.0;
                fR[Z][X] = 1.0;
            }
            return;
        }

        // compute moduli for the general case
        fmodGyz       = std::sqrt(fmodGyz);
        fmodGxyz      = std::sqrt(fmodGxyz);
        frecipmodGxyz = 1.0 / fmodGxyz;
        ftmp          = fmodGxyz / fmodGyz;

        // normalize the accelerometer reading into the z column
        for (i = X; i <= Z; i++) {
            fR[i][Z] = fGp[i] * frecipmodGxyz;
        }

        // construct x column of orientation matrix
        fR[X][X] = fmodGyz * frecipmodGxyz;
        fR[Y][X] = -fR[X][Z] * fR[Y][Z] * ftmp;
        fR[Z][X] = -fR[X][Z] * fR[Z][Z] * ftmp;

        // // construct y column of orientation matrix
        fR[X][Y] = 0.0;
        fR[Y][Y] = fR[Z][Z] * ftmp;
        fR[Z][Y] = -fR[Y][Z] * ftmp;

        return;
    }

    // Android accelerometer 3DOF tilt function computing rotation matrix fR
    template <typename Scalar>
    inline void f3DOFTiltAndroid(Scalar fR[][3], Scalar fGp[]) {
        // the Android tilt matrix is mathematically identical to the NED tilt matrix
        // the Android self-consistency twist occurs at 90 deg roll
        f3DOFTiltNED(fR, fGp);
        return;
    }

    // Windows 8 accelerometer 3DOF tilt function computing rotation matrix fR
    template <typename Scalar>
    inline void f3DOFTiltWin8(Scalar fR[][3], Scalar fGp[]) {
        // the Win8 self-consistency twist occurs at 90 deg roll

        // local variables
        Scalar fmodGxyz;       // modulus of the x, y, z accelerometer readings
        Scalar fmodGxz;        // modulus of the x, z accelerometer readings
        Scalar frecipmodGxyz;  // reciprocal of modulus
        Scalar ftmp;           // scratch variable
        int i;                 // counter

        // compute the accelerometer squared magnitudes
        fmodGxz  = fGp[X] * fGp[X] + fGp[Z] * fGp[Z];
        fmodGxyz = fmodGxz + fGp[Y] * fGp[Y];

        // check for freefall special case where no solution is possible
        if (fmodGxyz == 0.0) {
            f3x3matrixAeqI(fR);
            return;
        }

        // check for vertical up or down gimbal lock case
        if (fmodGxz == 0.0) {
            f3x3matrixAeqScalar(fR, 0.0);
            fR[X][X] = 1.0;
            if (fGp[Y] >= 0.0) {
                fR[Y][Z] = -1.0;
                fR[Z][Y] = 1.0;
            }
            else {
                fR[Y][Z] = 1.0;
                fR[Z][Y] = -1.0;
            }
            return;
        }

        // compute moduli for the general case
        fmodGxz       = std::sqrt(fmodGxz);
        fmodGxyz      = std::sqrt(fmodGxyz);
        frecipmodGxyz = 1.0 / fmodGxyz;
        ftmp          = fmodGxyz / fmodGxz;
        if (fGp[Z] < 0.0) {
            ftmp = -ftmp;
        }

        // normalize the negated accelerometer reading into the z column
        for (i = X; i <= Z; i++) {
            fR[i][Z] = -fGp[i] * frecipmodGxyz;
        }

        // construct x column of orientation matrix
        fR[X][X] = -fR[Z][Z] * ftmp;
        fR[Y][X] = 0.0;
        fR[Z][X] = fR[X][Z] * ftmp;
        ;

        // // construct y column of orientation matrix
        fR[X][Y] = fR[X][Z] * fR[Y][Z] * ftmp;
        fR[Y][Y] = -fmodGxz * frecipmodGxyz;
        if (fGp[Z] < 0.0) {
            fR[Y][Y] = -fR[Y][Y];
        }
        fR[Z][Y] = fR[Y][Z] * fR[Z][Z] * ftmp;

        return;
    }

    // computes normalized rotation quaternion from a rotation vector (deg)
    template <typename Scalar>
    inline void fQuaternionFromRotationVectorDeg(struct fquaternion* pq, const Scalar rvecdeg[], Scalar fscaling) {
        Scalar fetadeg;     // rotation angle (deg)
        Scalar fetarad;     // rotation angle (rad)
        Scalar fetarad2;    // eta (rad)^2
        Scalar fetarad4;    // eta (rad)^4
        Scalar sinhalfeta;  // sin(eta/2)
        Scalar fvecsq;      // q1^2+q2^2+q3^2
        Scalar ftmp;        // scratch variable

        // compute the scaled rotation angle eta (deg) which can be both positve or negative
        fetadeg  = fscaling * std::sqrt(rvecdeg[X] * rvecdeg[X] + rvecdeg[Y] * rvecdeg[Y] + rvecdeg[Z] * rvecdeg[Z]);
        fetarad  = fetadeg * FDEGTORAD;
        fetarad2 = fetarad * fetarad;

        // calculate the sine and cosine using small angle approximations or exact
        // angles under sqrt(0.02)=0.141 rad is 8.1 deg and 1620 deg/s (=936deg/s in 3 axes) at 200Hz and 405 deg/s at
        // 50Hz
        if (fetarad2 <= 0.02) {
            // use MacLaurin series up to and including third order
            sinhalfeta = fetarad * (0.5 - ONEOVER48 * fetarad2);
        }
        else if (fetarad2 <= 0.06) {
            // use MacLaurin series up to and including fifth order
            // angles under sqrt(0.06)=0.245 rad is 14.0 deg and 2807 deg/s (=1623deg/s in 3 axes) at 200Hz and 703
            // deg/s at 50Hz
            fetarad4   = fetarad2 * fetarad2;
            sinhalfeta = fetarad * (0.5 - ONEOVER48 * fetarad2 + ONEOVER3840 * fetarad4);
        }
        else {
            // use exact calculation
            sinhalfeta = (Scalar) std::sin(0.5 * fetarad);
        }

        // compute the vector quaternion components q1, q2, q3
        if (fetadeg != 0.0) {
            // general case with non-zero rotation angle
            ftmp   = fscaling * sinhalfeta / fetadeg;
            pq->q1 = rvecdeg[X] * ftmp;  // q1 = nx * sin(eta/2)
            pq->q2 = rvecdeg[Y] * ftmp;  // q2 = ny * sin(eta/2)
            pq->q3 = rvecdeg[Z] * ftmp;  // q3 = nz * sin(eta/2)
        }
        else {
            // zero rotation angle giving zero vector component
            pq->q1 = pq->q2 = pq->q3 = 0.0;
        }

        // compute the scalar quaternion component q0 by explicit normalization
        // taking care to avoid rounding errors giving negative operand to sqrt
        fvecsq = pq->q1 * pq->q1 + pq->q2 * pq->q2 + pq->q3 * pq->q3;
        if (fvecsq <= 1.0) {
            // normal case
            pq->q0 = std::sqrt(1.0 - fvecsq);
        }
        else {
            // rounding errors are present
            pq->q0 = 0.0;
        }

        return;
    }

    // compute the orientation quaternion from a 3x3 rotation matrix
    template <typename Scalar>
    inline void fQuaternionFromRotationMatrix(Scalar R[][3], struct fquaternion* pq) {
        Scalar fq0sq;     // q0^2
        Scalar recip4q0;  // 1/4q0

        // the quaternion is not explicitly normalized in this function on the assumption that it
        // is supplied with a normalized rotation matrix. if the rotation matrix is normalized then
        // the quaternion will also be normalized even if the case of small q0

        // get q0^2 and q0
        fq0sq  = 0.25 * (1.0 + R[X][X] + R[Y][Y] + R[Z][Z]);
        pq->q0 = std::sqrt(std::fabs(fq0sq));

        // normal case when q0 is not small meaning rotation angle not near 180 deg
        if (pq->q0 > SMALLQ0) {
            // calculate q1 to q3
            recip4q0 = 0.25 / pq->q0;
            pq->q1   = recip4q0 * (R[Y][Z] - R[Z][Y]);
            pq->q2   = recip4q0 * (R[Z][X] - R[X][Z]);
            pq->q3   = recip4q0 * (R[X][Y] - R[Y][X]);
        }  // end of general case
        else {
            // special case of near 180 deg corresponds to nearly symmetric matrix
            // which is not numerically well conditioned for division by small q0
            // instead get absolute values of q1 to q3 from leading diagonal
            pq->q1 = std::sqrt(std::fabs(0.5 * (1.0 + R[X][X]) - fq0sq));
            pq->q2 = std::sqrt(std::fabs(0.5 * (1.0 + R[Y][Y]) - fq0sq));
            pq->q3 = std::sqrt(std::fabs(0.5 * (1.0 + R[Z][Z]) - fq0sq));

            // correct the signs of q1 to q3 by examining the signs of differenced off-diagonal terms
            if ((R[Y][Z] - R[Z][Y]) < 0.0)
                pq->q1 = -pq->q1;
            if ((R[Z][X] - R[X][Z]) < 0.0)
                pq->q2 = -pq->q2;
            if ((R[X][Y] - R[Y][X]) < 0.0)
                pq->q3 = -pq->q3;
        }  // end of special case

        return;
    }

    // compute the rotation matrix from an orientation quaternion
    template <typename Scalar>
    inline void fRotationMatrixFromQuaternion(Scalar R[][3], const struct fquaternion* pq) {
        Scalar f2q;
        Scalar f2q0q0, f2q0q1, f2q0q2, f2q0q3;
        Scalar f2q1q1, f2q1q2, f2q1q3;
        Scalar f2q2q2, f2q2q3;
        Scalar f2q3q3;

        // calculate products
        f2q    = 2.0 * pq->q0;
        f2q0q0 = f2q * pq->q0;
        f2q0q1 = f2q * pq->q1;
        f2q0q2 = f2q * pq->q2;
        f2q0q3 = f2q * pq->q3;
        f2q    = 2.0 * pq->q1;
        f2q1q1 = f2q * pq->q1;
        f2q1q2 = f2q * pq->q2;
        f2q1q3 = f2q * pq->q3;
        f2q    = 2.0 * pq->q2;
        f2q2q2 = f2q * pq->q2;
        f2q2q3 = f2q * pq->q3;
        f2q3q3 = 2.0 * pq->q3 * pq->q3;

        // calculate the rotation matrix assuming the quaternion is normalized
        R[X][X] = f2q0q0 + f2q1q1 - 1.0;
        R[X][Y] = f2q1q2 + f2q0q3;
        R[X][Z] = f2q1q3 - f2q0q2;
        R[Y][X] = f2q1q2 - f2q0q3;
        R[Y][Y] = f2q0q0 + f2q2q2 - 1.0;
        R[Y][Z] = f2q2q3 + f2q0q1;
        R[Z][X] = f2q1q3 + f2q0q2;
        R[Z][Y] = f2q2q3 - f2q0q1;
        R[Z][Z] = f2q0q0 + f2q3q3 - 1.0;

        return;
    }

    // function calculate the rotation vector from a rotation matrix
    template <typename Scalar>
    inline void fRotationVectorDegFromRotationMatrix(Scalar R[][3], Scalar rvecdeg[]) {
        Scalar ftrace;    // trace of the rotation matrix
        Scalar fetadeg;   // rotation angle eta (deg)
        Scalar fmodulus;  // modulus of axis * angle vector = 2|sin(eta)|
        Scalar ftmp;      // scratch variable

        // calculate the trace of the rotation matrix = 1+2cos(eta) in range -1 to +3
        // and eta (deg) in range 0 to 180 deg inclusive
        // checking for rounding errors that might take the trace outside this range
        ftrace = R[X][X] + R[Y][Y] + R[Z][Z];
        if (ftrace >= 3.0) {
            fetadeg = 0.0;
        }
        else if (ftrace <= -1.0) {
            fetadeg = 180.0;
        }
        else {
            fetadeg = std::acos(0.5 * (ftrace - 1.0)) * FRADTODEG;
        }

        // set the rvecdeg vector to differences across the diagonal = 2*n*sin(eta)
        // and calculate its modulus equal to 2|sin(eta)|
        // the modulus approaches zero near 0 and 180 deg (when sin(eta) approaches zero)
        rvecdeg[X] = R[Y][Z] - R[Z][Y];
        rvecdeg[Y] = R[Z][X] - R[X][Z];
        rvecdeg[Z] = R[X][Y] - R[Y][X];
        fmodulus   = std::sqrt(rvecdeg[X] * rvecdeg[X] + rvecdeg[Y] * rvecdeg[Y] + rvecdeg[Z] * rvecdeg[Z]);

        // normalize the rotation vector for general, 0 deg and 180 deg rotation cases
        if (fmodulus > SMALLMODULUS) {
            // general case away from 0 and 180 deg rotation
            ftmp = fetadeg / fmodulus;
            rvecdeg[X] *= ftmp;  // set x component to eta(deg) * nx
            rvecdeg[Y] *= ftmp;  // set y component to eta(deg) * ny
            rvecdeg[Z] *= ftmp;  // set z component to eta(deg) * nz
        }                        // end of general case
        else if (ftrace >= 0.0) {
            // near 0 deg rotation (trace = 3): matrix is nearly identity matrix
            // R[Y][Z]-R[Z][Y]=2*nx*eta(rad) and similarly for other components
            ftmp = 0.5 * FRADTODEG;
            rvecdeg[X] *= ftmp;
            rvecdeg[Y] *= ftmp;
            rvecdeg[Z] *= ftmp;
        }  // end of zero deg case
        else {
            // near 180 deg (trace = -1): matrix is nearly symmetric
            // calculate the absolute value of the components of the axis-angle vector
            rvecdeg[X] = 180.0 * std::sqrt(std::abs(0.5 * (R[X][X] + 1.0)));
            rvecdeg[Y] = 180.0 * std::sqrt(std::abs(0.5 * (R[Y][Y] + 1.0)));
            rvecdeg[Z] = 180.0 * std::sqrt(std::abs(0.5 * (R[Z][Z] + 1.0)));

            // correct the signs of the three components by examining the signs of differenced off-diagonal terms
            if ((R[Y][Z] - R[Z][Y]) < 0.0)
                rvecdeg[X] = -rvecdeg[X];
            if ((R[Z][X] - R[X][Z]) < 0.0)
                rvecdeg[Y] = -rvecdeg[Y];
            if ((R[X][Y] - R[Y][X]) < 0.0)
                rvecdeg[Z] = -rvecdeg[Z];

        }  // end of 180 deg case

        return;
    }

    // computes rotation vector (deg) from rotation quaternion
    template <typename Scalar>
    inline void fRotationVectorDegFromQuaternion(struct fquaternion* pq, Scalar rvecdeg[]) {
        Scalar fetarad;     // rotation angle (rad)
        Scalar fetadeg;     // rotation angle (deg)
        Scalar sinhalfeta;  // sin(eta/2)
        Scalar ftmp;        // scratch variable

        // calculate the rotation angle in the range 0 <= eta < 360 deg
        if ((pq->q0 >= 1.0) || (pq->q0 <= -1.0)) {
            // rotation angle is 0 deg or 2*180 deg = 360 deg = 0 deg
            fetarad = 0.0;
            fetadeg = 0.0;
        }
        else {
            // general case returning 0 < eta < 360 deg
            fetarad = 2.0 * acosf(pq->q0);
            fetadeg = fetarad * FRADTODEG;
        }

        // map the rotation angle onto the range -180 deg <= eta < 180 deg
        if (fetadeg >= 180.0) {
            fetadeg -= 360.0;
            fetarad = fetadeg * FDEGTORAD;
        }

        // calculate sin(eta/2) which will be in the range -1 to +1
        sinhalfeta = (Scalar) sinf(0.5 * fetarad);

        // calculate the rotation vector (deg)
        if (sinhalfeta == 0.0) {
            // the rotation angle eta is zero and the axis is irrelevant
            rvecdeg[X] = rvecdeg[Y] = rvecdeg[Z] = 0.0;
        }
        else {
            // general case with non-zero rotation angle
            ftmp       = fetadeg / sinhalfeta;
            rvecdeg[X] = pq->q1 * ftmp;
            rvecdeg[Y] = pq->q2 * ftmp;
            rvecdeg[Z] = pq->q3 * ftmp;
        }

        return;
    }

    // function low pass filters an orientation quaternion and computes virtual gyro rotation rate
    template <typename Scalar>
    inline void fLPFOrientationQuaternion(struct fquaternion* pq,
                                          struct fquaternion* pLPq,
                                          Scalar flpf,
                                          Scalar fdeltat,
                                          int loopcounter) {
        // local variables
        struct fquaternion fdeltaq;  // delta rotation quaternion
        Scalar rvecdeg[3];           // rotation vector (deg)
        Scalar ftmp;                 // scratch variable

        // initialize delay line on first pass: LPq[n]=q[n]
        if (loopcounter == 0) {
            *pLPq = *pq;
        }

        // set fdeltaqn to the delta rotation quaternion conjg(fLPq[n-1) . fqn
        fdeltaq = filter::utilities::qconjgAxB(pLPq, pq);
        if (fdeltaq.q0 < 0.0) {
            fdeltaq.q0 = -fdeltaq.q0;
            fdeltaq.q1 = -fdeltaq.q1;
            fdeltaq.q2 = -fdeltaq.q2;
            fdeltaq.q3 = -fdeltaq.q3;
        }

        // set ftmp to a scaled lpf value which equals flpf in the limit of small rotations (q0=1)
        // but which rises as the delta rotation angle increases (q0 tends to zero)
        ftmp = flpf + 0.75 * (1.0 - fdeltaq.q0);
        if (ftmp > 1.0) {
            ftmp = 1.0;
        }

        // scale the delta rotation by the corrected lpf value
        fdeltaq.q1 *= ftmp;
        fdeltaq.q2 *= ftmp;
        fdeltaq.q3 *= ftmp;

        // compute the scalar component q0
        ftmp = fdeltaq.q1 * fdeltaq.q1 + fdeltaq.q2 * fdeltaq.q2 + fdeltaq.q3 * fdeltaq.q3;
        if (ftmp <= 1.0) {
            // normal case
            fdeltaq.q0 = std::sqrt(1.0 - ftmp);
        }
        else {
            // rounding errors present so simply set scalar component to 0
            fdeltaq.q0 = 0.0;
        }

        // calculate the delta rotation vector from fdeltaqn and the virtual gyro angular velocity (deg/s)
        fRotationVectorDegFromQuaternion(&fdeltaq, rvecdeg);
        ftmp = 1.0 / fdeltat;

        // set LPq[n] = LPq[n-1] . deltaq[n]
        qAeqAxB(pLPq, &fdeltaq);

        // renormalize the low pass filtered quaternion to prevent error accumulation
        // the renormalization function ensures that q0 is non-negative
        filter::utilities::fqAeqNormqA<Scalar>(pLPq);

        return;
    }

    // function low pass filters a scalar
    template <typename Scalar>
    inline void fLPFScalar(Scalar* pfS, Scalar* pfLPS, Scalar flpf, int loopcounter) {
        // set S[LP,n]=S[n] on first pass
        if (loopcounter == 0) {
            *pfLPS = *pfS;
        }

        // apply the exponential low pass filter
        *pfLPS += flpf * (*pfS - *pfLPS);

        return;
    }
}  // namespace filter::utilities

#endif  // UTILITIES_HPP
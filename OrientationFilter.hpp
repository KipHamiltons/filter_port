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
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL FREESCALE SEMICONDUCTOR, INC. BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef ORIENTATION_FILTER_HPP
#define ORIENTATION_FILTER_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "build.hpp"
#include "utilities.hpp"
namespace filter::kalman {

    using filter::utilities::f3DOFTiltAndroid;
    using filter::utilities::f3DOFTiltNED;
    using filter::utilities::f3DOFTiltWin8;
    using filter::utilities::quat_from_rot_mat;
    using filter::utilities::quat_from_rot_vec;
    using filter::utilities::rot_mat_from_quat;
    using filter::utilities::set_identity3x3;
    using filter::utilities::set_inverse_inplace;

    // *********************************************************************************
    // COMPUTE_6DOF_GY_KALMAN constants
    // *********************************************************************************
    // kalman filter noise variances
    static constexpr double FQVA_6DOF_GY_KALMAN = 2e-6;  // accelerometer noise g^2 so 1.4mg RMS
    // gyro noise (deg/s)^2
    static constexpr double FQVG_6DOF_GY_KALMAN = 0.3;
    // gyro offset drift (deg/s)^2: 1E-9 implies 0.09deg/s max at 50Hz
    static constexpr double FQWB_6DOF_GY_KALMAN = 1e-9;
    // linear acceleration drift g^2 (increase slows convergence to g but reduces sensitivity to shake)
    static constexpr double FQWA_6DOF_GY_KALMAN = 1e-4;
    // initialization of Qw covariance matrix
    static constexpr double FQWINITTHTH_6DOF_GY_KALMAN = 2000e-5;  // th_e * th_e terms
    static constexpr double FQWINITBB_6DOF_GY_KALMAN   = 250e-3;   // for FXAS21000: b_e * b_e terms
    static constexpr double FQWINITTHB_6DOF_GY_KALMAN  = 0.0;      // th_e * b_e terms
    // a_e * a_e terms (increase slows convergence to g but reduces sensitivity to shake)
    static constexpr double FQWINITAA_6DOF_GY_KALMAN = 10e-5;
    // linear acceleration time constant
    static constexpr double FCA_6DOF_GY_KALMAN = 0.5;  // linear acceleration decay factor
    // 6DOF Kalman filter accelerometer and gyroscope state vector structure

    template <typename Scalar>
    class OrientationFilter {
    public:
        // orientation matrix, quaternion and rotation vector
        Scalar posterior_orientation_mat[3][3];                // a posteriori  rotation matrix
        Eigen::Quaternion<Scalar> posterior_orientation_quat;  // a posteriori orientation quaternion
        Scalar fbPl[3];                                        // gyro offset (deg/s)
        Scalar fThErrPl[3];                                    // orientation error (deg)
        Scalar fbErrPl[3];                                     // gyro offset error (deg/s)

        Scalar fRMi[3][3];               // a priori rotation matrix
        Eigen::Quaternion<Scalar> fqMi;  // a priori orientation quaternion
        Eigen::Quaternion<Scalar> fDeltaq =
            Eigen::Quaternion<Scalar>::Identity();  // delta a priori or a posteriori quaternion
        Scalar faSePl[3];                           // linear acceleration (g, sensor frame)
        Scalar faErrSePl[3];                        // linear acceleration error (g, sensor frame)
        Scalar fgErrSeMi[3];                        // difference (g, sensor frame) of gravity vector (accel)
                                                    // and gravity vector (gyro)
        Scalar fgSeGyMi[3];                         // gravity vector (g, sensor frame) measurement from gyro
        Scalar faSeMi[3];                           // linear acceleration (g, sensor frame)
        Scalar fQvAA;                               // accelerometer terms of Qv
        Scalar fPPlus9x9[9][9];                     // covariance matrix P+
        Scalar fK9x3[9][3];                         // kalman filter gain matrix K
        Scalar fQw9x9[9][9];                        // covariance matrix Qw
        Scalar fC3x9[3][9];                         // measurement matrix C
        Scalar fcasq;                               // FCA * FCA;
        Scalar fFastdeltat;                         // sensor sampling interval (s) = 1 / SENSORFS
        Scalar fdeltat;                             // kalman filter sampling interval (s) = OVERSAMPLE_RATIO /
                                                    // SENSORFS
        Scalar fdeltatsq;                           // fdeltat * fdeltat;
        Scalar fQwbplusQvG;                         // FQWB + FQVG;
        int iFirstOrientationLock;                  // denotes that 6DOF orientation has locked to 3DOF
        int resetflag;                              // flag to request re-initialization on next pass

        void init_filter(int iSensorFS, int iOverSampleRatio) {
            int i, j;  // loop counters

            // reset the flag denoting that a first 6DOpthisSVtion lock has been achieved
            iFirstOrientationLock = 0;

            // compute and store useful product terms to save Scalaring point calculations later
            fFastdeltat = 1.0 / (Scalar) iSensorFS;
            fdeltat     = (Scalar) iOverSampleRatio * fFastdeltat;
            fdeltatsq   = fdeltat * fdeltat;
            fcasq       = FCA_6DOF_GY_KALMAN * FCA_6DOF_GY_KALMAN;
            fQwbplusQvG = FQWB_6DOF_GY_KALMAN + FQVG_6DOF_GY_KALMAN;

            // initialize the fixed entries in the measurement matrix C
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 9; j++) {
                    fC3x9[i][j] = 0.0;
                }
            }
            fC3x9[0][6] = fC3x9[1][7] = fC3x9[2][8] = 1.0;

            // zero a posteriori orientation, error vector xe+ (thetae+, be+, ae+) and b+
            set_identity3x3(posterior_orientation_mat);
            posterior_orientation_quat = Eigen::Quaternion<Scalar>::Identity();
            for (i = X; i <= Z; i++) {
                fThErrPl[i] = fbErrPl[i] = faErrSePl[i] = fbPl[i] = 0.0;
            }

            // initialize noise variance for Qv and Qw matrix updates
            fQvAA = FQVA_6DOF_GY_KALMAN + FQWA_6DOF_GY_KALMAN
                    + FDEGTORAD * FDEGTORAD * fdeltatsq * (FQWB_6DOF_GY_KALMAN + FQVG_6DOF_GY_KALMAN);

            // initialize the 6x6 noise covariance matrix Qw of the a priori error vector xe-
            // Qw is then recursively updated as P+ = (1 - K * C) * P- = (1 - K * C) * Qw  and Qw updated from P+
            // zero the matrix Qw
            for (i = 0; i < 9; i++) {
                for (j = 0; j < 9; j++) {
                    fQw9x9[i][j] = 0.0;
                }
            }
            // loop over non-zero values
            for (i = 0; i < 3; i++) {
                // theta_e * theta_e terms
                fQw9x9[i][i] = FQWINITTHTH_6DOF_GY_KALMAN;
                // b_e * b_e terms
                fQw9x9[i + 3][i + 3] = FQWINITBB_6DOF_GY_KALMAN;
                // th_e * b_e terms
                fQw9x9[i][i + 3] = fQw9x9[i + 3][i] = FQWINITTHB_6DOF_GY_KALMAN;
                // a_e * a_e terms
                fQw9x9[i + 6][i + 6] = FQWINITAA_6DOF_GY_KALMAN;
            }

            // clear the reset flag
            resetflag = false;

            return;
        }

        void run_filter(Scalar accel_reading[3], Scalar gyro_reading[3], int ithisCoordSystem, int iOverSampleRatio) {

            // local arrays and scalars
            Scalar rvec[3];         // rotation vector
            Scalar ftmpA9x3[9][3];  // scratch array

            // assorted array pointers
            Scalar* pfPPlus9x9kj;
            Scalar* pfPPlus9x9ij;
            Scalar* pfK9x3ij;
            Scalar* pfK9x3ik;
            Scalar* pftmpA9x3ik;
            Scalar* pftmpA9x3ij;
            Scalar* pftmpA9x3kj;
            Scalar* pfQw9x9ij;
            Scalar* pfQw9x9ik;
            Scalar* pfQw9x9kj;
            Scalar* pfC3x9ik;
            Scalar* pfC3x9jk;

            int i, j, k;  // loop counters

            // working arrays for 3x3 matrix inversion
            Scalar* pfRows[3];
            int iColInd[3];
            int iRowInd[3];
            int iPivot[3];

            // do a reset and return if requested
            if (resetflag) {
                init_filter(SENSORFS, OVERSAMPLE_RATIO);
                return;
            }

            // do a once-only orientation lock to accelerometer tilt
            if (!iFirstOrientationLock) {
                // get the 3DOF orientation matrix and initial inclination angle
                if (ithisCoordSystem == NED) {
                    // call NED tilt function
                    // f3DOFTiltNED(posterior_orientation_mat, pthisAccel->fGpFast);
                    f3DOFTiltNED(posterior_orientation_mat, accel_reading);
                }
                else if (ithisCoordSystem == ANDROID) {
                    // call Android tilt function
                    // f3DOFTiltAndroid(posterior_orientation_mat, pthisAccel->fGpFast);
                    f3DOFTiltAndroid(posterior_orientation_mat, accel_reading);
                }
                else {
                    // call Windows 8 tilt function
                    // f3DOFTiltWin8(posterior_orientation_mat, pthisAccel->fGpFast);
                    f3DOFTiltWin8(posterior_orientation_mat, accel_reading);
                }

                // get the orientation quaternion from the orientation matrix
                quat_from_rot_mat(posterior_orientation_mat, posterior_orientation_quat);

                // set the orientation lock flag so this initial alignment is only performed once
                iFirstOrientationLock = 1;
            }

            // *********************************************************************************
            // calculate a priori rotation matrix
            // *********************************************************************************

            // initialize the a priori orientation quaternion to the a posteriori orientation estimate
            fqMi = posterior_orientation_quat;

            // integrate the buffered high frequency(typically 200Hz) gyro readings for (j = 0; j < iOverSampleRatio;
            // j++) { compute the incremental fast (typically 200Hz) rotation vector rvec (deg)
            for (i = X; i <= Z; i++) {
                // rvec[i] = (((Scalar) pthisGyro->iYpFast[j][i] * pthisGyro->fDegPerSecPerCount) - fbPl[i]) *
                // fFastdeltat;
            }

            // compute the incremental quaternion fDeltaq from the rotation vector
            quat_from_rot_vec(fDeltaq, rvec, 1.0);

            // incrementally rotate the a priori orientation quaternion fqMi
            // the a posteriori orientation is re-normalized later so this update is stable
            fqMi = fqMi * fDeltaq;

            // get the a priori rotation matrix from the a priori quaternion
            rot_mat_from_quat(fRMi, fqMi);

            // *********************************************************************************
            // calculate a priori gyro and accelerometer estimates of the gravity vector
            // and the error between the two
            // *********************************************************************************

            // compute the a priori **gyro** estimate of the gravitational vector (g, sensor frame)
            // using an absolute rotation of the global frame gravity vector (with magnitude 1g)
            for (i = X; i <= Z; i++) {
                if (ithisCoordSystem == NED) {
                    // NED gravity is along positive z axis
                    fgSeGyMi[i] = fRMi[i][Z];
                }
                else {
                    // Android and Win8 gravity are along negative z axis
                    fgSeGyMi[i] = -fRMi[i][Z];
                }

                // compute a priori acceleration (a-) (g, sensor frame) from a posteriori estimate (g, sensor frame)
                faSeMi[i] = FCA_6DOF_GY_KALMAN * faSePl[i];

                // compute the a priori gravity error vector (accelerometer minus gyro estimates) (g, sensor frame)
                if ((ithisCoordSystem == NED) || (ithisCoordSystem == WIN8)) {
                    // NED and Windows 8 have positive sign for gravity: y = g - a and g = y + a
                    // fgErrSeMi[i] = pthisAccel->fGpFast[i] + faSeMi[i] - fgSeGyMi[i];
                    fgErrSeMi[i] = accel_reading[i] + faSeMi[i] - fgSeGyMi[i];
                }
                else {
                    // Android has negative sign for gravity: y = a - g, g = -y + a
                    // fgErrSeMi[i] = -pthisAccel->fGpFast[i] + faSeMi[i] - fgSeGyMi[i];
                    fgErrSeMi[i] = -accel_reading[i] + faSeMi[i] - fgSeGyMi[i];
                }
            }

            // *********************************************************************************
            // update variable elements of measurement matrix C
            // *********************************************************************************

            // update measurement matrix C (3x9) with -alpha(g-)x from gyro (g, sensor frame)
            fC3x9[0][1] = FDEGTORAD * fgSeGyMi[Z];
            fC3x9[0][2] = -FDEGTORAD * fgSeGyMi[Y];
            fC3x9[1][2] = FDEGTORAD * fgSeGyMi[X];
            fC3x9[1][0] = -fC3x9[0][1];
            fC3x9[2][0] = -fC3x9[0][2];
            fC3x9[2][1] = -fC3x9[1][2];
            fC3x9[0][4] = -fdeltat * fC3x9[0][1];
            fC3x9[0][5] = -fdeltat * fC3x9[0][2];
            fC3x9[1][5] = -fdeltat * fC3x9[1][2];
            fC3x9[1][3] = -fC3x9[0][4];
            fC3x9[2][3] = -fC3x9[0][5];
            fC3x9[2][4] = -fC3x9[1][5];

            // *********************************************************************************
            // calculate the Kalman gain matrix K (9x3)
            // K = P- * C^T * inv(C * P- * C^T + Qv) = Qw * C^T * inv(C * Qw * C^T + Qv)
            // Qw is used as a proxy for P- throughout the code
            // P+ is used here as a working array to reduce RAM usage and is re-computed later
            // *********************************************************************************

            // set ftmpA9x3 = P- * C^T = Qw * C^T where Qw and C are both sparse
            // C also has a significant number of +1 and -1 entries
            // ftmpA9x3 is also sparse but not symmetric
            for (i = 0; i < 9; i++)  // loop over rows of ftmpA9x3
            {
                // initialize pftmpA9x3ij for current i, j=0
                pftmpA9x3ij = ftmpA9x3[i];

                for (j = 0; j < 3; j++)  // loop over columns of ftmpA9x3
                {
                    // zero ftmpA9x3[i][j]
                    *pftmpA9x3ij = 0.0;

                    // initialize pfC3x9jk for current j, k=0
                    pfC3x9jk = fC3x9[j];

                    // initialize pfQw9x9ik for current i, k=0
                    pfQw9x9ik = fQw9x9[i];

                    // sum matrix products over inner loop over k
                    for (k = 0; k < 9; k++) {
                        if ((*pfQw9x9ik != 0.0) && (*pfC3x9jk != 0.0)) {
                            if (*pfC3x9jk == 1.0)
                                *pftmpA9x3ij += *pfQw9x9ik;
                            else if (*pfC3x9jk == -1.0)
                                *pftmpA9x3ij -= *pfQw9x9ik;
                            else
                                *pftmpA9x3ij += *pfQw9x9ik * *pfC3x9jk;
                        }

                        // increment pfC3x9jk and pfQw9x9ik for next iteration of k
                        pfC3x9jk++;
                        pfQw9x9ik++;

                    }  // end of loop over k

                    // increment pftmpA9x3ij for next iteration of j
                    pftmpA9x3ij++;

                }  // end of loop over j
            }      // end of loop over i

            // set symmetric P+ (3x3 scratch sub-matrix) to C * P- * C^T + Qv
            // = C * (Qw * C^T) + Qv = C * ftmpA9x3 + Qv
            // both C and ftmpA9x3 are sparse but not symmetric
            for (i = 0; i < 3; i++)  // loop over rows of P+
            {
                // initialize pfPPlus9x9ij for current i, j=i
                pfPPlus9x9ij = fPPlus9x9[i] + i;

                for (j = i; j < 3; j++)  // loop over above diagonal columns of P+
                {
                    // zero P+[i][j]
                    *pfPPlus9x9ij = 0.0;

                    // initialize pfC3x9ik for current i, k=0
                    pfC3x9ik = fC3x9[i];

                    // initialize pftmpA9x3kj for current j, k=0
                    pftmpA9x3kj = *ftmpA9x3 + j;

                    // sum matrix products over inner loop over k
                    for (k = 0; k < 9; k++) {
                        if ((*pfC3x9ik != 0.0) && (*pftmpA9x3kj != 0.0)) {
                            if (*pfC3x9ik == 1.0)
                                *pfPPlus9x9ij += *pftmpA9x3kj;
                            else if (*pfC3x9ik == -1.0)
                                *pfPPlus9x9ij -= *pftmpA9x3kj;
                            else
                                *pfPPlus9x9ij += *pfC3x9ik * *pftmpA9x3kj;
                        }

                        // update pfC3x9ik and pftmpA9x3kj for next iteration of k
                        pfC3x9ik++;
                        pftmpA9x3kj += 3;

                    }  // end of loop over k

                    // increment pfPPlus9x9ij for next iteration of j
                    pfPPlus9x9ij++;

                }  // end of loop over j
            }      // end of loop over i

            // add in noise covariance terms to the diagonal
            fPPlus9x9[0][0] += fQvAA;
            fPPlus9x9[1][1] += fQvAA;
            fPPlus9x9[2][2] += fQvAA;

            // copy above diagonal terms of P+ (3x3 scratch sub-matrix) to below diagonal terms
            fPPlus9x9[1][0] = fPPlus9x9[0][1];
            fPPlus9x9[2][0] = fPPlus9x9[0][2];
            fPPlus9x9[2][1] = fPPlus9x9[1][2];

            // calculate inverse of P+ (3x3 scratch sub-matrix) = inv(C * P- * C^T + Qv) = inv(C * Qw * C^T + Qv)
            for (i = 0; i < 3; i++) {
                pfRows[i] = fPPlus9x9[i];
            }
            set_inverse_inplace(pfRows, iColInd, iRowInd, iPivot, 3);

            // set K = P- * C^T * inv(C * P- * C^T + Qv) = Qw * C^T * inv(C * Qw * C^T + Qv)
            // = ftmpA9x3 * P+ (3x3 sub-matrix)
            // ftmpA9x3 = Qw * C^T is sparse but P+ (3x3 sub-matrix) is not
            // K is not symmetric because C is not symmetric
            for (i = 0; i < 9; i++)  // loop over rows of K9x3
            {
                // initialize pfK9x3ij for i, j=0
                pfK9x3ij = fK9x3[i];

                for (j = 0; j < 3; j++)  // loop over columns of K9x3
                {
                    // zero the matrix element fK9x3[i][j]
                    *pfK9x3ij = 0.0;

                    // initialize pftmpA9x3ik for i, k=0
                    pftmpA9x3ik = ftmpA9x3[i];

                    // initialize pfPPlus9x9kj for j, k=0
                    pfPPlus9x9kj = *fPPlus9x9 + j;

                    // sum matrix products over inner loop over k
                    for (k = 0; k < 3; k++) {
                        if (*pftmpA9x3ik != 0.0) {
                            *pfK9x3ij += *pftmpA9x3ik * *pfPPlus9x9kj;
                        }

                        // increment pftmpA9x3ik and pfPPlus9x9kj for next iteration of k
                        pftmpA9x3ik++;
                        pfPPlus9x9kj += 9;

                    }  // end of loop over k

                    // increment pfK9x3ij for the next iteration of j
                    pfK9x3ij++;

                }  // end of loop over j
            }      // end of loop over i

            // *********************************************************************************
            // calculate a posteriori error estimate: xe+ = K * ze-
            // *********************************************************************************

            // update the a posteriori state vector
            for (i = X; i <= Z; i++) {
                // zero a posteriori error terms
                fThErrPl[i] = fbErrPl[i] = faErrSePl[i] = 0.0;

                // accumulate the error vector terms K * ze-
                for (k = 0; k < 3; k++) {
                    fThErrPl[i] += fK9x3[i][k] * fgErrSeMi[k];
                    fbErrPl[i] += fK9x3[i + 3][k] * fgErrSeMi[k];
                    faErrSePl[i] += fK9x3[i + 6][k] * fgErrSeMi[k];
                }
            }

            // *********************************************************************************
            // apply the a posteriori error corrections to the a posteriori state vector
            // *********************************************************************************

            // get the a posteriori delta quaternion
            quat_from_rot_vec(fDeltaq, fThErrPl, -1.0);

            // compute the a posteriori orientation quaternion posterior_orientation_quat = fqMi * Deltaq(-thetae+)
            // the resulting quaternion may have negative scalar component q0
            posterior_orientation_quat = fqMi * fDeltaq;

            // normalize the a posteriori orientation quaternion to stop error propagation
            // the renormalization function ensures that the scalar component q0 is non-negative
            posterior_orientation_quat.normalize();

            // compute the a posteriori rotation matrix from the a posteriori quaternion
            rot_mat_from_quat(posterior_orientation_mat, posterior_orientation_quat);

            // update the a posteriori gyro offset vector b+ and linear acceleration vector a+ (sensor frame)
            for (i = X; i <= Z; i++) {
                // b+[k] = b-[k] - be+[k] = b+[k] - be+[k] (deg/s)
                fbPl[i] -= fbErrPl[i];
                // a+ = a- - ae+ (g, sensor frame)
                faSePl[i] = faSeMi[i] - faErrSePl[i];
            }

            // ***********************************************************************************
            // calculate (symmetric) a posteriori error covariance matrix P+
            // P+ = (I12 - K * C) * P- = (I12 - K * C) * Qw = Qw - K * (C * Qw)
            // both Qw and P+ are used as working arrays in this section
            // at the end of this section, P+ is valid but Qw is over-written
            // ***********************************************************************************

            // set P+ (3x9 scratch sub-matrix) to the product C (3x9) * Qw (9x9)
            // where both C and Qw are sparse and C has a significant number of +1 entries
            // the resulting matrix is sparse but not symmetric
            for (i = 0; i < 3; i++)  // loop over the rows of P+
            {
                // initialize pfPPlus9x9ij for current i, j=0
                pfPPlus9x9ij = fPPlus9x9[i];

                for (j = 0; j < 9; j++)  // loop over the columns of P+
                {
                    // zero P+[i][j]
                    *pfPPlus9x9ij = 0.0;

                    // initialize pfC3x9ik for current i, k=0
                    pfC3x9ik = fC3x9[i];

                    // initialize pfQw9x9kj for current j, k=0
                    pfQw9x9kj = &fQw9x9[0][j];

                    // sum matrix products over inner loop over k
                    for (k = 0; k < 9; k++) {
                        if ((*pfC3x9ik != 0.0) && (*pfQw9x9kj != 0.0)) {
                            if (*pfC3x9ik == 1.0)
                                *pfPPlus9x9ij += *pfQw9x9kj;
                            else if (*pfC3x9ik == -1.0)
                                *pfPPlus9x9ij -= *pfQw9x9kj;
                            else
                                *pfPPlus9x9ij += *pfC3x9ik * *pfQw9x9kj;
                        }

                        // update pfC3x9ik and pfQw9x9kj for next iteration of k
                        pfC3x9ik++;
                        pfQw9x9kj += 9;

                    }  // end of loop over k

                    // increment pfPPlus9x9ij for next iteration of j
                    pfPPlus9x9ij++;

                }  // end of loop over j
            }      // end of loop over i

            // compute P+ = (I9 - K * C) * Qw = Qw - K * (C * Qw) = Qw - K * P+ (3x9 sub-matrix)
            // storing result P+ in Qw and over-writing Qw which is OK since Qw is later computed from P+
            // where working array P+ (6x12 sub-matrix) is sparse but K is not sparse
            // only on and above diagonal terms of P+ are computed since P+ is symmetric
            for (i = 0; i < 9; i++) {
                // initialize pfQw9x9ij for i, j=i
                pfQw9x9ij = fQw9x9[i] + i;

                for (j = i; j < 9; j++) {
                    // initialize pfK9x3ik for i, k=0
                    pfK9x3ik = fK9x3[i];

                    // initialize pfPPlus9x9kj for j, k=0
                    pfPPlus9x9kj = *fPPlus9x9 + j;

                    // compute on and above diagonal matrix entry
                    for (k = 0; k < 3; k++) {
                        // check for non-zero values since P+ (3x9 scratch sub-matrix) is sparse
                        if (*pfPPlus9x9kj != 0.0) {
                            *pfQw9x9ij -= *pfK9x3ik * *pfPPlus9x9kj;
                        }

                        // increment pfK9x3ik and pfPPlus9x9kj for next iteration of k
                        pfK9x3ik++;
                        pfPPlus9x9kj += 9;

                    }  // end of loop over k

                    // increment pfQw9x9ij for next iteration of j
                    pfQw9x9ij++;

                }  // end of loop over j
            }      // end of loop over i

            // Qw now holds the on and above diagonal elements of P+ (9x9)
            // so perform a simple copy to the all elements of P+
            // after execution of this code P+ is valid but Qw remains invalid
            for (i = 0; i < 9; i++) {
                // initialize pfPPlus9x9ij and pfQw9x9ij for i, j=i
                pfPPlus9x9ij = fPPlus9x9[i] + i;
                pfQw9x9ij    = fQw9x9[i] + i;

                // copy the on-diagonal elements and increment pointers to enter loop at j=i+1
                *(pfPPlus9x9ij++) = *(pfQw9x9ij++);

                // loop over above diagonal columns j copying to below-diagonal elements
                for (j = i + 1; j < 9; j++) {
                    *(pfPPlus9x9ij++) = fPPlus9x9[j][i] = *(pfQw9x9ij++);
                }
            }

            // *********************************************************************************
            // re-create the noise covariance matrix Qw=fn(P+) for the next iteration
            // using the elements of P+ which are now valid
            // Qw was over-written earlier but is here recomputed (all elements)
            // *********************************************************************************

            // zero the matrix Qw (9x9)
            for (i = 0; i < 9; i++) {
                for (j = 0; j < 9; j++) {
                    fQw9x9[i][j] = 0.0;
                }
            }

            // update the covariance matrix components
            for (i = 0; i < 3; i++) {
                // Qw[th-th-] = Qw[0-2][0-2] = E[th-(th-)^T] = Q[th+th+] + deltat^2 * (Q[b+b+] + (Qwb + QvG) * I)
                fQw9x9[i][i] = fPPlus9x9[i][i] + fdeltatsq * (fPPlus9x9[i + 3][i + 3] + fQwbplusQvG);

                // Qw[b-b-] = Qw[3-5][3-5] = E[b-(b-)^T] = Q[b+b+] + Qwb * I
                fQw9x9[i + 3][i + 3] = fPPlus9x9[i + 3][i + 3] + FQWB_6DOF_GY_KALMAN;

                // Qw[th-b-] = Qw[0-2][3-5] = E[th-(b-)^T] = -deltat * (Q[b+b+] + Qwb * I) = -deltat * Qw[b-b-]
                fQw9x9[i][i + 3] = fQw9x9[i + 3][i] = -fdeltat * fQw9x9[i + 3][i + 3];

                // Qw[a-a-] = Qw[6-8][6-8] = E[a-(a-)^T] = ca^2 * Q[a+a+] + Qwa * I
                fQw9x9[i + 6][i + 6] = fcasq * fPPlus9x9[i + 6][i + 6] + FQWA_6DOF_GY_KALMAN;
            }

            return;
        }
    };

}  // namespace filter::kalman
#endif  // #ifndef ORIENTATION_FILTER_HPP
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

#include "build.hpp"
namespace filter::kalman {
    // *********************************************************************************
    // COMPUTE_6DOF_GY_KALMAN constants
    // *********************************************************************************
    // kalman filter noise variances
    static constexpr double FQVA_6DOF_GY_KALMAN = 2e-6;  // accelerometer noise g^2 so 1.4mg RMS
    // gyro noise (deg/s)^2
    static constexpr double FQVG_6DOF_GY_KALMAN = 0.3F;
    // gyro offset drift (deg/s)^2: 1E-9 implies 0.09deg/s max at 50Hz
    static constexpr double FQWB_6DOF_GY_KALMAN = 1e-9;
    // linear acceleration drift g^2 (increase slows convergence to g but reduces sensitivity to shake)
    static constexpr double FQWA_6DOF_GY_KALMAN = 1e-4;
    // initialization of Qw covariance matrix
    static constexpr double FQWINITTHTH_6DOF_GY_KALMAN = 2000e-5;  // th_e * th_e terms
    static constexpr double FQWINITBB_6DOF_GY_KALMAN   = 250e-3;   // for FXAS21000: b_e * b_e terms
    static constexpr double FQWINITTHB_6DOF_GY_KALMAN  = 0.0;      // th_e * b_e terms
    // a_e * a_e terms (increase slows convergence to g but reduces sensitivity to shake)
    static constexpr double FQWINITAA_6DOF_GY_KALMAN = 10e-5F;
    // linear acceleration time constant
    static constexpr double FCA_6DOF_GY_KALMAN = 0.5;  // linear acceleration decay factor
    // 6DOF Kalman filter accelerometer and gyroscope state vector structure

    class OrientationFilter {
    public:
        // start: elements common to all motion state vectors
        // Euler angles
        float fPhiPl;  // roll (deg)
        float fThePl;  // pitch (deg)
        float fPsiPl;  // yaw (deg)
        float fRhoPl;  // compass (deg)
        float fChiPl;  // tilt from vertical (deg)
        // orientation matrix, quaternion and rotation vector
        float fRPl[3][3];         // a posteriori  rotation matrix
        struct fquaternion fqPl;  // a posteriori orientation quaternion
        float fRVecPl[3];         // rotation vector
        // angular velocity
        float fOmega[3];  // angular velocity (deg/s)
        // systick timer for benchmarking
        int systick;  // systick timer
        // end: elements common to all motion state vectors

        // elements transmitted over bluetooth in kalman packet
        float fbPl[3];      // gyro offset (deg/s)
        float fThErrPl[3];  // orientation error (deg)
        float fbErrPl[3];   // gyro offset error (deg/s)
        // end elements transmitted in kalman packet

        float fzErrMi[3];            // angular error (deg) between a priori and eCompass
                                     // orientations
        float fRMi[3][3];            // a priori rotation matrix
        struct fquaternion fqMi;     // a priori orientation quaternion
        struct fquaternion fDeltaq;  // delta a priori or a posteriori quaternion
        float faSePl[3];             // linear acceleration (g, sensor frame)
        float faErrSePl[3];          // linear acceleration error (g, sensor frame)
        float fgErrSeMi[3];          // difference (g, sensor frame) of gravity vector (accel)
                                     // and gravity vector (gyro)
        float fgSeGyMi[3];           // gravity vector (g, sensor frame) measurement from gyro
        float faSeMi[3];             // linear acceleration (g, sensor frame)
        float fQvAA;                 // accelerometer terms of Qv
        float fPPlus9x9[9][9];       // covariance matrix P+
        float fK9x3[9][3];           // kalman filter gain matrix K
        float fQw9x9[9][9];          // covariance matrix Qw
        float fC3x9[3][9];           // measurement matrix C
        float fcasq;                 // FCA * FCA;
        float fFastdeltat;           // sensor sampling interval (s) = 1 / SENSORFS
        float fdeltat;               // kalman filter sampling interval (s) = OVERSAMPLE_RATIO /
                                     // SENSORFS
        float fdeltatsq;             // fdeltat * fdeltat;
        float fQwbplusQvG;           // FQWB + FQVG;
        int iFirstOrientationLock;   // denotes that 6DOF orientation has locked to 3DOF
        int resetflag;               // flag to request re-initialization on next pass

        void init_filter(int iSensorFS, int iOverSampleRatio);
        void run_filter(float accel_reading[3], float gyro_reading[3], int ithisCoordSystem, int iOverSampleRatio);
    };

}  // namespace filter::kalman
#endif  // #ifndef ORIENTATION_FILTER_HPP
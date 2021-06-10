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
#ifndef TASKS_HPP
#define TASKS_HPP

#include "build.hpp"
namespace filter::tasks {
    // 6DOF Kalman filter accelerometer and gyroscope state vector structure
    struct SV_6DOF_GY_KALMAN {
        // start: elements common to all motion state vectors
        // Euler angles
        double fPhiPl;  // roll (deg)
        double fThePl;  // pitch (deg)
        double fPsiPl;  // yaw (deg)
        double fRhoPl;  // compass (deg)
        double fChiPl;  // tilt from vertical (deg)
        // orientation matrix, quaternion and rotation vector
        double fRPl[3][3];        // a posteriori  rotation matrix
        struct fquaternion fqPl;  // a posteriori orientation quaternion
        double fRVecPl[3];        // rotation vector
        // angular velocity
        double fOmega[3];  // angular velocity (deg/s)
        // systick timer for benchmarking
        int32 systick;  // systick timer
        // end: elements common to all motion state vectors

        // elements transmitted over bluetooth in kalman packet
        double fbPl[3];      // gyro offset (deg/s)
        double fThErrPl[3];  // orientation error (deg)
        double fbErrPl[3];   // gyro offset error (deg/s)
        // end elements transmitted in kalman packet

        double fzErrMi[3];            // angular error (deg) between a priori and eCompass
                                      // orientations
        double fRMi[3][3];            // a priori rotation matrix
        struct fquaternion fqMi;      // a priori orientation quaternion
        struct fquaternion fDeltaq;   // delta a priori or a posteriori quaternion
        double faSePl[3];             // linear acceleration (g, sensor frame)
        double faErrSePl[3];          // linear acceleration error (g, sensor frame)
        double fgErrSeMi[3];          // difference (g, sensor frame) of gravity vector (accel)
                                      // and gravity vector (gyro)
        double fgSeGyMi[3];           // gravity vector (g, sensor frame) measurement from gyro
        double faSeMi[3];             // linear acceleration (g, sensor frame)
        double fQvAA;                 // accelerometer terms of Qv
        double fPPlus9x9[9][9];       // covariance matrix P+
        double fK9x3[9][3];           // kalman filter gain matrix K
        double fQw9x9[9][9];          // covariance matrix Qw
        double fC3x9[3][9];           // measurement matrix C
        double fcasq;                 // FCA * FCA;
        double fFastdeltat;           // sensor sampling interval (s) = 1 / SENSORFS
        double fdeltat;               // kalman filter sampling interval (s) = OVERSAMPLE_RATIO /
                                      // SENSORFS
        double fdeltatsq;             // fdeltat * fdeltat;
        double fQwbplusQvG;           // FQWB + FQVG;
        int16 iFirstOrientationLock;  // denotes that 6DOF orientation has locked to 3DOF
        int8 resetflag;               // flag to request re-initialization on next pass
    };

    // globals defined in tasks_func.c declared here for use elsewhere
    extern struct AccelSensor thisAccel;
    extern struct GyroSensor thisGyro;
    extern struct SV_6DOF_GY_KALMAN thisSV_6DOF_GY_KALMAN;

    // function prototypes for functions in tasks_func.c
    void ApplyAccelHAL(struct AccelSensor* pthisAccel);
    void ApplyMagHAL(struct MagSensor* pthisMag);
    void ApplyGyroHAL(struct GyroSensor* pthisGyro, int16 irow);
    void RdSensData_Init();
    void RdSensData_Run();
    void Fusion_Init();
    void Fusion_Run();
}  // namespace filter::tasks
#endif  // #ifndef TASKS_HPP
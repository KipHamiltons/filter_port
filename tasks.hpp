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

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "build.hpp"
namespace filter::tasks {
    // 6DOF Kalman filter accelerometer and gyroscope state vector structure
    // template <typename Scalar>
    struct SV_6DOF_GY_KALMAN {
        using Scalar = double;
        // start: elements common to all motion state vectors
        // Euler angles
        Scalar fPhiPl;  // roll (deg)
        Scalar fThePl;  // pitch (deg)
        Scalar fPsiPl;  // yaw (deg)
        Scalar fRhoPl;  // compass (deg)
        Scalar fChiPl;  // tilt from vertical (deg)
        // orientation matrix, quaternion and rotation vector
        Eigen::Matrix<Scalar, 3, 3> fRPl = Eigen::Matrix<Scalar, 3, 3>::Identity();  // a posteriori  rotation matrix
        Eigen::Quaternion<Scalar> fqPl = Eigen::Quaternion<Scalar>::Identity();  // a posteriori orientation quaternion
        Eigen::Matrix<Scalar, 3, 1> fRVecPl;                                     // rotation vector
        // angular velocity
        Eigen::Matrix<Scalar, 3, 1> fOmega;  // angular velocity (deg/s)
        // systick timer for benchmarking
        int systick;  // systick timer
        // end: elements common to all motion state vectors

        // elements transmitted over bluetooth in kalman packet
        Eigen::Matrix<Scalar, 3, 1> fbPl     = Eigen::Matrix<Scalar, 3, 1>::Zero();  // gyro offset (deg/s)
        Eigen::Matrix<Scalar, 3, 1> fThErrPl = Eigen::Matrix<Scalar, 3, 1>::Zero();  // orientation error (deg)
        Eigen::Matrix<Scalar, 3, 1> fbErrPl  = Eigen::Matrix<Scalar, 3, 1>::Zero();  // gyro offset error (deg/s)
        // end elements transmitted in kalman packet

        Eigen::Matrix<Scalar, 3, 1> fzErrMi;  // angular error (deg) between a priori and eCompass
                                              // orientations
        Eigen::Matrix<Scalar, 3, 3> fRMi;     // a priori rotation matrix
        Eigen::Quaternion<Scalar> fqMi;       // a priori orientation quaternion
        Eigen::Quaternion<Scalar> fDeltaq;    // delta a priori or a posteriori quaternion
        Eigen::Matrix<Scalar, 3, 1> faSePl;   // linear acceleration (g, sensor frame)
        Eigen::Matrix<Scalar, 3, 1> faErrSePl =
            Eigen::Matrix<Scalar, 3, 1>::Zero();  // linear acceleration error (g, sensor frame)
        Eigen::Matrix<Scalar, 3, 1> fgErrSeMi;    // difference (g, sensor frame) of gravity vector (accel)
                                                  // and gravity vector (gyro)
        Eigen::Matrix<Scalar, 3, 1> fgSeGyMi;     // gravity vector (g, sensor frame) measurement from gyro
        Eigen::Matrix<Scalar, 3, 1> faSeMi;       // linear acceleration (g, sensor frame)
        Scalar fQvAA;                             // accelerometer terms of Qv
        Eigen::Matrix<Scalar, 9, 9> fPPlus9x9;    // covariance matrix P+
        Eigen::Matrix<Scalar, 9, 3> fK9x3;        // kalman filter gain matrix K
        Eigen::Matrix<Scalar, 9, 9> fQw9x9 = Eigen::Matrix<Scalar, 9, 9>::Zero();  // covariance matrix Qw
        Eigen::Matrix<Scalar, 3, 9> fC3x9  = Eigen::Matrix<Scalar, 3, 9>::Zero();  // measurement matrix C
        Scalar fcasq;                                                              // FCA * FCA;
        Scalar fFastdeltat;         // sensor sampling interval (s) = 1 / SENSORFS
        Scalar fdeltat;             // kalman filter sampling interval (s) = OVERSAMPLE_RATIO /
                                    // SENSORFS
        Scalar fdeltatsq;           // fdeltat * fdeltat;
        Scalar fQwbplusQvG;         // FQWB + FQVG;
        int iFirstOrientationLock;  // denotes that 6DOF orientation has locked to 3DOF
        bool resetflag = false;     // flag to request re-initialization on next pass
    };

    // globals defined in tasks_func.c declared here for use elsewhere
    extern struct AccelSensor thisAccel;
    extern struct GyroSensor thisGyro;
    // extern struct SV_6DOF_GY_KALMAN<double> thisSV_6DOF_GY_KALMAN;
    extern struct SV_6DOF_GY_KALMAM thisSV_6DOF_GY_KALMAN;

    // function prototypes for functions in tasks_func.c
    // void ApplyAccelHAL(struct AccelSensor* pthisAccel);
    // void ApplyMagHAL(struct MagSensor* pthisMag);
    // void ApplyGyroHAL(struct GyroSensor* pthisGyro, int irow);
    // void RdSensData_Init();
    // void RdSensData_Run();
    // void Fusion_Init();
    // void Fusion_Run();
}  // namespace filter::tasks
#endif  // #ifndef TASKS_HPP
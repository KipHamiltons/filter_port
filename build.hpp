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
// This file contains build controls for a sensor fusion project.
// Board and MCU customization is done via Processor Expert.  The fusion
// code sits above that software layer, and can usually be ported from
// one board environment to another with no changes.

#ifndef BUILD_HPP
#define BUILD_HPP

// coordinate system for the build
#define NED             0        // identifier for NED angle output
#define ANDROID         1        // identifier for Android angle output
#define WIN8            2        // identifier for Windows 8 angle output
#define THISCOORDSYSTEM ANDROID  // the coordinate system to be used

#define COMPUTE_6DOF_GY_KALMAN  // 6DOF accel and gyro (Kalman): (1x accel + 1x gyro)

// sampling rate and kalman filter timing
#define FTM_INCLK_HZ     1000000  // int32: 1MHz FTM timer frequency set in PE: do not change
#define SENSORFS         200      // int32: 200Hz: frequency (Hz) of sensor sampling process
#define OVERSAMPLE_RATIO 1        // int32: 8x: 3DOF, 6DOF, 9DOF run at SENSORFS / OVERSAMPLE_RATIO Hz

// geomagnetic model parameters
#define DEFAULTB 50.0F  // default geomagnetic field (uT)

// useful multiplicative conversion constants
static constexpr double PI 3.141592654;               // Pi
static constexpr double FDEGTORAD 0.01745329251994;   // degrees to radians conversion = pi / 180
static constexpr double FRADTODEG 57.2957795130823;   // radians to degrees conversion = 180 / pi
static constexpr double FRECIP180 0.0055555555555;    // multiplicative factor 1/180
static constexpr double ONETHIRD 0.33333333;          // one third
static constexpr double ONESIXTH 0.166666667;         // one sixth
static constexpr double ONETWELFTH 0.0833333333;      // one twelfth
static constexpr double ONEOVER48 0.02083333333;      // 1 / 48
static constexpr double ONEOVER120 0.0083333333;      // 1 / 120
static constexpr double ONEOVER3840 0.0002604166667;  // 1 / 3840
static constexpr double ONEOVERROOT2 0.707106781;     // 1/sqrt(2)
static constexpr double ROOT3OVER2 0.866025403784;    // sqrt(3)/2

// the quaternion type to be transmitted
using quaternion_type = enum quaternion { Q3, Q3M, Q3G, Q6MA, Q6AG, Q9 };

// quaternion structure definition
// struct fquaternion {
//     double q0;  // scalar component
//     double q1;  // x vector component
//     double q2;  // y vector component
//     double q3;  // z vector component
// };

// We only care about fGpFast[3], so we'll just use that instead of this struct

// // accelerometer sensor structure definition
// struct AccelSensor {
//     int32 iSumGpFast[3];  // sum of fast measurements
//     double fGpFast[3];     // fast (typically 200Hz) readings (g)
//     double fGp[3];         // slow (typically 25Hz) averaged readings (g)
//     double fgPerCount;     // initialized to FGPERCOUNT
//     int16 iGpFast[3];     // fast (typically 200Hz) readings
//     int16 iGp[3];         // slow (typically 25Hz) averaged readings (counts)
// };

// For OVERSAMPLE_RATIO A.K.A. DECIMATION_FACTOR=1, we only care about fYp[3].
// For higher oversample ratio, we care about iYpFast.

// // gyro sensor structure definition
// struct GyroSensor {
//     int32 iSumYpFast[3];                 // sum of fast measurements
//     double fYp[3];                        // raw gyro sensor output (deg/s)
//     double fDegPerSecPerCount;            // initialized to FDEGPERSECPERCOUNT
//     int16 iYpFast[OVERSAMPLE_RATIO][3];  // fast (typically 200Hz) readings
//     int16 iYp[3];                        // averaged gyro sensor output (counts)
// };

#endif  // BUILD_HPP

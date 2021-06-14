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
static constexpr int NED     = 0;  // identifier for NED angle output
static constexpr int ANDROID = 1;  // identifier for Android angle output
static constexpr int WIN8    = 2;  // identifier for Windows 8 angle output

// sampling rate and kalman filter timing
static constexpr int FTM_INCLK_HZ     = 1000000;  // int: 1MHz FTM timer frequency set in PE: do not change
static constexpr int SENSORFS         = 200;      // int: 200Hz: frequency (Hz) of sensor sampling process
static constexpr int OVERSAMPLE_RATIO = 8;        // int: 8x: 3DOF, 6DOF, 9DOF run at SENSORFS / OVERSAMPLE_RATIO Hz

// vector components
static constexpr int X = 0;
static constexpr int Y = 1;
static constexpr int Z = 2;

// useful multiplicative conversion constants
static constexpr double PI           = 3.141592654;       // Pi
static constexpr double FDEGTORAD    = 0.01745329251994;  // degrees to radians conversion = pi / 180
static constexpr double FRADTODEG    = 57.2957795130823;  // radians to degrees conversion = 180 / pi
static constexpr double FRECIP180    = 0.0055555555555;   // multiplicative factor 1/180
static constexpr double ONETHIRD     = 0.33333333;        // one third
static constexpr double ONESIXTH     = 0.166666667;       // one sixth
static constexpr double ONETWELFTH   = 0.0833333333;      // one twelfth
static constexpr double ONEOVER48    = 0.02083333333;     // 1 / 48
static constexpr double ONEOVER120   = 0.0083333333;      // 1 / 120
static constexpr double ONEOVER3840  = 0.0002604166667;   // 1 / 3840
static constexpr double ONEOVERROOT2 = 0.707106781;       // 1/sqrt(2)
static constexpr double ROOT3OVER2   = 0.866025403784;    // sqrt(3)/2

// quaternion structure definition
struct fquaternion {
    double q0 = 1.0;  // scalar component
    double q1 = 0.0;  // x vector component
    double q2 = 0.0;  // y vector component
    double q3 = 0.0;  // z vector component
};

#endif  // BUILD_HPP

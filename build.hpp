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
static constexpr float PI           = 3.141592654F;       // Pi
static constexpr float FDEGTORAD    = 0.01745329251994F;  // degrees to radians conversion = pi / 180
static constexpr float FRADTODEG    = 57.2957795130823F;  // radians to degrees conversion = 180 / pi
static constexpr float FRECIP180    = 0.0055555555555F;   // multiplicative factor 1/180
static constexpr float ONETHIRD     = 0.33333333F;        // one third
static constexpr float ONESIXTH     = 0.166666667F;       // one sixth
static constexpr float ONETWELFTH   = 0.0833333333F;      // one twelfth
static constexpr float ONEOVER48    = 0.02083333333F;     // 1 / 48
static constexpr float ONEOVER120   = 0.0083333333F;      // 1 / 120
static constexpr float ONEOVER3840  = 0.0002604166667F;   // 1 / 3840
static constexpr float ONEOVERROOT2 = 0.707106781F;       // 1/sqrt(2)
static constexpr float ROOT3OVER2   = 0.866025403784F;    // sqrt(3)/2

// quaternion structure definition
struct fquaternion {
    float q0 = 1.0f;  // scalar component
    float q1 = 0.0f;  // x vector component
    float q2 = 0.0f;  // y vector component
    float q3 = 0.0f;  // z vector component
};

#endif  // BUILD_HPP

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
    struct SV_6DOF_GY_KALMAN {
        // start: elements common to all motion state vectors
        // Euler angles
        double roll_deg;     // roll (deg)
        double pitch_deg;    // pitch (deg)
        double yaw_deg;      // yaw (deg)
        double compass_deg;  // compass (deg)
        double tilt_deg;     // tilt from vertical (deg)
        // orientation matrix, quaternion and rotation vector
        Eigen::Matrix<double, 3, 3> posterior_rot_matrix;      // a posteriori  rotation matrix
        Eigen::Quaternion<double> posterior_orientation_quat;  // a posteriori orientation quaternion
        Eigen::Matrix<double, 3, 1> rot_vec;                   // rotation vector
        // angular velocity
        Eigen::Matrix<double, 3, 1> angular_velocity_vec;  // angular velocity (deg/s)
        // end: elements common to all motion state vectors

        // elements transmitted over bluetooth in kalman packet
        Eigen::Matrix<double, 3, 1> gyro_offset;            // gyro offset (deg/s)
        Eigen::Matrix<double, 3, 1> orientation_error_deg;  // orientation error (deg)
        Eigen::Matrix<double, 3, 1> gyro_offset_error;      // gyro offset error (deg/s)
        // end elements transmitted in kalman packet

        Eigen::Matrix<double, 3, 3> prior_rotation_matrix;             // a priori rotation matrix
        Eigen::Quaternion<double> prior_rotation_quat;                 // a priori orientation quaternion
        Eigen::Quaternion<double> delta_quaternion;                    // delta a priori or a posteriori quaternion
        Eigen::Matrix<double, 3, 1> linear_accel_g1;                   // linear acceleration (g, sensor frame)
        Eigen::Matrix<double, 3, 1> linear_accel_error_g1;             // linear acceleration error (g, sensor frame)
        Eigen::Matrix<double, 3, 1> gravity_accel_minus_gravity_gyro;  // difference (g, sensor frame) of gravity vector
                                                                       // (accel) and gravity vector (gyro)
        Eigen::Matrix<double, 3, 1> gyro_gravity_g;   // gravity vector (g, sensor frame) measurement from gyro
        Eigen::Matrix<double, 3, 1> linear_accel_g2;  // linear acceleration (g, sensor frame)
        double accel_term_from_Qv;                    // accelerometer terms of Qv
        Eigen::Matrix<double, 9, 9> P_plus;           // covariance matrix P+
        Eigen::Matrix<double, 9, 3> K;                // kalman filter gain matrix K
        Eigen::Matrix<double, 9, 9> Qw;               // covariance matrix Qw
        Eigen::Matrix<double, 3, 9> C;                // measurement matrix C
        static constexpr double FCA_squared     = FCA_6DOF_GY_KALMAN * FCA_6DOF_GY_KALMAN;  // FCA * FCA;
        static constexpr double DELTA_T         = OVERSAMPLE_RATIO / SENSORFS;  // kalman filter sampling interval (s)
        static constexpr double DELTA_T_SQUARED = DELTA_T * DELTA_T;
        static constexpr double FQWB_plus_FQVG  = FQWB_6DOF_GY_KALMAN + FQVG_6DOF_GY_KALMAN;
        bool iFirstOrientationLock;  // denotes that 6DOF orientation has locked to 3DOF
        bool resetflag;              // flag to request re-initialization on next pass
    };

    // globals defined in tasks_func.c declared here for use elsewhere
    extern struct AccelSensor thisAccel;
    extern struct GyroSensor thisGyro;
    extern struct SV_6DOF_GY_KALMAN thisSV_6DOF_GY_KALMAN;

}  // namespace filter::tasks
#endif  // #ifndef TASKS_HPP
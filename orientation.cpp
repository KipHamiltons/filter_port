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
#include "orientation.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

#include "build.hpp"
#include "utilities.hpp"


namespace filter::orientation {
    // compile time constants that are private to this file
    static constexpr double SMALLQ0      = 0.01;   // limit of quaternion scalar component requiring special algorithm
    static constexpr double CORRUPTQUAT  = 0.001;  // threshold for deciding rotation quaternion is corrupt
    static constexpr double SMALLMODULUS = 0.01;   // limit where rounding errors may appear

    using filter::utilities::acos_deg;
    using filter::utilities::asin_deg;
    using filter::utilities::atan2_deg;
    using filter::utilities::atan_deg;

    // Aerospace NED accelerometer 3DOF tilt function computing rotation matrix fR
    template <typename Scalar>
    Eigen::Matrix<Scalar, 3, 3> f3DOFTiltNED(Eigen::Matrix<Scalar, 3, 3> fR,
                                             const Eigen::Matrix<Scalar, 3, 1> measurement) {
        // the NED self-consistency twist occurs at 90 deg pitch

        // local variables
        int i                = 0;    // counter
        double fmodGxyz      = 0.0;  // modulus of the x, y, z accelerometer readings
        double fmodGyz       = 0.0;  // modulus of the y, z accelerometer readings
        double frecipmodGxyz = 0.0;  // reciprocal of modulus
        double ftmp          = 0.0;  // scratch variable

        // compute the accelerometer squared magnitudes
        fmodGyz  = measurement.y() * measurement.y() + measurement.z() * measurement.z();
        fmodGxyz = fmodGyz + measurement.x() * measurement.x();

        // check for freefall special case where no solution is possible
        if (fmodGxyz == 0.0) {
            // f3x3matrixAeqI(fR);
            fR = Eigen::Matrix<Scalar, 3, 3>::Identity();
            return;
        }

        // check for vertical up or down gimbal lock case
        if (fmodGyz == 0.0) {
            // f3x3matrixAeqScalar(fR, 0.0);
            fR       = Eigen::Matrix<Scalar, 3, 3>::Zero();
            fR[1][1] = 1.0;
            if (fGp[0] >= 0.0) {
                fR[0][2] = 1.0;
                fR[2][0] = -1.0;
            }
            else {
                fR[0][2] = -1.0;
                fR[2][0] = 1.0;
            }
            return;
        }

        // compute moduli for the general case
        fmodGyz       = std::sqrt(fmodGyz);
        fmodGxyz      = std::sqrt(fmodGxyz);
        frecipmodGxyz = 1.0 / fmodGxyz;
        ftmp          = fmodGxyz / fmodGyz;

        // normalize the accelerometer reading into the z column
        for (i = 0; i <= 2; i++) {
            fR[i][2] = fGp[i] * frecipmodGxyz;
        }

        // construct x column of orientation matrix
        fR[0][0] = fmodGyz * frecipmodGxyz;
        fR[1][0] = -fR[0][2] * fR[1][2] * ftmp;
        fR[2][0] = -fR[0][2] * fR[2][2] * ftmp;

        // // construct y column of orientation matrix
        fR[0][1] = 0.0;
        fR[1][1] = fR[2][2] * ftmp;
        fR[2][1] = -fR[1][2] * ftmp;
        return fR;
    }

    // Android accelerometer 3DOF tilt function computing rotation matrix fR
    template <typename Scalar>
    Eigen::Matrix<Scalar, 3, 3> f3DOFTiltAndroid(Eigen::Matrix<Scalar, 3, 3> fR,
                                                 const Eigen::Matrix<Scalar, 3, 1> measurement) {
        // the Android tilt matrix is mathematically identical to the NED tilt matrix
        // the Android self-consistency twist occurs at 90 deg roll
        return f3DOFTiltNED(fR, measurement);
    }

    // Windows 8 accelerometer 3DOF tilt function computing rotation matrix fR
    Eigen::Matrix<Scalar, 3, 3> f3DOFTiltWin8(Eigen::Matrix<Scalar, 3, 3> fR,
                                              const Eigen::Matrix<Scalar, 3, 1> measurement) {
        // the Win8 self-consistency twist occurs at 90 deg roll

        // local variables
        double fmodGxyz      = 0.0;  // modulus of the x, y, z accelerometer readings
        double fmodGxz       = 0.0;  // modulus of the x, z accelerometer readings
        double frecipmodGxyz = 0.0;  // reciprocal of modulus
        double ftmp          = 0.0;  // scratch variable
        int i                = 0;    // counter

        // compute the accelerometer squared magnitudes
        fmodGxz  = fGp[0] * fGp[0] + fGp[2] * fGp[2];
        fmodGxyz = fmodGxz + fGp[1] * fGp[1];

        // check for freefall special case where no solution is possible
        if (fmodGxyz == 0.0) {
            // f3x3matrixAeqI(fR);
            return Eigen::Matrix<Scalar, 3, 3>::Identity();
        }

        // check for vertical up or down gimbal lock case
        if (fmodGxz == 0.0) {
            // f3x3matrixAeqScalar(fR, 0.0);
            fR       = Eigen::Matrix<Scalar, 3, 3>::Zero();
            fR[0][0] = 1.0;
            if (fGp[1] >= 0.0) {
                fR[1][2] = -1.0;
                fR[2][1] = 1.0;
            }
            else {
                fR[1][2] = 1.0;
                fR[2][1] = -1.0;
            }
            return fR;
        }

        // compute moduli for the general case
        fmodGxz       = std::sqrt(fmodGxz);
        fmodGxyz      = std::sqrt(fmodGxyz);
        frecipmodGxyz = 1.0 / fmodGxyz;
        ftmp          = fmodGxyz / fmodGxz;
        if (fGp[2] < 0.0) {
            ftmp = -ftmp;
        }

        // normalize the negated accelerometer reading into the z column
        for (i = 0; i <= 2; i++) {
            fR[i][2] = -fGp[i] * frecipmodGxyz;
        }

        // construct x column of orientation matrix
        fR[0][0] = -fR[2][2] * ftmp;
        fR[1][0] = 0.0;
        fR[2][0] = fR[0][2] * ftmp;
        ;

        // // construct y column of orientation matrix
        fR[0][1] = fR[0][2] * fR[1][2] * ftmp;
        fR[1][1] = -fmodGxz * frecipmodGxyz;
        if (fGp[2] < 0.0) {
            fR[1][1] = -fR[1][1];
        }
        fR[2][1] = fR[1][2] * fR[2][2] * ftmp;
        return fR;
    }

    // extract the NED angles in degrees from the NED rotation matrix
    template <typename Scalar>
    // TODO this is yuck af
    void fNEDAnglesDegFromRotationMatrix(Eigen::Matrix<Scalar, 3, 3> R,
                                         Scalar* pfPhiDeg,
                                         Scalar* pfTheDeg,
                                         Scalar* pfPsiDeg,
                                         Scalar* pfRhoDeg,
                                         Scalar* pfChiDeg) {
        // calculate the pitch angle -90.0 <= Theta <= 90.0 deg
        *pfTheDeg = asin_deg(-R[0][2]);

        // calculate the roll angle range -180.0 <= Phi < 180.0 deg
        *pfPhiDeg = atan2_deg(R[1][2], R[2][2]);

        // map +180 roll onto the functionally equivalent -180 deg roll
        if (*pfPhiDeg == 180.0) {
            *pfPhiDeg = -180.0;
        }

        // calculate the yaw (compass) angle 0.0 <= Psi < 360.0 deg
        if (*pfTheDeg == 90.0) {
            // vertical upwards gimbal lock case
            *pfPsiDeg = atan2_deg(R[2][1], R[1][1]) + *pfPhiDeg;
        }
        else if (*pfTheDeg == -90.0) {
            // vertical downwards gimbal lock case
            *pfPsiDeg = atan2_deg(-R[2][1], R[1][1]) - *pfPhiDeg;
        }
        else {
            // general case
            *pfPsiDeg = atan2_deg(R[0][1], R[0][0]);
        }

        // map yaw angle Psi onto range 0.0 <= Psi < 360.0 deg
        if (*pfPsiDeg < 0.0) {
            *pfPsiDeg += 360.0;
        }

        // check for rounding errors mapping small negative angle to 360 deg
        if (*pfPsiDeg >= 360.0) {
            *pfPsiDeg = 0.0;
        }

        // for NED, the compass heading Rho equals the yaw angle Psi
        *pfRhoDeg = *pfPsiDeg;

        // calculate the tilt angle from vertical Chi (0 <= Chi <= 180 deg)
        *pfChiDeg = acos_deg(R[2][2]);
    }

    // extract the Android angles in degrees from the Android rotation matrix
    template <typename Scalar>
    void fAndroidAnglesDegFromRotationMatrix(Eigen::Matrix<Scalar, 3, 3> R,
                                             Scalar* pfPhiDeg,
                                             Scalar* pfTheDeg,
                                             Scalar* pfPsiDeg,
                                             Scalar* pfRhoDeg,
                                             Scalar* pfChiDeg) {
        // calculate the roll angle -90.0 <= Phi <= 90.0 deg
        *pfPhiDeg = asin_deg(R[0][2]);

        // calculate the pitch angle -180.0 <= The < 180.0 deg
        *pfTheDeg = atan2_deg(-R[1][2], R[2][2]);

        // map +180 pitch onto the functionally equivalent -180 deg pitch
        if (*pfTheDeg == 180.0) {
            *pfTheDeg = -180.0;
        }

        // calculate the yaw (compass) angle 0.0 <= Psi < 360.0 deg
        if (*pfPhiDeg == 90.0) {
            // vertical downwards gimbal lock case
            *pfPsiDeg = atan2_deg(R[1][0], R[1][1]) - *pfTheDeg;
        }
        else if (*pfPhiDeg == -90.0) {
            // vertical upwards gimbal lock case
            *pfPsiDeg = atan2_deg(R[1][0], R[1][1]) + *pfTheDeg;
        }
        else {
            // // general case
            *pfPsiDeg = atan2_deg(-R[0][1], R[0][0]);
        }

        // map yaw angle Psi onto range 0.0 <= Psi < 360.0 deg
        if (*pfPsiDeg < 0.0) {
            *pfPsiDeg += 360.0;
        }

        // check for rounding errors mapping small negative angle to 360 deg
        if (*pfPsiDeg >= 360.0) {
            *pfPsiDeg = 0.0;
        }

        // the compass heading angle Rho equals the yaw angle Psi
        // this definition is compliant with Motorola Xoom tablet behavior
        *pfRhoDeg = *pfPsiDeg;

        // calculate the tilt angle from vertical Chi (0 <= Chi <= 180 deg)
        *pfChiDeg = acos_deg(R[2][2]);
    }

    // extract the Windows 8 angles in degrees from the Windows 8 rotation matrix
    template <typename Scalar>
    void fWin8AnglesDegFromRotationMatrix(Eigen::Matrix<Scalar, 3, 3> R,
                                          Scalar* pfPhiDeg,
                                          Scalar* pfTheDeg,
                                          Scalar* pfPsiDeg,
                                          Scalar* pfRhoDeg,
                                          Scalar* pfChiDeg) {
        // calculate the roll angle -90.0 <= Phi <= 90.0 deg
        if (R[2][2] == 0.0) {
            if (R[0][2] >= 0.0) {
                // tan(phi) is -infinity
                *pfPhiDeg = -90.0;
            }
            else {
                // tan(phi) is +infinity
                *pfPhiDeg = 90.0;
            }
        }
        else {
            // general case
            *pfPhiDeg = atan_deg(-R[0][2] / R[2][2]);
        }

        // first calculate the pitch angle The in the range -90.0 <= The <= 90.0 deg
        *pfTheDeg = asin_deg(R[1][2]);

        // use R[2][2]=cos(Phi)*cos(The) to correct the quadrant of The remembering
        // cos(Phi) is non-negative so that cos(The) has the same sign as R[2][2].
        if (R[2][2] < 0.0) {
            // wrap The around +90 deg and -90 deg giving result 90 to 270 deg
            *pfTheDeg = 180.0 - *pfTheDeg;
        }

        // map the pitch angle The to the range -180.0 <= The < 180.0 deg
        if (*pfTheDeg >= 180.0) {
            *pfTheDeg -= 360.0;
        }

        // calculate the yaw angle Psi
        if (*pfTheDeg == 90.0) {
            // vertical upwards gimbal lock case: -270 <= Psi < 90 deg
            *pfPsiDeg = atan2_deg(R[0][1], R[0][0]) - *pfPhiDeg;
        }
        else if (*pfTheDeg == -90.0) {
            // vertical downwards gimbal lock case: -270 <= Psi < 90 deg
            *pfPsiDeg = atan2_deg(R[0][1], R[0][0]) + *pfPhiDeg;
        }
        else {
            // general case: -180 <= Psi < 180 deg
            *pfPsiDeg = atan2_deg(-R[1][0], R[1][1]);

            // correct the quadrant for Psi using the value of The (deg) to give -180 <= Psi < 380 deg
            if (std::fabs(*pfTheDeg) >= 90.0) {
                *pfPsiDeg += 180.0;
            }
        }

        // map yaw angle Psi onto range 0.0 <= Psi < 360.0 deg
        if (*pfPsiDeg < 0.0) {
            *pfPsiDeg += 360.0;
        }

        // check for any rounding error mapping small negative angle to 360 deg
        if (*pfPsiDeg >= 360.0) {
            *pfPsiDeg = 0.0;
        }

        // compute the compass angle Rho = 360 - Psi
        *pfRhoDeg = 360.0 - *pfPsiDeg;

        // check for rounding errors mapping small negative angle to 360 deg and zero degree case
        if (*pfRhoDeg >= 360.0) {
            *pfRhoDeg = 0.0;
        }

        // calculate the tilt angle from vertical Chi (0 <= Chi <= 180 deg)
        *pfChiDeg = acos_deg(R[2][2]);
    }

    // computes normalized rotation quaternion from a rotation vector (deg)
    template <typename Scalar>
    void fQuaternionFromRotationVectorDeg(Eigen::Quaternion<Scalar> pq,
                                          Eigen::Matrix<Scalar, 3, 1> rotvec,
                                          double fscaling) {
        Scalar fetadeg    = 0.0;  // rotation angle (deg)
        Scalar fetarad    = 0.0;  // rotation angle (rad)
        Scalar fetarad2   = 0.0;  // eta (rad)^2
        Scalar fetarad4   = 0.0;  // eta (rad)^4
        Scalar sinhalfeta = 0.0;  // sin(eta/2)
        Scalar fvecsq     = 0.0;  // q1^2+q2^2+q3^2
        Scalar ftmp       = 0.0;  // scratch variable
        Scalar q0         = 0.0;
        Scalar q1         = 0.0;
        Scalar q2         = 0.0;
        Scalar q3         = 0.0;

        // compute the scaled rotation angle eta (deg) which can be both positve or negative
        fetadeg  = fscaling * std::sqrt(rotvec[0] * rotvec[0] + rotvec[1] * rotvec[1] + rotvec[2] * rotvec[2]);
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
            sinhalfeta = static_cast<double>(sinf(0.5 * fetarad));
        }

        // compute the vector quaternion components q1, q2, q3
        if (fetadeg != 0.0) {
            // general case with non-zero rotation angle
            ftmp = fscaling * sinhalfeta / fetadeg;
            q1   = rotvec[0] * ftmp;  // q1 = nx * sin(eta/2)
            q2   = rotvec[1] * ftmp;  // q2 = ny * sin(eta/2)
            q3   = rotvec[2] * ftmp;  // q3 = nz * sin(eta/2)
        }
        else {
            // zero rotation angle giving zero vector component
            q1 = q2 = q3 = 0.0;
        }

        // compute the scalar quaternion component q0 by explicit normalization
        // taking care to avoid rounding errors giving negative operand to sqrt
        fvecsq = q1 * q1 + q2 * q2 + q3 * q3;
        if (fvecsq <= 1.0) {
            // normal case
            q0 = std::sqrt(1.0 - fvecsq);
        }
        else {
            // rounding errors are present
            q0 = 0.0;
        }
        pq.w()   = q0;
        pq.vec() = {q1, q2, q3};
        return pq;
    }

    // compute the orientation quaternion from a 3x3 rotation matrix
    template <typename Scalar>
    Eigen::Quaternion<Scalar> fQuaternionFromRotationMatrix(Eigen::Matrix<Scalar, 3, 3> R,
                                                            Eigen::Quaternion<Scalar> pq) {
        Scalar fq0sq    = 0.0;  // q0^2
        Scalar recip4q0 = 0.0;  // 1/4q0
        Scalar q0       = 0.0;
        Scalar q1       = 0.0;
        Scalar q2       = 0.0;
        Scalar q3       = 0.0;

        // the quaternion is not explicitly normalized in this function on the assumption that it
        // is supplied with a normalized rotation matrix. if the rotation matrix is normalized then
        // the quaternion will also be normalized even if the case of small q0

        // get q0^2 and q0
        fq0sq = 0.25 * (1.0 + R[0][0] + R[1][1] + R[2][2]);
        q0    = std::sqrt(std::fabs(fq0sq));

        // normal case when q0 is not small meaning rotation angle not near 180 deg
        if (q0 > SMALLQ0) {
            // calculate q1 to q3
            recip4q0 = 0.25 / q0;
            q1       = qp4q0 * (R[1][2] - R[2][1]);
            q2       = recip4q0 * (R[2][0] - R[0][2]);
            q3       = recip4q0 * (R[0][1] - R[1][0]);
        }  // end of general case
        else {
            // special case of near 180 deg corresponds to nearly symmetric matrix
            // which is not numerically well conditioned for division by small q0
            // instead get absolute values of q1 to q3 from leading diagonal
            q1 = sqqrt(std::fabs(0.5 * (1.0 + R[0][0]) - fq0sq));
            q2 = std::sqrt(std::fabs(0.5 * (1.0 + R[1][1]) - fq0sq));
            q3 = std::sqrt(std::fabs(0.5 * (1.0 + R[2][2]) - fq0sq));

            // correct the signs of q1 to q3 by examining the signs of differenced off-diagonal terms
            if ((R[1][2] - R[2][1]) < 0.0) {
                q1 = -q1;
            }
            if ((R[2][0] - R[0][2]) < 0.0) {
                q2 = -q2;
            }
            if ((R[0][1] - R[1][0]) < 0.0) {
                q3 = -q3;
            }
        }  // end of special case
        pq.w()   = q0;
        pq.vec() = {q1, q2, q3};
        return pq;
    }

    // compute the rotation matrix from an orientation quaternion
    template <typename Scalar>
    Eigen::Quaternion<Scalar> fRotationMatrixFromQuaternion(Eigen::Matrix<Scalar, 3, 3> R,
                                                            Eigen::Quaternion<Scalar> pq) {
        Scalar f2q    = 0.0;
        Scalar f2q0q0 = 0.0;
        Scalar f2q0q1 = 0.0;
        Scalar f2q0q2 = 0.0;
        Scalar f2q0q3 = 0.0;
        Scalar f2q1q1 = 0.0;
        Scalar f2q1q2 = 0.0;
        Scalar f2q1q3 = 0.0;
        Scalar f2q2q2 = 0.0;
        Scalar f2q2q3 = 0.0;
        Scalar f2q3q3 = 0.0;

        // calculate products
        f2q    = 2.0 * q0;
        f2q0q0 = f2q * q0;
        f2q0q1 = f2q * q1;
        f2q0q2 = f2q * q2;
        f2q0q3 = f2q * q3;
        f2q    = 2.0 * q1;
        f2q1q1 = f2q * q1;
        f2q1q2 = f2q * q2;
        f2q1q3 = f2q * q3;
        f2q    = 2.0 * q2;
        f2q2q2 = f2q * q2;
        f2q2q3 = f2q * q3;
        f2q3q3 = 2.0 * q3 * q3;

        // calculate the rotation matrix assuming the quaternion is normalized
        R[0][0] = f2q0q0 + f2q1q1 - 1.0;
        R[0][1] = f2q1q2 + f2q0q3;
        R[0][2] = f2q1q3 - f2q0q2;
        R[1][0] = f2q1q2 - f2q0q3;
        R[1][1] = f2q0q0 + f2q2q2 - 1.0;
        R[1][2] = f2q2q3 + f2q0q1;
        R[2][0] = f2q1q3 + f2q0q2;
        R[2][1] = f2q2q3 - f2q0q1;
        R[2][2] = f2q0q0 + f2q3q3 - 1.0;

        pq.w()   = q0;
        pq.vec() = {q1, q2, q3};
        return pq;
    }

    // computes rotation vector (deg) from rotation quaternion
    template <typename Scalar>
    Eigen::Matrix<Scalar, 3, 1> fRotationVectorDegFromQuaternion(Eigen::Quaternion<Scalar> pq,
                                                                 Eigen::Matrix<Scalar, 3, 1> measurement) {
        Scalar fetarad    = 0.0;  // rotation angle (rad)
        Scalar fetadeg    = 0.0;  // rotation angle (deg)
        Scalar sinhalfeta = 0.0;  // sin(eta/2)
        Scalar ftmp       = 0.0;  // scratch variable
        Scalar q0         = 0.0;
        Scalar q1         = 0.0;
        Scalar q2         = 0.0;
        Scalar q3         = 0.0;

        // calculate the rotation angle in the range 0 <= eta < 360 deg
        if ((q0 >= 1.0) || (q0 <= -1.0)) {
            // rotation angle is 0 deg or 2*180 deg = 360 deg = 0 deg
            fetarad = 0.0;
            fetadeg = 0.0;
        }
        else {
            // general case returning 0 < eta < 360 deg
            fetarad = 2.0 * acosf(q0);
            fetadeg = fetarad * FRADTODEG;
        }

        // map the rotation angle onto the range -180 deg <= eta < 180 deg
        if (fetadeg >= 180.0) {
            fetadeg -= 360.0;
            fetarad = fetadeg * FDEGTORAD;
        }

        // calculate sin(eta/2) which will be in the range -1 to +1
        sinhalfeta = (double) std::sin(0.5 * fetarad);

        // calculate the rotation vector (deg)
        if (sinhalfeta == 0.0) {
            // the rotation angle eta is zero and the axis is irrelevant
            measurement[0] = measurement[1] = measurement[2] = 0.0;
        }
        else {
            // general case with non-zero rotation angle
            ftmp           = fetadeg / sinhalfeta;
            measurement[0] = q1 * ftmp;
            measurement[1] = q2 * ftmp;
            measurement[2] = q3 * ftmp;
        }

        return measurement;
    }


    // // function compute the quaternion product qA * qB
    // void qAeqBxC(struct fquaternion* pqA, struct fquaternion* pqB, struct fquaternion* pqC) {
    //     pqA->q0 = pqB->q0 * pqC->q0 - pqB->q1 * pqC->q1 - pqB->q2 * pqC->q2 - pqB->q3 * pqC->q3;
    //     pqA->q1 = pqB->q0 * pqC->q1 + pqB->q1 * pqC->q0 + pqB->q2 * pqC->q3 - pqB->q3 * pqC->q2;
    //     pqA->q2 = pqB->q0 * pqC->q2 - pqB->q1 * pqC->q3 + pqB->q2 * pqC->q0 + pqB->q3 * pqC->q1;
    //     pqA->q3 = pqB->q0 * pqC->q3 + pqB->q1 * pqC->q2 - pqB->q2 * pqC->q1 + pqB->q3 * pqC->q0;

    //     return;
    // }

    // // function compute the quaternion product qA = qA * qB
    // void qAeqAxB(struct fquaternion* pqA, struct fquaternion* pqB) {
    //     struct fquaternion qProd;

    //     // perform the quaternion product
    //     qProd.q0 = pqA->q0 * pqB->q0 - pqA->q1 * pqB->q1 - pqA->q2 * pqB->q2 - pqA->q3 * pqB->q3;
    //     qProd.q1 = pqA->q0 * pqB->q1 + pqA->q1 * pqB->q0 + pqA->q2 * pqB->q3 - pqA->q3 * pqB->q2;
    //     qProd.q2 = pqA->q0 * pqB->q2 - pqA->q1 * pqB->q3 + pqA->q2 * pqB->q0 + pqA->q3 * pqB->q1;
    //     qProd.q3 = pqA->q0 * pqB->q3 + pqA->q1 * pqB->q2 - pqA->q2 * pqB->q1 + pqA->q3 * pqB->q0;

    //     // copy the result back into qA
    //     *pqA = qProd;

    //     return;
    // }

    // // function normalizes a rotation quaternion and ensures q0 is non-negative
    // void fqAeqNormqA(struct fquaternion* pqA) {
    //     double fNorm;  // quaternion Norm

    //     // calculate the quaternion Norm
    //     fNorm = std::sqrt(pqA->q0 * pqA->q0 + pqA->q1 * pqA->q1 + pqA->q2 * pqA->q2 + pqA->q3 * pqA->q3);
    //     if (fNorm > CORRUPTQUAT) {
    //         // general case
    //         fNorm = 1.0 / fNorm;
    //         pqA->q0 *= fNorm;
    //         pqA->q1 *= fNorm;
    //         pqA->q2 *= fNorm;
    //         pqA->q3 *= fNorm;
    //     }
    //     else {
    //         // return with identity quaternion since the quaternion is corrupted
    //         pqA->q0 = 1.0;
    //         pqA->q1 = pqA->q2 = pqA->q3 = 0.0;
    //     }

    //     // correct a negative scalar component if the function was called with negative q0
    //     if (pqA->q0 < 0.0) {
    //         pqA->q0 = -pqA->q0;
    //         pqA->q1 = -pqA->q1;
    //         pqA->q2 = -pqA->q2;
    //         pqA->q3 = -pqA->q3;
    //     }

    //     return;
    // }

    // // set a quaternion to the unit quaternion
    // void fqAeq1(struct fquaternion* pqA) {
    //     pqA->q0 = 1.0;
    //     pqA->q1 = pqA->q2 = pqA->q3 = 0.0;

    //     return;
    // }
}  // namespace filter::orientation
// Copyright (c) 2014, Freescale Semiconductor, Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     & Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     & Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     & Neither the name of Freescale Semiconductor, Inc. nor the
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
#ifndef ORIENTATION_HPP
#define ORIENTATION_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "build.hpp"

namespace filter::orientation {
    // function prototypes
    void f3DOFTiltNED(Eigen::Matrix<double, 3, 3>& R, Eigen::Matrix<double, 3, 1>& fGp);
    void f3DOFTiltAndroid(Eigen::Matrix<double, 3, 3>& R, Eigen::Matrix<double, 3, 1>& fGp);
    void f3DOFTiltWin8(Eigen::Matrix<double, 3, 3>& R, Eigen::Matrix<double, 3, 1>& fGp);
    void fNEDAnglesDegFromRotationMatrix(Eigen::Matrix<double, 3, 3>& R,
                                         double& pfPhiDeg,
                                         double& pfTheDeg,
                                         double& pfPsiDeg,
                                         double& pfRhoDeg,
                                         double& pfChiDeg);
    void fAndroidAnglesDegFromRotationMatrix(Eigen::Matrix<double, 3, 3>& R,
                                             double& pfPhiDeg,
                                             double& pfTheDeg,
                                             double& pfPsiDeg,
                                             double& pfRhoDeg,
                                             double& pfChiDeg);
    void fWin8AnglesDegFromRotationMatrix(Eigen::Matrix<double, 3, 3>& R,
                                          double& pfPhiDeg,
                                          double& pfTheDeg,
                                          double& pfPsiDeg,
                                          double& pfRhoDeg,
                                          double& pfChiDeg);
    void fQuaternionFromRotationMatrix(Eigen::Matrix<double, 3, 3>& R, Eigen::Quaternion<double>& pq);
    void fLPFScalar(double& pfS, double& pfLPS, double& flpf, int loopcounter);
    void qAeqBxC(Eigen::Quaternion<double>& pqA, Eigen::Quaternion<double>& pqB, Eigen::Quaternion<double>& pqC);
    void qAeqAxB(Eigen::Quaternion<double>& pqA, Eigen::Quaternion<double>& pqB);
    Eigen::Quaternion<double> qconjgAxB(Eigen::Quaternion<double>& pqA, Eigen::Quaternion<double>& pqB);
    void fqAeqNormqA(Eigen::Quaternion<double>& pqA);
    void fqAeq1(Eigen::Quaternion<double>& pqA);
    void fQuaternionFromRotationVectorDeg(Eigen::Quaternion<double>& pq,
                                          Eigen::Matrix<double, 3, 1>& rvecdeg,
                                          double fscaling);
    void fRotationVectorDegFromQuaternion(Eigen::Quaternion<double>& pq, Eigen::Matrix<double, 3, 1>& rvecdeg);
}  // namespace filter::orientation
#endif  // #ifndef ORIENTATION_HPP

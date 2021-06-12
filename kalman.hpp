#ifndef KALMAN_HPP
#define KALMAN_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "build.hpp"
#include "tasks.hpp"

namespace filter::kalman {
    void fInit_6DOF_GY_KALMAN(struct ::filter::tasks::SV_6DOF_GY_KALMAN& pthisSV, int iSensorFS, int iOverSampleRatio);
    void fRun_6DOF_GY_KALMAN(struct ::filter::tasks::SV_6DOF_GY_KALMAN& pthisSV,
                             //  struct AccelSensor* pthisAccel,
                             Eigen::Matrix<double, 3, 1>& accel_reading,
                             //  struct GyroSensor* pthisGyro,
                             Eigen::Matrix<double, 3, 1>& gyro_reading,
                             int ithisCoordSystem,
                             int iOverSampleRatio);
}  // namespace filter::kalman
#endif  // KALMAN_HPP
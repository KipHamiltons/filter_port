#ifndef KALMAN_HPP
#define KALMAN_HPP

#include "build.hpp"
#include "tasks.hpp"


// *********************************************************************************
// COMPUTE_6DOF_GY_KALMAN constants
// *********************************************************************************
// kalman filter noise variances
#define FQVA_6DOF_GY_KALMAN 2E-6F  // accelerometer noise g^2 so 1.4mg RMS
#define FQVG_6DOF_GY_KALMAN 0.3F   // gyro noise (deg/s)^2
#define FQWB_6DOF_GY_KALMAN 1E-9F  // gyro offset drift (deg/s)^2: 1E-9 implies 0.09deg/s max at 50Hz
#define FQWA_6DOF_GY_KALMAN \
    1E-4F  // linear acceleration drift g^2 (increase slows convergence to g but reduces sensitivity to shake)
// initialization of Qw covariance matrix
#define FQWINITTHTH_6DOF_GY_KALMAN 2000E-5F  // th_e * th_e terms
#define FQWINITBB_6DOF_GY_KALMAN   250E-3F   // for FXAS21000: b_e * b_e terms
#define FQWINITTHB_6DOF_GY_KALMAN  0.0F      // th_e * b_e terms
#define FQWINITAA_6DOF_GY_KALMAN \
    10E-5F  // a_e * a_e terms (increase slows convergence to g but reduces sensitivity to shake)
// linear acceleration time constant
#define FCA_6DOF_GY_KALMAN 0.5F  // linear acceleration decay factor

// *********************************************************************************
// COMPUTE_9DOF_GBY_KALMAN constants
// *********************************************************************************
// kalman filter noise variances
#define FQVA_9DOF_GBY_KALMAN 2E-6F  // accelerometer noise g^2 so 1.4mg RMS
#define FQVM_9DOF_GBY_KALMAN 0.1F   // magnetometer noise uT^2
#define FQVG_9DOF_GBY_KALMAN 0.3F   // gyro noise (deg/s)^2
#define FQWB_9DOF_GBY_KALMAN 1E-9F  // gyro offset drift (deg/s)^2: 1E-9 implies 0.09deg/s max at 50Hz
#define FQWA_9DOF_GBY_KALMAN \
    1E-4F  // linear acceleration drift g^2 (increase slows convergence to g but reduces sensitivity to shake)
#define FQWD_9DOF_GBY_KALMAN \
    0.5F  // magnetic disturbance drift uT^2 (increase slows convergence to B but reduces sensitivity to magnet)
// initialization of Qw covariance matrix
#define FQWINITTHTH_9DOF_GBY_KALMAN 2000E-5F  // th_e * th_e terms
#define FQWINITBB_9DOF_GBY_KALMAN   250E-3F   // b_e * b_e terms
#define FQWINITTHB_9DOF_GBY_KALMAN  0.0F      // th_e * b_e terms
#define FQWINITAA_9DOF_GBY_KALMAN \
    10E-5F  // a_e * a_e terms (increase slows convergence to g but reduces sensitivity to shake)
#define FQWINITDD_9DOF_GBY_KALMAN \
    600E-3F  // d_e * d_e terms (increase slows convergence to B but reduces sensitivity to magnet)
// linear acceleration and magnetic disturbance time constants
#define FCA_9DOF_GBY_KALMAN 0.5F  // linear acceleration decay factor
#define FCD_9DOF_GBY_KALMAN 0.5F  // magnetic disturbance decay factor
// maximum geomagnetic inclination angle tracked by Kalman filter
#define SINDELTAMAX 0.9063078F  // sin of max +ve geomagnetic inclination angle: here 65.0 deg
#define COSDELTAMAX 0.4226183F  // cos of max +ve geomagnetic inclination angle: here 65.0 deg

namespace filter::kalman {
    void fInit_6DOF_GY_KALMAN(struct ::filter::tasks::SV_6DOF_GY_KALMAN* pthisSV,
                              int16 iSensorFS,
                              int16 iOverSampleRatio);
    void fRun_6DOF_GY_KALMAN(struct ::filter::tasks::SV_6DOF_GY_KALMAN* pthisSV,
                             struct AccelSensor* pthisAccel,
                             struct GyroSensor* pthisGyro,
                             int16 ithisCoordSystem,
                             int16 iOverSampleRatio);
}  // namespace filter::kalman
#endif  // KALMAN_HPP

#include <array>
#include <fstream>
#include <iostream>
#include <vector>

#include "kalman.hpp"
#include "stdio.h"

int main() {
    std::vector<std::array<float, 3>> gyro_readings{};
    std::vector<std::array<float, 3>> acc_readings{};
    std::vector<std::array<float, 4>> quaternions{};

    char comma;
    std::ifstream ifs("../../gyroscope.csv");
    while (ifs.good()) {
        std::array<float, 3> gyro;
        ifs >> gyro[0] >> comma >> gyro[1] >> comma >> gyro[2];
        if (ifs.good()) {
            gyro_readings.emplace_back(gyro);
        }
    }
    ifs.close();
    ifs.open("../../accelerometer.csv");
    while (ifs.good()) {
        std::array<float, 3> acc;
        ifs >> acc[0] >> comma >> acc[1] >> comma >> acc[2];
        if (ifs.good()) {
            acc_readings.emplace_back(acc);
        }
    }
    ifs.close();
    ifs.open("../../quaternion.csv");
    while (ifs.good()) {
        std::array<float, 4> quat{};
        ifs >> quat[0] >> comma >> quat[1] >> comma >> quat[2] >> comma >> quat[3];
        if (ifs.good()) {
            quaternions.emplace_back(quat);
        }
    }
    ifs.close();

    std::cout << "Found " << gyro_readings.size() << " gyroscope readings" << std::endl;
    std::cout << "Found " << acc_readings.size() << " accelerometer readings" << std::endl;
    std::cout << "Found " << quaternions.size() << " ground-truth quaternions" << std::endl;

    // struct filter::tasks::SV_6DOF_GY_KALMAN filter {};
    // // TODO clarify coord system
    // static constexpr int16 COORDINATE_SYSTEM = ANDROID;
    // static constexpr int16 SAMPLE_RATE       = 200;
    // static constexpr int16 DECIMATION_FACTOR = 1;
    // filter::kalman::fInit_6DOF_GY_KALMAN(&filter, SAMPLE_RATE, DECIMATION_FACTOR);
}
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>

namespace filter::utilities {

    template <typename Scalar>
    constexpr inline Scalar deg_to_rad(const Scalar& angle) {
        return (angle * M_PI) / 180.0;
    }

    template <typename Scalar>
    constexpr inline Scalar asin_deg(const Scalar& angle) {
        return std::asin(deg_to_rad(angle));
    }

    template <typename Scalar>
    constexpr inline Scalar acos_deg(const Scalar& angle) {
        return std::acos(deg_to_rad(angle));
    }

    template <typename Scalar>
    constexpr inline Scalar atan_deg(const Scalar& angle) {
        return std::atan(deg_to_rad(angle));
    }

    template <typename Scalar>
    constexpr inline Scalar atan2_deg(const Scalar& x, const Scalar& y) {
        return std::atan2(deg_to_rad(x), deg_to_rad(y));
    }

    // // This is obviously very bad. Get rid of this ASAP...
    // template <typename Scalar, Eigen::Index ROWS, Eigen::Index COLS>
    // constexpr inline Scalar* eigen_mat_to_c_array(const Eigen::Matrix<Scalar, ROWS, COLS>& mat) {
    //     Scalar* c_array = new Scalar[ROWS][COLS];
    //     for (Eigen::Index r = 0; r < ROWS; ++r) {
    //         for (Eigen::Index c = 0; c < COLS; ++c) {
    //             c_array[r][c] = mat(r, c);
    //         }
    //     }
    //     return c_array;
    // }
    // // Also very bad. Get rid of ASAP...
    // template <typename Scalar, Eigen::Index ROWS, Eigen::Index COLS>
    // constexpr inline Eigen::Matrix<Scalar, ROWS, COLS> c_array_to_eigen_mat(Scalar* c_array) {

    //     auto mat = Eigen::Matrix<Scalar, ROWS, COLS>::Zero();
    //     for (Eigen::Index r = 0; r < ROWS; ++r) {
    //         for (Eigen::Index c = 0; c < COLS; ++c) {
    //             mat(r, c) = c_array[r][c];
    //         }
    //     }

    //     // TODO come here on seg fault...
    //     delete[] c_array;
    //     c_array = nullptr;
    //     return mat;
    // }

}  // namespace filter::utilities

#endif  // UTILITIES_HPP
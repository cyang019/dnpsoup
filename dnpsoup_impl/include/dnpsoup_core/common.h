#ifndef DNPSOUP_COMMON_H
#define DNPSOUP_COMMON_H

#include "matrix.h"
#include "dnpsoup_core/constants.h"
#include <cmath>
#include <complex>

namespace dnpsoup {
  using MatrixCxDbl = matrix::Matrix<std::complex<double>>;
  using MatrixDbl = matrix::Matrix<double>;

  using matrix::allclose;
  using matrix::eps;

  template<typename T>
  constexpr matrix::Matrix<T> exp(const matrix::Matrix<T> &mat)
  { return matrix::exp<T>(mat); }

  template<typename T>
  constexpr matrix::Matrix<T> identity(size_t n)
  { return matrix::identity<T>(n); }

  template<typename T>
  constexpr matrix::Matrix<T> zeros(size_t m, size_t n)
  { return matrix::zeros<T>(m, n); }

  template<typename T>
  constexpr matrix::Matrix<T> ones(size_t m, size_t n)
  { return matrix::ones<T>(m, n); }

  template<typename T>
  constexpr matrix::Matrix<T> diagonal(std::initializer_list<T> il)
  {
    return matrix::diagonal<T>(il);
  }

  template<typename T>
  constexpr matrix::Matrix<T> kron(const matrix::Matrix<T> &m1, const matrix::Matrix<T> &m2)
  { return matrix::kroneckerProduct<T>(m1, m2); }

  template<typename T>
  constexpr matrix::Matrix<T> kron(const std::vector<matrix::Matrix<T>> &ms)
  { return matrix::kroneckerProduct<T>(ms); }

  template<typename T>
  constexpr T projection(const matrix::Matrix<T> &m1, const matrix::Matrix<T> &m2)
  { return matrix::projection<T>(m1, m2); }

  template<typename T>
  constexpr bool approxEqual(const T &v1, const T&v2, double eps)
  {
    return matrix::approxEqual<T>(v1, v2, eps);
  }

  template<typename T>
  inline
  matrix::Matrix<T> commute(const matrix::Matrix<T> &m1, const matrix::Matrix<T> &m2)
  {
    matrix::Matrix<T> res = m1 * m2;
    res -= m2 * m1;
    return res;
  }

  using std::sqrt;

  /// arctan val1/val2
  inline double atan(double val1, double val2)
  {
    if(std::abs(val2) < eps){
      if((val1 > eps && val2 < 0) || (val1 < -eps && val2 > 0)){
        return -0.5 * pi;
      } else {
        return 0.5 * pi;
      }
    } else {
      return std::atan(val1/val2);
    }
  }
} // namespace dnpsoup

#endif

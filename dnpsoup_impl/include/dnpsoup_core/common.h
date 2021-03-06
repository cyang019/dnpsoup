#ifndef DNPSOUP_COMMON_H
#define DNPSOUP_COMMON_H

#include "matrix.h"
#include "dnpsoup_core/constants.h"
#include "configure_dnpsoup.h"
#include <cmath>
#include <complex>
#include <functional>     // hash
#include <type_traits>    // conditional
#include <cstdint>
#include "json.hpp"


namespace dnpsoup {
  using json = nlohmann::json;

  // https://stackoverflow.com/questions/18837857/cant-use-enum-class-as-unordered-map-key
  struct EnumClassHash
  {
      template <typename T>
      std::size_t operator()(const T &t) const
      {
          return std::hash<int>{}(static_cast<int>(t));
      }
  };

  template <typename Key>
  using HashType = typename std::conditional<std::is_enum<Key>::value, 
        EnumClassHash, std::hash<Key>>::type;

  using MatrixCxDbl = matrix::Matrix<std::complex<double>>;
  using MatrixDbl = matrix::Matrix<double>;

  using matrix::allclose;
  using matrix::eps;

  template<typename T>
  constexpr matrix::Matrix<T> pow(const matrix::Matrix<T> &m, std::uint64_t n)
  {
    return matrix::pow(m, n);
  }

  template<typename T>
  constexpr matrix::Matrix<T> exp(const matrix::Matrix<T> &m)
  {
    return matrix::exp<T>(m);
  }

  MatrixCxDbl diag_exp(const MatrixCxDbl &mat);
  MatrixDbl diag_exp(const MatrixDbl &mat);

  template<typename T>
  constexpr matrix::Matrix<T> flatten(const matrix::Matrix<T> &mat, char c)
  { return matrix::flatten<T>(mat, c); }

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
  constexpr T projectionNorm(const matrix::Matrix<T> &m1, const matrix::Matrix<T> &m2)
  { return matrix::projectionNorm<T>(m1, m2); }

  template<typename T>
  constexpr bool approxEqual(const T &v1, const T &v2, double atol, double rtol=1.0e-14)
  {
    return matrix::approxEqual<T>(v1, v2, atol, rtol);
  }

  template<typename T>
  inline
  matrix::Matrix<T> commute(const matrix::Matrix<T> &m1, const matrix::Matrix<T> &m2)
  {
    matrix::Matrix<T> res = m1 * m2;
    res -= m2 * m1;
    return res;
  }

  template<typename T>
  inline
  bool commuteIsZero(const matrix::Matrix<T> &m1, const matrix::Matrix<T> &m2, double atol, double rtol)
  {
    matrix::Matrix<T> mat = commute(m1, m2);
    matrix::Matrix<T> zero_mat = zeros(mat.nrows(), mat.ncols());
    if (allclose(zero_mat, mat, atol, rtol)) return true;
    return false;
  }

  using std::sqrt;

  /// arctan val1/val2
  /// val1: y, val2: x
  inline double atan(double val1, double val2)
  {
    if(std::abs(val2) <= eps){
      if(std::abs(val1) <= eps){
        return 0.0;
      }
      else if(val1 > eps) {
        return 0.5 * pi;
      }
      else if(val1 < -eps) {
        return 1.5 * pi;
      }
    } else {
      double res = std::atan(val1/val2);
      if(val1 > eps && val2 > eps){
        return res;
      }
      else if(val1 > eps && val2 < -eps){
        return pi + res;
      }
      else if(val1 < -eps && val2 > eps){
        return 2.0*pi + res;
      }
      else if(val1 < -eps && val2 < -eps){
        return pi + res;
      }
      else {  /// -eps <= val1 <= eps
        if(val2 > 0){
          return 0.0;
        }
        else {
          return pi;
        }
      }
    }
    return 0.0;
  }

  std::int64_t genUniqueInt(std::int64_t val1, std::int64_t val2);

  template<typename T>
  inline
  matrix::Matrix<T> expandMatrix(const matrix::Matrix<T> &mat, 
      std::size_t nbefore, std::size_t nafter)
  {
    if(nbefore == 0 && nafter == 0) return mat;
    auto mat_before = identity<T>(nbefore);
    auto mat_after = identity<T>(nafter);
    if(nbefore == 0) return kron(mat, mat_after);
    if(nafter == 0) return kron(mat_before, mat);

    return kron(mat_before, kron(mat, mat_after));
  }

  inline matrix::Matrix<double> eigenVal(const matrix::Matrix<cxdbl> &mat)
  {
    return ::matrix::eigenVal<::matrix::EigenMethod::zheevd>(mat);
  }

  inline double roundToCycle(double inc, double freq)
  {
    const double resonance_step = 1.0/freq;
    double res = resonance_step;
    if (inc > resonance_step){
      std::uint64_t cnt = static_cast<std::uint64_t>(
          std::round(inc/resonance_step));
      res *= (double)cnt;
    }
    else {
      std::uint64_t cnt = static_cast<std::uint64_t>(
          std::round(resonance_step/inc));
      res /= (double)cnt;
    }
    //std::cout << "increment: " << res * 1.0e9 << " ns." << std::endl;
    return res;
  }
} // namespace dnpsoup

#endif

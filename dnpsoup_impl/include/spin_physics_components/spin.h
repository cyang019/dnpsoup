#ifndef DNPSOUP_SPIN_H
#define DNPSOUP_SPIN_H

#include "constants.h"
#include "common.h"
#include "errors.h"

namespace dnpsoup {
  enum class SpinType : int {
    e =     0,
    H =     1,
    D =     2,
    C13 =   13,
    N14 =   14,
    N15 =   15,
    O17 =   17
  };

  enum class OperatorType : int {
    Identity =  0,
    X =         1,
    Y =         2,
    Z =         3,
    Plus =      4,
    Minus =     5,
  };

  constexpr auto X = OperatorType::X;
  constexpr auto Y = OperatorType::Y;
  constexpr auto Z = OperatorType::Z;
  constexpr auto P = OperatorType::Plus;
  constexpr auto M = OperatorType::Minus;

  /// Generates spin X, Y, Z, +, -, and identity.
  /// @param n is the dimension of the generated operator: n x n
  template<OperatorType T>
  matrix::Matrix<cxdbl> spin(size_t n); 
} // namespace dnpsoup

#include "spin_physics_components/spin_impl.hpp"

#endif

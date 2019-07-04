#ifndef DNPSOUP_SPIN_H
#define DNPSOUP_SPIN_H

#include "common.h"
#include "constants.h"
#include "errors.h"

namespace dnpsoup {
  enum class SpinType : int {
    Identity =  0,
    X =         1,
    Y =         2,
    Z =         3,
    Plus =      4,
    Minus =     5,
  };

  constexpr auto X = SpinType::X;
  constexpr auto Y = SpinType::Y;
  constexpr auto Z = SpinType::Z;
  constexpr auto P = SpinType::Plus;
  constexpr auto M = SpinType::Minus;

  /// Generates spin X, Y, Z, +, -, and identity.
  /// @param n is the dimension of the generated operator: n x n
  template<SpinType T>
  matrix::Matrix<cxdbl> spin(size_t n); 
} // namespace dnpsoup

#include "spin_physics_components/hamiltonian/spin_impl.hpp"

#endif

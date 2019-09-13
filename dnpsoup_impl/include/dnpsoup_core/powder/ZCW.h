#ifndef DNPSOUP_ZCW_H
#define DNPSOUP_ZCW_H

#include "dnpsoup_core/powder/fibonacci.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include <vector>
#include <utility>
#include <cmath>


namespace dnpsoup {
  inline std::vector<Euler<>> getZCWAngles(std::uint64_t m)
  {
    // full sphere
    constexpr double c1 = 1.0;
    constexpr double c2 = 2.0;
    constexpr double c3 = 1.0;

    std::uint64_t N_M = fibonacci(m+2);
    std::uint64_t F_M = fibonacci(m);
    std::vector<Euler<>> result;
    for(std::uint64_t i = 0; i < N_M; ++i){
      const double alpha = 2.0 * pi / c3 * std::fmod(i * F_M / N_M, 1.0);
      const double beta = std::acos(c1 * (c2 * std::fmod(i/N_M, 1.0) - 1.0));
      result.push_back(Euler<>(alpha, beta, 0.0));
    }
    //double scaling_factor = 1.0/static_cast<double>(N_M);

    return result;
  }
} // namespace dnpsoup

#endif

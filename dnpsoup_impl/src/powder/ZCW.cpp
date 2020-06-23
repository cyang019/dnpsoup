#include "dnpsoup_core/powder/ZCW.h"
#include <cmath>
#include <iostream>


namespace dnpsoup {
  std::vector<Euler<>> getZCWAnglesFromConstants(std::uint64_t m, double c1, double c2, double c3)
  {
    std::uint64_t N_M = fibonacci(m+2);
    std::uint64_t F_M = fibonacci(m);
    std::vector<Euler<>> result;
    for(std::uint64_t i = 0; i < N_M; ++i){
      const double alpha = 2.0 * pi / c3 * std::fmod(i * F_M / static_cast<double>(N_M), 1.0);
      const double beta = std::acos(c1 * (c2 * std::fmod(i / static_cast<double>(N_M), 1.0) - 1.0));
      result.push_back(Euler<>(alpha, beta, 0.0));
    }
    //double scaling_factor = 1.0/static_cast<double>(N_M);
#ifndef NDEBUG
    std::cout << "ZCW " << m << " gives " << result.size() << " powder angles." << std::endl;
#endif
    return result;
  }

  std::vector<Euler<>> getZCWAngles(std::uint64_t m)
  {
    // full sphere
    constexpr double c1 = 1.0;
    constexpr double c2 = 2.0;
    constexpr double c3 = 1.0;
    return getZCWAnglesFromConstants(m, c1, c2, c3);
  }

  std::vector<Euler<>> getZCWAnglesSTEP(std::uint64_t m, PowderSphere choice)
  {
    switch(choice) {
      case PowderSphere::full:
        return getZCWAnglesFromConstants(m, 1.0, 2.0, 1.0);
        break;
      case PowderSphere::hemi:
        return getZCWAnglesFromConstants(m, -1.0, 1.0, 1.0);
        break;
      case PowderSphere::octant:
        return getZCWAnglesFromConstants(m, 2.0, 1.0, 8.0);
        break;
      default:
        {
          std::vector<Euler<>> result;
          return result;
        }
    }
    std::vector<Euler<>> result;
    return result;
  }
} // namespace dnpsoup

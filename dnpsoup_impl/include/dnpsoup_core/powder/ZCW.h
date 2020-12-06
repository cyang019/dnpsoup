#ifndef DNPSOUP_ZCW_H
#define DNPSOUP_ZCW_H

#include "dnpsoup_core/powder/fibonacci.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include <vector>
#include <utility>
#include <cmath>


namespace dnpsoup {
  enum class PowderSphere : int
  {
    full=1,
    hemi=2,
    octant=3
  };

  std::vector<Euler<>> getZCWAngles(std::uint64_t m);

  std::vector<Euler<>> getZCWAnglesFromConstants(std::uint64_t m, double c1, double c2, double c3);

  std::vector<Euler<>> getZCWAnglesSTEP(std::uint64_t m, PowderSphere choice);
  std::vector<Euler<>> getZCWAnglesSTEP(std::uint64_t m, std::uint64_t cnt_gamma, PowderSphere choice);
} // namespace dnpsoup

#endif

#ifndef DNPSOUP_ZCW_H
#define DNPSOUP_ZCW_H

#include "dnpsoup_core/powder/fibonacci.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include <vector>
#include <utility>
#include <cmath>


namespace dnpsoup {
  std::vector<Euler<>> getZCWAngles(std::uint64_t m);
} // namespace dnpsoup

#endif

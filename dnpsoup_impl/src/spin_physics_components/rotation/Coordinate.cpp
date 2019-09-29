#include "dnpsoup_core/spin_physics_components/rotation/Coordinate.h"
#include "dnpsoup_core/common.h"
#include <cmath>

namespace dnpsoup {
  std::tuple<double, double> calcAnglesWithZ(const Coordinate &c)
  {
    double xy_projection = std::sqrt(c.x * c.x + c.y * c.y);
    double theta = dnpsoup::atan(xy_projection, c.z); /// angle with z
    double phi = dnpsoup::atan(c.y, c.x); /// angle with x

    return std::make_tuple(phi, theta);
  }
} // namespace dnpsoup

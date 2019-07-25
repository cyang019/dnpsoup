#include "dnpsoup_core/spin_physics_components/rotation/Coordinate.h"
#include "dnpsoup_core/common.h"
#include <cmath>

namespace dnpsoup {
  std::tuple<double, double> calcAnglesWithZ(const Coordinate &c)
  {
    double theta = 0; /// angle with z
    double phi = 0; /// angle with x
    if(std::abs(c.z) < eps){
      theta = 0.0;
    } else {
      double xy_projection = std::sqrt(c.x * c.x + c.y * c.y);
      theta = dnpsoup::atan(xy_projection, c.z);
    }

    if (std::abs(c.x) < eps){
      if(c.y > -eps) phi = 0.5 * pi;
      else phi = 1.5 * pi;
    } else {
      if(c.x > 0 && c.y >= 0){
        // 0 ~ pi/2
        phi = dnpsoup::atan(c.y, c.x);
      } else if(c.x < 0 && c.y >= 0){
        // pi/2 ~ pi
        phi = dnpsoup::atan(c.y, c.x);
      } else if(c.x < 0 && c.y < 0){
        // pi ~ 3/2pi
        phi = dnpsoup::atan(c.y, c.x) + pi;
      } else {
        // 3/2pi ~ 2pi
        phi = dnpsoup::atan(c.y, c.x) - pi;
      }
    }

    return std::make_tuple(phi, theta);
  }
} // namespace dnpsoup

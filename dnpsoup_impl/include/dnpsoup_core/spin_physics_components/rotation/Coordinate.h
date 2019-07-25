#ifndef DNPSOUP_COORDINATE_H
#define DNPSOUP_COORDINATE_H

#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "tuple"
#include <cmath>

namespace dnpsoup {
  struct Coordinate {
    Coordinate(double v1, double v2, double v3)
      : x(v1), y(v2), z(v3) {}
    Coordinate() : x(0.0), y(0.0), z(0.0) {}
    Coordinate(const Coordinate &) = default;
    Coordinate(Coordinate &&) noexcept = default;
    Coordinate& operator=(const Coordinate &) = default;
    Coordinate& operator=(Coordinate &&) noexcept = default;
    ~Coordinate() {}

    Coordinate& operator+=(const Coordinate &rhs)
    {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;
      return *this;
    }

    Coordinate& operator-=(const Coordinate &rhs)
    {
      x -= rhs.x;
      y -= rhs.y;
      z -= rhs.z;
      return *this;
    }

    double x;
    double y;
    double z;
  };

  inline Coordinate operator+(const Coordinate &c1, const Coordinate &c2)
  {
    Coordinate c = c1;
    c += c2;
    return c;
  }

  inline Coordinate operator-(const Coordinate &c1, const Coordinate &c2)
  {
    Coordinate c = c1;
    c -= c2;
    return c;
  }

  inline bool operator==(const Coordinate &c1, const Coordinate &c2)
  {
    if(std::abs(c1.x - c2.x) > eps) return false;
    if(std::abs(c1.y - c2.y) > eps) return false;
    if(std::abs(c1.z - c2.z) > eps) return false;
    return true;
  }

  std::tuple<double, double> calcAnglesWithZ(const Coordinate &c);  

} // namespace dnpsoup


#endif

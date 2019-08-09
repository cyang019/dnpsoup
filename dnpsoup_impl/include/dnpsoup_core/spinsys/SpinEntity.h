#ifndef DNPSOUP_SPINENTITY_H
#define DNPSOUP_SPINENTITY_H

#include "dnpsoup_core/spin_physics_components/rotation/Coordinate.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spinsys/SpinId.h"

namespace dnpsoup {
  class SpinEntity {
  public:
    SpinEntity();
    SpinEntity(const SpinType &t, const Coordinate &c)
      : m_type(t), m_position(c), m_t1(inf), m_t2(inf)
    {}

    SpinEntity(const SpinType &t, double x, double y, double z)
      : m_type(t), m_position(Coordinate(x,y,z)), m_t1(inf), m_t2(inf)
    {}

    SpinEntity(const SpinEntity &) = default;
    SpinEntity(SpinEntity &&) noexcept = default;
    SpinEntity& operator=(const SpinEntity &) = default;
    SpinEntity& operator=(SpinEntity &&) noexcept = default;
    ~SpinEntity() {}

    const Coordinate& getLocation() const { return m_position; }
    SpinEntity& setLocation(const Coordinate &c) 
    { m_position = c; return *this; }

    const SpinType& getSpinType() const { return m_type; }
    SpinEntity& setSpinType(const SpinType &t)
    { m_type = t; return *this; }

    const double &getT1() const { return m_t1; }
    SpinEntity& setT1(double t1);
    const double &getT2() const { return m_t2; }
    SpinEntity& setT2(double t2);
  private:
    SpinType m_type;
    Coordinate m_position;
    double m_t1;
    double m_t2;
  };
} // namespace dnpsoup

#endif

#include "dnpsoup_core/spin_physics_components/rotation/Coordinate.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spinsys/SpinId.h"

namespace dnpsoup {
  class SpinEntity {
  public:
    SpinEntity(const SpinType &t, const Coordinate &c)
      : m_type(t), m_position(c)
    {}

    SpinEntity(const SpinType &t, double x, double y, double z)
      : m_type(t), m_position(Coordinate(x,y,z))
    {}

    SpinEntity(const SpinEntity &) = default;
    SpinEntity(SpinEntity &&) noexcept = default;
    SpinEntity& operator=(const SpinEntity &) = default;
    SpinEntity& operator=(SpinEntity &&) noexcept = default;
    ~SpinEntity() {}

    Coordinate getLocation() const { return m_position; }
    SpinEntity& setLocation(const Coordinate &c) 
    { m_position = c; return *this; }

    SpinType getSpinType() const { return m_type; }
    SpinEntity& setSpinType(const SpinType &t)
    { m_type = t; return *this; }
  private:
    SpinType m_type;
    Coordinate m_position;
  };
} // namespace dnpsoup

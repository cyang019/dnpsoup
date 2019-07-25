#include "dnpsoup_core/spin_physics_components/rotation/Coordinate.h"
#include "dnpsoup_core/sph_physics_components/spin.h"
#include "dnpsoup_core/spinsys/SpinId.h"

namespace dnpsoup {
  class SpinEntity {
  public:
    SpinEntity(const SpinType &, const Coordinate &);
    SpinEntity(const SpinType &, double x, double y, double z);
  private:
    Coordinate m_position;
    SpinType m_type;
    double m_gamma;
  }
} // namespace dnpsoup

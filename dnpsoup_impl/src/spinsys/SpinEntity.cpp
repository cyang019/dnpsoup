#include "dnpsoup_core/spinsys/SpinEntity.h"


namespace dnpsoup {
  SpinEntity::SpinEntity()
    : m_type(SpinType::Null), m_position(Coordinate()),
    m_t1(inf), m_t2(inf)
  {}


  SpinEntity& SpinEntity::setT1(double t1)
  {
    m_t1 = t1;
    return *this;
  }

  SpinEntity& SpinEntity::setT2(double t2)
  {
    m_t2 = t2;
    return *this;
  }
} // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/DipolarInteraction.h"


namespace dnpsoup {
  template<DipoleType T>
  DipolarInteraction::DipolarInteraction(size_t n1, size_t n2)
    : m_n1(n1), m_n2(n2), m_nbefore(0), m_nbetween(0), m_nafter(0)
  {
  }
} // namespace dnpsoup

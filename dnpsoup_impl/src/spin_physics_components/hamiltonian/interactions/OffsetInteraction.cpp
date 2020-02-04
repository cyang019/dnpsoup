#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/OffsetInteraction.h"


namespace dnpsoup {

  OffsetInteraction::OffsetInteraction([[maybe_unused]] double gamma, size_t n)
    : InteractionInterface(), m_n(n), m_nbefore(0), m_nafter(0)
  { 
    m_z = spin<Z>(m_n);
  }

  OffsetInteraction::OffsetInteraction([[maybe_unused]] double gamma, size_t n, size_t nbefore, size_t nafter)
    : InteractionInterface(), m_n(n), m_nbefore(nbefore), m_nafter(nafter)
  { 
    m_z = kron(kron(identity<cxdbl>(nbefore), spin<Z>(n)), identity<cxdbl>(nafter));
  }

  // internally convert angles to active form for calculation
  MatrixCxDbl OffsetInteraction::genMatrix(
      const Property &p_freq,
      [[maybe_unused]] const Euler<ActiveRotation> &e) const
  {
    const double freq0 = p_freq.get(ValueName::offset);
    return freq0 * m_z;
  }

  MatrixCxDbl OffsetInteraction::genMatrix(
      const Property &p_freq,
      [[maybe_unused]] const Euler<PassiveRotation> &e) const
  {
    const double freq0 = p_freq.get(ValueName::offset);
    return freq0 * m_z;
  }

  size_t OffsetInteraction::dimension() const
  { if (m_n == 0 && m_nbefore == 0 && m_nafter == 0) 
      return 0;
    size_t dim_before = m_nbefore > 0 ? m_nbefore : 1;
    size_t dim_after = m_nafter > 0 ? m_nafter : 1;
    size_t m_mat = m_n > 0 ? m_n : 1;
    return dim_before * m_mat * dim_after;
  }
} // namespace dnpsoup


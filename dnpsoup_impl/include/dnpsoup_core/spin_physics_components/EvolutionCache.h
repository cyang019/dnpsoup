#ifndef DNPSOUP_EVOLUTIONCACHE_H
#define DNPSOUP_EVOLUTIONCACHE_H

#include "dnpsoup_core/common.h"
#include <cstdint>
#include <map>


namespace dnpsoup {
  /// to store intermediate rho_inf_eq, exp(-L dt)
  class EvolutionCache {
  public:
  private:
    std::uint64_t m_n_in_rotor_period;
    std::uint64_t m_capacity;
  };
} // namespace dnpsoup


#endif

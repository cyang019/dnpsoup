#ifndef DNPSOUP_EVOLUTIONCACHESTATIC_H
#define DNPSOUP_EVOLUTIONCACHESTATIC_H

#include "configure_dnpsoup.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/pulseseq/seq_common.h"
#include "dnpsoup_core/common.h"
#include <cstddef>
#include <cstdint>
#include <vector>
#include <map>
#include <utility>


namespace dnpsoup {
  /// rho_dt = exp(-L dt) * (rho - rho_inf_eq) + rho_inf_eq
  /// exp(-L dt) --> scaling_factor

  /// to store intermediate L, rho_inf_eq
  class EvolutionCacheStatic {
  public:
    EvolutionCacheStatic(std::uint64_t capacity);

    std::uint64_t capacity() const { return m_capacity; }

    const std::pair<MatrixCxDbl, MatrixCxDbl>& getCache(
        const pulseseq::Component &comp) const;
    /// better performance
    const std::pair<MatrixCxDbl, MatrixCxDbl>& getCache(
        int key) const;

    // L, rho_inf_eq
    EvolutionCacheStatic& saveCache(
        const pulseseq::Component &, 
        std::pair<MatrixCxDbl, MatrixCxDbl> &&);
    int getCacheIdentity(const pulseseq::Component &) const;
  private:
    std::uint64_t m_capacity;
    int m_key;
    std::vector<pulseseq::Component> m_cache_identities;
    /// L, rho_inf_eq
    std::vector<std::pair<MatrixCxDbl, MatrixCxDbl>> m_cache;
  };
} // namespace dnpsoup


#endif


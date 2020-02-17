#ifndef DNPSOUP_EVOLUTIONCACHE_H
#define DNPSOUP_EVOLUTIONCACHE_H

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
  struct EvolutionCacheElement {
    EvolutionCacheElement(
        const MatrixCxDbl &factor, const MatrixCxDbl &rho);
    EvolutionCacheElement(
        const MatrixCxDbl &factor, const MatrixCxDbl &rho,
        const MatrixCxDbl &rotate_mat_super, 
        const MatrixCxDbl &rotate_mat_super_inv);
    EvolutionCacheElement(
        MatrixCxDbl &&factor, const MatrixCxDbl &&rho,
        MatrixCxDbl &&rotate_mat_super,
        MatrixCxDbl &&rotate_mat_super_inv);
    EvolutionCacheElement();
    EvolutionCacheElement(const EvolutionCacheElement &);
    EvolutionCacheElement& operator=(const EvolutionCacheElement &);
    EvolutionCacheElement(EvolutionCacheElement &&) noexcept;
    EvolutionCacheElement& operator=(EvolutionCacheElement &&) noexcept;
    ~EvolutionCacheElement();

    MatrixCxDbl scaling_factor;
    MatrixCxDbl rho_inf_eq;
    MatrixCxDbl rotate_mat_super;
    MatrixCxDbl rotate_mat_super_inv;
  };

  /// to store intermediate rho_inf_eq, exp(-L dt)
  class EvolutionCache {
  public:
    EvolutionCache(std::uint64_t cnt, std::uint64_t capacity);
    EvolutionCache();

    std::uint64_t capacity() const { return m_capacity; }
    std::uint64_t cntPerRotorPeriod() const { return m_n_in_rotor_period; }

    EvolutionCacheElement getCache(
        const pulseseq::Component &comp,
        std::uint64_t idx) const;
    /// better performance
    EvolutionCacheElement getCache(
        int key, std::uint64_t idx) const;
    EvolutionCache& saveCache(const pulseseq::Component &, EvolutionCacheElement &&);
    int getCacheIdentity(const pulseseq::Component &) const;
    std::size_t getLength(int key) const;
  private:
    std::uint64_t m_n_in_rotor_period;
    std::uint64_t m_capacity;
    std::vector<pulseseq::Component> m_cache_identities;
    std::vector<std::vector<EvolutionCacheElement>> m_cache;
    int m_key;
  };
} // namespace dnpsoup


#endif

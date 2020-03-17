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

  /// to store intermediate exp(-Ldt), rho_inf_eq
  class EvolutionCacheStatic {
  public:
    EvolutionCacheStatic(std::uint64_t capacity);

    std::uint64_t capacity() const { return m_capacity; }

    const std::pair<MatrixCxDbl, MatrixCxDbl>& getCache(
        const pulseseq::Component &comp, std::uint64_t cnt) const;
    /// better performance
    const std::pair<MatrixCxDbl, MatrixCxDbl>& getCache(
        int key, std::uint64_t cnt) const;

    // exp(-Ldt), rho_inf_eq
    EvolutionCacheStatic& saveCache(
        const pulseseq::Component &, 
        const MatrixCxDbl &super_op,
        const MatrixCxDbl &rho_eq_super,
        std::uint64_t cnt,
        double dt);
    EvolutionCacheStatic& saveCache(
        const pulseseq::Component &, 
        MatrixCxDbl &&super_op,
        MatrixCxDbl &&rho_eq_super,
        std::uint64_t cnt,
        double dt);
    EvolutionCacheStatic& saveCache(
        const pulseseq::Component &comp,
        const MatrixCxDbl &super_op,
        const MatrixCxDbl &rho_eq_super,
        const MatrixCxDbl &scaling_factor,
        std::uint64_t cnt
        );
    EvolutionCacheStatic& saveCache(
        const pulseseq::Component &comp,
        MatrixCxDbl &&super_op,
        MatrixCxDbl &&rho_eq_super,
        MatrixCxDbl &&scaling_factor,
        std::uint64_t cnt
        );
    int getCacheIdentity(const pulseseq::Component &) const;

    /// @param comp: component
    /// @param comp_size: component size
    std::pair<int, bool>
      getCacheIdentity(const pulseseq::Component &comp, std::uint64_t comp_size) const;
  private:
    std::uint64_t m_capacity;
    int m_key;
    std::vector<pulseseq::Component> m_cache_identities;
    /// exp(-Ldt), rho_inf_eq
    std::vector<std::map<std::uint64_t, std::pair<MatrixCxDbl, MatrixCxDbl>>> m_cache;
    /// L
    std::vector<MatrixCxDbl> m_super_op_cache;
  };
} // namespace dnpsoup


#endif


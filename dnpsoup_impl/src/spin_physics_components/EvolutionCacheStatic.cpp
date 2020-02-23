#include "dnpsoup_core/spin_physics_components/EvolutionCacheStatic.h"
#include "dnpsoup_core/spin_physics_components/evolve.h"
#include "dnpsoup_core/pulseseq/seq_common.h"
#include "dnpsoup_core/errors.h"
#include <string>
#include <cstdint>

using namespace std;

namespace dnpsoup {
  // ==========================
  // class EvolutionCacheStatic
  // ==========================
  EvolutionCacheStatic::EvolutionCacheStatic(std::uint64_t capacity) 
      : m_capacity(capacity),
      m_key(0), 
      m_cache_identities(capacity, pulseseq::Component()),
      m_cache(capacity, map<uint64_t, pair<MatrixCxDbl, MatrixCxDbl>>()),
      m_super_op_cache(capacity, MatrixCxDbl())
  {}

  const std::pair<MatrixCxDbl, MatrixCxDbl>& 
    EvolutionCacheStatic::getCache(
      const pulseseq::Component &comp, uint64_t cnt) const
  {
    int key_idx = getCacheIdentity(comp);
#ifndef NDEBUG
    if(key_idx < 0) {
      throw CacheNotFoundError("Cache of the input component is missing.");
    }

    if(m_cache[key_idx].find(cnt) == m_cache[key_idx].end()){
      throw CacheNotFoundError(
          "Cache of the input component with " 
          + to_string(cnt) + " is missing.");
    }
#endif

    return m_cache[key_idx].at(cnt);
  }

  const std::pair<MatrixCxDbl, MatrixCxDbl>& 
    EvolutionCacheStatic::getCache(int key, uint64_t cnt) const
  {
#ifndef NDEBUG
    if (key < 0 || (std::uint64_t)key >= m_capacity){
      std::string err_msg = "Cache key not in range: "
        + to_string(key) + " not in [0, " + to_string(m_capacity) + "].";
      throw CacheError(err_msg);
    }

    if(m_cache[key].find(cnt) == m_cache[key].end()){
      throw CacheNotFoundError(
          "Cache of the input component with " 
          + to_string(cnt) + " is missing.");
    }
#endif

    return m_cache[key].at(cnt);
  }

  EvolutionCacheStatic& EvolutionCacheStatic::saveCache(
      const pulseseq::Component &comp, 
      const MatrixCxDbl &super_op,
      const MatrixCxDbl &rho_eq_super,
      std::uint64_t cnt,
      double dt)
  {
    int cache_idx = getCacheIdentity(comp);
    if(cache_idx < 0){
      m_cache_identities[m_key] = comp;
      m_super_op_cache[m_key] = super_op;
      cache_idx = m_key;
      ++m_key;
      m_key %= m_capacity;
    }

    const auto scaling_factor = calcExpEvolve(super_op, dt, cnt);
    m_cache[cache_idx].try_emplace(cnt, make_pair(scaling_factor, rho_eq_super));
    return *this;
  }

  EvolutionCacheStatic& EvolutionCacheStatic::saveCache(
      const pulseseq::Component &comp,
      const MatrixCxDbl &super_op,
      const MatrixCxDbl &rho_eq_super,
      const MatrixCxDbl &scaling_factor,
      std::uint64_t cnt
      )
  {
    int cache_idx = getCacheIdentity(comp);
    if(cache_idx < 0){
      m_cache_identities[m_key] = comp;
      m_super_op_cache[m_key] = super_op;
      cache_idx = m_key;
      ++m_key;
      m_key %= m_capacity;
    }
    m_cache[cache_idx].try_emplace(
        cnt, make_pair(scaling_factor, rho_eq_super));
    return *this;
  }
  
  int EvolutionCacheStatic::getCacheIdentity(const pulseseq::Component &comp) const
  {
    int idx = -1;
    for(int i = 0; i < (int)m_cache_identities.size(); ++i){
      if(sameValue(comp, m_cache_identities[i], eps)){
        return i;
      }
    }
    return idx;
  }

  std::pair<int, bool>
    EvolutionCacheStatic::getCacheIdentity(
        const pulseseq::Component &comp, std::uint64_t cnt) const
  {
    int idx = getCacheIdentity(comp);
    if(idx < 0) return make_pair(idx, false);
    bool found = (m_cache[idx].find(cnt) != m_cache[idx].end());
    return make_pair(idx, found);
  }
}   // namespace dnpsoup



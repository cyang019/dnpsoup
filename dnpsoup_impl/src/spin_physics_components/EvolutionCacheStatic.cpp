#include "dnpsoup_core/spin_physics_components/EvolutionCacheStatic.h"
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
      m_key(0)
  {
    m_cache_identities = std::vector<pulseseq::Component>(
        m_capacity, pulseseq::Component());
    m_cache = std::vector<std::pair<MatrixCxDbl, MatrixCxDbl>>(
        capacity, std::make_pair(MatrixCxDbl(), MatrixCxDbl()));
  }

  const std::pair<MatrixCxDbl, MatrixCxDbl>& 
    EvolutionCacheStatic::getCache(
      const pulseseq::Component &comp) const
  {
    int key_idx = getCacheIdentity(comp);
    if(key_idx < 0) {
      throw CacheNotFoundError("Cache of the input component is missing.");
    }

    return m_cache[key_idx];
  }

  const std::pair<MatrixCxDbl, MatrixCxDbl>& 
    EvolutionCacheStatic::getCache(int key) const
  {
    if (key < 0 || (std::uint64_t)key >= m_capacity){
      std::string err_msg = "Cache key not in range: "
        + to_string(key) + " not in [0, " + to_string(m_capacity) + "].";
      throw CacheError(err_msg);
    }

    return m_cache[key];
  }

  EvolutionCacheStatic& EvolutionCacheStatic::saveCache(
      const pulseseq::Component &comp, 
      std::pair<MatrixCxDbl, MatrixCxDbl> &&elem)
  {
    int cache_idx = getCacheIdentity(comp);
    if(cache_idx < 0){
      m_cache_identities[m_key] = comp;
      m_cache[m_key] = std::move(elem);
      ++m_key;
      m_key %= m_capacity;
    }
    else {
      m_cache[cache_idx] = std::move(elem);      
    }
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
}   // namespace dnpsoup



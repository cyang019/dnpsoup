#include "dnpsoup_core/spin_physics_components/EvolutionCache.h"
#include "dnpsoup_core/pulseseq/seq_common.h"
#include "dnpsoup_core/errors.h"
#include <string>
#include <cstdint>

using namespace std;

namespace dnpsoup {
  EvolutionCacheElement::EvolutionCacheElement(
      const MatrixCxDbl &factor, const MatrixCxDbl &rho)
    : scaling_factor(factor), rho_inf_eq(rho)
  {}

  EvolutionCacheElement::EvolutionCacheElement(
      const MatrixCxDbl &factor, const MatrixCxDbl &rho,
      const MatrixCxDbl &rotate_mat_super,
      const MatrixCxDbl &rotate_mat_super_inv)
    : scaling_factor(factor), rho_inf_eq(rho), 
    rotate_mat_super(rotate_mat_super), 
    rotate_mat_super_inv(rotate_mat_super_inv)
  {}

  EvolutionCacheElement::EvolutionCacheElement(
      MatrixCxDbl &&factor, const MatrixCxDbl &&rho,
      MatrixCxDbl &&rotate_mat_super,
      MatrixCxDbl &&rotate_mat_super_inv)
    : scaling_factor(std::move(factor)), rho_inf_eq(std::move(rho)),
    rotate_mat_super(std::move(rotate_mat_super)),
    rotate_mat_super_inv(std::move(rotate_mat_super_inv))
  {}

  EvolutionCacheElement::EvolutionCacheElement()
    : scaling_factor(), rho_inf_eq(), rotate_mat_super(), rotate_mat_super_inv()
  {}

  EvolutionCacheElement::EvolutionCacheElement(const EvolutionCacheElement &rhs)
    : scaling_factor(rhs.scaling_factor), rho_inf_eq(rhs.rho_inf_eq),
    rotate_mat_super(rhs.rotate_mat_super),
    rotate_mat_super_inv(rhs.rotate_mat_super_inv)
  {}

  EvolutionCacheElement& EvolutionCacheElement::operator=(const EvolutionCacheElement &rhs)
  {
    scaling_factor = rhs.scaling_factor;
    rho_inf_eq = rhs.rho_inf_eq;
    rotate_mat_super = rhs.rotate_mat_super;
    rotate_mat_super_inv = rhs.rotate_mat_super_inv;
    return *this;
  }

  EvolutionCacheElement::EvolutionCacheElement(EvolutionCacheElement &&rhs) noexcept
    : scaling_factor(std::move(rhs.scaling_factor)),
    rho_inf_eq(std::move(rhs.rho_inf_eq)),
    rotate_mat_super(std::move(rhs.rotate_mat_super)),
    rotate_mat_super_inv(std::move(rhs.rotate_mat_super_inv))
  {}

  EvolutionCacheElement& EvolutionCacheElement::operator=(
      EvolutionCacheElement &&rhs) noexcept
  {
    scaling_factor = std::move(rhs.scaling_factor);
    rho_inf_eq = std::move(rhs.rho_inf_eq);
    rotate_mat_super = std::move(rhs.rotate_mat_super);
    rotate_mat_super_inv = std::move(rhs.rotate_mat_super_inv);
    return *this;
  }

  EvolutionCacheElement::~EvolutionCacheElement()
  {}

  // ====================
  // class EvolutionCache
  // ====================
  EvolutionCache::EvolutionCache(std::uint64_t cnt, 
      std::uint64_t capacity)
      : m_n_in_rotor_period(cnt), m_capacity(capacity),
      m_key(0)
  {
    const pulseseq::Component default_comp = {{SpinType::Null, pulseseq::EMRadiation()}};
    m_cache_identities = std::vector<pulseseq::Component>(
        m_capacity, default_comp);
    m_cache = std::vector<std::vector<EvolutionCacheElement>>(
        capacity, std::vector<EvolutionCacheElement>());
  }

  EvolutionCache::EvolutionCache()
    : m_n_in_rotor_period(1u), m_capacity(1u), m_key(0u)
  {
    const pulseseq::Component default_comp = {{SpinType::Null, pulseseq::EMRadiation()}};
    m_cache_identities = std::vector<pulseseq::Component>(
        m_capacity, default_comp);
    m_cache = std::vector<std::vector<EvolutionCacheElement>>(
        m_capacity, std::vector<EvolutionCacheElement>());
  }

  EvolutionCacheElement EvolutionCache::getCache(
      const pulseseq::Component &comp,
      std::uint64_t idx) const
  {
    int key_idx = getCacheIdentity(comp);
    if(key_idx < 0) {
      throw CacheNotFoundError("Cache of the input component is missing.");
    }

    idx %= m_n_in_rotor_period;
    return m_cache[key_idx][idx];
  }

  EvolutionCacheElement EvolutionCache::getCache(
      int key, std::uint64_t idx) const
  {
    if (key < 0 || (std::uint64_t)key >= m_capacity){
      std::string err_msg = "Cache key not in range: "
        + to_string(key) + " not in [0, " + to_string(m_capacity) + "].";
      throw CacheError(err_msg);
    }

    return m_cache[key][idx];
  }

  EvolutionCache& EvolutionCache::saveCache(
      const pulseseq::Component &comp, EvolutionCacheElement &&elem)
  {
    int cache_idx = getCacheIdentity(comp);
    if(cache_idx < 0){
      m_cache_identities[m_key] = comp;
      m_cache[m_key] = std::vector<EvolutionCacheElement>();
      m_cache[m_key].push_back(std::move(elem));
      ++m_key;
      m_key %= m_capacity;
    }
    else {
      if(m_cache[cache_idx].size() >= m_n_in_rotor_period){
        throw CacheError("Cache cannot save over the length limit.");
      }

      m_cache[cache_idx].push_back(std::move(elem));      
    }
    return *this;
  }
  
  int EvolutionCache::getCacheIdentity(const pulseseq::Component &comp) const
  {
    int idx = -1;
    for(int i = 0; i < (int)m_cache_identities.size(); ++i){
      if(sameValue(comp, m_cache_identities[i], eps)){
        return i;
      }
    }
    return idx;
  }
  
  std::size_t EvolutionCache::getLength(int key) const
  {
    if(key < 0 || (uint64_t)key >= m_capacity){
      return 0;
    }
    return m_cache[key].size();
  }
}   // namespace dnpsoup



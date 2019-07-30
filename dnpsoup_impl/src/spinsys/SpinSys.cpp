#include "dnpsoup_core/spinsys/SpinSys.h"
#include "dnpsoup_core/errors.h"
#include <string>
#include <sstream>
#include <iterator>


namespace dnpsoup {
  std::size_t calcDimBeforeId(
      const std::map<SpinId, SpinEntity> &spins, const SpinId &sid)
  {
    auto sid_iter = spins.find(sid);
    std::size_t res = 0;
    for(auto iter = spins.begin(); iter != sid_iter; ++iter){
      const auto t = (iter->second).getSpinType();
      std::size_t dim = getMatrixDimension(t);
      if(dim > 0){
        res = (res > 0) res * dim : dim;
      }
    }
    return res;
  }

  std::size_t calcDimAfterId(
      const std::map<SpinId, SpinEntity> &spins, const SpinId &sid)
  {
    std::size_t res = 0;
    auto sid_iter = spins.find(sid);
    std::advance(sid_iter, 1);
    for(auto iter = sid_iter; iter != spins.end(); ++iter){
      const auto t = (iter->second).getSpinType();
      std::size_t dim = getMatrixDimension(t);
      if(dim > 0){
        res = (res > 0) res * dim : dim;
      }
    }
    return res;
  }

  std::size_t calcDimBetweenIds(
      const std::map<SpinId, SpinEntity> &spins, 
      const SpinId &sid1, const SpinId &sid2)
  {
    std::size_t res = 0;
    auto iter1 = spins.find(sid1);
    auto iter2 = spins.find(sid2);
    if(sid1.get() > sid2.get()){
      auto iter_temp = iter1;
      iter1 = iter2;
      iter2 = iter_temp;
    }
    std::advance(iter1, 1);
    for(auto iter = iter1; iter != iter2; ++iter){
      const auto t = (iter->second).getSpinType();
      std::size_t dim = getMatrixDimension(t);
      if(dim > 0){
        res = (res > 0) res * dim : dim;
      }
    }
    return res;
  }

  SpinSys::SpinSys()
    : m_e(Euler<>(0.0,0.0,0.0))
  {}

  SpinSys& SpinSys::addSpin(const SpinId &id_name, const SpinEntity &s) 
  {
    if(m_spins.find(id_name) == m_spins.end()){
      m_spins[id_name] = s;
    } else {
      throw DuplicationError("SpinId already in the SpinSys.");
    }
    m_dimensions.clear();
    std::size_t dim = getMatrixDimension(s.getSpinType()); 
    if(dim > 0){
      m_ntotal = (m_ntotal > 0) ? m_ntotal * dim : dim;
    }
    return *this;
  }

  SpinSys& SpinSys::addCsa(const SpinId &sid, 
      double xx, double yy, double zz, const Euler<> &e,
      double t1, double t2)
  {
    if(m_spins.find(sid) == m_spins.end()){
      const string id_str = std::to_string(sid.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    }
    auto csa = Observable(InteractionType::Csa, sid);
    auto p = Property();
    p.set(ValueName::xx, xx);
    p.set(ValueName::yy, yy);
    p.set(ValueName::zz, zz);

    csa.setProperty(p);
    csa.setEuler(e);
    csa.setT1(t1);
    csa.setT2(t2);

    const auto oid_name = ObservableId(InteractionType::Csa, sid);
    m_observables[oid_name] = csa;

    return *this;
  }

  SpinSys& SpinSys::addDipole(const SpinId &s1, const SpinId &s2, 
      double dist, double t1, double t2)
  {
    if(m_spins.find(id1) == m_spins.end()){
      const string id_str = std::to_string(s1.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    } else if(m_spins.find(id2) == m_spins.end()) {
      const string id_str = std::to_string(s2.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    }

    auto dipole = Observable(InteractionType::Dipole, s1, s2);
    auto p = Property();
    p.set(ValueName::distance, dist);

    dipole.setProperty(p);
    dipole.setT1(t1);
    dipole.setT2(t2);

    const auto oid_name = ObservableId(InteractionType::Dipole, s1, s2);
    m_observables[oid_name] = dipole;
    return *this;
  }

  SpinSys& SpinSys::addScalar(const SpinId &s1, const SpinId &s2,
      double val, double t1, double t2)
  {
    if(m_spins.find(id1) == m_spins.end()){
      const string id_str = std::to_string(s1.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    } else if(m_spins.find(id2) == m_spins.end()) {
      const string id_str = std::to_string(s2.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    }

    auto scalar = Observable(InteractionType::Scalar, s1, s2);
    auto p = Property();
    p.set(ValueName::scalar, val);
    
    scalar.setProperty(scalar);
    scalar.setT1(t1);
    scalar.setT2(t2);

    const auto oid_name = ObservableId(InteractionType::Scalar, s1, s2);
    m_observables[oid_name] = scalar;
    return *this;
  }

  SpinSys& SpinSys::addShielding(const SpinId &sid, 
      double gxx, double gyy, double gzz, const Euler<> &e,
      double t1, double t2)
  {
    if(m_spins.find(sid) == m_spins.end()){
      const string id_str = std::to_string(sid.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    }
    auto shielding = Observable(InteractionType::Shielding, sid);
    auto p = Property();
    p.set(ValueName::xx, gxx);
    p.set(ValueName::yy, gyy);
    p.set(ValueName::zz, gzz);

    shielding.setProperty(p);
    shielding.setEuler(e);
    shielding.setT1(t1);
    shielding.setT2(t2);

    const auto oid_name = ObservableId(InteractionType::Shielding, sid);
    m_observables[oid_name] = shielding;
    return *this;
  }

  SpinSys& SpinSys::rotate(const Euler<> &e)
  {
    m_e = m_e * e;    ///< first rotate e, then m_e
    return *this;
  }

  SpinSys& SpinSys::clearObservables()
  { m_observables.clear(); return *this; }

  SpinSys& SpinSys::clear()
  {
    m_observables.clear();
    m_spins.clear();
    m_dimensions.clear();
    m_ntotal = 0;
    return *this;
  }

  std::size_t SpinSys::calcTotalDimension() const 
  {
    std::size_t total_dim = 0;
    for(const auto &spin : m_spins){
      auto t = spin.second.getSpinType();
      std::size_t dim = getMatrixDimension(t);
      if(dim > 0){
        total_dim = (total_dim > 0) ? total_dim * dim : dim;
      }
    }
    return total_dim;
  }

  std::vector<std::size_t> SpinSys::calcDimensions() const
  {
    std::vector<std::size_t> dims;
    for(const auto &s : m_spins){
      const std::size_t d = getMatrixDimension(s.second.getSpinType());
      dims.push_back((d);
    }
    return dims;
  }

} // namespace dnpsoup

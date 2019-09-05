#include "dnpsoup_core/spinsys/SpinSys.h"
#include "dnpsoup_core/errors.h"
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <limits>

using namespace std;


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
        res = (res > 0) ? res * dim : dim;
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
        res = (res > 0) ? res * dim : dim;
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
        res = (res > 0) ? res * dim : dim;
      }
    }
    return res;
  }

  SpinSys::SpinSys()
    : m_e(Euler<>(0.0,0.0,0.0)), m_ntotal(0)
  {}

  SpinSys& SpinSys::addSpin(const SpinId &id_name, const SpinEntity &s, bool t_add_dipole) 
  {
    if(m_spins.find(id_name) != m_spins.end()){
      throw DuplicationError("SpinId already in the SpinSys.");
    }
    m_spins.insert({id_name, s});
    if(m_spin_types.find(s.getSpinType()) == m_spin_types.end()){
      m_spin_types.insert({s.getSpinType(), {id_name}});
    } else{
      m_spin_types[s.getSpinType()].push_back(id_name);
    }
    std::size_t dim = getMatrixDimension(s.getSpinType()); 
    if(dim > 0){
      m_ntotal = (m_ntotal > 0) ? m_ntotal * dim : dim;
    }

    // automatically add dipole interactions.
    if(t_add_dipole){
      for(const auto &s : m_spins){
        if(s.first != id_name){
          this->addDipole(s.first, id_name);
        }
      }
    }
    return *this;
  }

  SpinSys& SpinSys::addSpin(int id_val, const SpinEntity &s, bool t_add_dipole) 
  {
    SpinId id_name(id_val);
    return addSpin(id_name, s, t_add_dipole);
  }

  SpinSys& SpinSys::addSpin(int id_val, SpinType t, double x, double y, double z,
      bool t_add_dipole)
  {
    SpinId id_name(id_val);
    auto s = SpinEntity(t, x, y, z);
    return addSpin(id_name, s, t_add_dipole);
  }

  SpinSys& SpinSys::removeSpin(const SpinId &sid)
  {
    const auto t = m_spins.at(sid).getSpinType();
    m_spins.erase(sid);
    auto vec = m_spin_types.at(t);
    m_spin_types.at(t).erase(std::remove(vec.begin(), vec.end(), sid), vec.end());

    std::vector<ObservableId> obs_to_remove;
    for(const auto &ob_pair : m_observables){
      if(ob_pair.second.hasSpinId(sid)){
        obs_to_remove.push_back(ob_pair.first);
      }
    }

    for(const auto &oid : obs_to_remove){
      m_observables.erase(oid);
    }

    m_ntotal = this->calcTotalDimension();
    return *this;
  }

  SpinSys& SpinSys::removeSpin(int id_val)
  {
    auto id_name = SpinId(id_val);
    return this->removeSpin(id_name);
  }

  SpinSys& SpinSys::addCsa(const SpinId &sid, 
      double xx, double yy, double zz, const Euler<> &e)
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

    const auto oid_name = ObservableId(InteractionType::Csa, sid);
    m_observables.insert({oid_name, csa});

    return *this;
  }

  SpinSys& SpinSys::addDipole(const SpinId &s1, const SpinId &s2) 
  {
    if(m_spins.find(s1) == m_spins.end()){
      const string id_str = std::to_string(s1.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    } else if(m_spins.find(s2) == m_spins.end()) {
      const string id_str = std::to_string(s2.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    }

    auto dipole = Observable(InteractionType::Dipole, s1, s2);
    auto p = Property();
    Coordinate c1 = m_spins[s1].getLocation();
    Coordinate c2 = m_spins[s2].getLocation();
    double dist = calcDistance(c1, c2);
    p.set(ValueName::distance, dist);

    double phi, theta;
    std::tie(phi, theta) = calcAnglesWithZ(c2 - c1);
    auto e = Euler<>(phi, theta, 0);

    dipole.setProperty(p);
    dipole.setEuler(e);

    const auto oid_name = ObservableId(InteractionType::Dipole, s1, s2);
    m_observables.insert({oid_name, dipole});
    return *this;
  }

  SpinSys& SpinSys::addScalar(const SpinId &s1, const SpinId &s2,
      double val)
  {
    if(m_spins.find(s1) == m_spins.end()){
      const string id_str = std::to_string(s1.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    } else if(m_spins.find(s2) == m_spins.end()) {
      const string id_str = std::to_string(s2.get());
      std::ostringstream oss;
      oss << "SpinId " << id_str << " not found in the SpinSys.";
      throw IndexError(oss.str());
    }

    auto scalar = Observable(InteractionType::Scalar, s1, s2);
    auto p = Property();
    p.set(ValueName::scalar, val);
    
    scalar.setProperty(p);

    const auto oid_name = ObservableId(InteractionType::Scalar, s1, s2);
    m_observables.insert({oid_name, scalar});
    return *this;
  }

  SpinSys& SpinSys::addShielding(const SpinId &sid, 
      double gxx, double gyy, double gzz, const Euler<> &e)
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

    const auto oid_name = ObservableId(InteractionType::Shielding, sid);
    m_observables.insert({oid_name, shielding});
    return *this;
  }

  SpinSys& SpinSys::irradiateOn(const SpinType &t)
  {
    const auto irradiated_ids = m_spin_types.at(t);
    auto irradiation = Observable(InteractionType::EMR, irradiated_ids);
    auto p = Property();
    p.set(ValueName::freq, 0.0);
    p.set(ValueName::phase, 0.0);
    p.set(ValueName::offset, 0.0);
    irradiation.setProperty(p);
    auto oid_name = ObservableId(InteractionType::EMR, t);
    m_observables.insert({oid_name, irradiation});
    return *this;
  }

  SpinSys& SpinSys::acquireOn(const SpinType &t)
  {
    const auto acq_ids = m_spin_types.at(t);
    auto acq = Observable(InteractionType::Acquisition, acq_ids);
    auto p = Property();
    p.set(ValueName::freq, 0.0);
    p.set(ValueName::phase, 0.0);
    p.set(ValueName::offset, 0.0);
    acq.setProperty(p);
    const auto oid_name = ObservableId(InteractionType::Acquisition, t);
    m_observables.insert({oid_name, acq});
    return *this;
  }

  std::vector<RelaxationPacket> SpinSys::summarizeRelaxation() const
  {
    std::vector<RelaxationPacket> result;
    constexpr double dmax = std::numeric_limits<double>::max();
    for(const auto &s : m_spins){
      auto t1 = s.second.getT1();
      auto t2 = s.second.getT2();
      if(t1 > dmax && t2 > dmax) continue;  // if inf, skip
      auto nbefore = calcDimBeforeId(m_spins, s.first);
      auto nafter = calcDimAfterId(m_spins, s.first);
      result.emplace_back(s.first, s.second, nbefore, nafter);
    }
    return result;
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
      dims.push_back(d);
    }
    return dims;
  }

  std::vector<SpinType> SpinSys::getSpinTypes() const
  {
    std::vector<SpinType> types;

    for(const auto &s_pair : m_spins){
      types.push_back(s_pair.second.getSpinType());
    }
    return types;
  }
} // namespace dnpsoup

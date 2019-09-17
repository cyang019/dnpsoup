#include "dnpsoup_core/spinsys/SpinSys.h"
#include "dnpsoup_core/errors.h"
#include "json.hpp"
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

  SpinSys& SpinSys::addSpin(const SpinId &id_name, const SpinEntity &s, bool t_auto_add) 
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
    if(t_auto_add){
      for(const auto &s : m_spins){
        if(s.first != id_name){
          this->addDipole(s.first, id_name);
        }
      }

      auto t = s.getSpinType();
      switch(t){
        case SpinType::e:
          this->addShielding(id_name, 0.0, 0.0, 0.0, Euler<>(0.0, 0.0, 0.0));
          break;
        case SpinType::Null:
          break;
        default:
          this->addCsa(id_name, 0.0, 0.0, 0.0, Euler<>(0.0, 0.0, 0.0));
          break;
      }
    }
    return *this;
  }

  SpinSys& SpinSys::addSpin(int id_val, const SpinEntity &s, bool t_auto_add) 
  {
    SpinId id_name(id_val);
    return addSpin(id_name, s, t_auto_add);
  }

  SpinSys& SpinSys::addSpin(int id_val, SpinType t, double x, double y, double z,
      bool t_auto_add)
  {
    SpinId id_name(id_val);
    auto s = SpinEntity(t, x, y, z);
    return addSpin(id_name, s, t_auto_add);
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

  SpinSys& SpinSys::removeObservable(const ObservableId &oid)
  {
    m_observables.erase(oid);
    return *this;
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
    p.set(ValueName::b0, 0.0);

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
    p.set(ValueName::b0, 0.0);
    p.set(ValueName::offset, 0.0);

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

  MatrixCxDbl SpinSys::acquireOn(const SpinType &t) const
  {
    auto acq = AcquisitionInteraction<>(this->getSpinTypes(), t);
    Property p; ///< placeholder
    Euler<> e;  ///< placeholder
    return acq.genMatrix(p,  e);
  }

  MatrixCxDbl SpinSys::acquireOn(const std::vector<SpinId> &sids) const
  {
    auto acq = AcquisitionInteraction<>(m_spins, sids);
    Property p;
    Euler<> e;
    return acq.genMatrix(p, e);
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

  std::vector<SpinId> SpinSys::getSpinIds(const SpinType &t) const
  {
    return m_spin_types.at(t);
  }

  std::ostream& operator<<(std::ostream &os, const SpinSys &spin_sys)
  {
    json j;

    json euler_js;
    auto e = spin_sys.getEuler();
    euler_js["alpha"] = e.alpha();
    euler_js["beta"] = e.beta();
    euler_js["gamma"] = e.gamma();
    j["euler"] = euler_js;

    json spins_js;
    for(const auto &[sid, spin_entity] : spin_sys.m_spins){
      json spin_entity_js;
      spin_entity_js["type"] = toString(spin_entity.getSpinType());
      auto coord = spin_entity.getLocation();
      spin_entity_js["x"] = coord.x;
      spin_entity_js["y"] = coord.y;
      spin_entity_js["z"] = coord.z;
      auto t1 = spin_entity.getT1();
      if(t1 < std::numeric_limits<double>::max()){
        spin_entity_js["T1"] = spin_entity.getT1();
      }
      auto t2 = spin_entity.getT2();
      if(t2 < std::numeric_limits<double>::max()){
        spin_entity_js["T2"] = spin_entity.getT2();
      }
      spins_js[toString(sid)] = spin_entity_js;
    }
    j["spins"] = spins_js;

    json obs_js;
    json emr_js = json::array();
    for(const auto &[oid, obs] : spin_sys.m_observables){
      json named_js;
      InteractionType t = obs.getType();
      switch(t){
        case InteractionType::Scalar:
          {
            json scalar_js;
            auto ids = obs.getSpinIds();
            scalar_js["id1"] = ids[0].get();
            scalar_js["id2"] = ids[1].get();
            scalar_js["value"] = obs.getPropertyValue(ValueName::scalar);
            named_js["name"] = "scalar";
            named_js["entries"] = scalar_js;
          }
          obs_js.push_back(named_js);
          break;
        case InteractionType::Csa:
          {
            json csa_js;
            auto ids = obs.getSpinIds();
            csa_js["id"] = ids[0].get();
            csa_js["x"] = obs.getPropertyValue(ValueName::xx);
            csa_js["y"] = obs.getPropertyValue(ValueName::yy);
            csa_js["z"] = obs.getPropertyValue(ValueName::zz);
            json csa_euler_js;
            auto csa_euler = obs.getEuler();
            csa_euler_js["alpha"] = csa_euler.alpha();
            csa_euler_js["beta"] = csa_euler.beta();
            csa_euler_js["gamma"] = csa_euler.gamma();
            csa_js["euler"] = csa_euler_js;
            named_js["name"] = "csa";
            named_js["entries"] = csa_js;
          }
          obs_js.push_back(named_js);
          break;
        case InteractionType::Shielding:
          {
            json shielding_js;
            auto ids = obs.getSpinIds();
            shielding_js["id"] = ids[0].get();
            shielding_js["x"] = obs.getPropertyValue(ValueName::xx);
            shielding_js["y"] = obs.getPropertyValue(ValueName::yy);
            shielding_js["z"] = obs.getPropertyValue(ValueName::zz);
            json shielding_euler_js;
            auto shielding_euler = obs.getEuler();
            shielding_euler_js["alpha"] = shielding_euler.alpha();
            shielding_euler_js["beta"] =  shielding_euler.beta();
            shielding_euler_js["gamma"] = shielding_euler.gamma();
            shielding_js["euler"] = shielding_euler_js;
            named_js["name"] = "shielding";
            named_js["entries"] = shielding_js;
          }
          obs_js.push_back(named_js);
          break;
        case InteractionType::Dipole:
          {
            json dipole_js;
            auto ids = obs.getSpinIds();
            dipole_js["id1"] = ids[0].get();
            dipole_js["id2"] = ids[1].get();
            named_js["name"] = "dipole";
            named_js["entries"] = dipole_js;
          }
          obs_js.push_back(named_js);
          break;
        case InteractionType::EMR:
          {
            auto ids = obs.getSpinIds();
            auto s = spin_sys.m_spins.at(ids[0]);
            auto t_str = toString(s.getSpinType());

            if(emr_js.find(t_str) == emr_js.end()){
              emr_js.push_back(t_str);
            }
          }
          break;
        default:
          break;
      } // switch
    } // for
    j["interactions"] = obs_js;
    j["irradiation"] = emr_js;

    os << j.dump(4) << std::endl;
    return os;
  }

  std::istream& operator>>(std::istream &is, SpinSys &spin_sys)
  {
    json j;
    is >> j;
    for(auto &[sid_js, spin_js] : j["spins"].items()){
      auto t = toSpinType(spin_js["type"]);
      auto x = spin_js["x"].get<double>();
      auto y = spin_js["y"].get<double>();
      auto z = spin_js["z"].get<double>();
      auto spin_entity = SpinEntity(t, x, y, z);
      if(spin_js.find("T1") != spin_js.end()){
        spin_entity.setT1(spin_js["T1"].get<double>());
      }
      if(spin_js.find("T2") != spin_js.end()){
        spin_entity.setT2(spin_js["T2"].get<double>());
      }

      SpinId sid(stoi(sid_js));
      spin_sys.addSpin(sid, spin_entity, false);
    }
    auto euler_js = j["euler"];
    double alpha = euler_js["alpha"].get<double>();
    double beta = euler_js["beta"].get<double>();
    double gamma = euler_js["gamma"].get<double>();
    spin_sys.setEuler(Euler<>(alpha, beta, gamma));
    for(auto &interaction_js : j["interactions"]){
      auto interaction_name = interaction_js["name"].get<string>();
      auto interaction = interaction_js["entries"];
      if(interaction_name == "scalar"){
        int id1 = interaction["id1"].get<int>();
        int id2 = interaction["id2"].get<int>();
        double val = interaction["value"].get<double>();
        spin_sys.addScalar(id1, id2, val);
      }
      else if(interaction_name == "dipole"
          || interaction_name == "hyperfine"
          || interaction_name == "throughspace"){
        auto id1 = interaction["id1"].get<int>();
        auto id2 = interaction["id2"].get<int>();
        spin_sys.addDipole(SpinId(id1), SpinId(id2));
      }
      else if(interaction_name == "csa"
          || interaction_name == "shielding"){
        int sid = interaction["id"].get<int>();
        const double x = interaction["x"].get<double>();
        const double y = interaction["y"].get<double>();
        const double z = interaction["z"].get<double>();
        auto euler_js = interaction["euler"];
        const auto a = euler_js["alpha"].get<double>();
        const auto b = euler_js["beta"].get<double>();
        const auto g = euler_js["gamma"].get<double>();
        const auto angle = Euler<>(a,b,g);
        if(interaction_name == "csa"){
          spin_sys.addCsa(sid, x, y, z, angle);
        }
        else{
          spin_sys.addShielding(sid, x, y, z, angle);
        }
      }
    }

    if(j.find("irradiation") != j.end()){
      for(auto stype : j["irradiation"]){
        auto t = toSpinType(stype);
        spin_sys.irradiateOn(t);
      }
    }

    return is;
  }
} // namespace dnpsoup

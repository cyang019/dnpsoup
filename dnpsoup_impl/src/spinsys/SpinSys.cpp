#include "dnpsoup_core/spinsys/SpinSys.h"
#include "dnpsoup_core/errors.h"
#include <string>
#include <sstream>


namespace dnpsoup {
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
    return *this;
  }

  SpinSys& SpinSys::addCsa(const SpinId &sid, 
      double xx, double yy, double zz, const Euler<> &e)
  {
    if(m_spins.find(id_name) == m_spins.end()){
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

    return *this;
  }
  SpinSys& SpinSys::addCsa(const SpinId &, 
      double xx, double yy, double zz, double t1, double t2,
      const Euler<> &e)
  {
    return *this;
  }

  SpinSys& SpinSys::addDipole(const SpinId&, const SpinId&, double dist)
  {
    return *this;
  }

  SpinSys& SpinSys::addDipole(const SpinId&, const SpinId&, double dist, double t1, double t2)
  {
    return *this;
  }

  SpinSys& SpinSys::addScalar(const SpinId&, cosnt SpinId&, double val)
  {
    return *this;
  }

  SpinSys& SpinSys::addScalar(const SpinId&, cosnt SpinId&, double val, double t1, double t2)
  {
    return *this;
  }

  SpinSys& SpinSys::addShielding(const SpinId&, 
      double xx, double yy, double zz, double offset, const Euler<> &e)
  {
    return *this;
  }

  SpinSys& SpinSys::addShielding(const SpinId&, 
      double xx, double yy, double zz, double offset, double t1, double t2,
      const Euler<> &e)
  {
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
    return *this;
  }
} // namespace dnpsoup

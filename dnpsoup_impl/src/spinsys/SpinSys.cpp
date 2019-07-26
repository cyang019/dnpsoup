#include "dnpsoup_core/spinsys/SpinSys.h"
#include "dnpsoup_core/errors.h"


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

  SpinSys& SpinSys::addObservable(const InteractionType &t, const SpinId &s_id)
  {
    return *this;
  }

  SpinSys& SpinSys::addObservable(const InteractionType &t, 
      const SpinId &s_id1, const SpinId &s_id2)
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

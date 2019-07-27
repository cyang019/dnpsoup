#ifndef DNPSOUP_SPINSYS_H
#define DNPSOUP_SPINSYS_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/interactions.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/SpinEntity.h"
#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/spinsys/SpinPacket.h"
#include <vector>
#include <unordered_map>
#include <tuple>
#include <unique_ptr>
#include <string>


namespace dnpsoup {
  class SpinSys {
  public:
    SpinSys();

    SpinSys& addSpin(const SpinId &, const SpinEntity &); 
    SpinSys& addObservable(const InteractionType &, const SpinId &);
    SpinSys& addObservable(const InteractionType &, const SpinId &, const SpinId &);

    SpinPacketCollection Summarize() const;

    SpinSys& setEuler(const Euler<> &e) { m_e = e; return *this; }
    Euler<> getEuler() const { return m_e; }
    SpinSys& rotate(const Euler<> &);

    SpinSys& clearObservables();
    SpinSys& clear();
  private:
    std::unordered_map<SpinId, SpinEntity, SpinIdHash> m_spins;
    std::unordered_map<ObservableId, Observable, ObservableIdHash> m_observables;
    Euler<> m_e;
  };  // class SpinSys
} // namespace dnpsoup

#endif

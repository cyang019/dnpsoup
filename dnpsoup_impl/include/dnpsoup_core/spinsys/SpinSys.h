#ifndef DNPSOUP_SPINSYS_H
#define DNPSOUP_SPINSYS_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/interactions.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/SpinEntity.h"
#include "dnpsoup_core/spinsys/Observable.h"
#include <vector>
#include <unordered_map>
#include <tuple>
#include <unique_ptr>


namespace dnpsoup {
  class SpinSys {
  public:
    SpinSys();

    std::tuple<std::unique_ptr<InteractionInterface>, Property, Euler<>> Summarize() const;
  private:
    std::unordered_map<SpinId, SpinEntity, SpinIdHash> m_spins;
    std::unordered_map<ObservableId, Observable> m_observables;
    Euler<> m_e;
  };  // class SpinSys
} // namespace dnpsoup

#endif

#ifndef DNPSOUP_SPINSYS_H
#define DNPSOUP_SPINSYS_H

#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/interactions.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/SpinEntity.h"
#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/spinsys/SpinPacket.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <unique_ptr>
#include <string>
#include <type_traits>


namespace dnpsoup {
  class SpinSys {
  public:
    SpinSys();

    SpinSys& addSpin(const SpinId &, const SpinEntity &); 
    SpinSys& addCsa(const SpinId &, 
        double xx, double yy, double zz, const Euler<> &e);
    SpinSys& addCsa(const SpinId &, 
        double xx, double yy, double zz, double t1, double t2,
        const Euler<> &e);

    SpinSys& addDipole(const SpinId&, const SpinId&, double dist);
    SpinSys& addDipole(const SpinId&, const SpinId&, double dist, double t1, double t2);

    SpinSys& addScalar(const SpinId&, cosnt SpinId&, double val);
    SpinSys& addScalar(const SpinId&, cosnt SpinId&, double val, double t1, double t2);

    SpinSys& addShielding(const SpinId&, 
        double xx, double yy, double zz, double offset, const Euler<> &e);
    SpinSys& addShielding(const SpinId&, 
        double xx, double yy, double zz, double offset, double t1, double t2,
        const Euler<> &e);

    /// @param T: either DnpExperiment or Nmrexperiment
    /// If DnpExperiment only e in rotating frame, everything else in the lab frame.
    /// If nmr experiment, everything in the rotating frame
    template<typename T>
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

#include "dnpsoup_core/spinsys/SpinSysImpl.hpp"

#endif

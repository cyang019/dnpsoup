#ifndef DNPSOUP_SPINSYS_H
#define DNPSOUP_SPINSYS_H

#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/spin_physics_components/rotation/FrameType.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/interactions.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/SpinEntity.h"
#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/spinsys/RelaxationPacket.h"
#include <vector>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <tuple>
#include <memory>   // unique_ptr
#include <string>
#include <type_traits>
#include <iterator>


namespace dnpsoup {
  std::size_t calcDimBeforeId(
      const std::map<SpinId, SpinEntity> &spins, const SpinId &sid);
  std::size_t calcDimAfterId(
      const std::map<SpinId, SpinEntity> &spins, const SpinId &sid);
  std::size_t calcDimBetweenIds(
      const std::map<SpinId, SpinEntity> &spins, 
      const SpinId &sid1, const SpinId &sid2);

  class SpinSys {
  public:
    SpinSys();

    SpinSys& addSpin(const SpinId &, const SpinEntity &, bool t_add_dipole=true); 
    SpinSys& addSpin(int, const SpinEntity &, bool t_add_dipole=true); 
    SpinSys& addSpin(int, SpinType, double x, double y, double z, bool t_add_dipole=true);
    const std::map<SpinId, SpinEntity>& getSpins() const { return m_spins; }

    SpinSys& removeSpin(const SpinId &);
    SpinSys& removeSpin(int);

    // ===============================================
    // add observables
    // ===============================================
    SpinSys& addCsa(const SpinId &, 
        double xx, double yy, double zz, const Euler<> &e);

    SpinSys& addDipole(const SpinId&, const SpinId&);

    SpinSys& addScalar(const SpinId&, const SpinId&, double val);

    SpinSys& addShielding(const SpinId&, 
        double gxx, double gyy, double gzz, const Euler<> &e);
    // ===============================================
    
    SpinSys& irradiateOn(const SpinType &);

    /// @returns matrix for I+ operator
    MatrixCxDbl acquireOn(const SpinType &);
    MatrixCxDbl acquireOn(const std::vector<SpinId> &);

    /// @param T: either DnpExperiment or Nmrexperiment
    /// If DnpExperiment only e in rotating frame, everything else in the lab frame.
    /// If nmr experiment, everything in the rotating frame
    template<typename T>
    PacketCollection summarize() const;

    /// collect the OffsetInteraction operators (needed for relaxation, emradiation)
    template<typename T>
    PacketCollection summarizeOffset() const;

    std::vector<RelaxationPacket> summarizeRelaxation() const;

    SpinSys& setEuler(const Euler<> &e) { m_e = e; return *this; }
    const Euler<>& getEuler() const { return m_e; }
    SpinSys& rotate(const Euler<> &e);

    SpinSys& clearObservables();
    SpinSys& clear();

    std::size_t calcTotalDimension() const;
    std::vector<std::size_t> calcDimensions() const;

    std::vector<SpinType> getSpinTypes() const;
    std::vector<SpinId> getSpinIds(const SpinType &) const;
  private:
    std::map<SpinId, SpinEntity> m_spins;
    std::unordered_map<SpinType, std::vector<SpinId>> m_spin_types;
    std::unordered_map<ObservableId, Observable, ObservableIdHash> m_observables;
    Euler<> m_e;
    std::size_t m_ntotal;

    /// need to use position info from SpinSys
    template<typename T>
    std::unique_ptr<InteractionInterface> genInteractionFromObservable(const Observable&) const;
  };  // class SpinSys
} // namespace dnpsoup

#include "dnpsoup_core/spinsys/SpinSysImpl.hpp"

#endif

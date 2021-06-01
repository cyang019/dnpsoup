#ifndef DNPSOUP_ACQUISITIONINTERACTION_H
#define DNPSOUP_ACQUISITIONINTERACTION_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/SpinEntity.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/rotation/FrameType.h"
#include <vector>


namespace dnpsoup {
  /// Microwave or Radio Frequency Interaction
  /// enable only if T is a frame type (i.e. RotatingFrame or LabFrame)
  template<typename T = RotatingFrame, typename E = DnpExperiment>
  class AcquisitionInteraction
  : public InteractionInterface {
  public:
    AcquisitionInteraction(const std::vector<SpinType> &spins, const SpinType &irradiated);
    AcquisitionInteraction(const std::map<SpinId, SpinEntity> &spins,
                           const std::vector<SpinId> &irradiated);
    ~AcquisitionInteraction() {}

    MatrixCxDbl genMatrix(
        [[maybe_unused]] const Property &,
        [[maybe_unused]] const Euler<ActiveRotation> &,
        [[maybe_unused]] const Euler<ActiveRotation> &e2=default_euler_a,
        [[maybe_unused]] const Euler<ActiveRotation> &e3=default_euler_a) const override;

    MatrixCxDbl genMatrix(
        [[maybe_unused]] const Property &,
        [[maybe_unused]] const Euler<PassiveRotation> &,
        [[maybe_unused]] const Euler<PassiveRotation> &e2=default_euler_p,
        [[maybe_unused]] const Euler<PassiveRotation> &e3=default_euler_p) const override;

    std::size_t dimension() const;
  private:
    std::size_t m_ntotal;
    matrix::Matrix<cxdbl> m_op;
  };  // class ChemicalShiftInteraction
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/AcquisitionInteractionImpl.hpp"

#endif




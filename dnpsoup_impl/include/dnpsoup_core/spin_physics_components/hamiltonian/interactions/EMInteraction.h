#ifndef DNPSOUP_EMINTERACTION_H
#define DNPSOUP_EMINTERACTION_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/rotation/FrameType.h"
#include <vector>


namespace dnpsoup {
  /// Microwave or Radio Frequency Interaction
  /// enable only if T is a frame type (i.e. RotatingFrame or LabFrame)
  template<typename T>
  class EMInteraction
  : public InteractionInterface {
  public:
    EMInteraction(const std::vector<SpinType> &spins, const SpinType &irradiated);
    ~EMInteraction() {}

    matrix::Matrix<cxdbl> genMatrix(
        const Property &,   // freq = 0.5 * gamma * B1 in Hz; phase in degrees. phase0 in degrees
        [[maybe_unused]] const Euler<ActiveRotation> &) const override;
    matrix::Matrix<cxdbl> genMatrix(
        const Property &,   // freq = 0.5 * gamma * B1 in Hz; phase in degrees. phase0 in degrees
        [[maybe_unused]] const Euler<PassiveRotation> &) const override;

    std::size_t dimension() const;
  private:
    std::size_t m_ntotal;

    matrix::Matrix<cxdbl> m_x;
    matrix::Matrix<cxdbl> m_y;

    matrix::Matrix<cxdbl> genMatrixInternal(const Property &) const;
  };  // class ChemicalShiftInteraction
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/EMInteractionImpl.hpp"

#endif



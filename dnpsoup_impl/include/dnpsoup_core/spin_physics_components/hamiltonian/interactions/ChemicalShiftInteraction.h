#ifndef DNPSOUP_CHEMICALSHIFTINTERACTION_H
#define DNPSOUP_CHEMICALSHIFTINTERACTION_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/rotation/FrameType.h"


namespace dnpsoup {
  /// enable only if T is a frame type (i.e. RotatingFrame or LabFrame)
  template<typename T>
  class ChemicalShiftInteraction
  : public InteractionInterface {
  public:
    ChemicalShiftInteraction(double gamma, size_t n);
    ChemicalShiftInteraction(double gamma, size_t n, size_t nbefore, size_t nafter);
    ~ChemicalShiftInteraction() {}

    MatrixCxDbl genMatrix(
        const Property &,
        const Euler<ActiveRotation> &,
        [[maybe_unused]] const Euler<ActiveRotation> &e2=default_euler_a,
        [[maybe_unused]] const Euler<ActiveRotation> &e3=default_euler_a) const override;

    MatrixCxDbl genMatrix(
        const Property &,
        const Euler<PassiveRotation> &,
        [[maybe_unused]] const Euler<PassiveRotation> &e2=default_euler_p,
        [[maybe_unused]] const Euler<PassiveRotation> &e3=default_euler_p) const override;

    size_t dimension() const;
  private:
    matrix::Matrix<cxdbl> m_x;
    matrix::Matrix<cxdbl> m_y;
    matrix::Matrix<cxdbl> m_z;

    double m_gamma;
    size_t m_n;
    size_t m_nbefore;
    size_t m_nafter;
  };  // class ChemicalShiftInteraction
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/ChemicalShiftInteractionImpl.hpp"

#endif


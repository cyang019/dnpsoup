#ifndef DNPSOUP_SCALARINTERACTION_H
#define DNPSOUP_SCALARINTERACTION_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/rotation/FrameType.h"


namespace dnpsoup {
  /// enable only if T is a frame type (i.e. RotatingFrame or LabFrame)
  class ScalarInteraction
  : public InteractionInterface {
  public:
    ScalarInteraction(double g1, double g2, size_t n1, size_t n2);
    ScalarInteraction(double g1, double g2, size_t n1, size_t n2, 
        size_t nbefore, size_t nbetween, size_t nafter);
    ~ScalarInteraction() {}

    MatrixCxDbl genMatrix(
        const Property &,
        [[maybe_unused]] const Euler<ActiveRotation> &,
        [[maybe_unused]] const Euler<ActiveRotation> &e2=default_euler_a,
        [[maybe_unused]] const Euler<ActiveRotation> &e3=default_euler_a) const override;

    MatrixCxDbl genMatrix(
        const Property &,
        [[maybe_unused]] const Euler<PassiveRotation> &,
        [[maybe_unused]] const Euler<PassiveRotation> &e2=default_euler_p,
        [[maybe_unused]] const Euler<PassiveRotation> &e3=default_euler_p) const override;


    size_t dimension() const;
  private:
    matrix::Matrix<cxdbl> m_zz;

    size_t m_n1;
    size_t m_n2;
    size_t m_nbefore;
    size_t m_nbetween;
    size_t m_nafter;

    double m_gamma1;
    double m_gamma2;
  };  // class ChemicalShiftInteraction
}   // namespace dnpsoup

#endif



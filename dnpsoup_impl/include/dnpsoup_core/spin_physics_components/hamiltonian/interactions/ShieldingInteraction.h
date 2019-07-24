#ifndef DNPSOUP_SHIELDINGINTERACTION_H
#define DNPSOUP_SHIELDINGINTERACTION_H

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
  class ShieldingInteraction
  : public InteractionInterface {
  public:
    ShieldingInteraction(double beta, size_t n);
    ShieldingInteraction(double beta, size_t n, size_t nbefore, size_t nafter);
    ~ShieldingInteraction() {}

    template<typename R>
    matrix::Matrix<cxdbl> genMatrix(
        const Property &,
        const Euler<R> &) const override;

    size_t dimension() const;
  private:
    matrix::Matrix<cxdbl> m_x;
    matrix::Matrix<cxdbl> m_y;
    matrix::Matrix<cxdbl> m_z;

    size_t m_n;
    size_t m_nbefore;
    size_t m_nafter;

    double m_beta;
  };  // class ShieldingInteraction
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/ShieldingInteractionImpl.hpp"

#endif



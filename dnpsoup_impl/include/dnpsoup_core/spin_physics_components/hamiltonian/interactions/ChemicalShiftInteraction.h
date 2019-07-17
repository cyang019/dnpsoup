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
  template<FrameType T>
  class ChemicalShiftInteraction : public InteractionInterface {
  public:
    ChemicalShiftInteraction(size_t n);
    ChemicalShiftInteraction(size_t n, size_t nbefore, size_t nafter);
    ~ChemicalShiftInteraction() {}

    // active rotation
    virtual matrix::Matrix<cxdbl> genMatrix(
        const PropertyValue *,
        const Euler &) const override;

    size_t dimension() const;
  private:
    matrix::Matrix<cxdbl> m_iz;

    size_t m_n;
    size_t m_nbefore;
    size_t m_nafter;
  };  // class Shift
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/ChemicalShiftInteractionImpl.hpp"

#endif


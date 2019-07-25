#ifndef DNPSOUP_INTERACTIONINTERFACE_H
#define DNPSOUP_INTERACTIONINTERFACE_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"

namespace dnpsoup {
  class InteractionInterface {
  public:
    InteractionInterface() {}
    virtual ~InteractionInterface() {}

    virtual MatrixCxDbl genMatrix(
        const Property &, 
        [[maybe_unused]] const Euler<ActiveRotation> &) const = 0;

    virtual MatrixCxDbl genMatrix(
        const Property &, 
        [[maybe_unused]] const Euler<PassiveRotation> &) const = 0;
  };
} // namespace dnpsoup

#endif

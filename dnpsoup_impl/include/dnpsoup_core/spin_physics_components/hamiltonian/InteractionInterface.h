#ifndef DNPSOUP_INTERACTIONINTERFACE_H
#define DNPSOUP_INTERACTIONINTERFACE_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/PropertyValue.h"

namespace dnpsoup {
  class InteractionInterface {
  public:
    InteractionInterface() {}
    virtual ~InteractionInterface() {}

    virtual MatrixCxDbl genMatrix(
        const PropertyValue *, 
        const Euler &) const = 0;
  };
} // namespace dnpsoup

#endif

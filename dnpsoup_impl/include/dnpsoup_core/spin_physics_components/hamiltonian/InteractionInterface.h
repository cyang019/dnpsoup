#ifndef DNPSOUP_INTERACTIONINTERFACE_H
#define DNPSOUP_INTERACTIONINTERFACE_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotations/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/PropertyValueInterface.h"

namespace dnpsoup {
  class InteractionInterface {
    InteractionInterface() {}
    virtual ~InteractionInterface() {}

    virtual MatrixCxDbl genMatrix(
        const PropertyValueInterface *, 
        const Euler &) const = 0;
  }
} // namespace dnpsoup

#endif

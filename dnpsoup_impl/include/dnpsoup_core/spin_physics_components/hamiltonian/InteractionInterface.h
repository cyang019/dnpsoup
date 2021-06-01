#ifndef DNPSOUP_INTERACTIONINTERFACE_H
#define DNPSOUP_INTERACTIONINTERFACE_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"

namespace dnpsoup {
  inline const Euler<ActiveRotation> default_euler_a;
  inline const Euler<PassiveRotation> default_euler_p;

  class InteractionInterface {
  public:
    InteractionInterface() {}
    virtual ~InteractionInterface() {}

    virtual MatrixCxDbl genMatrix(
        const Property &, 
        [[maybe_unused]] const Euler<ActiveRotation> &,
        [[maybe_unused]] const Euler<ActiveRotation> &e2 = default_euler_a,
        [[maybe_unused]] const Euler<ActiveRotation> &e3 = default_euler_a
        ) const = 0;

    virtual MatrixCxDbl genMatrix(
        const Property &, 
        [[maybe_unused]] const Euler<PassiveRotation> &,
        [[maybe_unused]] const Euler<PassiveRotation> &e2 = default_euler_p,
        [[maybe_unused]] const Euler<PassiveRotation> &e3 = default_euler_p
        ) const = 0;
  };
} // namespace dnpsoup

#endif

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/constants.h"

namespace dnpsoup {
  double getGyromagneticRatio(const SpinType &t)
  {
    switch(t){
      case SpinType::e:
        return beta_e;
      case SpinType::H:
        return gamma_H1;
      case SpinType::D:
        return gamma_D2;
      case SpinType::C13:
        return gamma_C13;
      case SpinType::N14:
        return gamma_N14;
      case SpinType::N15:
        return gamma_N15;
      case SpinType::O17:
        return gamma_O17;
      default: return 0.0;
    }
  }
} // namespace dnpsoup

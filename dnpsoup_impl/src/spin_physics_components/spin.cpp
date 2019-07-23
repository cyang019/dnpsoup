#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/constants.h"
#include <utility>


namespace dnpsoup {
  std::size_t SpinTypeHash::operator()(const SpinType &spin_t) const
  {
    return std::hash<int>()(static_cast<int>(spin_t));
  }

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

  size_t getMatrixDimension(const SpinType &t)
  {
    switch(t){
      case SpinType::e:
        return 2;
      case SpinType::H:
        return 2;
      case SpinType::D:
        return 3;
      case SpinType::C13:
        return 2;
      case SpinType::N14:
        return 3;
      case SpinType::N15:
        return 2;
      case SpinType::O17:
        return 10;
      default: return 2;
    }
  }
} // namespace dnpsoup

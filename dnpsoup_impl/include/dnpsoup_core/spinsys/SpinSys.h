#ifndef DNPSOUP_SPINSYS_H
#define DNPSOUP_SPINSYS_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include <vector>


namespace dnpsoup {
  class SpinSys {
  public:
    SpinSys();
  private:
    std::vector<SpinType> spins;
    std::vector<
  };  // class SpinSys
} // namespace dnpsoup

#endif

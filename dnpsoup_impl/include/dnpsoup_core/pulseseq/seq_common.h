#ifndef DNPSOUP_SEQ_COMMON_H
#define DNPSOUP_SEQ_COMMON_H

#include "dnpsoup_core/pulseseq/EMRadiation.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include <string>
#include <iostream>
#include <map>

namespace dnpsoup {
  namespace pulseseq {
    using Name = std::string;
    using Component = std::map<SpinType, EMRadiation>;

    std::istream& operator>>(std::istream &is, Component &comp);
    std::ostream& operator<<(std::ostream &os, const Component &comp);
  } // namespace pulseseq
  bool sameValue(const pulseseq::Component &comp1, const pulseseq::Component &comp2, double eps);
}
  
  
#endif

#ifndef DNPSOUP_SEQ_COMMON_H
#define DNPSOUP_SEQ_COMMON_H

#include "dnpsoup_core/pulseseq/EMRadiation.h"
#include <string>
#include <map>

namespace dnpsoup {
  namespace pulseseq {
    using Name = std::string;
    using Component = std::map<SpinType, EMRadiation>;
  }
}
  
  
#endif

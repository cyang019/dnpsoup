#ifndef DNPSOUP_PULSEINTERFACE_H
#define DNPSOUP_PULSEINTERFACE_H

#include "dnpsoup_core/pulseseq/PulseComponent.h"
#include <map>
#include <utility>
#include <cstdint>


namespace dnpsoup {
  class PulseInterface {
  public:
    PulseInterface() {}
    virtual ~PulseInterface() {}

    virtual std::pair<std::map<SpinType, PulseComponent>, std::size_t> next() = 0;
    std::uint64_t size() const = 0; 
  };  // class PulseInterface
} // namespace dnpsoup


#endif

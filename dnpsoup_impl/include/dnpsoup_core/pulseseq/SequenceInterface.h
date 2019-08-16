#ifndef DNPSOUP_SEQUENCEINTERFACE_H
#define DNPSOUP_SEQUENCEINTERFACE_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/pulseseq/PulsePacket.h"
#include <map>
#include <utility>  // pair
#include <cstdint>


namespace dnpsoup {
  enum class SequenceType : int
  {
    Pulse         = 0,
    Sequence      = 1,
    Shape         = 10,
    Chirp         = 20
  };

  class SequenceInterface {
  public:
    virtual std::pair<std::map<SpinType, PulsePacket>, std::uint64_t> next() = 0;
    virtual std::uint64_t size() const = 0;
    virtual void setSize(std::uint64_t) = 0;
    virtual SequenceType type() const = 0;   // type of interaction

    virtual std::uint64_t getIndex() const = 0;
    virtual void resetIndex() = 0;
  };
} // namespace dnpsoup

#endif

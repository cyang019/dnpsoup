#ifndef DNPSOUP_PULSESEQUENCE_H
#define DNPSOUP_PULSESEQUENCE_H

#include "dnpsoup_core/pulseseq/PulseInterface.h"
#include <memory>
#include <vector>


namespace dnpsoup {
  class PulseSequence {
  public:
    PulseSequence();

    PulseSequence& push_back(std::unique_ptr<PulseInterface>);

    virtual std::pair<std::map<SpinType, PulseComponent>, std::size_t> next();
    std::uint64_t size() const;
  private:
    std::vector<std::unique_ptr<PulseInterface>> m_seq;
    std::uint64_t m_idx;
  };
} // namespace dnpsoup

#endif

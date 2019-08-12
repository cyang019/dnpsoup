#ifndef DNPSOUP_PULSESEQ_H
#define DNPSOUP_PULSESEQ_H

#include "dnpsoup_core/pulseseq/Segment.h"


namespace dnpsoup {
  class PulseSeq {
  public:
  private:
    std::vector<Segment> m_segments;
  };
} // namespace dnpsoup

#endif

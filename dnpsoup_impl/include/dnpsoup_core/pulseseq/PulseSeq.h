#ifndef DNPSOUP_PULSESEQ_H
#define DNPSOUP_PULSESEQ_H

#include "dnpsoup_core/pulseseq/Segment.h"
#include <vector>
#include <iostream>


namespace dnpsoup {
  class PulseSeq {
    friend std::istream& operator>>(std::istream&, PulseSeq &);
    friend std::ostream& operator<<(std::ostream&, const PulseSeq &);
  public:
    PulseSeq();
    PulseSeq(std::istream&);
    PulseSeq(const PulseSeq&) = default;
    PulseSeq(PulseSeq &&) noexcept = default;
    PulseSeq& operator=(const PulseSeq &) = default;
    PulseSeq& operator=(PulseSeq &&) noexcept = default;
    ~PulseSeq() {}

    PulseSeq& addSegment(const Segment &);
    std::tuple<std::map<SpinType, PulsePacket>, std::uint64_t> next();

    const Segment& operator[](std::size_t) const;
    Segment& operator[](std::size_t);

    std::size_t size() const { return m_segments.size(); }
  private:
    std::vector<Segment> m_segments;
    std::uint64_t m_seg_index;
  };

  std::istream& operator>>(std::istream&, PulseSeq &);
  std::ostream& operator<<(std::ostream&, const PulseSeq &);
} // namespace dnpsoup

#endif

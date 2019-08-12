#ifndef DNPSOUP_SEGMENT_H
#define DNPSOUP_SEGMENT_H

#include "dnpsoup_core/pulseseq/PulseComponent.h"
#include <vector>
#include <tuple>
#include <cstdint>


namespace dnpsoup {
  class Segment {
  public:
    Segment(double inc);
    Segment(const Segment&) = default;
    Segment(Segment &&) noexcept = default;
    Segment& operator=(const Segment &) = default;
    Segment& operator=(Segment &&) noexcept = default;
    ~Segment() {}

    double getIncrement() const { return m_inc; }
    Segment& setIncrement(double inc) 
    { m_inc = inc; return *this; }

    std::uint64_t getRepetition() const { return m_repetition; }
    Segment& setRepetition(std::uint64_t rep)
    { m_repetition = rep; return *this; }

    Segment& addComponent(const PulseComponent &);
    std::tuple<std::unordered_map<SpinType, PulsePacket>, std::uint64_t> next();

    const PulseComponent& operator[](std::size_t) const;
    PulseComponent& operator[](std::size_t);

    Segment& resetIndex();

    std::uint64_t getRepetitionIndex() const { return m_rep_index; }
    std::uint64_t getComponentIndex() const { return m_comp_index; }
  private:
    std::vector<PulseComponent> m_components; 
    double m_inc;   // default increment
    std::uint64_t m_repetition;
    std::uint64_t m_comp_index;
    std::uint64_t m_rep_index;

    std::tuple<std::unordered_map<SpinType, PulsePacket>, std::uint64_t> nextInOneRep();
  };  // class Segment
} // namespace dnpsoup

#endif

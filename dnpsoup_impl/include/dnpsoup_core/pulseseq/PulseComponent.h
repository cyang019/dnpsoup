#ifndef DNPSOUP_PULSECOMPONENT_H
#define DNPSOUP_PULSECOMPONENT_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/pulseseq/PulsePacket.h"
#include <unordered_map>
#include <utility>
#include <tuple>
#include <cstdint>


namespace dnpsoup {
  class PulseComponent {
  public:
    PulseComponent(std::uint64_t cnt, double inc);
    PulseComponent(const PulseComponent&) = default;
    PulseComponent(PulseComponent &&) noexcept = default;
    PulseComponent& operator=(const PulseComponent &) = default;
    PulseComponent& operator=(PulseComponent &&) noexcept = default;
    ~PulseComponent() {}

    double getIncrement() const { return m_inc; }
    std::uint64_t getCount() const { return m_count; }

    PulseComponent& setIncrement(double);
    PulseComponent& setCount(std::uint64_t);

    const PulsePacket& getChannel(const SpinType &) const;
    PulseComponent& setChannel(const SpinType &, const PulsePacket &);
    PulseComponent& removeChannel(const SpinType &);

    std::tuple<std::unordered_map<SpinType, PulsePacket>, std::uint64_t> next();
    PulseComponent& resetIndex();
    std::uint64_t getIndex() const { return m_index; }
  private:
    std::unordered_map<SpinType, PulsePacket> m_channels;
    std::uint64_t m_count; // number of steps
    double m_inc; // size of each step
    std::uint64_t m_index;
  };
} // namespace dnpsoup

#endif

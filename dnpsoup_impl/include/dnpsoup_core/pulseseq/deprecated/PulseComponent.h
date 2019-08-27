#ifndef DNPSOUP_PULSECOMPONENT_H
#define DNPSOUP_PULSECOMPONENT_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/pulseseq/PulsePacket.h"
#include "dnpsoup_core/pulseseq/SequenceInterface.h"
#include <map>
#include <utility>
#include <tuple>
#include <cstdint>
#include <string>


namespace dnpsoup {
  class PulseComponent : public SequenceInterface 
  {
  public:
    friend std::istream& operator>>(std::istream&, PulseComponent &);
    friend std::ostream& operator<<(std::ostream&, const PulseComponent &);
    PulseComponent(const std::string &);
    PulseComponent(const std::string &, std::uint64_t cnt);
    PulseComponent(const PulseComponent&) = default;
    PulseComponent(PulseComponent &&) noexcept = default;
    PulseComponent& operator=(const PulseComponent &) = default;
    PulseComponent& operator=(PulseComponent &&) noexcept = default;
    ~PulseComponent() {}

    std::pair<std::map<SpinType, PulsePacket>, std::uint64_t> next() override;
    std::uint64_t size() const override { return m_count; }
    void setSize(std::uint64_t) override;
    SequenceType type() const override;   // type of interaction

    PulseComponent& resetIndex() override;
    std::uint64_t getIndex() const override { return m_index; }

    const PulsePacket& getChannel(const SpinType &) const;
    PulseComponent& setChannel(const SpinType &, const PulsePacket &);
    PulseComponent& removeChannel(const SpinType &);

    std::string name;
  private:
    std::map<SpinType, PulsePacket> m_channels;
    std::uint64_t m_count; // number of steps
    std::uint64_t m_index;
  };

  std::istream& operator>>(std::istream&, PulseComponent &);
  std::ostream& operator<<(std::ostream&, const PulseComponent &);
} // namespace dnpsoup

#endif

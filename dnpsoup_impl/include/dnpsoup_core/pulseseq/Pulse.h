#ifndef DNPSOUP_PULSE_H
#define DNPSOUP_PULSE_H

#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"

namespace dnpsoup {
  namespace pulseseq {
    /// @brief constant pulse
    class Pulse : public SubSequenceInterface {
    public:
      Pulse();
      Pulse(const Name &);
      Pulse(std::uint64_t, const Name &);
      Pulse(const Pulse &) = default;
      Pulse(Pulse &&) noexcept = default;
      Pulse& operator=(const Pulse &) = default;
      Pulse& operator=(Pulse &&) noexcept = default;
      virtual ~Pulse() {};

      virtual std::unique_ptr<SubSequenceInterface> copy() const override;

      virtual std::tuple<Component, std::uint64_t, std::uint64_t> next(
          std::map<Name, Component> *components,
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *m_sections
          ) override;
      virtual SequenceType type() const override;
      virtual std::vector<Name> getNames() const override;
      virtual void resetIndex(
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *m_sections) override;
    private:
      Name m_component_name;
    };
  } // namespace pulseseq
} // namespace dnpsoup

#endif

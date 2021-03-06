#ifndef DNPSOUP_DELAY_H
#define DNPSOUP_DELAY_H

#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"

namespace dnpsoup {
  namespace pulseseq {
    class Delay : public SubSequenceInterface {
    public:
      Delay();
      Delay(std::uint64_t);
      Delay(const Delay &) = default;
      Delay(Delay &&) noexcept = default;
      Delay& operator=(const Delay &) = default;
      Delay& operator=(Delay &&) noexcept = default;
      virtual ~Delay() {};

      virtual std::unique_ptr<SubSequenceInterface> copy() const override;
      virtual std::tuple<Component, std::uint64_t, std::uint64_t> next(
          std::map<Name, Component> *components,
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
          ) override;
      virtual SequenceType type() const override;
      virtual std::vector<Name> getNames() const override;
      virtual void resetIndex(
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *m_sections) override;
    };
  } // namespace pulseseq
} // namespace dnpsoup

#endif


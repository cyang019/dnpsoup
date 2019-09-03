#ifndef DNPSOUP_PULSESEQUENCE_H
#define DNPSOUP_PULSESEQUENCE_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/pulseseq/seq_common.h"
#include "dnpsoup_core/pulseseq/EMRadiation.h"
#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <iostream>


namespace dnpsoup {
  namespace pulseseq{
    class PulseSequence {
      friend std::istream& operator>>(std::istream &is, PulseSequence &);
      friend std::ostream& operator<<(std::ostream &os, const PulseSequence &);
    public:
      PulseSequence();
      PulseSequence(double inc);
      PulseSequence(const std::string &, double inc);
      PulseSequence(const std::string &);
      PulseSequence(const PulseSequence &) = delete;
      PulseSequence(PulseSequence &&) noexcept = default;
      PulseSequence& operator=(const PulseSequence &) = delete;
      PulseSequence& operator=(PulseSequence &&) noexcept = default;
      ~PulseSequence();

      const Component& getComponent(const Name &) const;
      PulseSequence& set(const Name &, const Component &);
      PulseSequence& set(const Name &, std::unique_ptr<SubSequenceInterface>);

      std::vector<Name> getNames(const Name &) const;
      std::uint64_t getIdx() const { return m_idx; }

      double getParam(const Name &seq_name, const Name &param_name) const;
      PulseSequence& setParam(const Name &seq_name, const Name &param_name, double value);

      PulseSequence& set(const std::vector<Name> &seq_names);
      // reset sequence index to 0
      PulseSequence& remove(const Name &name);

      std::pair<Component, std::uint64_t> next();
      std::size_t size() const { return m_sections_in_order.size(); }

      Name name;
    private:
      std::map<Name, Component> m_components;
      std::map<Name, std::unique_ptr<SubSequenceInterface>> m_sections;
      std::vector<Name> m_sections_in_order;
      double m_inc;
      std::uint64_t m_idx;
    };

    std::istream& operator>>(std::istream &is, PulseSequence &);
    std::ostream& operator<<(std::ostream &os, const PulseSequence &);
  } // namespace pulseseq
} // namespace dnpsoup

#endif

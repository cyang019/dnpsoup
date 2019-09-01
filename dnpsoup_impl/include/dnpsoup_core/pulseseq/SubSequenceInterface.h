#ifndef DNPSOUP_SUBSEQUENCEINTERFACE_H
#define DNPSOUP_SUBSEQUENCEINTERFACE_H

#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/pulseseq/seq_common.h"
#include "dnpsoup_core/pulseseq/EMRadiation.h"
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <cstdint>
#include <iostream>


namespace dnpsoup {
  namespace pulseseq {
    enum class SequenceType : int
    {
      PulseType         = 0,
      DelayType         = 1,
      ChirpType         = 2,
      SectionType       = 10,
      DefaultType       = 99
    };

    std::string toString(const SequenceType &);
    SequenceType toSequenceType(const std::string &);

    std::ostream& operator<<(std::ostream &os, const SequenceType &t);
    std::istream& operator>>(std::istream &is, SequenceType &t);

    class SubSequenceInterface {
    public:
      SubSequenceInterface();
      SubSequenceInterface(std::uint64_t);
      SubSequenceInterface(const SubSequenceInterface &) = default;
      SubSequenceInterface(SubSequenceInterface &&) noexcept = default;
      SubSequenceInterface& operator=(const SubSequenceInterface &) = default;
      SubSequenceInterface& operator=(SubSequenceInterface &&) noexcept = default;
      virtual ~SubSequenceInterface() {};

      virtual std::pair<Component, std::uint64_t> next(
          std::map<Name, Component> *components,
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *m_sections
          ) = 0;
      virtual SequenceType type() const = 0;   // type of interaction

      virtual std::vector<Name> getNames() const = 0;

      virtual std::uint64_t size() const { return m_sz; }
      virtual void setSize(std::uint64_t n) { m_sz = n; }
      virtual std::uint64_t getIndex() const { return m_idx; }
      virtual void resetIndex() { m_idx = 0; }

      double getParam(const Name &n) const { return m_params.at(n); }
      void setParam(const Name &, double);
      std::vector<Name> getParamNames() const;

      Name name;
    protected:
      std::uint64_t m_sz;
      std::uint64_t m_idx;
    private:
      std::map<Name, double> m_params;
    };

  } // namespace pulseseq
} // namespace dnpsoup

#endif

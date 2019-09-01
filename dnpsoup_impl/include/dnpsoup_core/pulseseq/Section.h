#ifndef DNPSOUP_SECTION_H
#define DNPSOUP_SECTION_H

#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"
#include <vector>

namespace dnpsoup {
  namespace pulseseq {
    class Section : public SubSequenceInterface {
    public:
      Section();
      Section(const std::vector<Name> &);
      Section(std::uint64_t, const std::vector<Name> &);
      Section(const Section &) = default;
      Section(Section &&) noexcept = default;
      Section& operator=(const Section &) = default;
      Section& operator=(Section &&) noexcept = default;
      virtual ~Section() {};

      /// names in components should not overlap with names in sections
      virtual std::pair<Component, std::uint64_t> next(
          std::map<Name, Component> *components,
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
          ) override;

      virtual SequenceType type() const override;
      virtual std::vector<Name> getNames() const override;
      virtual void resetIndex() override;
    private:
      /// names in sections
      std::vector<Name> m_names;
      std::uint64_t m_names_idx;
    };
  } // namespace pulseseq
} // namespace dnpsoup

#endif


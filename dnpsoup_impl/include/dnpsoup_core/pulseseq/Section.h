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
      Section(std::uint64_t, const std::vector<Name> &, bool flag_reset);
      Section(std::uint64_t, const std::vector<Name> &, 
          bool flag_reset, std::uint64_t seed);
      Section(const Section &) = default;
      Section(Section &&) noexcept = default;
      Section& operator=(const Section &) = default;
      Section& operator=(Section &&) noexcept = default;
      virtual ~Section() {};

      /// names in components should not overlap with names in sections
      virtual std::tuple<Component, std::uint64_t, std::uint64_t> next(
          std::map<Name, Component> *components,
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
          ) override;

      virtual SequenceType type() const override;
      virtual std::vector<Name> getNames() const override;
      virtual void resetIndex() override;
    private:
      /// names in sections
      std::vector<Name> names_;
      std::uint64_t names_idx_;  ///< current index in m_names

      bool use_random_phase0_;
      std::uint_fast64_t seed_;    
      double phase0_;
      
      std::uniform_real_distribution<> dist_;
      std::mt19937 gen_;
    };
  } // namespace pulseseq
} // namespace dnpsoup

#endif


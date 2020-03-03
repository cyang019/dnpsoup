#include "configure_dnpsoup.h"
#include "dnpsoup_core/pulseseq/Section.h"
#include "dnpsoup_core/errors.h"


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    Section::Section()
      : SubSequenceInterface(), names_idx_(0), phase0_(0.0)
    {}

    Section::Section(const vector<Name> &component_names)
      : SubSequenceInterface(), names_(component_names), names_idx_(0),
      use_random_phase0_(false), seed_(0u), phase0_(0.0)
    {}

    Section::Section(std::uint64_t sz, const vector<Name> &component_names)
      : SubSequenceInterface(sz), names_(component_names), names_idx_(0),
      use_random_phase0_(false), seed_(0u), phase0_(0.0)
    {}

    Section::Section(
        std::uint64_t sz, const std::vector<Name> &comp_names, 
        bool flag_phase0)
      : SubSequenceInterface(sz), names_(comp_names), names_idx_(0),
      use_random_phase0_(flag_phase0), seed_(0u), phase0_(0.0),
      dist_(0.0, 360.0), gen_(seed_)
    {
      phase0_ = dist_(gen_) * use_random_phase0_;
    }

    Section::Section(std::uint64_t sz, const std::vector<Name> &comp_names, 
        bool flag_phase0, std::uint64_t seed)
      : SubSequenceInterface(sz), names_(comp_names), names_idx_(0),
      use_random_phase0_(flag_phase0), seed_(seed),
      dist_(0.0, 360.0), gen_(seed)
    {
      phase0_ = dist_(gen_) * use_random_phase0_;
    }

    std::unique_ptr<SubSequenceInterface> Section::copy() const
    {
      return make_unique<Section>(*this);
    }

    void Section::resetIndex(
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections)
    {
      names_idx_ = 0;
      idx_ = 0;
      for(const auto &name: names_){
        sections->at(name)->resetIndex(sections);
      }
    }

    std::tuple<Component, std::uint64_t, std::uint64_t> Section::next(
        std::map<Name, Component> *components,
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(idx_ >= sz_){    
        idx_ = 0;
        names_idx_ = 0;
        return make_tuple(Component(), 0, sz_);
      }

      auto [comp, comp_size, sub_idx] = 
        increment_within_section_(components, sections);
      /// new iteration
      while(sub_idx >= names_.size()
          && idx_ < sz_){
        ++idx_;
        if(idx_ >= sz_) break;
        names_idx_ = 0;
        for(const auto &name : names_){
          if(sections->find(name) == sections->end()){
            const auto err_msg = name + " not found in sections.";
            throw ::dnpsoup::PulseSequenceError(err_msg);
          }
          (sections->at(name))->resetIndex(sections);
        }
        std::tie(comp, comp_size, sub_idx) = 
          increment_within_section_(components, sections);
      }
#ifdef DNPSOUP_VERBOSE
          cout << "comp size: " << comp_size 
               << "  names_idx_: " << sub_idx 
               << "  idx_: " << idx_ << endl;
#endif
      return make_tuple(comp, comp_size, idx_);
    }

    SequenceType Section::type() const 
    { return SequenceType::SectionType; }

    // return subsequence names
    std::vector<Name> Section::getNames() const
    {
      return names_;
    }

    std::tuple<Component, std::uint64_t, std::uint64_t> 
      Section::increment_within_section_(
        std::map<Name, Component> *components,
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections)
    {
      auto sub_name = names_[names_idx_];
      if(sections->find(sub_name) == sections->end()){
        const auto err_msg = sub_name + " not found in sections.";
        throw ::dnpsoup::PulseSequenceError(err_msg);
      }
      auto [comp, comp_size, comp_idx] = (sections->at(sub_name))->next(
          components, sections);
      while(comp_idx >= sections->at(sub_name)->size()
          && names_idx_ < names_.size()){
        ++names_idx_;
        if(names_idx_ >= names_.size()) break;
        sub_name = names_[names_idx_];
        std::tie(comp, comp_size, comp_idx) = (sections->at(sub_name))->next(
          components, sections);
      }
      for(auto &[spin_t, emr]: comp){
        emr.phase += phase0_ * use_random_phase0_;
      }
      return make_tuple(comp, comp_size, names_idx_);
    }
  } // namespace pulseseq
} // namespace dnpsoup


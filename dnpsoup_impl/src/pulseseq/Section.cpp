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
      use_random_phase0_(flag_phase0), seed_(0u),
      dist_(0.0, 360.0), gen_(seed_)
    {
      phase0_ = dist_(gen_) * use_random_phase0_;
    }

    Section::Section(std::uint64_t sz, const std::vector<Name> &comp_names, 
        bool flag_phase0, std::uint64_t seed)
      : SubSequenceInterface(sz), names_(comp_names),
      use_random_phase0_(flag_phase0), seed_(seed),
      dist_(0.0, 360.0), gen_(seed)
    {
      phase0_ = dist_(gen_) * use_random_phase0_;
    }

    void Section::resetIndex()
    {
      names_idx_ = 0;
      idx_ = 0;
    }

    std::tuple<Component, std::uint64_t, std::uint64_t> Section::next(
        std::map<Name, Component> *components,
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      /// to track iterations
      if(idx_ >= sz_){    // protected members of SubSequenceInterface (parent)
        idx_ = 0;
        names_idx_ = 0;
        return make_tuple(Component(), 0, sz_);
      }

      /// iterate through Section
      while(idx_ < sz_){
        // per iteration
        if(names_idx_ >= names_.size()){
          ++idx_;
          phase0_ = dist_(gen_);
          // next iteration, reset all indices
          names_idx_ = 0;
          for(const auto &name : names_){
            (sections->at(name))->resetIndex();
          }
          // increment idx_
          continue;
        }

        // within the iteration
        std::uint64_t component_idx = 0;
        Component p;
        while(names_idx_ < names_.size()){  // iterate through a vector
          auto sub_name = names_[names_idx_];
          if(sections->find(sub_name) == sections->end()){
            const std::string error_str = "cannot find section " + sub_name + ".";
            throw PulseSequenceError(error_str);
          }

          if ((sections->at(sub_name))->type() == SequenceType::DefaultType)
          {
            throw NotImplementedError("cannot get default sub sequence type.");
          }
          uint64_t comp_size = 0;
          std::tie(p, comp_size, component_idx) = (sections->at(sub_name))->next(components, sections);
          if(component_idx >= sections->at(sub_name)->size()){
            ++names_idx_;
            // go to the next loop to return result
            continue;
          }
          else {
            for(auto &[spin_t, emr] : p){
              emr.phase += phase0_ * use_random_phase0_;
            }
            return make_tuple(p, comp_size, idx_);
          }
        } // inner while loop
      } // outer while loop for iteration

      // finished iteration
      idx_ = 0;
      names_idx_ = 0;
      //cout << "[Termination] section " << this->name << ": " << " idx_: " << idx_ << " " << " names_idx_: " << names_idx_ << endl;
      return make_tuple(Component(), 0, sz_);
    }

    SequenceType Section::type() const 
    { return SequenceType::SectionType; }

    // return subsequence names
    std::vector<Name> Section::getNames() const
    {
      return names_;
    }
  } // namespace pulseseq
} // namespace dnpsoup


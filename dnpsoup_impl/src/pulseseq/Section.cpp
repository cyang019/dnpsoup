#include "dnpsoup_core/pulseseq/Section.h"
#include "dnpsoup_core/errors.h"


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    Section::Section()
      : SubSequenceInterface(), m_names_idx(0)
    {}

    Section::Section(const vector<Name> &component_names)
      : SubSequenceInterface(), m_names(component_names), m_names_idx(0)
    {}

    Section::Section(std::uint64_t sz, const vector<Name> &component_names)
      : SubSequenceInterface(sz), m_names(component_names), m_names_idx(0)
    {}

    void Section::resetIndex()
    {
      m_names_idx = 0;
      m_idx = 0;
    }

    std::tuple<Component, std::uint64_t, std::uint64_t> Section::next(
        std::map<Name, Component> *components,
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      /// to track iterations
      if(m_idx >= m_sz){    // protected members of SubSequenceInterface (parent)
        m_idx = 0;
        m_names_idx = 0;
        //cout << "[Termination] section " << this->name << ": " << " m_idx: " << m_idx << " " << " m_names_idx: " << m_names_idx << endl;
        return make_tuple(Component(), 0, m_sz);
      }

      while(m_idx < m_sz){
        // per iteration
        if(m_names_idx >= m_names.size()){
          ++m_idx;
          // next iteration, reset all indices
          m_names_idx = 0;
          for(const auto &name : m_names){
            (sections->at(name))->resetIndex();
          }
          continue;
        }

        // within the iteration
        std::uint64_t component_idx = 0;
        Component p;
        while(m_names_idx < m_names.size()){  // iterate through a vector
          auto sub_name = m_names[m_names_idx];
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
            ++m_names_idx;
            // go to the next loop to return result
            continue;
          }
          else {
            return make_tuple(p, comp_size, m_idx);
          }
        } // inner while loop
      } // outer while loop for iteration

      // finished iteration
      m_idx = 0;
      m_names_idx = 0;
      //cout << "[Termination] section " << this->name << ": " << " m_idx: " << m_idx << " " << " m_names_idx: " << m_names_idx << endl;
      return make_tuple(Component(), 0, m_sz);
    }

    SequenceType Section::type() const 
    { return SequenceType::SectionType; }

    // return subsequence names
    std::vector<Name> Section::getNames() const
    {
      return m_names;
    }
  } // namespace pulseseq
} // namespace dnpsoup


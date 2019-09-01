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

    std::pair<Component, std::uint64_t> Section::next(
        std::map<Name, Component> *components,
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(m_idx >= m_sz){    // parent protected member
        return make_pair(Component(), m_sz);
      }

      while(m_idx < m_sz){
        if(m_names_idx >= m_names.size()){
          ++m_idx;
          m_names_idx = 0;
          for(const auto &name : m_names){
            (sections->at(name))->resetIndex();
          }
          continue;
        }

        std::uint64_t component_idx = 0;
        Component p;
        while(m_names_idx < m_names.size()){
          auto name = m_names[m_names_idx];
          if(sections->find(name) == sections->end()){
            const std::string error_str = "cannot find section " + name + ".";
            throw PulseSequenceError(error_str);
          }

          std::tie(p, component_idx) = (sections->at(name))->next(components, sections);
          if(component_idx == (sections->at(name)->size())){
            ++m_names_idx;
          }
          else{
            break;
          }
        }

        return make_pair(p, m_idx);
      }

      return make_pair(Component(), m_sz);
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


#include "dnpsoup_core/pulseseq/Section.h"
#include "dnpsoup_core/errors.h"


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    Section::Section()
      : SubSequenceInterface(), m_sub_idx(0)
    {}

    Section::Section(const Name &name, const vector<Name> &component_names)
      : SubSequenceInterface(name), m_names(component_names), m_sub_idx(0)
    {}

    Section::Section(const Name &name, std::uint64_t sz, const vector<Name> &component_name)
      : SubSequenceInterface(name, sz), m_names(component_names), m_sub_idx(0)
    {}

    std::pair<Component, std::uint64_t> Section::next(
        std::map<Name, Component> *components,
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(m_sub_idx >= m_names.size()){
        return make_pair(Component(), m_names.size());
      }
      std::uint64_t component_idx = 0;
      Component p;
      while(m_sub_idx < m_names.size()){
        auto name = m_names[m_sub_idx];
        if(sections->find(name) == sections->end()){
          const std::string error_str = "cannot find section " + name + ".";
          throw PulseSequenceError(error_str);
        }

        std::tie(p, component_idx) = (sections->at(name)).next(components, sections);
        if(component_idx == (sections->at(name).size())){
          ++m_sub_idx;
        }
        else{
          break;
        }
      }
      if(m_sub_idx < m_names.size()){
        return make_pair(p, m_sub_idx);
      }

      return make_pair(Component(), m_names.size());
    }

    SequenceType Section::type() const 
    { return SequenceType::Section; }

    // return subsequence names
    std::vector<Name> Section::getNames() const override
    {
      return m_names;
    }
  } // namespace pulseseq
} // namespace dnpsoup


#include "dnpsoup_core/pulseseq/Pulse.h"


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    Pulse::Pulse()
      : SubSequenceInterface()
    {}

    Pulse::Pulse(const Name &component_name)
      : SubSequenceInterface(), m_component_name(component_name)
    {}

    Pulse::Pulse(std::uint64_t sz, const Name &component_name)
      : SubSequenceInterface(sz), m_component_name(component_name)
    {}

    std::pair<Component, std::uint64_t> Pulse::next(
        std::map<Name, Component> *components,
        [[maybe_unused]] std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(m_idx >= m_sz){
        m_idx = 0;
        return make_pair(Component(), m_sz);
      }

      if(components->find(this->name) != components->end()){
        auto idx = m_idx;
        ++m_idx;
        return make_pair(components->at(name), idx);
      }

      return make_pair(Component(), m_sz);
    }

    SequenceType Pulse::type() const 
    { return SequenceType::PulseType; }

    std::vector<Name> Pulse::getNames() const
    {
      std::vector<Name> res = {m_component_name};
      return res;
    }
  } // namespace pulseseq
} // namespace dnpsoup

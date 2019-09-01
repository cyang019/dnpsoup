#include "dnpsoup_core/pulseseq/Delay.h"


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    Delay::Delay()
      : SubSequenceInterface()
    {}

    Delay::Delay(std::uint64_t sz)
      : SubSequenceInterface(sz)
    {}

    std::pair<Component, std::uint64_t> Delay::next(
        [[maybe_unused]] std::map<Name, Component> *components,
        [[maybe_unused]] std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(m_idx >= m_sz){
        m_idx = 0;
        return make_pair(Component(), m_sz);
      }
      else {
        auto idx = m_idx;
        ++m_idx;
        return make_pair(Component(), idx);
      }
    }

    SequenceType Delay::type() const 
    { return SequenceType::DelayType; }

    std::vector<Name> Delay::getNames() const
    {
      auto res = std::vector<Name>();
      return res;
    }
  } // namespace pulseseq
} // namespace dnpsoup


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

    std::tuple<Component, std::uint64_t, std::uint64_t> Delay::next(
        [[maybe_unused]] std::map<Name, Component> *components,
        [[maybe_unused]] std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(m_idx >= m_sz){
        m_idx = 0;
        return make_tuple(Component(), 0, m_sz);
      }
      else {
        m_idx = m_sz;
        return make_tuple(Component(), m_sz, 0);
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


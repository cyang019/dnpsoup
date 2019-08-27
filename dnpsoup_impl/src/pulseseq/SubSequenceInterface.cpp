#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"


namespace dnpsoup {
  namespace pulseseq {
    SubSequenceInterface::SubSequenceInterface()
      : name("default"), m_sz(0), m_idx(0)
    {}

    SubSequenceInterface::SubSequenceInterface(const Name &name)
      : name(name), m_sz(0), m_idx(0)
    {}

    SubSequenceInterface::SubSequenceInterface(const Name &name, std::uint64_t sz)
      : name(name), m_sz(sz), m_idx(0)
    {}

    void SubSequenceInterface::setParam(const Name &name, double val)
    {
      if(m_params.find(name) == m_params.end()){
        m_params.insert({name, val});
      } else {
        m_params[name] = val;
      }
    }
  } // namespace pulseseq
} // namespace dnpsoup

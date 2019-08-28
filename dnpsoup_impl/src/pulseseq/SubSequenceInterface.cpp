#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"
#include "dnpsoup_core/pulseseq/Delay.h"
#include "dnpsoup_core/pulseseq/Pulse.h"
#include "dnpsoup_core/pulseseq/Section.h"
#include <iostream>
#include <sstream>

using namespace std;


namespace dnpsoup {
  namespace pulseseq {
    std::ostream& operator<<(std::ostream &os, const SequenceType &t)
    {
      switch(t){
        case SequenceType::Pulse:
          os << "Pulse";
          break;
        case SequenceType::Delay:
          os << "Delay";
          break;
        case SequenceType::Section:
          os << "Section";
          break;
        default:
          os << "Default";
          break;
      }
      return os;
    }

    std::istream& operator>>(std::istream &is, SequenceType &t)
    {
      string name;
      is >> name;
      if(name == "Pulse"){
        t = SequenceType::Pulse;
      } 
      else if(name == "Delay"){
        t = SequenceType::Delay;
      }
      else if(name == "Section"){
        t = SequenceType::Section;
      }
      else {
        t = SequenceType::Default;
      }
      return is;
    }

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

    std::unique_ptr<SubSequenceInterface> genPtrSequence(std::istream &is)
    {
    }

    std::ostream operator<<(std::ostream &os, const std::unique_ptr<SubSequenceInterface> &uptr_seq)
    {
      os << "  " << uptr_seq->type() << " " uptr_seq->size() << "\n";
    }
  } // namespace pulseseq
} // namespace dnpsoup

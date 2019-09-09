#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/pulseseq/Delay.h"
#include "dnpsoup_core/pulseseq/Pulse.h"
#include "dnpsoup_core/pulseseq/Section.h"
#include <string>
#include <limits>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <utility>

using namespace std;


namespace dnpsoup {
  namespace pulseseq {
    std::string toString(const SequenceType &t)
    {
      switch(t){
      case SequenceType::PulseType:
        return "Pulse";
      case SequenceType::DelayType:
        return "Delay";
      case SequenceType::ChirpType:
        return "Chirp";
      case SequenceType::SectionType:
        return "Section";
      default:
        return "Default";
      };
      return "NotDefined";
    }

    SequenceType toSequenceType(const std::string &s)
    {
      if(s == "Pulse"){
        return SequenceType::PulseType;
      }
      else if(s == "Delay"){
        return SequenceType::DelayType;
      }
      else if(s == "Chirp"){
        return SequenceType::ChirpType;
      }
      else if(s == "Section"){
        return SequenceType::SectionType;
      }
      return SequenceType::DefaultType;
    }

    std::ostream& operator<<(std::ostream &os, const SequenceType &t)
    {
      switch(t){
        case SequenceType::PulseType:
          os << "Pulse";
          break;
        case SequenceType::DelayType:
          os << "Delay";
          break;
        case SequenceType::SectionType:
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
        t = SequenceType::PulseType;
      } 
      else if(name == "Delay"){
        t = SequenceType::DelayType;
      }
      else if(name == "Section"){
        t = SequenceType::SectionType;
      }
      else {
        t = SequenceType::DefaultType;
      }
      return is;
    }

    SubSequenceInterface::SubSequenceInterface()
      : m_sz(1), m_idx(0)
    {}

    SubSequenceInterface::SubSequenceInterface(std::uint64_t sz)
      : m_sz(sz), m_idx(0)
    {}

    void SubSequenceInterface::setParam(const Name &name, double val)
    {
      if(m_params.find(name) == m_params.end()){
        m_params.insert({name, val});
      } else {
        m_params[name] = val;
      }
    }

    std::vector<Name> SubSequenceInterface::getParamNames() const
    {
      vector<Name> res;
      for(const auto &p_pair : m_params){
        res.push_back(p_pair.first);
      }
      return res;
    }
  } // namespace pulseseq
} // namespace dnpsoup
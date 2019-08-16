#include "dnpsoup_core/pulseseq/Segment.h"
#include "dnpsoup_core/errors.h"
#include <sstream>
#include <iomanip>
#include <tuple>

using namespace std;


namespace dnpsoup {
  Segment::Segment()
    : m_inc(1.0e-9), m_repetition(0), m_comp_index(0), m_rep_index(0)
  {}

  Segment::Segment(double inc)
    : m_inc(inc), m_repetition(1), m_comp_index(0), m_rep_index(0)
  {}

  Segment& Segment::addComponent(const PulseComponent &p)
  {
    m_components.push_back(p);
    return *this;
  }

  std::tuple<
    std::map<SpinType, PulsePacket>, 
    std::uint64_t> Segment::next()
  {
    auto res = this->nextInOneRep();
    while(m_comp_index >= m_components.size() && m_rep_index < m_repetition){
      auto temp_rep_idx = m_rep_index;
      this->resetIndex();
      m_rep_index = temp_rep_idx + 1;
      res = this->nextInOneRep();
    }
    return make_tuple(std::get<0>(res), m_rep_index);
  }

  const PulseComponent& Segment::operator[](std::size_t idx) const
  {
    if(idx > m_components.size()){
      throw IndexError("vector of PulseComponent index out of range.");
    }
    return m_components[idx];
  }

  PulseComponent& Segment::operator[](std::size_t idx)
  {
    if(idx > m_components.size()){
      throw IndexError("vector of PulseComponent index out of range.");
    }
    return m_components[idx];
  }

  Segment& Segment::resetIndex()
  {
    for(auto &c : m_components){
      c.resetIndex();
    }
    m_comp_index = 0;
    m_rep_index = 0;
    return *this;
  }

  std::tuple<
    std::map<SpinType, PulsePacket>, 
    std::uint64_t> Segment::nextInOneRep()
  {
    if(m_comp_index >= m_components.size()){
      auto placeholder = std::map<SpinType, PulsePacket>();
      return make_tuple(placeholder, 0);
    }
    auto res = m_components[m_comp_index].next();
    while(std::get<1>(res) >= m_components[m_comp_index].getCount() 
        && m_comp_index < m_components.size()){
      ++m_comp_index;
      res = m_components[m_comp_index].next();
    }
    return res;
  }

  std::istream& operator>>(std::istream &is, Segment &s)
  {
    string line;
    string pulse_components_str = "";
    while(getline(is, line)){
      if(line.find_first_not_of(" \t") == std::string::npos)
        continue;
      istringstream iss(line);
      string word;
      iss >> word;
      if(word == "SegmentEnd" 
          || word == "segmentend"
          || word == "SEGMENTEND")
      {
        break;
      }
      else if(word == "Segment" 
          || word == "segment" || word == "SEGMENT")
      {
        while(iss >> word){
          if(word == "repeat" || word == "Repeat" || word == "REPEAT"){
            iss >> s.m_repetition;
          }
          else if(word == "increment" || word == "Increment" 
              || word == "INCREMENT")
          {
            iss >> s.m_inc;
          }
        }
      }
      else{
        iss.clear();
        pulse_components_str += line + "\n";
      }
    }
    PulseComponent p;
    istringstream comp_ss(pulse_components_str);
    while(comp_ss >> p){
      s.addComponent(p);
    }

    return is;
  }

  std::ostream& operator<<(std::ostream &os, const Segment &s)
  {
    auto ss = os.precision();
    os << setprecision(std::numeric_limits<double>::max_digits10);
    os << "  Segment " << " repeat " << s.m_repetition << " times, "
       << " increment " << s.m_inc << " seconds\n";
    for(const auto &c : s.m_components)
    {
      os << c << "\n";
    }
      
    os << "  SegmentEnd";
    os << setprecision(ss);
    return os;
  }
} // namespace dnpsoup

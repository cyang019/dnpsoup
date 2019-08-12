#include "dnpsoup_core/pulseseq/Segment.h"
#include "dnpsoup_core/errors.h"
#include <tuple>

using namespace std;


namespace dnpsoup {
  Segment::Segment(double inc)
    : m_inc(inc), m_repetition(1), m_comp_index(0), m_rep_index(0)
  {}

  Segment& Segment::addComponent(const PulseComponent &p)
  {
    m_components.push_back(p);
    return *this;
  }

  std::tuple<
    std::unordered_map<SpinType, PulsePacket>, 
    std::uint64_t> Segment::next()
  {
    auto res = this->nextInOneRep();
    while(m_comp_index >= m_components.size() && m_rep_index < m_repetition){
      auto temp_rep_idx = m_rep_index;
      this->resetIndex();
      m_rep_index = temp_rep_idx + 1;
      res = this->nextInOneRep();
    }
    return res;
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
    std::unordered_map<SpinType, PulsePacket>, 
    std::uint64_t> Segment::nextInOneRep()
  {
    if(m_comp_index >= m_components.size()){
      auto placeholder = std::unordered_map<SpinType, PulsePacket>();
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
} // namespace dnpsoup

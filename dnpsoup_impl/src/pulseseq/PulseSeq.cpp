#include "dnpsoup_core/pulseseq/PulseSeq.h"
#include <sstream>
#include <iomanip>
#include <iostream>

using namespace std;


namespace dnpsoup {
  PulseSeq::PulseSeq()
    : m_seg_index(0)
  {}

  PulseSeq& PulseSeq::addSegment(const Segment &s)
  {
    m_segments.push_back(s);
    return *this;
  }

  std::tuple<std::map<SpinType, PulsePacket>, 
             std::uint64_t>
  PulseSeq::next()
  {
    auto p_default = std::map<SpinType, PulsePacket>();
    if(m_seg_index >= m_segments.size()){
      return make_tuple(p_default, m_segments.size());
    }

    auto res = m_segments[m_seg_index].next();
    while(std::get<1>(res) >= m_segments[m_seg_index].getRepetition()
        && m_seg_index < m_segments.size()){
      ++m_seg_index;
      res = m_segments[m_seg_index].next();
    }
    if(m_seg_index < m_segments.size()){
      return make_tuple(std::get<0>(res), m_seg_index);
    }

    return make_tuple(std::get<0>(res), m_segments.size());
  }

  const Segment& PulseSeq::operator[](std::size_t pos) const
  {
    if(pos >= m_segments.size()){
      throw IndexError("PulseSeq index out of range for segments.");
    }
    return m_segments[pos];
  }

  Segment& PulseSeq::operator[](std::size_t pos)
  {
    if(pos >= m_segments.size()){
      throw IndexError("PulseSeq index out of range for segments.");
    }
    return m_segments[pos];
  }

  std::istream& operator>>(std::istream &is, PulseSeq &pseq)
  {
    string line;
    bool in_sequence = false;
    string sequence_content = "";
    while(getline(is, line)){
      if(line.find_first_not_of(" \t") == std::string::npos)
        continue;
      istringstream iss(line);
      string word;
      iss >> word;
      if(word == "PulseSequenceEnd"
          || word == "pulsesequenceend"
          || word == "PULSESEQUENCEEND"){
        break;
      }
      else if(word == "PulseSequence"
          || word == "pulsesequence"
          || word == "PULSESEQUENCE"){
        in_sequence = true;
      }

      if(in_sequence){
        sequence_content += line + "\n";
      }
    }
    if(sequence_content.find_first_not_of(" \t\n") != std::string::npos)
    {
      istringstream segment_ss(sequence_content);
      Segment s;
      while(segment_ss >> s){
        pseq.m_segments.push_back(s);
      }
    }

    return is;
  }

  std::ostream& operator<<(std::ostream &os, const PulseSeq &pseq)
  {
    auto ss = os.precision();
    os << setprecision(std::numeric_limits<double>::max_digits10);

    os << "PulseSequence\n";
    for(const auto &s : pseq.m_segments){
      os << s << "\n";
    }
    os << "PulseSequenceEnd\n";

    os << setprecision(ss);
    return os;
  }
} // namespace dnpsoup

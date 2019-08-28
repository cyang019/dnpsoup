#include "dnpsoup_core/pulseseq/PulseSequence.h"
#include "dnpsoup_core/errors.h"
#include <limits>
#include <sstream>
#include <iomanip>
#include <memory>
#include <iostream>

using namespace std;


namespace dnpsoup {
namespace pulseseq{
  PulseSequence::PulseSequence()
    : name("default"), m_seg_index(0), m_idx(0)
  {}
  PulseSequence::PulseSequence(double inc)
    : name("default"), m_inc(1.0e-9), m_idx(0)
  {}

  PulseSequence::PulseSequence(const std::string &name, double inc)
    : name(name), m_inc(inc), m_idx(0)
  {}

  PulseSequence::PulseSequence(const std::string &name)
    : name(name), m_inc(1.0e-9), m_idx(0)
  {}

  PulseSequence::~PulseSequence()
  {}
  
  const Component& PulseSequence::getComponent(const Name &n) const
  {
    return m_components.at(n);
  }

  PulseSequence& PulseSequence::set(const Name &n, const Component &c)
  {
    m_components[n] = c;
    return *this;
  }

  PulseSequence& PulseSequence::set(
      const Name &n, 
      unique_ptr<SubSequenceInterface> uptr_seq)
  {
    m_sections[n] = std::move(uptr_seq);
    return *this;
  }

  double PulseSequence::getParam(const Name &seq_name, const Name &param_name) const
  {
    return m_sections.at(seq_name)->getParam(param_name);
  }

  PulseSequence& PulseSequence::setParam(const Name &seq_name, const Name &param_name, double value)
  {
    return m_sections.at(seq_name)->setParam(param_name, value);
  }

  PulseSequence& PulseSequence::set(const vector<Name> &seq_names)
  {
#ifndef NDEBUG
    for(const auto &name : seq_names){
      if(m_sections.find(name) == m_sections.end()){
        throw ::dnpsoup::PulseSequenceError("sub sequence not defined yet.");
      }
    }
#endif
    m_sections_in_order = seq_names;
    return *this;
  }

  PulseSequence& PulseSequence::remove(const Name &seq_name)
  {
    auto it = m_sections.find(seq_name);
    if(it != m_sections.end()){
      m_sections.erase(it);
    }
    m_idx = 0;
    return *this;
  }

  Component PulseSequence::next()
  {
    if (m_idx == m_sections_in_order.size()){
      return Component();
    }

    auto name = m_sections_in_order[m_idx];
    auto [comp, idx] = m_sections[name]->next(&m_components, &m_sections);
    while(idx == m_sections[name]->size()){
      ++m_idx;
      if(m_idx == m_sections_in_order.size()) break;
      name = m_sections_in_order[m_idx];
      std::tie(comp, idx) = m_sections[name]->next(&m_components, &m_sections);
    }

    return comp;
  }

  std::istream& operator>>(std::istream &is, PulseSequence &pseq)
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
      if(word == "PulseSequenceuenceEnd"
          || word == "pulsesequenceend"
          || word == "PULSESEQUENCEEND"){
        break;
      }
      else if(word == "PulseSequenceuence"
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

  std::ostream& operator<<(std::ostream &os, const PulseSequence &pseq)
  {
    auto ss = os.precision();
    os << setprecision(std::numeric_limits<double>::max_digits10);
    os << pseq.name << "\n";
    os << "Increment " << pseq.m_inc << "\n\n";
    for(const auto &comp_pair : pseq.m_components){
      os << "Component " << comp_pair.first << "\n";
      os << comp_pair.second;
      os << "ComponentEnd\n";
    }
    os << "\n";

    for(const auto &section_pair : pseq.m_sections){
      os << "Section " << section_pair.first
         << section_pair.second->size() << "\n";
      os << section_pair.second;
      os << "SectionEnd\n";
    }
    os << "\n";

    os << "PulseSequence\n";
    for(const auto &name : pseq.m_sections_in_order){
      os << "  " << name << "\n";
    }
    os << "PulseSequenceEnd\n";

    os << setprecision(ss);
    return os;
  }
} // namespace pulseseq
} // namespace dnpsoup

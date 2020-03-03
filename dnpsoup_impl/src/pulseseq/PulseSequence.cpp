#include "dnpsoup_core/pulseseq/PulseSequence.h"
#include "dnpsoup_core/pulseseq/seq_common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/pulseseq/Pulse.h"
#include "dnpsoup_core/pulseseq/Delay.h"
#include "dnpsoup_core/pulseseq/Section.h"
#include <limits>
#include <sstream>
#include <iomanip>
#include <memory>
#include <iostream>
#include <limits>
#include <cstdint>
#include "json.hpp"

using namespace std;


namespace dnpsoup {
namespace pulseseq{
  PulseSequence::PulseSequence()
    : name("default"), m_inc(1.0e-9), m_idx(0)
  {
  }

  PulseSequence::PulseSequence(double inc)
    : name("default"), m_inc(inc), m_idx(0)
  {
  }

  PulseSequence::PulseSequence(const std::string &name, double inc)
    : name(name), m_inc(inc), m_idx(0)
  {
  }

  PulseSequence::PulseSequence(const std::string &name)
    : name(name), m_inc(1.0e-9), m_idx(0)
  {
  }

  PulseSequence::PulseSequence(const PulseSequence &seq)
    : name(seq.name), m_components(seq.m_components),
    m_sections_in_order(seq.m_sections_in_order),
    m_inc(seq.m_inc), m_idx(seq.m_idx)
  {
    for(const auto &[name, uptr_seq] : seq.m_sections){
      m_sections[name] = uptr_seq->copy();
    }
  }

  PulseSequence& PulseSequence::operator=(const PulseSequence &seq)
  {
    name = seq.name;
    m_components = seq.m_components;
    m_sections_in_order = seq.m_sections_in_order;
    m_inc = seq.m_inc;
    m_idx = seq.m_idx;
    for(const auto &[name, uptr_seq] : seq.m_sections){
      m_sections[name] = uptr_seq->copy();
    }

    return *this;
  }

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
    m_sections[n]->name = n;
    return *this;
  }

  std::vector<Name> PulseSequence::getNames(const Name &name) const
  {
    return m_sections.at(name)->getNames();
  }

  double PulseSequence::getParam(const Name &seq_name, const Name &param_name) const
  {
    return m_sections.at(seq_name)->getParam(param_name);
  }

  PulseSequence& PulseSequence::setParam(const Name &seq_name, const Name &param_name, double value)
  {
    m_sections.at(seq_name)->setParam(param_name, value);
    return *this;
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

  std::tuple<Component, std::uint64_t, std::uint64_t> PulseSequence::next()
  {
    if (m_idx >= m_sections_in_order.size()){
      m_idx = m_sections_in_order.size();
      return std::make_tuple(Component(), 0, m_sections_in_order.size());
    }

    auto name = m_sections_in_order[m_idx];
    auto [comp, comp_size, idx] = m_sections[name]->next(&m_components, &m_sections);
    while(idx == m_sections[name]->size()){
      ++m_idx;
      if(m_idx == m_sections_in_order.size()) break;
      name = m_sections_in_order[m_idx];
      std::tie(comp, comp_size, idx) = m_sections[name]->next(&m_components, &m_sections);
    }

    return std::make_tuple(comp, comp_size, m_idx);
  }

  std::istream& operator>>(std::istream &is, PulseSequence &pseq)
  {
    json j;
    is >> j;
    if(j.empty()){
      pseq = PulseSequence();
      return is;
    }
    if(j.find("pulse_sequence") != j.end()){
      j = j["pulse_sequence"];
    }
    if(j.find("increment") != j.end()){
      pseq.m_inc = j["increment"].get<double>();      
    }
    if(j.find("components") != j.end()){
      json components_js = j["components"];
      for(auto& [comp_name, comp_emrs_js] : components_js.items()){
        auto comp = Component();
        for(auto &[emr_name_str, emr_js] : comp_emrs_js.items()){
          auto emr_name = toSpinType(emr_name_str);
#ifndef NDEBUG
          const std::string err_str_postfix = 
            " not found for " + emr_name_str + " of " + comp_name;
          
          if(emr_js.find("frequency") == emr_js.end()){
            throw PulseSequenceError("frequency" + err_str_postfix);
          }
          if(emr_js.find("phase") == emr_js.end()){
            throw PulseSequenceError("phase" + err_str_postfix);
          }
          if(emr_js.find("offset") == emr_js.end()){
            throw PulseSequenceError("offset" + err_str_postfix);
          }
#endif
          const double freq = emr_js["frequency"].get<double>();
          const double phase = emr_js["phase"].get<double>();
          const double offset = emr_js["offset"].get<double>();
          comp[emr_name] = EMRadiation(freq, phase, offset);
        }
        pseq.m_components[comp_name] = comp;
      }
    }

    if(j.find("sections") != j.end()){
      for(auto& [section_name, section_info_js] : j["sections"].items()){
        const std::string name_str = section_name;
#ifndef NDEBUG
        if(section_info_js.find("type") == section_info_js.end()){
          throw PulseSequenceError("type not found in section");
        }
        if(section_info_js.find("size") == section_info_js.end()){
          throw PulseSequenceError("size not found in section");
        }
#endif
        const std::string seq_type_str = section_info_js["type"].get<std::string>();
        const std::uint64_t sz = section_info_js["size"].get<std::uint64_t>();
        std::vector<std::string> member_names;
        for(string name : section_info_js["names"]){
          member_names.push_back(name);
        }
        std::map<std::string, double> params;
        if(section_info_js.find("params") != section_info_js.end()
            && (!section_info_js["params"].empty())){
          for(auto& [key, value_js] : section_info_js["params"].items()){
            params[key] = value_js.get<double>();
          }
        }
        auto seq_type = toSequenceType(seq_type_str);
        unique_ptr<SubSequenceInterface> ptr_member;
        if(seq_type == SequenceType::PulseType){
          if(member_names.size() < 1){
#ifndef NDEBUG
            throw PulseSequenceError("Pulse EMR component missing.");
#endif
            continue;
          }
          ptr_member = make_unique<Pulse>(sz, member_names[0]);
        }
        else if(seq_type == SequenceType::DelayType){
          ptr_member = make_unique<Delay>(sz);
        }
        else if(seq_type == SequenceType::SectionType){
          bool use_phase0 = false;
          std::uint_fast64_t seed = 0u;
          if(section_info_js.find("phase0") != section_info_js.end()){
            if(section_info_js["phase0"].find("reset") != section_info_js["phase0"].end()){
              use_phase0 = section_info_js["phase0"]["reset"].get<bool>();
            }
            if(section_info_js["phase0"].find("seed") != section_info_js["phase0"].end()){
              seed = section_info_js["phase0"]["seed"].get<std::uint_fast64_t>();
            }
          }
          ptr_member = make_unique<Section>(sz, member_names, use_phase0, seed);
        }
        else{   // not implemented.
#ifndef NDEBUG
          throw PulseSequenceError("Such Sequence Type is not implemented.");
#endif
          continue;
        }
        pseq.m_sections[name_str] = std::move(ptr_member);
        pseq.m_sections[name_str]->name = name_str;
        for(auto& [param_name, val] : params){
          pseq.m_sections[name_str]->setParam(param_name, val);
        } // fill parameters
      } // for section in sections
    } // sections

    pseq.m_sections_in_order.clear();
#ifndef NDEBUG
    if(j.find("sequence") == j.end()){
      throw PulseSequenceError("sequence of components missing.");
    }
#endif
    for(auto &name : j["sequence"]){
      pseq.m_sections_in_order.push_back(name.get<std::string>());
    }

    return is;
  }

  std::ostream& operator<<(std::ostream &os, const PulseSequence &pseq)
  {
    auto ss = os.precision();
    os << setprecision(std::numeric_limits<double>::max_digits10);
    json seq_json = json::object();
    seq_json["increment"] = pseq.m_inc;
    json components_js = json::object();

    for(const auto &comp_pair : pseq.m_components){
      components_js[comp_pair.first] = {};
      for(const auto &emr_pair : comp_pair.second){
        components_js[comp_pair.first][toString(emr_pair.first)] = 
        {
          {"frequency", emr_pair.second.freq},
          {"phase", emr_pair.second.phase},
          {"offset", emr_pair.second.offset}
        };
      }
    }
    seq_json["components"] = components_js;

    json sections_js = json::object();
    for(const auto &section_pair : pseq.m_sections){
      auto names = section_pair.second->getNames();
      json names_js(names);
      json params_js = json::object();
      auto param_names = section_pair.second->getParamNames();
      for(const auto &pname : param_names){
        params_js[pname] = section_pair.second->getParam(pname);
      }
      sections_js[section_pair.first] = {
          {"type", toString(section_pair.second->type())},
          {"size", section_pair.second->size()},
          {"names", names_js},
          {"params", params_js}
      };
    }
    seq_json["sections"] = sections_js;

    json sequence(pseq.m_sections_in_order);
    seq_json["sequence"] = sequence;

    os << seq_json.dump(4) << "\n";
    os << setprecision(ss);
    return os;
  }
} // namespace pulseseq
} // namespace dnpsoup

#include "dnpsoup_core/pulseseq/PulseComponent.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>

using namespace std;


namespace dnpsoup {
    PulseComponent::PulseComponent()
      : m_count(0), m_index(0)
    {}

    PulseComponent::PulseComponent(std::uint64_t cnt)
      : m_count(cnt), m_index(0)
    {}

    PulseComponent& PulseComponent::setSize(std::uint64_t cnt)
    {
      m_count = cnt;
      return *this;
    }

    SequenceType PulseComponent::type() const override
    {
      return SequenceType::Pulse;
    }

    const PulsePacket& PulseComponent::getChannel(const SpinType &t) const
    {
      return m_channels.at(t);
    }

    PulseComponent& PulseComponent::setChannel(const SpinType &t, const PulsePacket &p)
    {
      m_channels[t] = p;
      return *this;
    }

    PulseComponent& PulseComponent::removeChannel(const SpinType &t)
    {
      auto it = m_channels.find(t);
      if(it != m_channels.end()){
        m_channels.erase(it);
      }
      return *this;
    }

    std::pair<
      std::map<SpinType, PulsePacket>,
      std::uint64_t> PulseComponent::next()
    {
      if(m_index < m_count){
        auto res = std::make_pair(m_channels, m_index);
        ++m_index;
        return res;
      } else {    // beyond range
        auto res = std::make_pair(
            std::map<SpinType, PulsePacket>(),
            m_count);
        return res;
      }
    }

    PulseComponent& PulseComponent::resetIndex()
    {
      m_index = 0;
      return *this;
    }

    std::istream& operator>>(std::istream &is, PulseComponent &p)
    {
      string line;
      while(getline(is, line)){
        if(line.find_first_not_of(" \t") == std::string::npos)
          continue;
        istringstream iss(line);
        string word;
        iss >> word;
        if(word == "PulseComponentEnd" 
            || word == "pulsecomponentend" 
            || word == "PULSECOMPONENTEND")
          break;
        if(word == "PulseComponent"
            || word == "pulsecomponent"
            || word == "PULSECOMPONENT")
        {
          while(iss >> word){
            if(word == "steps" || word == "STEPS" || word == "Steps"){
              iss >> p.m_count;
            }
            else if(word == "increment" 
                || word == "INCREMENT" 
                || word == "Increment"){
              iss >> p.m_inc;
            }
          } // while
        } // if
        else if(isSpinTypeStr(word)) {
          auto t = toSpinType(word);
          PulsePacket packet;
          iss >> packet;
          p.setChannel(t, packet);
        }
      }

      return is;
    }

    std::ostream& operator<<(std::ostream &os, const PulseComponent &p)
    {
      auto ss = os.precision();
      os << setprecision(std::numeric_limits<double>::max_digits10);
      os << "    PulseComponent " 
         << "steps " << p.m_count 
         << " increment " << p.m_inc << "\n";
      for(const auto &spin_info : p.m_channels){
        os << "      " << toString(spin_info.first) 
           << " " << spin_info.second << "\n";
      }
      os << "    PulseComponentEnd";
      os << setprecision(ss);
      return os;
    }
} // namespace dnpsoup

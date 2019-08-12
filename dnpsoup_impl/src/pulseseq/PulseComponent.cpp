#include "dnpsoup_core/pulseseq/PulseComponent.h"
#include <limits>
#include <iostream>


namespace dnpsoup {
  PulseComponent::PulseComponent(std::uint64_t cnt, double inc)
      : m_count(cnt), m_inc(inc), m_index(0)
    {}

    PulseComponent& PulseComponent::setIncrement(double inc)
    {
#ifndef NDEBUG
      if(inc < std::numeric_limits<double>::epsilon()){
        std::cerr << "Increment too small." << std::endl;
      }
#endif
      m_inc = inc;
      return *this;
    }

    PulseComponent& PulseComponent::setCount(std::uint64_t cnt)
    {
      m_count = cnt;
      return *this;
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

    std::tuple<
      std::unordered_map<SpinType, PulsePacket>,
      std::uint64_t> PulseComponent::next()
    {
      if(m_index != m_count){
        auto res = std::make_tuple(m_channels, m_index);
        ++m_index;
        return res;
      } else {    // beyond range
        auto res = std::make_tuple(
            std::unordered_map<SpinType, PulsePacket>(),
            m_index);
        return res;
      }
    }

    PulseComponent& PulseComponent::resetIndex()
    {
      m_index = 0;
      return *this;
    }
} // namespace dnpsoup

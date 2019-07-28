#include "dnpsoup_core/spinsys/SpinPacket.h"


namespace dnpsoup {
  SpinPacket::SpinPacket(
      std::unique_ptr<InteractionInterface> ptr_interface, 
      const Property &p, const Euler<> &e)
    : m_ptr_interface(std::move(ptr_interface)), m_property(p), m_e(e)
  {}

  SpinPacket& SpinPacket::rotate(const Euler<> &e)
  {
    m_e = m_e * e;
    return *this;
  }

  SpinPacket& SpinPacket::setPropertyValue(const ValueName &vname, double val)
  {
    m_property.set(vname, val);
    return *this;
  }

  double SpinPacket::getPropertyValue(const ValueName &vname) const
  {
    return m_property.get(vname);
  }

  MatrixCxDbl SpinPacket::genMatrix() const
  {
    return m_ptr_interface->genMatrix(m_property, m_e);
  }

  SpinPacketCollection& SpinPacketCollection::add(SpinPacket &&sp)
  {
    m_packets.push_back(std::move(sp));
    return *this;
  }

  SpinPacketCollection& SpinPacketCollection::rotate(const Euler<> &e)
  {
    for(auto &sp : m_packets){
      sp.rotate(e);
    }
    return *this;
  }

  SpinPacketCollection& SpinPacketCollection::setPropertyValue(const ValueName &vname, double val)
  {;;
    for(auto &sp : m_packets){
      sp.setPropertyValue(vname, val);
    }
    return *this;
  }

  MatrixCxDbl SpinPacketCollection::genMatrix() const
  {
    if(m_packets.size() == 0) return MatrixCxDbl();
    MatrixCxDbl res = m_packets[0].genMatrix();
    for(size_t i = 1; i < m_packets.size(); ++i){
      res += m_packets[i].genMatrix();
    }
    return res;
  }
} // namespace dnpsoup

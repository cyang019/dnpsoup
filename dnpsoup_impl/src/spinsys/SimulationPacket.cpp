#include "dnpsoup_core/spinsys/SimulationPacket.h"


namespace dnpsoup {
  SimulationPacket::SimulationPacket(
      std::unique_ptr<InteractionInterface> ptr_interface, 
      const Property &p, const Euler<> &e)
    : m_ptr_interface(std::move(ptr_interface)), m_property(p), m_e(e)
  {}

  SimulationPacket::SimulationPacket(SimulationPacket &&sim_pack) noexcept
    : m_ptr_interface(std::move(sim_pack.m_ptr_interface)),
    m_property(std::move(sim_pack.m_property)),
    m_e(std::move(sim_pack.m_e))
  {}

  SimulationPacket& SimulationPacket::operator=(SimulationPacket &&rhs) noexcept
  {
    m_ptr_interface = std::move(rhs.m_ptr_interface);
    m_property = std::move(rhs.m_property);
    m_e = std::move(rhs.m_e);
    return *this;
  }

  SimulationPacket& SimulationPacket::rotate(const Euler<> &e)
  {
    m_e = m_e * e;
    return *this;
  }

  SimulationPacket& SimulationPacket::setPropertyValue(const ValueName &vname, double val)
  {
    m_property.set(vname, val);
    return *this;
  }

  double SimulationPacket::getPropertyValue(const ValueName &vname) const
  {
    return m_property.get(vname);
  }

  MatrixCxDbl SimulationPacket::genMatrix() const
  {
    return m_ptr_interface->genMatrix(m_property, m_e);
  }

  // class PacketCollection
  PacketCollection& PacketCollection::add(SimulationPacket &&sp)
  {
    m_packets.push_back(std::move(sp));
    return *this;
  }

  PacketCollection& PacketCollection::rotate(const Euler<> &e)
  {
    for(auto &sp : m_packets){
      sp.rotate(e);
    }
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(const ValueName &vname, double val)
  {;;
    for(auto &sp : m_packets){
      sp.setPropertyValue(vname, val);
    }
    return *this;
  }

  MatrixCxDbl PacketCollection::genMatrix() const
  {
    if(m_packets.size() == 0) return MatrixCxDbl();
    MatrixCxDbl res = m_packets[0].genMatrix();
    for(size_t i = 1; i < m_packets.size(); ++i){
      res += m_packets[i].genMatrix();
    }
    return res;
  }
} // namespace dnpsoup

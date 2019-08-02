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
  PacketCollection& PacketCollection::add(const ObservableId &oid, SimulationPacket &&sp)
  {
    m_packets.insert({oid, std::move(sp)});
    return *this;
  }

  PacketCollection& PacketCollection::rotate(const Euler<> &e)
  {
    for(auto &sp_pair : m_packets){
      sp_pair.second.rotate(e);
    }
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(const ValueName &vname, double val)
  {;;
    for(auto &sp_pair : m_packets){
      sp_pair.second.setPropertyValue(vname, val);
    }
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(
        const ObservableId &oid, const ValueName &vname, double val)
  {
    m_packets.at(oid).setPropertyValue(vname, val);
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(
      const InteractionType &t, const SpinId &sid,
      const ValueName &vname, double val)
  {
    auto oid = ObservableId(t, sid);
    m_packets.at(oid).setPropertyValue(vname, val);
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(
      const InteractionType &t, const SpinId &id1, const SpinId &id2,
      const ValueName &vname, double val)
  {
    auto oid = ObservableId(t, id1, id2);
    m_packets.at(oid).setPropertyValue(vname, val);
    return *this;
  }

  double PacketCollection::getPropertyValue(
      const ObservableId &oid, const ValueName &vname) const
  {
    return m_packets.at(oid).getPropertyValue(vname);
  }

  double PacketCollection::getPropertyValue(
      const InteractionType &t, const SpinId &sid, 
      const ValueName &vname) const
  {
    auto oid = ObservableId(t, sid);
    return m_packets.at(oid).getPropertyValue(vname);
  }

  double PacketCollection::getPropertyValue(
      const InteractionType &t, const SpinId &id1, const SpinId &id2,
      const ValueName &vname) const
  {
    auto oid = ObservableId(t, id1, id2);
    return m_packets.at(oid).getPropertyValue(vname);
  }

  MatrixCxDbl PacketCollection::genMatrix() const
  {
    if(m_packets.size() == 0) return MatrixCxDbl();
    auto res = MatrixCxDbl();
    for(const auto &obs : m_packets){
      if(res.nrows() == 0){
        res = obs.second.genMatrix();
      } else{
        res += obs.second.genMatrix();
      }
    }
    return res;
  }

  std::size_t PacketCollection::getNumOfPackets() const
  {
    return m_packets.size();
  }
} // namespace dnpsoup

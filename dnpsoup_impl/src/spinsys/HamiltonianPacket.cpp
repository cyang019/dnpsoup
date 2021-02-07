#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/common.h"
#include <iostream>

using namespace std;

namespace dnpsoup {
  HamiltonianPacket::HamiltonianPacket(
      std::unique_ptr<InteractionInterface> ptr_interface, 
      const Property &p, const Euler<> &e)
    : m_ptr_interface(std::move(ptr_interface)), m_property(p), m_e(e)
  {}

  HamiltonianPacket::HamiltonianPacket(HamiltonianPacket &&sim_pack) noexcept
    : m_ptr_interface(std::move(sim_pack.m_ptr_interface)),
    m_property(std::move(sim_pack.m_property)),
    m_e(std::move(sim_pack.m_e))
  {}

  HamiltonianPacket& HamiltonianPacket::operator=(HamiltonianPacket &&rhs) noexcept
  {
    m_ptr_interface = std::move(rhs.m_ptr_interface);
    m_property = std::move(rhs.m_property);
    m_e = std::move(rhs.m_e);
    return *this;
  }

  HamiltonianPacket& HamiltonianPacket::rotate(const Euler<> &e)
  {
    m_e = m_e * e;
    return *this;
  }

  HamiltonianPacket& HamiltonianPacket::setPropertyValue(const ValueName &vname, double val)
  {
    m_property.set(vname, val);
    return *this;
  }

  double HamiltonianPacket::getPropertyValue(const ValueName &vname) const
  {
    return m_property.get(vname);
  }

  MatrixCxDbl HamiltonianPacket::genMatrix() const
  {
    return m_ptr_interface->genMatrix(m_property, m_e);
  }

  MatrixCxDbl HamiltonianPacket::genMatrix(const Euler<> &e) const
  {
    /// active rotation
    /// The order: [magic angle] * [sample angle] * [spin sys angle] * m_e
    ///                                                                 ^
    ///                                                           Observation Angle
    auto angle = e * m_e;
		//cout << "euler1: " << e << "\n"
    //     << "euler2: " << m_e << "\n"
    //     << "euler1 * euler2: " << angle << "\n" << endl;
    return m_ptr_interface->genMatrix(m_property, angle);
  }

  // class PacketCollection
  PacketCollection& PacketCollection::add(const ObservableId &oid, HamiltonianPacket &&sp)
  {
    m_packets.insert({oid, std::move(sp)});
    return *this;
  }

  PacketCollection& PacketCollection::addBulk(const ObservableId &oid, HamiltonianPacket &&sp)
  {
    m_bulk_packets.insert({oid, std::move(sp)});
    return *this;
  }

  PacketCollection& PacketCollection::rotate(const Euler<> &e)
  {
    for(auto &sp_pair : m_packets){
      sp_pair.second.rotate(e);
    }
    for(auto &sp_pair : m_bulk_packets){
      sp_pair.second.rotate(e);
    }
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(const ValueName &vname, double val)
  {;;
    for(auto &sp_pair : m_packets){
      sp_pair.second.setPropertyValue(vname, val);
    }
    for(auto &sp_pair : m_bulk_packets){
      sp_pair.second.setPropertyValue(vname, val);
    }
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(
        const ObservableId &oid, const ValueName &vname, double val)
  {
    m_packets.at(oid).setPropertyValue(vname, val);
    if(m_bulk_packets.find(oid) != m_bulk_packets.end()){
      m_bulk_packets.at(oid).setPropertyValue(vname, val);
    }
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(
      const InteractionType &t, const SpinId &sid,
      const ValueName &vname, double val)
  {
    auto oid = ObservableId(t, sid);
    m_packets.at(oid).setPropertyValue(vname, val);
    if(m_bulk_packets.find(oid) != m_bulk_packets.end()){
      m_bulk_packets.at(oid).setPropertyValue(vname, val);
    }
    return *this;
  }

  PacketCollection& PacketCollection::setPropertyValue(
      const InteractionType &t, const SpinId &id1, const SpinId &id2,
      const ValueName &vname, double val)
  {
    auto oid = ObservableId(t, id1, id2);
    m_packets.at(oid).setPropertyValue(vname, val);
    if(m_bulk_packets.find(oid) != m_bulk_packets.end()){
      m_bulk_packets.at(oid).setPropertyValue(vname, val);
    }
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

  MatrixCxDbl PacketCollection::genMatrix(const Euler<> &e) const
  {
    if(m_packets.size() == 0) return MatrixCxDbl();
    auto res = MatrixCxDbl();
    for(const auto &obs : m_packets){
      if(res.nrows() == 0){
        res = obs.second.genMatrix(e);
      } else{
        res += obs.second.genMatrix(e);
      }
    }
    return res;
  }

  bool PacketCollection::hasBulk() const
  {
    return m_bulk_packets.size() > 0;
  }

  MatrixCxDbl PacketCollection::genMatrixBulkRemainder() const
  {
    if(m_bulk_packets.size() == 0) return MatrixCxDbl();
    auto res = MatrixCxDbl();
    for(const auto &obs : m_bulk_packets){
      if(res.nrows() == 0){
        res = obs.second.genMatrix();
      } else{
        res += obs.second.genMatrix();
      }
    }
    return res;
  }

  MatrixCxDbl PacketCollection::genMatrixBulkRemainder(const Euler<> &e) const
  {
    if(m_bulk_packets.size() == 0) return MatrixCxDbl();
    auto res = MatrixCxDbl();
    for(const auto &obs : m_bulk_packets){
      if(res.nrows() == 0){
        res = obs.second.genMatrix(e);
      } else{
        res += obs.second.genMatrix(e);
      }
    }
    return res;
  }


  bool PacketCollection::hasPulseSeqComponent(const pulseseq::Component &comp) const
  {
    for(const auto &[spin_t, emr] : comp){
      auto obs_id = ObservableId(InteractionType::EMR, spin_t);
      double prev_freq = this->getPropertyValue(obs_id, ValueName::freq);
      double prev_phase = this->getPropertyValue(obs_id, ValueName::phase);
      double prev_offset = this->getPropertyValue(obs_id, ValueName::offset);
      if(!approxEqual<double>(prev_freq, emr.freq, eps)
          || !approxEqual<double>(prev_phase, emr.phase, eps)
          || !approxEqual<double>(prev_offset, emr.offset, eps)){
        return false;
      }
    }
    return true;
  }

  void PacketCollection::updatePulseSeqComponent(const pulseseq::Component &comp)
  {
    for(const auto &[spin_t, emr] : comp){
      auto obs_id = ObservableId(InteractionType::EMR, spin_t);
      this->setPropertyValue(
          obs_id, ValueName::freq, emr.freq);
      this->setPropertyValue(
          obs_id, ValueName::phase, emr.phase);
      this->setPropertyValue(
          obs_id, ValueName::offset, emr.offset);
    }
  }

  std::vector<ObservableId> PacketCollection::getObservableIds() const
  {
    std::vector<ObservableId> result;
    for(const auto &obs : m_packets){
      result.push_back(obs.first);
    }
    return result;
  }

  std::size_t PacketCollection::getNumOfPackets() const
  {
    return m_packets.size();
  }
} // namespace dnpsoup

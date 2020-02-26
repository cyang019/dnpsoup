#ifndef DNPSOUP_HAMILTONIANPACKET_H
#define DNPSOUP_HAMILTONIANPACKET_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/Observable.h"
#include "dnpsoup_core/pulseseq/seq_common.h"
#include <vector>
#include <tuple>
#include <memory>   // unique_ptr
#include <map>
#include <set>
#include <utility>


namespace dnpsoup {
  class HamiltonianPacket {
  public:
    HamiltonianPacket(std::unique_ptr<InteractionInterface>, 
        const Property &, const Euler<> &);
    HamiltonianPacket(const HamiltonianPacket &) = delete;
    HamiltonianPacket(HamiltonianPacket &&) noexcept;
    HamiltonianPacket& operator=(const HamiltonianPacket &) = delete;
    HamiltonianPacket& operator=(HamiltonianPacket &&) noexcept;
    ~HamiltonianPacket() {}

    // rotate e, then m_e
    HamiltonianPacket& rotate(const Euler<> &e);
    const Euler<>& getEuler() const { return m_e; }

    HamiltonianPacket& setPropertyValue(const ValueName &, double val);
    double getPropertyValue(const ValueName &) const;

    MatrixCxDbl genMatrix() const;
    MatrixCxDbl genMatrix(const Euler<> &e) const;
  private:
    std::unique_ptr<InteractionInterface> m_ptr_interface;
    Property m_property;
    Euler<> m_e;
  };

  class PacketCollection {
  public:
    PacketCollection() {}
    PacketCollection(const PacketCollection &) = delete;
    PacketCollection(PacketCollection &&) noexcept = default;
    PacketCollection& operator=(const PacketCollection &) = delete;
    PacketCollection& operator=(PacketCollection &&) noexcept = default;
    ~PacketCollection() {}

    PacketCollection& add(const ObservableId &, HamiltonianPacket &&);
    PacketCollection& rotate(const Euler<> &e);

    // set all values
    PacketCollection& setPropertyValue(const ValueName &, double val);
    // set individual value
    PacketCollection& setPropertyValue(const ObservableId &, const ValueName &, double val);
    PacketCollection& setPropertyValue(
        const InteractionType &, const SpinId &,
        const ValueName &, double val);
    PacketCollection& setPropertyValue(
        const InteractionType &, const SpinId &, const SpinId &,
        const ValueName &, double val);
    double getPropertyValue(const ObservableId &, const ValueName &) const;
    double getPropertyValue(
        const InteractionType &, const SpinId &, 
        const ValueName &) const;
    double getPropertyValue(
        const InteractionType &, const SpinId &, const SpinId &,
        const ValueName &) const;

    MatrixCxDbl genMatrix() const;
    MatrixCxDbl genMatrix(const Euler<> &e) const;

    bool hasPulseSeqComponent(const pulseseq::Component &comp) const;

    void updatePulseSeqComponent(const pulseseq::Component &comp);

    std::size_t getNumOfPackets() const;
    std::vector<ObservableId> getObservableIds() const;
    const HamiltonianPacket& getPacket(const ObservableId &oid) const { return m_packets.at(oid); }
  private:
    std::map<ObservableId, HamiltonianPacket> m_packets;
  };
} // namespace dnpsoup

#endif

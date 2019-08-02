#ifndef DNPSOUP_SIMULATIONPACKET_H
#define DNPSOUP_SIMULATIONPACKET_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/Observable.h"
#include <vector>
#include <tuple>
#include <memory>   // unique_ptr
#include <map>
#include <utility>


namespace dnpsoup {
  class SimulationPacket {
  public:
    SimulationPacket(std::unique_ptr<InteractionInterface>, 
        const Property &, const Euler<> &);
    SimulationPacket(const SimulationPacket &) = delete;
    SimulationPacket(SimulationPacket &&) noexcept;
    SimulationPacket& operator=(const SimulationPacket &) = delete;
    SimulationPacket& operator=(SimulationPacket &&) noexcept;
    ~SimulationPacket() {}

    // rotate e, then m_e
    SimulationPacket& rotate(const Euler<> &e);
    Euler<> getEuler() const { return m_e; }

    SimulationPacket& setPropertyValue(const ValueName &, double val);
    double getPropertyValue(const ValueName &) const;

    MatrixCxDbl genMatrix() const;
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

    PacketCollection& add(const ObservableId &, SimulationPacket &&);
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

    std::size_t getNumOfPackets() const;
  private:
    std::map<ObservableId, SimulationPacket> m_packets;
  };
} // namespace dnpsoup

#endif

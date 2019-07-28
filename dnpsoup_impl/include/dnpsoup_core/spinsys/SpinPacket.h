#ifndef DNPSOUP_SPINPACKET_H
#define DNPSOUP_SPINPACKET_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/Property.h"
#include <vector>
#include <tuple>
#include <unique_ptr>
#include <utility>


namespace dnpsoup {
  class SpinPacket {
  public:
    SpinPacket(std::unique_ptr<InteractionInterface>, 
        const Property &, const Euler<> &);
    SpinPacket(const SpinPacket &) = delete;
    SpinPacket(SpinPacket &&) noexcept = default;
    SpinPacket& operator=(const SpinPacket &) = delete;
    SpinPacket& operator=(SpinPacket &&) noexcept = default;
    ~SpinPacket() {}

    // rotate e, then m_e
    SpinPacket& rotate(const Euler<> &e);
    Euler<> getEuler() const { return m_e; }

    SpinPacket& setPropertyValue(const ValueName &, double val);
    double getPropertyValue(const ValueName &) const;

    MatrixCxDbl genMatrix() const;
  private:
    std::unique_ptr<InteractionInterface> m_ptr_interface;
    Property m_property;
    Euler<> m_e;
  };

  class SpinPacketCollection {
  public:
    SpinPacketCollection() {}
    SpinPacketCollection(const SpinPacketCollection &) = delete;
    SpinPacketCollection(SpinPacketCollection &&) noexcept = default;
    SpinPacketCollection& operator=(const SpinPacketCollection &) = delete;
    SpinPacketCollection& operator=(SpinPacketCollection &&) noexcept = default;
    ~SpinPacketCollection() {}

    SpinPacketCollection& add(SpinPacket &&);
    SpinPacketCollection& rotate(const Euler<> &e);
    SpinPacketCollection& setPropertyValue(const ValueName &, double val);

    MatrixCxDbl genMatrix() const;

    std::vector<SpinType> getSpinTypes() const;
    std::vector<SpinId> getSpinIds() const;
  private:
    std::vector<SpinPacket> m_packets;
    std::vector<std::pair<SpinId, SpinType>> m_ordered;
  };
} // namespace dnpsoup

#endif

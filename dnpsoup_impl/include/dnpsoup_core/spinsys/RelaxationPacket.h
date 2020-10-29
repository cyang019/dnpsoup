#ifndef DNPSOUP_RELAXATIONPACKET_H
#define DNPSOUP_RELAXATIONPACKET_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/SpinEntity.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/relaxation.h"
#include <vector>
#include <utility>


namespace dnpsoup {
  class RelaxationPacket {
  public:
    RelaxationPacket(const SpinId &, const SpinEntity &,
        std::size_t nbefore, std::size_t nafter);
    RelaxationPacket(const RelaxationPacket &) = default;
    RelaxationPacket(RelaxationPacket &&) noexcept = default;
    RelaxationPacket& operator=(const RelaxationPacket &) = default;
    RelaxationPacket& operator=(RelaxationPacket &&) noexcept = default;
    ~RelaxationPacket() {}

    void setT1(double t1) { m_t1 = t1; }
    void setT2(double t2) { m_t2 = t2; }

    double getT1() const { return m_t1; }
    double getT2() const { return m_t2; }
    SpinId getSpinId() const { return m_id; }

    // rotate from eigen frame to current frame, then [A_{-q}, [A_q, . ] ]
    //MatrixCxDbl genSuperOpT1(const MatrixCxDbl &eigenvec) const;
    //MatrixCxDbl genSuperOpT2(const MatrixCxDbl &eigenvec) const;
    MatrixCxDbl genSuperOpT1() const;
    MatrixCxDbl genSuperOpT2() const;
  private:
    SpinId m_id;
    MatrixCxDbl m_x;
    MatrixCxDbl m_y;
    MatrixCxDbl m_z;
    double m_t1;
    double m_t2;
  };

  class CustomRelaxationPacket {
  public:
    CustomRelaxationPacket(
        const std::vector<std::pair<SpinType, OperatorType>> &ops, double t, 
        double scale=1.0);
    CustomRelaxationPacket(const CustomRelaxationPacket &) = default;
    CustomRelaxationPacket(CustomRelaxationPacket &&) noexcept = default;
    CustomRelaxationPacket& operator=(const CustomRelaxationPacket &) = default;
    CustomRelaxationPacket& operator=(CustomRelaxationPacket &&) noexcept = default;
    ~CustomRelaxationPacket() {}

    void setT(double t) { m_t = t; }
    double getT() const { return m_t; }

    MatrixCxDbl genSuperOp() const;
  private:
    std::vector<std::pair<SpinType, OperatorType>> m_ops;
    MatrixCxDbl m_mat;
    double m_t;
  };

  class RelaxationPacketCollection {
  public:
    MatrixCxDbl genMatrix() const;
    void addRelaxationPacket(const SpinId &sid, const SpinEntity &sinfo,
        std::size_t nbefore, std::size_t nafter);
    void addCustomRelaxationPacket(
        const std::vector<std::pair<SpinType, OperatorType>> &ops,
        double t, double scale=1.0);
  private:
    std::vector<RelaxationPacket> m_rpackets;
    std::vector<CustomRelaxationPacket> m_crpackets;
  };
} // namespace dnpsoup

#endif


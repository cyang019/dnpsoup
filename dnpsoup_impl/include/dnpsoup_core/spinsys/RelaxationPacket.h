#ifndef DNPSOUP_RELAXATIONPACKET_H
#define DNPSOUP_RELAXATIONPACKET_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spinsys/SpinId.h"
#include "dnpsoup_core/spinsys/SpinEntity.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/relaxation.h"


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
    MatrixCxDbl genSuperOpT1(const MatrixCxDbl &eigenvec) const;
    MatrixCxDbl genSuperOpT2(const MatrixCxDbl &eigenvec) const;
  private:
    SpinId m_id;
    MatrixCxDbl m_x;
    MatrixCxDbl m_y;
    MatrixCxDbl m_z;
    double m_t1;
    double m_t2;
  };

} // namespace dnpsoup

#endif


#ifndef DNPSOUP_DIPOLARINTERACTION_H
#define DNPSOUP_DIPOLARINTERACTION_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/rotation/FrameType.h"
#include <type_traits>


namespace dnpsoup {
  /// enable only if T1, T2 are frame types (i.e. RotatingFrame or LabFrame)
  template<typename T1, typename T2>
  class DipolarInteraction
  : public InteractionInterface {
  public:
    DipolarInteraction(double g1, double g2, size_t n1, size_t n2);
    DipolarInteraction(double g1, double g2, size_t n1, size_t n2, 
        size_t nbefore, size_t nbetween, size_t nafter);
    ~DipolarInteraction() {}

    MatrixCxDbl genMatrix(
        const Property &,
        const Euler<ActiveRotation> &) const override;

    MatrixCxDbl genMatrix(
        const Property &,
        const Euler<PassiveRotation> &) const override;

    size_t dimension() const;
  private:
    /// A, B, C, D, E, F alphabets tensor correspondence
    MatrixCxDbl m_a2n2;
    MatrixCxDbl m_a2n1;
    MatrixCxDbl m_a20;
    MatrixCxDbl m_a21;
    MatrixCxDbl m_a22;

    size_t m_n1;
    size_t m_n2;
    size_t m_nbefore;
    size_t m_nbetween;
    size_t m_nafter;

    double m_gamma1;
    double m_gamma2;
  };  // class DipolarInteraction

  // Principles of Nuclear Magnetic Resonance in One and Two Dimensions Page 47
  double calcF20(double theta);
  cxdbl calcF21(double phi, double theta);
  cxdbl calcF2n1(double phi, double theta);
  cxdbl calcF22(double phi, double theta);
  cxdbl calcF2n2(double phi, double theta);

  inline double calcF20(double theta)
  {
    const double ct = std::cos(theta);
    return 1.0 - 3.0 * ct * ct;
  }

  inline cxdbl calcF21(double phi, double theta)
  {
    const double s2t = std::sin(2.0 * theta);
    const double sp = std::sin(phi);
    const double cp = std::cos(phi);
    return 0.5 * s2t * (cp - cxdbl(0,1) * sp);
  }

  inline cxdbl calcF2n1(double phi, double theta)
  {
    const double s2t = std::sin(2.0 * theta);
    const double sp = std::sin(phi);
    const double cp = std::cos(phi);
    return 0.5 * s2t * (cp + cxdbl(0,1) * sp);
  }

  inline cxdbl calcF22(double phi, double theta)
  {
    const double st = std::sin(theta);
    const double s2p = std::sin(2.0*phi);
    const double c2p = std::cos(2.0*phi);
    return st * st * (c2p - cxdbl(0,1) * s2p);
  }

  inline cxdbl calcF2n2(double phi, double theta)
  {
    const double st = std::sin(theta);
    const double s2p = std::sin(2.0*phi);
    const double c2p = std::cos(2.0*phi);
    return st * st * (c2p + cxdbl(0,1) * s2p);
  }
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/DipolarInteractionImpl.hpp"

#endif



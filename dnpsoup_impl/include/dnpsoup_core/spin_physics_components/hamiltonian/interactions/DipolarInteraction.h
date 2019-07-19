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
  template<typename T1, typename T2, typename Enable1 = void, typename Enable2 = void>
  class DipolarInteraction {};

  /// enable only if T1, T2 are frame types (i.e. RotatingFrame or LabFrame)
  template<typename T1, typename T2>
  class DipolarInteraction<T1, T2, 
        typename std::enable_if<is_frame_type<T1>::value>::type, 
        typename std::enable_if<is_frame_type<T2>::value>::type>
  : public InteractionInterface {
  public:
    DipolarInteraction(double g1, double g2, size_t n1, size_t n2);
    DipolarInteraction(double g1, double g2, size_t n1, size_t n2, 
        size_t nbefore, size_t nbetween, size_t nafter);
    ~DipolarInteraction() {}

    // active rotation
    MatrixCxDbl genMatrix(
        const Property *,
        const Euler &) const override;

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
  };  // class Shift
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/DipolarInteractionImpl.hpp"

#endif



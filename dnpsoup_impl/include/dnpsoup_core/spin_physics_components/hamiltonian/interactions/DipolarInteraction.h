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
  // enable only if T1, T2 are frame types
  template<typename T1, typename T2>
  template<typename U1=T1, typename U2=T2,
           std::enable_if_t<is_frame_type<U1>::value, int> = 0,
           std::enable_if_t<is_frame_type<U2>::value, int> = 0>
  class DipolarInteraction : public InteractionInterface {
  public:
    DipolarInteraction(double g1, double g2, size_t n1, size_t n2);
    DipolarInteraction(double g1, double g2, size_t n1, size_t n2, 
        size_t nbefore, size_t nbetween, size_t nafter);
    ~DipolarInteraction() {}

    // active rotation
    virtual matrix::Matrix<cxdbl> genMatrix(
        const PropertyValue *,
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

#include "dnpsoup_core/spin_physics_components/hamiltonian/DipolarInteractionImpl.hpp"

#endif



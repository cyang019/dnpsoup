#ifndef DNPSOUP_CHEMICALSHIFTINTERACTION_H
#define DNPSOUP_CHEMICALSHIFTINTERACTION_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/InteractionInterface.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include "dnpsoup_core/spin_physics_components/rotation/FrameType.h"


namespace dnpsoup {
  template<typename T, typename Enable = void>
  class ChemicalShiftInteraction {};

  /// enable only if T is a frame type (i.e. RotatingFrame or LabFrame)
  template<typename T>
  class ChemicalShiftInteraction<T, typename std::enable_if<is_frame_type<T>::value>::type>
  : public InteractionInterface {
  public:
    ChemicalShiftInteraction(double gamma, size_t n);
    ChemicalShiftInteraction(double gamma, size_t n, size_t nbefore, size_t nafter);
    ~ChemicalShiftInteraction() {}

    // active rotation
    matrix::Matrix<cxdbl> genMatrix(
        const Property *,
        const Euler &) const override;

    size_t dimension() const;
  private:
    matrix::Matrix<cxdbl> m_x;
    matrix::Matrix<cxdbl> m_y;
    matrix::Matrix<cxdbl> m_z;

    size_t m_n;
    size_t m_nbefore;
    size_t m_nafter;

    double m_gamma;
  };  // class Shift
}   // namespace dnpsoup

#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/ChemicalShiftInteractionImpl.hpp"

#endif


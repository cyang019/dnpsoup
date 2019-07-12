#ifndef DNPSOUP_CHEMICALSHIFTINTERACTION_H
#define DNPSOUP_CHEMICALSHIFTINTERACTION_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/spin_physics_components/InteractionInterface.h"


namespace dnpsoup {
  // abstract base class
  class ChemicalShiftInteraction : public InteractionInterface {
  public:
    ChemicalShiftInteraction(size_t n);
    ChemicalShiftInteraction(size_t n, size_t nbefore, size_t nafter);

    // active rotation
    virtual matrix::Matrix<cxdbl> genMatrix(
        const PropertyValueInterface *,
        const Euler &) const override;

    /// @param sxx: in Hz
    /// @param syy: in Hz
    /// @param szz: in Hz
    /// @param beta: angle in rads
    /// @param gamma: angle in rads
    double calcCoeffZ(double sxx, double syy, double szz,
        double beta, double gamma) const; 

    size_t dimension() const;
  private:
    matrix::Matrix<cxdbl> m_iz;

    size_t m_n;
    size_t m_nbefore;
    size_t m_nafter;
  };  // class Shift
}   // namespace dnpsoup

#endif


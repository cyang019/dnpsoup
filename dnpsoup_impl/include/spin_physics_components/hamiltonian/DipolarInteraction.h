#ifndef DNPSOUP_DIPOLARINTERACTION_H
#define DNPSOUP_DIPOLARINTERACTION_H

#include "constants.h"
#include "common.h"
#include "spin_physics_components/hamiltonian/irreducible_tensor_op.h"
#include "spin_physics_components/hamiltonian/spin.h"

namespace dnpsoup {
  template<int l, int m>
  constexpr double calcDipoleF();

  class DipolarInteraction {
  public:
    DipolarInteraction();
    ~DipolarInteraction();

    /// first calculate m_b * g1 * g2/r^3, then use theta to calculate F<l, m>'s., then combine
    /// @param g1: gamma1/2pi
    /// @param g2: gamma2/2pi
    /// @param r: distance in Anstrom (1e-10)
    MatrixCxDbl genOperator(double g1, double g2, double r, double theta, double phi) const;
  private:
    MatrixCxDbl m_A2n2;
    MatrixCxDbl m_A2n1;
    MatrixCxDbl m_A20;
    MatrixCxDbl m_A21;
    MatrixCxDbl m_A22;
    double m_b;   // mu_0 * hbar / 4 * 1e30
  };
} // namespace dnpsoup

#include "spin_physics_components/hamiltonian/DipolarInteraction.hpp"

#endif

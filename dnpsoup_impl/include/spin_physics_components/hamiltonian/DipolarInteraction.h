#ifndef DNPSOUP_DIPOLARINTERACTION_H
#define DNPSOUP_DIPOLARINTERACTION_H

#include "constants.h"
#include "common.h"
#include "spin_physics_components/spin.h"
#include <cmath>


namespace dnpsoup {
  // Abstract base class
  class DipolarInteraction {
  public:
    DipolarInteraction(size_t n1, size_t n2);
    DipolarInteraction(size_t n1, size_t n2, 
        size_t nbefore, size_t nbetween, size_t nafter);
    DipolarInteraction(const DipolarInteraction &) = default;
    DipolarInteraction(DipolarInteraction &&) noexcept = default;
    DipolarInteraction& operator=(const DipolarInteraction &) = default;
    DipolarInteraction& operator=(DipolarInteraction &&) noexcept = default;
    ~DipolarInteraction();

    /// first calculate m_b * g1 * g2/r^3, then use theta to calculate F<l, m>'s., then combine
    /// @param g1: gamma1/2pi
    /// @param g2: gamma2/2pi
    /// @param r: distance in Anstrom (1e-10)
    MatrixCxDbl genMatrix(double g1, double g2, double r, double theta, double phi) const = 0;

    double calcCoeff(int j, int m) const;
  private:
    /// A, B, C, D, E, F alphabets
    MatrixCxDbl m_a2n2;
    MatrixCxDbl m_a2n1;
    MatrixCxDbl m_a20;
    MatrixCxDbl m_a21;
    MatrixCxDbl m_a22;

    double m_b;   // mu_0 * hbar / 4 * 1e30

    size_t m_n1;
    size_t m_n2;
    size_t m_nbefore;
    size_t m_nbetween;
    size_t m_nafter;
  };
} // namespace dnpsoup

#endif

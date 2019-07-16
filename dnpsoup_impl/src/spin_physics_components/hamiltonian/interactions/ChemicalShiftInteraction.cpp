#include "dnpsoup_core/spin_physics_components/hamiltonian/interactions/ChemicalShiftInteraction.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include <cmath>
#include <complex>

using std::cos;
using std::sin;
using cxdbl = std::complex<double>;

namespace dnpsoup {
  ChemicalShiftInteraction::ChemicalShiftInteraction(size_t n)
    : InteractionInterface(), m_n(n), m_nbefore(0), m_nafter(0)
  { m_iz = spin<X>(n); }

  ChemicalShiftInteraction::ChemicalShiftInteraction(size_t n, size_t nbefore, size_t nafter)
    : InteractionInterface(), m_n(n), m_nbefore(nbefore), m_nafter(nafter)
  { m_iz = kron(kron(identity<cxdbl>(nbefore), spin<X>(n)), identity<cxdbl>(nafter)); }

  MatrixCxDbl ChemicalShiftInteraction::genMatrix(
      const PropertyValue *ptr_csa,
      const Euler &e) const
  {
    const double &sxx = ptr_csa->get("sxx");
    const double &syy = ptr_csa->get("syy");
    const double &szz = ptr_csa->get("szz");
    const double &gyro = ptr_csa->get("gyromagnetic ratio");
    const double &bz = ptr_csa->get("Bz");

    const double sb = sin(e.beta());
    const double cb = cos(e.beta());
    const double sg = sin(e.gamma());
    const double cg = cos(e.gamma());
    double coeff3 = szz * cb * cb + sb * sb * (sxx * cg * cg + syy * sg * sg);
    const double coeff_field = gyro * bz;
    MatrixCxDbl res = (coeff_field * coeff3) * m_iz;
    return res;
  }

  size_t ChemicalShiftInteraction::dimension() const
  { if (m_n == 0 && m_nbefore == 0 && m_nafter == 0) 
      return 0;
    size_t dim_before = m_nbefore > 0 ? m_nbefore : 1;
    size_t dim_after = m_nafter > 0 ? m_nafter : 1;
    size_t m_mat = m_n > 0 ? m_n : 1;
    return dim_before * m_mat * dim_after;
  }
} // namespace dnpsoup


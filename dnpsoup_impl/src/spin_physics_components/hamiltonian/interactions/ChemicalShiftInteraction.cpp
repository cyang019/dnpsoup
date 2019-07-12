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
    : m_n(n), m_nbefore(0), m_nafter(0)
  { m_iz = spin<X>(n); }

  ChemicalShiftInteraction::ChemicalShiftInteraction(size_t n, size_t nbefore, size_t nafter)
    : m_n(n), m_nbefore(nbefore), m_nafter(nafter)
  { m_iz = kron(kron(identity<cxdbl>(nbefore), spin<X>(n)), identity<cxdbl>(nafter)); }

  MatrixCxDbl ChemicalShiftInteraction::genMatrix(
      const PropertyValueInterface *ptr_csa,
      const Euler &e) const
  {
    std::vector<double> vals = ptr_csa->get();
    double &gyro = vals[0];
    double &bz = vals[1];
    double &sxx = vals[2];
    double &syy = vals[3];
    double &szz = vals[4];

    const double sb = sin(e.beta());
    const double cb = cos(e.beta());
    const double sg = sin(e.gamma());
    const double cg = cos(e.gamma());
    double coeff3 = szz * cb * cb + sb * sb * (sxx * cg * cg + syy * sg * sg);
    const double coeff_field = gyro * bz;
    MatrixCxDbl res = (coeff_field * coeff3) * m_iz;
    return res;
  }
} // namespace dnpsoup


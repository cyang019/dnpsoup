#include "spin_physics_components/hamiltonian/Shift.h"
#include "common.h"
#include "errors.h"
#include <cmath>
#include <complex>

using std::cos;
using std::sin;
using cxdbl = std::complex<double>;

namespace dnpsoup {
  Shift::Shift(size_t n)
    : m_ix(n, n), m_iy(n, n), m_iz(n, n), m_n(n), m_nbefore(0), m_nafter(0)
  {
    m_ix = spin<X>(m_n);
    m_iy = spin<Y>(m_n);
    m_iz = spin<Z>(m_n);
  }

  Shift::Shift(size_t n, size_t nbefore, size_t nafter)
    : m_n(n), m_nbefore(nbefore), m_nafter(nafter)
  {
    m_ix = kroneckerProduct(kroneckerProduct(matrix::identity<cxdbl>(nbefore), spin<X>(n)),
                       matrix::identity<cxdbl>(nafter));
    m_iy = kroneckerProduct(kroneckerProduct(matrix::identity<cxdbl>(nbefore), spin<Y>(n)),
                       matrix::identity<cxdbl>(nafter));
    m_iz = kroneckerProduct(kroneckerProduct(matrix::identity<cxdbl>(nbefore), spin<Z>(n)),
                       matrix::identity<cxdbl>(nafter));
  }

  // not implemented yet. placeholder
  matrix::Matrix<cxdbl> Shift::genMatrix(double gyro, double bz, 
        double sxx, double syy, double szz,
        double alpha, double beta, double gamma) const
  {
    const double sa = sin(alpha);
    const double ca = cos(alpha);
    const double sb = sin(beta);
    const double s2b = sin(2.0*beta);
    const double cb = cos(beta);
    const double sg = sin(gamma);
    const double s2g = sin(2.0*gamma);
    const double cg = cos(gamma);
    double coeff1 = (sxx - syy) * 0.5 * sa * sb * s2g + ca * s2b * (szz - sxx * (cg * cg) -syy * sg * sg);
    double coeff2 = (syy - sxx) * 0.5 * ca * sb * s2g + sa * s2b * (szz - sxx * (cg * cg) -syy * sg * sg);
    double coeff3 = szz * cb * cb + sb * sb * (sxx * cg * cg + syy * sg * sg);
    const double coeff_field = gyro * bz;
    matrix::Matrix<cxdbl> res = coeff_field * (coeff1 * m_ix + coeff2 * m_iy + coeff3 * m_iz);
    return res;
  }

  double Shift::calcCoeffX(double sxx, double syy, double szz,
      double alpha, double beta, double gamma) const
  {
    const double sa = sin(alpha);
    const double ca = cos(alpha);
    const double sb = sin(beta);
    const double s2b = sin(2.0*beta);
    const double s2g = sin(2.0*gamma);
    const double cg = cos(gamma);
    const double sg = sin(gamma);

    double res = sa * sb * s2g * 0.5 * (sxx - syy)
      + ca * 0.5 * s2b * (szz - sxx * cg * cg - syy * sg * sg);
    return res;
  }

  double Shift::calcCoeffY(double sxx, double syy, double szz,
      double alpha, double beta, double gamma) const
  {
    const double sa = sin(alpha);
    const double ca = cos(alpha);
    const double sb = sin(beta);
    const double s2b = sin(2.0*beta);
    const double s2g = sin(2.0*gamma);
    const double cg = cos(gamma);
    const double sg = sin(gamma);

    double res = ca * sb * s2g * 0.5 * (syy - sxx)
      + sa * 0.5 * s2b * (szz - sxx * cg * cg - syy * sg * sg);
    return res;
  }

  double Shift::calcCoeffZ(double sxx, double syy, double szz,
      double beta, double gamma) const 
  {
    const double sb = sin(beta);
    const double cb = cos(beta);
    const double cg = cos(gamma);
    const double sg = sin(gamma);
    double res = szz * cb * cb + sb * sb * (sxx * cg * cg + syy * sg * sg);
    return res;
  }

  size_t Shift::dimension() const
  { if (m_n == 0 && m_nbefore == 0 && m_nafter == 0) 
      return 0;
    size_t dim_before = m_nbefore > 0 ? m_nbefore : 1;
    size_t dim_after = m_nafter > 0 ? m_nafter : 1;
    size_t m_mat = m_n > 0 ? m_n : 1;
    return dim_before * m_mat * dim_after;
  }
} // namespace dnpsoup

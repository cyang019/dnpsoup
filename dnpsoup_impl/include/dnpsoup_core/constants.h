#ifndef DNPSOUP_CONSTANTS_H
#define DNPSOUP_CONSTANTS_H

#include <complex>
#include <limits>
#include <cmath>

namespace dnpsoup {
  using cxdbl = std::complex<double>;
  constexpr cxdbl img = cxdbl(0,1);
  /// pi from std::M_PI
  constexpr double pi = 3.14159265358979323846264338327950288;
  constexpr double mu0 = 4 * pi * 1.0e-7;

  constexpr inf = std::numeric_limits<double>::infinity();

  /// From NIST
  constexpr double kb = 1.38064852e-23;   ///< \f$m^2 kg s^{-2} K^{-1} \f$
  constexpr double kb_inv = 7.242973e22;  ///< inverse of kb

  constexpr double hbar = 1.0545718e-34; ///< \f$ m^2 kg /s \f$
  constexpr double h = 6.626070040e-34; ///< J s  

  constexpr double mu_b = 9.274009994e-24;  ///< \f$J T^{-1}\f$
  constexpr double mu_e = -9.284764620e-24;   ///< \f$J T^{-1}\f$
  constexpr double hbar_over_kb = 7.6382351e-12;

  // =======================================================
  // gyromagnetic ratios
  /// From Nist
  constexpr double beta_e = -mu_b/h; ///< in Hz/T
  constexpr double gamma_H1 = 42.57747892e6;  // in Hz/T, actually gamma_H1/2pi
  constexpr double gamma_N14 = 3.077e6; // in Hz/T, actually gamma_N14/2pi
  /// From wikipedia
  constexpr double gamma_D2 = 6.536e6;    // wiki, actually gamma_D2/2pi
  constexpr double gamma_O17 = -5.772e6;  // wiki, actually gamma_O17/2pi
  constexpr double gamma_C13 = 10.705e6;  // in Hz/T, actually gamma_C13/2pi
  constexpr double gamma_N15 = -4.316e6;  // in Hz/T, actually gamma_N15/2pi
  // =======================================================
  
  //acos(1.0/sqrt(3.0));  ///< in Rad
  constexpr double magic_angle = 0.9553166181245092;  ///< in Radius
} // namespace dnpsoup

#endif

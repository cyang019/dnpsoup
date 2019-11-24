#ifndef DNPSOUP_EMRADIATION_H
#define DNPSOUP_EMRADIATION_H

#include <iostream>

namespace dnpsoup {
  namespace pulseseq{
  struct EMRadiation {
      EMRadiation() : freq(0.0), phase(0.0), offset(0.0) {}
      EMRadiation(double f, double p)
        : freq(f), phase(p) {}
      EMRadiation(double f, double p, double offset)
        : freq(f), phase(p), offset(offset) {}
      EMRadiation(const EMRadiation &) = default;
      EMRadiation(EMRadiation &&) noexcept = default;
      EMRadiation& operator=(const EMRadiation &) = default;
      EMRadiation& operator=(EMRadiation &&) noexcept = default;
      ~EMRadiation() {}

      EMRadiation& reset();

      double freq;    ///< in Hz  (gamma B1)
      double phase;   ///< in rad (phase)
      double offset;  ///< in Hz  (carrier)
    };  // class EMRadiation

    std::istream& operator>>(std::istream&, EMRadiation &);
    std::ostream& operator<<(std::ostream&, const EMRadiation &);
  }
} // namespace dnpsoup

#endif

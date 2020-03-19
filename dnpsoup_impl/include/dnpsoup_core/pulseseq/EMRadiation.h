#ifndef DNPSOUP_EMRADIATION_H
#define DNPSOUP_EMRADIATION_H

#include <iostream>
#include "dnpsoup_core/common.h"

namespace dnpsoup {
  namespace pulseseq{
  struct EMRadiation {
      EMRadiation() 
        : freq(0.0), phase(0.0), offset(0.0), id_(-1) {}
      EMRadiation(double f, double p)
        : freq(f), phase(p), id_(-1) {}
      EMRadiation(double f, double p, double offset)
        : freq(f), phase(p), offset(offset), id_(-1) {}
      EMRadiation(const EMRadiation &) = default;
      EMRadiation(EMRadiation &&) noexcept = default;
      EMRadiation& operator=(const EMRadiation &) = default;
      EMRadiation& operator=(EMRadiation &&) noexcept = default;
      ~EMRadiation() {}

      EMRadiation& reset();

      double freq;    ///< in Hz  (gamma B1)
      double phase;   ///< in rad (phase)
      double offset;  ///< in Hz  (carrier)

      void setId(int value) { id_ = value; }
      int getId() const { return id_; }
    private:
      int id_;   ///< to lookup in cache
    };  // class EMRadiation

    std::istream& operator>>(std::istream&, EMRadiation &);
    std::ostream& operator<<(std::ostream&, const EMRadiation &);

  } // namespace pulseseq
  inline bool sameValue(const pulseseq::EMRadiation &emr1, const pulseseq::EMRadiation &emr2, double eps)
  {
    auto res = approxEqual<double>(emr1.freq, emr2.freq, eps, eps) 
      && approxEqual<double>(emr1.offset, emr2.offset, eps, eps)
      && approxEqual<double>(emr1.phase, emr2.phase, eps, eps);
    return res;
  }

  inline bool sameId(
      const pulseseq::EMRadiation &emr1,
      const pulseseq::EMRadiation &emr2)
  {
    return emr1.getId() == emr2.getId();
  }
  
} // namespace dnpsoup

#endif

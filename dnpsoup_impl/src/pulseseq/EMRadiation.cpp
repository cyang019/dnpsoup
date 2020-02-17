#include "dnpsoup_core/pulseseq/EMRadiation.h"
#include "dnpsoup_core/common.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace std;


namespace dnpsoup {
namespace pulseseq {
  EMRadiation& EMRadiation::reset()
  {
    freq = 0.0;
    phase = 0.0;
    offset = 0.0;
    return *this;
  }

  std::istream& operator>>(std::istream &is, EMRadiation &p)
  {
    p.reset();
    string word;
    while(is >> word){
      if(word == "freq" || word == "Freq" || word == "FREQ"){
        is >> p.freq;
      } else if(word == "phase" || word == "Phase" || word == "PHASE"){
        is >> p.phase;
      } else if(word == "offset" || word == "Offset" || word == "OFFSET"){
        is >> p.offset;
      }
    }
    return is;
  }

  std::ostream& operator<<(std::ostream &os, const EMRadiation &p)
  {
    std::streamsize ss = os.precision();
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "freq " << p.freq << " Hz, phase " << p.phase 
       << " rad, offset " << p.offset << " Hz";
    os << std::setprecision(ss);
    return os;
  }
} // namespace pulseseq
  bool sameValue(const pulseseq::EMRadiation &emr1, const pulseseq::EMRadiation &emr2, double eps)
  {
    auto res = approxEqual<double>(emr1.freq, emr2.freq, eps) 
      && approxEqual<double>(emr1.offset, emr2.offset, eps)
      && approxEqual<double>(emr1.phase, emr2.phase, eps);
    return res;
  }
} // namespace dnpsoup

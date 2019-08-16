#include "dnpsoup_core/pulseseq/PulsePacket.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>

using namespace std;


namespace dnpsoup {
  PulsePacket& PulsePacket::reset()
  {
    freq = 0.0;
    phase = 0.0;
    offset = 0.0;
    return *this;
  }

  std::istream& operator>>(std::istream &is, PulsePacket &p)
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

  std::ostream& operator<<(std::ostream &os, const PulsePacket &p)
  {
    std::streamsize ss = os.precision();
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "freq " << p.freq << " Hz, phase " << p.phase 
       << " rad, offset " << p.offset << " Hz";
    os << std::setprecision(ss);
    return os;
  }
} // namespace dnpsoup

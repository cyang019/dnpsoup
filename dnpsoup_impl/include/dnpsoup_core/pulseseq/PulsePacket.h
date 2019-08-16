#ifndef DNPSOUP_PULSEPACKET_H
#define DNPSOUP_PULSEPACKET_H

#include <iostream>

namespace dnpsoup {
  struct PulsePacket {
    PulsePacket() : freq(0.0), phase(0.0), offset(0.0) {}
    PulsePacket(double f, double p)
      : freq(f), phase(p) {}
    PulsePacket(double f, double p, double offset)
      : freq(f), phase(p), offset(offset) {}
    PulsePacket(const PulsePacket &) = default;
    PulsePacket(PulsePacket &&) noexcept = default;
    PulsePacket& operator=(const PulsePacket &) = default;
    PulsePacket& operator=(PulsePacket &&) noexcept = default;
    ~PulsePacket() {}

    PulsePacket& reset();

    double freq;    // in Hz
    double phase;       // in rad
    double offset;      // in Hz
  };  // class PulsePacket

  std::istream& operator>>(std::istream&, PulsePacket &);
  std::ostream& operator<<(std::ostream&, const PulsePacket &);
} // namespace dnpsoup

#endif

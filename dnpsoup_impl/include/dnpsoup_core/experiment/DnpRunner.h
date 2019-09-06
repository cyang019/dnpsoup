#ifndef DNPSOUP_DNPRUNNER_H
#define DNPSOUP_DNPRUNNER_H

#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/pulseseq/pulse_sequence.h"
#include "dnpsoup_core/spinsys/SpinSys.h"
#include "dnpsoup_core/hardware.h"
#include <vector>
#include <sstream>
#include <string>
#include <cstdint>

namespace dnpsoup {
  class DnpRunner {
    double calcEnhancement(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler);
  };
} // namespace dnpsoup

#endif

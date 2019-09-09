#ifndef DNPSOUP_DNPRUNNER_H
#define DNPSOUP_DNPRUNNER_H

#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/pulseseq/pulse_sequence.h"
#include "dnpsoup_core/spinsys/SpinSys.h"
#include "dnpsoup_core/experiment/hardware.h"
#include <vector>
#include <sstream>
#include <string>
#include <cstdint>
#include <utility>

namespace dnpsoup {
  class DnpRunner {
  public:
    double calcIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler) const;

    double calcPowderIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers) const;

    MatrixCxDbl evolve(
        const MatrixCxDbl &rho_prev, 
        const MatrixCxDbl &hamiltonian,
        const RelaxationPacket &rpacket,
        double dt) const;
  };

} // namespace dnpsoup

#endif

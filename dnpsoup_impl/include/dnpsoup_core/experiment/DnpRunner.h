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
    std::vector<std::vector<double>> calcEigenValues(
        const Magnet &m,
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const Euler<> &spin_sys_euler
        ) const;

    double calcIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler) const;

    std::vector<std::pair<double, double>> calcBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler) const;
    //std::vector<double> calcIntensityBuildUp(
    //    const Magnet &m, 
    //    const Gyrotron &g,
    //    const Probe &p,
    //    const SpinSys &spin_sys,
    //    const std::string &pulse_seq_str,
    //    const SpinType &acq_spin,
    //    const Euler<> &spin_sys_euler) const;

    std::vector<double> calcFieldProfile(
        const std::vector<Magnet> &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1) const;

    std::vector<double> calcFieldProfile(
        const Magnet &fields, 
        const std::vector<Gyrotron> &emrs,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1) const;

    double calcPowderIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1) const;

    MatrixCxDbl evolve(
        const MatrixCxDbl &rho_prev, 
        const MatrixCxDbl &hamiltonian,
        const std::vector<RelaxationPacket> &rpackets,
        double dt, double temperature) const;
  };

} // namespace dnpsoup

#endif

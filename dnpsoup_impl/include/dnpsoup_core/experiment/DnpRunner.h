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
  namespace DnpRunner {
    std::vector<std::vector<double>> calcEigenValues(
        const Magnet &m,
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const Euler<> &spin_sys_euler
        );

    double calcIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler);

    std::vector<std::pair<double, double>> calcBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler);

    std::vector<std::pair<double, double>> 
      calcBuildUpEnhancement(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler);

    std::vector<std::pair<double, double>> calcPowderBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1);

    std::vector<std::pair<double, double>> calcPowderBuildUpEnhancement(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1);

    std::vector<std::pair<double, double>> calcFieldProfile(
        const std::vector<Magnet> &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1);

    std::vector<std::pair<double, double>> calcFieldProfile(
        const Magnet &fields, 
        const std::vector<Gyrotron> &emrs,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1);

    double calcPowderIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores=1);

    MatrixCxDbl propagate(
        const MatrixCxDbl &rho_prev_super, 
        const MatrixCxDbl &hamiltonian,
        const MatrixCxDbl &hamiltonian_lab,
        const MatrixCxDbl &rotate_mat_super,
        const MatrixCxDbl &rotate_mat_super_inv,
        const std::vector<RelaxationPacket> &rpackets,
        double dt, std::uint64_t cnt, double temperature);

    MatrixCxDbl propagate(
        const MatrixCxDbl &rho_prev_super, 
        const PacketCollection &packets,
        const MatrixCxDbl &hamiltonian_offset,
        const std::vector<RelaxationPacket> &rpackets,
        const Gyrotron &g,
        const Euler<> &euler,
        double dt,
        std::uint64_t cnt,
        double temperature
        );
  }

} // namespace dnpsoup

#endif

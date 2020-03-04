#ifndef DNPSOUP_SCAN_H
#define DNPSOUP_SCAN_H

#include "dnpsoup_core/experiment/DnpRunner.h"
#include <vector>
#include <tuple>
#include <cstdint>
#include <string>


namespace dnpsoup {
  enum class ScanType : int
  {
    //FieldType = 0,
    //EmOffsetType = 1,

    GammaB1Type = 10,
    PhaseType = 11,
    LengthType = 12,

    EulerType = 100,

    T1Type = 200,
    T2Type = 201
  };

  /// x, y pairs
  using ScanResults1D = std::vector<std::pair<double, double>>;

  /// x, y, values
  using ScanResults2D = std::vector<std::tuple<double, double, double>>;

  struct Parameters {
    Parameters(const Magnet &m, const Gyrotron &g, const Probe &p,
        const SpinSys &spin_sys, const PulseSequence &seq,
        const SpinType &acq_spin, const std::vector<Euler<>> &spin_sys_eulers);
    Parameters(const Parameters &) = default;
    Parameters(Parameters &&) noexcept = default;
    Parameters& operator=(const Parameters &) = default;
    Parameters& operator=(Parameters &&) noexcept = default;
    ~Parameters() {}

    Magnet magnet;
    Gyrotron gyrotron;
    Probe probe;
    SpinSys spin_sys;
    PulseSequence seq;
    SpinType acq_spin;
    std::vector<Euler<>> spin_sys_eulers;
  };

  ScanResults1D scan1dEmrFreq(
      const Parameters &params,
      const std::string &name,
      const SpinType t,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores=1);

  ScanResults1D scan1dEmrPhase(
      const Parameters &params,
      const std::string &name,
      const SpinType t,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores=1);

  ScanResults1D scan1dEmrLength(
      const Parameters &params,
      const std::string &name,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores=1);

  std::vector<double> scan1d(
      const Parameters &params,
      const std::vector<PulseSequence> &seqs,
      int ncores=1);

  std::vector<double> populateValues(double val_beg, double val_end, std::uint64_t cnt);
}


#endif

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

  ScanResults1D scan1d(
      const Parameters &params,
      const ScanType &scan_t,
      const std::string &scan_name,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores=1);

  ScanResults2D scan2d(
      const Parameters &params,
      const ScanType &scan_t_1,
      const std::string &scan_name_1,
      double value_begin_1,
      double value_end_1,
      std::uint64_t cnt1,
      const ScanType &scan_t_2,
      const std::string &scan_name_2,
      double value_begin_2,
      double value_end_2,
      std::uint64_t cnt2,
      int ncores=1);

  std::vector<double> populateValues(double val_beg, double val_end, std::uint64_t cnt);
}


#endif

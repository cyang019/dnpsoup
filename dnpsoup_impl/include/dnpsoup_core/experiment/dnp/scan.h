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
    DefaultType = 0,

    EmrGammaB1Type = 10,
    EmrPhaseType = 11,
    EmrLengthType = 12,
  };

  ScanType getScanType(const std::string &type_str);

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

  struct Range {
    Range() : beg(0.0), end(0.0), cnt(1u) {}
    Range(double beg, double end, std::uint64_t cnt)
      : beg(beg), end(end), cnt(cnt)
    {}
    Range(const Range &) = default;
    Range(Range &&) noexcept = default;
    Range& operator=(const Range &) = default;
    Range& operator=(Range &&) noexcept = default;
    ~Range() {}

    std::vector<double> values() const;

    double beg;
    double end;
    std::uint64_t cnt;
  };

  class Selector {
  public:
    Selector();
    Selector(const ScanType &scan_t, 
        const std::string &name, 
        const SpinType &spin_t);

    Selector(const ScanType &scan_t, const std::string &name); 
    Selector(const Selector &) = default;
    Selector(Selector &&) noexcept = default;
    Selector& operator=(const Selector &) = default;
    Selector& operator=(Selector &&) noexcept = default;

    PulseSequence modify(const PulseSequence &seq, double value) const;
    PulseSequence modify(PulseSequence &&seq, double value) const;
  private:
    ScanType scan_t_;
    std::string name_;
    SpinType spin_t_;
  };

  ScanResults1D scan1d(
      const Parameters &params,
      const Selector &selector,
      const Range &range,
      int ncores=1);

  ScanResults2D scan2d(
      const Parameters &params,
      const Selector &selector1,
      const Range &range1,
      const Selector &selector2,
      const Range &range2,
      int ncores=1);

  std::vector<double> populateValues(double val_beg, double val_end, std::uint64_t cnt);
}


#endif

#ifndef DNPSOUP_SCAN_H
#define DNPSOUP_SCAN_H

#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/common.h"
#include "json.hpp"
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

  union ValueContainer 
  {
    double d_value;
    std::uint64_t sz_value;
  };

  class ScanValueType {
    friend ScanValueType operator+(const ScanValueType &, const ScanValueType &);
    friend bool operator<(const ScanValueType &, const ScanValueType &);
  public:
    ScanValueType();
    ScanValueType(double);
    ScanValueType(std::uint64_t);
    ScanValueType(const ScanValueType &);
    ScanValueType(ScanValueType &&) noexcept;
    ScanValueType& operator=(const ScanValueType &);
    ScanValueType& operator=(ScanValueType &&) noexcept;
    ScanValueType& operator=(const double &);
    ScanValueType& operator=(double &&) noexcept;
    ScanValueType& operator=(const std::uint64_t &);
    ScanValueType& operator=(std::uint64_t &&) noexcept;
    ScanValueType& operator+=(const ScanValueType &);
    ScanValueType& operator+=(const double &);
    ScanValueType& operator+=(double &&) noexcept;
    ScanValueType& operator+=(const std::uint64_t &);
    ScanValueType& operator+=(std::uint64_t &&) noexcept;
    ~ScanValueType() {}

    double getRealValue() const;
    std::uint64_t getSizeValue() const;

    ScanValueType& setRealValue(const double &val);
    ScanValueType& setIntValue(const std::uint64_t &val);

    bool isReal() const { return is_real_; }
  private:
    ValueContainer value_;
    bool is_real_;
  };

  ScanValueType operator+(const ScanValueType &lhs, const ScanValueType &rhs);
  bool operator<(const ScanValueType &, const ScanValueType &);

  class Range {
  public:
    Range();
    Range(double start, double stop, std::uint64_t cnt);
    Range(double start, double stop, double step);
    Range(std::uint64_t start, std::uint64_t stop, std::uint64_t step);
    Range(const Range &) = default;
    Range(Range &&) noexcept = default;
    Range& operator=(const Range &) = default;
    Range& operator=(Range &&) noexcept = default;
    ~Range() {}

    std::vector<ScanValueType> values() const;
  private:
    ScanValueType start_;
    ScanValueType stop_;
    ScanValueType step_;
  };

  template<typename T>
  Range genRangeFromJs(const json &js)
  {
    if(js.find("begin") == js.end()
        || js.find("end") == js.end()
        || js.find("step") == js.end()){
      throw PropertyNameNotFound("A Range needs 'begin', 'end' and 'step' properties.");
    }
    if constexpr(std::is_same<T, double>::value){
      auto start = js["begin"].get<double>();
      auto stop = js["end"].get<double>();
      auto step = js["step"].get<double>();
      return Range(start, stop, step);
    }
    else if constexpr(std::is_same<T, std::uint64_t>::value){
      auto start = js["begin"].get<std::uint64_t>();
      auto stop = js["end"].get<std::uint64_t>();
      auto step = js["step"].get<std::uint64_t>();
      return Range(start, stop, step);
    }
    else {
      throw NotImplementedError("Unknown Range Type.");
    }
  }

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

    PulseSequence modify(const PulseSequence &seq, const ScanValueType &value) const;
    PulseSequence modify(PulseSequence &&seq, const ScanValueType &value) const;
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
  std::vector<std::uint64_t> populateValues(std::uint64_t val_beg, std::uint64_t val_end, std::uint64_t cnt);
}


#endif

#include "dnpsoup_core/experiment/dnp/scan.h"
#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"
#include <cmath>
#include <limits>

using namespace std;


namespace dnpsoup {
  ScanType getScanType(const std::string &type_str)
  {
    if(type_str == "GammaB1"
        || type_str == "GammaB1Type"
        || type_str == "EmrGammaB1"
        || type_str == "EmrGammaB1Type"
        || type_str == "EMRGammaB1"
        || type_str == "EMRGammaB1Type"
        || type_str == "gamma-b1"
        || type_str == "gamma-b1-type"
        || type_str == "emr-gamma-b1"
        || type_str == "emr-gamma-b1-type"){
      return ScanType::EmrGammaB1Type;
    }
    else if(type_str == "EmrPhase"
        || type_str == "EmrPhaseType"
        || type_str == "EMRPhase"
        || type_str == "EMRPhaseType"
        || type_str == "emr-phase"
        || type_str == "emr-phase-type"){
      return ScanType::EmrPhaseType;
    }
    else if(type_str == "EmrLength"
        || type_str == "EmrLengthType"
        || type_str == "EMRLength"
        || type_str == "EMRLengthType"
        || type_str == "emr-length"
        || type_str == "emr-length-type"){
      return ScanType::EmrLengthType;
    }
    else{
      const string err_msg = type_str + " unknown for scan1d type.";
      throw NotImplementedError(err_msg);
    }
    return ScanType::DefaultType;
  }

  // class ScanValueType
  ScanValueType::ScanValueType()
    : value_({0.0}), is_real_(true)
  {}

  ScanValueType::ScanValueType(double val)
    : value_({val}), is_real_(true)
  {}

  ScanValueType::ScanValueType(std::uint64_t val)
    : is_real_(false)
  {
    value_.sz_value = val;
  }

  ScanValueType::ScanValueType(const ScanValueType &rhs)
    : value_(rhs.value_), is_real_(rhs.is_real_)
  {}

  ScanValueType::ScanValueType(ScanValueType &&rhs) noexcept
    : value_(std::move(rhs.value_)), is_real_(rhs.is_real_)
  {}

  ScanValueType& ScanValueType::operator=(const ScanValueType &rhs)
  {
    value_ = rhs.value_;
    is_real_ = rhs.is_real_;
    return *this;
  }

  ScanValueType& ScanValueType::operator=(ScanValueType &&rhs) noexcept
  {
    value_ = std::move(rhs.value_);
    is_real_ = std::move(rhs.is_real_);
    return *this;
  }

  ScanValueType& ScanValueType::operator=(const double &rhs)
  {
    value_.d_value = rhs;
    is_real_ = true;
    return *this;
  }
  ScanValueType& ScanValueType::operator=(double && rhs) noexcept
  {
    value_.d_value = std::move(rhs);
    is_real_ = true;
    return *this;
  }
  ScanValueType& ScanValueType::operator=(const std::uint64_t &rhs)
  {
    value_.sz_value = rhs;
    is_real_ = false;
    return *this;
  }

  ScanValueType& ScanValueType::operator=(std::uint64_t && rhs) noexcept
  {
    value_.sz_value = std::move(rhs);
    is_real_ = false;
    return *this;
  }

  ScanValueType& ScanValueType::operator+=(const double &rhs)
  {
    value_.d_value *= is_real_;
    value_.d_value += rhs;
    is_real_ = true;
    return *this;
  }

  ScanValueType& ScanValueType::operator+=(double && rhs) noexcept
  {
    value_.d_value *= is_real_;
    value_.d_value += std::move(rhs);
    is_real_ = true;
    return *this;
  }

  ScanValueType& ScanValueType::operator+=(const std::uint64_t &rhs)
  {
    value_.sz_value *= (!is_real_);
    value_.sz_value += rhs;
    is_real_ = false;
    return *this;
  }

  ScanValueType& ScanValueType::operator+=(std::uint64_t && rhs) noexcept
  {
    value_.sz_value *= (!is_real_);
    value_.sz_value += std::move(rhs);
    is_real_ = false;
    return *this;
  }

  double ScanValueType::getRealValue() const
  {
    return (is_real_) * value_.d_value;
  }

  std::uint64_t ScanValueType::getSizeValue() const
  {
    return (!is_real_) * value_.sz_value;
  }

  ScanValueType& ScanValueType::setRealValue(const double &val)
  {
    value_.d_value = val;
    is_real_ = true;
    return *this;
  }

  ScanValueType& ScanValueType::setIntValue(const std::uint64_t &val)
  {
    value_.sz_value = val;
    is_real_ = false;
    return *this;
  }
  // class ScanValueType
  
  // class Range
  Range::Range()
    : start_(0.0), stop_(0.0), step_(1.0)
  {}

  Range::Range(double start, double stop, std::uint64_t cnt)
    : start_(start), stop_(stop)
  {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double diff = stop - start;
    if(cnt <= 1) step_ = diff;
    else {
      step_ = diff / static_cast<double>(cnt - 1);
      step_ += 1.0 * (::dnpsoup::approxEqual(step_.getRealValue(), 0.0, eps, eps));
    }
  }

  Range::Range(double start, double stop, double step)
    : start_(start), stop_(stop), step_(step)
  {
      step_ += 1.0 * (::dnpsoup::approxEqual(step, 0.0, eps, eps));
  }

  Range::Range(std::uint64_t start, std::uint64_t stop, std::uint64_t step)
    : start_(start), stop_(stop), step_(step)
  {
    step_ += static_cast<std::uint64_t>(step_.getSizeValue() == 0u);
  }

  vector<double> Range::values() const
  {
    vector<double> results;
    double val = start_.getRealValue();
    const double inc = step_.getRealValue();
    while(val < stop_.getRealValue()){
      results.push_back(val);
      val += inc;
    }
    return results;
  }

  std::vector<std::uint64_t> Range::sz_values() const
  {
    vector<std::uint64_t> results;
    std::uint64_t val = start_.getSizeValue();
    const std::uint64_t inc = step_.getSizeValue();
    while(val < stop_.getSizeValue()){
      results.push_back(val);
      val += inc;
    }
    return results;
  }
  // class Range

  Selector::Selector()
    : scan_t_(ScanType::DefaultType), name_(""), spin_t_(SpinType::Null)
  {}

  Selector::Selector(const ScanType &scan_t, const std::string &name,
      const SpinType &spin_t)
    : scan_t_(scan_t), name_(name), spin_t_(spin_t)
  {}

  Selector::Selector(const ScanType &scan_t, const std::string &name)
    : scan_t_(scan_t), name_(name), spin_t_(SpinType::Null)
  {}

  PulseSequence Selector::modify(
      const PulseSequence &seq,
      const ScanValueType &value) const
  {
    PulseSequence result = seq;
    switch(scan_t_){
      case ScanType::EmrGammaB1Type:
        {
          result.setEmrFreq(name_, spin_t_, value.getRealValue());
          return result;
        }
        break;
      case ScanType::EmrPhaseType:
        {
          result.setEmrPhase(name_, spin_t_, value.getRealValue());
          return result;
        }
        break;
      case ScanType::EmrLengthType:
        {
          result.setSize(name_, value.getSizeValue());
          return result;
        }
        break;
      default:
        throw NotImplementedError("Unknown ScanType.");
    }
    return result;
  }
  
  PulseSequence Selector::modify(
      PulseSequence &&seq, const ScanValueType &value) const
  {
    PulseSequence pseq = std::move(seq);
    switch(scan_t_){
      case ScanType::EmrGammaB1Type:
        {
          pseq.setEmrFreq(name_, spin_t_, value.getRealValue());
          return pseq;
        }
        break;
      case ScanType::EmrPhaseType:
        {
          pseq.setEmrPhase(name_, spin_t_, value.getRealValue());
          return pseq;
        }
        break;
      case ScanType::EmrLengthType:
        {
          pseq.setSize(name_, value.getSizeValue());
          return pseq;
        }
        break;
      default:
        throw NotImplementedError("Unknown ScanType.");
    }
    return pseq;
  }

  Parameters::Parameters(const Magnet &m, const Gyrotron &g, const Probe &p,
        const SpinSys &spin_sys, const PulseSequence &seq,
        const SpinType &acq_spin, const std::vector<Euler<>> &spin_sys_eulers)
    : magnet(m), gyrotron(g), probe(p), spin_sys(spin_sys),
    seq(seq), acq_spin(acq_spin), spin_sys_eulers(spin_sys_eulers)
  {}

  ScanResults1D scan1d(
      const Parameters &params,
      const Selector &selector,
      const Range &range,
      int ncores)
  {
    PulseSequence pseq_ref;
    const double intensity_ref = DnpRunner::calcPowderIntensity(
        params.magnet, params.gyrotron, params.probe,
        params.spin_sys, pseq_ref, params.acq_spin,
        params.spin_sys_eulers, ncores);
    ScanResults1D results;
    vector<double> values = range.values();
    for(const auto &value : values){
      PulseSequence pseq = selector.modify(params.seq, value);
      double intensity = DnpRunner::calcPowderIntensity(
          params.magnet, params.gyrotron, params.probe,
          params.spin_sys, pseq, params.acq_spin,
          params.spin_sys_eulers, ncores);
      results.push_back(make_pair(value, intensity/intensity_ref));
      std::cout << "." << std::flush;
    }
    
    return results;
  }

  ScanResults2D scan2d(
      const Parameters &params,
      const Selector &selector1,
      const Range &range1,
      const Selector &selector2,
      const Range &range2,
      int ncores)
  {
    PulseSequence pseq_ref;
    const double intensity_ref = DnpRunner::calcPowderIntensity(
        params.magnet, params.gyrotron, params.probe,
        params.spin_sys, pseq_ref, params.acq_spin,
        params.spin_sys_eulers, ncores);
    ScanResults2D results;
    vector<double> values1 = range1.values();
    vector<double> values2 = range2.values();
    for(const auto &v1 : values1){
      for(const auto &v2 : values2) {
        auto seq_temp = selector1.modify(params.seq, v1);
        auto pseq = selector2.modify(std::move(seq_temp), v2);
        double intensity = DnpRunner::calcPowderIntensity(
            params.magnet, params.gyrotron, params.probe,
            params.spin_sys, pseq, params.acq_spin,
            params.spin_sys_eulers, ncores);
        results.push_back(
            make_tuple(v1, v2, intensity/intensity_ref));
        std::cout << "." << std::flush;
      }
      std::cout << std::endl;
    }
    return results;
  }

  vector<double> populateValues(double val_beg, double val_end, uint64_t cnt)
  {
    const double diff = val_end - val_beg;
    uint64_t nsteps = (cnt == 0) + (cnt != 0) * (cnt - 1);
    const double step = diff/static_cast<double>(nsteps);
    vector<double> results;
    double elem = val_beg;
    for(uint64_t i = 0; i < cnt; ++i){
      results.push_back(elem);
      elem += step;
    }
    return results;
  }

  std::vector<std::uint64_t> populateValues(
      std::uint64_t val_beg, std::uint64_t val_end, std::uint64_t cnt)
  {
    vector<uint64_t> results;
    auto step = (val_end - val_beg)/((cnt == 0) + (cnt != 0) * (cnt - 1));
    for(auto val = val_beg; val < val_end; val += step){
      results.push_back(val);
    }
    return results;
  }
}   // namespace dnpsoup

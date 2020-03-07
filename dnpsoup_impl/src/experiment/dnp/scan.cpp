#include "dnpsoup_core/experiment/dnp/scan.h"
#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/errors.h"
#include <cmath>

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

  vector<double> Range::values() const
  {
    return populateValues(beg, end, cnt);
  }

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

  PulseSequence Selector::modify(const PulseSequence &seq, double value) const
  {
    PulseSequence result = seq;
    switch(scan_t_){
      case ScanType::EmrGammaB1Type:
        {
          result.setEmrFreq(name_, spin_t_, value);
          return result;
        }
        break;
      case ScanType::EmrPhaseType:
        {
          result.setEmrPhase(name_, spin_t_, value);
          return result;
        }
        break;
      case ScanType::EmrLengthType:
        {
          std::uint64_t sz = static_cast<std::uint64_t>(std::round(
                value/result.getIncrement()));

          result.setSize(name_, sz);
          return result;
        }
        break;
      default:
        throw NotImplementedError("Unknown ScanType.");
    }
    return result;
  }
  
  PulseSequence Selector::modify(PulseSequence &&seq, double value) const
  {
    PulseSequence pseq = std::move(seq);
    switch(scan_t_){
      case ScanType::EmrGammaB1Type:
        {
          pseq.setEmrFreq(name_, spin_t_, value);
          return pseq;
        }
        break;
      case ScanType::EmrPhaseType:
        {
          pseq.setEmrPhase(name_, spin_t_, value);
          return pseq;
        }
        break;
      case ScanType::EmrLengthType:
        {
          std::uint64_t sz = static_cast<std::uint64_t>(std::round(
                value/seq.getIncrement()));

          pseq.setSize(name_, sz);
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
    ScanResults1D results;
    vector<double> values = range.values();
    for(const auto &value : values){
      PulseSequence pseq = selector.modify(params.seq, value);
      double intensity = DnpRunner::calcPowderIntensity(
          params.magnet, params.gyrotron, params.probe,
          params.spin_sys, pseq, params.acq_spin,
          params.spin_sys_eulers, ncores);
      results.push_back(make_pair(value, intensity));
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
        results.push_back(make_tuple(v1, v2, intensity));
      }
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
}   // namespace dnpsoup

#include "dnpsoup_core/experiment/dnp/scan.h"
#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/errors.h"
#include <cmath>

using namespace std;


namespace dnpsoup {
  vector<double> Range::values() const
  {
    return populateValues(beg, end, cnt);
  }

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
      case ScanType::GammaB1Type:
        {
          result.setEmrFreq(name_, spin_t_, value);
          return result;
        }
        break;
      case ScanType::PhaseType:
        {
          result.setEmrPhase(name_, spin_t_, value);
          return result;
        }
        break;
      case ScanType::LengthType:
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
      case ScanType::GammaB1Type:
        {
          pseq.setEmrFreq(name_, spin_t_, value);
          return pseq;
        }
        break;
      case ScanType::PhaseType:
        {
          pseq.setEmrPhase(name_, spin_t_, value);
          return pseq;
        }
        break;
      case ScanType::LengthType:
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

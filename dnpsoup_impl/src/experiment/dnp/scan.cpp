#include "dnpsoup_core/experiment/dnp/scan.h"
#include "dnpsoup_core/experiment/DnpRunner.h"
#include <cmath>

using namespace std;


namespace dnpsoup {
  // enum class ScanType : int
  // {
  //   FieldType = 0,
  //   EmOffsetType = 1,

  //   GammaB1Type = 10,
  //   PhaseType = 11,
  //   LengthType = 12,

  //   EulerType = 100,

  //   T1Type = 200,
  //   T2Type = 201
  // };
  ScanResults1D scan1dEmrFreq(
      const Parameters &params,
      const std::string &name,
      const SpinType t,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores)
  {
    vector<double> gamma_b1s = populateValues(
        value_begin, value_end, cnt);

    ScanResults1D results;
    for(const auto &val : gamma_b1s){
      PulseSequence seq = params.seq;
      seq.setEmrFreq(name, t, val);
      double intensity = DnpRunner::calcPowderIntensity(
          params.magnet, params.gyrotron, params.probe,
          params.spin_sys, seq, params.acq_spin,
          params.spin_sys_eulers, ncores);
      results.push_back(make_pair(val, intensity));
    }

    return results;
  }

  ScanResults1D scan1dEmrPhase(
      const Parameters &params,
      const std::string &name,
      const SpinType t,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores)
  {
    vector<double> phases = populateValues(
        value_begin, value_end, cnt);

    ScanResults1D results;
    for(const auto &val : phases){
      PulseSequence seq = params.seq;
      seq.setEmrPhase(name, t, val);

      double intensity = DnpRunner::calcPowderIntensity(
          params.magnet, params.gyrotron, params.probe,
          params.spin_sys, seq, params.acq_spin,
          params.spin_sys_eulers, ncores);
      results.push_back(make_pair(val, intensity));
    }

    return results;
  }

  ScanResults1D scan1dEmrLength(
      const Parameters &params,
      const std::string &name,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores)
  {
    vector<double> lengths = populateValues(
        value_begin, value_end, cnt);

    ScanResults1D results;
    for(const auto &val : lengths){
      PulseSequence seq = params.seq;
      std::uint64_t sz = static_cast<std::uint64_t>(std::round(
            val/seq.getIncrement()));

      seq.setSize(name, sz);

      double intensity = DnpRunner::calcPowderIntensity(
          params.magnet, params.gyrotron, params.probe,
          params.spin_sys, seq, params.acq_spin,
          params.spin_sys_eulers, ncores);
      results.push_back(make_pair(val, intensity));
    }

    return results;
  }

  vector<double> scan1d(
      const Parameters &params,
      const std::vector<PulseSequence> &seqs,
      int ncores)
  {
    vector<double> results;
    for(const auto &seq : seqs){
      double intensity = DnpRunner::calcPowderIntensity(
          params.magnet, params.gyrotron, params.probe,
          params.spin_sys, seq, params.acq_spin,
          params.spin_sys_eulers, ncores);
      results.push_back(intensity);
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

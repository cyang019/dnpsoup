#include "dnpsoup_core/experiment/dnp/scan.h"
#include "dnpsoup_core/experiment/DnpRunner.h"

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
  
  ScanResults1D scan1d(
      const Parameters &params,
      const ScanType &scan_t,
      const std::string &scan_name,
      double value_begin,
      double value_end,
      std::uint64_t cnt,
      int ncores)
  {
    switch(scan_t){
      case ScanType::GammaB1Type:
        {
          vector<double> gamma_b1s = populateValues(
              value_begin, value_end, cnt);
          vector<PulseSequence> sequences;
          for(const auto &val : gamma_b1s){
            PulseSequence seq;
            seq.setParam
          }
          

        }
        break;
    }
  }

  ScanResults2D scan2d(
      const Parameters &params,
      const ScanType &scan_t_1,
      double value_begin_1,
      double value_end_1,
      double value_step_1,
      const ScanType &scan_t_2,
      double value_begin_2,
      double value_end_2,
      double value_step_2,
      int ncores)
  {
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

#include "dnpsoup_api.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/experiment/hardware.h"
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

extern "C" {
  void eigenValues(const char *json_str, 
      double *values)
  {
    dnpsoup::json j;
    istringstream iss(json_str);
    iss >> j;

    auto magnet = dnpsoup::Magnet(j);
    auto gyrotron = dnpsoup::Gyrotron(j);
    auto probe = dnpsoup::Probe(j);
    
    if(j.find("spinsys") == j.end()){
      throw dnpsoup::InputError("Cannot find 'spinsys' in input json.");
    }
    const string spins_str = j["spinsys"].dump(4);
    istringstream spins_iss(spins_str);
    dnpsoup::SpinSys spins;
    spins_iss >> spins;

    const string pulse_seq_str = j["pulse_sequence"].dump(4);
    istringstream pulse_seq_iss(pulse_seq_str);
    dnpsoup::PulseSequence seq;
    pulse_seq_iss >> seq;
    
    auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
    if(j.find("euler") != j.end()){
      if(j["euler"].find("alpha") != j["euler"].end()){
        euler.alpha(j["euler"]["alpha"].get<double>());
      }
      if(j["euler"].find("beta") != j["euler"].end()){
        euler.beta(j["euler"]["beta"].get<double>());
      }
      if(j["euler"].find("gamma") != j["euler"].end()){
        euler.gamma(j["euler"]["gamma"].get<double>());
      }
    }

    auto results = dnpsoup::DnpRunner::calcEigenValues(
        magnet, gyrotron, probe, spins, seq,
        euler);

    if(results.size() == 0) return;
    unsigned idx = 0;
    for(const auto &row : results){
      for(const auto &val : row){
        values[idx++] = val;
      }
    }
    
    return;
  }

  double calculateIntensity(const char *json_str)
  {
    dnpsoup::json j;
    istringstream iss(json_str);
    iss >> j;

    auto magnet = dnpsoup::Magnet(j);
    auto gyrotron = dnpsoup::Gyrotron(j);
    auto probe = dnpsoup::Probe(j);
    
    if(j.find("spinsys") == j.end()){
      throw dnpsoup::InputError("Cannot find 'spinsys' in input json.");
    }
    const string spins_str = j["spinsys"].dump(4);
    istringstream spins_iss(spins_str);
    dnpsoup::SpinSys spins;
    spins_iss >> spins;

    const string pulse_seq_str = j["pulse_sequence"].dump(4);
    istringstream pulse_seq_iss(pulse_seq_str);
    dnpsoup::PulseSequence seq;
    pulse_seq_iss >> seq;
    
    auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
    if(j.find("euler") != j.end()){
      if(j["euler"].find("alpha") != j["euler"].end()){
        euler.alpha(j["euler"]["alpha"].get<double>());
      }
      if(j["euler"].find("beta") != j["euler"].end()){
        euler.beta(j["euler"]["beta"].get<double>());
      }
      if(j["euler"].find("gamma") != j["euler"].end()){
        euler.gamma(j["euler"]["gamma"].get<double>());
      }
    }

    if(j.find("acq") == j.end()){
      throw dnpsoup::InputError("Cannot find 'acq' in input json.");
    }
    dnpsoup::SpinType t_acq = dnpsoup::toSpinType(j["acq"].get<string>());

    double result = dnpsoup::DnpRunner::calcIntensity(
        magnet, gyrotron, probe, spins, seq,
        t_acq, euler);

    return result;
  }

  double powderIntensity(const char *json_str,
      int ncores)
  {
    dnpsoup::json j;
    istringstream iss(json_str);
    iss >> j;

    auto magnet = dnpsoup::Magnet(j);
    auto gyrotron = dnpsoup::Gyrotron(j);
    auto probe = dnpsoup::Probe(j);
    
    if(j.find("spinsys") == j.end()){
      throw dnpsoup::InputError("Cannot find 'spinsys' in input json.");
    }
    const string spins_str = j["spinsys"].dump(4);
    istringstream spins_iss(spins_str);
    dnpsoup::SpinSys spins;
    spins_iss >> spins;

    const string pulse_seq_str = j["pulse_sequence"].dump(4);
    istringstream pulse_seq_iss(pulse_seq_str);
    dnpsoup::PulseSequence seq;
    pulse_seq_iss >> seq;
    
    auto eulers = vector<dnpsoup::Euler<>>();
    if(j.find("eulers") != j.end()){
      for(const auto &euler_js : j["eulers"]){
        auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
        if(euler_js["euler"].find("alpha") != euler_js["euler"].end()){
          euler.alpha(euler_js["euler"]["alpha"].get<double>());
        }
        if(euler_js["euler"].find("beta") != euler_js["euler"].end()){
          euler.beta(euler_js["euler"]["beta"].get<double>());
        }
        if(euler_js["euler"].find("gamma") != euler_js["euler"].end()){
          euler.gamma(euler_js["euler"]["gamma"].get<double>());
        }
      }
    }

    if(j.find("acq") == j.end()){
      throw dnpsoup::InputError("Cannot find 'acq' in input json.");
    }
    dnpsoup::SpinType t_acq = dnpsoup::toSpinType(j["acq"].get<string>());

    double result = dnpsoup::DnpRunner::calcPowderIntensity(
        magnet, gyrotron, probe, spins, seq,
        t_acq, eulers, ncores);

    return result;
  }

  void fieldProfile(const char *json_str, double *values, 
      int ncores)
  {
    dnpsoup::json j;
    istringstream iss(json_str);
    iss >> j;

    auto probe = dnpsoup::Probe(j);
    
    if(j.find("spinsys") == j.end()){
      throw dnpsoup::InputError("Cannot find 'spinsys' in input json.");
    }
    const string spins_str = j["spinsys"].dump(4);
    istringstream spins_iss(spins_str);
    dnpsoup::SpinSys spins;
    spins_iss >> spins;

    const string pulse_seq_str = j["pulse_sequence"].dump(4);
    istringstream pulse_seq_iss(pulse_seq_str);
    dnpsoup::PulseSequence seq;
    pulse_seq_iss >> seq;
    
    auto eulers = vector<dnpsoup::Euler<>>();
    if(j.find("eulers") != j.end()){
      for(const auto &euler_js : j["eulers"]){
        auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
        if(euler_js["euler"].find("alpha") != euler_js["euler"].end()){
          euler.alpha(euler_js["euler"]["alpha"].get<double>());
        }
        if(euler_js["euler"].find("beta") != euler_js["euler"].end()){
          euler.beta(euler_js["euler"]["beta"].get<double>());
        }
        if(euler_js["euler"].find("gamma") != euler_js["euler"].end()){
          euler.gamma(euler_js["euler"]["gamma"].get<double>());
        }
        eulers.push_back(euler);
      }
    }

    if(j.find("acq") == j.end()){
      throw dnpsoup::InputError("Cannot find 'acq' in input json.");
    }
    dnpsoup::SpinType t_acq = dnpsoup::toSpinType(j["acq"].get<string>());

    if(j.find("fields") != j.end()){
      vector<dnpsoup::Magnet> fields;
      for(const auto &val : j["fields"]){
        double f_val = val.get<double>();
        fields.emplace_back(f_val);
      }
      auto gyrotron = dnpsoup::Gyrotron(j);

      auto results = dnpsoup::DnpRunner::calcFieldProfile(
          fields, gyrotron, probe, spins, seq, 
          t_acq, eulers, ncores);

      unsigned idx = 0;
      for(const auto &val : results){
        values[idx++] = val.second;
      }
    }
    else if(j.find("em_frequencies") != j.end()){
      auto magnet = dnpsoup::Magnet(j);
      vector<dnpsoup::Gyrotron> ems;
      for(const auto &val : j["em_frequencies"]){
        double em_val = val.get<double>();
        ems.emplace_back(em_val);
      }

      auto results = dnpsoup::DnpRunner::calcFieldProfile(
          magnet, ems, probe, spins, seq,
          t_acq, eulers, ncores);
      unsigned idx = 0;
      for(const auto &val : results){
        values[idx++] = val.second;
      }
    }

    throw dnpsoup::InputError("Cannot find 'fields' or 'em_frequencies' in input json.");
  }
} // extern C



#include "dnpsoup_core/experiment/hardware.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"
#include <iostream>

using namespace std;

namespace dnpsoup {
  Magnet::Magnet(const json &j)
  {
    if(j.find("Magnet") != j.end()){
      if(j.at("Magnet").find("b0") != j.at("Magnet").end()){
        b0 = j.at("Magnet").at("b0").get<double>();
      }
      else{
        throw InputError("Cannot find 'b0' in Magnet.");
      }
    }
    else if(j.find("b0") != j.end()){
      b0 = j.at("b0").get<double>();
    }
    else{
      throw InputError("Neither 'Magnet' nor 'b0' is in the input json.");
    }
  }

  Gyrotron::Gyrotron(const json &j)
  {
    if(j.find("Gyrotron") != j.end()){
      if(j["Gyrotron"].find("em_frequency") != j["Gyrotron"].end()){
        em_frequency = j["Gyrotron"]["em_frequency"].get<double>();
      }
      else{
        throw InputError("Cannot find 'em_frequency' in Gyrotron.");
      }
    }
    else if(j.find("em_frequency") != j.end()){
      em_frequency = j["em_frequency"].get<double>();
    }
    else{
      throw InputError("Neither 'Gyrotron' nor 'em_frequency' is in the input json.");
    }
  }

  Probe::Probe(const json &j)
    : mas_frequency(0.0), 
    magic_angle(Euler<>(0.0, ::dnpsoup::magic_angle, 0.0)),
    mas_increment(-1.0), temperature(100.0)
  {
    if(j.find("Probe") != j.end()){
      auto probe_js = j["Probe"];
      if(probe_js.find("mas_frequency") != probe_js.end()){
        mas_frequency = probe_js["mas_frequency"].get<double>();
      }
      magic_angle = Euler<>(0.0, ::dnpsoup::magic_angle, 0.0);
      if(probe_js.find("magic_angle") != probe_js.end()){
        double beta = probe_js["magic_angle"].get<double>();
        magic_angle.beta(beta);
      }

      if(probe_js.find("temperature") != probe_js.end()){
        temperature = probe_js["temperature"].get<double>();
      }

      if(probe_js.find("mas_increment") != probe_js.end()){
        mas_increment = probe_js["mas_increment"].get<double>();
      }
    }
    else{
      throw InputError("Cannot find 'Probe' in the input json.");
    }
  }
} // namespace dnpsoup

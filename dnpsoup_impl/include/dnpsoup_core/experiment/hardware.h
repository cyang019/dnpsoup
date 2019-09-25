#ifndef DNPSOUP_HARDWARE_H
#define DNPSOUP_HARDWARE_H

#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/rotation/Euler.h"


namespace dnpsoup {
  struct Magnet {
    Magnet() : b0(0.0) {}
    Magnet(double val) : b0(val) {}
    Magnet(const json &j);
    Magnet(const Magnet &) = default;
    Magnet(Magnet &&) noexcept = default;
    Magnet& operator=(const Magnet &) = default;
    Magnet& operator=(Magnet &&) noexcept = default;
    ~Magnet() {}

    double b0;
  };  // class Magnet

  /// e in rotating frame. Nuclei's all in lab frame
  struct Gyrotron {
    Gyrotron() : em_frequency(0.0) {}
    Gyrotron(double val) : em_frequency(val) {}
    Gyrotron(const json &j);
    Gyrotron(const Gyrotron &) = default;
    Gyrotron(Gyrotron &&) noexcept = default;
    Gyrotron& operator=(const Gyrotron &) = default;
    Gyrotron& operator=(Gyrotron &&) noexcept = default;
    ~Gyrotron() {}

    double em_frequency;   ///< in Hz
  };

  struct Probe {
    Probe() 
      : mas_frequency(0.0), magic_angle(0.0, ::dnpsoup::magic_angle, 0.0),
        mas_increment(-1.0), temperature(77.0) {}
    Probe(double mas_freq, double t) 
      : mas_frequency(mas_freq), magic_angle(0.0, ::dnpsoup::magic_angle, 0.0),
        mas_increment(-1.0), temperature(t)
    {}
    Probe(const json &j);
    Probe(const Probe &) = default;
    Probe(Probe &&) noexcept = default;
    Probe& operator=(const Probe &) = default;
    Probe& operator=(Probe &&) noexcept = default;
    ~Probe() {}

    double mas_frequency;
    Euler<> magic_angle;
    double mas_increment;
    double temperature;
  };
} // namespace dnpsoup

#endif

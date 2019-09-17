#ifndef DNPSOUP_HARDWARE_H
#define DNPSOUP_HARDWARE_H



namespace dnpsoup {
  struct Magnet {
    Magnet() : b0(0.0) {}
    Magnet(double val) : b0(val) {}
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
    Gyrotron(const Gyrotron &) = default;
    Gyrotron(Gyrotron &&) noexcept = default;
    Gyrotron& operator=(const Gyrotron &) = default;
    Gyrotron& operator=(Gyrotron &&) noexcept = default;
    ~Gyrotron() {}

    double em_frequency;   ///< in Hz
  };

  struct Probe {
    Probe() : mas_frequency(0.0), temperature(77.0) {}
    Probe(double mas, double t) : mas_frequency(mas), temperature(t) {}
    Probe(const Probe &) = default;
    Probe(Probe &&) noexcept = default;
    Probe& operator=(const Probe &) = default;
    Probe& operator=(Probe &&) noexcept = default;
    ~Probe() {}

    double mas_frequency;
    double temperature;
  };
} // namespace dnpsoup

#endif

#ifndef DNPSOUP_HARDWARE_H
#define DNPSOUP_HARDWARE_H



namespace dnpsoup {
  struct Magnet {
    double b0;
  };  // class Magnet

  /// e in rotating frame. Nuclei's all in lab frame
  struct Gyrotron {
    double em_frequency;   ///< in Hz
  };

  class Probe {
    double mas_frequency;
  };
} // namespace dnpsoup

#endif

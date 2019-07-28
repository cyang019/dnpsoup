#ifndef DNPSOUP_EXPERIMENTTYPES_H
#define DNPSOUP_EXPERIMENTTYPES_H

namespace dnpsoup {
  class Experiment {
  public:
    Experiment() : m_b0(0.0) {}
    Experiment(const Experiment &) = default;
    Experiment(Experiment &&) noexcept = default;
    Experiment& operator=(const Experiment &) = default;
    Experiment& operator=(Experiment &&) noexcept = default;
    virtual ~Experiment() {}

    double getB0() const { return m_bo; }
    Experiment& setB0(double val) { m_b0 = val; return *this; }
  private:
    double m_b0;
  };  // class Experiment

  /// e in rotating frame. Nuclei's all in lab frame
  class DnpExperiment : public Experiment
  {
  private:
    double m_rf_offset;   ///< in Hz
  };

  /// everything in the rotating frame
  class NmrExperiment : public Experiment
  {
  };
} // namespace dnpsoup

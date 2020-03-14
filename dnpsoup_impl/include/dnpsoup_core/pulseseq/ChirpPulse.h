#ifndef DNPSOUP_CHIRPPULSE_H
#define DNPSOUP_CHIRPPULSE_H

#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"

namespace dnpsoup {
  namespace pulseseq {
    /// @brief chirped pulse
    class ChirpPulse : public SubSequenceInterface {
    public:
      ChirpPulse();
      ChirpPulse(const Name &);
      ChirpPulse(std::uint64_t, const Name &);
      ChirpPulse(std::uint64_t, const Name &, 
          SpinType t, double span, double inc);
      ChirpPulse(const ChirpPulse &) = default;
      ChirpPulse(ChirpPulse &&) noexcept = default;
      ChirpPulse& operator=(const ChirpPulse &) = default;
      ChirpPulse& operator=(ChirpPulse &&) noexcept = default;
      virtual ~ChirpPulse() {};

      virtual std::unique_ptr<SubSequenceInterface> copy() const override;

      virtual std::tuple<Component, std::uint64_t, std::uint64_t> next(
          std::map<Name, Component> *components,
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *m_sections
          ) override;
      virtual SequenceType type() const override;
      virtual std::vector<Name> getNames() const override;
      virtual void resetIndex(
          std::map<Name, std::unique_ptr<SubSequenceInterface>> *m_sections) override;
    private:
      Name component_name_;
      SpinType spin_t_;
      double inc_;
      double f0_; ///< chirp frequency begin
      double c_;  ///< chirp rate

      /// $phi(t) = phi_0 + 0.5 * c * t^2 + f_0 * t$
      /// f_0: chirp beginning frequency: -span/2
      /// c: chirp rate: span/duration
      /// @return f0, c
      std::pair<double, double> gen_chirp_(double span, double inc, std::uint64_t sz);
    // inherited
    //protected:
    //  std::uint64_t sz_;   ///< # of iterations
    //  std::uint64_t idx_;  ///< current location in chirp
    };
  } // namespace pulseseq
} // namespace dnpsoup

#endif


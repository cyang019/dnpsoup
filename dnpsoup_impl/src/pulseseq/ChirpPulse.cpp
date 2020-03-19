#include "dnpsoup_core/pulseseq/ChirpPulse.h"
#include "dnpsoup_core/errors.h"
#include <cmath>


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    ChirpPulse::ChirpPulse()
      : SubSequenceInterface(), spin_t_(SpinType::Null),
      inc_(1.0e-9)
    {
      std::tie(f0_, c_) = gen_chirp_(0.0, inc_, sz_);
    }

    ChirpPulse::ChirpPulse(const Name &component_name)
      : SubSequenceInterface(1), component_name_(component_name),
      spin_t_(SpinType::Null), inc_(1.0e-9)
    {
      std::tie(f0_, c_) = gen_chirp_(0.0, inc_, sz_);
    }

    ChirpPulse::ChirpPulse(std::uint64_t sz, const Name &component_name)
      : SubSequenceInterface(sz), component_name_(component_name),
      spin_t_(SpinType::Null), inc_(1.0e-9)
    {
      std::tie(f0_, c_) = gen_chirp_(0.0, inc_, sz_);
    }

    ChirpPulse::ChirpPulse(std::uint64_t sz, const Name &component_name,
        SpinType t, double span, double inc)
      : SubSequenceInterface(sz), component_name_(component_name),
      spin_t_(t), inc_(inc)
    {
      std::tie(f0_, c_) = gen_chirp_(span, inc_, sz_);
    }

    std::unique_ptr<SubSequenceInterface> ChirpPulse::copy() const
    {
      auto res = make_unique<ChirpPulse>(*this);
      res->idx_ = 0u;
      return res;
    }

    std::tuple<Component, std::uint64_t, std::uint64_t> ChirpPulse::next(
        std::map<Name, Component> *components,
        [[maybe_unused]] std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(idx_ >= sz_){
        idx_ = 0;
        return make_tuple(Component(), 0, sz_);
      }

      if(components->find(this->component_name_) == components->end()){
        const string err_str = 
          "Inside " + this->name + ": " + component_name_ + " not found.";
        throw PulseSequenceError(err_str);
      }
      if(approxEqual(c_, 0.0, eps, eps)){
        idx_ = sz_;
        return make_tuple(components->at(component_name_), sz_, 0);
      }
      const double t = static_cast<double>(idx_) * inc_;
      /// in degree
      const double phase_t = 360.0 * (0.5 * c_ * t + f0_) * t;
      Component comp = components->at(component_name_);
      comp[spin_t_].phase = fmod(comp[spin_t_].phase + phase_t, 360.0);
      auto this_idx_ = idx_;

      ++idx_;
      return make_tuple(comp, 1, this_idx_);
    }

    SequenceType ChirpPulse::type() const 
    { return SequenceType::ChirpType; }

    std::vector<Name> ChirpPulse::getNames() const
    {
      std::vector<Name> res = {component_name_};
      return res;
    }

    void ChirpPulse::resetIndex(
        [[maybe_unused]]
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections)
    {
      idx_ = 0;
    }

    std::pair<double, double> ChirpPulse::gen_chirp_(
        double span, double inc, std::uint64_t sz)
    {
      const double f0 = span/(-2.0);
      const double total_duration = inc * static_cast<double>(sz);
      const double rate = span / total_duration;
      return make_pair(f0, rate);
    }
  } // namespace pulseseq
} // namespace dnpsoup


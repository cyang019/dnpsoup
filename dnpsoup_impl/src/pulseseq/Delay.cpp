#include "dnpsoup_core/pulseseq/Delay.h"


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    Delay::Delay()
      : SubSequenceInterface()
    {}

    Delay::Delay(std::uint64_t sz)
      : SubSequenceInterface(sz)
    {}

    std::unique_ptr<SubSequenceInterface> Delay::copy() const
    {
      auto res = make_unique<Delay>(*this);
      res->idx_ = 0u;
      return res;
    }

    std::tuple<Component, std::uint64_t, std::uint64_t> Delay::next(
        [[maybe_unused]] std::map<Name, Component> *components,
        [[maybe_unused]] std::map<Name, std::unique_ptr<SubSequenceInterface>> *sections
        )
    {
      if(idx_ == 0){
        idx_ += sz_;
        return make_tuple(Component(), sz_, 0);
      } else {
        idx_ = 0;
        return make_tuple(Component(), 0, sz_);
      }
    }

    SequenceType Delay::type() const 
    { return SequenceType::DelayType; }

    std::vector<Name> Delay::getNames() const
    {
      auto res = std::vector<Name>();
      return res;
    }

    void Delay::resetIndex(
        [[maybe_unused]]
        std::map<Name, std::unique_ptr<SubSequenceInterface>> *m_sections)
    {
      idx_ = 0;
    }
  } // namespace pulseseq
} // namespace dnpsoup


#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>

namespace {
    using SpinSys = dnpsoup::SpinSys;
    using SpinType = dnpsoup::SpinType;

    TEST(TestDnpsoup, CWCrossEffectXtal){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 2.02, 2.06, 2.09);
      spins.addSpin(2, SpinType::H, 2.5, 0.0, 0.0);
      spins.addSpin(3, SpinType::H, 0.0, 2.0, 1.0);
      spins.irradiateOn(SpinType::e);
      spins.acquireOn(SpinType::H);

      dnpsoup::PulseSequence p;
      dnpsoup::pulseseq::EMRadiation emr(1.0e6, 0.0, 400.0e6);
      dnpsoup::pulseseq::Component c;
      c.insert({SpinType::e, emr});
      p.set("emr", c);
      auto uptr_sec = std::make_unique<dnpsoup::pulseseq::Pulse>("emr");
      p.set("cw", std::move(uptr_sec));
      std::vector<std::string> seq_names = { "cw" };
      p.set(seq_names);
      std::cout << "CW Pulse Sequence:\n" << p << std::endl;
      std::ostringstream buffer;
      buffer << p;

      auto magnet = dnpsoup::Magnet(9.4);
      auto gyrotron = dnpsoup::Gyrotron(300.0e9);
      auto probe = dnpsoup::Probe(0.0, 77.0);

      dnpsoup::DnpRunner runner;
      auto res = runner.calcIntensity(
          magnet, gyrotron, probe, spins, 
          buffer.str(), 
          SpinType::H, dnpsoup::Euler<>(0.0,0.0,0.0));
      std::cout << "Intensity with DNP: " << res << std::endl;

      c[SpinType::e].freq = 0.0;
      p.set("emr", c);
      std::ostringstream buffer2;
      buffer2 << p;
      auto res2 = runner.calcIntensity(
          magnet, gyrotron, probe, spins, buffer2.str(),
          SpinType::H, dnpsoup::Euler<>(0.0,0.0,0.0));
      std::cout << "Intensity without radiation: " << res2 << std::endl;
    }
} // namespace dnpsoup



#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <memory>

namespace {
    using SpinSys = dnpsoup::SpinSys;
    using SpinType = dnpsoup::SpinType;

    TEST(TestDnpsoup, TotapolMASEigen){
      constexpr double pi = ::dnpsoup::pi;
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 1.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::H, 0.0, 0.0, 0.0);
      spins.addSpin(3, SpinType::e, 9.3, 0.0, 0.0);
      spins.irradiateOn(SpinType::e);
      spins.acquireOn(SpinType::H);
      spins.setShielding(dnpsoup::SpinId(1), 2.0017, 2.006, 2.0094, 
          dnpsoup::Euler<>(0.0,0.0,0.0));
      spins.setShielding(dnpsoup::SpinId(3), 2.0017, 2.006, 2.0094, 
          dnpsoup::Euler<>(124.0/180.0 * pi,108.0/180.0 * pi,-107.0/180.0 * pi));
      spins.setT1(1, 1.0e-3);
      spins.setT1(2, 2.0);
      spins.setT1(3, 1.0e-3);
      spins.setT2(1, 2.0e-6);
      spins.setT2(2, 1.0e-3);
      spins.setT2(3, 2.0e-6);

      dnpsoup::PulseSequence p(1.0e-8);
      // [number] x increment duration
      auto uptr_sec = std::make_unique<dnpsoup::pulseseq::Delay>(800);
      p.set("delay", std::move(uptr_sec));
      std::vector<std::string> seq_names = { "delay" };
      p.set(seq_names);
      std::ostringstream buffer;
      buffer << p;

      auto magnet = dnpsoup::Magnet(9.402);
      // Hz
      auto gyrotron = dnpsoup::Gyrotron(263.00e9);
      // MAS, Temperature
      auto probe = dnpsoup::Probe(8e3, 77.0);

      auto res = dnpsoup::DnpRunner::calcEigenValues(
          magnet, gyrotron, probe, spins,
          buffer.str(), ::dnpsoup::Euler<>(320.0/180.0*pi, 141.0/180.0*pi, 80.0/180.0*pi));

      std::ofstream ofss("eigen_values.log");
      ofss << "EigenValues:\n";
      for(auto vals : res){
        for(auto val : vals){
          ofss << val << ",\t";
        }
        ofss << std::endl;
      }
      ofss << std::endl;
    }

    TEST(TestDnpsoup, SEXtal){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::H, 0.7, 0.7, 1.0);
      spins.irradiateOn(SpinType::e);
      spins.acquireOn(SpinType::H);
      spins.setShielding(dnpsoup::SpinId(1), 2.00263, 2.00259, 2.00234, dnpsoup::Euler<>(0.0,0.0,0.0));
      spins.setT1(dnpsoup::SpinId(1), 1.0e-3);
      spins.setT2(dnpsoup::SpinId(1), 2.0e-6);
      spins.setT1(dnpsoup::SpinId(2), 1.0);
      spins.setT2(dnpsoup::SpinId(2), 4.0e-3);

      auto magnet = dnpsoup::Magnet(9.36932070693522);
      auto gyrotron = dnpsoup::Gyrotron(263.0e9);
      auto probe = dnpsoup::Probe(0.0, 100.0);

      dnpsoup::PulseSequence p;
      p.setIncrement(1.0/gyrotron.em_frequency * 250.0);
      dnpsoup::pulseseq::EMRadiation emr(1.0e6, 0.0, 0.0e6);
      dnpsoup::pulseseq::Component c;
      c.insert({SpinType::e, emr});
      p.set("emr", c);
      // CW irradiation only has one component
      auto uptr_sec = std::make_unique<dnpsoup::pulseseq::Pulse>(200, "emr");
      p.set("cw", std::move(uptr_sec));
      std::vector<std::string> seq_names = { "cw" };
      p.set(seq_names);
      std::ostringstream buffer;
      buffer << p;

      auto res = dnpsoup::DnpRunner::calcIntensity(
          magnet, gyrotron, probe, spins, 
          buffer.str(), 
          SpinType::H, dnpsoup::Euler<>(0.0,0.0,0.0));
      std::cout << "Intensity with DNP: " << res << std::endl;

      c[SpinType::e].freq = 0.0;
      p.set("emr", c);
      std::ostringstream buffer2;
      buffer2 << p;
      auto res2 = dnpsoup::DnpRunner::calcIntensity(
          magnet, gyrotron, probe, spins, buffer2.str(),
          SpinType::H, dnpsoup::Euler<>(0.0,0.0,0.0));
      std::cout << "Intensity without radiation: " << res2 << std::endl;
    }

    TEST(TestDnpsoup, SEPowder){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0, true);
      spins.addSpin(2, SpinType::H, 0.7, 0.7, 1.0, true);
      spins.irradiateOn(SpinType::e);
      spins.acquireOn(SpinType::H);
      spins.setShielding(dnpsoup::SpinId(1), 2.00263, 2.00259, 2.00234, dnpsoup::Euler<>(0.0,0.0,0.0));
      spins.setT1(dnpsoup::SpinId(1), 1.0e-3);
      spins.setT2(dnpsoup::SpinId(1), 2.0e-6);
      spins.setT1(dnpsoup::SpinId(2), 1.0);
      spins.setT2(dnpsoup::SpinId(2), 4.0e-3);

      auto magnet = dnpsoup::Magnet(9.36932070693522);
      auto gyrotron = dnpsoup::Gyrotron(263.0e9);
      auto probe = dnpsoup::Probe(0.0, 77.0);

      dnpsoup::PulseSequence p;
      p.setIncrement(1.0/gyrotron.em_frequency * 250.0);
      dnpsoup::pulseseq::EMRadiation emr(1.0e6, 0.0, 0.0);
      dnpsoup::pulseseq::Component c;
      c.insert({SpinType::e, emr});
      p.set("emr", c);
      auto uptr_sec = std::make_unique<dnpsoup::pulseseq::Pulse>(200, "emr");
      p.set("cw", std::move(uptr_sec));
      std::vector<std::string> seq_names = { "cw" };
      p.set(seq_names);
      std::ostringstream buffer;
      buffer << p;

      auto eulers = dnpsoup::getZCWAngles(2);
      auto res = dnpsoup::DnpRunner::calcPowderIntensity(
          magnet, gyrotron, probe, spins, 
          buffer.str(), 
          SpinType::H, eulers);
      std::cout << "Intensity with DNP sequential: " << res << std::endl;
      auto res_concurrent = dnpsoup::DnpRunner::calcPowderIntensity(
          magnet, gyrotron, probe, spins, 
          buffer.str(), 
          SpinType::H, eulers, 4);
      ASSERT_DOUBLE_EQ(res, res_concurrent);

      std::cout << "Intensity with DNP using multiple cores: " << res_concurrent << std::endl;

      c[SpinType::e].freq = 0.0;
      p.set("emr", c);
      std::ostringstream buffer2;
      buffer2 << p;
      auto res2 = dnpsoup::DnpRunner::calcPowderIntensity(
          magnet, gyrotron, probe, spins, buffer2.str(),
          SpinType::H, eulers);
      std::cout << "Intensity without radiation: " << res2 << std::endl;
    }
    
#ifdef ENABLEFP
    TEST(TestDnpsoup, BDPASolidEffectFieldProfile){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::H, 2.0, 0.0, 0.0);
      spins.irradiateOn(SpinType::e);
      spins.acquireOn(SpinType::H);
      spins.setShielding(dnpsoup::SpinId(1), 2.00263, 2.00259, 2.00234, 
      spins.setT1(1, 1.0e-3);
      spins.setT1(2, 2.0);
      spins.setT2(1, 2.0e-6);
      spins.setT2(2, 1.0e-3);
      spins.setScalar(1, 2, 2.0e6);

      dnpsoup::PulseSequence p(5.0e-9);
      dnpsoup::pulseseq::EMRadiation emr(0.85e6, 0.0, 0.0);
      dnpsoup::pulseseq::Component c;
      c.insert({SpinType::e, emr});
      p.set("emr", c);
      // [number] x increment duration
      auto uptr_sec = std::make_unique<dnpsoup::pulseseq::Pulse>(10, "emr");
      p.set("cw", std::move(uptr_sec));
      std::vector<std::string> seq_names = { "cw" };
      p.set(seq_names);
      std::ostringstream buffer;
      buffer << p;

      std::vector<dnpsoup::Magnet> fields;
      for(int i = 0; i < 30; ++i){
        fields.push_back(dnpsoup::Magnet(9.397 + static_cast<double>(i) * 0.0002));
      }
      // Hz
      auto gyrotron = dnpsoup::Gyrotron(263.46e9);
      // MAS, Temperature
      auto probe = dnpsoup::Probe(8000.0, 77.0);

      auto eulers = dnpsoup::getZCWAngles(2); 
      auto res = dnpsoup::DnpRunner::calcFieldProfile(
          fields, gyrotron, probe, spins, 
          buffer.str(), 
          SpinType::H, eulers);
      std::cout << "Field Profile with DNP: ";
      for(auto f : res){
        std::cout << f << ", ";
      }
      std::cout << std::endl;

      c[SpinType::e].freq = 0.0;
      p.set("emr", c);
      std::ostringstream buffer2;
      buffer2 << p;
      auto res2 = dnpsoup::DnpRunner::calcFieldProfile(
          fields, gyrotron, probe, spins, buffer2.str(),
          SpinType::H, eulers);
      std::cout << "Field Profile without radiation: ";
      for(auto f : res2){
        std::cout << f << ", ";
      }
      std::cout << std::endl;
    }
#endif

    TEST(TestDnpsoup, CWCrossEffectXtal){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::e, 6.0, 0.0, 0.0);
      spins.addSpin(3, SpinType::H, 0.0, -2.0, 1.0);
      spins.setShielding(dnpsoup::SpinId(1), 2.09, 2.06, 2.02, dnpsoup::Euler<>(0.0,0.0,0.0));
      spins.setShielding(dnpsoup::SpinId(2), 2.09, 2.06, 2.02, dnpsoup::Euler<>(20.0,100.0,0.0));
      spins.setT1(dnpsoup::SpinId(1), 0.001);
      spins.setT1(dnpsoup::SpinId(2), 0.001);
      spins.setT1(dnpsoup::SpinId(3), 1.0);
      spins.setT2(dnpsoup::SpinId(1), 1.0e-6);
      spins.setT2(dnpsoup::SpinId(2), 1.0e-6);
      spins.setT2(dnpsoup::SpinId(3), 1.0e-3);
      spins.irradiateOn(SpinType::e);
      std::vector<dnpsoup::SpinId> acq_spins;
      acq_spins.push_back(dnpsoup::SpinId(3));
      spins.acquireOn(acq_spins);

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

      auto res = dnpsoup::DnpRunner::calcIntensity(
          magnet, gyrotron, probe, spins, 
          buffer.str(), 
          SpinType::H, dnpsoup::Euler<>(0.0,0.0,0.0));
      std::cout << "Intensity with DNP: " << res << std::endl;

      c[SpinType::e].freq = 0.0;
      p.set("emr", c);
      std::ostringstream buffer2;
      buffer2 << p;
      auto res2 = dnpsoup::DnpRunner::calcIntensity(
          magnet, gyrotron, probe, spins, buffer2.str(),
          SpinType::H, dnpsoup::Euler<>(0.0,0.0,0.0));
      std::cout << "Intensity without radiation: " << res2 << std::endl;
    }
} // namespace dnpsoup



#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using SpinSys = dnpsoup::SpinSys;
    using SpinType = dnpsoup::SpinType;

    TEST(TestDnpsoup, GenPackets){
      SpinSys spins;
      auto packets1 = spins.summarize<dnpsoup::DnpExperiment>();
      auto packets2 = spins.summarize<dnpsoup::NmrExperiment>();
      ASSERT_EQ(0u, packets1.getNumOfPackets());
      ASSERT_EQ(0u, packets2.getNumOfPackets());
    }

    TEST(TestDnpsoup, PacketsToMatrices){
      SpinSys spins;
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::H, 2.5, 0.0, 0.0);
      spins.addSpin(3, SpinType::H, 0.0, 2.0, 1.0);
      spins.irradiateOn(SpinType::e);

      auto packets1 = spins.summarize<dnpsoup::DnpExperiment>();
      auto packets2 = spins.summarize<dnpsoup::NmrExperiment>();
      ASSERT_EQ(4u, packets1.getNumOfPackets());
      ASSERT_EQ(4u, packets2.getNumOfPackets());

      auto mat = packets1.genMatrix();
      std::cout << "Matrix of a \'e H H\' spin system.\n";
      std::cout << mat << std::endl;
    }
} // namespace dnpsoup


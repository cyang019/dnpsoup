#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using SpinSys = dnpsoup::SpinSys;

    TEST(TestDnpsoup, SpinSysCtor){
      SpinSys spins;
      auto packets1 = spins.summarize<dnpsoup::DnpExperiment>();
      auto packets2 = spins.summarize<dnpsoup::NmrExperiment>();
      ASSERT_EQ(0u, packets1.getNumOfPackets());
      ASSERT_EQ(0u, packets2.getNumOfPackets());
    }
} // namespace dnpsoup


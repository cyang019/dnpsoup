#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>


namespace {
    using PulseSeq = dnpsoup::PulseSequence;

    TEST(TestDnpsoup, EmptyPulseSeq){
      auto p = PulseSeq();
      std::ostringstream oss;
      oss << p;
      std::istringstream iss(oss.str());
      PulseSeq p2;
      iss >> p2;
      ASSERT_EQ(0u, p2.size());
    }

    TEST(TestDnpsoup, PulseSeqUtopia){
      auto p = PulseSeq();

      std::cout << p << std::endl;
    }
} // namespace dnpsoup



#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>

namespace {
    TEST(TestDnpsoup, ZCW){
      auto angles = ::dnpsoup::getZCWAngles(4);
      for(const auto &angle : angles){
        std::cout << angle.alpha() << "\t" << angle.beta() << "\t" << angle.gamma() << "\n";
      }
      std::cout << std::endl;
      ASSERT_EQ(144u, angles.size());
    }
} // namespace



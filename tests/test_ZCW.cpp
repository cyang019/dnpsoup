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
      //for(const auto &angle : angles){
      //  std::cout << angle.alpha() << "\t" << angle.beta() << "\t" << angle.gamma() << "\n";
      //}
      std::cout << std::endl;
      ASSERT_EQ(144u, angles.size());

      auto angles2 = ::dnpsoup::getZCWAngles(15);
      std::cout << "Number of points in ZCW(15): " << angles2.size() << std::endl;
      double sum = 0.0;
      for(auto angle : angles2){
        sum += std::sin(angle.beta());
      }
      const double scaling_factor = static_cast<double>(angles2.size());
      sum = sum/scaling_factor/dnpsoup::pi * 4.0;
      ASSERT_NEAR(1.0, sum, 1.0e-5);
    }
} // namespace



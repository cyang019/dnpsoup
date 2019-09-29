#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using Euler = dnpsoup::Euler<>;

    TEST(TestDnpsoup, EulerMultiplication){
      Euler e1(0.0, 0.0, 0.0);
      Euler e2(0.0, dnpsoup::pi * 0.5, 0.0);
      Euler e3 = e2 * e1;
      std::cout << "e1: " << e1 << "\n"
                << "e2: " << e2 << "\n";
      std::cout << "e2 * e1:" << e3 << std::endl;

      ASSERT_DOUBLE_EQ(0.0, e3.alpha());
      ASSERT_DOUBLE_EQ(dnpsoup::pi * 0.5, e3.beta());
      ASSERT_DOUBLE_EQ(0.0, e3.gamma());
    }
} // namespace dnpsoup


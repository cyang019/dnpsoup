#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <vector>
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

    TEST(TestDnpsoup, EulerMultiplication2){
      constexpr size_t cnt = 100;
      constexpr double step_size1 = dnpsoup::pi * 2 / static_cast<double>(cnt);
      constexpr double step_size2 = (dnpsoup::pi - 0.11) / static_cast<double>(cnt);
      constexpr double eps = std::numeric_limits<double>::epsilon() * 60;
      for(size_t i = 0; i < cnt; ++i) {
        const double val1 = step_size1 * static_cast<double>(i);
        const double val2 = step_size2 * static_cast<double>(i);
        auto temp1 = Euler(val1, 0.1, 0.0);
        auto temp2 = Euler(0.0, val2, 0.0);
        auto temp3 = temp1 * temp2;
        //auto temp4 = temp2 * temp1;
        //std::cout << "e1:\n" << temp1 << "\n\n"
        //          << "e2:\n" << temp2 << "\n\n"
        //          << "e1 * e2:\n" << temp3 << "\n\n";

        ASSERT_NEAR(val1, temp3.alpha(), eps);
        ASSERT_NEAR(val2 + 0.1, temp3.beta(), eps);
        ASSERT_NEAR(0.0, temp3.gamma(), eps);
        //ASSERT_NEAR(val1, temp4.alpha(), eps);
        //ASSERT_NEAR(val2, temp4.beta(), eps);
        //ASSERT_NEAR(0.0, temp4.gamma(), eps);
      }
    }
} // namespace dnpsoup


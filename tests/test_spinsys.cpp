#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using SpinSys = dnpsoup::SpinSys;
    using SpinType = dnpsoup::SpinType;

    TEST(TestDnpsoup, SpinSysCtor){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::H, 2.5, 0.0, 0.0);
      spins.addSpin(3, SpinType::H, 0.0, 2.0, 1.0);

      auto dim = spins.calcTotalDimension();
      auto dims = spins.calcDimensions();
      ASSERT_EQ(8u, dim);

      auto dims_desired = std::vector<std::size_t>({2u, 2u, 2u});
      ASSERT_TRUE(dims_desired == dims);
    }
} // namespace dnpsoup

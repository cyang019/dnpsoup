#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using Coord = dnpsoup::Coordinate;

    TEST(TestDnpsoup, Coordinate){
      Coord c1(1.0, 0.0, 0.0);
      Coord c2(2.0, 0.0, 0.0);
      ASSERT_TRUE(c2 - c1 == Coord(1.0, 0.0, 0.0));
      auto [phi, theta] = calcAnglesWithZ(c1);
      ASSERT_DOUBLE_EQ(0.0, phi);
      ASSERT_DOUBLE_EQ(0.0, theta);
    }
} // namespace dnpsoup


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
      ASSERT_DOUBLE_EQ(dnpsoup::pi * 0.5, theta);

      auto c3 = Coord(1.0, 2.0, 3.0);
      auto c4 = Coord(3.0, 5.0, 9.0);
      double dist = dnpsoup::calcDistance(c3, c4);
      ASSERT_DOUBLE_EQ(7.0, dist);

      Coord c5(0.0, 0.0, 1.0);
      Coord c6(std::sqrt(0.5), std::sqrt(0.5), 2.0);
      ASSERT_TRUE(c6 - c5 == Coord(std::sqrt(0.5), std::sqrt(0.5), 1.0));
      auto [phi2, theta2] = calcAnglesWithZ(c6 - c5);
      ASSERT_DOUBLE_EQ(dnpsoup::pi * 0.25, phi2);
      ASSERT_DOUBLE_EQ(dnpsoup::pi * 0.25, theta2);

    }
} // namespace dnpsoup


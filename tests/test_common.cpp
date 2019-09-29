#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;
    using cxdbl = dnpsoup::cxdbl;

    TEST(TestDnpsoup, Commute){
      auto x = dnpsoup::spin<dnpsoup::X>(4);
      auto y = dnpsoup::spin<dnpsoup::Y>(4);
      auto z = dnpsoup::spin<dnpsoup::Z>(4);
      MatrixCD desired = cxdbl(0,1) * z;
      auto res = dnpsoup::commute(x, y);
      ASSERT_TRUE(dnpsoup::allclose(res, desired, 1e-14));
      res = dnpsoup::commute(y, z);
      ASSERT_TRUE(dnpsoup::allclose(res, cxdbl(0,1) * x, 1e-14));
      res = dnpsoup::commute(z, x);
      ASSERT_TRUE(dnpsoup::allclose(res, cxdbl(0,1) * y, 1e-14));
    }

    TEST(TestDnpsoup, ArcTan){
      double val1 = dnpsoup::atan(0.0, 0.0);
      double val2 = dnpsoup::atan(0.5, 0.5);
      double val3 = dnpsoup::atan(-0.5, 0.5);
      ASSERT_DOUBLE_EQ(0.0, val1);
      ASSERT_DOUBLE_EQ(dnpsoup::pi * 0.25, val2);
      ASSERT_DOUBLE_EQ(dnpsoup::pi * 1.75, val3);
    }
}


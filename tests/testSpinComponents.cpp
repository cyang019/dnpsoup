#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;

    TEST(TestDnpsoup, SpinX){
      MatrixCD desired2 = {{0.0, 0.5}, {0.5, 0.0}};
      MatrixCD spinX2 = dnpsoup::spin<dnpsoup::X>(2);
      ASSERT_TRUE(dnpsoup::allclose(desired2, spinX2, 1.0e-14));

      MatrixCD desired3 = {{0.0, 0.5 * std::sqrt(2), 0},
                           {0.5 * std::sqrt(2), 0, 0.5 * std::sqrt(2)},
                           {0.0, 0.5 * std::sqrt(2), 0}};
      MatrixCD spinX3 = dnpsoup::spin<dnpsoup::X>(3);
      ASSERT_TRUE(dnpsoup::allclose(desired3, spinX3, 1.0e-14));
    }

    TEST(TestDnpsoup, SpinY){
      MatrixCD desired2 = {{0.0, -dnpsoup::cxdbl(0, 0.5)}, {dnpsoup::cxdbl(0,0.5), 0.0}};
      MatrixCD spinY2 = dnpsoup::spin<dnpsoup::Y>(2);
      ASSERT_TRUE(dnpsoup::allclose(desired2, spinY2, 1.0e-14));

      MatrixCD desired3 = {{0.0, -dnpsoup::cxdbl(0, 0.5 * std::sqrt(2)), 0},
                           {dnpsoup::cxdbl(0, 0.5 * std::sqrt(2)), 0, -dnpsoup::cxdbl(0, 0.5 * std::sqrt(2))},
                           {0.0, dnpsoup::cxdbl(0, 0.5 * std::sqrt(2)), 0}};
      MatrixCD spinY3 = dnpsoup::spin<dnpsoup::Y>(3);
      ASSERT_TRUE(dnpsoup::allclose(desired3, spinY3, 1.0e-14));
    }

    TEST(TestDnpsoup, SpinZ){
      MatrixCD desired2 = {{0.5, 0.0}, {0.0, -0.5}};
      MatrixCD spinZ2 = dnpsoup::spin<dnpsoup::Z>(2);
      ASSERT_TRUE(dnpsoup::allclose(desired2, spinZ2, 1.0e-14));

      MatrixCD desired3 = dnpsoup::diagonal<dnpsoup::cxdbl>({1.0, 0.0, -1.0});
      MatrixCD spinZ3 = dnpsoup::spin<dnpsoup::Z>(3);
      ASSERT_TRUE(dnpsoup::allclose(desired3, spinZ3, 1.0e-14));
    }

    TEST(TestDnpsoup, SpinPlus){
      MatrixCD desired2 = {{0.0, 1.0}, {0.0, 0.0}};
      MatrixCD spinP2 = dnpsoup::spin<dnpsoup::P>(2);
      ASSERT_TRUE(dnpsoup::allclose(desired2, spinP2, 1.0e-14));

      MatrixCD desired3 = {{0.0, std::sqrt(2), 0},
                           {0.0, 0.0, std::sqrt(2)},
                           {0.0, 0.0, 0.0}};
      MatrixCD spinP3 = dnpsoup::spin<dnpsoup::P>(3);
      ASSERT_TRUE(dnpsoup::allclose(desired3, spinP3, 1.0e-14));
    }

    TEST(TestDnpsoup, SpinMinus){
      MatrixCD desired2 = {{0.0, 0.0}, {1.0, 0.0}};
      MatrixCD spinM2 = dnpsoup::spin<dnpsoup::M>(2);
      ASSERT_TRUE(dnpsoup::allclose(desired2, spinM2, 1.0e-14));

      MatrixCD desired3 = {{0.0, 0.0, 0.0},
                           {std::sqrt(2), 0.0, 0.0},
                           {0.0, std::sqrt(2), 0.0}};
      MatrixCD spinM3 = dnpsoup::spin<dnpsoup::M>(3);
      ASSERT_TRUE(dnpsoup::allclose(desired3, spinM3, 1.0e-14));
    }
}


#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;

    TEST(TestDnpsoup, SpinX){
      MatrixCD desired2 = {{0.0, 0.5}, {0.5, 0.0}};
      MatrixCD spinX2 = dnpsoup::spin<dnpsoup::X>(2);
      ASSERT_TRUE(dnpsoup::allclose(desired2, spinX2, 1.0e-14));
    }
}


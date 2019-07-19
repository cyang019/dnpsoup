#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;
    using cxdbl = dnpsoup::cxdbl;

    TEST(TestDnpsoup, DipolarInteraction){
      auto dipoleRR = DipolarInteraction<RotatingFrame, RotatingFrame>(1.0, 1,0, 2, 2);
      auto dipoleRL = DipolarInteraction<RotatingFrame, LabFrame>(1.0, 1,0, 2, 2);
    }
} // namespace


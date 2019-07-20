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
      auto dipoleRR = dnpsoup::DipolarInteraction<dnpsoup::RotatingFrame, dnpsoup::RotatingFrame>(
          dnpsoup::gamma_H1, dnpsoup::gamma_H1, 2, 2);
      auto dipoleRL = dnpsoup::DipolarInteraction<dnpsoup::RotatingFrame, dnpsoup::LabFrame>(
          dnpsoup::gamma_H1, dnpsoup::gamma_H1, 2, 2);

      auto p = dnpsoup::genDipoleProperty(1.0);
      auto e = dnpsoup::Euler(0,0,0);
      auto matRR = dipoleRR.genMatrix(p, e);
      std::cout << "matRR:\n" << matRR << std::endl;
      auto matRL = dipoleRL.genMatrix(p, e);
      std::cout << "matRL:\n" << matRL << std::endl;
    }
} // namespace


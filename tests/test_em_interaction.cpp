#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <vector>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;
    using cxdbl = dnpsoup::cxdbl;
    using namespace dnpsoup;
    using namespace std;

    TEST(TestDnpsoup, EMInteractionOnElectron){
      auto emr_interaction = dnpsoup::EMInteraction<dnpsoup::RotatingFrame>(
          std::vector<SpinType>({SpinType::H, SpinType::e}),
          SpinType::e);

      auto x2 = dnpsoup::kron(identity<cxdbl>(2), spin<X>(2));
      auto y2 = dnpsoup::kron(identity<cxdbl>(2), spin<Y>(2));
      double freq = 1.0e6;
      auto p90 = dnpsoup::genWaveProperty(1.0e6, 90.0, 0.0);
      auto p0 = dnpsoup::genWaveProperty(1.0e6, 0.0, 0.0);
      auto euler = dnpsoup::Euler<>(0,0,0);
      auto matY = emr_interaction.genMatrix(p90, euler);
      ASSERT_TRUE(matrix::allclose(y2*freq, matY, 1.0e-10, 1.0e-10));
      std::cout << "Irradiation on y:\n" << matY << std::endl;
      auto matX = emr_interaction.genMatrix(p0, euler);
      ASSERT_TRUE(matrix::allclose(x2*freq, matX, 1.0e-10, 1.0e-10));
      std::cout << "Irradiation on x:\n" << matX << std::endl;
    }
} // namespace



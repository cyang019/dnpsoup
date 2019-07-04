#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;
    using dnpsoup::sqrt;
    using dnpsoup::allclose;

    TEST(TestDnpsoup, ShowTensors){
      auto printTensorOneSpin = [](int n) {
        std::cout << "0, 0" << "\n";
        std::cout << dnpsoup::tensor<0,0>(n) << "\n";
        std::cout << "\n";
        std::cout << "1, 0" << "\n";
        std::cout << dnpsoup::tensor<1,0>(n) << "\n";
        std::cout << "1, -1" << "\n";
        std::cout << dnpsoup::tensor<1,-1>(n) << "\n";
        std::cout << "1, 1" << "\n";
        std::cout << dnpsoup::tensor<1,1>(n) << "\n";
        std::cout << "\n";
      };
      auto printTensorTwoSpins = [](int n1, int n2){
        std::cout << "0, 0" << "\n";
        std::cout << dnpsoup::tensor<0,0>(n1, n2) << "\n";
        std::cout << "\n";
        std::cout << "1, 0" << "\n";
        std::cout << dnpsoup::tensor<1,0>(n1, n2) << "\n";
        std::cout << "1, -1" << "\n";
        std::cout << dnpsoup::tensor<1,-1>(n1, n2) << "\n";
        std::cout << "1, 1" << "\n";
        std::cout << dnpsoup::tensor<1,1>(n1, n2) << "\n";
        std::cout << "\n";
        std::cout << "2, 0" << "\n";
        std::cout << dnpsoup::tensor<2,0>(n1, n2) << "\n";
        std::cout << "2, -1" << "\n";
        std::cout << dnpsoup::tensor<2,-1>(n1, n2) << "\n";
        std::cout << "2, 1" << "\n";
        std::cout << dnpsoup::tensor<2,1>(n1, n2) << "\n";
        std::cout << "2, -2" << "\n";
        std::cout << dnpsoup::tensor<2,-2>(n1, n2) << "\n";
        std::cout << "2, 2" << "\n";
        std::cout << dnpsoup::tensor<2,2>(n1, n2) << "\n";
        std::cout << "\n";
      };

      printTensorOneSpin(2);
      printTensorTwoSpins(2, 2);
    }

    TEST(TestDnpsoup, Tensor00){
      MatrixCD t00 = {{sqrt(2.0)*0.5, 0.0}, {0.0, sqrt(2.0)*0.5}};
      auto res = dnpsoup::tensor<0,0>(2);
      ASSERT_TRUE(allclose(t00, res, 1.0e-14));
    }

    TEST(TestDnpsoup, TensorNorm){
      MatrixCD t00 = dnpsoup::tensor<0,0>(2,2);
      MatrixCD t10 = dnpsoup::tensor<1,0>(2,2);
      MatrixCD t1n1 = dnpsoup::tensor<1,-1>(2,2);
      MatrixCD t11 = dnpsoup::tensor<1,1>(2,2);
      MatrixCD t20 = dnpsoup::tensor<2,0>(2,2);
      MatrixCD t2n1 = dnpsoup::tensor<2,-1>(2,2);
      MatrixCD t21 = dnpsoup::tensor<2,1>(2,2);
      MatrixCD t2n2 = dnpsoup::tensor<2,-2>(2,2);
      MatrixCD t22 = dnpsoup::tensor<2,2>(2,2);
      dnpsoup::cxdbl n = dnpsoup::projection(t00, t00);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t10, t10);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t11, t11);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t1n1, t1n1);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t20, t20);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t21, t21);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t2n1, t2n1);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t22, t22);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t2n2, t2n2);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(1,0), 1.0e-14));
      n = dnpsoup::projection(t10, t00);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(0,0), 1.0e-14));
      n = dnpsoup::projection(t10, t21);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(0,0), 1.0e-14));
      n = dnpsoup::projection(t22, t2n1);
      ASSERT_TRUE(dnpsoup::approxEqual(n, dnpsoup::cxdbl(0,0), 1.0e-14));
    }
}

#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>

namespace {
  using namespace std;
  using namespace dnpsoup;

  TEST(TestDnpsoup, ProjectionNorm){
    auto mat1 = kron(spin<Z>(2), identity<cxdbl>(2));
    auto norm1 = dnpsoup::projectionNorm(mat1, mat1);
    ASSERT_DOUBLE_EQ(1.0, norm1.real());
    ASSERT_DOUBLE_EQ(0.0, norm1.imag());

    mat1 = kron(spin<Z>(2), identity<cxdbl>(2));
    auto mat2 = kron(identity<cxdbl>(2), spin<Z>(2));
    auto norm2 = dnpsoup::projectionNorm(mat2, mat1);
    ASSERT_DOUBLE_EQ(0.0, norm2.real());
    ASSERT_DOUBLE_EQ(0.0, norm2.imag());
  }

}

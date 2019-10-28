#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>

namespace {
  using namespace dnpsoup;
  using namespace std;

  TEST(TestDnpsoup, SuperOperator){
    auto mat = spin<Z>(2);
    auto u = dnpsoup::exp(cxdbl(0,-1.0) * mat);
    auto u_inv = dnpsoup::exp(cxdbl(0,1.0) * mat);
    auto d_super = rotationSuperOp(u);
    auto d_super_inv = rotationSuperOp(u_inv);
    auto desired = identity<cxdbl>(d_super.nrows());
    auto result = d_super * d_super_inv;
    ASSERT_TRUE(allclose(desired, result, 1.0e-14));
  }

  TEST(TestDnpsoup, CommutationSuperOp){
    auto Iy = kron(identity<cxdbl>(2), spin<Y>(2));
    auto ham_Iy_super = commutationSuperOp(Iy);
    auto Iy_super = dnpsoup::flatten(Iy, 'c');
    auto Iy_evolved = ham_Iy_super * Iy_super;
    auto desired = zeros<cxdbl>(Iy_evolved.nrows(),1);
    ASSERT_TRUE(allclose(desired, Iy_evolved, 1.0e-14));
  }
}

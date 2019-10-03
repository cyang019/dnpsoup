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
    std::cout << d_super << std::endl;
    std::cout << d_super_inv << std::endl;
    std::cout << d_super * d_super_inv << std::endl;
  }
}

#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>

namespace {
    TEST(TestDnpsoup, fibonacci){
      auto n1 = ::dnpsoup::fibonacci(0);
      auto n2 = ::dnpsoup::fibonacci(1);
      auto n3 = ::dnpsoup::fibonacci(2);
      auto n4 = ::dnpsoup::fibonacci(3);
      auto n5 = ::dnpsoup::fibonacci(4);
      auto n6 = ::dnpsoup::fibonacci(5);
      auto n7 = ::dnpsoup::fibonacci(6);
      auto n8 = ::dnpsoup::fibonacci(7);

      ASSERT_EQ(8u, n1);
      ASSERT_EQ(13u, n2);
      ASSERT_EQ(21u, n3);
      ASSERT_EQ(34u, n4);
      ASSERT_EQ(55u, n5);
      ASSERT_EQ(89u, n6);
      ASSERT_EQ(144u, n7);
      ASSERT_EQ(233u, n8);
    }
} // namespace



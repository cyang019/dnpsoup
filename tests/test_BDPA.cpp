#include "dnpsoup.h"
#include "configure_dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;
    using cxdbl = dnpsoup::cxdbl;

    TEST(TestDnpsoup, BDPA){
      std::string bdpa_js_file = std::string(EXAMPLE_DIR) + "/BDPA_spinsys.json";
      std::ifstream spinsys_stream;
	    spinsys_stream.exceptions(std::ios::failbit | std::ios::badbit);
	    spinsys_stream.open(bdpa_js_file.c_str());
      dnpsoup::json spinsys_js;
	    spinsys_stream >> spinsys_js;
	    if(spinsys_js.find("spinsys") == spinsys_js.end()){
	    	throw std::runtime_error("Need 'spinsys' entry in the input json file.");
	    }
      std::istringstream spinsys_iss(spinsys_js["spinsys"].dump());
      dnpsoup::SpinSys spinsys;
	    spinsys_iss >> spinsys;

      ASSERT_EQ(2u, spinsys.spinCount());
      ASSERT_EQ(2u, spinsys.typeCount());
      ASSERT_EQ(4u, spinsys.observableCount());

      auto spins = spinsys.getSpins();
      auto id1 = dnpsoup::SpinId(1);
      auto id2 = dnpsoup::SpinId(2);
      ASSERT_DOUBLE_EQ(1.0e-3, spins[id1].getT1());
      ASSERT_DOUBLE_EQ(1.0e-6, spins[id1].getT2());
      ASSERT_DOUBLE_EQ(1.0, spins[id2].getT1());
      ASSERT_DOUBLE_EQ(1.0e-3, spins[id2].getT2());
    }
}

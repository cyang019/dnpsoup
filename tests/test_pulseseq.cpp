#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>


namespace {
    using PulseSeq = dnpsoup::PulseSeq;
    using Segment = dnpsoup::Segment;
    using PulseComponent = dnpsoup::PulseComponent;
    using PulsePacket = dnpsoup::PulsePacket;

    TEST(TestDnpsoup, EmptyPulseSeq){
      auto p = PulseSeq();
      std::ostringstream oss;
      oss << p;
      std::istringstream iss(oss.str());
      PulseSeq p2;
      iss >> p2;
      ASSERT_EQ(0u, p2.size());
    }

    TEST(TestDnpsoup, PulseSeqUtopia){
      auto packet = PulsePacket(50.0e6, 0);
      auto c = PulseComponent(500, 1.0e-9);
      c.setChannel(dnpsoup::SpinType::e, packet);
      auto c_delay = PulseComponent(500, 1.0e-9);
      auto s = Segment();
      s.addComponent(c);
      s.addComponent(c_delay);
      s.setRepetition(50);
      auto p = PulseSeq();
      p.addSegment(s);

      std::cout << p << std::endl;
    }
} // namespace dnpsoup



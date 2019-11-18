#include "dnpsoup.h"
#include "json.hpp"
#include "gtest/gtest.h"
#include <limits>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <sstream>


namespace {
    using PulseSeq = dnpsoup::PulseSequence;

    TEST(TestDnpsoup, EmptyPulseSeq){
      auto p = PulseSeq();
      std::ostringstream oss;
      oss << p;
      std::istringstream iss(oss.str());
      PulseSeq p2;
      iss >> p2;
      ASSERT_EQ(0u, p2.size());
    }
    TEST(TestDnpsoup, PulseSeqCW){
      dnpsoup::PulseSequence p;
      dnpsoup::pulseseq::EMRadiation emr(1.0e6, 0.0, 0.0);
      dnpsoup::pulseseq::Component c;
      c.insert({dnpsoup::SpinType::e, emr});
      p.set("emr", c);
      std::unique_ptr<dnpsoup::pulseseq::SubSequenceInterface>
        uptr_sec = std::make_unique<dnpsoup::pulseseq::Pulse>(10, "emr");
      p.set("cw", std::move(uptr_sec));
      std::vector<std::string> seq_names = { "cw" };
      p.set(seq_names);
      //std::cout << "CW Pulse Sequence:\n" << p << std::endl;
      std::ostringstream oss;
      oss << p;
      std::istringstream iss(oss.str());
      dnpsoup::PulseSequence p2;
      std::cout << "Initialize pulse sequence from stream..." << std::endl;
      iss >> p2;
      std::cout << "CW Pulse Sequence:\n" << p2 << std::endl;
    }

    TEST(TestDnpsoup, PulseSeqTopDnp){
      // Tan, Kong Ooi, Chen Yang, Ralph T. Weber, Guinevere Mathies, and Robert G. Griffin. "Time-optimized pulsed dynamic nuclear polarization." Science advances 5, no. 1 (2019): eaav6909.
      const std::string pulse_seq_str = 
        "{\n"
        "   \"increment\": 1.0e-9,"
        "   \"components\": " 
        "   {"
        "     \"emr1\": {"
        "            \"e\": { \"frequency\": 2.0e6, \"phase\": 0.0, \"offset\": 0.0 }"
        "           }"
        "   },"
        "   \"sections\": "
        "   {"
        "     \"loop\":"
        "     {"
        "       \"type\": \"Section\","
        "       \"size\": 3,"
        "       \"names\": [\"pulse_train\", \"d2\"],"
        "       \"params\": {}"
        "     },"
        "     \"pulse_train\":"
        "     {"
        "       \"type\": \"Section\","
        "       \"size\": 5,"
        "       \"names\": [\"p1\", \"d1\"],"
        "       \"params\": {}"
        "     },"
        "     \"p1\":"
        "     {"
        "       \"type\": \"Pulse\","
        "       \"size\": 50,"
        "       \"names\": [\"emr1\"],"
        "       \"params\": {}"
        "     },"
        "     \"d1\":"
        "     {"
        "       \"type\": \"Delay\","
        "       \"size\": 65,"
        "       \"names\": [],"
        "       \"params\": {}"
        "     },"
        "     \"d2\":"
        "     {"
        "       \"type\": \"Delay\","
        "       \"size\": 100,"
        "       \"names\": [],"
        "       \"params\": {}"
        "     }"
        "   },"            
        "   \"sequence\": [\"loop\"]"
        "}";
      PulseSeq top_dnp_seq;
      std::istringstream iss(pulse_seq_str);
      iss >> top_dnp_seq;
      std::cout << top_dnp_seq << std::endl;

      // call next on PulseSeq
      std::vector<dnpsoup::pulseseq::Component> emrs;
      dnpsoup::pulseseq::Component temp;
      std::uint64_t idx = 0;
      std::uint64_t sz = 0;
      while(idx < top_dnp_seq.size()){
        std::tie(temp, sz, idx) = top_dnp_seq.next();
        if(idx >= top_dnp_seq.size()) break;
        emrs.push_back(temp);
      }

      //std::size_t desired = ((50 + 65) * 5 + 100) * 3;
      std::size_t desired = ((1 + 1) * 5 + 1) * 3;
      ASSERT_EQ(desired, emrs.size());
      //std::cout << "Print EMRadiations: " << std::endl;
      //for(const auto &emr : emrs){
      //  std::cout << emr << "\n";
      //}
      //std::cout << "incremented " << emrs.size() << " components." << std::endl;
    }
} // namespace dnpsoup



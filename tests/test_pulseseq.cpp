#include "dnpsoup.h"
#include "configure_dnpsoup.h"
#include "json.hpp"
#include "gtest/gtest.h"
#include <limits>
#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <sstream>


namespace {
    using PulseSeq = dnpsoup::PulseSequence;
    TEST(TestDnpsoup, EmptyPulseSeq) {
      const std::string seq_str = "{}";
      std::istringstream iss(seq_str);
      auto seq = PulseSeq();
      iss >> seq;

      ASSERT_EQ(0u, seq.size());
      auto [comp, comp_size, idx] = seq.next();
      ASSERT_EQ(0u, comp_size);
      ASSERT_EQ(0u, idx);
      ASSERT_EQ(0u, comp.size());
    }

    TEST(TestPulseSeq, CWFromJs){
      std::string pulseseq_js_file = 
        std::string(EXAMPLE_DIR) + "/pulseseq/cw_pulse.json";

      PulseSeq pseq;
      std::ifstream pseq_ss;
      pseq_ss.open(pulseseq_js_file.c_str());
      
      pseq_ss >> pseq;
      pseq_ss.close();

      auto [comp, comp_size, idx] = pseq.next();
      std::cout << "max uint64_t:"
                << std::numeric_limits<std::uint64_t>::max()
                << std::endl;
      
      ASSERT_EQ(2000000000u, comp_size);
      ASSERT_EQ(1, comp.size());
      ASSERT_EQ(0u, idx);
    }

    TEST(TestPulseSeq, EmptyPulseSeq){
      auto p = PulseSeq();
      std::ostringstream oss;
      oss << p;
      std::istringstream iss(oss.str());
      PulseSeq p2;
      iss >> p2;
      ASSERT_EQ(0u, p2.size());
    }

    TEST(TestPulseSeq, DelayOnly){
      dnpsoup::PulseSequence seq(1.0e-8);
      // [number] x increment duration
      auto uptr_sec = std::make_unique<dnpsoup::pulseseq::Delay>(800);
      seq.set("delay", std::move(uptr_sec));
      std::vector<std::string> seq_names = { "delay" };
      seq.set(seq_names);

      std::vector<dnpsoup::pulseseq::Component> iterations;
      std::uint64_t sz = 0;
      size_t cnt = 0;
      while(sz < seq.size()){
        auto [comp, comp_sz, idx] = seq.next();
        if(comp_sz > 0){
          iterations.push_back(comp);
          ASSERT_EQ(800, comp_sz);
        }
        sz = idx;
        if (cnt > 100) {
          throw std::runtime_error("infinite loop.");
        }
        ++cnt;
      }
      ASSERT_EQ(1, iterations.size());
    }

    TEST(TestPulseSeq, PulseSeqCW){
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

    TEST(TestPulseSeq, PulseSeqTopDnp){
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
        "       \"size\": 1000,"
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
      std::vector<std::uint64_t> indices;
      auto names = top_dnp_seq.getNames("loop");
      ASSERT_EQ(2u, names.size());
      ASSERT_EQ("pulse_train", names[0]);
      ASSERT_EQ("d2", names[1]);

      names = top_dnp_seq.getNames("pulse_train");
      ASSERT_EQ(2u, names.size());
      ASSERT_EQ("p1", names[0]);
      ASSERT_EQ("d1", names[1]);

      std::cout << "\n";
      while(idx < top_dnp_seq.size()){
        std::tie(temp, sz, idx) = top_dnp_seq.next();
        if(idx >= top_dnp_seq.size()) break;
        //std::cout << "temp: " << temp << "\t\t"
        //          << "sz: " << sz << " idx: " << idx << std::endl;
        emrs.push_back(temp);
      }
      std::vector<std::uint64_t> desired_idxs = {

      };

      //std::size_t desired = ((50 + 65) * 5 + 100) * 3;
      std::size_t desired = ((1 + 1) * 5 + 1) * 3;
      ASSERT_EQ(desired, emrs.size());

      auto selector = dnpsoup::Selector(dnpsoup::ScanType::EmrLengthType, "loop");
      dnpsoup::ScanValueType new_val = dnpsoup::ScanValueType(std::uint64_t(5u));
      auto pseq2 = selector.modify(top_dnp_seq, new_val.getSizeValue());
      //std::cout << "pseq2:\n" << pseq2 << std::endl;
      std::size_t desired2 = ((1 + 1) * 5 + 1) * 5;
      idx = 0u;
      std::vector<dnpsoup::pulseseq::Component> emrs2;
      while(idx < pseq2.size()){
        std::tie(temp, sz, idx) = pseq2.next();
        if(idx >= pseq2.size()) break;
        //std::cout << "temp: " << temp << "\t\t"
        //          << "sz: " << sz << " idx: " << idx << std::endl;
        emrs2.push_back(temp);

      }
      ASSERT_EQ(desired2, emrs2.size());
      //std::cout << "Print EMRadiations: " << std::endl;
      //for(const auto &emr : emrs){
      //  std::cout << emr << "\n";
      //}
      //std::cout << "incremented " << emrs.size() << " components." << std::endl;
    }
} // namespace dnpsoup



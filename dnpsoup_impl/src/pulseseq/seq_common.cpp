#include "dnpsoup_core/pulseseq/seq_common.h"
#include "dnpsoup_core/spin_physics_components/spin.h"
#include <iostream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <string>


using namespace std;

namespace dnpsoup {
  namespace pulseseq {
    std::istream& operator>>(std::istream &is, Component &comp)
    {
      string line;
      while(getline(comp, line)){
        if(line.find_first_not_of(" \n\t") == string::npos) continue;
      
        istringstream iss(line);

        string name;
        iss >> name;
        auto t = toSpinType(name);
        auto emr = EMRadiation();
        iss >> emr.freq;
        iss >> emr.phase;
        iss >> emr.offset;
        comp[t] = emr;
      }
      return is;
    }

    std::ostream& operator<<(std::ostream &os, const Component &comp)
    {
      auto ss = os.precision();
      os << setprecision(std::numeric_limits<double>::max_digits10);
      for(const auto &emr_pair : comp){
        os << "  " << toString(emr_pair.first) << " "
           << emr_pair.second.freq << " "
           << emr_pair.second.phase << " "
           << emr_pair.second.offset << "\n";
      }
      os << setprecision(ss);
      return os;
    }
  }   // namespace pulseseq
}

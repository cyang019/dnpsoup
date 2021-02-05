#include "dnpsoup_core/spin_physics_components/spin.h"
#include "dnpsoup_core/constants.h"
#include <utility>


namespace dnpsoup {
  bool isSpinTypeStr(const std::string &s)
  {
    const std::vector<std::string> candidates = {
      "e", "E", 
      "H", "H1", "h", "h1",
      "D", "D2", "d", "d2",
      "C13", "c13",
      "N14", "n14",
      "N15", "n15",
      "O17", "o17",
      "BulkH", "HBulk"
    };
    for(const auto &name : candidates){
      if(s == name) return true;
    }
    return false;
  }

  SpinType toSpinType(const std::string &s)
  {
    if(s == "e" || s == "E"){
      return SpinType::e;
    }
    else if(s == "H" || s == "H1" || s == "h" || s == "h1"){
      return SpinType::H;
    }
    else if(s == "D" || s == "D1" || s == "d1" || s== "d"){
      return SpinType::D;
    }
    else if(s == "C13" || s == "c13"){
      return SpinType::C13;
    }
    else if(s == "N14" || s == "n14"){
      return SpinType::N14;
    }
    else if(s == "N15" || s == "n15"){
      return SpinType::N15;
    }
    else if(s == "O17" || s == "o17"){
      return SpinType::O17;
    }
    else if(s == "BulkH" || s == "HBulk"){
      return SpinType::BulkH;
    }
    return SpinType::Null;
  }

  std::string toString(const SpinType &t){
    switch(t){
      case SpinType::Null:
        return "Null";
        break;
      case SpinType::e:
        return "e";
        break;
      case SpinType::H:
        return "H";
        break;
      case SpinType::D:
        return "D";
        break;
      case SpinType::C13:
        return "C13";
        break;
      case SpinType::N14:
        return "N14";
        break;
      case SpinType::N15:
        return "N15";
        break;
      case SpinType::O17:
        return "O17";
        break;
      case SpinType::BulkH:
        return "BulkH";
        break;
      default:
        return "Unknown";
        break;
    }
    return "Unknown";
  }

  double getGyromagneticRatio(const SpinType &t)
  {
    switch(t){
      case SpinType::e:
        return beta_e;
      case SpinType::H: case SpinType::BulkH:
        return gamma_H1;
      case SpinType::D:
        return gamma_D2;
      case SpinType::C13:
        return gamma_C13;
      case SpinType::N14:
        return gamma_N14;
      case SpinType::N15:
        return gamma_N15;
      case SpinType::O17:
        return gamma_O17;
      default: return 0.0;
    }
  }

  size_t getMatrixDimension(const SpinType &t)
  {
    switch(t){
      case SpinType::e:
        return 2;
      case SpinType::H:
        return 2;
      case SpinType::BulkH:
        return 2;
      case SpinType::D:
        return 3;
      case SpinType::C13:
        return 2;
      case SpinType::N14:
        return 3;
      case SpinType::N15:
        return 2;
      case SpinType::O17:
        return 10;
      default: return 2;
    }
  }
} // namespace dnpsoup

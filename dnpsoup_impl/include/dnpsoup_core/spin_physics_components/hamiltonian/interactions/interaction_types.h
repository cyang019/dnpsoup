#ifndef DNPSOUP_INTERACTIONTYPES_H
#define DNPSOUP_INTERACTIONTYPES_H

#include "dnpsoup_core/common.h"
#include <string>
#include <unordered_set>


namespace dnpsoup {
  enum class InteractionType : int 
  {
    Scalar  = 0,
    Csa     = 1,
    Dipole  = 2,
    Shielding   = 3,
    EMR         = 10,     ///< RF or Microwave Irradiation
    Offset      = 20
  };

  inline bool operator>(const InteractionType &t1, const InteractionType &t2)
  {
    return (static_cast<int>(t1) > static_cast<int>(t2));
  }

  inline bool operator<(const InteractionType &t1, const InteractionType &t2)
  {
    return (static_cast<int>(t1) < static_cast<int>(t2));
  }

  const std::unordered_set<InteractionType, EnumClassHash> OneSpinInteractions({
      InteractionType::Csa,
      InteractionType::Shielding
      });
  const std::unordered_set<InteractionType, EnumClassHash> TwoSpinInteractions({
      InteractionType::Scalar,
      InteractionType::Dipole
      });

  inline std::string interactionName(const InteractionType &t)
  {
    std::string name;
    switch(t){
      case InteractionType::Scalar:
        name = "Scalar";
        break;
      case InteractionType::Csa:
        name = "CSA";
        break;
      case InteractionType::Dipole:
        name = "Dipole";
        break;
      case InteractionType::Shielding:
        name = "Shielding";
        break;
      case InteractionType::EMR:
        name = "EMR";
        break;
      case InteractionType::Offset:
        name = "Offset";
        break;
      default:
        name = "Default";
        break;
    }
    return name;
  }
} // namespace dnpsoup


#endif

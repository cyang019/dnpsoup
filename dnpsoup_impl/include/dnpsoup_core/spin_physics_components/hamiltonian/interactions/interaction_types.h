#ifndef DNPSOUP_INTERACTIONTYPES_H
#define DNPSOUP_INTERACTIONTYPES_H

#include <string>


namespace dnpsoup {
  enum class InteractionType {
    Scalar  = 0,
    Csa     = 1,
    Dipole  = 2,
    Shielding   = 3,
    Wave        = 10     ///< RF or Microwave Irradiation
  };

  inline std::string interactionName(const InteractionType &t)
  {
    string name;
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
      case InteractionType::Wave:
        name = "Wave";
        break;
      default:
        name = "Default";
        break;
    }
    return name;
  }
} // namespace dnpsoup


#endif

#ifndef DNPSOUP_INTERACTIONTYPES_H
#define DNPSOUP_INTERACTIONTYPES_H


namespace dnpsoup {
  enum class InteractionType {
    Scalar  = 0,
    Csa     = 1,
    Dipole  = 2,
    Shielding   = 3,
    Wave        = 10     ///< RF or Microwave Irradiation
  };
} // namespace dnpsoup


#endif

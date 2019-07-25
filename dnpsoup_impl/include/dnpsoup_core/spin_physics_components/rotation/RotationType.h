#ifndef DNPSOUP_ROTATIONTYPE_H
#define DNPSOUP_ROTATIONTYPE_H

#include <type_traits>

namespace dnpsoup {
  class ActiveRotation {};

  class PassiveRotation {};

  template<typename T> struct is_rotation_type : std::false_type {};
  template<> struct is_rotation_type<ActiveRotation> : std::true_type {};
  template<> struct is_rotation_type<PassiveRotation> : std::true_type {};
} // namespace dnpsoup


#endif

#ifndef DNPSOUP_FRAMETYPE_H
#define DNPSOUP_FRAMETYPE_H

#include <type_traits>

namespace dnpsoup {
  class RotatingFrame {};
  class LabFrame {};

  template<typename T> struct is_frame_type : std::false_type {};
  template<> struct is_frame_type<RotatingFrame> : std::true_type {};
  template<> struct is_frame_type<LabFrame> : std::true_type {};
} // namespace dnpsoup

#endif

#ifndef DNPSOUP_FRAMETYPE_H
#define DNPSOUP_FRAMETYPE_H

#include <type_traits>

namespace dnpsoup {
  class RotatingFrame : public FrameType {};
  class LabFrame : public FrameType {};

  template<typename T> struct is_frame_type : std::false_type {};
  template<typename T> struct is_frame_type<RotatingFrame> : std::true_type {};
  template<typename T> struct is_frame_type<LabFrame> : std::true_type {};
} // namespace dnpsoup

#endif

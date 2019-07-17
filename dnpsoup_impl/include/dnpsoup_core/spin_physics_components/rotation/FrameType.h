#ifndef DNPSOUP_FRAMETYPE_H
#define DNPSOUP_FRAMETYPE_H

namespace dnpsoup {
  /// rotated frame or lab frame
  enum class FrameType : int {
    Lab = 0,
    Rotate = 1
  };
} // namespace dnpsoup

#endif

#include "dnpsoup_core/common.h"
#include <algorithm>
#include <cmath>

namespace dnpsoup {
  std::int64_t genUniqueInt(std::int64_t val1, std::int64_t val2)
  {
    if(val1 < val2) std::swap(val1, val2);
    std::int64_t result = (val1 + val2) * (val1 + val2 + 1)/2 + val2;
    return result;
  }
} // namespace dnpsoup

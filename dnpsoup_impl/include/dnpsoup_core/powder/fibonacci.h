#ifndef DNPSOUP_FIBONACCI_H
#define DNPSOUP_FIBONACCI_H

#include <cstdint>


namespace dnpsoup {
  inline std::uint64_t fibonacci(std::uint64_t m)
  {
    switch(m){
      case 0:
        return 8u;
        break;
      case 1:
        return 13u;
        break;
      default:
        return fibonacci(m-1) + fibonacci(m-2);
        break;
    }
    return 0u;
  }
} // namespace dnpsoup

#endif

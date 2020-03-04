#ifndef DNPSOUP_SCAN_H
#define DNPSOUP_SCAN_H

#include "dnpsoup_core/experiment/DnpRunner.h"


namespace dnpsoup {
  enum class ScanType : int
  {
    FieldType = 0,
    EmOffsetType = 1,

    GammaB1Type = 10,
    PhaseType = 11,
    LengthType = 12,

    EulerType = 100,

    T1Type = 200,
    T2Type = 201
  };


}


#endif

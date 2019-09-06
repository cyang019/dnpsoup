#ifndef DNPSOUP_RHO_EQ_H
#define DNPSOUP_RHO_EQ_H

#include "dnpsoup_core/common.h"


namespace dnpsoup {
  MatrixCxDbl genRhoEq(const MatrixCxDbl &mat);

  std::tuple<MatrixDbl, MatrixCxDbl> diagonalizeMat(const MatrixCxDbl &mat);
} // namespace dnpsoup

#endif

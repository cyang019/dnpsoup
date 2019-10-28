#include "dnpsoup_core/common.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace std;

namespace dnpsoup {
  std::int64_t genUniqueInt(std::int64_t val1, std::int64_t val2)
  {
    if(val1 < val2) std::swap(val1, val2);
    std::int64_t result = (val1 + val2) * (val1 + val2 + 1)/2 + val2;
    return result;
  }

  MatrixCxDbl diag_exp(const MatrixCxDbl &mat)
  {
#ifndef NDEBUG
    // check if offdiagonal elements larger than eps
      for(size_t i = 0; i < mat.ncols(); ++i){
        for(size_t j = 0; j < mat.nrows(); ++j){
          if(i == j) continue; ///< diagonal element
          if(abs(mat(i,j)) > eps)
            throw runtime_error("Matrix contains non-diagonal element(s), cannot use diag_exp()");
        }
      }
#endif
      auto result = zeros<cxdbl>(mat.nrows(), mat.ncols());
      for(size_t i = 0; i < result.nrows(); ++i){
        result(i,i) = std::exp(mat(i,i));
      }
      return result;
  }

  MatrixDbl diag_exp(const MatrixDbl &mat)
  {
#ifndef NDEBUG
    // check if offdiagonal elements larger than eps
      for(size_t i = 0; i < mat.ncols(); ++i){
        for(size_t j = 0; j < mat.nrows(); ++j){
          if(i == j) continue; ///< diagonal element
          if(abs(mat(i,j)) > eps)
            throw runtime_error("Matrix contains non-diagonal element(s), cannot use diag_exp()");
        }
      }
#endif
      auto result = zeros<double>(mat.nrows(), mat.ncols());
      for(size_t i = 0; i < result.nrows(); ++i){
        result(i,i) = std::exp(mat(i,i));
      }
      return result;
  }
} // namespace dnpsoup

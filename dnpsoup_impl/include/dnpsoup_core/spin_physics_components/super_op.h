#ifndef DNPSOUP_SUPER_OP_H
#define DNPSOUP_SUPER_OP_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/errors.h"


namespace dnpsoup {
  template<typename T>
  matrix::Matrix<T> commutationSuperOp(const matrix::Matrix<T> &mat);

  template<typename T>
  matrix::Matrix<T> rotationSuperOp(const matrix::Matrix<T> &mat);

  template<typename T>
  matrix::Matrix<T> densitySuperOp(const matrix::Matrix<T> &rho);
} // namespace dnpsoup

namespace dnpsoup {
  template<typename T>
  inline
  matrix::Matrix<T> commutationSuperOp(const matrix::Matrix<T> &mat)
  {
    auto nrows = mat.nrows();
#ifndef NDEBUG
    auto ncols = mat.ncols();
    if(nrows != ncols)
      throw SizeMismatchError("need a square matrix for commutationSuperOp()");
#endif
    const auto i = identity<T>(nrows);
    matrix::Matrix<T> res = kron(mat, i) - kron(i, mat.t());
    return res;
  }

  template<typename T>
  inline
  matrix::Matrix<T> rotationSuperOp(const matrix::Matrix<T> &mat)
  {
    return kron(mat, mat.adjoint());
  }

  template<typename T>
  inline
  matrix::Matrix<T> densitySuperOp(const matrix::Matrix<T> &rho)
  {
    matrix::Matrix<T> res(rho.nelements(), 1);
    const auto nrows = rho.nrows();
    const auto ncols = rho.ncols();
    for(size_t i = 0; i < nrows; ++i){
      for(size_t j = 0; j < ncols; ++j){
        res(i * ncols + j, 1) = rho(i,j);
      }
    }
    return res;
  }
} // namespace dnpsoup

#endif

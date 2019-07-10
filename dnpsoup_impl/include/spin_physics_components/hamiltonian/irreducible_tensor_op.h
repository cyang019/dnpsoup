#ifndef DNPSOUP_IRREDUCIBLE_TENSOR_OP_H
#define DNPSOUP_IRREDUCIBLE_TENSOR_OP_H

#include "common.h"
#include "spin_physics_components/spin.h"

namespace dnpsoup {
  /// single spin tensors
  template<int l, int m>
  MatrixCxDbl tensor(size_t n);

  /// two spin tensors
  template<int l, int m>
  MatrixCxDbl tensor(size_t n1, size_t n2);
} // namespace dnpsoup

#include "spin_physics_components/hamiltonian/irreducible_tensor_op_impl.hpp"

#endif

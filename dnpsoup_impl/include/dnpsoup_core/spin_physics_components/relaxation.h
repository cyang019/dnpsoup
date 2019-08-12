#ifndef DNPSOUP_RELAXATION_H
#define DNPSOUP_RELAXATION_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/spin.h"


namespace dnpsoup {
  MatrixCxDbl secularRelaxationSuperOp(const MatrixCxDbl &op);
  MatrixCxDbl secularRelaxationSuperOp(const MatrixCxDbl &op, const MatrixCxDbl &eigen_vec);

  MatrixCxDbl t1SuperOp(
      double t1_inv,
      const MatrixCxDbl &eigen_vec,
      const MatrixCxDbl &x, 
      const MatrixCxDbl &y);

  MatrixCxDbl t2SuperOp(
      double t2_inv,
      const MatrixCxDbl &eigen_vec,
      const MatrixCxDbl &z);
} // namespace dnpsoup

#endif

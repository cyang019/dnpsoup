#include "dnpsoup_core/spin_physics_components/relaxation.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"


namespace dnpsoup {
  MatrixCxDbl secularRelaxationSuperOp(const MatrixCxDbl &op)
  {
    auto Aq = commutationSuperOp(op);
    auto Anq = commutationSuperOp(op.adjoint());
    return Anq * Aq;
  }

  MatrixCxDbl secularRelaxationSuperOp(const MatrixCxDbl &op, const MatrixCxDbl &eigen_vec)
  {
    auto Aq = commutationSuperOp(eigen_vec.adjoint() * op * eigen_vec);
    auto Anq = commutationSuperOp(eigen_vec.adjoint() * op.adjoint() * eigen_vec);
    return Anq * Aq;
  }

  MatrixCxDbl t1SuperOp(
      double t1_inv,
      const MatrixCxDbl &eigen_vec,
      const MatrixCxDbl &x, 
      const MatrixCxDbl &y)
  {
    auto x_rotated = eigen_vec.adjoint() * x * eigen_vec;
    auto y_rotated = eigen_vec.adjoint() * y * eigen_vec;

    auto part_x = secularRelaxationSuperOp(x_rotated);
    auto part_y = secularRelaxationSuperOp(y_rotated);

    return 0.5 * t1_inv * (part_x + part_y);
  }

  MatrixCxDbl t2SuperOp(
      double t2_inv,
      const MatrixCxDbl &eigen_vec,
      const MatrixCxDbl &z)
  {
    auto z_rotated = eigen_vec.adjoint() * z * eigen_vec;
    return t2_inv * secularRelaxationSuperOp(z_rotated);
  }
} // namespace dnpsoup

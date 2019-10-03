#include "dnpsoup_core/spin_physics_components/hamiltonian/rho_eq.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/constants.h"


namespace dnpsoup {
  /// @param hamiltonian in Hz
  MatrixCxDbl genRhoEq(const MatrixCxDbl &hamiltonian, double temperature)
  {
#ifndef NDEBUG
    if(hamiltonian.nrows() != hamiltonian.ncols()){
      throw SizeMismatchError("Need square matrices for matrix exponential.");
    }
#endif

    const double prefix = -dnpsoup::h / (dnpsoup::kb * temperature);
    auto energy_mat = hamiltonian * prefix;
    const double n1 = matrix::norm1(energy_mat);
    if(n1 > 1.0 - dnpsoup::eps){
      MatrixCxDbl rho_pre = ::dnpsoup::exp(energy_mat);
      return rho_pre * (1.0 / trace(rho_pre));
    }

    // trace(Identity) = n
    MatrixCxDbl rho_pre = expMinusIdentity(energy_mat);
    const double n = static_cast<double>(energy_mat.nrows());
    return rho_pre * (1.0 / (n + trace(rho_pre)));
  }

  std::tuple<MatrixDbl, MatrixCxDbl> diagonalizeMat(const MatrixCxDbl &mat)
  {
    auto [eigen_vals, eigen_vecs] = matrix::eigenSys<>(mat);
    return std::make_tuple(eigen_vals, eigen_vecs);
  }
} // namespace dnpsoup

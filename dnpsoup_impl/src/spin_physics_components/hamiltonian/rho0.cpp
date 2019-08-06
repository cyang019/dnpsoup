#include "dnpsoup_core/spin_physics_components/hamiltonian/rho0.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/constants.h"


namespace dnpsoup {
  /// @param hamiltonian in Hz
  MatrixCxDbl genRho0(const MatrixCxDbl &hamiltonian, double temperature)
  {
#ifndef NDEBUG
    if(hamiltonian.nrows() != hamiltonian.ncols()){
      throw SizeMismatchError("Need square matrices for matrix exponential.");
    }
#endif

    const double n1 = matrix::norm1(hamiltonian);
    const double prefix = -dnpsoup::h / (dnpsoup::kb * temperature);
    if(n1 > 1.0 - dnpsoup::eps){
      MatrixCxDbl rho_pre = exp(prefix * hamiltonian);
      return rho_pre * (1.0 / trace(rho_pre));
    }

    // trace(Identity) = n
    MatrixCxDbl rho_pre = expMinusIdentity(prefix * hamiltonian);
    const double n = static_cast<double>(hamiltonian.nrows());
    return rho_pre * (1.0 / (n + trace(rho_pre)));
  }
} // namespace dnpsoup

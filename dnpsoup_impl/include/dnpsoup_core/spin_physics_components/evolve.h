#ifndef DNPSOUP_EVOLVE_H
#define DNPSOUP_EVOLVE_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/rho_eq.h"
#include "dnpsoup_core/spinsys/RelaxationPacket.h"
#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/experiment/hardware.h"
#include <vector>
#include <tuple>
#include <cstdint>


namespace dnpsoup {
  /// (iH_super + Gamma_super) * rho_deq_super = Gamma * rho_ss_super
  /// @returns rho_deq_super
  MatrixCxDbl calcRhoDynamicEq(
      const MatrixCxDbl &h_super, 
      const MatrixCxDbl &gamma_super,
      const MatrixCxDbl &rho_ss_super);

  /// exp(-Lambda * dt * cnt)
  inline MatrixCxDbl calcExpEvolve(
      const MatrixCxDbl &super_op,
      double t_atomic,
      size_t cnt)
  {
#ifdef DNPSOUP_VERBOSE
    std::cout << "calcExpEvolve: " << "inc: " << t_atomic 
                                   << "  cnt " << cnt << std::endl;
#endif
    auto scaling_factor = ::dnpsoup::exp(cxdbl(-1.0 * t_atomic, 0) * super_op);
#ifdef DNPSOUP_VERBOSE
    std::cout << "scaling_factor shape: " 
              << scaling_factor.nrows() << ", " 
              << scaling_factor.ncols() << "\n";
#endif
    scaling_factor = ::dnpsoup::pow(scaling_factor, cnt);
#ifdef DNPSOUP_VERBOSE
    std::cout << "scaling_factor shape: " 
              << scaling_factor.nrows() << ", " 
              << scaling_factor.ncols() << "\n";
#endif
    return scaling_factor;
  }

  inline MatrixCxDbl evolveRho(
      const MatrixCxDbl &rho_prev_super,
      const MatrixCxDbl &rho_eq_super,
      const MatrixCxDbl &scaling_factor)
  {
    const auto c1 = rho_prev_super - rho_eq_super;
    return scaling_factor * c1 + rho_eq_super;
  }

  MatrixCxDbl calcGammaSuper(
      const MatrixCxDbl &rho_ss,
      const std::vector<RelaxationPacket> &rpackets);

  /// @returns h_super, gamma_super_internal, rho_eq_super
  std::tuple<MatrixCxDbl, MatrixCxDbl, MatrixCxDbl> calcSuperOpsForMasterEq(
      const MatrixCxDbl &ham,
      const MatrixCxDbl &ham_lab,
      const MatrixCxDbl &rotate_mat_super,
      const MatrixCxDbl &rotate_mat_super_inv,
      const std::vector<RelaxationPacket> &rpackets,
      double temperature);

  /// @returns rotate_mat_super, rotate_mat_super_inv
  std::pair<MatrixCxDbl, MatrixCxDbl> calcRotationSuperOps(
      const MatrixCxDbl &ham_offset,
      const Gyrotron &g,
      double dt,
      std::uint64_t cnt
      );


  inline MatrixCxDbl calcLambdaSuper(
      const MatrixCxDbl &h_super, const MatrixCxDbl &gamma_super)
  {
    return std::complex<double>(0, 1.0) * h_super + gamma_super;
  }

  MatrixCxDbl evolve(
      const MatrixCxDbl &rho_prev_super,
      const MatrixCxDbl &rho_eq_super,
      const MatrixCxDbl &scaling_factor,    ///< exp(-L*dt)
      const MatrixCxDbl &rotate_mat_super,
      const MatrixCxDbl &rotate_mat_super_inv
      );

  MatrixCxDbl evolveMASCnstEmr(
      const MatrixCxDbl &rho_prev_super,
      double mas_frequency,
      const pulseseq::Component &comp,
      const PacketCollection &packets,
      const std::vector<RelaxationPacket> &rpackets,
      const MatrixCxDbl ham_offset,
      const Euler<> &spin_sys_euler,
      Euler<> magic_angle,
      const Gyrotron &g,
      double inc, 
      std::uint64_t cnt,
      std::uint64_t mas_cnt,
      std::uint64_t total_rotor_cnt,
      double temperature);

  std::vector<std::pair<double, double>> evolveMASCnstEmr(
      const MatrixCxDbl &rho_prev_super,
      const MatrixCxDbl &detect_op_super,
      double t0,
      double result_ref,
      double mas_frequency,
      const pulseseq::Component &comp,
      const PacketCollection &packets,
      const std::vector<RelaxationPacket> &rpackets,
      const MatrixCxDbl ham_offset,
      const Euler<> &spin_sys_euler,
      Euler<> magic_angle,
      const Gyrotron &g,
      double inc, 
      std::uint64_t cnt,
      std::uint64_t mas_cnt,
      std::uint64_t total_rotor_cnt,
      double temperature);

} // namespace dnpsoup

#endif

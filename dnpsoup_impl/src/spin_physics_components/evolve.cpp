#include "dnpsoup_core/spin_physics_components/evolve.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
#include "dnpsoup_core/spin_physics_components/EvolutionCache.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spinsys/RelaxationPacket.h"
#include "dnpsoup_core/common.h"
#include <complex>

using namespace std;

namespace dnpsoup {
  MatrixCxDbl calcRhoDynamicEq(
      const MatrixCxDbl &h_super, 
      const MatrixCxDbl &gamma_super,
      const MatrixCxDbl &rho_ss_super)
  {
    const auto super_op = complex<double>(0, 1.0) * h_super + gamma_super;
    const auto rho_right = gamma_super * rho_ss_super;
    auto [rho_eq_super, status] = matrix::lstsq(super_op, rho_right);
#ifndef NDEBUG
      if(status != 0){
        string err_msg = "lstsq error: ";
        if(status < 0){
          err_msg = "The " + std::to_string(-status) + "-th argument had an illegal value.";
        } else {
          err_msg = "The algorithm for computing the SVD failed to converge;\n";
          err_msg += std::to_string(status) + "off-diagonal elements of an intermediate"
            + "bidiagonal form did not converge to zero.";
        }
        throw CalculationError(err_msg);
      }
#endif
    return rho_eq_super;
  }


  MatrixCxDbl calcGammaSuper(
      const MatrixCxDbl &rho_ss,
      const std::vector<RelaxationPacket> &rpackets)
  {
    auto [eigenvals, eigenvec] = diagonalizeMat(rho_ss);
    const size_t sz = rho_ss.nrows() * rho_ss.ncols();
    auto t1_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
    auto t2_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
    for(const auto &rpacket : rpackets){
      t1_superop += rpacket.genSuperOpT1(eigenvec);
      t2_superop += rpacket.genSuperOpT2(eigenvec);
    }
    auto gamma_super = t1_superop + t2_superop;
    return gamma_super;
  }

  std::tuple<MatrixCxDbl, MatrixCxDbl, MatrixCxDbl> calcSuperOpsForMasterEq(
      const MatrixCxDbl &ham,
      const MatrixCxDbl &ham_lab,
      const MatrixCxDbl &rotate_mat_super,
      const MatrixCxDbl &rotate_mat_super_inv,
      const std::vector<RelaxationPacket> &rpackets,
      double temperature)
  {
    // static state thermo equilibrium
      auto rho_ss = genRhoEq(ham_lab, temperature);

      auto rho_ss_super = ::dnpsoup::flatten(rho_ss, 'c');
      if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
        rho_ss_super = rotate_mat_super * rho_ss_super;
      }
      auto gamma_super = calcGammaSuper(rho_ss, rpackets);
      MatrixCxDbl gamma_super_int;
      if(rotate_mat_super.nrows() == 0 || rotate_mat_super_inv.nrows() == 0){
        gamma_super_int = gamma_super;
      } else {
        gamma_super_int = rotate_mat_super * gamma_super * rotate_mat_super_inv;
      }
      auto h_super = commutationSuperOp(ham);
      auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;
      auto rho_eq_super = calcRhoDynamicEq(h_super, gamma_super_int, rho_ss_super);
      return std::make_tuple(std::move(h_super), 
                             std::move(gamma_super_int),
                             std::move(rho_eq_super));
  }

  std::pair<MatrixCxDbl, MatrixCxDbl> calcRotationSuperOps(
      const MatrixCxDbl &ham_offset,
      const Gyrotron &g,
      double dt,
      std::uint64_t cnt
      )
  {
    // to determine if we need rotate_mat
    const double cycle = dt * static_cast<double>(cnt) * g.em_frequency;
    double cycle_int_part;
    double cycle_fraction_part = modf(cycle, &cycle_int_part);
    if (std::abs(cycle_fraction_part) < 1.0e-4) {
      // placeholders
      MatrixCxDbl rotate_mat_super, rotate_mat_super_inv;
      return make_pair(rotate_mat_super, rotate_mat_super_inv);
    } else {
      MatrixCxDbl rotate_mat = dnpsoup::exp((cxdbl(0,-2) * dt * pi) * ham_offset);
      MatrixCxDbl rotate_mat_inv = dnpsoup::exp((cxdbl(0, 2) * dt * pi) * ham_offset);
      MatrixCxDbl rotate_mat_super = rotationSuperOp(rotate_mat);
      MatrixCxDbl rotate_mat_super_inv = rotationSuperOp(rotate_mat_inv);
      return make_pair(rotate_mat_super, rotate_mat_super_inv);
    }
  }

  MatrixCxDbl evolve(
      const MatrixCxDbl &rho_prev_super,
      const MatrixCxDbl &rho_eq_super,
      const MatrixCxDbl &scaling_factor,
      const MatrixCxDbl &rotate_mat_super,
      const MatrixCxDbl &rotate_mat_super_inv
      )
  {
    auto rho_super = rho_prev_super;
    if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
      rho_super = rotate_mat_super * rho_prev_super;
    }

    rho_super = evolveRho(rho_super, rho_eq_super, scaling_factor);
    if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
      rho_super = rotate_mat_super_inv * rho_super;
    }
    //MatrixCxDbl rho_post(rho_prev.nrows(), rho_prev.ncols());

    //const auto nrows = rho_post.nrows();
    //const auto ncols = rho_post.ncols();
    //for(size_t i = 0; i < nrows; ++i){
    //  for(size_t j = 0; j < ncols; ++j){
    //    rho_post(i, j) = rho_super(j + i * ncols, 0);
    //  }
    //}
    return rho_super;
  }

  MatrixCxDbl evolveMASCnstEmr(
      const MatrixCxDbl &rho_prev_super,
      double mas_frequency,
      const pulseseq::Component &comp,
      const PacketCollection &packets,
      const std::vector<RelaxationPacket> &rpackets,
      const MatrixCxDbl ham_offset,
      const Euler<> &spin_sys_euler,
      Euler<> mas_angle,
      const Gyrotron &g,
      double inc, 
      std::uint64_t cnt,
      std::uint64_t mas_inc_cnt,
      std::uint64_t total_rotor_cnt,
      double temperature)
  {
    auto cache = EvolutionCache(total_rotor_cnt, 1);
    auto [rotate_mat_super, rotate_mat_super_inv] = 
      calcRotationSuperOps(ham_offset, g, inc, mas_inc_cnt);
    
    double t = 0.0;
    uint64_t rotor_cnt = 0u;

    MatrixCxDbl rho_evolve_super = rho_prev_super;
    
    while(cnt > 0){
      auto temp_euler = spin_sys_euler * mas_angle;
      if(cnt >= mas_inc_cnt){
        if(rotor_cnt < total_rotor_cnt) { // need to save to cache
          auto ham = packets.genMatrix(temp_euler);
          auto ham_lab = ham + ham_offset;
          auto [h_super, gamma_super_internal, rho_eq_super] =
            calcSuperOpsForMasterEq(ham, ham_lab, 
                rotate_mat_super, rotate_mat_super_inv, 
                rpackets, temperature);
          auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
          auto scaling_factor = calcExpEvolve(super_op, inc, mas_inc_cnt);
          rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
              scaling_factor, rotate_mat_super, rotate_mat_super_inv);
          cache.saveCache(comp, EvolutionCacheElement(scaling_factor, rho_eq_super)); 
        }
        else{ // retrieve from cache
          const auto idx = rotor_cnt % total_rotor_cnt;
          auto cache_elem = cache.getCache(0, idx);
          rho_evolve_super = evolve(rho_evolve_super, 
              cache_elem.rho_inf_eq, cache_elem.scaling_factor,
              rotate_mat_super, rotate_mat_super_inv);
        }
        
        ++rotor_cnt;
        cnt -= mas_inc_cnt;
        t += static_cast<double>(mas_inc_cnt) * inc;

        mas_angle.gamma(t * mas_frequency * 2.0 * pi);
      }
      else {    /// if pulse length is not a integer multiple of rotor period.
        auto [rotate_mat_super, rotate_at_super_inv] = 
          calcRotationSuperOps(ham_offset, g, inc, cnt);

        auto ham = packets.genMatrix(temp_euler);
        auto ham_lab = ham + ham_offset;
        auto [h_super, gamma_super_internal, rho_eq_super] =
          calcSuperOpsForMasterEq(ham, ham_lab, 
              rotate_mat_super, rotate_mat_super_inv, 
              rpackets, temperature);
        auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
        auto scaling_factor = calcExpEvolve(super_op, inc, cnt);
        rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
            scaling_factor, rotate_mat_super, rotate_mat_super_inv);
        cnt = 0;
      }
    }
    return rho_evolve_super;
  }

  std::vector<std::pair<double, double>> evolveMASCnstEmr(
      const MatrixCxDbl &rho_prev_super,
      const MatrixCxDbl &acq_mat_super,
      double t0,
      double result_ref,
      double mas_frequency,
      const pulseseq::Component &comp,
      const PacketCollection &packets,
      const std::vector<RelaxationPacket> &rpackets,
      const MatrixCxDbl ham_offset,
      const Euler<> &spin_sys_euler,
      Euler<> mas_angle,
      const Gyrotron &g,
      double inc, 
      std::uint64_t cnt,
      std::uint64_t mas_inc_cnt,
      std::uint64_t total_rotor_cnt,
      double temperature)
  {
    auto cache = EvolutionCache(total_rotor_cnt, 1);
    auto [rotate_mat_super, rotate_mat_super_inv] = 
      calcRotationSuperOps(ham_offset, g, inc, mas_inc_cnt);
    
    double t = 0.0;
    uint64_t rotor_cnt = 0u;

    MatrixCxDbl rho_evolve_super = rho_prev_super;
    
    std::vector<std::pair<double, double>> results;
    while(cnt > 0){
      auto temp_euler = spin_sys_euler * mas_angle;
      if(cnt >= mas_inc_cnt){
        if(rotor_cnt < total_rotor_cnt) { // need to save to cache
          auto ham = packets.genMatrix(temp_euler);
          auto ham_lab = ham + ham_offset;
          auto [h_super, gamma_super_internal, rho_eq_super] =
            calcSuperOpsForMasterEq(ham, ham_lab, 
                rotate_mat_super, rotate_mat_super_inv, 
                rpackets, temperature);
          auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
          auto scaling_factor = calcExpEvolve(super_op, inc, mas_inc_cnt);
          rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
              scaling_factor, rotate_mat_super, rotate_mat_super_inv);
          cache.saveCache(comp, EvolutionCacheElement(scaling_factor, rho_eq_super)); 
          double result = ::dnpsoup::projectionNorm(
              rho_evolve_super, acq_mat_super).real();
          results.push_back(make_pair(t0 + t, result/result_ref));
        }
        else{ // retrieve from cache
          const auto idx = rotor_cnt % total_rotor_cnt;
          auto cache_elem = cache.getCache(0, idx);
          rho_evolve_super = evolve(rho_evolve_super, 
              cache_elem.rho_inf_eq, cache_elem.scaling_factor,
              rotate_mat_super, rotate_mat_super_inv);
          double result = ::dnpsoup::projectionNorm(
              rho_evolve_super, acq_mat_super).real();
          results.push_back(make_pair(t0 + t, result/result_ref));
        }
        
        ++rotor_cnt;
        cnt -= mas_inc_cnt;
        t += static_cast<double>(mas_inc_cnt) * inc;

        mas_angle.gamma(t * mas_frequency * 2.0 * pi);
      }
      else {    /// if pulse length is not a integer multiple of rotor period.
        auto [rotate_mat_super, rotate_at_super_inv] = 
          calcRotationSuperOps(ham_offset, g, inc, cnt);

        auto ham = packets.genMatrix(temp_euler);
        auto ham_lab = ham + ham_offset;
        auto [h_super, gamma_super_internal, rho_eq_super] =
          calcSuperOpsForMasterEq(ham, ham_lab, 
              rotate_mat_super, rotate_mat_super_inv, 
              rpackets, temperature);
        auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
        auto scaling_factor = calcExpEvolve(super_op, inc, cnt);
        rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
            scaling_factor, rotate_mat_super, rotate_mat_super_inv);
        cnt = 0;
        double result = ::dnpsoup::projectionNorm(
            rho_evolve_super, acq_mat_super).real();
        results.push_back(make_pair(t0 + t, result/result_ref));
      }
    }
    return results;
  }
} // namespace dnpsoup

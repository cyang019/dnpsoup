#include "dnpsoup_core/spin_physics_components/evolve.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/spinsys/RelaxationPacket.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/constants.h"
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
      MatrixCxDbl h_super = commutationSuperOp(ham);
      //auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;
      MatrixCxDbl rho_eq_super = calcRhoDynamicEq(h_super, gamma_super_int, rho_ss_super);
      //return std::make_tuple(h_super, gamma_super_int, rho_eq_super);
      return std::make_tuple(std::move(h_super), 
                             std::move(gamma_super_int),
                             std::move(rho_eq_super));
  }

  /// @returns h_super, gamma_super_internal, rho_eq_super
  std::tuple<MatrixCxDbl, MatrixCxDbl, MatrixCxDbl> calcSuperOpsForMasterEq(
      const MatrixCxDbl &ham,
      const MatrixCxDbl &ham_lab,
      const std::vector<RelaxationPacket> &rpackets,
      double temperature)
  {
    // static state thermo equilibrium
      auto rho_ss = genRhoEq(ham_lab, temperature);

      auto rho_ss_super = ::dnpsoup::flatten(rho_ss, 'c');
      auto gamma_super = calcGammaSuper(rho_ss, rpackets);
      auto h_super = commutationSuperOp(ham);
      //auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;
      auto rho_eq_super = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
      //return std::make_tuple(h_super, gamma_super, rho_eq_super);
      return std::make_tuple(std::move(h_super), 
                             std::move(gamma_super),
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
    return rho_super;
  }

  MatrixCxDbl evolveMASCnstEmr(
      const MatrixCxDbl &rho_prev_super,
      double mas_frequency,
      //const pulseseq::Component &comp,
      const PacketCollection &packets,
      const std::vector<RelaxationPacket> &rpackets,
      const MatrixCxDbl &ham_offset,
      const Euler<> &spin_sys_euler,
      Euler<> mas_angle,
      //const Gyrotron &g,
      double inc, 
      std::uint64_t cnt,
      std::uint64_t mas_inc_cnt,
      std::uint64_t total_rotor_cnt,
      double temperature)
  {
    /// 2d cache
    //auto cache = EvolutionCache(total_rotor_cnt, 1);
    vector<MatrixCxDbl> cache_rho;
    vector<MatrixCxDbl> cache_scaling_factor;
    cache_rho.reserve(total_rotor_cnt);
    cache_scaling_factor.reserve(total_rotor_cnt);
    
    //auto [rotate_mat_super, rotate_mat_super_inv] = 
    //  calcRotationSuperOps(ham_offset, g, inc, mas_inc_cnt);
    const double mas_inc = inc * static_cast<double>(mas_inc_cnt);
    const double gamma0 = mas_angle.gamma();
    
    double t = 0.0;
    uint64_t rotor_cnt = 0u;

    MatrixCxDbl rho_evolve_super = rho_prev_super;
    size_t gamma_step_size = static_cast<size_t>(std::round(GAMMA_STEP_MIN/mas_inc));
    gamma_step_size += (gamma_step_size == 0);
    MatrixCxDbl gamma_super;
    
    while(cnt > 0){
      auto temp_euler = mas_angle * spin_sys_euler;
      if(cnt >= mas_inc_cnt){
        t += static_cast<double>(mas_inc_cnt) * inc;
        if(rotor_cnt < total_rotor_cnt) { // need to save to cache
          const auto ham = packets.genMatrix(temp_euler);
          const auto ham_lab = ham + ham_offset;
          const auto rho_ss = genRhoEq(ham_lab, temperature);
          auto rho_ss_super = ::dnpsoup::flatten(rho_ss, 'c');
          //if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
          //  rho_ss_super = rotate_mat_super * rho_ss_super;
          //}
          if(rotor_cnt % gamma_step_size == 0) {
            gamma_super = calcGammaSuper(rho_ss, rpackets);
            //if (rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0) {
            //  gamma_super = rotate_mat_super * gamma_super * rotate_mat_super_inv;
            //}
          }
          const auto h_super = commutationSuperOp(ham);
          //auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;
          const auto super_op = calcLambdaSuper(h_super, gamma_super);
          const auto rho_eq_super = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
          const auto scaling_factor = calcExpEvolve(super_op, inc, mas_inc_cnt);
          //rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
          //    scaling_factor, rotate_mat_super, rotate_mat_super_inv);
          rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, scaling_factor);
          cache_rho.push_back(rho_eq_super);
          cache_scaling_factor.push_back(scaling_factor);
          //cache.saveCache(
          //    comp, 
          //    EvolutionCacheElement(
          //      std::move(scaling_factor), std::move(rho_eq_super))); 
        }
        else{ // retrieve from cache
          t += static_cast<double>(cnt) * inc;
          const auto idx = rotor_cnt % total_rotor_cnt;
          rho_evolve_super = evolve(rho_evolve_super, 
              cache_rho[idx], cache_scaling_factor[idx]);
          //auto cache_elem = cache.getCache(0, idx);
          //rho_evolve_super = evolve(rho_evolve_super, cache_elem.rho_inf_eq, cache_elem.scaling_factor);
          //rho_evolve_super = evolve(rho_evolve_super, 
          //    cache_elem.rho_inf_eq, cache_elem.scaling_factor,
          //    rotate_mat_super, rotate_mat_super_inv);
        }
        
        ++rotor_cnt;
        cnt -= mas_inc_cnt;
        const double new_angle = gamma0 + t * mas_frequency * 2.0 * pi;
        mas_angle.gamma(new_angle);
      }
      else {    
        /// if pulse length is less than an integer multiple of rotor step size.
        //auto [rotate_mat_super, rotate_at_super_inv] = 
        //  calcRotationSuperOps(ham_offset, g, inc, cnt);

        const auto ham = packets.genMatrix(temp_euler);
        const auto ham_lab = ham + ham_offset;
        auto [h_super, gamma_super_internal, rho_eq_super] =
          calcSuperOpsForMasterEq(ham, ham_lab, 
              //rotate_mat_super, rotate_mat_super_inv, 
              rpackets, temperature);
        const auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
        const auto scaling_factor = calcExpEvolve(super_op, inc, cnt);
        //rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
        //    scaling_factor, rotate_mat_super, rotate_mat_super_inv);
        rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, scaling_factor);
        cnt = 0;
      }
    }
    return rho_evolve_super;
  }

  std::pair<std::vector<std::pair<double, double>>, MatrixCxDbl>
    evolveMASCnstEmr(
      const MatrixCxDbl &rho_prev_super,
      const MatrixCxDbl &acq_mat_super,
      double t0,
      double result_ref,
      double mas_frequency,
      //const pulseseq::Component &comp,
      const PacketCollection &packets,
      const std::vector<RelaxationPacket> &rpackets,
      const MatrixCxDbl &ham_offset,
      const Euler<> &spin_sys_euler,
      Euler<> mas_angle,
      //const Gyrotron &g,
      double inc, 
      std::uint64_t cnt,    ///< comp size
      std::uint64_t mas_inc_cnt,
      std::uint64_t total_rotor_cnt,  ///< # of steps per rotor period
      double temperature)
  {
    //auto cache = EvolutionCache(total_rotor_cnt, 1);
    vector<MatrixCxDbl> cache_rho;
    vector<MatrixCxDbl> cache_scaling_factor;
    cache_rho.reserve(total_rotor_cnt);
    cache_scaling_factor.reserve(total_rotor_cnt);
    // auto [rotate_mat_super, rotate_mat_super_inv] = 
    //   calcRotationSuperOps(ham_offset, g, inc, mas_inc_cnt);
#ifndef NDEBUG
    cout << "mas_inc_cnt: " << mas_inc_cnt << endl;
    cout << "total_rotor_cnt: " << total_rotor_cnt << endl;
#endif
    
    double t = 0.0;
    uint64_t rotor_cnt = 0u;  ///< number of mas_inc_cnt
    const double mas_inc = inc * static_cast<double>(mas_inc_cnt);

    MatrixCxDbl rho_evolve_super = rho_prev_super;
    size_t gamma_step_size = static_cast<size_t>(std::round(GAMMA_STEP_MIN/mas_inc));
    gamma_step_size += (gamma_step_size == 0);
    MatrixCxDbl gamma_super;
    
#ifndef NDEBUG
    auto print_mas_info = [&](const Euler<> &euler1, const Euler<> &euler2) {
      cout << "[" << t << "] "
           << "\n\t\tmas_angle: " << mas_angle 
           << "\n\t\tspin_sys_euler: " << euler1
           << "\n\t\ttemp_euler: " << euler2
           << "\n\t\t" << "rotor_cnt: " << rotor_cnt << " \t" << "cnt: " << cnt << endl;
    };
#endif

    std::vector<std::pair<double, double>> results;
    const double gamma0 = mas_angle.gamma();
    while(cnt > 0){
      auto temp_euler = mas_angle * spin_sys_euler;
#ifndef NDEBUG
      print_mas_info(spin_sys_euler, temp_euler);
#endif
      if(cnt >= mas_inc_cnt){
        // increment t here, otherwise will have two points on the 1st point zero
        t += static_cast<double>(mas_inc_cnt) * inc;
        if(rotor_cnt < total_rotor_cnt) { // need to save to cache
          const auto ham = packets.genMatrix(temp_euler);
          const auto ham_lab = ham + ham_offset;
          const auto rho_ss = genRhoEq(ham_lab, temperature);
          const auto rho_ss_super = ::dnpsoup::flatten(rho_ss, 'c');
          //if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
          //  rho_ss_super = rotate_mat_super * rho_ss_super;
          //}
          if(rotor_cnt % gamma_step_size == 0) {
            gamma_super = calcGammaSuper(rho_ss, rpackets);
          }
          //gamma_super = calcGammaSuper(rho_ss, rpackets);
          const auto h_super = commutationSuperOp(ham);
          const auto super_op = calcLambdaSuper(h_super, gamma_super);
          const auto rho_eq_super = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
          const auto scaling_factor = calcExpEvolve(super_op, inc, mas_inc_cnt);
          //rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
          //    scaling_factor, rotate_mat_super, rotate_mat_super_inv);
          rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, scaling_factor);
          cache_rho.push_back(rho_eq_super);
          cache_scaling_factor.push_back(scaling_factor);
          //cache.saveCache(comp, 
          //    EvolutionCacheElement(
          //      std::move(scaling_factor), std::move(rho_eq_super))); 
          //rho_evolve_super = rho_eq_super;
          double result = ::dnpsoup::projectionNorm(
              rho_evolve_super, acq_mat_super).real();
          results.push_back(make_pair(t0 + t, result/result_ref));
        }
        else{ // retrieve from cache
          // cache can be big, avoid using it
          const auto idx = rotor_cnt % total_rotor_cnt;
          if(idx >= cache_rho.size()) {
            cout << "overflow cache_row: " << idx << endl;
          }
          if(idx >= cache_scaling_factor.size()) {
            cout << "overflow cache_scaling_factor: " << idx << endl;
          }
          //auto cache_elem = cache.getCache(0, idx);
          //rho_evolve_super = evolve(rho_evolve_super, 
          //    cache_elem.rho_inf_eq, cache_elem.scaling_factor,
          //    rotate_mat_super, rotate_mat_super_inv);
          rho_evolve_super = evolve(rho_evolve_super, 
              cache_rho[idx], cache_scaling_factor[idx]);
          double result = ::dnpsoup::projectionNorm(
              rho_evolve_super, acq_mat_super).real();
          results.push_back(make_pair(t0 + t, result/result_ref));
        }
        
        ++rotor_cnt;
        cnt -= mas_inc_cnt;
        double new_angle = gamma0 + t * mas_frequency * 2.0 * pi;
        mas_angle.gamma(new_angle);
      }
      else {    /// if pulse length is not a integer multiple of rotor period step size.
        t += static_cast<double>(cnt) * inc;
        const auto ham = packets.genMatrix(temp_euler);
        const auto ham_lab = ham + ham_offset;
        auto [h_super, gamma_super_internal, rho_eq_super] =
          calcSuperOpsForMasterEq(ham, ham_lab, 
              //rotate_mat_super, rotate_mat_super_inv, 
              rpackets, temperature);
        const auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
        const auto scaling_factor = calcExpEvolve(super_op, inc, cnt);
        //rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, 
        //    scaling_factor, rotate_mat_super, rotate_mat_super_inv);
        rho_evolve_super = evolve(rho_evolve_super, rho_eq_super, scaling_factor);
        //rho_evolve_super = rho_eq_super;
        cnt = 0;
        double result = ::dnpsoup::projectionNorm(
            rho_evolve_super, acq_mat_super).real();
        results.push_back(make_pair(t0 + t, result/result_ref));
      }
    }
    return make_pair(results, rho_evolve_super);
  }
} // namespace dnpsoup

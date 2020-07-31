#include "configure_dnpsoup.h"
#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/rho_eq.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
#include "dnpsoup_core/spin_physics_components/evolve.h"
#include "dnpsoup_core/spin_physics_components/EvolutionCache.h"
#include "dnpsoup_core/spin_physics_components/EvolutionCacheStatic.h"
#include "dnpsoup_core/spin_physics_components/MasterEqTerms.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/constants.h"
#include "lean/threadpool.h"
#include <complex>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <memory>

using namespace std;
//#define TEST_LSTSQ


namespace dnpsoup {
namespace DnpRunner {
    std::vector<std::vector<double>> calcEigenValues(
        const Magnet &m,
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const Euler<> &spin_sys_euler
        )
    {
      auto irradiated_types = spin_sys.irradiated();
      pulseseq::Component default_comp;
      for(const auto &t : irradiated_types){
        default_comp.insert_or_assign(t, pulseseq::EMRadiation());
      }

      double inc = seq.getIncrement();
      inc = roundToCycle(inc, g.em_frequency);
      vector<vector<double>> results;
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();

      double t = 0.0;
      Euler<> mas_angle = p.magic_angle;
      double t_inc = p.mas_increment;
      if(p.mas_frequency > 1.0){
        double t_inc_limit = 1.0/p.mas_frequency * 0.02; ///< override mas_inc
        if(t_inc_limit < t_inc){
          t_inc = t_inc_limit;
          cout << "Overriding existing mac_increment with " << t_inc_limit << endl;
        }
      }
      if(t_inc < eps && t_inc < inc){
        t_inc = inc;
      }
      if(t_inc < eps){
        throw PulseSequenceError("Incrementation too small.");
      }
      const uint64_t n_inc = static_cast<uint64_t>(std::round(t_inc/inc));
//#ifndef NDEBUG
//      cout << "n_inc = t_inc/inc = " << n_inc << endl;
//#endif

      packets.setPropertyValue(ValueName::b0, m.b0);
      offset_packets.setPropertyValue(ValueName::b0, m.b0);
      auto spin_ids = spin_sys.getSpinIds(SpinType::e);
      for(const auto &sid : spin_ids){
        packets.setPropertyValue(
            InteractionType::Shielding, sid, ValueName::offset, g.em_frequency);
        offset_packets.setPropertyValue(
            InteractionType::Offset, sid, ValueName::offset, g.em_frequency);
      }

      uint64_t comp_size;
      uint64_t idx = 0;

      do {
        pulseseq::Component comp;
        std::tie(comp, comp_size, idx) = seq.next();
        if(idx == seq.size()) break;
        packets.updatePulseSeqComponent(default_comp);
        packets.updatePulseSeqComponent(comp);

        while(comp_size > 0){
          uint64_t step_sz = std::min(comp_size, n_inc);

          mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);

          const auto temp_angle = mas_angle * spin_sys_euler;
          MatrixCxDbl hamiltonian = packets.genMatrix(temp_angle);
          MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(temp_angle);
          //auto eigen_values = ::dnpsoup::eigenVal(hamiltonian);
          auto eigen_values = ::dnpsoup::eigenVal(hamiltonian + hamiltonian_offset);
          vector<double> eigen_values_temp;
          eigen_values_temp.push_back(t);
          for(size_t i = 0; i < eigen_values.nelements(); ++i){
            eigen_values_temp.push_back(eigen_values(i, 0));
          }
          results.push_back(eigen_values_temp);

          comp_size -= step_sz;
          t += static_cast<double>(step_sz) * inc;
        }
      } while(idx < seq.size());

      return results;
    }

    double calcIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const Euler<> &sample_euler,
        bool ignore_all_power)
    {
      constexpr double eps = std::numeric_limits<double>::epsilon();
      // Active rotation
      const Euler<> spin_sys_euler = sample_euler * spin_sys.getEuler();
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();
      auto rpackets = spin_sys.summarizeRelaxation();
      auto acq_mat = spin_sys.acquireOn(acq_spin);
      auto acq_mat_super = ::dnpsoup::flatten(acq_mat, 'c');
      double mas_inc = p.mas_increment;

      double rotor_period = 0.0;
      uint64_t total_rotor_cnt = 0u;
      const bool has_mas = (p.mas_frequency > eps);
      if(has_mas) {
        rotor_period = 1.0 / p.mas_frequency;
        if(p.mas_increment > eps){
          total_rotor_cnt = static_cast<uint64_t>(
              round(rotor_period/p.mas_increment));
          // total_rotor_cnt >= 1
          total_rotor_cnt += (total_rotor_cnt == 0u);
          mas_inc = rotor_period / static_cast<double>(total_rotor_cnt);
        }
        else {
          total_rotor_cnt = 100;
          mas_inc = rotor_period / static_cast<double>(total_rotor_cnt);
          ///< 1% of MAS
        }
      }
      double t = 0.0;
      auto mas_angle = p.magic_angle;

      auto irradiated_types = spin_sys.irradiated();
      pulseseq::Component default_comp;
      for(const auto &t : irradiated_types){
        default_comp.insert_or_assign(t, pulseseq::EMRadiation());
      }
      double inc = seq.getIncrement();
      //inc = roundToCycle(inc, g.em_frequency);

      packets.setPropertyValue(ValueName::b0, m.b0);
      offset_packets.setPropertyValue(ValueName::b0, m.b0);
      auto e_spin_ids = spin_sys.getSpinIds(SpinType::e);
      for(const auto &sid : e_spin_ids){
        offset_packets.setPropertyValue(
            InteractionType::Offset, sid, ValueName::offset, g.em_frequency);
        packets.setPropertyValue(
            InteractionType::Shielding, sid, ValueName::offset, g.em_frequency);
      }

      unique_ptr<EvolutionCacheStatic> uptr_cache = nullptr;
      if (!has_mas) {
#ifndef NDEBUG
        cout << "EvolutionCacheStatic capacity: " << seq.uniqueEmrsCount() << endl;
#endif
        uptr_cache = make_unique<EvolutionCacheStatic>(
            seq.uniqueEmrsCount());
      }
//#ifndef NDEBUG
//      std::cout << "mas_inc_cnt: " << mas_inc_cnt << "\n";
//#endif

      const auto angle = mas_angle * spin_sys_euler;
      MatrixCxDbl hamiltonian = packets.genMatrix(angle);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(angle);
      const auto hamiltonian_lab = hamiltonian + hamiltonian_offset;
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab, p.temperature);
      auto rho0_evolve = rho0_lab;
      auto rho0_evolve_super = ::dnpsoup::flatten(rho0_evolve, 'c');
      if(seq.size() == 0){
        double result = 
          ::dnpsoup::projectionNorm(rho0_evolve_super, acq_mat_super).real();
        return result;
      }

      // --------------------------------------------------------
      /// via MasterEqTerms static case (MAS=0)
      // --------------------------------------------------------
      if(!has_mas) {
        auto [terms, ptr_packets] = genMasterEqTerms(&packets, rpackets, 
            hamiltonian_offset, seq, irradiated_types, g, angle, p.temperature, inc);
        rho0_evolve_super = evolve(rho0_evolve_super, terms);
      } else {  // MAS
        std::uint64_t comp_size = 0u;
        uint64_t idx = 0;
        pulseseq::Component comp;
//#ifndef NDEBUG
//        cout << "seq.size() = " << seq.size() << endl;
//#endif
        if(ignore_all_power) {
          std::uint64_t total_seq_size = 0u;
          while(idx < seq.size() || comp_size > 0) {
            std::tie(comp, comp_size, idx) = seq.next();
            total_seq_size += comp_size;
            if(idx >= seq.size()) break;
          }
          packets.updatePulseSeqComponent(default_comp);
          MasterEqTerms terms = genMasterEqTermsMAS(
              &packets, rpackets, hamiltonian_offset,// g,
              spin_sys_euler, mas_angle, total_seq_size,
              p.temperature, p.mas_frequency, inc, mas_inc); 
          rho0_evolve_super = evolve(rho0_evolve_super, terms);
        } else {
          while(idx < seq.size() || comp_size > 0) {
            /// step-wise consume pulse sequence
            std::tie(comp, comp_size, idx) = seq.next();
            if(idx >= seq.size()) break;
            // reset
            packets.updatePulseSeqComponent(default_comp);
            // update with new value
            packets.updatePulseSeqComponent(comp);

            MasterEqTerms terms = genMasterEqTermsMAS(
                &packets, rpackets, hamiltonian_offset,// g,
                spin_sys_euler, mas_angle, comp_size,
                p.temperature, p.mas_frequency, inc, mas_inc); 
            rho0_evolve_super = evolve(rho0_evolve_super, terms);

            //uint64_t mas_inc_cnt = static_cast<uint64_t>(round(mas_inc/inc));
            // EvolutionCache is used per comp
            //rho0_evolve_super = evolveMASCnstEmr(
            //    rho0_evolve_super, 
            //    p.mas_frequency, //comp, 
            //    packets, rpackets, hamiltonian_offset,
            //    spin_sys_euler, mas_angle,// g,
            //    inc, comp_size, mas_inc_cnt, total_rotor_cnt, p.temperature);
            t += inc * static_cast<double>(comp_size);
            mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
          } // while
        } // if not ignore_all_power
      } // if has_mas
      
      double result = ::dnpsoup::projectionNorm(
          rho0_evolve_super, acq_mat_super).real();
      return result;
    }

    std::vector<std::pair<double, double>> calcPowderBuildUpEnhancement(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores)
    {
      auto results = calcPowderBuildUp(
          m, g, p, spin_sys, seq, 
          acq_spin, spin_sys_eulers, ncores);
      const bool has_mas = std::abs(p.mas_frequency) > eps;
      if(has_mas) {
        const auto ref_intensities = calcPowderBuildUp(
            m, g, p, spin_sys, seq, acq_spin, spin_sys_eulers, ncores, true);
      }
      double ref_intensity = calcPowderIntensity(
            m, g, p, spin_sys, PulseSequence(), 
            acq_spin, spin_sys_eulers, ncores);
      /// to avoid divide by zero error
      ref_intensity += std::abs(ref_intensity) < eps;
#ifndef NDEBUG
      cout << "ref intensity: " << ref_intensity << endl;
#endif
      for(size_t i = 0; i < results.size(); ++i){
        results[i].second /= ref_intensity;
      }
      return results;
    }

    std::vector<std::pair<double, double>> calcPowderBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &eulers,
        int ncores,
        bool ignore_all_power)
    {
      std::vector<std::pair<double, double>> results;
      const double scaling_factor = 1.0 / static_cast<double>(eulers.size());
#ifndef NDEBUG
      cout << "powder scaling factor: " << scaling_factor << endl;
#endif
      if(ncores == 1) {
        for(const auto &euler : eulers){
          auto xtal_results = calcBuildUp(m, g, p, spin_sys, seq,
              acq_spin, euler, ignore_all_power);
          if(results.size() == 0){  // initial crystal point
            for(const auto &pt_res : xtal_results){
              const double temp = pt_res.second * std::sin(euler.beta());
              results.push_back(make_pair(pt_res.first, temp));
            }
          }
          else {  // additional points
#ifndef NDEBUG
            if(results.size() != xtal_results.size()){
              throw SizeMismatchError("buildup number of points mismatch."); 
            }
#endif
            for(size_t i = 0; i < results.size(); ++i){
              const double temp = xtal_results[i].second * std::sin(euler.beta());
              results[i].second += temp;
            }
          }
        } // for loop
      }
      else {  // multithreading
        ::lean::ThreadPool<std::vector<std::pair<double, double>>> tpool(ncores);
        for(const auto &euler : eulers){
          auto task = [=](){
            auto xtal_results = 
              calcBuildUp(m, g, p, spin_sys, seq, acq_spin, euler, ignore_all_power);
            const double factor = std::sin(euler.beta());
            for(auto &xtal_result : xtal_results){
              xtal_result.second *= factor;
            }

            return xtal_results;
          };
          tpool.add_task(std::move(task));
        }
        tpool.run();
        auto intensities = tpool.get_results();
        for(auto &xtal_results : intensities){
          if(results.size() == 0) {
            for(const auto &pt_res : xtal_results) {
              results.push_back(pt_res);
            }
          }
          else {
#ifndef NDEBUG
            if(results.size() != xtal_results.size()){
              throw SizeMismatchError("buildup number of points mismatch."); 
            }
#endif
            for(size_t i = 0; i < results.size(); ++i){
              results[i].second += xtal_results[i].second;
            }
          }
        }
      }
      for(auto &pt : results){
        pt.second *= (scaling_factor/pi * 4.0);
      }
      std::cout << std::endl;
      return results;
    }

    std::vector<pair<double, double>> calcBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const Euler<> &sample_euler,
        bool ignore_all_power)
    {
      // Active Rotation
      const Euler<> spin_sys_euler = sample_euler * spin_sys.getEuler();

      constexpr double eps = std::numeric_limits<double>::epsilon();
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();
      auto rpackets = spin_sys.summarizeRelaxation();
      auto acq_mat = spin_sys.acquireOn(acq_spin);
      auto acq_mat_super = ::dnpsoup::flatten(acq_mat, 'c');
      double mas_inc = p.mas_increment;

      double rotor_period = 0.0;
      uint64_t total_rotor_cnt = 0u;
      const bool has_mas = (p.mas_frequency > eps);
      if(has_mas) {
        rotor_period = 1.0 / p.mas_frequency;
        if(p.mas_increment > eps){
          total_rotor_cnt = static_cast<uint64_t>(
              round(rotor_period/p.mas_increment));
          // total_rotor_cnt >= 1
          total_rotor_cnt += (total_rotor_cnt == 0u);
          mas_inc = rotor_period / static_cast<double>(total_rotor_cnt);
        }
        else {
          total_rotor_cnt = 100;
          mas_inc = rotor_period / static_cast<double>(total_rotor_cnt);
          ///< 1% of MAS
        }
      }
      else {
        total_rotor_cnt = 1u;
      }
      double t = 0.0;
      auto mas_angle = p.magic_angle;

      auto irradiated_types = spin_sys.irradiated();
      pulseseq::Component default_comp;
      for(const auto &t : irradiated_types){
        default_comp.insert_or_assign(t, pulseseq::EMRadiation());
      }
      double inc = seq.getIncrement();
      //inc = roundToCycle(inc, g.em_frequency);
#ifndef NDEBUG
      if(mas_inc < inc - eps) {
        throw InputError("For BuildUp, MAS increment is needed. (>= inc of a pulse sequence)");
      }
#endif
      uint64_t mas_inc_cnt = static_cast<uint64_t>(round(mas_inc/inc));

      packets.setPropertyValue(ValueName::b0, m.b0);
      offset_packets.setPropertyValue(ValueName::b0, m.b0);
      auto e_spin_ids = spin_sys.getSpinIds(SpinType::e);
      for(const auto &sid : e_spin_ids){
        offset_packets.setPropertyValue(
            InteractionType::Offset, sid, ValueName::offset, g.em_frequency);
        packets.setPropertyValue(
            InteractionType::Shielding, sid, ValueName::offset, g.em_frequency);
      }

      unique_ptr<EvolutionCacheStatic> uptr_cache = nullptr;
      if (!has_mas) {
#ifndef NDEBUG
        cout << "EvolutionCacheStatic capacity: " << seq.uniqueEmrsCount() << endl;
#endif
        uptr_cache = make_unique<EvolutionCacheStatic>(
            seq.uniqueEmrsCount());
      }

      const auto temp_angle = mas_angle * spin_sys_euler;
      MatrixCxDbl hamiltonian = packets.genMatrix(temp_angle);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(temp_angle);
      auto hamiltonian_lab = hamiltonian + hamiltonian_offset;
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab, p.temperature);
      const auto rho_ss_super = ::dnpsoup::flatten(rho0_lab, 'c');
      auto rho0_evolve = rho0_lab;
      auto rho0_evolve_super = ::dnpsoup::flatten(rho0_evolve, 'c');
      //MatrixCxDbl gamma_super = calcGammaSuper(rho0_lab, rpackets);
      //const auto h_super = commutationSuperOp(hamiltonian);
      //const auto super_op = calcLambdaSuper(h_super, gamma_super);
      //const auto rho_eq_super = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
      //MatrixCxDbl rho0_evolve_super = rho_eq_super;

      double val = ::dnpsoup::projectionNorm(rho0_evolve_super, acq_mat_super).real();

      constexpr double result_ref = 1.0;  ///< intensity, not enhancement
      vector<pair<double, double>> results;
      results.push_back(make_pair(0.0, val));

      std::uint64_t comp_size = 0u;
      uint64_t idx = 0;
      pulseseq::Component comp;

      while(idx < seq.size() || comp_size > 0) {
        /// step-wise consume pulse sequence
        std::tie(comp, comp_size, idx) = seq.next();

        if(idx >= seq.size()) break;

        packets.updatePulseSeqComponent(default_comp);
        if(!ignore_all_power) {
          packets.updatePulseSeqComponent(comp);
        } else {
          comp = default_comp;
        }
        if(has_mas) {   ///< mas
#ifdef DNPSOUP_VERBOSE
          cout << "call evolveMASCnstEmr()..." << "\n";
#endif
          ///< with MAS or if need to recalculate super operators
          /// result_ref = 1.0 for intensity (not enhancement)
          vector<pair<double, double>> temp_results;
          std::tie(temp_results, rho0_evolve_super) = evolveMASCnstEmr(
                rho0_evolve_super, acq_mat_super, t, result_ref,
                p.mas_frequency, //comp, 
                packets, rpackets, hamiltonian_offset,
                spin_sys_euler, mas_angle,// g, 
                inc, comp_size, 
                mas_inc_cnt, total_rotor_cnt, p.temperature);
          t += static_cast<double>(comp_size) * inc;
          mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
          std::copy(temp_results.begin(), temp_results.end(), 
              std::back_inserter(results));
        }
        else {    ///< static
#ifdef DNPSOUP_VERBOSE
          cout << "static evolve..." << "\n";
#endif
          auto [rotate_mat_super, rotate_mat_super_inv] = calcRotationSuperOps(
              hamiltonian_offset, g, inc, mas_inc_cnt);
          const auto temp_euler = mas_angle * spin_sys_euler;
          auto [cache_idx, has_super_op_section] 
            = uptr_cache->getCacheIdentity(comp, mas_inc_cnt);
          MatrixCxDbl h_super, gamma_super_internal, rho_eq_super;
          MatrixCxDbl scaling_factor;
          const auto ham = packets.genMatrix(temp_euler);
          const auto ham_lab = ham + hamiltonian_offset; 
          if(!has_super_op_section){
            std::tie(h_super, gamma_super_internal, rho_eq_super) = 
              calcSuperOpsForMasterEq(ham, ham_lab,
                  rotate_mat_super, rotate_mat_super_inv,
                  rpackets, p.temperature);
            auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
            scaling_factor = calcExpEvolve(super_op, inc, mas_inc_cnt);
            uptr_cache->saveCache(
                comp, 
                super_op, 
                rho_eq_super, 
                scaling_factor, 
                mas_inc_cnt);
          }
          else {
            std::tie(scaling_factor, rho_eq_super) = uptr_cache->getCache(
                cache_idx, mas_inc_cnt);
          }
#ifdef DNPSOUP_VERBOSE
          //cout << setprecision(12);
          //cout << "emr comp: " << comp << endl;
          //cout << "\nrho0_evolve_super:\n" << rho0_evolve_super
          //     << "\nrho_eq_super:\n" << rho_eq_super
          //     << "\nscaling_factor:\n" << scaling_factor
          //     << "\nrotate_mat_super:\n" << rotate_mat_super
          //     << "\nrotate_mat_super_inv:\n" << rotate_mat_super_inv << endl;
#endif
          while(comp_size >= mas_inc_cnt){
#ifdef DNPSOUP_VERBOSE
            cout << "in static case, call evolve() with comp_size: " << comp_size << endl;
#endif

#ifndef NDEBUG
#endif
            rho0_evolve_super = evolve(rho0_evolve_super, rho_eq_super,
                scaling_factor,
                rotate_mat_super, rotate_mat_super_inv);
            t += inc * static_cast<double>(mas_inc_cnt);
            double val = ::dnpsoup::projectionNorm(
                rho0_evolve_super, acq_mat_super).real();
            results.push_back(make_pair(t, val));
            comp_size -= mas_inc_cnt;
          } // while
          if(comp_size > 0){
            std::tie(rotate_mat_super, rotate_mat_super_inv) = calcRotationSuperOps(
                hamiltonian_offset, g, inc, comp_size);
            std::tie(h_super, gamma_super_internal, rho_eq_super) = 
              calcSuperOpsForMasterEq(ham, ham_lab,
                  rotate_mat_super, rotate_mat_super_inv,
                  rpackets, p.temperature);
            auto super_op = calcLambdaSuper(h_super, gamma_super_internal);
            scaling_factor = calcExpEvolve(
                super_op, inc, comp_size);
            rho0_evolve_super = evolve(rho0_evolve_super, rho_eq_super, 
                scaling_factor,
                rotate_mat_super, rotate_mat_super_inv);
            t += inc * static_cast<double>(comp_size);
            double val = ::dnpsoup::projectionNorm(
                rho0_evolve_super, acq_mat_super).real();
            results.push_back(make_pair(t, val));
          }
        }
      } // while
      std::cout << "." << std::flush;

      return results;
    }

    std::vector<std::pair<double, double>> 
      calcBuildUpEnhancement(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const Euler<> &sample_euler)
    {
      auto intensities = calcBuildUp(
          m, g, p, spin_sys, seq, acq_spin, sample_euler);
      const bool has_mas = std::abs(p.mas_frequency) > eps;
      if(has_mas) {
        const auto ref_intensities = calcBuildUp(
          m, g, p, spin_sys, seq, acq_spin, sample_euler, true);
        for(size_t i = 0; i < intensities.size(); ++i) {
          const double val = ref_intensities[i].second + (std::abs(ref_intensities[i].second) < eps);
          intensities[i].second /= val;
        }
      } else {
        double ref_intensity = calcIntensity(
            m, g, p, spin_sys, PulseSequence(), acq_spin, sample_euler);
        //ref_intensity += (ref_intensity < eps);
        for(auto &value_pair : intensities){
          value_pair.second /= ref_intensity;
        }
      }
      return intensities;
    }


    std::vector<std::pair<double, double>> calcFieldProfile(
        const std::vector<Magnet> &fields, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &eulers,
        [[maybe_unused]] int ncores)
    {
      std::vector<std::pair<double, double>> result;
      const bool has_mas = std::abs(p.mas_frequency) > eps;
      if(eulers.size() == 1){
        const auto &euler = eulers[0];
        ::lean::ThreadPool<pair<double, double>> tpool(ncores);
        for(const auto &field : fields){
          auto task = [=](){
            double ref = 0.0;
            if(has_mas) {
              ref = 
                calcIntensity(field, g, p, spin_sys, seq, acq_spin, euler, true);
            } else {
              ref = 
                calcIntensity(field, g, p, spin_sys, PulseSequence(), acq_spin, euler);
            }
            ref += (std::abs(ref) < eps);
            auto intensity = 
              calcIntensity(field, g, p, spin_sys, seq, acq_spin, euler);
            std::cout << "." << std::flush;
            return make_pair(field.b0, intensity/ref);
          };
          tpool.add_task(std::move(task));
        }
        tpool.run();
        result = tpool.get_results();
        std::sort(result.begin(), result.end(),
            [](const pair<double, double> &val1, 
              const pair<double, double> &val2){
              return val1.first < val2.first;
        });
      }
      else {
        for(const auto &field : fields){
          //constexpr double ref = 1.0;
          double ref = 0.0;
          if(has_mas) {
            ref = calcPowderIntensity(
                field, g, p, spin_sys, seq,
                acq_spin, eulers, ncores, true);
          } else {
            ref = calcPowderIntensity(
                field, g, p, spin_sys, PulseSequence(),
                acq_spin, eulers, ncores);
          }
          ref += std::abs(ref) < eps;
          const double res = calcPowderIntensity(
              field, g, p, spin_sys, seq, acq_spin, eulers, ncores);
          const double ratio = res/ref;
          result.push_back(std::make_pair(field.b0, ratio));
          std::cout << "." << std::flush;
        }
      }
      std::cout << std::endl;
      return result;
    }

    std::vector<std::pair<double, double>> calcFieldProfile(
        const Magnet &m, 
        const std::vector<Gyrotron> &emrs,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &eulers,
        [[maybe_unused]] int ncores)
    {
      std::vector<std::pair<double, double>> result;
      const bool has_mas = std::abs(p.mas_frequency) > eps;
      if(eulers.size() == 1){
        const auto &euler = eulers[0];
        ::lean::ThreadPool<pair<double, double>> tpool(ncores);
        for(const auto &emr : emrs){
          auto task = [=](){
            double ref = 0.0;
            if(has_mas) {
              ref = calcIntensity(m, emr, p, spin_sys, seq, acq_spin, euler, true);
            } else {
              ref = calcIntensity(m, emr, p, spin_sys, PulseSequence(), acq_spin, euler);
            }
            ref += (std::abs(ref) < eps);
            auto intensity = 
              calcIntensity(m, emr, p, spin_sys, seq, acq_spin, euler);
            return make_pair(emr.em_frequency, intensity/ref);
          };
          tpool.add_task(std::move(task));
        }
        tpool.run();
        result = tpool.get_results();
        std::sort(result.begin(), result.end(),
            [](const pair<double, double> &val1, 
              const pair<double, double> &val2){
              return val1.first < val2.first;
        });
      }
      else {
        for(const auto &g : emrs) {
          //constexpr double ref = 1.0;
          double ref = 0.0;
          if (has_mas) {
            ref = calcPowderIntensity(
                m, g, p, spin_sys, seq,
                acq_spin, eulers, ncores, true);
          } else {
            ref = calcPowderIntensity(
                m, g, p, spin_sys, PulseSequence(),
                acq_spin, eulers, ncores);
          }
          ref += (std::abs(ref) < eps);
          const double res = calcPowderIntensity(
              m, g, p, spin_sys, seq, acq_spin, eulers, ncores);
          const double ratio = res/ref;
          result.emplace_back(make_pair(g.em_frequency, ratio));
          std::cout << "." << std::flush;
        }
      }
      std::cout << std::endl;
      return result;
    }

    double calcPowderIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        PulseSequence seq,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &eulers,
        [[maybe_unused]] int ncores,
        bool ignore_all_power)
    {
      double result = 0.0;
      const double scaling_factor = 1.0 / static_cast<double>(eulers.size());
      if(ncores == 1) {
        for(const auto &e : eulers){
          auto xtal_intensity = 
            calcIntensity(m, g, p, spin_sys, seq, acq_spin, e, ignore_all_power);
          result += xtal_intensity * std::sin(e.beta());
        }
      }
      else {
        ::lean::ThreadPool<double> tpool(ncores);
        for(const auto &e : eulers){
          auto task = [=](){
            auto xtal_intensity = 
              calcIntensity(m, g, p, spin_sys, seq, acq_spin, e, ignore_all_power);
            return xtal_intensity * std::sin(e.beta());
          };
          tpool.add_task(std::move(task));
        }
        tpool.run();
        auto intensities = tpool.get_results();
        for(auto &intensity : intensities){
          result += intensity;
        }
      }
      result = result * scaling_factor / pi * 4.0;
      //std::cout << "." << std::flush;

      return result;
    }

    MatrixCxDbl propagate(
        const MatrixCxDbl &rho_prev_super, 
        const MatrixCxDbl &hamiltonian,
        const MatrixCxDbl &hamiltonian_lab,
        const MatrixCxDbl &rotate_mat_super,
        const MatrixCxDbl &rotate_mat_super_inv,
        const std::vector<RelaxationPacket> &rpackets,
        double dt, std::uint64_t cnt, double temperature)
    {
      auto [h_super, gamma_super_int, rho_eq_super] = calcSuperOpsForMasterEq(
          hamiltonian, hamiltonian_lab, rotate_mat_super, rotate_mat_super_inv,
          rpackets, temperature);
      auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;

      //auto rho_prev_super = ::dnpsoup::flatten(rho_prev, 'c');
      auto rho_super = rho_prev_super;
      if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
        rho_super = rotate_mat_super * rho_prev_super;
      }

      auto scaling_factor = calcExpEvolve(super_op, dt, cnt);
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

    MatrixCxDbl propagate(
        const MatrixCxDbl &rho_prev_super, 
        const PacketCollection &packets,
        const MatrixCxDbl &hamiltonian_offset,
        const std::vector<RelaxationPacket> &rpackets,
        const Gyrotron &g,
        const Euler<> &euler,
        double dt,
        std::uint64_t cnt,
        double temperature
        )
    {
      auto [rotate_mat_super, rotate_mat_super_inv] = calcRotationSuperOps(
          hamiltonian_offset, g, dt, cnt);

      // auto rho_prev_super = ::dnpsoup::flatten(rho_prev, 'c');
      auto rho_super = rho_prev_super;
      if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
        rho_super = rotate_mat_super * rho_prev_super;
      }
      auto ham = packets.genMatrix(euler);
      auto ham_lab = ham + hamiltonian_offset;

      auto [h_super, gamma_super_int, rho_eq_super] = calcSuperOpsForMasterEq(
          ham, ham_lab, rotate_mat_super, rotate_mat_super_inv,
          rpackets, temperature);

      // iH + Gamma
      auto super_op = calcLambdaSuper(h_super, gamma_super_int);

      auto scaling_factor = calcExpEvolve(super_op, dt, cnt);

      rho_super = evolve(rho_super, rho_eq_super, scaling_factor,
          rotate_mat_super, rotate_mat_super_inv);
      return rho_super;
    }

} // namespace DnpRunner
} // namespace dnpsoup

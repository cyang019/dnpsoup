#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/rho_eq.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
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
#include <string>
#include <cmath>

using namespace std;
//#define TEST_LSTSQ


namespace dnpsoup {
namespace DnpRunner {
    std::vector<std::vector<double>> calcEigenValues(
        const Magnet &m,
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const Euler<> &spin_sys_euler
        )
    {
      PulseSequence seq;
      istringstream iss(pulse_seq_str);
      try{
        iss >> seq;
      }
      catch(const exception &e){
        throw PulseSequenceError(e.what());
      }

      const double inc = seq.getIncrement();
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
#ifndef NDEBUG
      cout << "n_inc = t_inc/inc = " << n_inc << endl;
#endif

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
        packets.updatePulseSeqComponent(comp);

        while(comp_size > 0){
          uint64_t step_sz = std::min(comp_size, n_inc);

          mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
// #ifndef NDEBUG
//           cout << "[mas] t: " << t << " \t" << "mas_angle:"
//                << " a: " << mas_angle.alpha() 
//                << " b: " << mas_angle.beta()
//                << " g: " << mas_angle.gamma() << endl;
// #endif

          MatrixCxDbl hamiltonian = packets.genMatrix(spin_sys_euler * mas_angle);
          MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler * mas_angle);
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
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler)
    {
      constexpr double eps = std::numeric_limits<double>::epsilon();
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();
      auto rpackets = spin_sys.summarizeRelaxation();
      auto acq_mat = spin_sys.acquireOn(acq_spin);
      double mas_inc = -1;
      if(p.mas_frequency > eps){
        if(p.mas_increment > eps){
          mas_inc = p.mas_increment;
        }
        else {
          mas_inc = 1.0/p.mas_frequency * 0.001;   /// 0.1% of MAS
        }
      }
      double t = 0.0;
      auto mas_angle = p.magic_angle;

      packets.setPropertyValue(ValueName::b0, m.b0);
      offset_packets.setPropertyValue(ValueName::b0, m.b0);
      auto e_spin_ids = spin_sys.getSpinIds(SpinType::e);
      for(const auto &sid : e_spin_ids){
        offset_packets.setPropertyValue(
            InteractionType::Offset, sid, ValueName::offset, g.em_frequency);
        packets.setPropertyValue(
            InteractionType::Shielding, sid, ValueName::offset, g.em_frequency);
      }

      PulseSequence seq;
      istringstream iss(pulse_seq_str);
      try{
        iss >> seq;
      }
      catch(const exception &e){
        throw PulseSequenceError(e.what());
      }
      const double inc = seq.getIncrement();
      uint64_t mas_inc_cnt = 0;
      if (mas_inc > 0){
        mas_inc_cnt = static_cast<uint64_t>(round(mas_inc/inc));
      }
//#ifndef NDEBUG
//      std::cout << "mas_inc_cnt: " << mas_inc_cnt << "\n";
//#endif
      uint64_t idx = 0;

      MatrixCxDbl hamiltonian = packets.genMatrix(spin_sys_euler * mas_angle);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler * mas_angle);
      auto hamiltonian_lab = hamiltonian + hamiltonian_offset;
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab, p.temperature);
      auto rho0_evolve = rho0_lab;
      if(seq.size() == 0){
        double result = ::dnpsoup::projectionNorm(rho0_evolve, acq_mat).real();
        return result;
      }
      std::uint64_t cnt = 0u;    /// keep track of identical hamiltonians
      std::uint64_t comp_size = 0u;
      do {
        pulseseq::Component comp;
        comp_size = 0u;
        if(idx < seq.size()){
          std::tie(comp, comp_size, idx) = seq.next();
        }
//#ifndef NDEBUG
//        std::cout << "comp_size: " << comp_size << " \tidx: " << idx << std::endl;
//#endif

        bool same_comp = packets.hasPulseSeqComponent(comp);
        if(idx >= seq.size()){    // end of sequence
          same_comp = false;
        }

        if(same_comp){
          cnt += comp_size;
          continue;
        }
        else {
          if(mas_inc_cnt > 0){
            while(cnt > 0){
              auto temp_euler = spin_sys_euler * mas_angle;
              if(cnt >= mas_inc_cnt){
                rho0_evolve = evolve(rho0_evolve,
                    packets, hamiltonian_offset, rpackets,
                    g, temp_euler,
                    inc, mas_inc_cnt, p.temperature); 
                t += static_cast<double>(mas_inc_cnt) * inc;
                cnt -= mas_inc_cnt;
              }
              else {
                rho0_evolve = evolve(rho0_evolve,
                    packets, hamiltonian_offset, rpackets,
                    g, temp_euler,
                    inc, cnt, p.temperature); 
                t += static_cast<double>(cnt) * inc;
                cnt = 0;
              }
              mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
            }
          }
          else {    // no MAS
            if(cnt > 0){
              auto temp_euler = spin_sys_euler * mas_angle;
              rho0_evolve = evolve(rho0_evolve,
                  packets, hamiltonian_offset, rpackets,
                  g, temp_euler,
                  inc, cnt, p.temperature); 
              t += static_cast<double>(cnt) * inc;
              cnt = 0;
            }
          }
          packets.updatePulseSeqComponent(comp);
          cnt = comp_size;
        }
      } while(idx < seq.size() || cnt > 0);

      double result = ::dnpsoup::projectionNorm(rho0_evolve, acq_mat).real();
      return result;
    }

    std::vector<std::pair<double, double>> calcPowderBuildUpEnhancement(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores)
    {
      auto raw_results = calcPowderBuildUp(
          m, g, p, spin_sys, pulse_seq_str, acq_spin, spin_sys_eulers, ncores);
      const double ref_intensity = calcPowderIntensity(
            m.b0, g, p, spin_sys, "{}", acq_spin, spin_sys_eulers, ncores);
      auto res = raw_results;
      for(size_t i = 0; i < raw_results.size(); ++i){
        res[i].second /= ref_intensity;
      }
      return res;
    }

    std::vector<std::pair<double, double>> calcPowderBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        int ncores)
    {
      std::vector<std::pair<double, double>> results;
      const double scaling_factor = 1.0 / static_cast<double>(spin_sys_eulers.size());
      if(ncores == 1) {
        for(const auto &e : spin_sys_eulers){
          auto xtal_results = calcBuildUp(m, g, p, spin_sys, pulse_seq_str,
              acq_spin, e, false);
          if(results.size() == 0){  // initial crystal point
            for(const auto &pt_res : xtal_results){
              const double temp = pt_res.second * std::sin(e.beta());
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
              const double temp = xtal_results[i].second * std::sin(e.beta());
              results[i].second += temp;
            }
          }
        }
      }
      else {  // multithreading
        ::lean::ThreadPool<std::vector<std::pair<double, double>>> tpool(ncores);
        for(const auto &e : spin_sys_eulers){
          auto task = [=](){
            auto xtal_results = 
              calcBuildUp(m, g, p, spin_sys, pulse_seq_str, acq_spin, e, false);
            const double factor = std::sin(e.beta());
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
      return results;
    }

    std::vector<pair<double, double>> calcBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler,
        bool enhancement)
    {
      constexpr double eps = std::numeric_limits<double>::epsilon();
      PulseSequence seq;
      istringstream iss(pulse_seq_str);
      try{
        iss >> seq;
      }
      catch(const exception &e){
        throw PulseSequenceError(e.what());
      }
      const double inc = seq.getIncrement();

      vector<pair<double, double>> results;
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();
      auto rpackets = spin_sys.summarizeRelaxation();
      auto acq_mat = spin_sys.acquireOn(acq_spin);
      double mas_inc = p.mas_increment;
      if(mas_inc < inc - eps) {
        throw InputError("For BuildUp, MAS increment is needed. (>= inc of a pulse sequence)");
      }

      double t = 0.0;
      auto mas_angle = p.magic_angle;

      packets.setPropertyValue(ValueName::b0, m.b0);
      offset_packets.setPropertyValue(ValueName::b0, m.b0);
      auto e_spin_ids = spin_sys.getSpinIds(SpinType::e);
      for(const auto &sid : e_spin_ids){
        offset_packets.setPropertyValue(
            InteractionType::Offset, sid, ValueName::offset, g.em_frequency);
        packets.setPropertyValue(
            InteractionType::Shielding, sid, ValueName::offset, g.em_frequency);
      }

      uint64_t mas_inc_cnt = 0;
      if (mas_inc > 0){
        mas_inc_cnt = static_cast<uint64_t>(round(mas_inc/inc));
      }
      uint64_t idx = 0;

      MatrixCxDbl hamiltonian = packets.genMatrix(spin_sys_euler * mas_angle);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler * mas_angle);
      auto hamiltonian_lab = hamiltonian + hamiltonian_offset;
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab, p.temperature);
      const double result_ref = 
        enhancement ? ::dnpsoup::projectionNorm(rho0_lab, acq_mat).real() : 1.0;
#ifndef NDEBUG
      //std::cout << "\nresult_ref: " << result_ref << std::endl;
      std::cout << "." << std::flush;
#endif
      auto rho0_evolve = rho0_lab;
      std::uint64_t cnt = 0u;    /// keep track of identical hamiltonians
      std::uint64_t comp_size = 0u;
      do {
        pulseseq::Component comp;
        comp_size = 0u;
        if(idx < seq.size()){
          std::tie(comp, comp_size, idx) = seq.next();
        }

        bool same_comp = packets.hasPulseSeqComponent(comp);
        if(idx >= seq.size()){    // end of sequence
          same_comp = false;
        }

        if(same_comp){
          cnt += comp_size;
          continue;
        }
        else {
          while(cnt > 0){
            auto temp_euler = spin_sys_euler * mas_angle;
            if(cnt >= mas_inc_cnt){
              rho0_evolve = evolve(rho0_evolve,
                  packets, hamiltonian_offset, rpackets,
                  g, temp_euler,
                  inc, mas_inc_cnt, p.temperature); 
              double result = ::dnpsoup::projectionNorm(rho0_evolve, acq_mat).real();
              const double ratio = result/result_ref;
              results.push_back(make_pair(t, ratio));
              t += static_cast<double>(mas_inc_cnt) * inc;
              cnt -= mas_inc_cnt;
            }
            else {
              rho0_evolve = evolve(rho0_evolve,
                  packets, hamiltonian_offset, rpackets,
                  g, temp_euler,
                  inc, cnt, p.temperature); 
              double result = ::dnpsoup::projectionNorm(rho0_evolve, acq_mat).real();
              const double ratio = result/result_ref;
              results.push_back(make_pair(t, ratio));
              t += static_cast<double>(cnt) * inc;
              cnt = 0;
            }
            mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
          }
          packets.updatePulseSeqComponent(comp);
          cnt = comp_size;
        }
      } while(idx < seq.size() || cnt > 0);

      return results;
    }

    std::vector<std::pair<double, double>> calcFieldProfile(
        const std::vector<Magnet> &fields, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        [[maybe_unused]] int ncores)
    {
      std::vector<std::pair<double, double>> result;
      for(const auto &field : fields){
        //constexpr double ref = 1.0;
        const double ref = calcPowderIntensity(
            field, g, p, spin_sys, "{}", acq_spin, spin_sys_eulers, ncores);
        const double res = calcPowderIntensity(
            field, g, p, spin_sys, pulse_seq_str, acq_spin, spin_sys_eulers, ncores);
        const double ratio = res/ref;
        result.push_back(std::make_pair(field.b0, ratio));
#ifndef NDEBUG
        std::cout << "." << std::flush;
#endif
      }
#ifndef NDEBUG
        std::cout << std::endl;
#endif
      return result;
    }

    std::vector<std::pair<double, double>> calcFieldProfile(
        const Magnet &m, 
        const std::vector<Gyrotron> &emrs,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        [[maybe_unused]] int ncores)
    {
      std::vector<std::pair<double, double>> result;
      for(const auto g : emrs) {
        //constexpr double ref = 1.0;
        const double ref = calcPowderIntensity(
            m, g, p, spin_sys, "{}", acq_spin, spin_sys_eulers, ncores);
        const double res = calcPowderIntensity(
            m, g, p, spin_sys, pulse_seq_str, acq_spin, spin_sys_eulers, ncores);
        const double ratio = res/ref;
        result.emplace_back(make_pair(g.em_frequency, ratio));
      }
      return result;
    }

    double calcPowderIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        [[maybe_unused]] int ncores)
    {
      double result = 0.0;
      const double scaling_factor = 1.0 / static_cast<double>(spin_sys_eulers.size());
      if(ncores == 1) {
        for(const auto &e : spin_sys_eulers){
          auto xtal_intensity = 
            calcIntensity(m, g, p, spin_sys, pulse_seq_str, acq_spin, e);
          result += xtal_intensity * std::sin(e.beta());
        }
      }
      else {
        ::lean::ThreadPool<double> tpool(ncores);
        for(const auto &e : spin_sys_eulers){
          auto task = [=](){
            auto xtal_intensity = 
              calcIntensity(m, g, p, spin_sys, pulse_seq_str, acq_spin, e);
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
      return result;
    }

    MatrixCxDbl evolve(
        const MatrixCxDbl &rho_prev, 
        const MatrixCxDbl &hamiltonian,
        const MatrixCxDbl &hamiltonian_lab,
        const MatrixCxDbl &rotate_mat_super,
        const MatrixCxDbl &rotate_mat_super_inv,
        const std::vector<RelaxationPacket> &rpackets,
        double dt, std::uint64_t cnt, double temperature)
    {
      // static state thermo equilibrium
      auto rho_ss = genRhoEq(hamiltonian_lab, temperature);

      auto rho_ss_super = ::dnpsoup::flatten(rho_ss, 'c');
      if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
        rho_ss_super = rotate_mat_super * rho_ss_super;
      }
      auto [eigenvals, eigenvec] = diagonalizeMat(rho_ss);
      const size_t sz = hamiltonian.nrows() * hamiltonian.ncols();
      auto t1_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
      auto t2_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
      for(const auto &rpacket : rpackets){
        t1_superop += rpacket.genSuperOpT1(eigenvec);
        t2_superop += rpacket.genSuperOpT2(eigenvec);
      }
      auto gamma_super = t1_superop + t2_superop;
#ifdef TEST_LSTSQ
      cout << "gamma_super:\n" << gamma_super << endl;
#endif
      MatrixCxDbl gamma_super_int;
      if(rotate_mat_super.nrows() == 0 || rotate_mat_super_inv.nrows() == 0){
        gamma_super_int = gamma_super;
      } else {
        gamma_super_int = rotate_mat_super * gamma_super * rotate_mat_super_inv;
      }
#ifdef TEST_LSTSQ
      cout << "gamma_super_int:\n" << gamma_super_int << endl;
#endif
      auto h_super = commutationSuperOp(hamiltonian);
      auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;
//#ifndef NDEBUG
//      cout << "[evolve] " << "super_op calculated..." << "\n";
//#endif
#ifdef TEST_LSTSQ
      cout << "rho_ss_super:\n" << rho_ss_super << endl;
#endif
      auto rho_right = gamma_super_int * rho_ss_super;
#ifdef TEST_LSTSQ
      cout << "super_op:\n" << super_op << endl;
      cout << "rho_right:\n" << rho_right << endl;
#endif
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
//#ifndef NDEBUG
//      cout << "[evolve] " << "lstsq() calculated..." << endl;
//#endif
      auto rho_prev_super = ::dnpsoup::flatten(rho_prev, 'c');
      if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
        rho_prev_super = rotate_mat_super * rho_prev_super;
      }

      auto c1 = rho_prev_super - rho_eq_super;
      auto scaling_factor = ::dnpsoup::exp(cxdbl(-1.0 * dt, 0) * super_op);
//#ifndef NDEBUG
//      cout << "[evolve] " << "exp(-Lambda * dt) calculated..." << endl;
//#endif
      scaling_factor = ::dnpsoup::pow(scaling_factor, cnt);
//#ifndef NDEBUG
//      cout << "[evolve] " << "exp(-Lambda * dt * cnt) calculated..." << endl;
//#endif
      auto rho_super = scaling_factor * c1 + rho_eq_super;
      if(rotate_mat_super.nrows() != 0 && rotate_mat_super_inv.nrows() != 0){
        rho_super = rotate_mat_super * rho_super;
      }
      MatrixCxDbl rho_post(rho_prev.nrows(), rho_prev.ncols());

      const auto nrows = rho_post.nrows();
      const auto ncols = rho_post.ncols();
      for(size_t i = 0; i < nrows; ++i){
        for(size_t j = 0; j < ncols; ++j){
          rho_post(i, j) = rho_super(j + i * ncols, 0);
        }
      }
      return rho_post;
    }

    MatrixCxDbl evolve(
        const MatrixCxDbl &rho_prev, 
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
        MatrixCxDbl hamiltonian = packets.genMatrix(euler);
        const double cycle = dt * static_cast<double>(cnt) * g.em_frequency;
        double cycle_int_part;
        double cycle_fraction_part = modf(cycle, &cycle_int_part);
        if (std::abs(cycle_fraction_part) < eps) {
          // placeholders
          MatrixCxDbl rotate_mat, rotate_mat_inv, rotate_mat_super, rotate_mat_super_inv;
          auto rho0_evolve = evolve(rho_prev, 
              hamiltonian, hamiltonian + hamiltonian_offset, 
              rotate_mat_super, rotate_mat_super_inv,  
              rpackets, dt, cnt, temperature);
          return rho0_evolve;
        } else {
          MatrixCxDbl rotate_mat = dnpsoup::exp((cxdbl(0,-2) * dt * pi) * hamiltonian_offset);
          MatrixCxDbl rotate_mat_inv = dnpsoup::exp((cxdbl(0, 2) * dt * pi) * hamiltonian_offset);
          MatrixCxDbl rotate_mat_super = rotationSuperOp(rotate_mat);
          MatrixCxDbl rotate_mat_super_inv = rotationSuperOp(rotate_mat_inv);
          auto rho0_evolve = evolve(rho_prev, 
              hamiltonian, hamiltonian + hamiltonian_offset, 
              rotate_mat_super, rotate_mat_super_inv,  
              rpackets, dt, cnt, temperature);
          return rho0_evolve;
        }
    }
} // namespace DnpRunner
} // namespace dnpsoup

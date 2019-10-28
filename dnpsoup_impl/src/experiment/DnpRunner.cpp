#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/rho_eq.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/common.h"
#include "dnpsoup_core/constants.h"
#include <complex>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;


namespace dnpsoup {
    std::vector<std::vector<double>> DnpRunner::calcEigenValues(
        const Magnet &m,
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const Euler<> &spin_sys_euler
        ) const
    {
      vector<vector<double>> results;
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();

      double t = 0.0;
      Euler<> mas_angle = p.magic_angle;

      packets.setPropertyValue(ValueName::b0, m.b0);
      offset_packets.setPropertyValue(ValueName::b0, m.b0);
      auto spin_ids = spin_sys.getSpinIds(SpinType::e);
      for(const auto &sid : spin_ids){
        packets.setPropertyValue(
            InteractionType::Shielding, sid, ValueName::offset, g.em_frequency);
        offset_packets.setPropertyValue(
            InteractionType::Offset, sid, ValueName::offset, g.em_frequency);
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
      uint64_t idx = 0;

      do {
        pulseseq::Component comp;
        std::tie(comp, idx) = seq.next();
        if(idx == seq.size()) break;

        for(auto &[t_spin, emr] : comp){
          auto obs_id = ObservableId(InteractionType::EMR, t_spin);
          packets.setPropertyValue(
              obs_id, ValueName::freq, emr.freq);
          packets.setPropertyValue(
              obs_id, ValueName::phase, emr.phase);
          packets.setPropertyValue(
              obs_id, ValueName::offset, emr.offset);
        }
        mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);

        MatrixCxDbl hamiltonian = packets.genMatrix(spin_sys_euler * mas_angle);
        MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler * mas_angle);
        //auto eigen_values = ::dnpsoup::eigenVal(hamiltonian);
        auto eigen_values = ::dnpsoup::eigenVal(hamiltonian + hamiltonian_offset);
        vector<double> eigen_values_temp;
        for(size_t i = 0; i < eigen_values.nelements(); ++i){
          eigen_values_temp.push_back(eigen_values(i, 0));
        }
        results.push_back(eigen_values_temp);

        t = t + inc;
      } while(idx < seq.size());

      return results;
    }

    double DnpRunner::calcIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler) const
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
      double t_prev = 0.0;
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
      uint64_t idx = 0;

      MatrixCxDbl hamiltonian = packets.genMatrix(spin_sys_euler * mas_angle);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler * mas_angle);
      auto hamiltonian_lab = hamiltonian + hamiltonian_offset;
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab, p.temperature);
      auto rho0_evolve = rho0_lab;
      std::uint64_t cnt = 1;    /// keep track of identical hamiltonians
      do {
        pulseseq::Component comp;
        std::tie(comp, idx) = seq.next();
        if(idx >= seq.size()) break;

        bool if_hamiltonian_changed = false;
        for(auto &[spin_t, emr] : comp){
          auto obs_id = ObservableId(InteractionType::EMR, spin_t);
          double prev_freq = packets.getPropertyValue(obs_id, ValueName::freq);
          double prev_phase = packets.getPropertyValue(obs_id, ValueName::phase);
          double prev_offset = packets.getPropertyValue(obs_id, ValueName::offset);
          if(!approxEqual<double>(prev_freq, emr.freq, eps)
              || !approxEqual<double>(prev_phase, emr.phase, eps)
              || !approxEqual<double>(prev_offset, emr.offset, eps)){
            packets.setPropertyValue(
                obs_id, ValueName::freq, emr.freq);
            packets.setPropertyValue(
                obs_id, ValueName::phase, emr.phase);
            packets.setPropertyValue(
                obs_id, ValueName::offset, emr.offset);
            if_hamiltonian_changed = true;
          }
          else {
            if_hamiltonian_changed = false;
          }
        }
        t += inc;
        if((t - t_prev) > mas_inc - eps && mas_inc > 0){
          // active rotation
          mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
          t_prev = t;
          if_hamiltonian_changed = true;
        }

        if(!if_hamiltonian_changed){
          ++cnt;
          continue;
        }
        else {
          auto hamiltonian_temp = packets.genMatrix(spin_sys_euler * mas_angle);
          if(cnt > 0){
            const double dt = inc * static_cast<double>(cnt);
            MatrixCxDbl rotate_mat = dnpsoup::exp((cxdbl(0,-2) * dt * pi) * hamiltonian_offset);
            MatrixCxDbl rotate_mat_inv = dnpsoup::exp((cxdbl(0, 2) * dt * pi) * hamiltonian_offset);
            MatrixCxDbl rotate_mat_super = rotationSuperOp(rotate_mat);
            MatrixCxDbl rotate_mat_super_inv = rotationSuperOp(rotate_mat_inv);
            rho0_evolve = evolve(rho0_evolve, 
                hamiltonian, hamiltonian + hamiltonian_offset, 
                rotate_mat_super, rotate_mat_super_inv,  
                rpackets, dt, p.temperature);
          }
          hamiltonian = hamiltonian_temp;
          cnt = 1;
        }
      } while(idx < seq.size());
      if(cnt > 0){
        const double dt = inc * static_cast<double>(cnt);
        MatrixCxDbl rotate_mat 
          = dnpsoup::exp((cxdbl(0,-2) * dt * pi) * hamiltonian_offset);
        MatrixCxDbl rotate_mat_inv 
          = dnpsoup::exp((cxdbl(0, 2) * dt * pi) * hamiltonian_offset);
        MatrixCxDbl rotate_mat_super = rotationSuperOp(rotate_mat);
        MatrixCxDbl rotate_mat_super_inv = rotationSuperOp(rotate_mat_inv);
        rho0_evolve = evolve(rho0_evolve, 
            hamiltonian, hamiltonian + hamiltonian_offset, 
            rotate_mat_super, rotate_mat_super_inv,  
            rpackets, dt, p.temperature);
      }

      double result = ::dnpsoup::projectionNorm(rho0_evolve, acq_mat).real();
      return result;
    }

    std::vector<pair<double, double>> DnpRunner::calcBuildUp(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler) const
    {
      constexpr double eps = std::numeric_limits<double>::epsilon();
      vector<pair<double, double>> results;
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();
      auto rpackets = spin_sys.summarizeRelaxation();
      auto acq_mat = spin_sys.acquireOn(acq_spin);
      double mas_inc = -1;
      if(p.mas_frequency > eps){
        if(p.mas_increment > 0.0){
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
      uint64_t idx = 0;

      MatrixCxDbl hamiltonian = packets.genMatrix(spin_sys_euler * mas_angle);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler * mas_angle);
      MatrixCxDbl rotate_mat 
        = dnpsoup::exp((cxdbl(0,-2.0) * inc * pi) * hamiltonian_offset);
      MatrixCxDbl rotate_mat_inv 
        = dnpsoup::exp((cxdbl(0, 2.0) * inc * pi) * hamiltonian_offset);
      MatrixCxDbl rotate_mat_super = rotationSuperOp(rotate_mat);
      MatrixCxDbl rotate_mat_super_inv = rotationSuperOp(rotate_mat_inv);
      auto hamiltonian_lab = hamiltonian + hamiltonian_offset;
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab, p.temperature);
      auto rho0_evolve = rho0_lab;
      double result = ::dnpsoup::projectionNorm(rho0_evolve, acq_mat).real();
      results.push_back(make_pair(t, result));
      pulseseq::Component comp;
      std::tie(comp, idx) = seq.next();
      while(idx < seq.size()){
        for(auto &[spin_t, emr] : comp){
          auto obs_id = ObservableId(InteractionType::EMR, spin_t);
          packets.setPropertyValue(
              obs_id, ValueName::freq, emr.freq);
          packets.setPropertyValue(
              obs_id, ValueName::phase, emr.phase);
          packets.setPropertyValue(
              obs_id, ValueName::offset, emr.offset);
        }
        t += inc;
        // active rotation
        mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
        hamiltonian = packets.genMatrix(spin_sys_euler * mas_angle);
        rho0_evolve = evolve(rho0_evolve, 
            hamiltonian, hamiltonian + hamiltonian_offset,
            rotate_mat_super, rotate_mat_super_inv,
            rpackets, inc, p.temperature);
        result = ::dnpsoup::projectionNorm(rho0_evolve, acq_mat).real();
        results.push_back(make_pair(t, result));
        std::tie(comp, idx) = seq.next();
      }

      return results;
    }

    std::vector<double> DnpRunner::calcFieldProfile(
        const std::vector<Magnet> &fields, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        [[maybe_unused]] int ncores) const
    {
      std::vector<double> result;
      for(const auto field : fields){
        const double res = this->calcPowderIntensity(
            field, g, p, spin_sys, pulse_seq_str, acq_spin, spin_sys_eulers);
        result.push_back(res);
#ifndef NDEBUG
        std::cout << "." << std::flush;
#endif
      }
#ifndef NDEBUG
        std::cout << std::endl;
#endif
      return result;
    }

    std::vector<double> DnpRunner::calcFieldProfile(
        const Magnet &m, 
        const std::vector<Gyrotron> &emrs,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        [[maybe_unused]] int ncores) const
    {
      std::vector<double> result;
      for(const auto g : emrs) {
        const double res = this->calcPowderIntensity(
            m, g, p, spin_sys, pulse_seq_str, acq_spin, spin_sys_eulers);
        result.push_back(res);
      }
      return result;
    }

    double DnpRunner::calcPowderIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers,
        [[maybe_unused]] int ncores) const
    {
      double result = 0.0;
      const double scaling_factor = 1.0 / static_cast<double>(spin_sys_eulers.size());
      for(const auto &e : spin_sys_eulers){
        auto xtal_intensity = 
          this->calcIntensity(m, g, p, spin_sys, pulse_seq_str, acq_spin, e);
        result += xtal_intensity * std::sin(e.beta());
      }
      result = result * scaling_factor / pi * 4.0;
      return result;
    }

    MatrixCxDbl DnpRunner::evolve(
        const MatrixCxDbl &rho_prev, 
        const MatrixCxDbl &hamiltonian,
        const MatrixCxDbl &hamiltonian_lab,
        const MatrixCxDbl &rotate_mat_super,
        const MatrixCxDbl &rotate_mat_super_inv,
        const std::vector<RelaxationPacket> &rpackets,
        double dt, double temperature) const
    {
      auto rho_eq = genRhoEq(hamiltonian_lab, temperature);
      auto delta_rho = rho_prev - rho_eq;
      auto delta_rho_super = ::dnpsoup::flatten(delta_rho, 'c');
      auto delta_rho_int_super = rotate_mat_super * delta_rho_super;
      auto [eigenvals, eigenvec] = diagonalizeMat(rho_eq);
      const size_t sz = hamiltonian.nrows() * hamiltonian.ncols();
      auto t1_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
      auto t2_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
      for(const auto &rpacket : rpackets){
        t1_superop += rpacket.genSuperOpT1(eigenvec);
        t2_superop += rpacket.genSuperOpT2(eigenvec);
      }
      auto gamma_super = t1_superop + t2_superop;
      auto gamma_super_int = rotate_mat_super * gamma_super * rotate_mat_super_inv;
      auto h_super = commutationSuperOp(hamiltonian);
      auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;
      delta_rho_int_super = dnpsoup::exp((-2.0 * pi * dt) * super_op) * delta_rho_int_super;
      delta_rho_super = rotate_mat_super_inv * delta_rho_int_super;

      //auto d_delta_rho = MatrixCxDbl(rho_eq.nrows(), rho_eq.ncols());
      const auto nrows = delta_rho.nrows();
      const auto ncols = delta_rho.ncols();
      for(size_t i = 0; i < nrows; ++i){
        for(size_t j = 0; j < ncols; ++j){
          rho_eq(i, j) += delta_rho_super(j + i * ncols, 0);
        }
      }
      //return d_delta_rho + rho_eq;
      return rho_eq;
    }
} // namespace dnpsoup

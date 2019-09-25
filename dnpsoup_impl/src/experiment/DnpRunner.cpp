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
      if(p.mas_frequency > 1.0){
        if(p.mas_increment > 0.0){
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

      MatrixCxDbl hamiltonian = packets.genMatrix(spin_sys_euler);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler);
      auto hamiltonian_lab = hamiltonian + hamiltonian_offset;
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab, p.temperature);
      MatrixCxDbl rho0_relative = genRhoEq(hamiltonian, p.temperature);
      auto rho0_offset = rho0_lab - rho0_relative;
      std::uint64_t cnt = 0;    /// keep track of identical hamiltonians
      do {
        pulseseq::Component comp;
        std::tie(comp, idx) = seq.next();

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
        if((t - t_prev) > mas_inc - eps && mas_inc > 0){
          // active rotation
          mas_angle.gamma(t * p.mas_frequency * 2.0 * pi);
          t_prev = t;
        }

        auto hamiltonian_temp = packets.genMatrix(spin_sys_euler * mas_angle);
        if(allclose(hamiltonian_temp, hamiltonian, eps)){
          ++cnt;
          continue;
        } else {
          if(cnt > 0){
            const double dt = inc * static_cast<double>(cnt);
            rho0_relative = evolve(rho0_relative, hamiltonian, rpackets, dt, p.temperature);
          }
          hamiltonian = hamiltonian_temp;
          cnt = 1;
        }
      } while(idx < seq.size());
      if(cnt > 0){
        const double dt = inc * static_cast<double>(cnt);
        rho0_relative = evolve(rho0_relative, hamiltonian, rpackets, dt, p.temperature);
      }

      double result = ::dnpsoup::projectionNorm(acq_mat, rho0_relative + rho0_offset).real();
      return result;
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
        const std::vector<RelaxationPacket> &rpackets,
        double dt, double temperature) const
    {
      auto rho_eq = genRhoEq(hamiltonian, temperature);
      auto delta_rho = ::dnpsoup::flatten(rho_prev - rho_eq, 'c');
      auto [eigenvals, eigenvec] = diagonalizeMat(hamiltonian);
      const size_t sz = hamiltonian.nrows() * hamiltonian.ncols();
      auto t1_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
      auto t2_superop = zeros<cxdbl>(sz, sz); ///< sz: 'size'
      for(const auto &rpacket : rpackets){
        t1_superop += rpacket.genSuperOpT1(eigenvec);
        t2_superop += rpacket.genSuperOpT2(eigenvec);
      }
      auto h_super = commutationSuperOp(hamiltonian);
      auto super_op = complex<double>(0,-1) * h_super + t1_superop + t2_superop;
      auto d_delta_rho_super = exp((2.0 * pi * dt) * super_op) * delta_rho;

      auto d_delta_rho = MatrixCxDbl(rho_eq.nrows(), rho_eq.ncols());
      const size_t nrows = d_delta_rho.nrows();
      for(size_t i = 0; i < d_delta_rho.ncols(); ++i){
        for(size_t j = 0; j < nrows; ++j){
          d_delta_rho(j, i) = d_delta_rho_super(i * nrows + j, 0);
        }
      }
      return d_delta_rho + rho_eq;
    }
} // namespace dnpsoup

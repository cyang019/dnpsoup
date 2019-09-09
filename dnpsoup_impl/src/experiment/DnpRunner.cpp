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
    double DnpRunner::calcIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler) const
    {
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();
      auto rpackets = spin_sys.summarizeRelaxation();
      auto acq_mat = spin_sys.acquireOn(acq_spin);

      packets.setPropertyValue(ValueName::b0, m.b0);
      offset_packets.setPropertyValue(ValueName::b0, m.b0);
      auto spin_ids = spin_sys.getSpinIds(SpinType::e);
      for(const auto &sid : spin_ids){
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
      MatrixCxDbl rho0_lab = genRhoEq(hamiltonian_lab);
      MatrixCxDbl rho0_relative = genRhoEq(hamiltonian);
      auto rho0_offset = rho0_lab - rho0_relative;
      std::uint64_t cnt = 0;    /// keep track of identical hamiltonians
      constexpr double eps = std::numeric_limits<double>::epsilon();
      do {
        pulseseq::Component comp;
        std::tie(comp, idx) = seq.next();

        for(auto &[t, emr] : comp){
          auto sids = spin_sys.getSpinIds(t);
          for(const auto &sid : sids){
            packets.setPropertyValue(
                InteractionType::EMR, sid, ValueName::freq, emr.freq);
            packets.setPropertyValue(
                InteractionType::EMR, sid, ValueName::phase, emr.phase);
            //packets.setPropertyValue(
            //    InteractionType::EMR, sid, ValueName::offset, emr.offset);
          }
        }
        auto hamiltonian_temp = packets.genMatrix(spin_sys_euler);
        if(allclose(hamiltonian_temp, hamiltonian, eps)){
          ++cnt;
          continue;
        } else {
          if(cnt > 0){
            const double dt = inc * static_cast<double>(cnt);
            rho0_relative = evolve(rho0_relative, hamiltonian, rpackets, dt);
          }
          hamiltonian = hamiltonian_temp;
          cnt = 1;
        }
      } while(idx < seq.size());
      if(cnt > 0){
        const double dt = inc * static_cast<double>(cnt);
        rho0_relative = evolve(rho0_relative, hamiltonian, rpackets, dt);
      }

      double result = projection(acq_mat, rho0_relative + rho0_offset);
      return result;
    }

    double DnpRunner::calcPowderIntensity(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const std::vector<Euler<>> &spin_sys_eulers) const
    {
      double result = 0.0;
      const double scaling_factor = 1.0 / static_cast<double>(spin_sys_eulers.size());
      for(const auto &e : spin_sys_eulers){
        auto xtal_intensity = 
          this->calcIntensity(m, g, p, spin_sys, pulse_seq_str, acq_spin, e);
        result += xtal_intensity * std::sin(e.beta());
      }
      result = result * scaling_factor / (4.0 * pi);
      return result;
    }

    MatrixCxDbl DnpRunner::evolve(
        const MatrixCxDbl &rho_prev, 
        const MatrixCxDbl &hamiltonian,
        const std::vector<RelaxationPacket> &rpackets,
        double dt) const
    {
      auto rho_eq = genRhoEq(hamiltonian);
      auto delta_rho = rho_prev - rho_eq;
      auto [eigenvals, eigenvec] = diagonalizeMat(hamiltonian);
      MatrixCxDbl t1_superop = rpacket.genSuperOpT1(eigenvec);
      MatrixCxDbl t2_superop = rpacket.genSuperOpT2(eigenvec);
      auto h_super = commutationSuperOp(hamiltonian);
      auto super_op = complex<double>(0,-1) * h_super - t1_super + t2_super;
      auto d_delta_rho = ::dnpsoup::exp(-super_op * dt) * delta_rho;
      return d_delta_rho + rho_eq;
    }
} // namespace dnpsoup

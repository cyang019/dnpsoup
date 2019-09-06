#include "dnpsoup_core/experiment/DnpRunner.h"
#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/rho_eq.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
#include "dnpsoup_core/errors.h"
#include "dnpsoup_core/common.h"
#include <complex>
#include <sstream>
#include <string>

using namespace std;


namespace dnpsoup {
    double DnpRunner::calcEnhancement(
        const Magnet &m, 
        const Gyrotron &g,
        const Probe &p,
        const SpinSys &spin_sys,
        const std::string &pulse_seq_str,
        const SpinType &acq_spin,
        const Euler<> &spin_sys_euler)
    {
      double result = 0.0;
      auto packets = spin_sys.summarize<DnpExperiment>();
      auto offset_packets = spin_sys.summarizeOffset<DnpExperiment>();
      auto relaxation_packets = spin_sys.summarizeRelaxation<DnpExperiment>();
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
      istringstream iss(spin_sys_euler);
      try{
        iss >> seq;
      }
      catch(const exception &e){
        throw PulseSequenceError(e.what());
      }
      const double inc = seq.getIncrement();
      uint64_t idx = 0;

      MatrixCxDbl hamiltonian_relative = packets.genMatrix(spin_sys_euler);
      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spin_sys_euler);
      auto hamiltonian_lab = hamiltonian_relative + hamiltonian_offset;
      MatrixCxDbl rho0 = genRhoEq(hamiltonian_lab);
      MatrixCxDbl rho0_relative = genRhoEq(hamiltonian_relative);
      auto rho0_offset = rho0 - rho0_relative;
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
            packets.setPropertyValue(
                InteractionType::EMR, sid, ValueName::offset, emr.offset);
          }
        }
        hamiltonian_relative = packets.genMatrix(spin_sys_euler);

        MatrixCxDbl rho1_relative = genRhoEq(hamiltonian_relative);

        // deviation from current equilibrium
        MatrixCxDbl delta_rho = rho0_relative - rho1_relative;

        auto [eigenvals, eigenvec] = diagonalizeMat(hamiltonian_relative);
        MatrixCxDbl t1_superop = genSuperOpT1(eigenvec);
        MatrixCxDbl t2_superop = genSuperOpT2(eigenvec);
        auto h_super = commutationSuperOp(hamiltonian_relative);
        auto super_op = - complex<double>(0,1) * h_super - t1_super + t2_super;
        delta_rho_evolved = ::dnpsoup::exp(-super_op * inc) * delta_rho;

        rho0_relative = delta_rho_evolved + rho1_relative;
      } while(idx < seq.size());

      result = projection(acq_mat, rho0_relative + rho0_offset);
      return result;
    }
} // namespace dnpsoup

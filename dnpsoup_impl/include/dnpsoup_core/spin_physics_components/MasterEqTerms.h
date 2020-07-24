#ifndef DNPSOUP_MASTEREQTERMS_H
#define DNPSOUP_MASTEREQTERMS_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/experiment/experiment_types.h"
#include "dnpsoup_core/experiment/hardware.h"
#include "dnpsoup_core/spin_physics_components/hamiltonian/rho_eq.h"
#include "dnpsoup_core/spinsys/HamiltonianPacket.h"
#include "dnpsoup_core/spinsys/RelaxationPacket.h"
#include "dnpsoup_core/pulseseq/pulse_sequence.h"
#include <vector>
#include <tuple>
#include <cstdint>


namespace dnpsoup {
  /// rho_{t+dt} = E(rho_t - c1) + c1prime
  /// E = sum_i{exp(-L_idt_i)}
  /// for static case
  struct MasterEqTerms {
    MatrixCxDbl E;
    MatrixCxDbl c1;
    MatrixCxDbl c1prime;

    MasterEqTerms() {}
    MasterEqTerms(const MasterEqTerms &);
    MasterEqTerms(MasterEqTerms &&) noexcept;
    MasterEqTerms& operator=(const MasterEqTerms &);
    MasterEqTerms& operator=(MasterEqTerms &&) noexcept;
    MasterEqTerms(MatrixCxDbl &&exp_term, MatrixCxDbl &&c1)
      : E(std::move(exp_term)), c1(std::move(c1)) {}
    MasterEqTerms(MatrixCxDbl &&exp_term, MatrixCxDbl &&c1, MatrixCxDbl &&c1prime)
      : E(std::move(exp_term)), c1(std::move(c1)), c1prime(std::move(c1prime)) {}
    ~MasterEqTerms() {}
  };
  std::ostream& operator<<(std::ostream &os, const MasterEqTerms &);

  MasterEqTerms combineMasterEqTerms(const std::vector<MasterEqTerms> &terms, size_t n);

  /// find an order to calculate sections so that 
  /// its children would have results first.
  std::vector<std::string> findSectionOrder(
    const std::map<std::string, std::unique_ptr<pulseseq::SubSequenceInterface>> *ptr_sections  ///< all nodes
  );

  inline MatrixCxDbl evolve(const MatrixCxDbl &rho_super, const MasterEqTerms &m)
  {
    if (m.c1prime.nelements() > 0) {
      return m.E * (rho_super - m.c1) + m.c1prime;
    } else {
      return m.E * (rho_super - m.c1) + m.c1;
    }
  }

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTermsFromSingletonSeq(
    PacketCollection *packets,
    const std::vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const pulseseq::SubSequenceInterface* ptr_section,
    const std::map<std::string, pulseseq::Component>* ptr_components,
    const std::vector<SpinType> &irradiated,
    const Gyrotron &g,
    const Euler<> &euler,
    double temperature,
    double inc
  );

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTerms(
    PacketCollection *packets,
    const std::vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const PulseSequence &pseq,
    const std::vector<SpinType> &irradiated,
    const Gyrotron &g,
    const Euler<> &euler,
    double temperature,
    double inc
  );

  MasterEqTerms genMasterEqTermsMAS(
    PacketCollection *packets,
    const std::vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    //const Gyrotron &g,  // em freq
    const Euler<> &sample_euler,
    const Euler<> &magic_angle,
    size_t comp_size,
    double temperature,
    double mas_freq,
    double inc,
    double mas_inc,
    size_t cache_size = 256
  );
}

#endif

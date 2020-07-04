#ifndef DNPSOUP_MASTEREQTERMS_H
#define DNPSOUP_MASTEREQTERMS_H

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/experiment/experiment_types.h"
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
    MasterEqTerms(
        double inc);
    ~MasterEqTerms() {}
  };

  MasterEqTerms combineMasterEqTerms(const std::vector<MasterEqTerms> &terms, size_t n);

  /// find an order to calculate sections so that 
  /// its children would have results first.
  std::vector<std::string> findSectionOrder(
    const std::map<std::string, std::unique_ptr<pulseseq::SubSequenceInterface>> *ptr_sections  ///< all nodes
  );

  inline MatrixCxDbl evolve(const MatrixCxDbl &rho_super, const MasterEqTerms &m)
  {
    return m.E * (rho_super - m.c1) + m.c1prime;
  }

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTermsFromSingletonSeq(
    PacketCollection *packets,
    const std::vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const pulseseq::SubSequenceInterface* ptr_section,
    const std::map<std::string, pulseseq::Component>* ptr_components,
    const std::vector<SpinType> &irradiated,
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
    const Euler<> &euler,
    double temperature,
    double inc
  );
}

#endif

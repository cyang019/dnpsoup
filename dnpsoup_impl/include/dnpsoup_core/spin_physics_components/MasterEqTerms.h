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
#include <string>


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
    MasterEqTerms(MatrixCxDbl &&exp_term, MatrixCxDbl &&c1_rhs)
      : E(std::move(exp_term)), c1(std::move(c1_rhs)) {}
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
#ifndef NDEBUG
    if(rho_super.nrows() != m.c1.nrows()) {
      std::string err_msg = "evolve(rho, terms) rho_super and c1 # of rows mismatch: " 
        + std::to_string(rho_super.nrows()) + " vs " + std::to_string(m.c1.nrows());
      std::cout << err_msg << std::endl;
      throw SizeMismatchError(err_msg);
    }
    if(rho_super.ncols() != m.c1.ncols()) {
      std::string err_msg = "evolve(rho, terms) rho_super and c1 # of columns mismatch: " 
        + std::to_string(rho_super.ncols()) + " vs " + std::to_string(m.c1.ncols());
      std::cout << err_msg << std::endl;
      throw SizeMismatchError(err_msg);
    }
#endif
    const auto rho_diff = rho_super - m.c1;
#ifndef NDEBUG
    if (m.E.ncols() != rho_diff.nrows()) {
      std::string err_msg = "evolve(rho, terms) E and rho_diff size mismatch: E "
        + std::to_string(m.E.nrows()) + "x" + std::to_string(m.E.ncols()) 
        + " vs rho_diff " 
        + std::to_string(rho_diff.nrows()) + "x" + std::to_string(rho_diff.ncols());
      std::cout << err_msg << std::endl;
      throw SizeMismatchError(err_msg);
    }
#endif
    if (m.c1prime.nelements() > 0) {
      return m.E * rho_diff + m.c1prime;
    } else {
      return m.E * rho_diff + m.c1;
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
    size_t cache_size = 128
  );
  
  // including customized relaxation
  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTermsFromSingletonSeq(
    PacketCollection *packets,
    const RelaxationPacketCollection &rpackets,
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
    const RelaxationPacketCollection &rpackets,
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
    const RelaxationPacketCollection &rpackets,
    const MatrixCxDbl &ham_offset,
    //const Gyrotron &g,  // em freq
    const Euler<> &sample_euler,
    const Euler<> &magic_angle,
    size_t comp_size,
    double temperature,
    double mas_freq,
    double inc,
    double mas_inc,
    size_t cache_size = 128
  );

  // if the entire mat is nan, just need to check the first element
  inline bool hasNan(const MatrixCxDbl &mat)
  {
    if(mat.nelements() > 0) {
      if(std::isnan(mat(0,0).real()) || std::isnan(mat(0,0).imag())) {
        return true;
      }
    }
    return false;
  }

  inline bool hasNan(const MasterEqTerms &term)
  {
    return hasNan(term.E) || hasNan(term.c1) || hasNan(term.c1prime);
  }

  inline bool hasNanInC1prime(const MasterEqTerms &term)
  {
    return hasNan(term.c1prime);
  }
}

#endif

#include "dnpsoup_core/common.h"
#include "dnpsoup_core/spin_physics_components/MasterEqTerms.h"
#include "dnpsoup_core/pulseseq/ChirpPulse.h"
#include "dnpsoup_core/spin_physics_components/evolve.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

namespace dnpsoup {
  MasterEqTerms::MasterEqTerms(const MasterEqTerms &rhs)
    : E(rhs.E), c1(rhs.c1), c1prime(rhs.c1prime)
  {}

  MasterEqTerms& MasterEqTerms::operator=(const MasterEqTerms &rhs)
  {
    E = rhs.E;
    c1 = rhs.c1;
    c1prime = rhs.c1prime;
    return *this;
  }

  MasterEqTerms::MasterEqTerms(MasterEqTerms &&rhs) noexcept
    : E(std::move(rhs.E)), c1(std::move(rhs.c1)), 
    c1prime(std::move(rhs.c1prime))
  {}

  MasterEqTerms& MasterEqTerms::operator=(MasterEqTerms &&rhs) noexcept
  {
    E = std::move(rhs.E);
    c1 = std::move(rhs.c1);
    c1prime = std::move(rhs.c1prime);
    return *this;
  }

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTermsFromSingletonSeq(
    PacketCollection *packets,
    const std::vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const pulseseq::SubSequenceInterface* ptr_section,
    const map<string, pulseseq::Component>* ptr_components,
    const vector<SpinType> &irradiated,
    const Euler<> &euler,
    double temperature,
    double inc
  )
  {
    pulseseq::Component default_comp;
    MasterEqTerms result;
    for(const auto &t : irradiated){
      default_comp.insert_or_assign(t, pulseseq::EMRadiation());
    }

    using SequenceType = pulseseq::SequenceType;
    const SequenceType sec_type = ptr_section->type();
    switch (sec_type) {
      case SequenceType::DelayType:
        {
          packets->updatePulseSeqComponent(default_comp);
          MatrixCxDbl ham = packets->genMatrix(euler);
          MatrixCxDbl ham_lab = ham + ham_offset;
          auto [h_super, gamma_super, rho_eq_super] = 
            calcSuperOpsForMasterEq(ham, ham_lab,
                rpackets, temperature);
          const auto l_super = calcLambdaSuper(h_super, gamma_super);
          result.E = calcExpEvolve(l_super, inc, ptr_section->size());
          result.c1 = std::move(rho_eq_super);
        }
        break;
      case SequenceType::PulseType:
        {
          packets->updatePulseSeqComponent(default_comp);
          const string comp_name = ptr_section->getNames()[0];
          const auto comp = ptr_components->at(comp_name);
          packets->updatePulseSeqComponent(comp);

          MatrixCxDbl ham = packets->genMatrix(euler);
          MatrixCxDbl ham_lab = ham + ham_offset;
          auto [h_super, gamma_super, rho_eq_super] = 
            calcSuperOpsForMasterEq(ham, ham_lab,
                rpackets, temperature);
          const auto l_super = calcLambdaSuper(h_super, gamma_super);
          result.E = calcExpEvolve(l_super, inc, ptr_section->size());
          result.c1 = std::move(rho_eq_super);
        }
        break;
      case SequenceType::ChirpType:
        {
          auto ptr_chirp = ptr_section->copy();
          const auto seq_size = ptr_chirp->size();
          ptr_chirp->resetSelfIndex();
          std::map<string, std::unique_ptr<pulseseq::SubSequenceInterface>> sections_placeholder;
          std::map<string, pulseseq::Component> components = *ptr_components;
          packets->updatePulseSeqComponent(default_comp);
          bool first_iter = true;
          MatrixCxDbl ham, ham_lab, rho_ss, rho_ss_super, gamma_super, h_super;
          while(true) {
            auto [comp, comp_size, idx] = ptr_chirp->next(
                &components,
                &sections_placeholder
                );
            if (idx >= seq_size) break;
            packets->updatePulseSeqComponent(comp);

            ham = packets->genMatrix(euler);
            ham_lab = ham + ham_offset;
            // static state thermo equilibrium
            rho_ss = genRhoEq(ham_lab, temperature);

            rho_ss_super = ::dnpsoup::flatten(rho_ss, 'c');
            gamma_super = calcGammaSuper(rho_ss, rpackets);
            h_super = commutationSuperOp(ham);
            const auto l_super = calcLambdaSuper(h_super, gamma_super);
            if (first_iter) {
              result.E = calcExpEvolve(l_super, inc, ptr_section->size());
              result.c1 = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
              first_iter = false;
            } else {
              result.E = calcExpEvolve(l_super, inc, ptr_section->size()) * result.E;
            }
          }
          // if at least one iteration was executed
          if (! first_iter) {
            // calculate the last iter
            result.c1prime = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
          }
        }
        break;
      default:
#ifdef VERBOSE
        cout << "Expecting a Pulse, Delay or Chirp, but saw " << toString(sec_type) << endl;
#endif
        break;
    }
    packets->updatePulseSeqComponent(default_comp);
    return make_tuple(result, packets);
  }

  MasterEqTerms combineMasterEqTerms(const std::vector<MasterEqTerms> &terms, size_t n)
  {
    // time complexity: log(n) * (m+1)
    MasterEqTerms result;
    if (terms.size() == 0) return result;

    const size_t cnt = terms.size();
    // in total terms.size() multiples + 0th term
    vector<MatrixCxDbl> term_residuals;
    term_residuals.reserve(cnt);
    // edcba, edcb, edc, ed, e
    for(size_t i = 0; i < cnt; ++i){
      for(size_t j = 0; j < term_residuals.size(); ++j){
        term_residuals[j] = terms[i].E * term_residuals[j];
      }
      term_residuals.push_back(terms[i].E);
    }
    result.E = dnpsoup::pow(term_residuals[0], n);
    result.c1 = terms[0].c1;

    const size_t nrows = result.E.nrows();
    //----------
    //calculate c1prime
    //----------
    // edcba, aedcb, baedc, cbaed, dcbae
    vector<MatrixCxDbl> term_ratios;
    term_ratios.reserve(cnt);
    for(size_t i = 0; i < cnt; ++i) {
      MatrixCxDbl temp = term_residuals[i];
      for(size_t j = 0; j < i; ++j) {
        temp = terms[i].E * temp;
      }
      term_ratios.push_back(temp);
    }
    // geometric sequences
    vector<MatrixCxDbl> geometric_seqs;
    geometric_seqs.reserve(cnt);
#ifndef NDEBUG
    auto handleStatus = [](int status) {
      if(status != 0){
        string err_msg = "lstsq error during combineMasterEqTerms(): ";
        if(status < 0){
          err_msg = "The " + std::to_string(-status) + "-th argument had an illegal value.";
        } else {
          err_msg = "The algorithm for computing the SVD failed to converge;\n";
          err_msg += std::to_string(status) + "off-diagonal elements of an intermediate"
            + "bidiagonal form did not converge to zero.";
        }
        throw CalculationError(err_msg);
      }
    };
#endif
    //first one has n-1 terms
    if (n == 1) {
      geometric_seqs.push_back(zeros<cxdbl>(nrows, nrows));
    } else {
      // c sum = y
      // c = (1 - r)
      // y = residual
      const MatrixCxDbl i0 = identity<cxdbl>(nrows);
      MatrixCxDbl c = i0 - term_ratios[0];
      MatrixCxDbl y = term_residuals[0] * (i0 - dnpsoup::pow(term_ratios[0], n-1));

      auto [sum, status] = matrix::lstsq(c, y);
#ifndef NDEBUG
      handleStatus(status);
#endif
      geometric_seqs.push_back(sum);
    }
    // 1st has n-1 terms, the other have n terms
    for(size_t i = 1; i < cnt; ++i){
      const MatrixCxDbl i0 = identity<cxdbl>(nrows);
      MatrixCxDbl c = i0 - term_ratios[i];
      MatrixCxDbl y = term_residuals[i] * (i0 - dnpsoup::pow(term_ratios[i], n-1));

      auto [sum, status] = matrix::lstsq(c, y);
#ifndef NDEBUG
      handleStatus(status);
#endif
      geometric_seqs.push_back(sum);
    }

    // dynamic equilibrium constants c1prime differences
    vector<MatrixCxDbl> deq_diffs;
    deq_diffs.reserve(cnt);
    size_t idx1 = cnt - 1;
    size_t idx2 = 0;
    for(size_t i = 0; i < cnt; ++i){
      const MasterEqTerms &term_l = terms[idx1];
      const MasterEqTerms &term_r = terms[idx2];
      const MatrixCxDbl& c_l = term_l.c1prime.nelements() == 0 ? term_l.c1 : term_l.c1prime;
      const MatrixCxDbl& c_r = term_r.c1prime.nelements() == 0 ? term_r.c1 : term_r.c1prime;
      deq_diffs.push_back(c_l - c_r);
      idx1 = (idx1 + 1) % cnt;
      idx2 = (idx2 + 1) % cnt;
    }
    result.c1prime = geometric_seqs[0] * deq_diffs[0];
    for(size_t i = 1; i < cnt; ++i){
      result.c1prime += geometric_seqs[i] * deq_diffs[i];
    }

    return result;
  }

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTerms(
    PacketCollection *packets,
    const vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const pulseseq::SubSequenceInterface* ptr_section,
    const map<string, pulseseq::Component>* ptr_components,
    const vector<SpinType> &irradiated,
    const Euler<> &euler,
    double temperature,
    double inc
  )
  {
    MasterEqTerms result;
    using SequenceType = pulseseq::SequenceType;
    const SequenceType sec_type = ptr_section->type();
    switch (sec_type) {
      case SequenceType::DelayType:
      case SequenceType::PulseType:
      case SequenceType::ChirpType:
        {
          return genMasterEqTermsFromSingletonSeq(
              packets, rpackets, ham_offset, ptr_section, ptr_components,
              irradiated, euler, temperature, inc);
        }
        break;
      case SequenceType::SectionType:
        {
          // dfs look for singleton section (pulse/delay/chirp)
        }
        break;
      default:
        break;
    }
    return make_tuple(result, packets);
  }
} // namespace dnpsoup

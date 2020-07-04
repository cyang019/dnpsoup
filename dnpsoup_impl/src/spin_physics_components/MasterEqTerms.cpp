#include "dnpsoup_core/common.h"
#include "dnpsoup_core/pulseseq/SubSequenceInterface.h"
#include "dnpsoup_core/spin_physics_components/MasterEqTerms.h"
#include "dnpsoup_core/pulseseq/ChirpPulse.h"
#include "dnpsoup_core/spin_physics_components/evolve.h"
#include "dnpsoup_core/spin_physics_components/super_op.h"
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <algorithm>

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
#ifndef NDEBUG
        cout << "Expecting a Pulse, Delay or Chirp, but saw " << toString(sec_type) << endl;
        throw CalculationError("Unexpected SequenceType deriving MasterEqTerms for singleton section");
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

  std::vector<std::string> findSectionOrder(
    const std::map<std::string, std::unique_ptr<pulseseq::SubSequenceInterface>> *ptr_sections  ///< all nodes
    )
  {
    using SubSeqInterface = pulseseq::SubSequenceInterface;
    vector<pair<string, int>> ranked_sections;

    for(const auto &[name, ptr_section]: *ptr_sections) {
      int score = 0;
      if(ptr_section->isPure()) {
        score = -1;
        ranked_sections.push_back(make_pair(name, score));
      } else {
        vector<string> children = ptr_section->getNames();
        for(const auto &child : children) {
          const SubSeqInterface* ptr_child = ptr_sections->at(child).get();
          if(!ptr_child->isPure()) {
            ++score;
          }
        }
        ranked_sections.push_back(make_pair(name, score));
      }
    }
    std::sort(ranked_sections.begin(), ranked_sections.end(), 
        [](const pair<string, int> &a, const pair<string, int>&b) {
        return a.second < b.second;
    });
    vector<string> result;
    for(const auto &ranked_s : ranked_sections) {
      result.push_back(ranked_s.first);
    }

    return result;
  }

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTerms(
    PacketCollection *packets,
    const vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const PulseSequence &pseq,
    const vector<SpinType> &irradiated,
    const Euler<> &euler,
    double temperature,
    double inc
  )
  {
    using SubSeqInterface = pulseseq::SubSequenceInterface;
    const map<string, pulseseq::Component>* ptr_components = pseq.getPtrComponents();
    const map<string, std::unique_ptr<SubSeqInterface>>* ptr_sections = pseq.getPtrSections();

    vector<string> sections = pseq.getNames();  ///< in pulseseq order
    vector<string> sections_in_calc_order = findSectionOrder(ptr_sections);
    vector<MasterEqTerms> term_per_section;
    for(const string &name : sections_in_calc_order) {
      const SubSeqInterface *ptr_sec = ptr_sections->at(name).get();
      if(ptr_sec->isPure()){
        
        auto [term, packets_temp] = genMasterEqTermsFromSingletonSeq(
              packets, rpackets, ham_offset, ptr_sec, ptr_components,
              irradiated, euler, temperature, inc);
        packets = packets_temp;
        term_per_section.push_back(std::move(term));
      } else {
        vector<string> children = ptr_sec->getNames();
        vector<MasterEqTerms> child_terms;
        for(const string &child : children) {
          const SubSeqInterface *ptr_child = ptr_sections->at(child).get();
          if(ptr_child->isPure()){
            auto [term, packets_temp] = genMasterEqTermsFromSingletonSeq(
                  packets, rpackets, ham_offset, ptr_child, ptr_components,
                  irradiated, euler, temperature, inc);
            packets = packets_temp;
            child_terms.push_back(std::move(term));
          } else {
            auto found = std::find(sections_in_calc_order.begin(), sections_in_calc_order.end(), child); 
            size_t idx = std::distance(sections_in_calc_order.begin(), found);
#ifndef NDEBUG
            if (idx >= term_per_section.size()) {
              throw CalculationError("MasterEqTerms referenced before calculated.");
            }
#endif
            child_terms.push_back(term_per_section[idx]);
          }
        }
        MasterEqTerms term = combineMasterEqTerms(child_terms, ptr_sec->size());
        term_per_section.push_back(std::move(term));
      }
    }
    vector<MasterEqTerms> finalTerms;
    for(const string &name : sections) {
      auto found = std::find(sections_in_calc_order.begin(), sections_in_calc_order.end(), name);
      size_t idx = std::distance(sections_in_calc_order.begin(), found);
      finalTerms.push_back(term_per_section[idx]);
    }
    MasterEqTerms result = combineMasterEqTerms(finalTerms, 1);
    
    return make_tuple(result, packets);
  }
} // namespace dnpsoup

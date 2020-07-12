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

  std::ostream& operator<<(std::ostream &os, const MasterEqTerms &t)
  {
    os << "E:\n" << t.E << "\n"
       << "c1:\n" << t.c1 << "\n";
    if(t.c1prime.nelements() > 0) {
      os << "c1prime:\n" << t.c1prime << "\n";
    }
    return os;
  }

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTermsFromSingletonSeq(
    PacketCollection *packets,
    const std::vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const pulseseq::SubSequenceInterface* ptr_section,
    const map<string, pulseseq::Component>* ptr_components,
    const vector<SpinType> &irradiated,
    const Gyrotron &g,
    const Euler<> &euler,
    double temperature,
    double inc
  )
  {
    pulseseq::Component default_comp;
    for(const auto &t : irradiated){
      default_comp.insert_or_assign(t, pulseseq::EMRadiation());
    }
#ifndef NDEBUG
    cout << ptr_section->name 
         << " ptr_section->size(): " << ptr_section->size() << endl;
#endif

    MasterEqTerms result;
    using SequenceType = pulseseq::SequenceType;
    const SequenceType sec_type = ptr_section->type();
    switch (sec_type) {
      case SequenceType::DelayType:
        {
#ifndef NDEBUG
          cout << "\tdelay" << "\n";
#endif
          packets->updatePulseSeqComponent(default_comp);
          const MatrixCxDbl ham = packets->genMatrix(euler);
          const MatrixCxDbl ham_lab = ham + ham_offset;
          auto [rotate_mat_super, rotate_mat_super_inv] = calcRotationSuperOps(
              ham_offset, g, inc, ptr_section->size());
          auto [h_super, gamma_super, rho_eq_super] = 
            calcSuperOpsForMasterEq(ham, ham_lab,
                rotate_mat_super, rotate_mat_super_inv,
                rpackets, temperature);
          const auto l_super = calcLambdaSuper(h_super, gamma_super);
          result.E = calcExpEvolve(l_super, inc, ptr_section->size());
          result.c1 = std::move(rho_eq_super);
        }
        break;
      case SequenceType::PulseType:
        {
#ifndef NDEBUG
          cout << "\tpulse" << "\n";
#endif
          packets->updatePulseSeqComponent(default_comp);
          const string comp_name = ptr_section->getNames()[0];
          const auto comp = ptr_components->at(comp_name);
          packets->updatePulseSeqComponent(comp);

          const MatrixCxDbl ham = packets->genMatrix(euler);
          const MatrixCxDbl ham_lab = ham + ham_offset;
          auto [rotate_mat_super, rotate_mat_super_inv] = calcRotationSuperOps(
              ham_offset, g, inc, ptr_section->size());
          auto [h_super, gamma_super, rho_eq_super] = 
            calcSuperOpsForMasterEq(ham, ham_lab,
                rotate_mat_super, rotate_mat_super_inv,
                rpackets, temperature);
          const auto l_super = calcLambdaSuper(h_super, gamma_super);
          result.E = calcExpEvolve(l_super, inc, ptr_section->size());
          result.c1 = std::move(rho_eq_super);
        }
        break;
      case SequenceType::ChirpType:
        {
#ifndef NDEBUG
          cout << "\tchirp " << "\n";
#endif
          auto ptr_chirp = ptr_section->copy();
          const auto seq_size = ptr_chirp->size();
          ptr_chirp->resetSelfIndex();
          std::map<string, std::unique_ptr<pulseseq::SubSequenceInterface>> sections_placeholder;
          std::map<string, pulseseq::Component> components = *ptr_components;
          packets->updatePulseSeqComponent(default_comp);
          bool first_iter = true;
          MatrixCxDbl ham, ham_lab, rho_ss, rho_ss_super, gamma_super, h_super;
          MatrixCxDbl rotate_mat_super, rotate_mat_super_inv;
          while(true) {
            auto [comp, comp_size, idx] = ptr_chirp->next(
                &components,
                &sections_placeholder
                );
            if (idx >= seq_size) break;
            if(rotate_mat_super.nelements() == 0){
              std::tie(rotate_mat_super, rotate_mat_super_inv) = calcRotationSuperOps(
                  ham_offset, g, inc, comp_size);
            }
            packets->updatePulseSeqComponent(comp);

            ham = packets->genMatrix(euler);
            ham_lab = ham + ham_offset;
            std::tie(h_super, gamma_super, rho_ss_super) = 
              calcSuperOpsForMasterEq(ham, ham_lab,
                  rotate_mat_super, rotate_mat_super_inv,
                  rpackets, temperature);
            const auto l_super = calcLambdaSuper(h_super, gamma_super);
            if (first_iter) {
              result.E = calcExpEvolve(l_super, inc, comp_size);
              result.c1 = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
              first_iter = false;
            } else {
              result.E = calcExpEvolve(l_super, inc, comp_size) * result.E;
            }
          }
          // if at least one iteration was executed
          if (!first_iter) {
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
#ifndef NDEBUG
    cout << "combineMasterEqTerms " << terms.size() << " terms, of size " << n << "\n";
#endif
    MasterEqTerms result;
    if (terms.size() == 0) return result;
    else if(terms.size() == 1 && n == 1) return terms[0];

    const size_t cnt = terms.size();
    // in total terms.size() multiples + 0th term
    // matrices
    vector<MatrixCxDbl> term_residuals;
    term_residuals.reserve(cnt);
    size_t idx = cnt - 1;
    MatrixCxDbl temp_term = terms[idx].E;
    term_residuals.push_back(temp_term);
    while(idx != 0) {
      temp_term = temp_term * terms[idx-1].E;
      term_residuals.push_back(temp_term);
      --idx;
    }
    result.E = dnpsoup::pow(temp_term, n);
    result.c1 = terms[0].c1;
    result.c1prime = terms.back().c1prime.nelements() > 0 ? terms.back().c1prime : terms.back().c1;

    const size_t nrows = result.E.nrows();
    //----------
    //calculate c1prime
    //----------
    vector<MatrixCxDbl> c_diffs;
    c_diffs.reserve(cnt);
    for(size_t i = cnt; i > 1; --i) {
      const MatrixCxDbl &c_m = terms[i-1].c1;
      const MatrixCxDbl &c_mn1prime = terms[i-2].c1prime.nelements() > 0 ? terms[i-2].c1prime : terms[i-2].c1;
      c_diffs.push_back(c_mn1prime - c_m);
    }
    for(size_t i = 0; i < c_diffs.size(); ++i){
      result.c1prime += term_residuals[i] * c_diffs[i];
    }
    if (n == 1) return result;
    // result.E = temp_term^n
    const MatrixCxDbl i0 = identity<cxdbl>(nrows);
    const MatrixCxDbl c = i0 - temp_term;
    const MatrixCxDbl c1_diff = result.c1prime - result.c1;
    const MatrixCxDbl sum = ::matrix::geometricSum(temp_term, n-1) * c1_diff;
    result.c1prime += sum;

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
#ifndef NDEBUG
    cout << "sections in calculation order: \n\t";
    for(const auto &name: result) {
      cout << name << " \t";
    }
    cout << endl;
#endif

    return result;
  }

  std::tuple<MasterEqTerms, PacketCollection*>
  genMasterEqTerms(
    PacketCollection *packets,
    const vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    const PulseSequence &pseq,
    const vector<SpinType> &irradiated,
    const Gyrotron &g,
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
              irradiated, g, euler, temperature, inc);
        packets = packets_temp;
        term_per_section.push_back(std::move(term));
      } else {
        vector<string> children = ptr_sec->getNames();
        vector<MasterEqTerms> child_terms;
        for(const string &child : children) {
          const SubSeqInterface *ptr_child = ptr_sections->at(child).get();
          if(ptr_child->isPure()){
            auto found = std::find(sections_in_calc_order.begin(), sections_in_calc_order.end(), child); 
            size_t idx = std::distance(sections_in_calc_order.begin(), found);
            if(idx < term_per_section.size()) {
              child_terms.push_back(term_per_section[idx]);
            } else {
              auto [term, packets_temp] = genMasterEqTermsFromSingletonSeq(
                    packets, rpackets, ham_offset, ptr_child, ptr_components,
                    irradiated, g, euler, temperature, inc);
              packets = packets_temp;
              child_terms.push_back(std::move(term));
            }
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
#ifndef NDEBUG
        cout << "terms to combine: ";
        for(const string &child : children) {
          cout << child << " ";
        }
        cout << "\n";
#endif
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
#ifndef NDEBUG
        cout << "terms to combine: ";
        for(const string &child : sections) {
          cout << child << " ";
        }
        cout << "\n";
#endif
    MasterEqTerms result = combineMasterEqTerms(finalTerms, 1);
    
    return make_tuple(result, packets);
  }
} // namespace dnpsoup

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
#include <cmath>

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
#ifndef NDEBUG
            if (isnan(result.c1prime(0,0).real())) {
              cout << "calcRhoDynamicEq() resulted in NAN.\ninputs:\n"
                   << "h_super:\n" << h_super
                   << "gamma_super:\n" << gamma_super
                   << "rho_ss_super:\n" << rho_ss_super << endl;
            }
#endif
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
    if (terms.size() == 0 || n == 0) return result;
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
#ifndef NDEBUG
    if(hasNanInC1prime(result)) {
      cout << "[WARNING] combineMasterEqTerms() c1prime is nan." << endl;
    }
#endif

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
#ifndef NDEBUG
      if(hasNanInC1prime(result)) {
        cout << "[WARNING] combineMasterEqTerms() adding residuals saw nan at [" << i << "]:\n"
             << "c_diffs[i]:\n" << c_diffs[i] << endl;
      }
#endif
    }
    if (n == 1) return result;
    // result.E = temp_term^n
    const MatrixCxDbl i0 = identity<cxdbl>(nrows);
    const MatrixCxDbl c = i0 - temp_term;
    const MatrixCxDbl c1_diff = result.c1prime - result.c1;
#ifndef NDEBUG
    if (isnan(c1_diff(0,0).real()) || isnan(c1_diff(0,0).imag())) {
      cout << "[Warning] combineMasterEqTerms() c1_diff is nan." << endl;
    }
#endif
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

#ifndef NDEBUG
    if(std::isnan(result.c1(0,0).real())) {
      cout << "[WARNING] genMasterEqTerms result.c1 realpart has nan..." << endl;
    }
    if(std::isnan(result.c1(0,0).imag())) {
      cout << "[WARNING] genMasterEqTerms result.c1 imagpart has nan..." << endl;
    }
    if(std::isnan(result.c1prime(0,0).real())) {
      cout << "[WARNING] genMasterEqTerms result.c1prime realpart has nan..." << endl;
    }
    if(std::isnan(result.c1prime(0,0).imag())) {
      cout << "[WARNING] genMasterEqTerms result.c1prime imagpart has nan..." << endl;
    }
    if(std::isnan(result.E(0,0).real())) {
      cout << "[WARNING] genMasterEqTerms result.E realpart has nan..." << endl;
    }
    if(std::isnan(result.E(0,0).imag())) {
      cout << "[WARNING] genMasterEqTerms result.E imagpart has nan..." << endl;
    }
#endif
    
    return make_tuple(result, packets);
  }

  MasterEqTerms genMasterEqTermsMAS(
    PacketCollection *packets,
    const std::vector<RelaxationPacket> &rpackets,
    const MatrixCxDbl &ham_offset,
    //const Gyrotron &g,  // em freq
    const Euler<> &sample_euler,
    const Euler<> &magic_angle,
    size_t comp_size, // component size
    double temperature,
    double mas_freq,
    double inc,
    double mas_inc, // step size
    size_t cache_size // chunk size to combine
  )
  {
    size_t step_size = static_cast<size_t>(std::round(mas_inc/inc));
    step_size += (step_size == 0);
    // ------------------------------------
    constexpr double gamma_step_min = GAMMA_STEP_MIN;    ///< 0.1us
    // ------------------------------------
    size_t gamma_step_size = static_cast<size_t>(std::round(gamma_step_min/mas_inc));
    gamma_step_size += (gamma_step_size == 0);
#ifndef NDEBUG
    cout << "gamma_step_size: " << gamma_step_size << endl;
#endif
    
    size_t cnt_in_one_cycle = 
      static_cast<size_t>(std::round(1.0/(mas_freq * inc)));
    cnt_in_one_cycle += (cnt_in_one_cycle == 0);
    const size_t cnt_part = comp_size % cnt_in_one_cycle;
    const size_t n_rotor_period = comp_size / cnt_in_one_cycle;
#ifndef NDEBUG
    cout << "component size: " << comp_size << ", which is "
         << n_rotor_period << " rotor periods." << endl;
#endif
    //size_t gamma_step_size = static_cast<size_t>(static_cast<double>(cnt_in_one_cycle) * MAS_L_RATIO);
    //gamma_step_size += (gamma_step_size == 0);

    MatrixCxDbl gamma_super;
    //MatrixCxDbl rotate_mat_super, rotate_mat_super_inv;
    //std::tie(rotate_mat_super, rotate_mat_super_inv) = calcRotationSuperOps(
    //    ham_offset, g, inc, step_size);
    auto calcTerms = [&](const Euler<> &euler, size_t sz, bool calc_gamma) {
      const MatrixCxDbl ham = packets->genMatrix(euler);
      const MatrixCxDbl ham_lab = ham + ham_offset;
      const auto rho_ss = genRhoEq(ham_lab, temperature);
      const auto rho_ss_super = ::dnpsoup::flatten(rho_ss, 'c');
      //gamma_super = calcGammaSuper(rho_ss, rpackets);
      if (calc_gamma) {
        gamma_super = calcGammaSuper(rho_ss, rpackets);
      }
      const auto h_super = commutationSuperOp(ham);
      //auto super_op = complex<double>(0,1.0) * h_super + gamma_super_int;
      auto rho_eq_super = calcRhoDynamicEq(h_super, gamma_super, rho_ss_super);
      const auto l_super = calcLambdaSuper(h_super, gamma_super);
      MatrixCxDbl exp_term = calcExpEvolve(l_super, inc, sz);
      return MasterEqTerms(std::move(exp_term), std::move(rho_eq_super));
    };
    std::vector<MasterEqTerms> terms;
    terms.reserve(cache_size + 1);  ///< slot zero reserved for final result
    MasterEqTerms last_term;
    Euler<> temp_magic_angle = magic_angle;
#ifndef NDEBUG
    cout << "rotor period: " << cnt_in_one_cycle << endl;
    cout << "\nentering loop" << endl;
    cout << "cnt_part: " << cnt_part << endl;
#endif

    // euler angle gamma
    auto calcGamma = [](double dt, double mas_freq, double gamma0){
      constexpr double TwoPi = 2.0 * dnpsoup::pi;
      double new_gamma = TwoPi * dt * mas_freq + gamma0;
      while (new_gamma > TwoPi) {
        new_gamma -= TwoPi;
      }
      return new_gamma;
    };

    size_t part1 = cnt_part;
    double temp_t = 0;
    size_t mas_idx = 0;
    while(part1 > 0) {
			//cout << "[euler]: " << temp_magic_angle << endl;
      if(terms.size() >= cache_size) {
        const auto combined_terms = combineMasterEqTerms(terms, 1);
#ifndef NDEBUG
        if(hasNan(combined_terms)) {
          cout << "[WARNING] genMasterEqTermsMAS() saw nan when reducing cache." << endl;
        }
#endif
        terms.clear();
        terms.push_back(std::move(combined_terms));
      }
      if(part1 >= step_size) {
        const auto euler = temp_magic_angle * sample_euler;
        const bool need_to_calc_gamma = (mas_idx % gamma_step_size == 0);
        terms.push_back(calcTerms(euler, step_size, need_to_calc_gamma));

        ++mas_idx;
        part1 -= step_size;
        temp_t += static_cast<double>(step_size) * inc;
        const double new_gamma = calcGamma(temp_t, mas_freq, magic_angle.gamma());
        //const double new_gamma = TwoPi * temp_t * mas_freq + magic_angle.gamma();
        temp_magic_angle.gamma(new_gamma);
      } else {  // if(part1 < step_size)
        const auto euler = temp_magic_angle * sample_euler;
        constexpr bool need_to_calc_gamma = true;
        terms.push_back(calcTerms(euler, part1, need_to_calc_gamma));

        temp_t += static_cast<double>(part1) * inc;
        //const double new_gamma = TwoPi * temp_t * mas_freq + magic_angle.gamma();
        const double new_gamma = calcGamma(temp_t, mas_freq, magic_angle.gamma());
        temp_magic_angle.gamma(new_gamma);
        part1 = 0;
      }
    }
    if (terms.size() > 0) {
      last_term = combineMasterEqTerms(terms, 1);
#ifndef NDEBUG
      if(hasNan(last_term)) {
        cout << "[WARNING] genMasterEqTermsMAS() last_term saw nan." << endl;
      }
#endif
      if(terms.size() > 1 && n_rotor_period > 0) {
        terms.clear();
        terms.push_back(last_term);
      }
    }
#ifndef NDEBUG
    cout << "\nfinish last term\n";
#endif
    if (n_rotor_period > 0) {
      size_t part2 = cnt_in_one_cycle;
      mas_idx = 0;
      while(part2 > cnt_part) {
				//cout << "[euler]: " << temp_magic_angle << endl;
        if(terms.size() >= cache_size) {
          const auto combined_terms = combineMasterEqTerms(terms, 1);
          terms.clear();
          terms.push_back(std::move(combined_terms));
        }
        if(part2 >= step_size + cnt_part) {
          const auto euler = temp_magic_angle * sample_euler;
          const bool need_to_calc_gamma = (mas_idx % gamma_step_size == 0);
          terms.push_back(calcTerms(euler, step_size, need_to_calc_gamma));

          ++mas_idx;
          part2 -= step_size;
          temp_t += static_cast<double>(step_size) * inc;
          const double new_gamma = calcGamma(temp_t, mas_freq, magic_angle.gamma());
          //const double new_gamma = TwoPi * temp_t * mas_freq + magic_angle.gamma();
          temp_magic_angle.gamma(new_gamma);
        } else {  // if(part2 >= step_size + cnt_part)
          const auto euler = temp_magic_angle * sample_euler;
          constexpr bool need_to_calc_gamma = true;
          terms.push_back(calcTerms(euler, part2 - cnt_part, need_to_calc_gamma));

          temp_t += static_cast<double>(part2 - cnt_part) * inc;
          const double new_gamma = calcGamma(temp_t, mas_freq, magic_angle.gamma());
          //const double new_gamma = TwoPi * temp_t * mas_freq + magic_angle.gamma();
          temp_magic_angle.gamma(new_gamma);
          part2 = cnt_part;
        }
      }
    }
#ifndef NDEBUG
    cout << "\nfinish first term" << endl;
#endif
    if(n_rotor_period == 0) {
      return last_term;
    } else if(cnt_part == 0) {
      return combineMasterEqTerms(terms, n_rotor_period);
    } else {
      MasterEqTerms first_term = combineMasterEqTerms(terms, n_rotor_period);
      vector<MasterEqTerms> final_terms;
      final_terms.push_back(std::move(first_term));
      final_terms.push_back(std::move(last_term));
      return combineMasterEqTerms(final_terms, 1);
    }
  }
} // namespace dnpsoup

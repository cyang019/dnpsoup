#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <vector>
#include <complex>
#include <iostream>
#include <sstream>
#include <cmath>

namespace {
    using SpinSys = dnpsoup::SpinSys;
    using SpinType = dnpsoup::SpinType;

    TEST(TestDnpsoup, EvolutionOfIz){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0, true);
      spins.addSpin(2, SpinType::H, 0.7, 0.0, 0.7, true);
      auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      spins.setShielding(dnpsoup::SpinId(1), 2.00252, 2.00252, 2.00252, euler);
      spins.irradiateOn(SpinType::e);
      //spins.setT1(dnpsoup::SpinId(1), 1.0e-3);
      //spins.setT1(dnpsoup::SpinId(2), 1.0);
      //spins.setT2(dnpsoup::SpinId(1), 1.0e-6);
      //spins.setT2(dnpsoup::SpinId(2), 2.0e-3);
      //spins.setCsa(dnpsoup::SpinId(2), 0.0, 0.0, 0.0, euler);

      auto rpackets = spins.summarizeRelaxation();
      double b0 = 9.36932070693522;
      double em_freq = 263.0e9;
      constexpr double dt = 4.0e-8;
      constexpr size_t cnt = 25;
      constexpr double temperature = 100.0;

      auto collection = spins.summarize<dnpsoup::DnpExperiment>();
      auto offset_collection = spins.summarizeOffset<dnpsoup::DnpExperiment>();
      collection.setPropertyValue(dnpsoup::ValueName::b0, b0);
      auto emr_id = dnpsoup::ObservableId(dnpsoup::InteractionType::EMR, dnpsoup::SpinType::e);

      auto spin_ids = spins.getSpinIds(SpinType::e);
      for(const auto &sid : spin_ids){
        collection.setPropertyValue(
            dnpsoup::InteractionType::Shielding, sid, dnpsoup::ValueName::offset, em_freq);
        offset_collection.setPropertyValue(
            dnpsoup::InteractionType::Offset, sid, dnpsoup::ValueName::offset, em_freq);
      }
      auto euler_spins = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      auto hamiltonian0 = collection.genMatrix(euler_spins);
      auto hamiltonian0_offset = offset_collection.genMatrix(euler_spins);

      collection.setPropertyValue(emr_id, dnpsoup::ValueName::freq, 1.0e6);
      auto hamiltonian = collection.genMatrix(euler_spins);
      std::cout << "Hamiltonian0:\n" << hamiltonian0 << "\n";
      std::cout << "Hamiltonian:\n" << hamiltonian << "\n";
      dnpsoup::MatrixCxDbl rho_eq = dnpsoup::kron(
          dnpsoup::spin<dnpsoup::Z>(2), dnpsoup::spin<dnpsoup::OperatorType::Identity>(2));
      auto detection = dnpsoup::kron(
          dnpsoup::spin<dnpsoup::OperatorType::Identity>(2), dnpsoup::spin<dnpsoup::Z>(2));
      auto detection0 = dnpsoup::kron(
          dnpsoup::spin<dnpsoup::Z>(2),
          dnpsoup::spin<dnpsoup::OperatorType::Identity>(2));
      std::cout << "Detection:\n" << detection << "\n";
      std::cout << "Remaining:\n" << detection0 << "\n";

      // ================================
      // manual evolution
      std::vector<std::complex<double>> results;
      std::vector<std::complex<double>> results0;
      auto rho_flatten = dnpsoup::flatten(rho_eq, 'c');
      auto hamiltonian_super = dnpsoup::commutationSuperOp(hamiltonian);
      //std::cout << "Hamiltonian:\n" << hamiltonian << std::endl;
      //std::cout << "SuperOperator:\n" << hamiltonian_super << std::endl;
      const auto mat_evolve = dnpsoup::exp((2 * dnpsoup::pi * dnpsoup::cxdbl(0,-1) * dt) * hamiltonian_super);
      //std::cout << "Evolution Super Operator:\n" << mat_evolve << std::endl;
      for(size_t i = 0; i < cnt; ++i){
        auto result = dnpsoup::projectionNorm(rho_eq, detection);
        auto result0 = dnpsoup::projectionNorm(rho_eq, detection0);
        results.push_back(result);
        results0.push_back(result0);

        rho_flatten = mat_evolve * rho_flatten;
        for(size_t i = 0; i < rho_eq.nrows(); ++i){
          for(size_t j = 0; j < rho_eq.ncols(); ++j){
            rho_eq(i, j) = rho_flatten(i * rho_eq.ncols() + j, 0);
          }
        }
      }
      for(size_t i = 0; i < results.size(); ++i){
        std::cout << results[i].real() << ", " << results[i].imag() << "  \t"
                  << results0[i].real() << ", " << results0[i].imag() << "\n";
      }
      std::cout << std::endl;
      // ================================
      std::cout << " ================================\n";
      
      dnpsoup::DnpRunner r;
      std::vector<std::complex<double>> results_r;
      std::vector<std::complex<double>> results0_r;
      rho_eq = dnpsoup::kron(
          dnpsoup::spin<dnpsoup::Z>(2), dnpsoup::spin<dnpsoup::OperatorType::Identity>(2));
      auto d = dnpsoup::exp(dnpsoup::cxdbl(0, -1.0 * dt) * hamiltonian0_offset);
      auto d_inv = dnpsoup::exp(dnpsoup::cxdbl(0, 1.0 * dt) * hamiltonian0_offset);
      auto d_super = dnpsoup::rotationSuperOp(d);
      auto d_super_inv = dnpsoup::rotationSuperOp(d_inv);
      for(size_t i = 0; i < cnt; ++i){
        rho_eq = r.evolve(rho_eq, 
            hamiltonian, hamiltonian0_offset + hamiltonian,
            d_super, d_super_inv,
            rpackets, dt, temperature); 
        auto result = dnpsoup::projectionNorm(rho_eq, detection);
        auto result0 = dnpsoup::projectionNorm(rho_eq, detection0);
        results_r.push_back(result);
        results0_r.push_back(result0);
      }
      for(size_t i = 0; i < results_r.size(); ++i){
        std::cout << results_r[i].real() << ", " << results_r[i].imag() << "  \t"
                  << results0_r[i].real() << ", " << results0_r[i].imag() << "\n";
      }
      std::cout << std::endl;

      // ================================
      std::cout << " ================================\n";
      auto rho_eq0_lab = dnpsoup::genRhoEq(hamiltonian0 + hamiltonian0_offset, temperature);
      std::vector<dnpsoup::cxdbl> results_2;
      std::vector<dnpsoup::cxdbl> results0_2;
      rho_eq = rho_eq0_lab;
      std::cout << "Rho_eq\n" << rho_eq << std::endl;
      
      for(size_t i = 0; i < cnt; ++i){
        auto result = dnpsoup::projectionNorm(rho_eq, detection);
        auto result0 = dnpsoup::projectionNorm(rho_eq, detection0);
        rho_eq = r.evolve(rho_eq, 
            hamiltonian, hamiltonian0_offset + hamiltonian,
            d_super, d_super_inv,
            rpackets, dt, temperature); 
        results_2.push_back(result);
        results0_2.push_back(result0);
      }
      std::cout << std::setprecision(15);
      for(size_t i = 0; i < results_2.size(); ++i){
        std::cout << results_2[i].real() << ", " << results_2[i].imag() << "  \t"
                  << results0_2[i].real() << ", " << results0_2[i].imag() << "\n";
      }
      std::cout << std::endl;

    }
} // namespace dnpsoup



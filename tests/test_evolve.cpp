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
    using namespace std;
    using namespace dnpsoup;

    TEST(TestDnpsoup, EvolutionFromSzToIz){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0, true);
      spins.addSpin(2, SpinType::H, 1.2, 0.0, 1.2, true);
      double distance = sqrt(1.2*1.2 + 1.2*1.2);
      auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      spins.setShielding(dnpsoup::SpinId(1), 2.00252, 2.00252, 2.00252, euler);
      spins.irradiateOn(SpinType::e);
      //spins.setT1(dnpsoup::SpinId(1), 1.0e-3);
      //spins.setT1(dnpsoup::SpinId(2), 1.0);
      //spins.setT2(dnpsoup::SpinId(1), 1.0e-6);
      //spins.setT2(dnpsoup::SpinId(2), 2.0e-3);

      auto rpackets = spins.summarizeRelaxation();
      double b0 = 9.36932070693522;
      double em_freq = 263.0e9;
      double dt = 1000.0/em_freq;
      constexpr size_t cnt = 1e4;
      constexpr double temperature = 100.0;

      auto collection = spins.summarize<dnpsoup::DnpExperiment>();
      auto offset_collection = spins.summarizeOffset<dnpsoup::DnpExperiment>();
      collection.setPropertyValue(dnpsoup::ValueName::b0, b0);
      offset_collection.setPropertyValue(dnpsoup::ValueName::b0, b0);
      auto b0_val = collection.getPropertyValue(InteractionType::Shielding, SpinId(1), ValueName::b0);
      auto b0_val2 = offset_collection.getPropertyValue(InteractionType::Offset, SpinId(1), ValueName::b0);
      auto b0_val3 = collection.getPropertyValue(InteractionType::Csa, SpinId(2), ValueName::b0);
      auto d_1_2 = collection.getPropertyValue(InteractionType::Dipole, SpinId(1), SpinId(2), ValueName::distance);
      ASSERT_DOUBLE_EQ(b0, b0_val);
      ASSERT_DOUBLE_EQ(b0, b0_val2);
      ASSERT_DOUBLE_EQ(b0, b0_val3);
      ASSERT_DOUBLE_EQ(distance, d_1_2);
      auto emr_id = dnpsoup::ObservableId(dnpsoup::InteractionType::EMR, dnpsoup::SpinType::e);

      auto spin_ids = spins.getSpinIds(SpinType::e);
      ASSERT_EQ(1u, spin_ids.size());
      ASSERT_EQ(SpinId(1), spin_ids[0]);
      for(const auto &sid : spin_ids){
        collection.setPropertyValue(
            dnpsoup::InteractionType::Shielding, sid, dnpsoup::ValueName::offset, em_freq);
        offset_collection.setPropertyValue(
            dnpsoup::InteractionType::Offset, sid, dnpsoup::ValueName::offset, em_freq);
      }
      auto euler_spins = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      auto hamiltonian0 = collection.genMatrix(euler_spins);
      auto hamiltonian0_offset = offset_collection.genMatrix(euler_spins);
      auto oids = collection.getObservableIds();
      for(const auto &oid : oids){
        std::cout << oid.get() << "\n";
        const auto& packet = collection.getPacket(oid);
        std::cout << packet.genMatrix() << "\n";
      }

      collection.setPropertyValue(emr_id, dnpsoup::ValueName::freq, 1.0e6);
      auto hamiltonian = collection.genMatrix(euler_spins);
      std::cout << "Hamiltonian0:\n" << hamiltonian0 << "\n";
      std::cout << "Hamiltonian:\n" << hamiltonian << "\n";
      std::cout << "Hamiltonian_offset:\n" << hamiltonian0_offset << "\n";
      auto SzIx = dnpsoup::kron(spin<Z>(2), spin<X>(2));
      std::cout << "SzIx:\n" << SzIx << "\n";
      auto Sx_e = dnpsoup::kron(
          dnpsoup::spin<dnpsoup::OperatorType::X>(2), dnpsoup::spin<dnpsoup::OperatorType::Identity>(2));
      std::cout << "Sx_e:\n" << Sx_e << "\n";

      auto detection = dnpsoup::kron(
          dnpsoup::spin<dnpsoup::OperatorType::Identity>(2), dnpsoup::spin<dnpsoup::Z>(2));
      auto detection0 = dnpsoup::kron(
          dnpsoup::spin<dnpsoup::Z>(2),
          dnpsoup::spin<dnpsoup::OperatorType::Identity>(2));
      auto acq_mat = spins.acquireOn(dnpsoup::SpinType::H);
      auto acq_mat0 = spins.acquireOn(dnpsoup::SpinType::e);
      ASSERT_TRUE(dnpsoup::allclose(acq_mat, detection, 1.0e-14));
      ASSERT_TRUE(dnpsoup::allclose(acq_mat0, detection0, 1.0e-14));
      std::cout << "Detection:\n" << detection << "\n";
      std::cout << "Remaining:\n" << detection0 << "\n";
      std::cout << std::endl;

      auto rho_eq0_lab = dnpsoup::genRhoEq(hamiltonian0 + hamiltonian0_offset, temperature);
      std::cout << "rho_eq0:\n" << rho_eq0_lab;
      std::vector<dnpsoup::cxdbl> results;
      std::vector<dnpsoup::cxdbl> results0;

      auto Ux = dnpsoup::exp((dnpsoup::cxdbl(0, -1.0e-1) * 2.0 * dnpsoup::pi) * Sx_e);
      auto Ux_inv = dnpsoup::exp((dnpsoup::cxdbl(0, 1.0e-1) * 2.0 * dnpsoup::pi) * Sx_e);

      auto rho = rho_eq0_lab;
      auto rho_lab = genRhoEq(hamiltonian + hamiltonian0_offset, temperature);
      
      //dnpsoup::DnpRunner r;
      auto d = dnpsoup::diag_exp((dnpsoup::cxdbl(0,-1.0) * dt * 2.0 * dnpsoup::pi) * hamiltonian0_offset);
      auto d_inv = dnpsoup::diag_exp((dnpsoup::cxdbl(0,1.0) * dt * 2.0 * dnpsoup::pi) * hamiltonian0_offset);
      std::cout << "D:\n" << d << "\n";
      std::cout << "D^-1:\n" << d_inv << "\n";
      auto d_super = dnpsoup::rotationSuperOp(d);
      auto d_super_inv = dnpsoup::rotationSuperOp(d_inv);
      for(size_t i = 0; i < cnt; ++i){
        auto result = dnpsoup::projectionNorm(rho, detection);
        auto result0 = dnpsoup::projectionNorm(rho, detection0);
        results.push_back(result);
        results0.push_back(result0);

        //rho_eq = r.evolve(rho_eq, 
        //    hamiltonian, hamiltonian0_offset + hamiltonian,
        //    d_super, d_super_inv,
        //    rpackets, dt, temperature); 
        auto rho_eq_temp = genRhoEq(hamiltonian + hamiltonian0_offset, temperature);
        auto rho_eq_r = d * rho_eq_temp * d_inv;
        auto rho_r = d * rho * d_inv;
        auto rho_r_super = dnpsoup::flatten<cxdbl>(rho_r, 'c');

        auto h_super = commutationSuperOp(hamiltonian);
        auto U_super = dnpsoup::exp((cxdbl(0,-2.0) * dt * dnpsoup::pi) * h_super);
        rho_r_super = U_super * rho_r_super;
        for(size_t i = 0; i < rho_eq_temp.nrows(); ++i){
          for(size_t j = 0; j < rho_eq_temp.ncols(); ++j){
            rho_r(i,j) = rho_r_super(i * rho_eq_temp.ncols() + j, 0);
          }
        }
        rho = d_inv * rho_r * d;
      }
      auto result = dnpsoup::projectionNorm(rho, detection);
      auto result0 = dnpsoup::projectionNorm(rho, detection0);
      results.push_back(result);
      results0.push_back(result0);
      auto ss = std::cout.precision();
      std::cout << std::setprecision(15);
      for(size_t i = 0; i < results.size(); ++i){
        if(i % 100 == 0){
          std::cout << results[i].real() << ", " << results[i].imag() << "  \t"
                    << results0[i].real() << ", " << results0[i].imag() << "\n";
        }
      }
      std::cout << std::setprecision(ss);
      std::cout << std::endl;
    }
} // namespace dnpsoup



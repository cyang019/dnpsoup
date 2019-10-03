#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;
    using cxdbl = dnpsoup::cxdbl;

    TEST(TestDnpsoup, flatten){
      MatrixD mat = {{1.0,2.0,3.0}, {4.0,5.0,6.0}, {7.0,8.0,9.0}};
      MatrixD flatted = dnpsoup::flatten(mat, 'c');
      ASSERT_DOUBLE_EQ(1.0, flatted(0,0));
      ASSERT_DOUBLE_EQ(2.0, flatted(1,0));
      ASSERT_DOUBLE_EQ(3.0, flatted(2,0));
      ASSERT_DOUBLE_EQ(4.0, flatted(3,0));
      ASSERT_DOUBLE_EQ(5.0, flatted(4,0));
      ASSERT_DOUBLE_EQ(6.0, flatted(5,0));
      ASSERT_DOUBLE_EQ(7.0, flatted(6,0));
      ASSERT_DOUBLE_EQ(8.0, flatted(7,0));
      ASSERT_DOUBLE_EQ(9.0, flatted(8,0));

      MatrixD mat2(3,3);
      for(size_t i = 0; i < mat2.nrows(); ++i){
        for(size_t j = 0; j < mat2.ncols(); ++j){
          mat2(i,j) = flatted(i * mat2.ncols() + j, 0);
        }
      }
      ASSERT_DOUBLE_EQ(mat(0,1), mat2(0,1));
      ASSERT_DOUBLE_EQ(mat(1,1), mat2(1,1));
      ASSERT_DOUBLE_EQ(mat(2,1), mat2(2,1));
    }

    TEST(TestDnpsoup, RhoEq){
      auto spins = dnpsoup::SpinSys();
      spins.addSpin(1, dnpsoup::SpinType::e, 0.0, 0.0, 0.0, true);
      spins.addSpin(2, dnpsoup::SpinType::H, 0.7, 0.0, 0.7, true);
      auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      spins.setShielding(dnpsoup::SpinId(1), 2.00252, 2.00252, 2.00252, euler);
      //spins.setCsa(dnpsoup::SpinId(2), 0.0, 0.0, 0.0, euler);

      auto collection = spins.summarize<dnpsoup::DnpExperiment>();
      auto offset_collection = spins.summarizeOffset<dnpsoup::DnpExperiment>();
      double b0 = 9.36932070693522;
      double em_freq = 263.0e9;
      collection.setPropertyValue(dnpsoup::ValueName::b0, b0);
      offset_collection.setPropertyValue(dnpsoup::ValueName::b0, b0);
      auto spin_ids = spins.getSpinIds(dnpsoup::SpinType::e);
      for(const auto &sid : spin_ids){
        collection.setPropertyValue(
            dnpsoup::InteractionType::Shielding, sid, dnpsoup::ValueName::offset, em_freq);
        offset_collection.setPropertyValue(
            dnpsoup::InteractionType::Offset, sid, dnpsoup::ValueName::offset, em_freq);
      }
      auto euler_spins = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      auto hamiltonian = collection.genMatrix(euler_spins);
      auto hamiltonian_offset = offset_collection.genMatrix(euler_spins);
      std::cout << "hamiltonian:\n" << hamiltonian << "\n"
           << "offset hamiltonian:\n" << hamiltonian_offset << "\n";

      const double temperature = 100.0;
      auto rho_eq_lab = dnpsoup::genRhoEq(hamiltonian + hamiltonian_offset, temperature);
      auto rho_eq_relative = dnpsoup::genRhoEq(hamiltonian, temperature);
      std::cout << "Rho_equilibrium_lab =\n" << rho_eq_lab << std::endl;
      std::cout << "Rho_equilibrium =\n" << rho_eq_relative << std::endl;

      //auto rho_lab_flatten = dnpsoup::flatten(rho_eq_lab, 'c');
      //std::cout << "Rho_eq_lab flatten:\n" << rho_lab_flatten << "\n";
      //auto rho_relative_flatten = dnpsoup::flatten(rho_eq_relative, 'c');
      //std::cout << "Rho_eq_relative flatten:\n" << rho_relative_flatten << "\n";

      auto [eigenvals, eigenvec] = dnpsoup::diagonalizeMat(rho_eq_relative);
      std::cout << "Eigenvalues:\n" << eigenvals << "\n";
      std::cout << "Eigenvector:\n" << eigenvec << std::endl;
    }
}


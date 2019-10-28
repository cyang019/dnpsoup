#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <cmath>

namespace {
    using MatrixD = dnpsoup::MatrixDbl;
    using MatrixCD = dnpsoup::MatrixCxDbl;
    using cxdbl = dnpsoup::cxdbl;
    using namespace dnpsoup;
    using namespace std;

    TEST(TestDnpsoup, DipolarInteractionHH){
      auto dipoleRR = dnpsoup::DipolarInteraction<dnpsoup::RotatingFrame, dnpsoup::RotatingFrame>(
          dnpsoup::gamma_H1, dnpsoup::gamma_H1, 2, 2);
      auto dipoleRL = dnpsoup::DipolarInteraction<dnpsoup::RotatingFrame, dnpsoup::LabFrame>(
          dnpsoup::gamma_H1, dnpsoup::gamma_H1, 2, 2);

      auto p = dnpsoup::genDipoleProperty(1.0);
      auto e = dnpsoup::Euler(0,0,0);
      auto matRR = dipoleRR.genMatrix(p, e);
      std::cout << "matRR:\n" << matRR << std::endl;
      auto matRL = dipoleRL.genMatrix(p, e);
      std::cout << "matRL:\n" << matRL << std::endl;
    }
    TEST(TestDnpsoup, DipolarInteraction_eH){
      auto dipole = DipolarInteraction<dnpsoup::RotatingFrame, dnpsoup::LabFrame>(
          dnpsoup::beta_e, dnpsoup::gamma_H1, 2, 2);

      auto p = dnpsoup::genDipoleProperty(2.0);
      auto e = dnpsoup::Euler(0,0,0);
      auto mat = dipole.genMatrix(p, e);
      std::cout << "mat_e_H:\n" << mat << std::endl;
    }

    TEST(TestDnpsoup, EvolveUnderHyperfine){
      auto hf = DipolarInteraction<dnpsoup::RotatingFrame, dnpsoup::LabFrame>(
          dnpsoup::beta_e, dnpsoup::gamma_H1, 2, 2);
      auto p = genDipoleProperty(3.0);  // 3 anstrom
      auto e = Euler<>(0.0,0.1,0.0);
      auto mat_hf = hf.genMatrix(p, e);
      auto sz = kron(spin<Z>(2), identity<cxdbl>(2));
      auto iz = kron(identity<cxdbl>(2), spin<Z>(2));
      auto mat_Sz = 400e6 * sz;
      auto mat_Iz = 400e6 * iz;
      cout << "mat_hf:\n" << mat_hf << "\n";

      auto rho0 = kron(spin<Z>(2), identity<cxdbl>(2));
      auto rho = rho0;
      auto acq = kron(identity<cxdbl>(2), spin<Z>(2));

      auto ham = mat_hf + mat_Sz + mat_Iz;
      cout << "hamiltonian:\n" << ham << "\n";
      constexpr double dt = 1.0e-4;
      auto U = dnpsoup::exp((cxdbl(0,-1) * 2.0 * pi * dt) * ham);
      auto U_inv = dnpsoup::exp((cxdbl(0,1) * 2.0 * pi * dt) * ham);
      
      ASSERT_TRUE(allclose(U.adjoint(), U_inv, 1.0e-14));
      auto intensity0 = dnpsoup::projectionNorm(rho0, acq);
      cout << "rho before evolve:\n" << rho0;
      rho = U * rho0 * U_inv;
      cout << "rho after evolve:\n" << rho;
      auto intensity1 = dnpsoup::projectionNorm(rho, acq);

      auto mat_emr = 1.0e6 * kron(spin<X>(2), identity<cxdbl>(2));
      ham += mat_emr;
      U = dnpsoup::exp((cxdbl(0,-1) * 2.0 * pi * dt) * ham);
      U_inv = dnpsoup::exp((cxdbl(0,1) * 2.0 * pi * dt) * ham);

      rho = U * rho0 * U_inv;
      cout << "rho after evolve under EMR:\n" << rho << "\n";
      auto intensity2 = dnpsoup::projectionNorm(rho, acq);
      auto h_super = commutationSuperOp(ham);
      auto u_super = dnpsoup::exp((cxdbl(0, -1.0) * 2.0 * pi * dt) * h_super);
      auto rho2_super = u_super * dnpsoup::flatten(rho0, 'c');
      auto rho2 = zeros<cxdbl>(rho.nrows(), rho.ncols());
      for(size_t i = 0; i < rho2.nrows(); ++i){
        for(size_t j = 0; j < rho2.ncols(); ++j){
          rho2(i,j) = rho2_super(i * rho2.ncols() + j, 0);
        }
      }
      cout << "rho2 under superoperator:\n" << rho2 << "\n";

      auto intensity22 = dnpsoup::projectionNorm(rho2, acq);
      cout << "intensities:\n" << "intensity0: " << intensity0
           << "\nintensity1: " << intensity1
           << "\nintensity2: " << intensity2 << "\n"
           << "\nintensity22: " << intensity22 << "\n";
      ASSERT_NEAR(rho2_super(0,0).real(), rho(0,0).real(), 1.0e-8);
      ASSERT_NEAR(rho2_super(1,0).real(), rho(0,1).real(), 1.0e-8);
      ASSERT_NEAR(rho2_super(2,0).real(), rho(0,2).real(), 1.0e-8);
      ASSERT_NEAR(rho2_super(3,0).real(), rho(0,3).real(), 1.0e-8);
      ASSERT_NEAR(rho2_super(4,0).real(), rho(1,0).real(), 1.0e-8);
    }
} // namespace


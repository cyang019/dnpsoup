#include "dnpsoup.h"
#include "configure_dnpsoup.h"
#include "dnpsoup_core/constants.h"
#include "dnpsoup_core/experiment/hardware.h"
#include "dnpsoup_core/pulseseq/pulse_sequence.h"
#include "gtest/gtest.h"
#include <limits>
#include <vector>
#include <complex>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

namespace {
    using SpinSys = dnpsoup::SpinSys;
    using SpinType = dnpsoup::SpinType;
    using PulseSeq = dnpsoup::PulseSequence;
    using dnpsoup::cxdbl;
    using namespace std;
    using namespace dnpsoup;

    TEST(TestDnpsoup, EvolutionSSI){
      std::string totapol_js_file = 
        std::string(EXAMPLE_DIR) + "/spinsys/TOTAPOL.json";
      std::ifstream spinsys_stream;
	    spinsys_stream.exceptions(std::ios::failbit | std::ios::badbit);
	    spinsys_stream.open(totapol_js_file.c_str());
      dnpsoup::json spinsys_js;
	    spinsys_stream >> spinsys_js;
	    if(spinsys_js.find("spinsys") == spinsys_js.end()){
	    	throw std::runtime_error("Need 'spinsys' entry in the input json file.");
	    }
      std::istringstream spinsys_iss(spinsys_js["spinsys"].dump());
      SpinSys spinsys;
	    spinsys_iss >> spinsys;

      std::string pulseseq_js_file = 
        std::string(EXAMPLE_DIR) + "/pulseseq/cw_pulse_short.json";

      PulseSeq pseq;
      std::ifstream pseq_ss;
      pseq_ss.open(pulseseq_js_file.c_str());
      
      pseq_ss >> pseq;
      pseq_ss.close();

      auto packets = spinsys.summarize<DnpExperiment>();
      auto offset_packets = spinsys.summarizeOffset<DnpExperiment>();
      auto rpackets = spinsys.summarizeRelaxation();
      constexpr double b0 = 9.3520000000000074;
      constexpr double em_freq = 263.0e9;
      packets.setPropertyValue(ValueName::b0, b0);
      offset_packets.setPropertyValue(ValueName::b0, b0);
      auto e_spin_ids = spinsys.getSpinIds(SpinType::e);
      for(const auto &sid : e_spin_ids){
        offset_packets.setPropertyValue(
            InteractionType::Offset, sid, ValueName::offset, em_freq);
        packets.setPropertyValue(
            InteractionType::Shielding, sid, ValueName::offset, em_freq);
      }

      auto mas = dnpsoup::Euler<>(0.0, dnpsoup::magic_angle, 0.0);
      const auto spinsys_euler = dnpsoup::Euler<>(1.2935969750075629, 1.9315090471584213, 0.0);      

      MatrixCxDbl hamiltonian_offset = offset_packets.genMatrix(spinsys_euler * mas);
      auto temp_euler = dnpsoup::Euler<>(1.2935969750075631, 2.8868256652829287, 4.3982297150257113);
      constexpr double inc = 0.0000000010000000000000001;
      constexpr size_t cnt = 300000;
      MatrixCxDbl rho0_evolve = {
        {cxdbl(-0.0146845,6.1943e-14), cxdbl(1.20867e-08,3.02729e-09), cxdbl(2.03182e-05,-5.63397e-08), cxdbl(1.32193e-10,3.43585e-11), cxdbl(7.76608e-06,-8.25245e-09), cxdbl(-2.87418e-09,-7.14417e-10), cxdbl(1.38251e-09,-5.16038e-12), cxdbl(-5.02045e-13,-1.23254e-13)},
        {cxdbl(1.20867e-08,-3.02729e-09), cxdbl(-0.0147044,-1.5646e-14), cxdbl(-2.49258e-09,5.99363e-10), cxdbl(2.03135e-05,-5.63281e-08), cxdbl(-7.05009e-09,1.78759e-09), cxdbl(7.80427e-06,-8.33556e-09), cxdbl(-1.80376e-12,4.43273e-13), cxdbl(1.39069e-09,-5.19857e-12)},
        {cxdbl(2.03182e-05,5.63397e-08), cxdbl(-2.49258e-09,-5.99363e-10), cxdbl(2.86812e-06,7.1262e-14), cxdbl(1.46452e-08,3.67061e-09), cxdbl(-1.16105e-08,2.23107e-12), cxdbl(1.63178e-11,4.10232e-12), cxdbl(8.77511e-06,-9.26167e-09), cxdbl(-3.24467e-09,-8.06542e-10)},
        {cxdbl(1.32193e-10,-3.43585e-11), cxdbl(2.03135e-05,5.63281e-08), cxdbl(1.46452e-08,-3.67061e-09), cxdbl(-2.06518e-05,-3.18736e-14), cxdbl(9.05255e-11,-2.26838e-11), cxdbl(-1.18308e-08,2.23699e-12), cxdbl(-7.94146e-09,2.01357e-09), cxdbl(8.81922e-06,-9.35611e-09)},
        {cxdbl(7.76608e-06,8.25245e-09), cxdbl(-7.05009e-09,-1.78759e-09), cxdbl(-1.16105e-08,-2.231e-12), cxdbl(9.05255e-11,2.26838e-11), cxdbl(2.41212e-05,1.12669e-14), cxdbl(-1.47062e-08,-3.67916e-09), cxdbl(2.29339e-05,-6.33163e-08), cxdbl(1.47604e-10,3.83805e-11)},
        {cxdbl(-2.87418e-09,7.14417e-10), cxdbl(7.80427e-06,8.33556e-09), cxdbl(1.63178e-11,-4.10232e-12), cxdbl(-1.18308e-08,-2.2369e-12), cxdbl(-1.47062e-08,3.67916e-09), cxdbl(7.96011e-07,-1.03428e-13), cxdbl(-3.53136e-09,8.28443e-10), cxdbl(2.29312e-05,-6.33101e-08)},
        {cxdbl(1.38251e-09,5.16038e-12), cxdbl(-1.80376e-12,-4.4327e-13), cxdbl(8.77511e-06,9.26167e-09), cxdbl(-7.94146e-09,-2.01357e-09), cxdbl(2.29339e-05,6.33163e-08), cxdbl(-3.53136e-09,-8.28443e-10), cxdbl(0.0166806,7.65708e-14), cxdbl(-1.5818e-08,-3.96221e-09)},
        {cxdbl(-5.02046e-13,1.23253e-13), cxdbl(1.39069e-09,5.19857e-12), cxdbl(-3.24467e-09,8.06542e-10), cxdbl(8.81922e-06,9.35611e-09), cxdbl(1.47604e-10,-3.83805e-11), cxdbl(2.29312e-05,6.33101e-08), cxdbl(-1.5818e-08,3.96221e-09), cxdbl(0.016655,-5.77801e-14)}
      };
      auto rho0_evolve_super = ::dnpsoup::flatten(rho0_evolve, 'c');
      rho0_evolve_super = dnpsoup::DnpRunner::propagate(rho0_evolve_super,
          packets, hamiltonian_offset, rpackets,
          dnpsoup::Gyrotron(em_freq), temp_euler,
          inc, cnt, 100); 
    
      auto acq_mat = spinsys.acquireOn(dnpsoup::SpinType::H);
      auto acq_mat_super = ::dnpsoup::flatten(acq_mat, 'c');
      double result = ::dnpsoup::projectionNorm(rho0_evolve_super, acq_mat_super).real();
      std::cout << result << std::endl;
    }

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

        //rho_eq = r.propagate(rho_eq, 
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



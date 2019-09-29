#include "dnpsoup.h"
#include "gtest/gtest.h"
#include <limits>
#include <iostream>
#include <sstream>
#include <cmath>

namespace {
    using SpinSys = dnpsoup::SpinSys;
    using SpinType = dnpsoup::SpinType;

    TEST(TestDnpsoup, SpinSysCtor){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::H, 2.5, 0.0, 0.0);
      spins.addSpin(3, SpinType::H, 0.0, 2.0, 1.0);

      auto dim = spins.calcTotalDimension();
      auto dims = spins.calcDimensions();
      ASSERT_EQ(8u, dim);

      auto dims_desired = std::vector<std::size_t>({2u, 2u, 2u});
      ASSERT_TRUE(dims_desired == dims);
    }

    TEST(TestDnpsoup, SpinSysIO){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
      spins.addSpin(2, SpinType::H, 2.5, 0.0, 0.0);
      spins.addSpin(3, SpinType::H, 0.0, 2.0, 1.0);
      spins.irradiateOn(SpinType::e);
      spins.acquireOn(SpinType::H);

      std::cout << spins << std::endl;
      std::ostringstream oss;
      oss << spins;
      std::istringstream iss(oss.str());

      SpinSys spins2;
      iss >> spins2;
      std::cout << spins2 << std::endl;
    }

    TEST(TestDnpsoup, SpinSysHamiltonianEigen){
      auto spins = SpinSys();
      spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0, true);
      spins.addSpin(2, SpinType::H, 0.7, 0.0, 0.7, true);
      //auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      //spins.setShielding(dnpsoup::SpinId(1), 2.00252, 2.00252, 2.00252, euler);
      //spins.setCsa(dnpsoup::SpinId(2), 0.0, 0.0, 0.0, euler);

      auto collection = spins.summarize<dnpsoup::DnpExperiment>();
      double b0 = 9.36932070693522;
      double em_freq = 263.0e9;
      collection.setPropertyValue(dnpsoup::ValueName::b0, b0);
      auto spin_ids = spins.getSpinIds(SpinType::e);
      for(const auto &sid : spin_ids){
        collection.setPropertyValue(
            dnpsoup::InteractionType::Shielding, sid, dnpsoup::ValueName::offset, em_freq);
      }
      auto euler_spins = dnpsoup::Euler<>(0.0, 0.0, 0.0);
      std::cout << "components in packets:\n";
      auto oids = collection.getObservableIds();
      for(const auto &oid : oids){
        std::cout << oid.get() << "\n";
        const auto& packet = collection.getPacket(oid);
        std::cout << packet.genMatrix() << "\n";
      }

      auto hamiltonian = collection.genMatrix(euler_spins);
      std::cout << "Hamiltonian of spins e-H pair:\n" << hamiltonian << std::endl;
    }
} // namespace dnpsoup


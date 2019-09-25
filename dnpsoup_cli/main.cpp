#include "engine.h"
#include <iostream>
#include <fstream>

using namespace std;

// inline void calcSolidEffect() {
//     // // 1.39962450412e10
//     // constexpr double beta_e = -mu_b/h; ///< in Hz/T
//     // constexpr double gamma_H1 = 42.57747892e6;  // in Hz/T, actually gamma_H1/2pi
//     auto spins = SpinSys();
//     spins.addSpin(1, SpinType::e, 0.0, 0.0, 0.0);
//     spins.addSpin(2, SpinType::H, 2.0, 0.0, 0.0);
//     spins.irradiateOn(SpinType::e);
//     spins.acquireOn(SpinType::H);
//     spins.setShielding(dnpsoup::SpinId(1), 2.00263, 2.00259, 2.00234, 
//         dnpsoup::Euler<>(0.0,0.0,0.0));
//     spins.setT1(1, 1.0e-3);
//     spins.setT1(2, 2.0);
//     spins.setT2(1, 2.0e-6);
//     spins.setT2(2, 1.0e-3);
// 
//     dnpsoup::PulseSequence p(1.0e-9);
//     dnpsoup::pulseseq::EMRadiation emr(0.85e6, 0.0, 0.0);
//     dnpsoup::pulseseq::Component c;
//     c.insert({SpinType::e, emr});
//     p.set("emr", c);
//     // [number] x increment duration
//     auto uptr_sec = std::make_unique<dnpsoup::pulseseq::Pulse>(1000, "emr");
//     p.set("cw", std::move(uptr_sec));
//     std::vector<std::string> seq_names = { "cw" };
//     p.set(seq_names);
//     std::ostringstream buffer;
//     buffer << p;
// 
//     std::vector<dnpsoup::Magnet> fields;
//     for(int i = 0; i < 200; ++i){
//       fields.push_back(dnpsoup::Magnet(9.385 + static_cast<double>(i) * 0.000145));
//     }
//     // Hz
//     auto gyrotron = dnpsoup::Gyrotron(263.46e9);
//     // MAS, Temperature
//     auto probe = dnpsoup::Probe(8000.0, 100.0);
// 
//     dnpsoup::DnpRunner runner;
//     auto eulers = dnpsoup::getZCWAngles(3); 
//     auto res = runner.calcFieldProfile(
//         fields, gyrotron, probe, spins, 
//         buffer.str(), 
//         SpinType::H, eulers);
//     std::cout << "Field Profile with DNP: ";
//     for(auto f : res){
//       std::cout << f << ", ";
//     }
//     std::cout << std::endl;
// 
//     c[SpinType::e].freq = 0.0;
//     p.set("emr", c);
//     std::ostringstream buffer2;
//     buffer2 << p;
//     auto res2 = runner.calcFieldProfile(
//         fields, gyrotron, probe, spins, buffer2.str(),
//         SpinType::H, eulers);
//     std::cout << "Field Profile without radiation: ";
//     for(auto f : res2){
//       std::cout << f << ", ";
//     }
//     std::cout << std::endl;
// }

int main(int argc, char **argv)
{
	if(argc != 5){
		std::cout << "Saw " << argc - 1 << " arguments." << std::endl;
		std::cout << "Need exactly 4 argumments: "
							<< "spinsys_filename, pulse_sequence_filename, "
							<< "experiment_filename, and result_filename." << std::endl;
		return 1;
	}
	try{
		dnpsoup_exec(argv[1], argv[2], argv[3], argv[4]);
	}
	catch(const exception &e){
		std::cout << e.what() << std::endl;
		return 1;
	}
  return 0;
}

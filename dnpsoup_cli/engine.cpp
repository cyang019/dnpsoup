#include "engine.h"
#include "dnpsoup.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace dnpsoup;


void dnpsoup_exec(const std::string &spinsys_filename,
                  const std::string &pulse_sequence_filename,
								  const std::string &experiment_filename,
									const std::string &result_filename)
{
  std::ifstream spinsys_stream;
	spinsys_stream.exceptions(std::ios::failbit | std::ios::badbit);
	spinsys_stream.open(spinsys_filename.c_str());
	json spinsys_js;
	spinsys_stream >> spinsys_js;
	if(spinsys_js.find("spinsys") == spinsys_js.end()){
		throw runtime_error("Need 'spinsys' entry in the input json file.");
	}
	istringstream spinsys_iss(spinsys_js["spinsys"].dump());

  SpinSys spinsys;
	spinsys_iss >> spinsys;
	cout << "spinsys file loaded: " << spinsys_filename << endl;

  std::ifstream pulseseq_stream;
	pulseseq_stream.exceptions(std::ios::failbit | std::ios::badbit);
	pulseseq_stream.open(pulse_sequence_filename.c_str());
	ostringstream oss;
	oss << pulseseq_stream.rdbuf();
	string pulse_seq_str = oss.str();
	cout << "pulse sequence file loaded: " << pulse_sequence_filename << endl;

	std::ifstream exp_stream;
	exp_stream.exceptions(std::ios::failbit | std::ios::badbit);
	exp_stream.open(experiment_filename.c_str());
	json j;
	exp_stream >> j;
	cout << "experiment file loaded: " << experiment_filename << endl;
	
	if(j.find("task") == j.end()){
		throw std::runtime_error("Cannot find task in the input json.");
	}
	string task_str = j["task"].get<string>();
	cout << "task: " << task_str << endl;
	auto probe = Probe(j);
	if(j.find("acq") == j.end()){
		throw std::runtime_error("Cannot find acq in the input json.");
	}
	cout << "Probe loaded..." << endl;
	
	DnpRunner runner;
	if(task_str == "EigenValues") {
		auto magnet = Magnet(j);
		cout << "Magnet loaded..." << endl;
		auto gyrotron = Gyrotron(j);
		cout << "Gyrotron loaded..." << endl;
		auto euler = Euler<>(0.0, 0.0, 0.0);
		if(j.find("euler") != j.end()){
			if(j["euler"].find("alpha") != j["euler"].end()){
				euler.alpha(j["euler"]["alpha"].get<double>());
			}
			if(j["euler"].find("beta") != j["euler"].end()){
				euler.beta(j["euler"]["beta"].get<double>());
			}
			if(j["euler"].find("gamma") != j["euler"].end()){
				euler.gamma(j["euler"]["gamma"].get<double>());
			}
			cout << "Euler angle loaded..." << endl;
		}

		auto result = runner.calcEigenValues(magnet, gyrotron, probe,
			spinsys, pulse_seq_str, euler);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "Eigen Values:\n";
		for(const auto &row : result){
			for(size_t i = 0; i+1 < row.size(); ++i){
				result_stream << row[i] << ", ";
			}
			if(row.size() > 0){
				result_stream << row.back() << "\n";
			}
		}
		return;
	}

	SpinType acq_t = toSpinType(j["acq"].get<string>());
	if(task_str == "Intensity"){
		auto magnet = Magnet(j);
		cout << "Magnet loaded..." << endl;
		auto gyrotron = Gyrotron(j);
		cout << "Gyrotron loaded..." << endl;
		auto euler = Euler<>(0.0, 0.0, 0.0);
		if(j.find("euler") != j.end()){
			if(j["euler"].find("alpha") != j["euler"].end()){
				euler.alpha(j["euler"]["alpha"].get<double>());
			}
			if(j["euler"].find("beta") != j["euler"].end()){
				euler.beta(j["euler"]["beta"].get<double>());
			}
			if(j["euler"].find("gamma") != j["euler"].end()){
				euler.gamma(j["euler"]["gamma"].get<double>());
			}
			cout << "Euler angle loaded..." << endl;
		}

		auto result = runner.calcIntensity(magnet, gyrotron, probe,
			spinsys, pulse_seq_str, acq_t, euler);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "Intensity: " << result;
		return;
	}

	auto eulers = vector<Euler<>>();
  if(j.find("euler_scheme") != j.end()){
    if(j["euler_scheme"].find("zcw") != j["euler_scheme"].end()){
      auto zcw_input = j["euler_scheme"]["zcw"].get<uint64_t>();
      eulers = getZCWAngles(zcw_input);
      cout << "ZCW " << eulers.size() << " angles loaded..." << endl;
    }
    else{
      throw runtime_error("only 'zcw' euler_scheme is supported at the momonet. Please try to put euler angles explicitly if using other schemes.");
    }
  }
  else if(j.find("eulers") != j.end()){
		for(const auto &euler_js : j["eulers"]){
        auto euler = dnpsoup::Euler<>(0.0, 0.0, 0.0);
        if(euler_js.find("alpha") != euler_js.end()){
          euler.alpha(euler_js["alpha"].get<double>());
        }
        if(euler_js.find("beta") != euler_js.end()){
          euler.beta(euler_js["beta"].get<double>());
        }
        if(euler_js.find("gamma") != euler_js.end()){
          euler.gamma(euler_js["gamma"].get<double>());
        }
        eulers.push_back(euler);
		}
	}
	int ncores = 1;
	if(j.find("ncores") != j.end()){
		ncores = j["ncores"].get<int>();
	}
	if(task_str == "PowderIntensity"){
		auto magnet = Magnet(j);		
		auto gyrotron = Gyrotron(j);
		auto result = runner.calcPowderIntensity(magnet, gyrotron, probe,
			spinsys, pulse_seq_str, acq_t, eulers, ncores);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "PowderIntensity: " << result;
		return;
	}

	if(task_str == "FieldProfile"){
		std::vector<double> result;
		if(j.find("fields") != j.end()){
			vector<Magnet> magnets;
			for(const auto &f_val : j["fields"]) {
				auto temp_val = f_val.get<double>();
				magnets.push_back(Magnet(temp_val));
			}
			auto gyrotron = Gyrotron(j);
		  result = runner.calcFieldProfile(magnets, gyrotron, probe,
		  	spinsys, pulse_seq_str, acq_t, eulers, ncores);
		}
		else if(j.find("emrs") != j.end()){
			vector<Gyrotron> emrs;
			for(const auto &freq_val : j["emrs"]){
				auto temp_val = freq_val.get<double>();
				emrs.push_back(Gyrotron(temp_val));
			}
			auto magnet = Magnet(j);
		  result = runner.calcFieldProfile(magnet, emrs, probe,
		  	spinsys, pulse_seq_str, acq_t, eulers, ncores);
		}
		else {
			throw runtime_error("Neither 'fields' nor 'emrs' was in the input json.");
		}
    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "Field Profile:\n";
		for(size_t i = 0; i + 1 < result.size(); ++i){
			result_stream << result[i] << ", ";
		}
		if(result.size() > 0){
			result_stream << result.back() << "\n";
		}
		return;
	}

	throw runtime_error("Accepted tasks: 'EigenValues', 'Intensity', 'PowderIntensity', 'FieldProfile'");
	return;
}
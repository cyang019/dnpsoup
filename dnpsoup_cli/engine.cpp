#include "engine.h"
#include "dnpsoup.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <limits>
#include <cmath>
#include <thread>

using namespace std;
using namespace dnpsoup;


extern "C"
void dnpsoup_exec(const std::string &spinsys_filename,
                  const std::string &pulse_sequence_filename,
								  const std::string &experiment_filename,
									const std::string &result_filename)
{
  auto start_time = chrono::high_resolution_clock::now();
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

		auto result = DnpRunner::calcEigenValues(magnet, gyrotron, probe,
			spinsys, pulse_seq_str, euler);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# Eigen Values:\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
		for(const auto &row : result){
			for(size_t i = 0; i+1 < row.size(); ++i){
				result_stream << row[i] << ",";
			}
			if(row.size() > 0){
				result_stream << row.back() << "\n";
			}
		}
    auto end_time = chrono::high_resolution_clock::now();
	  cout << "Total time: " 
         << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() 
         << " seconds." << endl;
		return;
	}

	SpinType acq_t = toSpinType(j["acq"].get<string>());
	if(task_str == "Intensity" || task_str == "BuildUp"){
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

    if(task_str == "Intensity"){
		  auto result = DnpRunner::calcIntensity(magnet, gyrotron, probe,
		  	spinsys, pulse_seq_str, acq_t, euler);

      std::ofstream result_stream;
	    result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	    result_stream.open(result_filename.c_str());
      result_stream << setprecision(numeric_limits<double>::max_digits10);
		  result_stream << "# Intensity:\n" << result;
      auto end_time = chrono::high_resolution_clock::now();
	    cout << "Total time: " 
           << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() 
           << " seconds." << endl;
		  return;
    } else {
		  auto results = DnpRunner::calcBuildUp(magnet, gyrotron, probe,
		  	spinsys, pulse_seq_str, acq_t, euler);

      std::ofstream result_stream;
	    result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	    result_stream.open(result_filename.c_str());
		  result_stream << "# BuildUp:\n";
		  result_stream << "# time,intensity\n";
      result_stream << setprecision(numeric_limits<double>::max_digits10);
      for(const auto &val_pair : results){
        result_stream << val_pair.first << "," << val_pair.second << "\n";
      }
      auto end_time = chrono::high_resolution_clock::now();
	    cout << "Total time: " 
           << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() 
           << " seconds." << endl;
		  return;
    }
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
    cout << std::thread::hardware_concurrency() << " cores detected on the current machine." << endl;
		ncores = j["ncores"].get<int>();
    if (ncores > (int)std::thread::hardware_concurrency()){
      ncores = std::thread::hardware_concurrency();
    }
    cout << "use " << ncores << " core(s)." << endl;
	}
	if(task_str == "PowderIntensity"){
		auto magnet = Magnet(j);		
		auto gyrotron = Gyrotron(j);
		auto result = DnpRunner::calcPowderIntensity(magnet, gyrotron, probe,
			spinsys, pulse_seq_str, acq_t, eulers, ncores);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
    result_stream << setprecision(numeric_limits<double>::max_digits10);
		result_stream << "# PowderIntensity:\n" << result;
    auto end_time = chrono::high_resolution_clock::now();
	  cout << "Total time: " 
         << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() 
         << " seconds." << endl;
		return;
	}

  if(task_str == "PowderBuildUp"){
		auto magnet = Magnet(j);		
		auto gyrotron = Gyrotron(j);
		auto results = DnpRunner::calcPowderBuildUp(magnet, gyrotron, probe,
			spinsys, pulse_seq_str, acq_t, eulers, ncores);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# PowderBuildUp:\n";
		result_stream << "# time,intensity\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
    for(const auto &val_pair : results){
      result_stream << val_pair.first << "," << val_pair.second << "\n";
    }
    auto end_time = chrono::high_resolution_clock::now();
	  cout << "Total time: " 
         << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() 
         << " seconds." << endl;
    return;
  }

  if(task_str == "PowderBuildUpEnhancement"){
		auto magnet = Magnet(j);		
		auto gyrotron = Gyrotron(j);
		auto results = DnpRunner::calcPowderBuildUpEnhancement(magnet, gyrotron, probe,
			spinsys, pulse_seq_str, acq_t, eulers, ncores);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# PowderBuildUpEnhancement:\n";
		result_stream << "# time,intensity\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
    for(const auto &val_pair : results){
      result_stream << val_pair.first << "," << val_pair.second << "\n";
    }
    auto end_time = chrono::high_resolution_clock::now();
	  cout << "Total time: " 
         << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() 
         << " seconds." << endl;
    return;
  }

	if(task_str == "FieldProfile"){
		std::vector<std::pair<double, double>> result;
		if(j.find("fields") != j.end()){
			vector<Magnet> magnets;
			for(const auto &f_val : j["fields"]) {
				auto temp_val = f_val.get<double>();
				magnets.push_back(Magnet(temp_val));
			}
			auto gyrotron = Gyrotron(j);
		  result = DnpRunner::calcFieldProfile(magnets, gyrotron, probe,
		  	spinsys, pulse_seq_str, acq_t, eulers, ncores);
		}
    else if(j.find("field range") != j.end()){
      vector<Magnet> magnets;
      double f_beg = j["field range"]["begin"].get<double>();
      double f_end = j["field range"]["end"].get<double>();
      double f_step = j["field range"]["step"].get<double>();
      if(std::abs(f_step) < eps){
        f_step = (f_end - f_beg)/99;
      }
      double f_temp = f_beg;
      while(f_temp < f_end + eps){
        magnets.push_back(Magnet(f_temp));
        f_temp += f_step;
      }
			auto gyrotron = Gyrotron(j);
		  result = DnpRunner::calcFieldProfile(magnets, gyrotron, probe,
		  	spinsys, pulse_seq_str, acq_t, eulers, ncores);
    }
		else if(j.find("emrs") != j.end()){
			vector<Gyrotron> emrs;
			for(const auto &freq_val : j["emrs"]){
				auto temp_val = freq_val.get<double>();
				emrs.push_back(Gyrotron(temp_val));
			}
			auto magnet = Magnet(j);
		  result = DnpRunner::calcFieldProfile(magnet, emrs, probe,
		  	spinsys, pulse_seq_str, acq_t, eulers, ncores);
		}
    else if(j.find("emr range") != j.end()){
      vector<Gyrotron> emrs;
      double emr_beg = j["emr range"]["begin"].get<double>();
      double emr_end = j["emr range"]["end"].get<double>();
      double emr_step = j["emr range"]["step"].get<double>();
      if(std::abs(emr_step) < eps){
        emr_step = (emr_end - emr_beg)/99;
      }
      double emr_temp = emr_beg;
      while(emr_temp < emr_end + eps){
        emrs.push_back(Gyrotron(emr_temp)); 
        emr_temp += emr_step;
      }
			auto magnet = Magnet(j);
		  result = DnpRunner::calcFieldProfile(magnet, emrs, probe,
		  	spinsys, pulse_seq_str, acq_t, eulers, ncores);
    }
		else {
			throw runtime_error("Neither 'fields' nor 'emrs' was in the input json.");
		}
    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# Field Profile:\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
		for(size_t i = 0; i < result.size(); ++i){
			result_stream << result[i].first << "," << result[i].second << "\n";
		}
    auto end_time = chrono::high_resolution_clock::now();
	  cout << "Total time: " 
         << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() 
         << " seconds." << endl;
		return;
	}

	throw runtime_error("Accepted tasks: 'EigenValues', 'Intensity', 'PowderIntensity', 'FieldProfile', 'BuildUp'.");
	return;
}

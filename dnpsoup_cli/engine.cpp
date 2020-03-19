#include "engine.h"
#include "dnpsoup.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
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
	//ostringstream oss;
	//oss << pulseseq_stream.rdbuf();
	//string pulse_seq_str = oss.str();
  PulseSequence seq;
  //istringstream iss(pulse_seq_str);
  try{
    //iss >> seq;
    pulseseq_stream >> seq;
  }
  catch(const exception &e){
    throw PulseSequenceError(e.what());
  }
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
			spinsys, seq, euler);

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
		return;
	}

	SpinType acq_t = toSpinType(j["acq"].get<string>());
	if(task_str == "Intensity" || task_str == "BuildUp"
      || task_str == "BuildUpEnhancement"){
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
		  	spinsys, seq, acq_t, euler);

      std::ofstream result_stream;
	    result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	    result_stream.open(result_filename.c_str());
      result_stream << setprecision(numeric_limits<double>::max_digits10);
		  result_stream << "# Intensity:\n" << result;
		  return;
    } else {  // "BuildUp" or "BuildUpEnhancement"
      vector<pair<double, double>> results;
      if(task_str == "BuildUp"){
		    results = DnpRunner::calcBuildUp(
            magnet, gyrotron, probe,
		  	    spinsys, seq, acq_t, euler);
      }
      else {  // BuildUpEnhancement
		    results = DnpRunner::calcBuildUpEnhancement(
            magnet, gyrotron, probe,
		  	    spinsys, seq, acq_t, euler);
      }

      std::ofstream result_stream;
	    result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	    result_stream.open(result_filename.c_str());
		  result_stream << "# BuildUp:\n";
		  result_stream << "# time,intensity\n";
      result_stream << setprecision(numeric_limits<double>::max_digits10);
      for(const auto &val_pair : results){
        result_stream << val_pair.first << "," << val_pair.second << "\n";
      }
		  return;
    }
	}

  /// otherwise simulating powder
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
    cout << std::thread::hardware_concurrency() << " cores detected on the current machine." << "\n";
		ncores = j["ncores"].get<int>();
    if (ncores > (int)std::thread::hardware_concurrency()){
      cout << ncores << " exceeds the number of cores on the current machine, overriding..." << endl;
      ncores = std::thread::hardware_concurrency();
    }
    cout << "use " << ncores << " core(s)." << endl;
	}
	if(task_str == "PowderIntensity"){
		auto magnet = Magnet(j);		
		auto gyrotron = Gyrotron(j);
		auto result = DnpRunner::calcPowderIntensity(magnet, gyrotron, probe,
			spinsys, seq, acq_t, eulers, ncores);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
    result_stream << setprecision(numeric_limits<double>::max_digits10);
		result_stream << "# PowderIntensity:\n" << result;
		return;
	}
  else if(task_str == "PowderBuildUp"){
		auto magnet = Magnet(j);		
		auto gyrotron = Gyrotron(j);
		auto results = DnpRunner::calcPowderBuildUp(magnet, gyrotron, probe,
			spinsys, seq, acq_t, eulers, ncores);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# PowderBuildUp:\n";
		result_stream << "# time,intensity\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
    for(const auto &val_pair : results){
      result_stream << val_pair.first << "," << val_pair.second << "\n";
    }
    return;
  }
  else if(task_str == "PowderBuildUpEnhancement"){
		auto magnet = Magnet(j);		
		auto gyrotron = Gyrotron(j);
		auto results = DnpRunner::calcPowderBuildUpEnhancement(magnet, gyrotron, probe,
			spinsys, seq, acq_t, eulers, ncores);

    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# PowderBuildUpEnhancement:\n";
		result_stream << "# time,intensity\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
    for(const auto &val_pair : results){
      result_stream << val_pair.first << "," << val_pair.second << "\n";
    }
    return;
  }
	else if(task_str == "FieldProfile"){
		std::vector<std::pair<double, double>> result;
		if(j.find("fields") != j.end()){
			vector<Magnet> magnets;
			for(const auto &f_val : j["fields"]) {
				auto temp_val = f_val.get<double>();
				magnets.push_back(Magnet(temp_val));
			}
			auto gyrotron = Gyrotron(j);
		  result = DnpRunner::calcFieldProfile(magnets, gyrotron, probe,
		  	spinsys, seq, acq_t, eulers, ncores);
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
		  	spinsys, seq, acq_t, eulers, ncores);
    }
		else if(j.find("emrs") != j.end()){
			vector<Gyrotron> emrs;
			for(const auto &freq_val : j["emrs"]){
				auto temp_val = freq_val.get<double>();
				emrs.push_back(Gyrotron(temp_val));
			}
			auto magnet = Magnet(j);
		  result = DnpRunner::calcFieldProfile(magnet, emrs, probe,
		  	spinsys, seq, acq_t, eulers, ncores);
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
		  	spinsys, seq, acq_t, eulers, ncores);
    }
		else {
			throw runtime_error("Missing 'fields', 'emrs' or 'field range', 'emr range' in the input json.");
		}
    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# Field Profile:\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
		for(size_t i = 0; i < result.size(); ++i){
			result_stream << result[i].first << "," << result[i].second << "\n";
		}
		return;
	} // field profile
  else if (task_str == "scan1d") {
		auto gyrotron = Gyrotron(j);
		auto magnet = Magnet(j);
    auto params = Parameters(magnet, gyrotron, probe,
        spinsys, seq, acq_t, eulers);
    auto j_task = j["task details"];
    if(j_task.find("name") == j_task.end()){
      throw runtime_error("Need a 'name' for scan1d task details.");
    }
    auto name = j_task["name"].get<string>();
    if(j_task.find("type") == j_task.end()){
      throw runtime_error("Need a type for scan1d tasks in 'task details'.");
    }
    auto selector = Selector();
    auto scan_type = getScanType(j_task["type"].get<string>());
    if(j_task.find("range") == j_task.end()){
      throw runtime_error("Need a 'range' for task details.");
    }
    auto task_range_js = j_task["range"];
    Range range;
    if(task_range_js.find("begin") == task_range_js.end()
        || task_range_js.find("end") == task_range_js.end()
        || task_range_js.find("step") == task_range_js.end()){
      throw runtime_error("Need 'begin', 'end', and 'step' for scan1d task 'range'.");
    }
    switch (scan_type){
      case ScanType::EmrGammaB1Type:
      case ScanType::EmrPhaseType:
        {
          if(j_task.find("spin") == j_task.end()){
            throw runtime_error("Need to set 'spin' for scan1d task details.");
          }
          auto spin_t = toSpinType(j_task["spin"].get<string>());
          selector = Selector(scan_type, name, spin_t);
          range = genRangeFromJs<double>(task_range_js);
        }
        break;
      case ScanType::EmrLengthType:
        {
          selector = Selector(scan_type, name);
          range = genRangeFromJs<std::uint64_t>(task_range_js);
        }
        break;
      case ScanType::DefaultType:
        throw runtime_error("Scan type unknown.");
    }

    auto result = ScanResults1D();
    
    result = scan1d(params, selector, range, ncores);
    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# scan1d:\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
		for(size_t i = 0; i < result.size(); ++i){
			result_stream << result[i].first << "," << result[i].second << "\n";
		}
		return;
  }
  else if (task_str == "scan2d") {
		auto gyrotron = Gyrotron(j);
		auto magnet = Magnet(j);
    auto params = Parameters(magnet, gyrotron, probe,
        spinsys, seq, acq_t, eulers);
    auto j_task = j["task details"];
    if(j_task.find("range1") == j_task.end()){
      throw runtime_error("Need a 'range1' for scan2d 'task details'.");
    }

    if(j_task.find("name1") == j_task.end()){
      throw runtime_error("Need a 'name1' for scan2d task details.");
    }
    auto name1 = j_task["name1"].get<string>();
    if(j_task.find("name2") == j_task.end()){
      throw runtime_error("Need a 'name2' for scan2d task details.");
    }
    auto name2 = j_task["name2"].get<string>();
    if(j_task.find("type1") == j_task.end()){
      throw runtime_error("Need a type1 for scan2d tasks in 'task details'.");
    }
    auto scan_type1 = getScanType(j_task["type1"].get<string>());
    if(j_task.find("range1") == j_task.end()){
      throw runtime_error("Need a 'range1' for scan2d task details.");
    }
    auto task_range1_js = j_task["range1"];
    if(task_range1_js.find("begin") == task_range1_js.end()
        || task_range1_js.find("end") == task_range1_js.end()
        || task_range1_js.find("step") == task_range1_js.end()){
      throw runtime_error("Need 'begin', 'end' and 'step' for scan2d task 'range1'.");
    }
    Selector selector1;
    Range range1;
    switch (scan_type1){
      case ScanType::EmrGammaB1Type:
      case ScanType::EmrPhaseType:
        {
          if(j_task.find("spin1") == j_task.end()){
            throw runtime_error("Need to set 'spin' for scan1d task details.");
          }
          auto spin_t = toSpinType(j_task["spin1"].get<string>());
          selector1 = Selector(scan_type1, name1, spin_t);
          range1 = genRangeFromJs<double>(task_range1_js);
        }
        break;
      case ScanType::EmrLengthType:
        {
          selector1 = Selector(scan_type1, name1);
          range1 = genRangeFromJs<uint64_t>(task_range1_js);
        }
        break;
      case ScanType::DefaultType:
        throw runtime_error("Scan type1 unknown.");
    }
    if(j_task.find("type2") == j_task.end()){
      throw runtime_error("Need a type2 for scan2d tasks in 'task details'.");
    }
    auto scan_type2 = getScanType(j_task["type2"].get<string>());
    if(j_task.find("range2") == j_task.end()){
      throw runtime_error("Need a 'range2' for scan2d 'task details'.");
    }
    auto task_range2_js = j_task["range2"];
    if(task_range2_js.find("begin") == task_range2_js.end()
        || task_range2_js.find("end") == task_range2_js.end()
        || task_range2_js.find("step") == task_range2_js.end()){
      throw runtime_error("Need 'begin', 'end' and 'step' for scan2d task 'range2'.");
    }
    Selector selector2;
    Range range2;
    switch (scan_type2){
      case ScanType::EmrGammaB1Type:
      case ScanType::EmrPhaseType:
        {
          if(j_task.find("spin2") == j_task.end()){
            throw runtime_error("Need to set 'spin' for scan1d task details.");
          }
          auto spin_t = toSpinType(j_task["spin2"].get<string>());
          selector2 = Selector(scan_type2, name2, spin_t);
          range2 = genRangeFromJs<double>(task_range2_js);
        }
        break;
      case ScanType::EmrLengthType:
        {
          selector2 = Selector(scan_type2, name2);
          range2 = genRangeFromJs<uint64_t>(task_range2_js);
        }
        break;
      case ScanType::DefaultType:
        throw runtime_error("Scan type2 unknown.");
    }

    auto result = ScanResults2D();
    result = scan2d(params, selector1, range1, selector2, range2, ncores);
    
    std::ofstream result_stream;
	  result_stream.exceptions(std::ios::failbit | std::ios::badbit);
	  result_stream.open(result_filename.c_str());
		result_stream << "# Scan2d:\n";
    result_stream << setprecision(numeric_limits<double>::max_digits10);
		for(size_t i = 0; i < result.size(); ++i){
      const auto [x, y, z] = result[i];
      result_stream << x << "," << y << "," << z << "\n";
		}
		return;
  }

  const string err_msg = "Accepted tasks: "
    "'EigenValues', 'Intensity', 'PowderIntensity', 'FieldProfile', "
    "'BuildUp', 'PowderBuildUpEnhancement', "
    "'scan1d', and 'scan2d'";
	throw runtime_error(err_msg);
	return;
}

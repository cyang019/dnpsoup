#include "engine.h"
#include "engine_impl.h"
#include "json.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

using namespace std;
using json = nlohmann::json;


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
	cout << "spinsys file loaded: " << spinsys_filename << endl;

  std::ifstream pulseseq_stream;
	pulseseq_stream.exceptions(std::ios::failbit | std::ios::badbit);
	pulseseq_stream.open(pulse_sequence_filename.c_str());
  json pulseseq_js;
  pulseseq_stream >> pulseseq_js;
	cout << "pulse sequence file loaded: " << pulse_sequence_filename << endl;

	std::ifstream exp_stream;
	exp_stream.exceptions(std::ios::failbit | std::ios::badbit);
	exp_stream.open(experiment_filename.c_str());
	json exp_js;
	exp_stream >> exp_js;
	cout << "experiment file loaded: " << experiment_filename << endl;
  dnpsoup_exec_internal(spinsys_js, pulseseq_js, exp_js, result_filename);
}


extern "C"
void dnpsoup_exec0(const std::string &config_filename,
                  const std::string &result_filename)
{
  std::ifstream config_stream;
	config_stream.exceptions(std::ios::failbit | std::ios::badbit);
	config_stream.open(config_filename.c_str());
	json config_js;
  config_stream >> config_js; 
  if(config_js.find("spinsys") == config_js.end()){
    throw std::runtime_error("Cannot find spinsys configuration");
  }
  if(config_js.find("pulseseq") == config_js.end()){
    throw std::runtime_error("Cannot find pulseseq configuration");
  }
  if(config_js.find("settings") == config_js.end()){
    throw std::runtime_error("Cannot find settings configuration");
  }
  json spinsys_js = config_js["spinsys"];
	cout << "spinsys loaded... " << endl;
  json pulseseq_js = config_js["pulseseq"];
  cout << "pulse sequence loaded..." << endl;
  json exp_js = config_js["settings"];
  cout << "experiment loaded..." << endl;
  dnpsoup_exec_internal(spinsys_js, pulseseq_js, exp_js, result_filename);
}

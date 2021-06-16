// dnpsoup, a dynamic nuclear polarization simulation program.
#include "configure_dnpsoup.h"
#include "engine.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>
#ifndef APPLE
  #include <filesystem>
#endif

using namespace std;

bool fileExists(const string &filename)
{
  ifstream f(filename.c_str());
  return f.good();
}

int main(int argc, char **argv)
{
  const string authors = "DNPSOUP authors";
  const string license_short_line = " Copyright (C) " + 
    to_string(LICENSE_START) + "-" + to_string(LICENSE_END) + " " +
    authors;
  std::cout << PROJECT_NAME << " v" << PROJECT_VER << license_short_line << std::endl;
  int ret_code = 0;
	if(argc != 5 && argc != 3 && argc != 2){
		std::cout << "Saw " << argc - 1 << " arguments." << std::endl;
    for(int i = 0; i < argc; ++i) {
      std::cout << argv[i] << " ";
    }
    std::cout << std::endl;
		std::cout << "Need 1, 2 or 4 argumments: "
              << "dnpsoup_exec [configure_file]\n"
              << "dnpsoup_exec [configure_file] [output file]\n"
              << "or\n"
              << "dnpsoup_exec [spinsys_file] [pulseseq_file] [settings_file] [output file]\n"
							<< std::endl;
		return 1;
	}

  if(argc == 2) {
	  try{
      if(!fileExists(argv[1])){
        std::cout << "cannot find " << argv[1] << std::endl;
        return 1;
      }
      
      auto start_time = chrono::high_resolution_clock::now();
      const string input_filename = argv[1];
      const string default_output_filename = "results/" + input_filename + ".result.log";
	  	dnpsoup_exec0(argv[1], default_output_filename);
      auto end_time = chrono::high_resolution_clock::now();
      auto millis = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
	    cout << "Total time: " 
           << millis.count()/1000 
           << "."
           << millis.count()%1000
           << " seconds." << endl;
	  }
	  catch(const exception &e){
	  	std::cout << e.what() << std::endl;
      ret_code = 1;
	  }
  } else if(argc == 3) {
	  try{
      if(!fileExists(argv[1])){
        std::cout << "cannot find " << argv[1] << std::endl;
        return 1;
      }
      
#ifndef APPLE
      const auto output_dir = std::filesystem::path(argv[2]).remove_filename();
      if(!output_dir.empty() && !std::filesystem::exists(output_dir)) {
        std::filesystem::create_directory(output_dir);
        cout << "create directory: " << output_dir << endl;
      }
#endif
      auto start_time = chrono::high_resolution_clock::now();
	  	dnpsoup_exec0(argv[1], argv[2]);
      auto end_time = chrono::high_resolution_clock::now();
      auto millis = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
	    cout << "Total time: " 
           << millis.count()/1000 
           << "."
           << millis.count()%1000
           << " seconds." << endl;
	  }
	  catch(const exception &e){
	  	std::cout << e.what() << std::endl;
      ret_code = 1;
	  }
  } else if(argc == 5) {
	  try{
      if(!fileExists(argv[1])){
        std::cout << "cannot find " << argv[1] << std::endl;
        return 1;
      }
      if(!fileExists(argv[2])){
        std::cout << "cannot find " << argv[2] << std::endl;
        return 1;
      }
      if(!fileExists(argv[3])){
        std::cout << "cannot find " << argv[3] << std::endl;
        return 1;
      }
      
#ifndef APPLE
      const auto output_dir = std::filesystem::path(argv[4]).remove_filename();
      if(!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directory(output_dir);
        cout << "create directory: " << output_dir << endl;
      }
#endif
      auto start_time = chrono::high_resolution_clock::now();
	  	dnpsoup_exec(argv[1], argv[2], argv[3], argv[4]);
      auto end_time = chrono::high_resolution_clock::now();
      auto millis = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
	    cout << "Total time: " 
           << millis.count()/1000 
           << "."
           << millis.count()%1000
           << " seconds." << endl;
	  }
	  catch(const exception &e){
	  	std::cout << e.what() << std::endl;
      ret_code = 1;
	  }
  }

  return ret_code;
}

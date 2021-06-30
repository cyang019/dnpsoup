// dnpsoup, a dynamic nuclear polarization simulation program.
#include "configure_dnpsoup.h"
#include "engine.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <ctime>

using namespace std;

std::string timestamp()
{
  auto current_time = chrono::system_clock::now();
  char buffer[80];
  auto transformed = current_time.time_since_epoch().count() / 1000000;
  auto millis = transformed % 1000;
  std::time_t tt;
  tt = chrono::system_clock::to_time_t(current_time);
  auto timeinfo = localtime(&tt);
  strftime(buffer, 80, "%F-%H-%M-%S", timeinfo);
  sprintf(buffer, "%s-%03d", buffer, (int)millis);
  return std::string(buffer);
}

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
      
      auto argv_2_str = std::string(argv[2]);
      auto output_filepath = std::filesystem::path(argv_2_str);
      auto output_dir = output_filepath.remove_filename();
      auto output_filename = output_filepath.filename().string();
      if (argv_2_str.size() > 0 && output_filename.size() == 0) {
        output_filename = argv_2_str;
      }
      if(std::filesystem::is_directory(output_filepath) || output_filename.size() == 0) {
        // generates a default filename if non provided
        output_dir = output_filepath;
        output_filename = "dnpsoup-" + timestamp() + ".result";
        output_filepath = output_dir / std::filesystem::path(output_filename);
      }
      if(!output_dir.empty() && !std::filesystem::exists(output_dir)) {
        // creates a directory if non existed
        std::filesystem::create_directory(output_dir);
        cout << "create directory: " << output_dir << endl;
      }
      auto start_time = chrono::high_resolution_clock::now();
      // test to see if can correctly write to a file
      {
        std::ofstream result_stream(output_filepath.string().c_str());
        result_stream.close();
        cout << "result filename: " << output_filepath << endl;
      }
	  	dnpsoup_exec0(argv[1], output_filepath.string().c_str());
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
      
      const auto output_dir = std::filesystem::path(argv[4]).remove_filename();
      if(!output_dir.empty() && !std::filesystem::exists(output_dir)) {
        std::filesystem::create_directory(output_dir);
        cout << "create directory: " << output_dir << endl;
      }
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

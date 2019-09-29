#include "engine.h"
#include <iostream>
#include <fstream>

using namespace std;


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

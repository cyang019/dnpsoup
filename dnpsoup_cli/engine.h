#ifndef DNPSOUP_ENGINE_H
#define DNPSOUP_ENGINE_H
#include <string>

void dnpsoup_exec(const std::string &spinsys_filename,
                  const std::string &pulse_sequence_filename,
								  const std::string &experiment_filename,
									const std::string &result_filename);

#endif

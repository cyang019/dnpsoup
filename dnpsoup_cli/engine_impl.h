#ifndef DNPSOUP_ENGINE_IMPL_H
#define DNPSOUP_ENGINE_IMPL_H
#include <string>
#include "json.hpp"

void dnpsoup_exec_internal(
    const nlohmann::json &spinsys_js,
    const nlohmann::json &pulseseq_js,
    const nlohmann::json &settings_js,
    const std::string &output_name);

#endif
